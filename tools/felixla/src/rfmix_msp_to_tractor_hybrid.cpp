// rfmix_msp_to_tractor_hybrid.cpp
// designed by Kai, implemented by codex
//
// Convert phased genotype VCF/BCF plus RFMix .msp.tsv local ancestry calls
// into ancestry-aware packed files for the tractor_hybrid backend.

#include <htslib/hts.h>
#include <htslib/vcf.h>

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>

#include <unistd.h>

struct RareCarrierPacked {
    uint32_t pos_index;
    uint32_t anc_hap;
};

struct CarrierHap {
    uint32_t hap_id;
    uint8_t ancestry;
};

struct AncestryState {
    std::vector<std::vector<uint64_t>> masks;
    std::vector<int8_t> hap_ancestry;
};

struct OpenAncestryBlock {
    bool active = false;
    uint32_t block_id = 0;
    int geno_rid = -1;
    std::string chr;
    int64_t start_pos = -1;
    int64_t end_pos = -1;
    AncestryState state;
};

struct MspRecord {
    bool valid = false;
    int geno_rid = -1;
    std::string chr;
    int64_t start_pos = 0;
    int64_t end_pos = 0;
    AncestryState state;
};

static constexpr uint32_t kMaxPackedHapId = (1u << 27) - 1u;
static constexpr uint64_t kProgressRecordInterval = 10000;
static constexpr int kProgressSecondsInterval = 2;

[[noreturn]] static void die(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    std::fputs("ERROR: ", stderr);
    std::vfprintf(stderr, fmt, args);
    std::fputc('\n', stderr);
    va_end(args);
    std::exit(1);
}

struct IndexRecordCount {
    bool available = false;
    uint64_t total = 0;
};

static hts_idx_t* load_any_vcf_index(const char* path) {
    hts_idx_t* idx = hts_idx_load3(path, nullptr, HTS_FMT_CSI, HTS_IDX_SILENT_FAIL);
    if (idx) return idx;
    return hts_idx_load3(path, nullptr, HTS_FMT_TBI, HTS_IDX_SILENT_FAIL);
}

static IndexRecordCount get_index_record_count(const char* path) {
    IndexRecordCount result;
    hts_idx_t* idx = load_any_vcf_index(path);
    if (!idx) return result;

    int nseq = hts_idx_nseq(idx);
    if (nseq <= 0) {
        hts_idx_destroy(idx);
        return result;
    }

    uint64_t total = 0;
    for (int tid = 0; tid < nseq; ++tid) {
        uint64_t mapped = 0;
        uint64_t unmapped = 0;
        if (hts_idx_get_stat(idx, tid, &mapped, &unmapped) != 0) {
            hts_idx_destroy(idx);
            return result;
        }

        if (mapped > std::numeric_limits<uint64_t>::max() - total) {
            total = std::numeric_limits<uint64_t>::max();
        } else {
            total += mapped;
        }
    }

    hts_idx_destroy(idx);
    result.available = true;
    result.total = total;
    return result;
}

class ProgressReporter {
public:
    ProgressReporter(const char* input_label, const char* input_path, const char* converted_label)
        : converted_label_(converted_label),
          total_(get_index_record_count(input_path)),
          last_report_(std::chrono::steady_clock::now()),
          stderr_is_tty_(isatty(fileno(stderr)) != 0) {
        if (total_.available) {
            std::fprintf(
                stderr,
                "Progress: %s index reports %llu records.\n",
                input_label,
                static_cast<unsigned long long>(total_.total)
            );
        } else {
            std::fprintf(
                stderr,
                "Progress: %s index record count unavailable; reporting scanned records only.\n",
                input_label
            );
        }
    }

    void record_scanned() {
        ++records_scanned_;
    }

    void maybe_report(
        uint64_t converted,
        uint64_t common,
        uint64_t rare,
        uint64_t msp_rows,
        bool force = false
    ) {
        auto now = std::chrono::steady_clock::now();
        bool count_due = records_scanned_ >= next_record_report_;
        bool time_due = now - last_report_ >= std::chrono::seconds(kProgressSecondsInterval);

        if (!force && !count_due && !time_due) {
            return;
        }

        print(converted, common, rare, msp_rows);
        last_report_ = now;
        printed_ = true;

        if (records_scanned_ >= next_record_report_) {
            if (records_scanned_ > std::numeric_limits<uint64_t>::max() - kProgressRecordInterval) {
                next_record_report_ = std::numeric_limits<uint64_t>::max();
            } else {
                next_record_report_ = records_scanned_ + kProgressRecordInterval;
            }
        }
    }

    void finish(uint64_t converted, uint64_t common, uint64_t rare, uint64_t msp_rows) {
        maybe_report(converted, common, rare, msp_rows, true);
        if (stderr_is_tty_ && printed_) {
            std::fputc('\n', stderr);
        }
    }

private:
    void print(uint64_t converted, uint64_t common, uint64_t rare, uint64_t msp_rows) const {
        const char* prefix = stderr_is_tty_ ? "\r" : "";
        const char* suffix = stderr_is_tty_ ? "" : "\n";

        if (total_.available && total_.total > 0) {
            uint64_t capped_records = std::min(records_scanned_, total_.total);
            double pct = 100.0 * static_cast<double>(capped_records) /
                         static_cast<double>(total_.total);
            std::fprintf(
                stderr,
                "%sProgress: records %llu/%llu (%.1f%%), MSP rows %llu, converted %llu %s, common %llu, rare %llu%s",
                prefix,
                static_cast<unsigned long long>(records_scanned_),
                static_cast<unsigned long long>(total_.total),
                pct,
                static_cast<unsigned long long>(msp_rows),
                static_cast<unsigned long long>(converted),
                converted_label_,
                static_cast<unsigned long long>(common),
                static_cast<unsigned long long>(rare),
                suffix
            );
        } else if (total_.available) {
            std::fprintf(
                stderr,
                "%sProgress: records %llu/0, MSP rows %llu, converted %llu %s, common %llu, rare %llu%s",
                prefix,
                static_cast<unsigned long long>(records_scanned_),
                static_cast<unsigned long long>(msp_rows),
                static_cast<unsigned long long>(converted),
                converted_label_,
                static_cast<unsigned long long>(common),
                static_cast<unsigned long long>(rare),
                suffix
            );
        } else {
            std::fprintf(
                stderr,
                "%sProgress: records %llu, MSP rows %llu, converted %llu %s, common %llu, rare %llu%s",
                prefix,
                static_cast<unsigned long long>(records_scanned_),
                static_cast<unsigned long long>(msp_rows),
                static_cast<unsigned long long>(converted),
                converted_label_,
                static_cast<unsigned long long>(common),
                static_cast<unsigned long long>(rare),
                suffix
            );
        }

        std::fflush(stderr);
    }

    const char* converted_label_;
    IndexRecordCount total_;
    uint64_t records_scanned_ = 0;
    uint64_t next_record_report_ = kProgressRecordInterval;
    std::chrono::steady_clock::time_point last_report_;
    bool stderr_is_tty_ = false;
    bool printed_ = false;
};

static inline uint32_t pack_anc_hap(uint32_t ancestry, uint32_t hap_id) {
    return ((ancestry & 31u) << 27) | (hap_id & 0x07ffffffu);
}

static inline void set_bit(std::vector<uint64_t>& bits, uint32_t h) {
    bits[h >> 6] |= (1ULL << (h & 63));
}

static inline bool is_gt_vector_end(int32_t value) {
    return value == bcf_int32_vector_end;
}

static bool same_masks(
    const std::vector<std::vector<uint64_t>>& a,
    const std::vector<std::vector<uint64_t>>& b
) {
    if (a.size() != b.size()) return false;

    for (size_t k = 0; k < a.size(); ++k) {
        if (a[k] != b[k]) return false;
    }

    return true;
}

static uint64_t tell_or_die(FILE* fp, const char* path) {
    off_t offset = ftello(fp);
    if (offset < 0) die("ftello failed for %s", path);
    return static_cast<uint64_t>(offset);
}

static FILE* open_output_or_die(const std::string& path, const char* mode) {
    FILE* fp = std::fopen(path.c_str(), mode);
    if (!fp) die("cannot open output file: %s", path.c_str());

    if (setvbuf(fp, nullptr, _IOFBF, 1 << 24) != 0) {
        die("setvbuf failed for output file: %s", path.c_str());
    }

    return fp;
}

static void write_exact(FILE* fp, const void* data, size_t bytes, const char* what) {
    if (bytes == 0) return;
    if (std::fwrite(data, 1, bytes, fp) != bytes) {
        die("failed writing %s", what);
    }
}

static void write_u32_le(FILE* fp, uint32_t value, const char* what) {
    uint8_t bytes[4] = {
        static_cast<uint8_t>(value),
        static_cast<uint8_t>(value >> 8),
        static_cast<uint8_t>(value >> 16),
        static_cast<uint8_t>(value >> 24)
    };
    write_exact(fp, bytes, sizeof(bytes), what);
}

static void write_u64_le(FILE* fp, uint64_t value, const char* what) {
    uint8_t bytes[8] = {
        static_cast<uint8_t>(value),
        static_cast<uint8_t>(value >> 8),
        static_cast<uint8_t>(value >> 16),
        static_cast<uint8_t>(value >> 24),
        static_cast<uint8_t>(value >> 32),
        static_cast<uint8_t>(value >> 40),
        static_cast<uint8_t>(value >> 48),
        static_cast<uint8_t>(value >> 56)
    };
    write_exact(fp, bytes, sizeof(bytes), what);
}

static void write_i64_le(FILE* fp, int64_t value, const char* what) {
    write_u64_le(fp, static_cast<uint64_t>(value), what);
}

static void write_string(FILE* fp, const std::string& value, const char* what) {
    if (value.size() > std::numeric_limits<uint32_t>::max()) {
        die("string too long while writing %s", what);
    }

    write_u32_le(fp, static_cast<uint32_t>(value.size()), what);
    write_exact(fp, value.data(), value.size(), what);
}

static void write_magic_header(FILE* fp, const char magic[8]) {
    write_exact(fp, magic, 8, "magic header");
}

static void write_common_offset_idx_record(
    FILE* fp,
    uint64_t common_index,
    uint32_t global_variant_index,
    uint64_t mks_offset,
    uint64_t geno_offset
) {
    write_u64_le(fp, common_index, "common offset index");
    write_u32_le(fp, global_variant_index, "common offset index");
    write_u64_le(fp, mks_offset, "common offset index");
    write_u64_le(fp, geno_offset, "common offset index");
}

static void write_rare_offset_idx_record(
    FILE* fp,
    uint64_t rare_index,
    uint32_t global_variant_index,
    uint64_t mks_offset,
    uint64_t carrier_offset,
    uint32_t n_carriers
) {
    write_u64_le(fp, rare_index, "rare offset index");
    write_u32_le(fp, global_variant_index, "rare offset index");
    write_u64_le(fp, mks_offset, "rare offset index");
    write_u64_le(fp, carrier_offset, "rare offset index");
    write_u32_le(fp, n_carriers, "rare offset index");
}

static void write_anc_offset_idx_record(
    FILE* fp,
    uint32_t block_id,
    uint64_t mks_offset,
    uint64_t anc_offset
) {
    write_u32_le(fp, block_id, "ancestry offset index");
    write_u64_le(fp, mks_offset, "ancestry offset index");
    write_u64_le(fp, anc_offset, "ancestry offset index");
}

static void write_common_mks_record(
    FILE* fp,
    uint64_t common_index,
    uint32_t global_variant_index,
    const char* chr,
    int64_t pos,
    const std::string& id,
    const char* ref,
    const char* alt,
    uint32_t alt_index,
    uint32_t block_id,
    uint64_t geno_offset,
    uint32_t mac
) {
    write_u64_le(fp, common_index, "common mks");
    write_u32_le(fp, global_variant_index, "common mks");
    write_string(fp, chr, "common mks chr");
    write_i64_le(fp, pos, "common mks pos");
    write_string(fp, id, "common mks id");
    write_string(fp, ref, "common mks ref");
    write_string(fp, alt, "common mks alt");
    write_u32_le(fp, alt_index, "common mks");
    write_u32_le(fp, block_id, "common mks");
    write_u64_le(fp, geno_offset, "common mks");
    write_u32_le(fp, mac, "common mks");
}

static void write_rare_mks_record(
    FILE* fp,
    uint64_t rare_index,
    uint32_t global_variant_index,
    const char* chr,
    int64_t pos,
    const std::string& id,
    const char* ref,
    const char* alt,
    uint32_t alt_index,
    uint64_t carrier_offset,
    uint32_t n_carriers,
    uint32_t mac
) {
    write_u64_le(fp, rare_index, "rare mks");
    write_u32_le(fp, global_variant_index, "rare mks");
    write_string(fp, chr, "rare mks chr");
    write_i64_le(fp, pos, "rare mks pos");
    write_string(fp, id, "rare mks id");
    write_string(fp, ref, "rare mks ref");
    write_string(fp, alt, "rare mks alt");
    write_u32_le(fp, alt_index, "rare mks");
    write_u64_le(fp, carrier_offset, "rare mks");
    write_u32_le(fp, n_carriers, "rare mks");
    write_u32_le(fp, mac, "rare mks");
}

static void write_anc_mks_record(
    FILE* fp,
    uint32_t block_id,
    const char* chr,
    int64_t start_pos,
    int64_t end_pos,
    uint64_t anc_offset
) {
    write_u32_le(fp, block_id, "ancestry mks");
    write_string(fp, chr, "ancestry mks chr");
    write_i64_le(fp, start_pos, "ancestry mks start");
    write_i64_le(fp, end_pos, "ancestry mks end");
    write_u64_le(fp, anc_offset, "ancestry mks offset");
}

static void write_sidecars(
    bcf_hdr_t* ghdr,
    const std::string& samples_path,
    const std::string& meta_path,
    const char* geno_vcf,
    const char* msp_path,
    const std::string& subpopulation_codes,
    int n_samples,
    uint64_t n_haps,
    int n_words,
    int n_ancestries,
    int rare_threshold
) {
    FILE* samples_fp = open_output_or_die(samples_path, "w");
    for (int i = 0; i < n_samples; ++i) {
        std::fprintf(samples_fp, "%s\n", ghdr->samples[i]);
    }
    std::fclose(samples_fp);

    FILE* meta_fp = open_output_or_die(meta_path, "w");
    std::fprintf(meta_fp, "format_version\t1\n");
    std::fprintf(meta_fp, "n_samples\t%d\n", n_samples);
    std::fprintf(meta_fp, "n_haps\t%llu\n", static_cast<unsigned long long>(n_haps));
    std::fprintf(meta_fp, "n_words\t%d\n", n_words);
    std::fprintf(meta_fp, "n_ancestries\t%d\n", n_ancestries);
    std::fprintf(meta_fp, "rare_threshold\t%d\n", rare_threshold);
    std::fprintf(meta_fp, "source_genotype\t%s\n", geno_vcf);
    std::fprintf(meta_fp, "source_rfmix_msp\t%s\n", msp_path);
    std::fprintf(
        meta_fp,
        "rfmix_msp_interval_note\tfirst interval includes spos; later shared epos/spos boundaries belong to the previous interval\n"
    );
    if (!subpopulation_codes.empty()) {
        std::string sanitized_codes = subpopulation_codes;
        for (char& c : sanitized_codes) {
            if (c == '\t') c = ' ';
        }
        std::fprintf(meta_fp, "rfmix_subpopulation_order_codes\t%s\n", sanitized_codes.c_str());
    }
    std::fclose(meta_fp);
}

static void write_anc_block(
    FILE* anc_fp,
    const char* anc_bin_path,
    FILE* anc_mks_fp,
    const char* anc_mks_path,
    FILE* anc_idx_fp,
    uint32_t block_id,
    const char* chr,
    int64_t start_pos,
    int64_t end_pos,
    const AncestryState& state,
    int n_ancestries,
    int n_words
) {
    uint64_t offset = tell_or_die(anc_fp, anc_bin_path);

    for (int k = 0; k < n_ancestries; ++k) {
        size_t wrote = std::fwrite(
            state.masks[k].data(),
            sizeof(uint64_t),
            static_cast<size_t>(n_words),
            anc_fp
        );

        if (wrote != static_cast<size_t>(n_words)) {
            die("failed writing ancestry block %u", block_id);
        }
    }

    uint64_t mks_offset = tell_or_die(anc_mks_fp, anc_mks_path);
    write_anc_mks_record(anc_mks_fp, block_id, chr, start_pos, end_pos, offset);
    write_anc_offset_idx_record(anc_idx_fp, block_id, mks_offset, offset);
}

static void close_open_block(
    FILE* anc_fp,
    const char* anc_bin_path,
    FILE* anc_mks_fp,
    const char* anc_mks_path,
    FILE* anc_idx_fp,
    OpenAncestryBlock& block,
    uint32_t& n_blocks_written,
    int n_ancestries,
    int n_words
) {
    if (!block.active) return;

    write_anc_block(
        anc_fp,
        anc_bin_path,
        anc_mks_fp,
        anc_mks_path,
        anc_idx_fp,
        block.block_id,
        block.chr.c_str(),
        block.start_pos,
        block.end_pos,
        block.state,
        n_ancestries,
        n_words
    );

    ++n_blocks_written;
}

static void update_open_block_from_msp(
    FILE* anc_fp,
    const char* anc_bin_path,
    FILE* anc_mks_fp,
    const char* anc_mks_path,
    FILE* anc_idx_fp,
    OpenAncestryBlock& block,
    uint32_t& n_blocks_written,
    int n_ancestries,
    int n_words,
    int geno_rid,
    const char* chr,
    int64_t interval_start,
    int64_t interval_end,
    AncestryState& new_state
) {
    if (interval_start > interval_end) return;

    if (!block.active) {
        block.active = true;
        block.block_id = n_blocks_written;
        block.geno_rid = geno_rid;
        block.chr = chr;
        block.start_pos = interval_start;
        block.end_pos = interval_end;
        block.state = std::move(new_state);
        return;
    }

    bool same_chromosome = block.geno_rid == geno_rid;
    bool touches = interval_start <= block.end_pos + 1;
    bool unchanged = same_chromosome && touches && same_masks(block.state.masks, new_state.masks);

    if (unchanged) {
        block.end_pos = std::max(block.end_pos, interval_end);
        return;
    }

    if (same_chromosome && interval_start <= block.end_pos) {
        die(
            "overlapping MSP intervals with different ancestry states at %s:%lld-%lld",
            chr,
            static_cast<long long>(interval_start),
            static_cast<long long>(interval_end)
        );
    }

    close_open_block(
        anc_fp,
        anc_bin_path,
        anc_mks_fp,
        anc_mks_path,
        anc_idx_fp,
        block,
        n_blocks_written,
        n_ancestries,
        n_words
    );

    block.active = true;
    block.block_id = n_blocks_written;
    block.geno_rid = geno_rid;
    block.chr = chr;
    block.start_pos = interval_start;
    block.end_pos = interval_end;
    block.state = std::move(new_state);
}

static std::string make_split_id(
    const char* raw_id,
    const char* chr,
    int64_t pos,
    const char* ref,
    const char* alt
) {
    if (raw_id == nullptr || std::strcmp(raw_id, ".") == 0) {
        return std::string(chr) + ":" + std::to_string(pos) + ":" + ref + ":" + alt;
    }

    return std::string(raw_id) + "_" + ref + "_" + alt;
}

static void add_alt_carrier(
    std::vector<std::vector<CarrierHap>>& carriers_by_alt,
    const AncestryState& state,
    int allele,
    uint32_t hap_id
) {
    if (allele <= 0 || allele >= static_cast<int>(carriers_by_alt.size())) return;

    int8_t ancestry = state.hap_ancestry[hap_id];
    if (ancestry < 0) return;

    carriers_by_alt[allele].push_back(
        CarrierHap{hap_id, static_cast<uint8_t>(ancestry)}
    );
}

static void validate_extra_ploidy(
    const int32_t* sample_gt,
    int ploidy,
    const char* chr,
    int64_t pos,
    int sample_index
) {
    for (int p = 2; p < ploidy; ++p) {
        int32_t gt = sample_gt[p];
        if (is_gt_vector_end(gt)) break;
        die(
            "non-diploid genotype at %s:%lld sample index %d",
            chr,
            static_cast<long long>(pos),
            sample_index
        );
    }
}

static void process_genotypes(
    bcf_hdr_t* ghdr,
    bcf1_t* grec,
    const OpenAncestryBlock& block,
    int n_samples,
    int32_t** gt_arr,
    int* ngt_arr,
    std::vector<std::vector<CarrierHap>>& carriers_by_alt
) {
    const char* chr = bcf_hdr_id2name(ghdr, grec->rid);
    int64_t pos = static_cast<int64_t>(grec->pos) + 1;

    int ngt = bcf_get_genotypes(ghdr, grec, gt_arr, ngt_arr);
    if (ngt <= 0) {
        die("missing FORMAT/GT at %s:%lld", chr, static_cast<long long>(pos));
    }

    if (ngt % n_samples != 0) {
        die("GT field length is not divisible by sample count at %s:%lld", chr, static_cast<long long>(pos));
    }

    int ploidy = ngt / n_samples;
    if (ploidy < 2) {
        die("expected diploid GT at %s:%lld", chr, static_cast<long long>(pos));
    }

    carriers_by_alt.assign(static_cast<size_t>(grec->n_allele), std::vector<CarrierHap>{});

    for (int i = 0; i < n_samples; ++i) {
        const int32_t* sample_gt = *gt_arr + static_cast<size_t>(i) * ploidy;
        int32_t g0 = sample_gt[0];
        int32_t g1 = sample_gt[1];

        validate_extra_ploidy(sample_gt, ploidy, chr, pos, i);

        if (is_gt_vector_end(g0) || is_gt_vector_end(g1) ||
            bcf_gt_is_missing(g0) || bcf_gt_is_missing(g1)) {
            die(
                "missing genotype at %s:%lld sample index %d",
                chr,
                static_cast<long long>(pos),
                i
            );
        }

        int allele0 = bcf_gt_allele(g0);
        int allele1 = bcf_gt_allele(g1);

        if (!bcf_gt_is_phased(g1)) {
            die(
                "unphased genotype at %s:%lld sample index %d",
                chr,
                static_cast<long long>(pos),
                i
            );
        }

        add_alt_carrier(
            carriers_by_alt,
            block.state,
            allele0,
            static_cast<uint32_t>(2 * i)
        );
        add_alt_carrier(
            carriers_by_alt,
            block.state,
            allele1,
            static_cast<uint32_t>(2 * i + 1)
        );
    }
}

static int read_next_record(htsFile* fp, bcf_hdr_t* hdr, bcf1_t* rec) {
    int ret = bcf_read(fp, hdr, rec);
    if (ret == 0) bcf_unpack(rec, BCF_UN_STR);
    return ret;
}

static std::string trim_copy(const std::string& value) {
    size_t first = 0;
    while (first < value.size() && std::isspace(static_cast<unsigned char>(value[first]))) {
        ++first;
    }

    size_t last = value.size();
    while (last > first && std::isspace(static_cast<unsigned char>(value[last - 1]))) {
        --last;
    }

    return value.substr(first, last - first);
}

static bool starts_with(const std::string& value, const char* prefix) {
    size_t n = std::strlen(prefix);
    return value.size() >= n && value.compare(0, n, prefix) == 0;
}

static std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> fields;
    size_t start = 0;

    while (true) {
        size_t tab = line.find('\t', start);
        if (tab == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }

        fields.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }

    return fields;
}

static std::vector<std::string> split_whitespace(const std::string& line) {
    std::vector<std::string> fields;
    size_t i = 0;

    while (i < line.size()) {
        while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
        if (i >= line.size()) break;

        size_t start = i;
        while (i < line.size() && !std::isspace(static_cast<unsigned char>(line[i]))) ++i;
        fields.push_back(line.substr(start, i - start));
    }

    return fields;
}

static std::vector<std::string> split_msp_fields(const std::string& line) {
    if (line.find('\t') != std::string::npos) {
        return split_tabs(line);
    }

    return split_whitespace(line);
}

static int64_t parse_i64(const std::string& value, const char* field_name, uint64_t line_no) {
    char* end = nullptr;
    long long parsed = std::strtoll(value.c_str(), &end, 10);
    if (!end || *end != '\0') {
        die(
            "invalid integer in MSP %s at line %llu: %s",
            field_name,
            static_cast<unsigned long long>(line_no),
            value.c_str()
        );
    }
    return static_cast<int64_t>(parsed);
}

static int parse_ancestry_code(
    const std::string& value,
    int n_ancestries,
    uint64_t line_no,
    size_t hap_index
) {
    if (value.empty() || value == "." || value == "NA" || value == "nan") {
        die(
            "missing ancestry in MSP at line %llu hap column %llu",
            static_cast<unsigned long long>(line_no),
            static_cast<unsigned long long>(hap_index)
        );
    }

    char* end = nullptr;
    long parsed = std::strtol(value.c_str(), &end, 10);
    if (!end || *end != '\0') {
        die(
            "invalid ancestry code in MSP at line %llu hap column %llu: %s",
            static_cast<unsigned long long>(line_no),
            static_cast<unsigned long long>(hap_index),
            value.c_str()
        );
    }

    if (parsed < 0 || parsed >= n_ancestries) {
        die(
            "ancestry code out of range in MSP at line %llu hap column %llu: %ld",
            static_cast<unsigned long long>(line_no),
            static_cast<unsigned long long>(hap_index),
            parsed
        );
    }

    return static_cast<int>(parsed);
}

static bool hap_label_matches(
    const std::string& label,
    const std::string& sample,
    int hap_index
) {
    std::string hap = std::to_string(hap_index);
    return label == sample + "." + hap ||
           label == sample + "_" + hap ||
           label == sample + ":" + hap;
}

static int resolve_msp_chrom_to_geno_rid(
    bcf_hdr_t* ghdr,
    const std::string& msp_chr,
    std::string& resolved_chr
) {
    int rid = bcf_hdr_name2id(ghdr, msp_chr.c_str());
    if (rid >= 0) {
        resolved_chr = bcf_hdr_id2name(ghdr, rid);
        return rid;
    }

    if (starts_with(msp_chr, "chr")) {
        std::string without_chr = msp_chr.substr(3);
        rid = bcf_hdr_name2id(ghdr, without_chr.c_str());
        if (rid >= 0) {
            resolved_chr = bcf_hdr_id2name(ghdr, rid);
            return rid;
        }
    } else {
        std::string with_chr = "chr" + msp_chr;
        rid = bcf_hdr_name2id(ghdr, with_chr.c_str());
        if (rid >= 0) {
            resolved_chr = bcf_hdr_id2name(ghdr, rid);
            return rid;
        }
    }

    resolved_chr.clear();
    return -1;
}

static int infer_n_ancestries_from_codes(const std::string& line) {
    size_t colon = line.find(':');
    if (colon == std::string::npos) return 0;

    std::string rest = line.substr(colon + 1);
    std::vector<std::string> fields = split_whitespace(rest);
    int max_code = -1;

    for (const std::string& field : fields) {
        size_t eq = field.find('=');
        if (eq == std::string::npos || eq + 1 >= field.size()) continue;

        char* end = nullptr;
        long code = std::strtol(field.c_str() + eq + 1, &end, 10);
        if (!end || *end != '\0') continue;
        if (code >= 0 && code <= std::numeric_limits<int>::max()) {
            max_code = std::max(max_code, static_cast<int>(code));
        }
    }

    return max_code + 1;
}

class MspReader {
public:
    MspReader(
        const char* path,
        bcf_hdr_t* ghdr,
        int n_samples,
        int n_ancestries,
        int n_words
    ) : path_(path),
        ghdr_(ghdr),
        n_samples_(n_samples),
        n_ancestries_(n_ancestries),
        n_words_(n_words) {
        fp_ = hts_open(path, "r");
        if (!fp_) die("cannot open RFMix MSP file: %s", path);
        read_header();
    }

    ~MspReader() {
        if (line_.s) std::free(line_.s);
        if (fp_) hts_close(fp_);
    }

    const std::string& subpopulation_codes() const {
        return subpopulation_codes_;
    }

    bool read_next(MspRecord& out) {
        while (hts_getline(fp_, '\n', &line_) >= 0) {
            ++line_no_;
            std::string raw(line_.s, line_.l);
            if (!raw.empty() && raw.back() == '\r') raw.pop_back();
            std::string line = trim_copy(raw);
            if (line.empty()) continue;
            if (line[0] == '#') continue;

            std::vector<std::string> fields = split_msp_fields(line);
            size_t expected = 6 + static_cast<size_t>(2 * n_samples_);
            if (fields.size() != expected) {
                die(
                    "MSP data line %llu has %llu fields; expected %llu (6 interval fields plus 2 haplotypes per sample)",
                    static_cast<unsigned long long>(line_no_),
                    static_cast<unsigned long long>(fields.size()),
                    static_cast<unsigned long long>(expected)
                );
            }

            std::string resolved_chr;
            int geno_rid = resolve_msp_chrom_to_geno_rid(ghdr_, fields[0], resolved_chr);
            if (geno_rid < 0) {
                if (warned_missing_contigs_.insert(fields[0]).second) {
                    std::fprintf(
                        stderr,
                        "WARNING: skipping MSP contig absent from genotype header: %s\n",
                        fields[0].c_str()
                    );
                }
                continue;
            }

            int64_t start_pos = parse_i64(fields[1], "spos", line_no_);
            int64_t end_pos = parse_i64(fields[2], "epos", line_no_);
            if (start_pos < 1 || end_pos < 1) {
                die(
                    "MSP positions must be positive at line %llu",
                    static_cast<unsigned long long>(line_no_)
                );
            }

            if (start_pos > end_pos) {
                die(
                    "invalid MSP interval at line %llu: spos %lld > epos %lld",
                    static_cast<unsigned long long>(line_no_),
                    static_cast<long long>(start_pos),
                    static_cast<long long>(end_pos)
                );
            }

            out.valid = true;
            out.geno_rid = geno_rid;
            out.chr = resolved_chr;
            out.start_pos = start_pos;
            out.end_pos = end_pos;
            out.state.masks.assign(
                static_cast<size_t>(n_ancestries_),
                std::vector<uint64_t>(static_cast<size_t>(n_words_), 0)
            );
            out.state.hap_ancestry.assign(static_cast<size_t>(2 * n_samples_), -1);

            for (size_t h = 0; h < static_cast<size_t>(2 * n_samples_); ++h) {
                int ancestry = parse_ancestry_code(fields[6 + h], n_ancestries_, line_no_, h);
                set_bit(out.state.masks[ancestry], static_cast<uint32_t>(h));
                out.state.hap_ancestry[h] = static_cast<int8_t>(ancestry);
            }

            ++records_read_;
            return true;
        }

        out = MspRecord{};
        return false;
    }

    uint64_t records_read() const {
        return records_read_;
    }

private:
    void read_header() {
        while (hts_getline(fp_, '\n', &line_) >= 0) {
            ++line_no_;
            std::string raw(line_.s, line_.l);
            if (!raw.empty() && raw.back() == '\r') raw.pop_back();
            std::string line = trim_copy(raw);
            if (line.empty()) continue;

            if (starts_with(line, "#Subpopulation order/codes:")) {
                subpopulation_codes_ = trim_copy(line.substr(1));
                int inferred = infer_n_ancestries_from_codes(line);
                if (inferred > 0 && inferred != n_ancestries_) {
                    die(
                        "MSP subpopulation code count implies %d ancestries, but n_ancestries is %d",
                        inferred,
                        n_ancestries_
                    );
                }
                continue;
            }

            if (starts_with(line, "#chm") || starts_with(line, "chm")) {
                parse_header_line(line);
                header_seen_ = true;
                return;
            }

            if (line[0] == '#') continue;

            die("encountered MSP data before #chm header at line %llu", static_cast<unsigned long long>(line_no_));
        }

        die("RFMix MSP header not found in %s", path_);
    }

    void parse_header_line(const std::string& line) {
        std::vector<std::string> fields = split_msp_fields(line);
        size_t hap_start = 0;

        if (fields.size() >= 6 && fields[5] == "n snps") {
            hap_start = 6;
        } else if (fields.size() >= 7 && fields[5] == "n" && fields[6] == "snps") {
            hap_start = 7;
        } else {
            die("MSP #chm header must contain columns: #chm spos epos sgpos egpos n snps");
        }

        size_t expected_haps = static_cast<size_t>(2 * n_samples_);
        if (fields.size() != hap_start + expected_haps) {
            die(
                "MSP header has %llu haplotype columns; expected %llu from genotype VCF samples",
                static_cast<unsigned long long>(fields.size() - hap_start),
                static_cast<unsigned long long>(expected_haps)
            );
        }

        for (int i = 0; i < n_samples_; ++i) {
            std::string sample = ghdr_->samples[i];
            const std::string& h0 = fields[hap_start + static_cast<size_t>(2 * i)];
            const std::string& h1 = fields[hap_start + static_cast<size_t>(2 * i + 1)];

            if (!hap_label_matches(h0, sample, 0) ||
                !hap_label_matches(h1, sample, 1)) {
                die(
                    "MSP haplotype columns must match VCF sample order as sample.0/sample.1; at sample index %d expected %s.0 and %s.1, got %s and %s",
                    i,
                    sample.c_str(),
                    sample.c_str(),
                    h0.c_str(),
                    h1.c_str()
                );
            }
        }
    }

    const char* path_;
    bcf_hdr_t* ghdr_;
    int n_samples_;
    int n_ancestries_;
    int n_words_;
    htsFile* fp_ = nullptr;
    kstring_t line_ = {0, 0, nullptr};
    uint64_t line_no_ = 0;
    uint64_t records_read_ = 0;
    bool header_seen_ = false;
    std::string subpopulation_codes_;
    std::unordered_set<std::string> warned_missing_contigs_;
};

static bool block_covers(const OpenAncestryBlock& block, int rid, int64_t pos) {
    return block.active &&
           block.geno_rid == rid &&
           pos >= block.start_pos &&
           pos <= block.end_pos;
}

static void print_usage(const char* prog) {
    std::fprintf(
        stderr,
        "Usage:\n"
        "  %s genotype.phased.vcf.gz rfmix.msp.tsv[.gz] n_ancestries rare_threshold|auto out_prefix\n\n"
        "Example:\n"
        "  %s chr22.phased.vcf.gz chr22.rfmix.msp.tsv 5 auto chr22\n",
        prog,
        prog
    );
}

static int parse_rare_threshold_arg(const char* text) {
    if (std::strcmp(text, "auto") == 0 || std::strcmp(text, "default") == 0) {
        return -1;
    }

    char* end = nullptr;
    long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0') {
        die("rare_threshold must be a non-negative integer or auto: %s", text);
    }
    if (value < 0 || value > std::numeric_limits<int>::max()) {
        die("rare_threshold out of range: %s", text);
    }
    return static_cast<int>(value);
}

static int default_rare_threshold_from_samples(int n_samples) {
    return (n_samples + 31) / 32;
}

int main(int argc, char** argv) {
    if (argc < 6) {
        print_usage(argv[0]);
        return 1;
    }

    const char* geno_vcf = argv[1];
    const char* msp_path = argv[2];
    int n_ancestries = std::atoi(argv[3]);
    int rare_threshold = parse_rare_threshold_arg(argv[4]);
    const char* out_prefix = argv[5];

    if (n_ancestries <= 0 || n_ancestries > 32) {
        die("n_ancestries must be in [1, 32]");
    }

    htsFile* gfp = bcf_open(geno_vcf, "r");
    if (!gfp) {
        die("cannot open genotype VCF/BCF: %s", geno_vcf);
    }

    bcf_hdr_t* ghdr = bcf_hdr_read(gfp);
    if (!ghdr) {
        die("cannot read genotype header");
    }

    ProgressReporter progress("genotype VCF", geno_vcf, "split variants");

    int n_samples = bcf_hdr_nsamples(ghdr);
    if (n_samples <= 0) {
        die("genotype VCF has no samples");
    }
    if (rare_threshold < 0) {
        rare_threshold = default_rare_threshold_from_samples(n_samples);
        std::fprintf(
            stderr,
            "Using rare_threshold=%d = ceil(%d / 32) from genotype VCF sample count.\n",
            rare_threshold,
            n_samples
        );
    }

    uint64_t n_haps = static_cast<uint64_t>(n_samples) * 2ULL;
    if (n_haps > static_cast<uint64_t>(kMaxPackedHapId) + 1ULL) {
        die("hap_id exceeds 27-bit packed limit");
    }

    int n_words = static_cast<int>((n_haps + 63ULL) / 64ULL);

    MspReader msp_reader(msp_path, ghdr, n_samples, n_ancestries, n_words);

    std::string common_bin = std::string(out_prefix) + ".common.geno.bin";
    std::string common_mks = std::string(out_prefix) + ".common.variant.mks";
    std::string common_idx = std::string(out_prefix) + ".common.variant.idx";
    std::string rare_bin = std::string(out_prefix) + ".rare.carrier.bin";
    std::string rare_mks = std::string(out_prefix) + ".rare.variant.mks";
    std::string rare_idx = std::string(out_prefix) + ".rare.variant.idx";
    std::string anc_bin = std::string(out_prefix) + ".ancblock.bin";
    std::string anc_mks = std::string(out_prefix) + ".ancblock.mks";
    std::string anc_idx = std::string(out_prefix) + ".ancblock.idx";
    std::string samples_path = std::string(out_prefix) + ".samples";
    std::string meta_path = std::string(out_prefix) + ".meta";

    FILE* common_fp = open_output_or_die(common_bin, "wb");
    FILE* common_mks_fp = open_output_or_die(common_mks, "wb");
    FILE* common_idx_fp = open_output_or_die(common_idx, "wb");
    FILE* rare_fp = open_output_or_die(rare_bin, "wb");
    FILE* rare_mks_fp = open_output_or_die(rare_mks, "wb");
    FILE* rare_idx_fp = open_output_or_die(rare_idx, "wb");
    FILE* anc_fp = open_output_or_die(anc_bin, "wb");
    FILE* anc_mks_fp = open_output_or_die(anc_mks, "wb");
    FILE* anc_idx_fp = open_output_or_die(anc_idx, "wb");

    const char common_mks_magic[8] = {'T', 'R', 'C', 'M', 'M', 'K', 'S', '1'};
    const char rare_mks_magic[8] = {'T', 'R', 'R', 'A', 'M', 'K', 'S', '1'};
    const char anc_mks_magic[8] = {'T', 'R', 'A', 'N', 'M', 'K', 'S', '1'};
    const char common_idx_magic[8] = {'T', 'R', 'C', 'M', 'I', 'D', 'X', '2'};
    const char rare_idx_magic[8] = {'T', 'R', 'R', 'A', 'I', 'D', 'X', '2'};
    const char anc_idx_magic[8] = {'T', 'R', 'A', 'N', 'I', 'D', 'X', '2'};
    write_magic_header(common_mks_fp, common_mks_magic);
    write_magic_header(rare_mks_fp, rare_mks_magic);
    write_magic_header(anc_mks_fp, anc_mks_magic);
    write_magic_header(common_idx_fp, common_idx_magic);
    write_magic_header(rare_idx_fp, rare_idx_magic);
    write_magic_header(anc_idx_fp, anc_idx_magic);

    write_sidecars(
        ghdr,
        samples_path,
        meta_path,
        geno_vcf,
        msp_path,
        msp_reader.subpopulation_codes(),
        n_samples,
        n_haps,
        n_words,
        n_ancestries,
        rare_threshold
    );

    bcf1_t* grec = bcf_init();
    int32_t* gt_arr = nullptr;
    int ngt_arr = 0;

    OpenAncestryBlock block;
    std::vector<std::vector<CarrierHap>> carriers_by_alt;

    uint32_t n_blocks_written = 0;
    uint64_t common_index = 0;
    uint64_t rare_index = 0;
    uint32_t global_variant_index = 0;
    uint64_t msp_rows_consumed = 0;

    MspRecord msp_record;
    bool has_msp_record = msp_reader.read_next(msp_record);
    int last_msp_geno_rid = -1;
    int64_t last_msp_end = 0;

    while (read_next_record(gfp, ghdr, grec) == 0) {
        progress.record_scanned();

        const char* g_chr = bcf_hdr_id2name(ghdr, grec->rid);
        int64_t g_pos = static_cast<int64_t>(grec->pos) + 1;

        while (has_msp_record) {
            bool current_block_covers = block_covers(block, grec->rid, g_pos);
            if (current_block_covers) {
                break;
            }

            if (msp_record.geno_rid > grec->rid ||
                (msp_record.geno_rid == grec->rid && msp_record.start_pos > g_pos)) {
                break;
            }

            MspRecord interval_record = std::move(msp_record);
            has_msp_record = msp_reader.read_next(msp_record);

            if (last_msp_geno_rid > interval_record.geno_rid) {
                die("MSP records are not sorted by genotype header contig order");
            }

            int64_t interval_start = interval_record.start_pos;
            if (last_msp_geno_rid == interval_record.geno_rid) {
                if (interval_record.end_pos <= last_msp_end) {
                    die(
                        "MSP records do not advance at %s:%lld-%lld",
                        interval_record.chr.c_str(),
                        static_cast<long long>(interval_record.start_pos),
                        static_cast<long long>(interval_record.end_pos)
                    );
                }

                if (interval_start <= last_msp_end) {
                    interval_start = last_msp_end + 1;
                }
            }

            update_open_block_from_msp(
                anc_fp,
                anc_bin.c_str(),
                anc_mks_fp,
                anc_mks.c_str(),
                anc_idx_fp,
                block,
                n_blocks_written,
                n_ancestries,
                n_words,
                interval_record.geno_rid,
                interval_record.chr.c_str(),
                interval_start,
                interval_record.end_pos,
                interval_record.state
            );

            last_msp_geno_rid = interval_record.geno_rid;
            last_msp_end = interval_record.end_pos;
            ++msp_rows_consumed;
        }

        if (!block_covers(block, grec->rid, g_pos)) {
            progress.maybe_report(global_variant_index, common_index, rare_index, msp_rows_consumed);
            continue;
        }

        if (grec->n_allele < 2) {
            progress.maybe_report(global_variant_index, common_index, rare_index, msp_rows_consumed);
            continue;
        }

        process_genotypes(
            ghdr,
            grec,
            block,
            n_samples,
            &gt_arr,
            &ngt_arr,
            carriers_by_alt
        );

        if (carriers_by_alt.empty()) {
            progress.maybe_report(global_variant_index, common_index, rare_index, msp_rows_consumed);
            continue;
        }

        const char* raw_id = grec->d.id;
        const char* ref = grec->d.allele[0];

        for (int alt_idx = 1; alt_idx < grec->n_allele; ++alt_idx) {
            if (global_variant_index == std::numeric_limits<uint32_t>::max()) {
                die("global split variant index exceeds uint32_t limit");
            }

            const char* alt = grec->d.allele[alt_idx];
            const auto& carriers = carriers_by_alt[alt_idx];
            uint32_t mac = static_cast<uint32_t>(carriers.size());
            std::string split_id = make_split_id(raw_id, g_chr, g_pos, ref, alt);

            if (mac <= static_cast<uint32_t>(rare_threshold)) {
                uint64_t carrier_offset = tell_or_die(rare_fp, rare_bin.c_str());

                for (const CarrierHap& carrier : carriers) {
                    RareCarrierPacked packed{
                        global_variant_index,
                        pack_anc_hap(carrier.ancestry, carrier.hap_id)
                    };

                    size_t wrote = std::fwrite(
                        &packed,
                        sizeof(RareCarrierPacked),
                        1,
                        rare_fp
                    );

                    if (wrote != 1) {
                        die("failed writing rare carrier for split variant %u", global_variant_index);
                    }
                }

                uint64_t rare_mks_offset = tell_or_die(rare_mks_fp, rare_mks.c_str());
                write_rare_mks_record(
                    rare_mks_fp,
                    rare_index,
                    global_variant_index,
                    g_chr,
                    g_pos,
                    split_id,
                    ref,
                    alt,
                    static_cast<uint32_t>(alt_idx),
                    carrier_offset,
                    mac,
                    mac
                );
                write_rare_offset_idx_record(
                    rare_idx_fp,
                    rare_index,
                    global_variant_index,
                    rare_mks_offset,
                    carrier_offset,
                    mac
                );

                ++rare_index;
            } else {
                std::vector<uint64_t> bits(static_cast<size_t>(n_words), 0);

                for (const CarrierHap& carrier : carriers) {
                    set_bit(bits, carrier.hap_id);
                }

                uint64_t geno_offset = tell_or_die(common_fp, common_bin.c_str());
                size_t wrote = std::fwrite(
                    bits.data(),
                    sizeof(uint64_t),
                    static_cast<size_t>(n_words),
                    common_fp
                );

                if (wrote != static_cast<size_t>(n_words)) {
                    die("failed writing common genotype for split variant %u", global_variant_index);
                }

                uint64_t common_mks_offset = tell_or_die(common_mks_fp, common_mks.c_str());
                write_common_mks_record(
                    common_mks_fp,
                    common_index,
                    global_variant_index,
                    g_chr,
                    g_pos,
                    split_id,
                    ref,
                    alt,
                    static_cast<uint32_t>(alt_idx),
                    block.block_id,
                    geno_offset,
                    mac
                );
                write_common_offset_idx_record(
                    common_idx_fp,
                    common_index,
                    global_variant_index,
                    common_mks_offset,
                    geno_offset
                );

                ++common_index;
            }

            ++global_variant_index;
        }

        progress.maybe_report(global_variant_index, common_index, rare_index, msp_rows_consumed);
    }

    close_open_block(
        anc_fp,
        anc_bin.c_str(),
        anc_mks_fp,
        anc_mks.c_str(),
        anc_idx_fp,
        block,
        n_blocks_written,
        n_ancestries,
        n_words
    );

    if (gt_arr) std::free(gt_arr);

    bcf_destroy(grec);
    bcf_hdr_destroy(ghdr);
    bcf_close(gfp);

    std::fclose(common_fp);
    std::fclose(common_mks_fp);
    std::fclose(common_idx_fp);
    std::fclose(rare_fp);
    std::fclose(rare_mks_fp);
    std::fclose(rare_idx_fp);
    std::fclose(anc_fp);
    std::fclose(anc_mks_fp);
    std::fclose(anc_idx_fp);

    progress.finish(global_variant_index, common_index, rare_index, msp_rows_consumed);

    std::fprintf(stderr, "Finished.\n");
    std::fprintf(stderr, "Global split variants: %u\n", global_variant_index);
    std::fprintf(stderr, "Common variants:       %llu\n", static_cast<unsigned long long>(common_index));
    std::fprintf(stderr, "Rare variants:         %llu\n", static_cast<unsigned long long>(rare_index));
    std::fprintf(stderr, "Ancestry blocks:       %u\n", n_blocks_written);
    std::fprintf(stderr, "MSP rows consumed:     %llu\n", static_cast<unsigned long long>(msp_rows_consumed));

    return 0;
}
