// tractor_dosage_vcf_to_hybrid.cpp
// designed by Kai, implemented by codex
//
// Convert a split-biallelic SAIGE-TRACTOR-style dosage VCF/BCF with scalar
// hardcall DS1..DSk and ANC1..ANCk FORMAT fields into ancestry-aware packed
// files for the tractor_hybrid backend.

#include <htslib/hts.h>
#include <htslib/vcf.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
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
    int rid = -1;
    std::string chr;
    int64_t start_pos = -1;
    int64_t end_pos = -1;
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
        : input_label_(input_label),
          converted_label_(converted_label),
          total_(get_index_record_count(input_path)),
          last_report_(std::chrono::steady_clock::now()),
          stderr_is_tty_(isatty(fileno(stderr)) != 0) {
        if (total_.available) {
            std::fprintf(
                stderr,
                "Progress: %s index reports %llu records.\n",
                input_label_,
                static_cast<unsigned long long>(total_.total)
            );
        } else {
            std::fprintf(
                stderr,
                "Progress: %s index record count unavailable; reporting scanned records only.\n",
                input_label_
            );
        }
    }

    void record_scanned() {
        ++records_scanned_;
    }

    void maybe_report(uint64_t converted, uint64_t common, uint64_t rare, bool force = false) {
        auto now = std::chrono::steady_clock::now();
        bool count_due = records_scanned_ >= next_record_report_;
        bool time_due = now - last_report_ >= std::chrono::seconds(kProgressSecondsInterval);

        if (!force && !count_due && !time_due) {
            return;
        }

        print(converted, common, rare);
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

    void finish(uint64_t converted, uint64_t common, uint64_t rare) {
        maybe_report(converted, common, rare, true);
        if (stderr_is_tty_ && printed_) {
            std::fputc('\n', stderr);
        }
    }

private:
    void print(uint64_t converted, uint64_t common, uint64_t rare) const {
        const char* prefix = stderr_is_tty_ ? "\r" : "";
        const char* suffix = stderr_is_tty_ ? "" : "\n";

        if (total_.available && total_.total > 0) {
            uint64_t capped_records = std::min(records_scanned_, total_.total);
            double pct = 100.0 * static_cast<double>(capped_records) /
                         static_cast<double>(total_.total);
            std::fprintf(
                stderr,
                "%sProgress: records %llu/%llu (%.1f%%), converted %llu %s, common %llu, rare %llu%s",
                prefix,
                static_cast<unsigned long long>(records_scanned_),
                static_cast<unsigned long long>(total_.total),
                pct,
                static_cast<unsigned long long>(converted),
                converted_label_,
                static_cast<unsigned long long>(common),
                static_cast<unsigned long long>(rare),
                suffix
            );
        } else if (total_.available) {
            std::fprintf(
                stderr,
                "%sProgress: records %llu/0, converted %llu %s, common %llu, rare %llu%s",
                prefix,
                static_cast<unsigned long long>(records_scanned_),
                static_cast<unsigned long long>(converted),
                converted_label_,
                static_cast<unsigned long long>(common),
                static_cast<unsigned long long>(rare),
                suffix
            );
        } else {
            std::fprintf(
                stderr,
                "%sProgress: records %llu, converted %llu %s, common %llu, rare %llu%s",
                prefix,
                static_cast<unsigned long long>(records_scanned_),
                static_cast<unsigned long long>(converted),
                converted_label_,
                static_cast<unsigned long long>(common),
                static_cast<unsigned long long>(rare),
                suffix
            );
        }

        std::fflush(stderr);
    }

    const char* input_label_;
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
    bcf_hdr_t* hdr,
    const std::string& samples_path,
    const std::string& meta_path,
    const char* dosage_vcf,
    int n_samples,
    uint64_t n_haps,
    int n_words,
    int n_ancestries,
    int rare_threshold
) {
    FILE* samples_fp = open_output_or_die(samples_path, "w");
    for (int i = 0; i < n_samples; ++i) {
        std::fprintf(samples_fp, "%s\n", hdr->samples[i]);
    }
    std::fclose(samples_fp);

    FILE* meta_fp = open_output_or_die(meta_path, "w");
    std::fprintf(meta_fp, "format_version\t1\n");
    std::fprintf(meta_fp, "n_samples\t%d\n", n_samples);
    std::fprintf(meta_fp, "n_haps\t%llu\n", static_cast<unsigned long long>(n_haps));
    std::fprintf(meta_fp, "n_words\t%d\n", n_words);
    std::fprintf(meta_fp, "n_ancestries\t%d\n", n_ancestries);
    std::fprintf(meta_fp, "rare_threshold\t%d\n", rare_threshold);
    std::fprintf(meta_fp, "source_tractor_dosage_vcf\t%s\n", dosage_vcf);
    std::fprintf(
        meta_fp,
        "dosage_source_note\tcanonical haplotypes reconstructed from hardcall DS#/ANC# counts\n"
    );
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

static uint32_t update_open_block_from_state(
    FILE* anc_fp,
    const char* anc_bin_path,
    FILE* anc_mks_fp,
    const char* anc_mks_path,
    FILE* anc_idx_fp,
    OpenAncestryBlock& block,
    uint32_t& n_blocks_written,
    int n_ancestries,
    int n_words,
    int rid,
    const char* chr,
    int64_t pos,
    AncestryState& new_state
) {
    if (!block.active) {
        block.active = true;
        block.block_id = n_blocks_written;
        block.rid = rid;
        block.chr = chr;
        block.start_pos = pos;
        block.end_pos = pos;
        block.state = std::move(new_state);
        return block.block_id;
    }

    bool same_chromosome = block.rid == rid;
    bool unchanged = same_chromosome && same_masks(block.state.masks, new_state.masks);

    if (unchanged) {
        block.end_pos = std::max(block.end_pos, pos);
        return block.block_id;
    }

    if (same_chromosome && pos <= block.end_pos) {
        die(
            "duplicate position %s:%lld has inconsistent ANC fields; hybrid rare markers need a unique ancestry block per position",
            chr,
            static_cast<long long>(pos)
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
    block.rid = rid;
    block.chr = chr;
    block.start_pos = pos;
    block.end_pos = pos;
    block.state = std::move(new_state);
    return block.block_id;
}

static std::string make_variant_id(
    const char* raw_id,
    const char* chr,
    int64_t pos,
    const char* ref,
    const char* alt
) {
    if (raw_id == nullptr || std::strcmp(raw_id, ".") == 0) {
        return std::string(chr) + ":" + std::to_string(pos) + ":" + ref + ":" + alt;
    }

    return raw_id;
}

static int read_next_record(htsFile* fp, bcf_hdr_t* hdr, bcf1_t* rec) {
    int ret = bcf_read(fp, hdr, rec);
    if (ret == 0) bcf_unpack(rec, BCF_UN_STR);
    return ret;
}

static int hardcall_from_float(
    float value,
    const char* field_name,
    const char* chr,
    int64_t pos,
    int sample_index
) {
    if (bcf_float_is_missing(value) || bcf_float_is_vector_end(value)) {
        die(
            "missing FORMAT/%s at %s:%lld sample index %d",
            field_name,
            chr,
            static_cast<long long>(pos),
            sample_index
        );
    }

    if (!std::isfinite(static_cast<double>(value))) {
        die(
            "non-finite FORMAT/%s at %s:%lld sample index %d",
            field_name,
            chr,
            static_cast<long long>(pos),
            sample_index
        );
    }

    double rounded = std::round(static_cast<double>(value));
    if (std::fabs(static_cast<double>(value) - rounded) > 1e-6) {
        die(
            "FORMAT/%s is fractional at %s:%lld sample index %d: %.8g. Current hybrid conversion requires hardcall integer DS/ANC values",
            field_name,
            chr,
            static_cast<long long>(pos),
            sample_index,
            static_cast<double>(value)
        );
    }

    if (rounded < 0.0 || rounded > 2.0) {
        die(
            "FORMAT/%s out of diploid range at %s:%lld sample index %d",
            field_name,
            chr,
            static_cast<long long>(pos),
            sample_index
        );
    }

    return static_cast<int>(rounded);
}

static void read_scalar_hardcall_field(
    bcf_hdr_t* hdr,
    bcf1_t* rec,
    const std::string& field_name,
    int n_samples,
    float** values,
    int* n_values,
    std::vector<int>& out
) {
    const char* chr = bcf_hdr_id2name(hdr, rec->rid);
    int64_t pos = static_cast<int64_t>(rec->pos) + 1;

    int n = bcf_get_format_float(hdr, rec, field_name.c_str(), values, n_values);
    if (n <= 0) {
        die("missing FORMAT/%s at %s:%lld", field_name.c_str(), chr, static_cast<long long>(pos));
    }

    if (n != n_samples) {
        die("FORMAT/%s must be scalar Number=1 at %s:%lld", field_name.c_str(), chr, static_cast<long long>(pos));
    }

    out.resize(static_cast<size_t>(n_samples));
    for (int i = 0; i < n_samples; ++i) {
        out[i] = hardcall_from_float((*values)[i], field_name.c_str(), chr, pos, i);
    }
}

static void build_state_and_carriers(
    const std::vector<std::vector<int>>& anc_counts,
    const std::vector<std::vector<int>>& ds_counts,
    const std::vector<std::string>& ds_names,
    int n_samples,
    int n_ancestries,
    int n_words,
    const char* chr,
    int64_t pos,
    AncestryState& state,
    std::vector<CarrierHap>& carriers
) {
    state.masks.assign(
        static_cast<size_t>(n_ancestries),
        std::vector<uint64_t>(static_cast<size_t>(n_words), 0)
    );
    state.hap_ancestry.assign(static_cast<size_t>(2 * n_samples), -1);
    carriers.clear();

    for (int i = 0; i < n_samples; ++i) {
        int anc_sum = 0;
        int ds_sum = 0;
        uint32_t next_hap = static_cast<uint32_t>(2 * i);

        for (int k = 0; k < n_ancestries; ++k) {
            int anc = anc_counts[k][i];
            int ds = ds_counts[k][i];

            anc_sum += anc;
            if (anc_sum > 2) {
                die(
                    "ANC fields sum above 2 at %s:%lld sample index %d",
                    chr,
                    static_cast<long long>(pos),
                    i
                );
            }

            if (ds > anc) {
                die(
                    "FORMAT/%s exceeds matching ANC count at %s:%lld sample index %d",
                    ds_names[k].c_str(),
                    chr,
                    static_cast<long long>(pos),
                    i
                );
            }

            ds_sum += ds;
            if (ds_sum > 2) {
                die(
                    "DS fields sum above 2 at %s:%lld sample index %d",
                    chr,
                    static_cast<long long>(pos),
                    i
                );
            }

            for (int copy = 0; copy < anc; ++copy) {
                uint32_t hap_id = next_hap++;
                set_bit(state.masks[k], hap_id);
                state.hap_ancestry[hap_id] = static_cast<int8_t>(k);

                if (copy < ds) {
                    carriers.push_back(CarrierHap{hap_id, static_cast<uint8_t>(k)});
                }
            }
        }

        if (anc_sum != 2) {
            die(
                "ANC fields must sum to 2 at %s:%lld sample index %d",
                chr,
                static_cast<long long>(pos),
                i
            );
        }
    }
}

static void print_usage(const char* prog) {
    std::fprintf(
        stderr,
        "Usage:\n"
        "  %s tractor_dosage.vcf[.gz|.bcf] n_ancestries rare_threshold|auto out_prefix\n\n"
        "Input FORMAT must contain scalar hardcall DS1..DSk and ANC1..ANCk values.\n",
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
    if (argc < 5) {
        print_usage(argv[0]);
        return 1;
    }

    const char* dosage_vcf = argv[1];
    int n_ancestries = std::atoi(argv[2]);
    int rare_threshold = parse_rare_threshold_arg(argv[3]);
    const char* out_prefix = argv[4];

    if (n_ancestries <= 0 || n_ancestries > 32) {
        die("n_ancestries must be in [1, 32]");
    }

    htsFile* fp = bcf_open(dosage_vcf, "r");
    if (!fp) {
        die("cannot open input VCF/BCF");
    }

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        die("cannot read input header");
    }

    ProgressReporter progress("dosage VCF", dosage_vcf, "variants");

    int n_samples = bcf_hdr_nsamples(hdr);
    if (n_samples <= 0) {
        die("input has no samples");
    }
    if (rare_threshold < 0) {
        rare_threshold = default_rare_threshold_from_samples(n_samples);
        std::fprintf(
            stderr,
            "Using rare_threshold=%d = ceil(%d / 32) from dosage VCF sample count.\n",
            rare_threshold,
            n_samples
        );
    }

    uint64_t n_haps = static_cast<uint64_t>(n_samples) * 2ULL;
    if (n_haps > static_cast<uint64_t>(kMaxPackedHapId) + 1ULL) {
        die("hap_id exceeds 27-bit packed limit");
    }

    int n_words = static_cast<int>((n_haps + 63ULL) / 64ULL);

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
        hdr,
        samples_path,
        meta_path,
        dosage_vcf,
        n_samples,
        n_haps,
        n_words,
        n_ancestries,
        rare_threshold
    );

    std::vector<std::string> ds_names;
    std::vector<std::string> anc_names;
    ds_names.reserve(static_cast<size_t>(n_ancestries));
    anc_names.reserve(static_cast<size_t>(n_ancestries));
    for (int k = 0; k < n_ancestries; ++k) {
        ds_names.push_back("DS" + std::to_string(k + 1));
        anc_names.push_back("ANC" + std::to_string(k + 1));
    }

    bcf1_t* rec = bcf_init();

    std::vector<float*> ds_arrays(static_cast<size_t>(n_ancestries), nullptr);
    std::vector<int> nds_arrays(static_cast<size_t>(n_ancestries), 0);
    std::vector<float*> anc_arrays(static_cast<size_t>(n_ancestries), nullptr);
    std::vector<int> nanc_arrays(static_cast<size_t>(n_ancestries), 0);
    std::vector<std::vector<int>> ds_counts(static_cast<size_t>(n_ancestries));
    std::vector<std::vector<int>> anc_counts(static_cast<size_t>(n_ancestries));

    OpenAncestryBlock block;
    AncestryState state;
    std::vector<CarrierHap> carriers;

    uint32_t n_blocks_written = 0;
    uint64_t common_index = 0;
    uint64_t rare_index = 0;
    uint32_t global_variant_index = 0;

    while (read_next_record(fp, hdr, rec) == 0) {
        progress.record_scanned();

        const char* chr = bcf_hdr_id2name(hdr, rec->rid);
        int64_t pos = static_cast<int64_t>(rec->pos) + 1;

        if (rec->n_allele != 2) {
            die(
                "input must be split biallelic because DS#/ANC# are scalar; found %d alleles at %s:%lld",
                rec->n_allele,
                chr,
                static_cast<long long>(pos)
            );
        }

        for (int k = 0; k < n_ancestries; ++k) {
            read_scalar_hardcall_field(
                hdr,
                rec,
                anc_names[k],
                n_samples,
                &anc_arrays[k],
                &nanc_arrays[k],
                anc_counts[k]
            );
            read_scalar_hardcall_field(
                hdr,
                rec,
                ds_names[k],
                n_samples,
                &ds_arrays[k],
                &nds_arrays[k],
                ds_counts[k]
            );
        }

        build_state_and_carriers(
            anc_counts,
            ds_counts,
            ds_names,
            n_samples,
            n_ancestries,
            n_words,
            chr,
            pos,
            state,
            carriers
        );

        uint32_t block_id = update_open_block_from_state(
            anc_fp,
            anc_bin.c_str(),
            anc_mks_fp,
            anc_mks.c_str(),
            anc_idx_fp,
            block,
            n_blocks_written,
            n_ancestries,
            n_words,
            rec->rid,
            chr,
            pos,
            state
        );

        if (global_variant_index == std::numeric_limits<uint32_t>::max()) {
            die("global split variant index exceeds uint32_t limit");
        }

        const char* ref = rec->d.allele[0];
        const char* alt = rec->d.allele[1];
        std::string variant_id = make_variant_id(rec->d.id, chr, pos, ref, alt);
        uint32_t mac = static_cast<uint32_t>(carriers.size());

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
                    die("failed writing rare carriers for split variant %u", global_variant_index);
                }
            }

            uint64_t rare_mks_offset = tell_or_die(rare_mks_fp, rare_mks.c_str());
            write_rare_mks_record(
                rare_mks_fp,
                rare_index,
                global_variant_index,
                chr,
                pos,
                variant_id,
                ref,
                alt,
                1,
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
                die("failed writing common genotype bitset for split variant %u", global_variant_index);
            }

            uint64_t common_mks_offset = tell_or_die(common_mks_fp, common_mks.c_str());
            write_common_mks_record(
                common_mks_fp,
                common_index,
                global_variant_index,
                chr,
                pos,
                variant_id,
                ref,
                alt,
                1,
                block_id,
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
        progress.maybe_report(global_variant_index, common_index, rare_index);
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

    for (float* ptr : ds_arrays) {
        if (ptr) std::free(ptr);
    }
    for (float* ptr : anc_arrays) {
        if (ptr) std::free(ptr);
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    std::fclose(common_fp);
    std::fclose(common_mks_fp);
    std::fclose(common_idx_fp);
    std::fclose(rare_fp);
    std::fclose(rare_mks_fp);
    std::fclose(rare_idx_fp);
    std::fclose(anc_fp);
    std::fclose(anc_mks_fp);
    std::fclose(anc_idx_fp);

    progress.finish(global_variant_index, common_index, rare_index);

    std::fprintf(stderr, "Finished.\n");
    std::fprintf(stderr, "Global variants:       %u\n", global_variant_index);
    std::fprintf(stderr, "Common variants:       %llu\n", static_cast<unsigned long long>(common_index));
    std::fprintf(stderr, "Rare variants:         %llu\n", static_cast<unsigned long long>(rare_index));
    std::fprintf(stderr, "Ancestry blocks:       %u\n", n_blocks_written);

    return 0;
}
