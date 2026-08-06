// flare_subset_to_tractor_hybrid.cpp
// designed by Kai, implemented by codex
//
// Convert phased genotype VCF/BCF plus FLARE local ancestry VCF/BCF into
// ancestry-aware packed files for a SAIGE-TRACTOR genotype backend.

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>

#include <algorithm>
#include <chrono>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <unistd.h>

struct RareCarrierPacked {
    uint32_t pos_index; // split biallelic variant ordinal
    uint32_t anc_hap;   // high 5 bits ancestry, low 27 bits hap_id
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

struct LaiRecord {
    bool valid = false;
    int geno_rid = -1;
    std::string chr;
    int64_t pos = 0;
    AncestryState state;
};

struct Region {
    bool active = false;
    std::string label;
    std::string chr;
    int64_t start = 0;
    int64_t end = 0;
    int geno_rid = -1;
};

struct KeepSamples {
    bool active = false;
    std::string path;
    std::unordered_set<std::string> ids;
};

struct SampleSelection {
    std::vector<int> genotype_raw_indices;
    std::vector<int> flare_raw_indices;
    std::vector<std::string> sample_ids;
    int genotype_raw_sample_count = 0;
    int flare_raw_sample_count = 0;
    bool keep_active = false;
    std::string keep_path;
};

struct ExtractPosition {
    std::string chr;
    int64_t pos = 0;
    std::string ref;
};

struct ExtractSites {
    bool active = false;
    std::string path;
    std::unordered_set<std::string> allele_keys;
    std::unordered_map<std::string, std::string> ref_by_position;
    std::unordered_set<std::string> contigs;
    std::vector<ExtractPosition> positions;
};

static constexpr uint32_t kMaxPackedHapId = (1u << 27) - 1u;
static constexpr uint64_t kProgressRecordInterval = 10000;
static constexpr int kProgressSecondsInterval = 2;
static constexpr int64_t kMaxVcfCoordinate = 2147483647LL;

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
    hts_idx_t* idx = bcf_index_load3(path, nullptr, HTS_IDX_SILENT_FAIL);
    if (idx) return idx;

    idx = hts_idx_load3(path, nullptr, HTS_FMT_CSI, HTS_IDX_SILENT_FAIL);
    if (idx) return idx;
    return hts_idx_load3(path, nullptr, HTS_FMT_TBI, HTS_IDX_SILENT_FAIL);
}

static int64_t parse_i64_string(const std::string& text, const char* label) {
    if (text.empty()) die("%s is empty", label);

    char* end = nullptr;
    long long value = std::strtoll(text.c_str(), &end, 10);
    if (end == text.c_str() || *end != '\0') {
        die("invalid integer for %s: %s", label, text.c_str());
    }
    if (value <= 0) {
        die("%s must be positive: %s", label, text.c_str());
    }
    return static_cast<int64_t>(value);
}

static Region parse_region_string(const std::string& region_text) {
    size_t colon = region_text.find(':');
    size_t dash = region_text.find('-', colon == std::string::npos ? 0 : colon + 1);
    if (colon == std::string::npos || dash == std::string::npos || dash <= colon + 1) {
        die("region must look like chr:start-end: %s", region_text.c_str());
    }

    Region region;
    region.active = true;
    region.chr = region_text.substr(0, colon);
    region.start = parse_i64_string(region_text.substr(colon + 1, dash - colon - 1), "region start");
    region.end = parse_i64_string(region_text.substr(dash + 1), "region end");
    if (region.chr.empty()) die("region chromosome is empty");
    if (region.start > region.end) {
        die("invalid region coordinates: %s", region_text.c_str());
    }
    region.label = region.chr + ":" + std::to_string(region.start) + "-" + std::to_string(region.end);
    return region;
}

class TextLineReader {
public:
    TextLineReader(const std::string& path, const char* label)
        : path_(path) {
        fp_ = hts_open(path.c_str(), "r");
        if (!fp_) die("cannot open %s: %s", label, path.c_str());
    }

    TextLineReader(const TextLineReader&) = delete;
    TextLineReader& operator=(const TextLineReader&) = delete;

    ~TextLineReader() {
        std::free(line_.s);
        if (fp_) hts_close(fp_);
    }

    bool next(std::string& out) {
        int ret = hts_getline(fp_, '\n', &line_);
        if (ret >= 0) {
            out.assign(line_.s, static_cast<size_t>(line_.l));
            return true;
        }
        if (ret == -1) return false;
        die("error reading text file: %s", path_.c_str());
    }

private:
    std::string path_;
    htsFile* fp_ = nullptr;
    kstring_t line_{0, 0, nullptr};
};

static std::string strip_trailing_cr(std::string value) {
    if (!value.empty() && value.back() == '\r') {
        value.pop_back();
    }
    return value;
}

static std::vector<std::string> split_commas(const std::string& value) {
    std::vector<std::string> parts;
    size_t begin = 0;
    while (begin <= value.size()) {
        size_t comma = value.find(',', begin);
        size_t end = comma == std::string::npos ? value.size() : comma;
        std::string part = value.substr(begin, end - begin);
        if (part.empty()) die("empty ALT allele in --extract list");
        parts.push_back(part);
        if (comma == std::string::npos) break;
        begin = comma + 1;
    }
    return parts;
}

static std::string position_key(const std::string& chr, int64_t pos) {
    return chr + "\t" + std::to_string(pos);
}

static std::string allele_key(
    const std::string& chr,
    int64_t pos,
    const std::string& ref,
    const std::string& alt
) {
    return position_key(chr, pos) + "\t" + ref + "\t" + alt;
}

static KeepSamples load_keep_samples(const std::string& path) {
    KeepSamples keep;
    if (path.empty()) return keep;

    keep.active = true;
    keep.path = path;

    TextLineReader reader(path, "--keep sample list");

    std::string line;
    uint64_t line_no = 0;
    while (reader.next(line)) {
        ++line_no;
        line = strip_trailing_cr(line);
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string sample_id;
        std::string extra;
        iss >> sample_id;
        if (sample_id.empty()) continue;
        if (iss >> extra) {
            die("--keep expects one sample ID per line at %s:%llu", path.c_str(),
                static_cast<unsigned long long>(line_no));
        }
        if (!keep.ids.insert(sample_id).second) {
            die("duplicate sample ID in --keep list at %s:%llu: %s", path.c_str(),
                static_cast<unsigned long long>(line_no), sample_id.c_str());
        }
    }

    if (keep.ids.empty()) {
        die("--keep sample list is empty: %s", path.c_str());
    }

    return keep;
}

static bool sample_headers_identical(bcf_hdr_t* ghdr, bcf_hdr_t* ahdr) {
    int n_genotype = bcf_hdr_nsamples(ghdr);
    int n_flare = bcf_hdr_nsamples(ahdr);
    if (n_genotype != n_flare) return false;

    for (int i = 0; i < n_genotype; ++i) {
        if (std::strcmp(ghdr->samples[i], ahdr->samples[i]) != 0) {
            return false;
        }
    }
    return true;
}

static std::unordered_map<std::string, int> build_sample_index(
    bcf_hdr_t* hdr,
    const char* label
) {
    int n_samples = bcf_hdr_nsamples(hdr);
    std::unordered_map<std::string, int> by_id;
    by_id.reserve(static_cast<size_t>(n_samples) * 2);

    for (int i = 0; i < n_samples; ++i) {
        std::string sample(hdr->samples[i]);
        if (!by_id.emplace(sample, i).second) {
            die("duplicate sample ID in %s VCF header: %s", label, sample.c_str());
        }
    }
    return by_id;
}

static SampleSelection build_sample_selection(
    bcf_hdr_t* ghdr,
    bcf_hdr_t* ahdr,
    const std::string& keep_path
) {
    SampleSelection selection;
    selection.genotype_raw_sample_count = bcf_hdr_nsamples(ghdr);
    selection.flare_raw_sample_count = bcf_hdr_nsamples(ahdr);
    if (selection.genotype_raw_sample_count <= 0) die("genotype VCF has no samples");
    if (selection.flare_raw_sample_count <= 0) die("FLARE VCF has no samples");

    KeepSamples keep = load_keep_samples(keep_path);
    selection.keep_active = keep.active;
    selection.keep_path = keep.path;
    std::unordered_set<std::string> keep_missing_from_genotype = keep.ids;

    std::unordered_set<std::string> genotype_seen;
    genotype_seen.reserve(static_cast<size_t>(selection.genotype_raw_sample_count) * 2);

    bool identical = sample_headers_identical(ghdr, ahdr);
    std::unordered_map<std::string, int> flare_index;
    if (!identical) {
        flare_index = build_sample_index(ahdr, "FLARE");
    }

    selection.genotype_raw_indices.reserve(static_cast<size_t>(selection.genotype_raw_sample_count));
    selection.flare_raw_indices.reserve(static_cast<size_t>(selection.genotype_raw_sample_count));
    selection.sample_ids.reserve(static_cast<size_t>(selection.genotype_raw_sample_count));

    for (int genotype_i = 0; genotype_i < selection.genotype_raw_sample_count; ++genotype_i) {
        std::string sample(ghdr->samples[genotype_i]);
        if (!genotype_seen.insert(sample).second) {
            die("duplicate sample ID in genotype VCF header: %s", sample.c_str());
        }

        if (keep.active) {
            auto keep_it = keep.ids.find(sample);
            if (keep_it == keep.ids.end()) continue;
            keep_missing_from_genotype.erase(sample);
        }

        int flare_i = genotype_i;
        if (!identical) {
            auto flare_it = flare_index.find(sample);
            if (flare_it == flare_index.end()) continue;
            flare_i = flare_it->second;
        }

        selection.genotype_raw_indices.push_back(genotype_i);
        selection.flare_raw_indices.push_back(flare_i);
        selection.sample_ids.push_back(sample);
    }

    if (keep.active && !keep_missing_from_genotype.empty()) {
        die("sample ID in --keep list is absent from genotype VCF: %s",
            keep_missing_from_genotype.begin()->c_str());
    }

    if (selection.sample_ids.empty()) {
        if (keep.active) {
            die("--keep retained zero samples after genotype/FLARE intersection: %s",
                keep.path.c_str());
        }
        die("genotype and FLARE sample intersection is empty");
    }

    if (keep.active) {
        std::fprintf(
            stderr,
            "Loaded --keep sample list: retaining %llu of %llu requested sample(s) after genotype/FLARE intersection.\n",
            static_cast<unsigned long long>(selection.sample_ids.size()),
            static_cast<unsigned long long>(keep.ids.size())
        );
    }

    if (!identical || static_cast<int>(selection.sample_ids.size()) != selection.genotype_raw_sample_count ||
        static_cast<int>(selection.sample_ids.size()) != selection.flare_raw_sample_count) {
        std::fprintf(
            stderr,
            "Using genotype/FLARE sample intersection: retaining %llu sample(s) in genotype VCF order (genotype %d, FLARE %d).\n",
            static_cast<unsigned long long>(selection.sample_ids.size()),
            selection.genotype_raw_sample_count,
            selection.flare_raw_sample_count
        );
    }

    return selection;
}

static ExtractSites load_extract_sites(const std::string& path) {
    ExtractSites sites;
    if (path.empty()) return sites;

    sites.active = true;
    sites.path = path;

    TextLineReader reader(path, "--extract site list");

    std::string line;
    uint64_t line_no = 0;
    while (reader.next(line)) {
        ++line_no;
        line = strip_trailing_cr(line);
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string chr;
        std::string pos_text;
        std::string ignored_id;
        std::string ref;
        std::string alts;
        // Consume the ID column only to preserve PVAR/VCF column alignment.
        // Site filtering is keyed exclusively by CHROM, POS, REF, and ALT.
        if (!(iss >> chr >> pos_text >> ignored_id >> ref >> alts)) {
            die("--extract expects PVAR/VCF columns CHROM POS ID REF ALT at %s:%llu",
                path.c_str(), static_cast<unsigned long long>(line_no));
        }
        if (chr == "CHROM" || chr == "#CHROM") continue;

        int64_t pos = parse_i64_string(pos_text, "--extract POS");
        if (ref.empty() || ref == ".") {
            die("--extract requires known REF at %s:%llu", path.c_str(),
                static_cast<unsigned long long>(line_no));
        }
        if (alts.empty() || alts == ".") {
            die("--extract requires known ALT at %s:%llu", path.c_str(),
                static_cast<unsigned long long>(line_no));
        }

        std::string pos_key = position_key(chr, pos);
        auto ref_it = sites.ref_by_position.find(pos_key);
        if (ref_it == sites.ref_by_position.end()) {
            sites.ref_by_position.emplace(pos_key, ref);
            sites.contigs.insert(chr);
            sites.positions.push_back(ExtractPosition{chr, pos, ref});
        } else if (ref_it->second != ref) {
            die("conflicting REF values in --extract list at %s:%llu for %s:%lld: %s vs %s",
                path.c_str(), static_cast<unsigned long long>(line_no), chr.c_str(),
                static_cast<long long>(pos), ref_it->second.c_str(), ref.c_str());
        }

        for (const std::string& alt : split_commas(alts)) {
            std::string key = allele_key(chr, pos, ref, alt);
            if (!sites.allele_keys.insert(key).second) {
                die("duplicate allele in --extract list at %s:%llu: %s:%lld %s>%s",
                    path.c_str(), static_cast<unsigned long long>(line_no), chr.c_str(),
                    static_cast<long long>(pos), ref.c_str(), alt.c_str());
            }
        }
    }

    if (sites.allele_keys.empty()) {
        die("--extract site list is empty: %s", path.c_str());
    }
    return sites;
}

static bool extract_record_position_matches(
    const ExtractSites& sites,
    const char* chr,
    int64_t pos,
    const char* ref
) {
    if (!sites.active) return true;

    std::string pos_key = position_key(chr, pos);
    auto it = sites.ref_by_position.find(pos_key);
    if (it == sites.ref_by_position.end()) return false;

    if (!ref || !*ref || std::strcmp(ref, ".") == 0) {
        die("genotype VCF has unknown REF at --extract position %s:%lld",
            chr, static_cast<long long>(pos));
    }
    if (it->second != ref) {
        die("REF mismatch for --extract site %s:%lld: list has %s, genotype VCF has %s",
            chr, static_cast<long long>(pos), it->second.c_str(), ref);
    }
    return true;
}

static bool extract_allele_matches(
    const ExtractSites& sites,
    const char* chr,
    int64_t pos,
    const char* ref,
    const char* alt
) {
    if (!sites.active) return true;
    return sites.allele_keys.count(allele_key(chr, pos, ref, alt)) > 0;
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
    ProgressReporter(
        const char* input_label,
        const char* input_path,
        const char* converted_label,
        const char* region_label = nullptr
    )
        : input_label_(input_label),
          converted_label_(converted_label),
          total_(region_label ? IndexRecordCount{} : get_index_record_count(input_path)),
          last_report_(std::chrono::steady_clock::now()),
          stderr_is_tty_(isatty(fileno(stderr)) != 0) {
        if (region_label) {
            std::fprintf(
                stderr,
                "Progress: %s region %s; reporting scanned region records only.\n",
                input_label_,
                region_label
            );
        } else if (total_.available) {
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
    write_exact(fp, magic, 8, "index header");
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
    const std::string& samples_path,
    const std::string& meta_path,
    const char* geno_vcf,
    const char* flare_vcf,
    const std::vector<std::string>& sample_ids,
    uint64_t n_haps,
    int n_words,
    int n_ancestries,
    int rare_threshold,
    const char* selected_region,
    const char* keep_path,
    const char* extract_path
) {
    FILE* samples_fp = open_output_or_die(samples_path, "w");
    for (const std::string& sample_id : sample_ids) {
        std::fprintf(samples_fp, "%s\n", sample_id.c_str());
    }
    std::fclose(samples_fp);

    FILE* meta_fp = open_output_or_die(meta_path, "w");
    std::fprintf(meta_fp, "format_version\t1\n");
    std::fprintf(meta_fp, "n_samples\t%llu\n",
        static_cast<unsigned long long>(sample_ids.size()));
    std::fprintf(meta_fp, "n_haps\t%llu\n", static_cast<unsigned long long>(n_haps));
    std::fprintf(meta_fp, "n_words\t%d\n", n_words);
    std::fprintf(meta_fp, "n_ancestries\t%d\n", n_ancestries);
    std::fprintf(meta_fp, "rare_threshold\t%d\n", rare_threshold);
    std::fprintf(meta_fp, "source_genotype\t%s\n", geno_vcf);
    std::fprintf(meta_fp, "source_flare\t%s\n", flare_vcf);
    if (selected_region) {
        std::fprintf(meta_fp, "selected_region\t%s\n", selected_region);
    }
    if (keep_path) {
        std::fprintf(meta_fp, "keep_samples\t%s\n", keep_path);
    }
    if (extract_path) {
        std::fprintf(meta_fp, "extract_sites\t%s\n", extract_path);
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

static void build_state_from_flare(
    bcf_hdr_t* ahdr,
    bcf1_t* arec,
    int flare_raw_n_samples,
    const std::vector<int>& flare_raw_indices,
    int n_ancestries,
    int n_words,
    int32_t** an1_arr,
    int* nan1_arr,
    int32_t** an2_arr,
    int* nan2_arr,
    AncestryState& state
) {
    int n_an1 = bcf_get_format_int32(ahdr, arec, "AN1", an1_arr, nan1_arr);
    int n_an2 = bcf_get_format_int32(ahdr, arec, "AN2", an2_arr, nan2_arr);

    if (n_an1 <= 0 || n_an2 <= 0) {
        die("FLARE VCF must contain FORMAT/AN1 and FORMAT/AN2");
    }

    if (n_an1 != flare_raw_n_samples || n_an2 != flare_raw_n_samples) {
        die("FLARE VCF must contain scalar FORMAT/AN1 and FORMAT/AN2");
    }

    int n_output_samples = static_cast<int>(flare_raw_indices.size());
    state.masks.assign(
        static_cast<size_t>(n_ancestries),
        std::vector<uint64_t>(static_cast<size_t>(n_words), 0)
    );
    state.hap_ancestry.assign(static_cast<size_t>(2 * n_output_samples), -1);

    for (int out_i = 0; out_i < n_output_samples; ++out_i) {
        int raw_i = flare_raw_indices[out_i];
        int32_t raw_a1 = (*an1_arr)[raw_i];
        int32_t raw_a2 = (*an2_arr)[raw_i];

        if (raw_a1 == bcf_int32_missing || raw_a1 == bcf_int32_vector_end) {
            die("missing FORMAT/AN1 in FLARE record at sample index %d", raw_i);
        }

        if (raw_a2 == bcf_int32_missing || raw_a2 == bcf_int32_vector_end) {
            die("missing FORMAT/AN2 in FLARE record at sample index %d", raw_i);
        }

        int a1 = static_cast<int>(raw_a1);
        int a2 = static_cast<int>(raw_a2);

        if (a1 < 0 || a1 >= n_ancestries) {
            die("FORMAT/AN1 value out of range at sample index %d: %d", raw_i, a1);
        }

        if (a2 < 0 || a2 >= n_ancestries) {
            die("FORMAT/AN2 value out of range at sample index %d: %d", raw_i, a2);
        }

        uint32_t hap0 = static_cast<uint32_t>(2 * out_i);
        uint32_t hap1 = static_cast<uint32_t>(2 * out_i + 1);
        set_bit(state.masks[a1], hap0);
        set_bit(state.masks[a2], hap1);
        state.hap_ancestry[hap0] = static_cast<int8_t>(a1);
        state.hap_ancestry[hap1] = static_cast<int8_t>(a2);
    }
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

static void update_open_block_from_flare(
    FILE* anc_fp,
    const char* anc_bin_path,
    FILE* anc_mks_fp,
    const char* anc_mks_path,
    FILE* anc_idx_fp,
    OpenAncestryBlock& block,
    uint32_t& n_blocks_written,
    int n_ancestries,
    int n_words,
    int flare_geno_rid,
    const char* flare_chr,
    int64_t interval_start,
    int64_t interval_end,
    AncestryState& new_state
) {
    if (interval_start > interval_end) {
        die(
            "invalid LAI interval for %s: %lld > %lld",
            flare_chr,
            static_cast<long long>(interval_start),
            static_cast<long long>(interval_end)
        );
    }

    if (!block.active) {
        block.active = true;
        block.block_id = n_blocks_written;
        block.geno_rid = flare_geno_rid;
        block.chr = flare_chr;
        block.start_pos = interval_start;
        block.end_pos = interval_end;
        block.state = std::move(new_state);
        return;
    }

    bool same_chromosome = block.geno_rid == flare_geno_rid;
    bool unchanged = same_chromosome && same_masks(block.state.masks, new_state.masks);

    if (unchanged) {
        block.end_pos = std::max(block.end_pos, interval_end);
        return;
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
    block.geno_rid = flare_geno_rid;
    block.chr = flare_chr;
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
    int genotype_raw_n_samples,
    const std::vector<int>& genotype_raw_indices,
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

    if (ngt % genotype_raw_n_samples != 0) {
        die("GT field length is not divisible by sample count at %s:%lld", chr, static_cast<long long>(pos));
    }

    int ploidy = ngt / genotype_raw_n_samples;
    if (ploidy < 2) {
        die("expected diploid GT at %s:%lld", chr, static_cast<long long>(pos));
    }

    carriers_by_alt.assign(static_cast<size_t>(grec->n_allele), std::vector<CarrierHap>{});

    int n_output_samples = static_cast<int>(genotype_raw_indices.size());
    for (int out_i = 0; out_i < n_output_samples; ++out_i) {
        int raw_i = genotype_raw_indices[out_i];
        const int32_t* sample_gt = *gt_arr + static_cast<size_t>(raw_i) * ploidy;
        int32_t g0 = sample_gt[0];
        int32_t g1 = sample_gt[1];

        validate_extra_ploidy(sample_gt, ploidy, chr, pos, raw_i);

        if (is_gt_vector_end(g0) || is_gt_vector_end(g1) ||
            bcf_gt_is_missing(g0) || bcf_gt_is_missing(g1)) {
            die(
                "missing genotype at %s:%lld sample index %d",
                chr,
                static_cast<long long>(pos),
                raw_i
            );
        }

        int allele0 = bcf_gt_allele(g0);
        int allele1 = bcf_gt_allele(g1);

        if (!bcf_gt_is_phased(g1)) {
            die(
                "unphased genotype at %s:%lld sample index %d",
                chr,
                static_cast<long long>(pos),
                raw_i
            );
        }

        add_alt_carrier(
            carriers_by_alt,
            block.state,
            allele0,
            static_cast<uint32_t>(2 * out_i)
        );
        add_alt_carrier(
            carriers_by_alt,
            block.state,
            allele1,
            static_cast<uint32_t>(2 * out_i + 1)
        );
    }
}

static int read_next_record(htsFile* fp, bcf_hdr_t* hdr, bcf1_t* rec, hts_itr_t* itr = nullptr) {
    int ret = itr ? bcf_itr_next(fp, itr, rec) : bcf_read(fp, hdr, rec);
    if (ret == 0) bcf_unpack(rec, BCF_UN_STR);
    return ret;
}

static bool add_contig_to_header_if_missing(
    bcf_hdr_t* hdr,
    const std::string& contig,
    const char* header_label
) {
    if (contig.empty()) return false;
    if (bcf_hdr_name2id(hdr, contig.c_str()) >= 0) return false;

    std::string line = "##contig=<ID=" + contig + ">";
    if (bcf_hdr_append(hdr, line.c_str()) != 0 || bcf_hdr_sync(hdr) != 0) {
        die(
            "failed to add inferred contig to %s header: %s",
            header_label,
            contig.c_str()
        );
    }
    return true;
}

static int add_contigs_from_header_if_missing(
    bcf_hdr_t* dst,
    bcf_hdr_t* src,
    const char* dst_label,
    const char* src_label
) {
    int added = 0;
    int n_contigs = src ? src->n[BCF_DT_CTG] : 0;
    for (int i = 0; i < n_contigs; ++i) {
        const char* contig = bcf_hdr_id2name(src, i);
        if (!contig || !*contig) continue;
        if (add_contig_to_header_if_missing(dst, contig, dst_label)) {
            ++added;
        }
    }
    if (added > 0) {
        std::fprintf(
            stderr,
            "Added %d %s contig(s) to %s header.\n",
            added,
            src_label,
            dst_label
        );
    }
    return added;
}

static int add_contigs_from_extract_if_missing(
    bcf_hdr_t* dst,
    const ExtractSites& extract_sites,
    const char* dst_label
) {
    if (!extract_sites.active) return 0;

    std::unordered_set<std::string> seen;
    int added = 0;
    for (const ExtractPosition& position : extract_sites.positions) {
        if (!seen.insert(position.chr).second) continue;
        if (add_contig_to_header_if_missing(dst, position.chr, dst_label)) {
            ++added;
        }
    }

    if (added > 0) {
        std::fprintf(
            stderr,
            "Added %d --extract contig(s) to %s header.\n",
            added,
            dst_label
        );
    }
    return added;
}

static bool extract_position_in_region(
    const ExtractPosition& position,
    const Region& region
) {
    if (!region.active) return true;
    return position.chr == region.chr &&
           position.pos >= region.start &&
           position.pos <= region.end;
}

static int extract_position_sort_rid(bcf_hdr_t* hdr, const std::string& chr) {
    int rid = bcf_hdr_name2id(hdr, chr.c_str());
    if (rid >= 0) return rid;
    return std::numeric_limits<int>::max();
}

static std::vector<ExtractPosition> sorted_extract_positions(
    const ExtractSites& extract_sites,
    bcf_hdr_t* ghdr,
    const Region& region
) {
    std::vector<ExtractPosition> positions;
    if (!extract_sites.active) return positions;

    positions.reserve(extract_sites.positions.size());
    for (const ExtractPosition& position : extract_sites.positions) {
        if (extract_position_in_region(position, region)) {
            positions.push_back(position);
        }
    }

    std::sort(
        positions.begin(),
        positions.end(),
        [ghdr](const ExtractPosition& a, const ExtractPosition& b) {
            int arid = extract_position_sort_rid(ghdr, a.chr);
            int brid = extract_position_sort_rid(ghdr, b.chr);
            if (arid != brid) return arid < brid;
            if (a.chr != b.chr) return a.chr < b.chr;
            return a.pos < b.pos;
        }
    );

    return positions;
}

static std::string region_contig_token(const std::string& chr) {
    if (chr.find_first_of(":-") == std::string::npos) {
        return chr;
    }
    return "{" + chr + "}";
}

static void append_region_query(
    std::string& query,
    const std::string& chr,
    int64_t start,
    int64_t end
) {
    if (!query.empty()) query.push_back(',');
    query += region_contig_token(chr);
    query.push_back(':');
    query += std::to_string(start);
    query.push_back('-');
    query += std::to_string(end);
}

static std::string build_extract_genotype_region_query(
    const ExtractSites& extract_sites,
    bcf_hdr_t* ghdr,
    const Region& region,
    uint64_t& n_positions,
    uint64_t& n_intervals
) {
    n_positions = 0;
    n_intervals = 0;
    std::vector<ExtractPosition> positions =
        sorted_extract_positions(extract_sites, ghdr, region);

    std::string query;
    int64_t last_pos = -1;
    std::string last_chr;
    std::string interval_chr;
    int64_t interval_start = 0;
    int64_t interval_end = 0;
    bool have_interval = false;

    auto flush_interval = [&]() {
        if (!have_interval) return;
        append_region_query(query, interval_chr, interval_start, interval_end);
        ++n_intervals;
    };

    for (const ExtractPosition& position : positions) {
        if (position.chr == last_chr && position.pos == last_pos) continue;

        if (!have_interval) {
            interval_chr = position.chr;
            interval_start = position.pos;
            interval_end = position.pos;
            have_interval = true;
        } else if (position.chr == interval_chr) {
            interval_end = position.pos;
        } else {
            flush_interval();
            interval_chr = position.chr;
            interval_start = position.pos;
            interval_end = position.pos;
            have_interval = true;
        }

        last_chr = position.chr;
        last_pos = position.pos;
        ++n_positions;
    }
    flush_interval();

    return query;
}

static std::string build_extract_flare_region_query(
    const ExtractSites& extract_sites,
    bcf_hdr_t* ghdr,
    const Region& region,
    uint64_t& n_contigs
) {
    n_contigs = 0;
    std::vector<ExtractPosition> positions =
        sorted_extract_positions(extract_sites, ghdr, region);

    std::string query;
    std::string last_chr;
    for (const ExtractPosition& position : positions) {
        if (position.chr == last_chr) continue;
        append_region_query(query, position.chr, 1, kMaxVcfCoordinate);
        last_chr = position.chr;
        ++n_contigs;
    }

    return query;
}

static int infer_contigs_from_records(
    const char* vcf_path,
    bcf_hdr_t* dst,
    const char* dst_label
) {
    htsFile* fp = bcf_open(vcf_path, "r");
    if (!fp) die("cannot reopen VCF for contig inference: %s", vcf_path);

    bcf_hdr_t* scan_hdr = bcf_hdr_read(fp);
    if (!scan_hdr) die("cannot read VCF header for contig inference");

    bcf1_t* rec = bcf_init();
    if (!rec) die("failed to allocate VCF record for contig inference");

    std::unordered_set<std::string> seen;
    int added = 0;
    int ret = 0;
    while ((ret = read_next_record(fp, scan_hdr, rec)) == 0) {
        const char* chr = bcf_hdr_id2name(scan_hdr, rec->rid);
        if (!chr || !*chr) {
            die("cannot resolve chromosome during contig inference");
        }
        std::string contig(chr);
        if (seen.insert(contig).second &&
            add_contig_to_header_if_missing(dst, contig, dst_label)) {
            ++added;
        }
    }

    if (ret < -1) {
        die("error while scanning VCF for contig inference");
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(scan_hdr);
    bcf_close(fp);

    std::fprintf(
        stderr,
        "Inferred %d contig(s) for %s header by scanning %s.\n",
        added,
        dst_label,
        vcf_path
    );
    return added;
}

static bool ends_with(const std::string& value, const char* suffix) {
    size_t n = std::strlen(suffix);
    return value.size() >= n &&
           value.compare(value.size() - n, n, suffix) == 0;
}

static bool looks_indexable_variant_path(const char* path) {
    std::string value(path);
    return ends_with(value, ".bcf") ||
           ends_with(value, ".vcf.gz") ||
           ends_with(value, ".vcf.bgz") ||
           ends_with(value, ".bcf.gz");
}

static bool init_synced_region_reader(
    bcf_srs_t** out,
    const char* path,
    const std::string& query,
    const char* label,
    const char* query_label = nullptr
) {
    const char* shown_query = query_label ? query_label : query.c_str();

    if (!looks_indexable_variant_path(path)) {
        std::fprintf(
            stderr,
            "WARNING: %s is not bgzip-compressed VCF/BCF; falling back to scan/filter.\n",
            label
        );
        return false;
    }

    bcf_srs_t* reader = bcf_sr_init();
    if (!reader) die("failed to initialize synced reader for %s", label);

    if (bcf_sr_set_regions(reader, query.c_str(), 0) != 0) {
        std::fprintf(
            stderr,
            "WARNING: cannot set %s region %s; falling back to scan/filter.\n",
            label,
            shown_query
        );
        bcf_sr_destroy(reader);
        return false;
    }

    if (!bcf_sr_add_reader(reader, path)) {
        std::fprintf(
            stderr,
            "WARNING: cannot use indexed %s reader for %s (%s); falling back to scan/filter.\n",
            label,
            shown_query,
            bcf_sr_strerror(reader->errnum)
        );
        bcf_sr_destroy(reader);
        return false;
    }

    *out = reader;
    return true;
}

static int read_next_synced_record(bcf_srs_t* reader, bcf1_t* rec) {
    while (bcf_sr_next_line(reader)) {
        bcf1_t* line = bcf_sr_get_line(reader, 0);
        if (!line) continue;
        if (!bcf_copy(rec, line)) die("failed to copy synced VCF record");
        bcf_unpack(rec, BCF_UN_STR);
        return 0;
    }
    return -1;
}

static bool record_in_region(bcf1_t* rec, const Region& region) {
    if (!region.active) return true;
    if (rec->rid != region.geno_rid) return false;
    int64_t pos = static_cast<int64_t>(rec->pos) + 1;
    return pos >= region.start && pos <= region.end;
}

static int read_next_genotype_record(
    htsFile* fp,
    bcf_hdr_t* hdr,
    bcf1_t* rec,
    bcf_srs_t* synced_reader,
    const Region& region
) {
    if (synced_reader) {
        return read_next_synced_record(synced_reader, rec);
    }

    int ret = 0;
    while ((ret = read_next_record(fp, hdr, rec)) == 0) {
        if (record_in_region(rec, region)) return 0;
    }
    return ret;
}

static bool read_next_lai_record(
    htsFile* fp,
    bcf_hdr_t* ahdr,
    bcf_hdr_t* ghdr,
    bcf1_t* rec,
    bcf_srs_t* synced_reader,
    int flare_raw_n_samples,
    const std::vector<int>& flare_raw_indices,
    int n_ancestries,
    int n_words,
    int32_t** an1_arr,
    int* nan1_arr,
    int32_t** an2_arr,
    int* nan2_arr,
    std::unordered_set<std::string>& warned_missing_flare_contigs,
    LaiRecord& out
) {
    while ((synced_reader ? read_next_synced_record(synced_reader, rec)
                          : read_next_record(fp, ahdr, rec)) == 0) {
        const char* a_chr = bcf_hdr_id2name(ahdr, rec->rid);
        int a_geno_rid = bcf_hdr_name2id(ghdr, a_chr);

        if (a_geno_rid < 0) {
            std::string contig(a_chr);
            if (warned_missing_flare_contigs.insert(contig).second) {
                std::fprintf(
                    stderr,
                    "WARNING: skipping FLARE contig absent from genotype header: %s\n",
                    a_chr
                );
            }
            continue;
        }

        out.valid = true;
        out.geno_rid = a_geno_rid;
        out.chr = a_chr;
        out.pos = static_cast<int64_t>(rec->pos) + 1;

        build_state_from_flare(
            ahdr,
            rec,
            flare_raw_n_samples,
            flare_raw_indices,
            n_ancestries,
            n_words,
            an1_arr,
            nan1_arr,
            an2_arr,
            nan2_arr,
            out.state
        );

        return true;
    }

    out = LaiRecord{};
    return false;
}

static void print_usage(const char* prog) {
    std::fprintf(
        stderr,
        "Usage:\n"
        "  %s genotype.phased.vcf.gz flare.anc.vcf.gz n_ancestries rare_threshold|auto out_prefix [chr:start-end] [--keep samples.txt] [--extract sites.pvar|sites.vcf]\n\n"
        "Example:\n"
        "  %s chr1.phased.vcf.gz chr1.flare.anc.vcf.gz 3 auto chr1\n"
        "  %s chr22.phased.vcf.gz chr22.flare.anc.vcf.gz 5 512 chr22.1 chr22:1-50000000\n",
        prog,
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
    const char* flare_vcf = argv[2];
    int n_ancestries = std::atoi(argv[3]);
    int rare_threshold = parse_rare_threshold_arg(argv[4]);
    const char* out_prefix = argv[5];
    Region region;
    std::string keep_path;
    std::string extract_path;
    for (int argi = 6; argi < argc; ++argi) {
        std::string arg = argv[argi];
        if (arg == "--keep") {
            if (argi + 1 >= argc) die("--keep requires a value");
            keep_path = argv[++argi];
        } else if (arg == "--extract") {
            if (argi + 1 >= argc) die("--extract requires a value");
            extract_path = argv[++argi];
        } else if (!arg.empty() && arg[0] == '-') {
            die("unknown option: %s", arg.c_str());
        } else {
            if (region.active) die("multiple region arguments supplied");
            region = parse_region_string(arg);
        }
    }

    if (n_ancestries <= 0 || n_ancestries > 32) {
        die("n_ancestries must be in [1, 32]");
    }

    htsFile* gfp = bcf_open(geno_vcf, "r");
    htsFile* afp = bcf_open(flare_vcf, "r");

    if (!gfp || !afp) {
        die("cannot open input VCF/BCF");
    }

    bcf_hdr_t* ghdr = bcf_hdr_read(gfp);
    bcf_hdr_t* ahdr = bcf_hdr_read(afp);

    if (!ghdr || !ahdr) {
        die("cannot read input headers");
    }

    SampleSelection sample_selection = build_sample_selection(ghdr, ahdr, keep_path);
    int genotype_raw_n_samples = sample_selection.genotype_raw_sample_count;
    int flare_raw_n_samples = sample_selection.flare_raw_sample_count;
    int n_samples = static_cast<int>(sample_selection.sample_ids.size());
    if (n_samples <= 0) {
        die("retained zero samples");
    }

    ExtractSites extract_sites = load_extract_sites(extract_path);
    if (extract_sites.active) {
        std::fprintf(
            stderr,
            "Loaded --extract site list: retaining up to %llu split allele(s).\n",
            static_cast<unsigned long long>(extract_sites.allele_keys.size())
        );
    }

    if (ghdr->n[BCF_DT_CTG] == 0) {
        std::fprintf(
            stderr,
            "WARNING: genotype VCF header has no ##contig lines; adding FLARE/header region contigs to the in-memory genotype header.\n"
        );
        add_contigs_from_header_if_missing(ghdr, ahdr, "genotype", "FLARE");
        add_contigs_from_extract_if_missing(ghdr, extract_sites, "genotype");
        if (!region.active && ghdr->n[BCF_DT_CTG] == 0) {
            infer_contigs_from_records(geno_vcf, ghdr, "genotype");
        }
    }

    bcf_srs_t* genotype_region_reader = nullptr;
    bcf_srs_t* flare_region_reader = nullptr;
    bool skip_genotype_loop = false;
    bool bounded_extract_reader = false;
    std::string progress_scope;

    if (region.active) {
        add_contig_to_header_if_missing(ghdr, region.chr, "genotype");
        region.geno_rid = bcf_hdr_name2id(ghdr, region.chr.c_str());
        if (region.geno_rid < 0) {
            die("failed to add region chromosome to genotype header: %s", region.chr.c_str());
        }
    }

    if (extract_sites.active) {
        uint64_t n_extract_positions = 0;
        uint64_t n_extract_intervals = 0;
        std::string genotype_query = build_extract_genotype_region_query(
            extract_sites,
            ghdr,
            region,
            n_extract_positions,
            n_extract_intervals
        );

        if (n_extract_positions == 0) {
            skip_genotype_loop = true;
            progress_scope = "--extract bounded intervals";
            std::fprintf(
                stderr,
                "No --extract positions overlap the selected conversion scope; skipping genotype scan.\n"
            );
        } else {
            bool have_genotype_extract_reader = init_synced_region_reader(
                &genotype_region_reader,
                geno_vcf,
                genotype_query,
                "genotype VCF",
                "--extract bounded intervals"
            );

            if (have_genotype_extract_reader) {
                bounded_extract_reader = true;
                progress_scope = "--extract bounded intervals";
                std::fprintf(
                    stderr,
                    "Using indexed --extract bounded genotype reader: %llu target position(s) across %llu interval(s), %llu split allele(s).\n",
                    static_cast<unsigned long long>(n_extract_positions),
                    static_cast<unsigned long long>(n_extract_intervals),
                    static_cast<unsigned long long>(extract_sites.allele_keys.size())
                );
            } else {
                std::fprintf(
                    stderr,
                    "WARNING: --extract could not use genotype VCF bounded random access; scanning records and filtering by CHROM/POS/REF/ALT.\n"
                );
            }
        }
    }

    if (!extract_sites.active && region.active) {
        std::string flare_query =
            region.chr + ":1-" + std::to_string(kMaxVcfCoordinate);
        bool have_genotype_region_reader = init_synced_region_reader(
            &genotype_region_reader,
            geno_vcf,
            region.label,
            "genotype VCF"
        );
        if (have_genotype_region_reader) {
            progress_scope = region.label;
        }
        bool have_flare_region_reader = init_synced_region_reader(
            &flare_region_reader,
            flare_vcf,
            flare_query,
            "FLARE VCF"
        );

        if (!have_genotype_region_reader || !have_flare_region_reader) {
            if (genotype_region_reader) {
                bcf_sr_destroy(genotype_region_reader);
                genotype_region_reader = nullptr;
            }
            if (flare_region_reader) {
                bcf_sr_destroy(flare_region_reader);
                flare_region_reader = nullptr;
            }
            std::fprintf(
                stderr,
                "WARNING: indexed region reader unavailable; scanning inputs and filtering %s.\n",
                region.label.c_str()
            );
        }
    } else if (region.active) {
        std::string flare_query =
            region.chr + ":1-" + std::to_string(kMaxVcfCoordinate);
        bool have_flare_region_reader = init_synced_region_reader(
            &flare_region_reader,
            flare_vcf,
            flare_query,
            "FLARE VCF"
        );
        if (!have_flare_region_reader) {
            std::fprintf(
                stderr,
                "WARNING: indexed FLARE reader unavailable; scanning FLARE records while using --extract filter.\n"
            );
        }
    } else if (extract_sites.active && !skip_genotype_loop) {
        uint64_t n_extract_contigs = 0;
        std::string flare_query = build_extract_flare_region_query(
            extract_sites,
            ghdr,
            region,
            n_extract_contigs
        );
        bool have_flare_extract_reader = init_synced_region_reader(
            &flare_region_reader,
            flare_vcf,
            flare_query,
            "FLARE VCF",
            "--extract contig list"
        );
        if (have_flare_extract_reader) {
            std::fprintf(
                stderr,
                "Using indexed FLARE reader for --extract: %llu contig(s).\n",
                static_cast<unsigned long long>(n_extract_contigs)
            );
        } else {
            std::fprintf(
                stderr,
                "WARNING: indexed FLARE reader unavailable for --extract contigs; scanning FLARE records.\n"
            );
        }
    }

    ProgressReporter progress(
        "genotype VCF",
        geno_vcf,
        "split variants",
        progress_scope.empty() ? (region.active ? region.label.c_str() : nullptr) : progress_scope.c_str()
    );

    if (rare_threshold < 0) {
        rare_threshold = default_rare_threshold_from_samples(n_samples);
        std::fprintf(
            stderr,
            "Using rare_threshold=%d = ceil(%d / 32) from retained sample count.\n",
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
        samples_path,
        meta_path,
        geno_vcf,
        flare_vcf,
        sample_selection.sample_ids,
        n_haps,
        n_words,
        n_ancestries,
        rare_threshold,
        region.active ? region.label.c_str() : nullptr,
        keep_path.empty() ? nullptr : keep_path.c_str(),
        extract_path.empty() ? nullptr : extract_path.c_str()
    );

    bcf1_t* grec = bcf_init();
    bcf1_t* arec = bcf_init();

    int32_t* gt_arr = nullptr;
    int ngt_arr = 0;
    int32_t* an1_arr = nullptr;
    int nan1_arr = 0;
    int32_t* an2_arr = nullptr;
    int nan2_arr = 0;

    OpenAncestryBlock block;
    std::vector<std::vector<CarrierHap>> carriers_by_alt;
    std::unordered_set<std::string> warned_missing_flare_contigs;

    uint32_t n_blocks_written = 0;
    uint64_t common_index = 0;
    uint64_t rare_index = 0;
    uint32_t global_variant_index = 0;

    LaiRecord lai_record;
    bool has_lai_record = false;
    if (!skip_genotype_loop) {
        has_lai_record = read_next_lai_record(
            afp,
            ahdr,
            ghdr,
            arec,
            flare_region_reader,
            flare_raw_n_samples,
            sample_selection.flare_raw_indices,
            n_ancestries,
            n_words,
            &an1_arr,
            &nan1_arr,
            &an2_arr,
            &nan2_arr,
            warned_missing_flare_contigs,
            lai_record
        );
    }
    int last_flare_geno_rid = -1;
    int64_t last_flare_pos = 0;

    AncestryState pre_region_state;
    int pre_region_rid = -1;
    std::string pre_region_chr;
    bool has_pre_region_state = false;

    while (!skip_genotype_loop &&
           read_next_genotype_record(gfp, ghdr, grec, genotype_region_reader, region) == 0) {
        progress.record_scanned();

        const char* g_chr = bcf_hdr_id2name(ghdr, grec->rid);
        int64_t g_pos = static_cast<int64_t>(grec->pos) + 1;

        while (has_lai_record) {
            bool current_block_covers =
                block.active &&
                block.geno_rid == grec->rid &&
                g_pos <= block.end_pos;

            if (lai_record.geno_rid > grec->rid ||
                (lai_record.geno_rid == grec->rid && current_block_covers)) {
                break;
            }

            LaiRecord interval_record = std::move(lai_record);
            has_lai_record = read_next_lai_record(
                afp,
                ahdr,
                ghdr,
                arec,
                flare_region_reader,
                flare_raw_n_samples,
                sample_selection.flare_raw_indices,
                n_ancestries,
                n_words,
                &an1_arr,
                &nan1_arr,
                &an2_arr,
                &nan2_arr,
                warned_missing_flare_contigs,
                lai_record
            );

            while (has_lai_record &&
                   lai_record.geno_rid == interval_record.geno_rid &&
                   lai_record.pos == interval_record.pos) {
                interval_record = std::move(lai_record);
                has_lai_record = read_next_lai_record(
                    afp,
                    ahdr,
                    ghdr,
                    arec,
                    flare_region_reader,
                    flare_raw_n_samples,
                    sample_selection.flare_raw_indices,
                    n_ancestries,
                    n_words,
                    &an1_arr,
                    &nan1_arr,
                    &an2_arr,
                    &nan2_arr,
                    warned_missing_flare_contigs,
                    lai_record
                );
            }

            if (last_flare_geno_rid > interval_record.geno_rid) {
                die("FLARE records are not sorted by genotype header contig order");
            }

            int64_t interval_start = 1;
            if (last_flare_geno_rid == interval_record.geno_rid) {
                if (interval_record.pos <= last_flare_pos) {
                    die(
                        "FLARE records are not strictly increasing at %s:%lld",
                        interval_record.chr.c_str(),
                        static_cast<long long>(interval_record.pos)
                    );
                }
                interval_start = last_flare_pos + 1;
            }

            int64_t interval_end = interval_record.pos;
            int64_t output_start = interval_start;
            int64_t output_end = interval_end;
            bool emit_interval = true;

            if (region.active) {
                if (interval_record.geno_rid != region.geno_rid) {
                    emit_interval = false;
                } else if (interval_end < region.start) {
                    pre_region_state = interval_record.state;
                    pre_region_rid = interval_record.geno_rid;
                    pre_region_chr = interval_record.chr;
                    has_pre_region_state = true;
                    emit_interval = false;
                } else if (interval_start > region.end) {
                    emit_interval = false;
                } else {
                    output_start = std::max(interval_start, region.start);
                    output_end = std::min(interval_end, region.end);
                    emit_interval = output_start <= output_end;
                }
            }

            if (!region.active && extract_sites.active && bounded_extract_reader &&
                extract_sites.contigs.count(interval_record.chr) == 0) {
                emit_interval = false;
            }

            if (emit_interval) {
                update_open_block_from_flare(
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
                    output_start,
                    output_end,
                    interval_record.state
                );
            }

            last_flare_geno_rid = interval_record.geno_rid;
            last_flare_pos = interval_record.pos;
        }

        if (region.active && !block.active && !has_lai_record &&
            has_pre_region_state && pre_region_rid == grec->rid) {
            AncestryState state = pre_region_state;
            update_open_block_from_flare(
                anc_fp,
                anc_bin.c_str(),
                anc_mks_fp,
                anc_mks.c_str(),
                anc_idx_fp,
                block,
                n_blocks_written,
                n_ancestries,
                n_words,
                pre_region_rid,
                pre_region_chr.c_str(),
                region.start,
                region.end,
                state
            );
        }

        if (block.active && block.geno_rid == grec->rid && g_pos > block.end_pos &&
            (!has_lai_record || lai_record.geno_rid > grec->rid)) {
            block.end_pos = g_pos;
        }

        if (!block.active || block.geno_rid != grec->rid || g_pos > block.end_pos) {
            progress.maybe_report(global_variant_index, common_index, rare_index);
            continue;
        }

        if (grec->n_allele < 2) {
            progress.maybe_report(global_variant_index, common_index, rare_index);
            continue;
        }

        const char* raw_id = grec->d.id;
        const char* ref = grec->d.allele[0];
        if (!extract_record_position_matches(extract_sites, g_chr, g_pos, ref)) {
            progress.maybe_report(global_variant_index, common_index, rare_index);
            continue;
        }

        process_genotypes(
            ghdr,
            grec,
            block,
            genotype_raw_n_samples,
            sample_selection.genotype_raw_indices,
            &gt_arr,
            &ngt_arr,
            carriers_by_alt
        );

        if (carriers_by_alt.empty()) {
            progress.maybe_report(global_variant_index, common_index, rare_index);
            continue;
        }

        for (int alt_idx = 1; alt_idx < grec->n_allele; ++alt_idx) {
            const char* alt = grec->d.allele[alt_idx];
            if (!extract_allele_matches(extract_sites, g_chr, g_pos, ref, alt)) {
                continue;
            }

            if (global_variant_index == std::numeric_limits<uint32_t>::max()) {
                die("global split variant index exceeds uint32_t limit");
            }

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

    if (gt_arr) std::free(gt_arr);
    if (an1_arr) std::free(an1_arr);
    if (an2_arr) std::free(an2_arr);

    bcf_destroy(grec);
    bcf_destroy(arec);
    if (genotype_region_reader) bcf_sr_destroy(genotype_region_reader);
    if (flare_region_reader) bcf_sr_destroy(flare_region_reader);
    bcf_hdr_destroy(ghdr);
    bcf_hdr_destroy(ahdr);
    bcf_close(gfp);
    bcf_close(afp);

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
    std::fprintf(stderr, "Global split variants: %u\n", global_variant_index);
    std::fprintf(stderr, "Common variants:       %llu\n", static_cast<unsigned long long>(common_index));
    std::fprintf(stderr, "Rare variants:         %llu\n", static_cast<unsigned long long>(rare_index));
    std::fprintf(stderr, "Ancestry blocks:       %u\n", n_blocks_written);

    return 0;
}
