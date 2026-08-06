// tractor_hybrid_to_vcf.cpp
// designed by Kai, implemented by codex
//
// Reconstruct a split-biallelic phased VCF from the ancestry-aware packed
// backend files produced by flare_subset_to_tractor_hybrid.

#include <htslib/bgzf.h>
#include <htslib/tbx.h>

#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

struct RareCarrierPacked {
    uint32_t pos_index;
    uint32_t anc_hap;
};

struct Meta {
    uint64_t n_samples = 0;
    uint64_t n_haps = 0;
    uint64_t n_words = 0;
    uint64_t n_ancestries = 0;
};

struct VariantRecord {
    bool valid = false;
    bool common = false;
    uint32_t global_variant_index = 0;
    std::string chr;
    int64_t pos = 0;
    std::string id;
    std::string ref;
    std::string alt;
    uint64_t offset = 0;
    uint32_t n_carriers = 0;
};

[[noreturn]] static void die(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    std::fputs("ERROR: ", stderr);
    std::vfprintf(stderr, fmt, args);
    std::fputc('\n', stderr);
    va_end(args);
    std::exit(1);
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

static uint64_t parse_u64(const std::string& value, const char* field) {
    char* end = nullptr;
    unsigned long long parsed = std::strtoull(value.c_str(), &end, 10);
    if (!end || *end != '\0') {
        die("invalid unsigned integer for %s: %s", field, value.c_str());
    }
    return static_cast<uint64_t>(parsed);
}

static FILE* open_file_or_die(const std::string& path, const char* mode) {
    FILE* fp = std::fopen(path.c_str(), mode);
    if (!fp) die("cannot open %s", path.c_str());

    if (setvbuf(fp, nullptr, _IOFBF, 1 << 24) != 0) {
        die("setvbuf failed for %s", path.c_str());
    }

    return fp;
}

static bool ends_with(const std::string& value, const std::string& suffix) {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static void bgzf_write_all(BGZF* fp, const std::string& text, const std::string& path) {
    ssize_t wrote = bgzf_write(fp, text.data(), text.size());
    if (wrote < 0 || static_cast<size_t>(wrote) != text.size()) {
        die("failed writing BGZF output %s", path.c_str());
    }
}

static void seek_or_die(FILE* fp, uint64_t offset, const char* path) {
    if (offset > static_cast<uint64_t>(std::numeric_limits<off_t>::max())) {
        die("offset too large for this platform in %s", path);
    }

    if (fseeko(fp, static_cast<off_t>(offset), SEEK_SET) != 0) {
        die("seek failed in %s", path);
    }
}

static void read_exact(FILE* fp, void* data, size_t bytes, const char* path) {
    if (bytes == 0) return;
    if (std::fread(data, 1, bytes, fp) != bytes) {
        die("unexpected EOF while reading %s", path);
    }
}

static uint32_t read_u32_le(FILE* fp, const char* path) {
    uint8_t bytes[4];
    read_exact(fp, bytes, sizeof(bytes), path);
    return static_cast<uint32_t>(bytes[0]) |
           (static_cast<uint32_t>(bytes[1]) << 8) |
           (static_cast<uint32_t>(bytes[2]) << 16) |
           (static_cast<uint32_t>(bytes[3]) << 24);
}

static uint64_t read_u64_le(FILE* fp, const char* path) {
    uint8_t bytes[8];
    read_exact(fp, bytes, sizeof(bytes), path);
    return static_cast<uint64_t>(bytes[0]) |
           (static_cast<uint64_t>(bytes[1]) << 8) |
           (static_cast<uint64_t>(bytes[2]) << 16) |
           (static_cast<uint64_t>(bytes[3]) << 24) |
           (static_cast<uint64_t>(bytes[4]) << 32) |
           (static_cast<uint64_t>(bytes[5]) << 40) |
           (static_cast<uint64_t>(bytes[6]) << 48) |
           (static_cast<uint64_t>(bytes[7]) << 56);
}

static int64_t read_i64_le(FILE* fp, const char* path) {
    return static_cast<int64_t>(read_u64_le(fp, path));
}

static std::string read_string(FILE* fp, const char* path) {
    uint32_t len = read_u32_le(fp, path);
    std::string value(len, '\0');
    if (len > 0) {
        read_exact(fp, value.data(), len, path);
    }
    return value;
}

static void read_idx_header(FILE* fp, const std::string& path, const char magic[8]) {
    char got[8];
    read_exact(fp, got, sizeof(got), path.c_str());
    if (std::memcmp(got, magic, sizeof(got)) != 0) {
        die("bad index magic in %s", path.c_str());
    }
}

static Meta read_meta(const std::string& path) {
    std::ifstream in(path);
    if (!in) die("cannot open %s", path.c_str());

    std::unordered_map<std::string, std::string> kv;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> fields = split_tabs(line);
        if (fields.size() >= 2) kv[fields[0]] = fields[1];
    }

    Meta meta;
    if (!kv.count("n_samples")) die("missing n_samples in %s", path.c_str());
    if (!kv.count("n_haps")) die("missing n_haps in %s", path.c_str());
    if (!kv.count("n_words")) die("missing n_words in %s", path.c_str());
    if (!kv.count("n_ancestries")) die("missing n_ancestries in %s", path.c_str());

    meta.n_samples = parse_u64(kv["n_samples"], "n_samples");
    meta.n_haps = parse_u64(kv["n_haps"], "n_haps");
    meta.n_words = parse_u64(kv["n_words"], "n_words");
    meta.n_ancestries = parse_u64(kv["n_ancestries"], "n_ancestries");

    if (meta.n_samples == 0) die("n_samples must be positive");
    if (meta.n_haps != meta.n_samples * 2) die("n_haps must equal 2 * n_samples");
    if (meta.n_words != (meta.n_haps + 63) / 64) die("n_words does not match n_haps");
    if (meta.n_ancestries == 0 || meta.n_ancestries > 32) die("n_ancestries must be in [1, 32]");

    return meta;
}

static std::vector<std::string> read_samples(const std::string& path, uint64_t expected_n) {
    std::ifstream in(path);
    if (!in) die("cannot open %s", path.c_str());

    std::vector<std::string> samples;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) samples.push_back(line);
    }

    if (samples.size() != expected_n) {
        die(
            "sample count mismatch in %s: got %llu, expected %llu",
            path.c_str(),
            static_cast<unsigned long long>(samples.size()),
            static_cast<unsigned long long>(expected_n)
        );
    }

    return samples;
}

class IndexReader {
public:
    IndexReader(const std::string& path, bool common)
        : path_(path), common_(common), fp_(open_file_or_die(path, "rb")) {
        const char common_magic[8] = {'T', 'R', 'C', 'M', 'M', 'K', 'S', '1'};
        const char rare_magic[8] = {'T', 'R', 'R', 'A', 'M', 'K', 'S', '1'};
        read_idx_header(fp_, path_, common_ ? common_magic : rare_magic);

        read_next();
    }

    ~IndexReader() {
        if (fp_) std::fclose(fp_);
    }

    const VariantRecord& current() const {
        return current_;
    }

    void advance() {
        read_next();
    }

private:
    void read_next() {
        current_ = VariantRecord{};

        int first = std::fgetc(fp_);
        if (first == EOF) {
            return;
        }
        if (std::ungetc(first, fp_) == EOF) {
            die("ungetc failed for %s", path_.c_str());
        }

        current_.valid = true;
        current_.common = common_;
        (void)read_u64_le(fp_, path_.c_str());
        current_.global_variant_index = read_u32_le(fp_, path_.c_str());
        current_.chr = read_string(fp_, path_.c_str());
        current_.pos = read_i64_le(fp_, path_.c_str());
        current_.id = read_string(fp_, path_.c_str());
        current_.ref = read_string(fp_, path_.c_str());
        current_.alt = read_string(fp_, path_.c_str());
        (void)read_u32_le(fp_, path_.c_str());

        if (common_) {
            (void)read_u32_le(fp_, path_.c_str());
            current_.offset = read_u64_le(fp_, path_.c_str());
            (void)read_u32_le(fp_, path_.c_str());
        } else {
            current_.offset = read_u64_le(fp_, path_.c_str());
            current_.n_carriers = read_u32_le(fp_, path_.c_str());
            (void)read_u32_le(fp_, path_.c_str());
        }
    }

    std::string path_;
    bool common_;
    FILE* fp_ = nullptr;
    VariantRecord current_;
};

static inline bool bit_at(const std::vector<uint64_t>& words, uint64_t hap_id) {
    return (words[hap_id >> 6] >> (hap_id & 63)) & 1ULL;
}

static void read_common_bits(
    FILE* fp,
    const std::string& path,
    const VariantRecord& record,
    uint64_t n_words,
    std::vector<uint64_t>& words
) {
    words.assign(static_cast<size_t>(n_words), 0);
    seek_or_die(fp, record.offset, path.c_str());

    size_t got = std::fread(
        words.data(),
        sizeof(uint64_t),
        static_cast<size_t>(n_words),
        fp
    );

    if (got != n_words) {
        die("failed reading common bitset for global variant %u", record.global_variant_index);
    }
}

static void read_rare_haps(
    FILE* fp,
    const std::string& path,
    const VariantRecord& record,
    uint64_t n_haps,
    std::vector<uint8_t>& alt_haps
) {
    alt_haps.assign(static_cast<size_t>(n_haps), 0);
    seek_or_die(fp, record.offset, path.c_str());

    for (uint32_t i = 0; i < record.n_carriers; ++i) {
        RareCarrierPacked carrier{};
        size_t got = std::fread(&carrier, sizeof(RareCarrierPacked), 1, fp);
        if (got != 1) {
            die("failed reading rare carrier for global variant %u", record.global_variant_index);
        }

        if (carrier.pos_index != record.global_variant_index) {
            die(
                "rare carrier pos_index mismatch: carrier has %u, index has %u",
                carrier.pos_index,
                record.global_variant_index
            );
        }

        uint32_t hap_id = carrier.anc_hap & 0x07ffffffu;
        if (hap_id >= n_haps) {
            die("rare carrier hap_id %u exceeds n_haps", hap_id);
        }

        alt_haps[hap_id] = 1;
    }
}

static void write_vcf_header(BGZF* out_fp, const std::string& out_path, const std::vector<std::string>& samples) {
    std::string line;
    line += "##fileformat=VCFv4.2\n";
    line += "##source=tractor_hybrid_to_vcf\n";
    line += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype\">\n";
    line += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

    for (const std::string& sample : samples) {
        line += '\t';
        line += sample;
    }

    line += '\n';
    bgzf_write_all(out_fp, line, out_path);
}

static void write_common_variant(
    BGZF* out_fp,
    const std::string& out_path,
    const VariantRecord& record,
    const std::vector<uint64_t>& words,
    uint64_t n_samples
) {
    std::ostringstream line;
    line << record.chr << '\t'
         << record.pos << '\t'
         << record.id << '\t'
         << record.ref << '\t'
         << record.alt << "\t.\tPASS\t.\tGT";

    for (uint64_t i = 0; i < n_samples; ++i) {
        int a0 = bit_at(words, 2 * i) ? 1 : 0;
        int a1 = bit_at(words, 2 * i + 1) ? 1 : 0;
        line << '\t' << a0 << '|' << a1;
    }

    line << '\n';
    bgzf_write_all(out_fp, line.str(), out_path);
}

static void write_rare_variant(
    BGZF* out_fp,
    const std::string& out_path,
    const VariantRecord& record,
    const std::vector<uint8_t>& alt_haps,
    uint64_t n_samples
) {
    std::ostringstream line;
    line << record.chr << '\t'
         << record.pos << '\t'
         << record.id << '\t'
         << record.ref << '\t'
         << record.alt << "\t.\tPASS\t.\tGT";

    for (uint64_t i = 0; i < n_samples; ++i) {
        int a0 = alt_haps[2 * i] ? 1 : 0;
        int a1 = alt_haps[2 * i + 1] ? 1 : 0;
        line << '\t' << a0 << '|' << a1;
    }

    line << '\n';
    bgzf_write_all(out_fp, line.str(), out_path);
}

static void print_usage(const char* prog) {
    std::fprintf(
        stderr,
        "Usage:\n"
        "  %s packed_prefix out.vcf.gz\n\n"
        "Example:\n"
        "  %s chr1 chr1.roundtrip.vcf.gz\n",
        prog,
        prog
    );
}

int main(int argc, char** argv) {
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }

    std::string prefix = argv[1];
    std::string out_vcf = argv[2];

    if (!ends_with(out_vcf, ".gz")) {
        die("output path must end with .gz so tabix index naming is unambiguous");
    }

    std::string meta_path = prefix + ".meta";
    std::string samples_path = prefix + ".samples";
    std::string common_mks_path = prefix + ".common.variant.mks";
    std::string rare_mks_path = prefix + ".rare.variant.mks";
    std::string common_bin_path = prefix + ".common.geno.bin";
    std::string rare_bin_path = prefix + ".rare.carrier.bin";

    Meta meta = read_meta(meta_path);
    std::vector<std::string> samples = read_samples(samples_path, meta.n_samples);

    IndexReader common_idx(common_mks_path, true);
    IndexReader rare_idx(rare_mks_path, false);

    FILE* common_fp = open_file_or_die(common_bin_path, "rb");
    FILE* rare_fp = open_file_or_die(rare_bin_path, "rb");
    BGZF* out_fp = bgzf_open(out_vcf.c_str(), "w");
    if (!out_fp) die("cannot open BGZF output %s", out_vcf.c_str());

    write_vcf_header(out_fp, out_vcf, samples);

    std::vector<uint64_t> common_words;
    std::vector<uint8_t> rare_haps;
    uint64_t variants_written = 0;

    while (common_idx.current().valid || rare_idx.current().valid) {
        bool take_common = false;

        if (common_idx.current().valid && !rare_idx.current().valid) {
            take_common = true;
        } else if (common_idx.current().valid && rare_idx.current().valid) {
            uint32_t common_global = common_idx.current().global_variant_index;
            uint32_t rare_global = rare_idx.current().global_variant_index;
            if (common_global == rare_global) {
                die("global variant %u exists in both common and rare indexes", common_global);
            }
            take_common = common_global < rare_global;
        }

        if (take_common) {
            const VariantRecord& record = common_idx.current();
            read_common_bits(common_fp, common_bin_path, record, meta.n_words, common_words);
            write_common_variant(out_fp, out_vcf, record, common_words, meta.n_samples);
            common_idx.advance();
        } else {
            const VariantRecord& record = rare_idx.current();
            read_rare_haps(rare_fp, rare_bin_path, record, meta.n_haps, rare_haps);
            write_rare_variant(out_fp, out_vcf, record, rare_haps, meta.n_samples);
            rare_idx.advance();
        }

        ++variants_written;
    }

    std::fclose(common_fp);
    std::fclose(rare_fp);
    if (bgzf_close(out_fp) != 0) {
        die("failed closing BGZF output %s", out_vcf.c_str());
    }

    std::string tbi_path = out_vcf + ".tbi";
    std::remove(tbi_path.c_str());
    int index_ret = tbx_index_build(out_vcf.c_str(), 0, &tbx_conf_vcf);
    if (index_ret != 0) {
        die("failed building tabix index for %s", out_vcf.c_str());
    }

    std::fprintf(stderr, "Finished.\n");
    std::fprintf(stderr, "Split VCF variants written: %llu\n", static_cast<unsigned long long>(variants_written));
    std::fprintf(stderr, "Tabix index: %s\n", tbi_path.c_str());

    return 0;
}
