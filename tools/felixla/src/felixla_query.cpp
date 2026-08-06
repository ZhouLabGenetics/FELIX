// felixla_query.cpp
//
// Query ancestry-specific dosages from the FELIXla v0 packed layout.
// FELIXla v0 is binary-compatible with the tractor_hybrid packed files.

#include <algorithm>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

struct Meta {
    uint32_t format_version = 0;
    uint64_t n_samples = 0;
    uint64_t n_haps = 0;
    uint64_t n_words = 0;
    uint64_t n_ancestries = 0;
};

struct VariantRecord {
    bool common = false;
    uint32_t global_variant_index = 0;
    std::string chr;
    int64_t pos = 0;
    std::string id;
    std::string ref;
    std::string alt;
    uint32_t alt_index = 0;
    uint32_t block_id = 0;
    uint64_t payload_offset = 0;
    uint32_t n_carriers = 0;
    uint32_t mac = 0;
};

struct Args {
    std::string prefix;
    bool by_global_index = false;
    uint32_t global_index = 0;
    std::string chr;
    int64_t pos = 0;
    bool has_ref_alt = false;
    std::string ref;
    std::string alt;
    bool nonzero_only = false;
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

static void usage(FILE* out, const char* argv0) {
    std::fprintf(
        out,
        "Usage:\n"
        "  %s PREFIX CHR POS [REF ALT] [--nonzero-only]\n"
        "  %s PREFIX CHR:POS [REF ALT] [--nonzero-only]\n"
        "  %s PREFIX --global-index N [--nonzero-only]\n\n"
        "Output columns are variant metadata, sample, total ALT dosage (DSALL),\n"
        "and ancestry-specific ALT dosages DS1..DSk. DS1 is ancestry code 0.\n",
        argv0,
        argv0,
        argv0
    );
}

static bool starts_with(const std::string& value, const char* prefix) {
    size_t n = std::strlen(prefix);
    return value.size() >= n && value.compare(0, n, prefix) == 0;
}

static uint64_t parse_u64(const std::string& text, const char* label) {
    if (text.empty()) die("%s is empty", label);
    char* end = nullptr;
    unsigned long long value = std::strtoull(text.c_str(), &end, 10);
    if (end == text.c_str() || *end != '\0') {
        die("invalid unsigned integer for %s: %s", label, text.c_str());
    }
    return static_cast<uint64_t>(value);
}

static int64_t parse_i64(const std::string& text, const char* label) {
    if (text.empty()) die("%s is empty", label);
    char* end = nullptr;
    long long value = std::strtoll(text.c_str(), &end, 10);
    if (end == text.c_str() || *end != '\0') {
        die("invalid integer for %s: %s", label, text.c_str());
    }
    return static_cast<int64_t>(value);
}

static void parse_chr_pos(const std::string& text, std::string& chr, int64_t& pos) {
    size_t colon = text.find(':');
    if (colon == std::string::npos || colon == 0 || colon + 1 >= text.size()) {
        die("position query must look like CHR:POS: %s", text.c_str());
    }
    chr = text.substr(0, colon);
    pos = parse_i64(text.substr(colon + 1), "position");
    if (pos <= 0) die("position must be positive: %s", text.c_str());
}

static Args parse_args(int argc, char** argv) {
    if (argc == 1) {
        usage(stderr, argv[0]);
        std::exit(2);
    }

    std::vector<std::string> positional;
    Args args;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            usage(stdout, argv[0]);
            std::exit(0);
        } else if (arg == "--nonzero-only") {
            args.nonzero_only = true;
        } else if (arg == "--global-index") {
            if (i + 1 >= argc) die("--global-index requires an argument");
            uint64_t value = parse_u64(argv[++i], "global index");
            if (value > std::numeric_limits<uint32_t>::max()) {
                die("global index exceeds uint32_t limit: %llu", static_cast<unsigned long long>(value));
            }
            args.by_global_index = true;
            args.global_index = static_cast<uint32_t>(value);
        } else if (starts_with(arg, "--")) {
            die("unknown option: %s", arg.c_str());
        } else {
            positional.push_back(arg);
        }
    }

    if (positional.empty()) die("missing PREFIX");
    args.prefix = positional[0];

    if (args.by_global_index) {
        if (positional.size() != 1) {
            die("--global-index cannot be combined with CHR/POS arguments");
        }
        return args;
    }

    if (positional.size() == 2 || positional.size() == 4) {
        parse_chr_pos(positional[1], args.chr, args.pos);
        if (positional.size() == 4) {
            args.has_ref_alt = true;
            args.ref = positional[2];
            args.alt = positional[3];
        }
    } else if (positional.size() == 3 || positional.size() == 5) {
        args.chr = positional[1];
        args.pos = parse_i64(positional[2], "position");
        if (args.chr.empty()) die("chromosome is empty");
        if (args.pos <= 0) die("position must be positive");
        if (positional.size() == 5) {
            args.has_ref_alt = true;
            args.ref = positional[3];
            args.alt = positional[4];
        }
    } else {
        usage(stderr, argv[0]);
        std::exit(2);
    }

    return args;
}

static FILE* open_file_or_die(const std::string& path, const char* mode) {
    FILE* fp = std::fopen(path.c_str(), mode);
    if (!fp) die("cannot open %s", path.c_str());
    if (setvbuf(fp, nullptr, _IOFBF, 1 << 24) != 0) {
        die("setvbuf failed for %s", path.c_str());
    }
    return fp;
}

static void seek_or_die(FILE* fp, uint64_t offset, const std::string& path) {
    if (offset > static_cast<uint64_t>(std::numeric_limits<off_t>::max())) {
        die("offset too large for this platform in %s", path.c_str());
    }
    if (fseeko(fp, static_cast<off_t>(offset), SEEK_SET) != 0) {
        die("seek failed in %s", path.c_str());
    }
}

static void read_exact(FILE* fp, void* data, size_t bytes, const std::string& path) {
    if (bytes == 0) return;
    if (std::fread(data, 1, bytes, fp) != bytes) {
        die("unexpected EOF while reading %s", path.c_str());
    }
}

static uint32_t read_u32_le(FILE* fp, const std::string& path) {
    uint8_t b[4];
    read_exact(fp, b, sizeof(b), path);
    return static_cast<uint32_t>(b[0]) |
           (static_cast<uint32_t>(b[1]) << 8) |
           (static_cast<uint32_t>(b[2]) << 16) |
           (static_cast<uint32_t>(b[3]) << 24);
}

static uint64_t read_u64_le(FILE* fp, const std::string& path) {
    uint8_t b[8];
    read_exact(fp, b, sizeof(b), path);
    return static_cast<uint64_t>(b[0]) |
           (static_cast<uint64_t>(b[1]) << 8) |
           (static_cast<uint64_t>(b[2]) << 16) |
           (static_cast<uint64_t>(b[3]) << 24) |
           (static_cast<uint64_t>(b[4]) << 32) |
           (static_cast<uint64_t>(b[5]) << 40) |
           (static_cast<uint64_t>(b[6]) << 48) |
           (static_cast<uint64_t>(b[7]) << 56);
}

static int64_t read_i64_le(FILE* fp, const std::string& path) {
    return static_cast<int64_t>(read_u64_le(fp, path));
}

static std::string read_string(FILE* fp, const std::string& path) {
    uint32_t len = read_u32_le(fp, path);
    std::string value(len, '\0');
    if (len > 0) read_exact(fp, value.data(), len, path);
    return value;
}

static void read_magic(FILE* fp, const std::string& path, const char magic[8]) {
    char got[8];
    read_exact(fp, got, sizeof(got), path);
    if (std::memcmp(got, magic, sizeof(got)) != 0) {
        die("bad magic in %s", path.c_str());
    }
}

static bool has_more(FILE* fp, const std::string& path) {
    int first = std::fgetc(fp);
    if (first == EOF) return false;
    if (std::ungetc(first, fp) == EOF) die("ungetc failed for %s", path.c_str());
    return true;
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

    const char* required[] = {"format_version", "n_samples", "n_haps", "n_words", "n_ancestries"};
    for (const char* key : required) {
        if (!kv.count(key)) die("missing %s in %s", key, path.c_str());
    }

    Meta meta;
    meta.format_version = static_cast<uint32_t>(parse_u64(kv["format_version"], "format_version"));
    meta.n_samples = parse_u64(kv["n_samples"], "n_samples");
    meta.n_haps = parse_u64(kv["n_haps"], "n_haps");
    meta.n_words = parse_u64(kv["n_words"], "n_words");
    meta.n_ancestries = parse_u64(kv["n_ancestries"], "n_ancestries");

    if (meta.format_version != 1) die("unsupported FELIXla format_version: %u", meta.format_version);
    if (meta.n_samples == 0) die("n_samples must be positive");
    if (meta.n_haps != 2 * meta.n_samples) die("n_haps must equal 2 * n_samples");
    if (meta.n_words != (meta.n_haps + 63) / 64) die("n_words does not match n_haps");
    if (meta.n_ancestries == 0 || meta.n_ancestries > 32) die("n_ancestries must be in [1, 32]");
    if (meta.n_samples > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
        die("n_samples is too large for this platform");
    }

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

static VariantRecord read_variant_record(FILE* fp, const std::string& path, bool common) {
    VariantRecord record;
    record.common = common;
    (void)read_u64_le(fp, path);
    record.global_variant_index = read_u32_le(fp, path);
    record.chr = read_string(fp, path);
    record.pos = read_i64_le(fp, path);
    record.id = read_string(fp, path);
    record.ref = read_string(fp, path);
    record.alt = read_string(fp, path);
    record.alt_index = read_u32_le(fp, path);
    if (common) {
        record.block_id = read_u32_le(fp, path);
        record.payload_offset = read_u64_le(fp, path);
        record.mac = read_u32_le(fp, path);
    } else {
        record.payload_offset = read_u64_le(fp, path);
        record.n_carriers = read_u32_le(fp, path);
        record.mac = read_u32_le(fp, path);
    }
    return record;
}

static bool matches_query(const VariantRecord& record, const Args& args) {
    if (args.by_global_index) {
        return record.global_variant_index == args.global_index;
    }
    if (record.chr != args.chr || record.pos != args.pos) return false;
    if (args.has_ref_alt && (record.ref != args.ref || record.alt != args.alt)) return false;
    return true;
}

static void collect_matches(
    const std::string& prefix,
    bool common,
    const Args& args,
    std::vector<VariantRecord>& out
) {
    std::string path = prefix + (common ? ".common.variant.mks" : ".rare.variant.mks");
    FILE* fp = open_file_or_die(path, "rb");
    const char common_magic[8] = {'T', 'R', 'C', 'M', 'M', 'K', 'S', '1'};
    const char rare_magic[8] = {'T', 'R', 'R', 'A', 'M', 'K', 'S', '1'};
    read_magic(fp, path, common ? common_magic : rare_magic);

    while (has_more(fp, path)) {
        VariantRecord record = read_variant_record(fp, path, common);
        if (matches_query(record, args)) out.push_back(std::move(record));
    }

    std::fclose(fp);
}

static bool bit_at(const std::vector<uint64_t>& words, uint64_t hap_id) {
    uint64_t word_i = hap_id >> 6;
    uint64_t bit_i = hap_id & 63u;
    if (word_i >= words.size()) return false;
    return ((words[static_cast<size_t>(word_i)] >> bit_i) & 1u) != 0;
}

static std::vector<uint64_t> read_words(FILE* fp, const std::string& path, uint64_t offset, uint64_t n_words) {
    std::vector<uint64_t> words(static_cast<size_t>(n_words), 0);
    seek_or_die(fp, offset, path);
    for (uint64_t i = 0; i < n_words; ++i) {
        words[static_cast<size_t>(i)] = read_u64_le(fp, path);
    }
    return words;
}

static std::vector<std::vector<uint64_t>> load_ancestry_masks(
    const std::string& prefix,
    FILE* anc_idx_fp,
    FILE* anc_bin_fp,
    const Meta& meta,
    uint32_t block_id
) {
    std::string idx_path = prefix + ".ancblock.idx";
    std::string bin_path = prefix + ".ancblock.bin";
    uint64_t idx_offset = 8ULL + static_cast<uint64_t>(block_id) * 20ULL;
    seek_or_die(anc_idx_fp, idx_offset, idx_path);
    uint32_t got_block_id = read_u32_le(anc_idx_fp, idx_path);
    (void)read_u64_le(anc_idx_fp, idx_path);
    uint64_t anc_offset = read_u64_le(anc_idx_fp, idx_path);
    if (got_block_id != block_id) {
        die("ancestry block index mismatch: got %u, expected %u", got_block_id, block_id);
    }

    std::vector<std::vector<uint64_t>> masks(
        static_cast<size_t>(meta.n_ancestries),
        std::vector<uint64_t>(static_cast<size_t>(meta.n_words), 0)
    );
    seek_or_die(anc_bin_fp, anc_offset, bin_path);
    for (uint64_t ancestry = 0; ancestry < meta.n_ancestries; ++ancestry) {
        for (uint64_t word_i = 0; word_i < meta.n_words; ++word_i) {
            masks[static_cast<size_t>(ancestry)][static_cast<size_t>(word_i)] = read_u64_le(anc_bin_fp, bin_path);
        }
    }
    return masks;
}

static void print_header(const Meta& meta) {
    std::printf("global_variant_index\tchr\tpos\tid\tref\talt\tsample\tDSALL");
    for (uint64_t ancestry = 0; ancestry < meta.n_ancestries; ++ancestry) {
        std::printf("\tDS%llu", static_cast<unsigned long long>(ancestry + 1));
    }
    std::printf("\n");
}

static void print_row(
    const VariantRecord& record,
    const std::string& sample,
    uint32_t ds_all,
    const std::vector<uint32_t>& ds_by_ancestry
) {
    std::printf(
        "%u\t%s\t%lld\t%s\t%s\t%s\t%s\t%u",
        record.global_variant_index,
        record.chr.c_str(),
        static_cast<long long>(record.pos),
        record.id.c_str(),
        record.ref.c_str(),
        record.alt.c_str(),
        sample.c_str(),
        ds_all
    );
    for (uint32_t value : ds_by_ancestry) std::printf("\t%u", value);
    std::printf("\n");
}

static void query_common(
    const VariantRecord& record,
    const Meta& meta,
    const std::vector<std::string>& samples,
    FILE* common_fp,
    FILE* anc_idx_fp,
    FILE* anc_bin_fp,
    const std::string& prefix,
    bool nonzero_only
) {
    std::string common_path = prefix + ".common.geno.bin";
    std::vector<uint64_t> genotype = read_words(common_fp, common_path, record.payload_offset, meta.n_words);
    std::vector<std::vector<uint64_t>> masks =
        load_ancestry_masks(prefix, anc_idx_fp, anc_bin_fp, meta, record.block_id);

    std::vector<uint32_t> ds_by_ancestry(static_cast<size_t>(meta.n_ancestries), 0);
    for (uint64_t sample_i = 0; sample_i < meta.n_samples; ++sample_i) {
        uint64_t hap0 = 2 * sample_i;
        uint64_t hap1 = hap0 + 1;
        uint32_t ds_all = static_cast<uint32_t>(bit_at(genotype, hap0)) +
                          static_cast<uint32_t>(bit_at(genotype, hap1));
        if (nonzero_only && ds_all == 0) continue;

        std::fill(ds_by_ancestry.begin(), ds_by_ancestry.end(), 0);
        for (uint64_t ancestry = 0; ancestry < meta.n_ancestries; ++ancestry) {
            const std::vector<uint64_t>& mask = masks[static_cast<size_t>(ancestry)];
            uint32_t ds = 0;
            if (bit_at(genotype, hap0) && bit_at(mask, hap0)) ++ds;
            if (bit_at(genotype, hap1) && bit_at(mask, hap1)) ++ds;
            ds_by_ancestry[static_cast<size_t>(ancestry)] = ds;
        }

        print_row(record, samples[static_cast<size_t>(sample_i)], ds_all, ds_by_ancestry);
    }
}

static void query_rare(
    const VariantRecord& record,
    const Meta& meta,
    const std::vector<std::string>& samples,
    FILE* rare_fp,
    const std::string& prefix,
    bool nonzero_only
) {
    std::string rare_path = prefix + ".rare.carrier.bin";
    std::vector<uint32_t> ds_all(static_cast<size_t>(meta.n_samples), 0);
    std::vector<uint32_t> ds_by_ancestry(
        static_cast<size_t>(meta.n_samples * meta.n_ancestries),
        0
    );

    seek_or_die(rare_fp, record.payload_offset, rare_path);
    for (uint32_t i = 0; i < record.n_carriers; ++i) {
        uint32_t pos_index = read_u32_le(rare_fp, rare_path);
        uint32_t anc_hap = read_u32_le(rare_fp, rare_path);
        if (pos_index != record.global_variant_index) {
            die("rare carrier pos_index mismatch: got %u, expected %u", pos_index, record.global_variant_index);
        }

        uint32_t ancestry = anc_hap >> 27;
        uint32_t hap_id = anc_hap & ((1u << 27) - 1u);
        if (hap_id >= meta.n_haps) die("rare carrier hap_id %u exceeds n_haps", hap_id);
        if (ancestry >= meta.n_ancestries) die("rare carrier ancestry %u exceeds n_ancestries", ancestry);

        uint64_t sample_i = hap_id >> 1;
        ++ds_all[static_cast<size_t>(sample_i)];
        ++ds_by_ancestry[static_cast<size_t>(sample_i * meta.n_ancestries + ancestry)];
    }

    std::vector<uint32_t> row(static_cast<size_t>(meta.n_ancestries), 0);
    for (uint64_t sample_i = 0; sample_i < meta.n_samples; ++sample_i) {
        uint32_t total = ds_all[static_cast<size_t>(sample_i)];
        if (nonzero_only && total == 0) continue;

        for (uint64_t ancestry = 0; ancestry < meta.n_ancestries; ++ancestry) {
            row[static_cast<size_t>(ancestry)] =
                ds_by_ancestry[static_cast<size_t>(sample_i * meta.n_ancestries + ancestry)];
        }
        print_row(record, samples[static_cast<size_t>(sample_i)], total, row);
    }
}

int main(int argc, char** argv) {
    Args args = parse_args(argc, argv);
    Meta meta = read_meta(args.prefix + ".meta");
    std::vector<std::string> samples = read_samples(args.prefix + ".samples", meta.n_samples);

    std::vector<VariantRecord> matches;
    collect_matches(args.prefix, true, args, matches);
    collect_matches(args.prefix, false, args, matches);
    std::sort(matches.begin(), matches.end(), [](const VariantRecord& a, const VariantRecord& b) {
        return a.global_variant_index < b.global_variant_index;
    });

    if (matches.empty()) die("no FELIXla variants matched the query");

    FILE* common_fp = open_file_or_die(args.prefix + ".common.geno.bin", "rb");
    FILE* rare_fp = open_file_or_die(args.prefix + ".rare.carrier.bin", "rb");
    FILE* anc_idx_fp = open_file_or_die(args.prefix + ".ancblock.idx", "rb");
    FILE* anc_bin_fp = open_file_or_die(args.prefix + ".ancblock.bin", "rb");
    const char anc_idx_magic[8] = {'T', 'R', 'A', 'N', 'I', 'D', 'X', '2'};
    read_magic(anc_idx_fp, args.prefix + ".ancblock.idx", anc_idx_magic);

    print_header(meta);
    for (const VariantRecord& record : matches) {
        if (record.common) {
            query_common(record, meta, samples, common_fp, anc_idx_fp, anc_bin_fp, args.prefix, args.nonzero_only);
        } else {
            query_rare(record, meta, samples, rare_fp, args.prefix, args.nonzero_only);
        }
    }

    std::fclose(common_fp);
    std::fclose(rare_fp);
    std::fclose(anc_idx_fp);
    std::fclose(anc_bin_fp);
    return 0;
}
