// tractor_hybrid_extract_region.cpp
// designed by Kai, implemented by codex
//
// Extract a 1-based inclusive genomic region from an existing tractor_hybrid
// packed prefix into a new tractor_hybrid packed prefix.

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
    uint64_t n_samples = 0;
    uint64_t n_haps = 0;
    uint64_t n_words = 0;
    uint64_t n_ancestries = 0;
    std::vector<std::string> lines;
};

struct Region {
    std::string chr;
    int64_t start = 1;
    int64_t end = 0;
    std::string label;
};

struct VariantRecord {
    bool common = false;
    uint64_t local_index = 0;
    uint32_t global_variant_index = 0;
    std::string chr;
    int64_t pos = 0;
    std::string id;
    std::string ref;
    std::string alt;
    uint32_t alt_index = 0;
    uint32_t block_id = UINT32_MAX;
    uint64_t payload_offset = 0;
    uint32_t n_carriers = 0;
    uint32_t mac = 0;
};

struct AncBlockRecord {
    uint32_t old_block_id = 0;
    uint32_t new_block_id = 0;
    std::string chr;
    int64_t start = 0;
    int64_t end = 0;
    uint64_t anc_offset = 0;
};

static constexpr uint32_t kPackedHapMask = 0x07ffffffu;

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

static std::string strip_commas(std::string value) {
    value.erase(std::remove(value.begin(), value.end(), ','), value.end());
    return value;
}

static int64_t parse_i64_string(const std::string& value, const char* field) {
    std::string clean = strip_commas(value);
    char* end = nullptr;
    long long parsed = std::strtoll(clean.c_str(), &end, 10);
    if (!end || *end != '\0') {
        die("invalid integer for %s: %s", field, value.c_str());
    }
    return static_cast<int64_t>(parsed);
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

static uint64_t tell_or_die(FILE* fp, const char* path) {
    off_t offset = ftello(fp);
    if (offset < 0) die("ftello failed for %s", path);
    return static_cast<uint64_t>(offset);
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

static void write_exact(FILE* fp, const void* data, size_t bytes, const char* what) {
    if (bytes == 0) return;
    if (std::fwrite(data, 1, bytes, fp) != bytes) {
        die("failed writing %s", what);
    }
}

static uint32_t read_u32_le(FILE* fp, const char* path) {
    uint8_t b[4];
    read_exact(fp, b, sizeof(b), path);
    return static_cast<uint32_t>(b[0]) |
           (static_cast<uint32_t>(b[1]) << 8) |
           (static_cast<uint32_t>(b[2]) << 16) |
           (static_cast<uint32_t>(b[3]) << 24);
}

static uint64_t read_u64_le(FILE* fp, const char* path) {
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

static int64_t read_i64_le(FILE* fp, const char* path) {
    return static_cast<int64_t>(read_u64_le(fp, path));
}

static std::string read_string(FILE* fp, const char* path) {
    uint32_t len = read_u32_le(fp, path);
    std::string value(len, '\0');
    if (len > 0) read_exact(fp, value.data(), len, path);
    return value;
}

static void write_u32_le(FILE* fp, uint32_t value, const char* what) {
    uint8_t b[4] = {
        static_cast<uint8_t>(value),
        static_cast<uint8_t>(value >> 8),
        static_cast<uint8_t>(value >> 16),
        static_cast<uint8_t>(value >> 24)
    };
    write_exact(fp, b, sizeof(b), what);
}

static void write_u64_le(FILE* fp, uint64_t value, const char* what) {
    uint8_t b[8] = {
        static_cast<uint8_t>(value),
        static_cast<uint8_t>(value >> 8),
        static_cast<uint8_t>(value >> 16),
        static_cast<uint8_t>(value >> 24),
        static_cast<uint8_t>(value >> 32),
        static_cast<uint8_t>(value >> 40),
        static_cast<uint8_t>(value >> 48),
        static_cast<uint8_t>(value >> 56)
    };
    write_exact(fp, b, sizeof(b), what);
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

static void read_magic(FILE* fp, const std::string& path, const char magic[8]) {
    char got[8];
    read_exact(fp, got, sizeof(got), path.c_str());
    if (std::memcmp(got, magic, sizeof(got)) != 0) {
        die("bad magic in %s", path.c_str());
    }
}

static void write_magic(FILE* fp, const char magic[8]) {
    write_exact(fp, magic, 8, "magic header");
}

static bool try_record_start(FILE* fp, const std::string& path) {
    int first = std::fgetc(fp);
    if (first == EOF) return false;
    if (std::ungetc(first, fp) == EOF) die("ungetc failed for %s", path.c_str());
    return true;
}

static bool chrom_matches(const std::string& observed, const std::string& requested) {
    if (observed == requested) return true;

    if (observed.size() > 3 && observed.compare(0, 3, "chr") == 0) {
        return observed.substr(3) == requested;
    }

    if (requested.size() > 3 && requested.compare(0, 3, "chr") == 0) {
        return requested.substr(3) == observed;
    }

    return false;
}

static bool in_region(const std::string& chr, int64_t pos, const Region& region) {
    return chrom_matches(chr, region.chr) && pos >= region.start && pos <= region.end;
}

static bool overlaps_region(const AncBlockRecord& block, const Region& region) {
    return chrom_matches(block.chr, region.chr) &&
           block.end >= region.start &&
           block.start <= region.end;
}

static bool ancestry_blocks_cover_position(
    const std::vector<AncBlockRecord>& blocks,
    const std::string& chr,
    int64_t pos
) {
    for (const AncBlockRecord& block : blocks) {
        if (chrom_matches(block.chr, chr) && pos >= block.start && pos <= block.end) {
            return true;
        }
    }
    return false;
}

static Region parse_region_string(const std::string& region_text) {
    size_t colon = region_text.find(':');
    size_t dash = region_text.find('-', colon == std::string::npos ? 0 : colon + 1);
    if (colon == std::string::npos || dash == std::string::npos || dash <= colon + 1) {
        die("region must look like chr:start-end: %s", region_text.c_str());
    }

    Region region;
    region.chr = region_text.substr(0, colon);
    region.start = parse_i64_string(region_text.substr(colon + 1, dash - colon - 1), "region start");
    region.end = parse_i64_string(region_text.substr(dash + 1), "region end");
    if (region.chr.empty()) die("region chromosome is empty");
    if (region.start <= 0 || region.end <= 0 || region.start > region.end) {
        die("invalid region coordinates: %s", region_text.c_str());
    }
    region.label = region.chr + ":" + std::to_string(region.start) + "-" + std::to_string(region.end);
    return region;
}

static Region parse_region_args(int argc, char** argv, int& out_prefix_arg) {
    if (argc == 4) {
        out_prefix_arg = 3;
        return parse_region_string(argv[2]);
    }

    if (argc == 6) {
        Region region;
        region.chr = argv[2];
        region.start = parse_i64_string(argv[3], "region start");
        region.end = parse_i64_string(argv[4], "region end");
        if (region.chr.empty()) die("region chromosome is empty");
        if (region.start <= 0 || region.end <= 0 || region.start > region.end) {
            die("invalid region coordinates: %s %s %s", argv[2], argv[3], argv[4]);
        }
        region.label = region.chr + ":" + std::to_string(region.start) + "-" + std::to_string(region.end);
        out_prefix_arg = 5;
        return region;
    }

    die("invalid argument count");
}

static Meta read_meta(const std::string& path) {
    std::ifstream in(path);
    if (!in) die("cannot open %s", path.c_str());

    std::unordered_map<std::string, std::string> kv;
    Meta meta;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        meta.lines.push_back(line);
        std::vector<std::string> fields = split_tabs(line);
        if (fields.size() >= 2) kv[fields[0]] = fields[1];
    }

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

static void copy_text_file(const std::string& in_path, const std::string& out_path) {
    std::ifstream in(in_path, std::ios::binary);
    if (!in) die("cannot open %s", in_path.c_str());
    std::ofstream out(out_path, std::ios::binary);
    if (!out) die("cannot open %s", out_path.c_str());
    out << in.rdbuf();
    if (!out) die("failed writing %s", out_path.c_str());
}

static void write_meta(
    const Meta& meta,
    const std::string& out_path,
    const std::string& in_prefix,
    const Region& region
) {
    std::ofstream out(out_path);
    if (!out) die("cannot open %s", out_path.c_str());

    for (const std::string& line : meta.lines) {
        out << line << '\n';
    }
    out << "source_hybrid_prefix\t" << in_prefix << '\n';
    out << "extracted_region\t" << region.label << '\n';
    if (!out) die("failed writing %s", out_path.c_str());
}

static VariantRecord read_common_record(FILE* fp, const std::string& path) {
    VariantRecord record;
    record.common = true;
    record.local_index = read_u64_le(fp, path.c_str());
    record.global_variant_index = read_u32_le(fp, path.c_str());
    record.chr = read_string(fp, path.c_str());
    record.pos = read_i64_le(fp, path.c_str());
    record.id = read_string(fp, path.c_str());
    record.ref = read_string(fp, path.c_str());
    record.alt = read_string(fp, path.c_str());
    record.alt_index = read_u32_le(fp, path.c_str());
    record.block_id = read_u32_le(fp, path.c_str());
    record.payload_offset = read_u64_le(fp, path.c_str());
    record.mac = read_u32_le(fp, path.c_str());
    return record;
}

static VariantRecord read_rare_record(FILE* fp, const std::string& path) {
    VariantRecord record;
    record.common = false;
    record.local_index = read_u64_le(fp, path.c_str());
    record.global_variant_index = read_u32_le(fp, path.c_str());
    record.chr = read_string(fp, path.c_str());
    record.pos = read_i64_le(fp, path.c_str());
    record.id = read_string(fp, path.c_str());
    record.ref = read_string(fp, path.c_str());
    record.alt = read_string(fp, path.c_str());
    record.alt_index = read_u32_le(fp, path.c_str());
    record.payload_offset = read_u64_le(fp, path.c_str());
    record.n_carriers = read_u32_le(fp, path.c_str());
    record.mac = read_u32_le(fp, path.c_str());
    return record;
}

static AncBlockRecord read_anc_record(FILE* fp, const std::string& path) {
    AncBlockRecord record;
    record.old_block_id = read_u32_le(fp, path.c_str());
    record.chr = read_string(fp, path.c_str());
    record.start = read_i64_le(fp, path.c_str());
    record.end = read_i64_le(fp, path.c_str());
    record.anc_offset = read_u64_le(fp, path.c_str());
    return record;
}

static void write_common_record(
    FILE* fp,
    uint64_t local_index,
    uint32_t global_variant_index,
    const VariantRecord& record,
    uint32_t new_block_id,
    uint64_t geno_offset
) {
    write_u64_le(fp, local_index, "common mks");
    write_u32_le(fp, global_variant_index, "common mks");
    write_string(fp, record.chr, "common mks chr");
    write_i64_le(fp, record.pos, "common mks pos");
    write_string(fp, record.id, "common mks id");
    write_string(fp, record.ref, "common mks ref");
    write_string(fp, record.alt, "common mks alt");
    write_u32_le(fp, record.alt_index, "common mks");
    write_u32_le(fp, new_block_id, "common mks");
    write_u64_le(fp, geno_offset, "common mks");
    write_u32_le(fp, record.mac, "common mks");
}

static void write_rare_record(
    FILE* fp,
    uint64_t local_index,
    uint32_t global_variant_index,
    const VariantRecord& record,
    uint64_t carrier_offset
) {
    write_u64_le(fp, local_index, "rare mks");
    write_u32_le(fp, global_variant_index, "rare mks");
    write_string(fp, record.chr, "rare mks chr");
    write_i64_le(fp, record.pos, "rare mks pos");
    write_string(fp, record.id, "rare mks id");
    write_string(fp, record.ref, "rare mks ref");
    write_string(fp, record.alt, "rare mks alt");
    write_u32_le(fp, record.alt_index, "rare mks");
    write_u64_le(fp, carrier_offset, "rare mks");
    write_u32_le(fp, record.n_carriers, "rare mks");
    write_u32_le(fp, record.mac, "rare mks");
}

static void write_anc_record(FILE* fp, const AncBlockRecord& record, uint64_t anc_offset) {
    write_u32_le(fp, record.new_block_id, "ancestry mks");
    write_string(fp, record.chr, "ancestry mks chr");
    write_i64_le(fp, record.start, "ancestry mks start");
    write_i64_le(fp, record.end, "ancestry mks end");
    write_u64_le(fp, anc_offset, "ancestry mks offset");
}

static void copy_fixed_bytes(
    FILE* in_fp,
    const std::string& in_path,
    uint64_t offset,
    FILE* out_fp,
    const char* out_what,
    uint64_t bytes
) {
    seek_or_die(in_fp, offset, in_path.c_str());
    std::vector<uint8_t> buffer(static_cast<size_t>(bytes));
    read_exact(in_fp, buffer.data(), buffer.size(), in_path.c_str());
    write_exact(out_fp, buffer.data(), buffer.size(), out_what);
}

static void read_selected_anc_blocks(
    const std::string& in_prefix,
    const Meta& meta,
    const Region& region,
    std::vector<AncBlockRecord>& selected_blocks,
    std::unordered_map<uint32_t, uint32_t>& block_id_map
) {
    const std::string anc_mks_path = in_prefix + ".ancblock.mks";
    FILE* fp = open_file_or_die(anc_mks_path, "rb");
    const char anc_mks_magic[8] = {'T', 'R', 'A', 'N', 'M', 'K', 'S', '1'};
    read_magic(fp, anc_mks_path, anc_mks_magic);

    while (try_record_start(fp, anc_mks_path)) {
        AncBlockRecord block = read_anc_record(fp, anc_mks_path);
        if (!overlaps_region(block, region)) continue;

        block.start = std::max(block.start, region.start);
        block.end = std::min(block.end, region.end);
        block.new_block_id = static_cast<uint32_t>(selected_blocks.size());

        if (block_id_map.count(block.old_block_id)) {
            die("duplicate ancestry block id in %s: %u", anc_mks_path.c_str(), block.old_block_id);
        }
        if (selected_blocks.size() > std::numeric_limits<uint32_t>::max()) {
            die("too many selected ancestry blocks");
        }
        block_id_map[block.old_block_id] = block.new_block_id;
        selected_blocks.push_back(block);
    }

    std::fclose(fp);

    if (meta.n_ancestries > 0 && meta.n_words > 0 && selected_blocks.empty()) {
        std::fprintf(
            stderr,
            "WARNING: no ancestry blocks overlap %s; output variants may be empty or unreadable if selected variants need ancestry blocks.\n",
            region.label.c_str()
        );
    }
}

static void write_selected_anc_blocks(
    const std::string& in_prefix,
    const std::string& out_prefix,
    const Meta& meta,
    const std::vector<AncBlockRecord>& selected_blocks
) {
    const std::string in_anc_bin_path = in_prefix + ".ancblock.bin";
    const std::string out_anc_bin_path = out_prefix + ".ancblock.bin";
    const std::string out_anc_mks_path = out_prefix + ".ancblock.mks";
    const std::string out_anc_idx_path = out_prefix + ".ancblock.idx";

    FILE* in_bin = open_file_or_die(in_anc_bin_path, "rb");
    FILE* out_bin = open_file_or_die(out_anc_bin_path, "wb");
    FILE* out_mks = open_file_or_die(out_anc_mks_path, "wb");
    FILE* out_idx = open_file_or_die(out_anc_idx_path, "wb");

    const char anc_mks_magic[8] = {'T', 'R', 'A', 'N', 'M', 'K', 'S', '1'};
    const char anc_idx_magic[8] = {'T', 'R', 'A', 'N', 'I', 'D', 'X', '2'};
    write_magic(out_mks, anc_mks_magic);
    write_magic(out_idx, anc_idx_magic);

    uint64_t block_bytes = meta.n_ancestries * meta.n_words * sizeof(uint64_t);
    for (const AncBlockRecord& block : selected_blocks) {
        uint64_t anc_offset = tell_or_die(out_bin, out_anc_bin_path.c_str());
        copy_fixed_bytes(in_bin, in_anc_bin_path, block.anc_offset, out_bin, "ancestry block", block_bytes);

        uint64_t mks_offset = tell_or_die(out_mks, out_anc_mks_path.c_str());
        write_anc_record(out_mks, block, anc_offset);

        write_u32_le(out_idx, block.new_block_id, "ancestry idx");
        write_u64_le(out_idx, mks_offset, "ancestry idx");
        write_u64_le(out_idx, anc_offset, "ancestry idx");
    }

    std::fclose(in_bin);
    std::fclose(out_bin);
    std::fclose(out_mks);
    std::fclose(out_idx);
}

static bool next_variant(
    FILE* common_mks_fp,
    const std::string& common_mks_path,
    FILE* rare_mks_fp,
    const std::string& rare_mks_path,
    bool& have_common,
    bool& have_rare,
    VariantRecord& common_record,
    VariantRecord& rare_record,
    VariantRecord& out
) {
    if (!have_common && try_record_start(common_mks_fp, common_mks_path)) {
        common_record = read_common_record(common_mks_fp, common_mks_path);
        have_common = true;
    }

    if (!have_rare && try_record_start(rare_mks_fp, rare_mks_path)) {
        rare_record = read_rare_record(rare_mks_fp, rare_mks_path);
        have_rare = true;
    }

    if (!have_common && !have_rare) return false;

    bool take_common = have_common && !have_rare;
    if (have_common && have_rare) {
        if (common_record.global_variant_index == rare_record.global_variant_index) {
            die("global variant %u exists in both common and rare streams", common_record.global_variant_index);
        }
        take_common = common_record.global_variant_index < rare_record.global_variant_index;
    }

    if (take_common) {
        out = common_record;
        have_common = false;
    } else {
        out = rare_record;
        have_rare = false;
    }
    return true;
}

static void write_selected_variants(
    const std::string& in_prefix,
    const std::string& out_prefix,
    const Meta& meta,
    const Region& region,
    const std::unordered_map<uint32_t, uint32_t>& block_id_map,
    const std::vector<AncBlockRecord>& selected_blocks,
    uint64_t& common_written,
    uint64_t& rare_written,
    uint64_t& total_written
) {
    const std::string in_common_mks_path = in_prefix + ".common.variant.mks";
    const std::string in_rare_mks_path = in_prefix + ".rare.variant.mks";
    const std::string in_common_bin_path = in_prefix + ".common.geno.bin";
    const std::string in_rare_bin_path = in_prefix + ".rare.carrier.bin";
    const std::string out_common_mks_path = out_prefix + ".common.variant.mks";
    const std::string out_rare_mks_path = out_prefix + ".rare.variant.mks";
    const std::string out_common_idx_path = out_prefix + ".common.variant.idx";
    const std::string out_rare_idx_path = out_prefix + ".rare.variant.idx";
    const std::string out_common_bin_path = out_prefix + ".common.geno.bin";
    const std::string out_rare_bin_path = out_prefix + ".rare.carrier.bin";

    FILE* in_common_mks = open_file_or_die(in_common_mks_path, "rb");
    FILE* in_rare_mks = open_file_or_die(in_rare_mks_path, "rb");
    FILE* in_common_bin = open_file_or_die(in_common_bin_path, "rb");
    FILE* in_rare_bin = open_file_or_die(in_rare_bin_path, "rb");
    FILE* out_common_mks = open_file_or_die(out_common_mks_path, "wb");
    FILE* out_rare_mks = open_file_or_die(out_rare_mks_path, "wb");
    FILE* out_common_idx = open_file_or_die(out_common_idx_path, "wb");
    FILE* out_rare_idx = open_file_or_die(out_rare_idx_path, "wb");
    FILE* out_common_bin = open_file_or_die(out_common_bin_path, "wb");
    FILE* out_rare_bin = open_file_or_die(out_rare_bin_path, "wb");

    const char common_mks_magic[8] = {'T', 'R', 'C', 'M', 'M', 'K', 'S', '1'};
    const char rare_mks_magic[8] = {'T', 'R', 'R', 'A', 'M', 'K', 'S', '1'};
    const char common_idx_magic[8] = {'T', 'R', 'C', 'M', 'I', 'D', 'X', '2'};
    const char rare_idx_magic[8] = {'T', 'R', 'R', 'A', 'I', 'D', 'X', '2'};

    read_magic(in_common_mks, in_common_mks_path, common_mks_magic);
    read_magic(in_rare_mks, in_rare_mks_path, rare_mks_magic);
    write_magic(out_common_mks, common_mks_magic);
    write_magic(out_rare_mks, rare_mks_magic);
    write_magic(out_common_idx, common_idx_magic);
    write_magic(out_rare_idx, rare_idx_magic);

    bool have_common = false;
    bool have_rare = false;
    VariantRecord common_record;
    VariantRecord rare_record;
    VariantRecord record;

    uint64_t common_bytes = meta.n_words * sizeof(uint64_t);
    std::vector<uint8_t> common_payload(static_cast<size_t>(common_bytes));

    while (next_variant(
        in_common_mks,
        in_common_mks_path,
        in_rare_mks,
        in_rare_mks_path,
        have_common,
        have_rare,
        common_record,
        rare_record,
        record
    )) {
        if (!in_region(record.chr, record.pos, region)) continue;

        if (total_written > std::numeric_limits<uint32_t>::max()) {
            die("selected variant count exceeds uint32_t limit");
        }
        uint32_t new_global_index = static_cast<uint32_t>(total_written);

        if (record.common) {
            auto block_it = block_id_map.find(record.block_id);
            if (block_it == block_id_map.end()) {
                die(
                    "common variant %s:%lld references ancestry block %u outside selected region",
                    record.chr.c_str(),
                    static_cast<long long>(record.pos),
                    record.block_id
                );
            }

            seek_or_die(in_common_bin, record.payload_offset, in_common_bin_path.c_str());
            read_exact(in_common_bin, common_payload.data(), common_payload.size(), in_common_bin_path.c_str());

            uint64_t geno_offset = tell_or_die(out_common_bin, out_common_bin_path.c_str());
            write_exact(out_common_bin, common_payload.data(), common_payload.size(), "common genotype payload");

            uint64_t mks_offset = tell_or_die(out_common_mks, out_common_mks_path.c_str());
            write_common_record(
                out_common_mks,
                common_written,
                new_global_index,
                record,
                block_it->second,
                geno_offset
            );
            write_u64_le(out_common_idx, common_written, "common idx");
            write_u32_le(out_common_idx, new_global_index, "common idx");
            write_u64_le(out_common_idx, mks_offset, "common idx");
            write_u64_le(out_common_idx, geno_offset, "common idx");

            ++common_written;
        } else {
            if (!ancestry_blocks_cover_position(selected_blocks, record.chr, record.pos)) {
                die(
                    "rare variant %s:%lld is not covered by any selected ancestry block",
                    record.chr.c_str(),
                    static_cast<long long>(record.pos)
                );
            }

            uint64_t carrier_offset = tell_or_die(out_rare_bin, out_rare_bin_path.c_str());
            seek_or_die(in_rare_bin, record.payload_offset, in_rare_bin_path.c_str());

            for (uint32_t i = 0; i < record.n_carriers; ++i) {
                uint32_t old_pos_index = read_u32_le(in_rare_bin, in_rare_bin_path.c_str());
                uint32_t anc_hap = read_u32_le(in_rare_bin, in_rare_bin_path.c_str());
                if (old_pos_index != record.global_variant_index) {
                    die(
                        "rare carrier pos_index mismatch: carrier has %u, marker has %u",
                        old_pos_index,
                        record.global_variant_index
                    );
                }
                uint32_t hap_id = anc_hap & kPackedHapMask;
                if (hap_id >= meta.n_haps) {
                    die("rare carrier hap_id %u exceeds n_haps", hap_id);
                }

                write_u32_le(out_rare_bin, new_global_index, "rare carrier");
                write_u32_le(out_rare_bin, anc_hap, "rare carrier");
            }

            uint64_t mks_offset = tell_or_die(out_rare_mks, out_rare_mks_path.c_str());
            write_rare_record(out_rare_mks, rare_written, new_global_index, record, carrier_offset);
            write_u64_le(out_rare_idx, rare_written, "rare idx");
            write_u32_le(out_rare_idx, new_global_index, "rare idx");
            write_u64_le(out_rare_idx, mks_offset, "rare idx");
            write_u64_le(out_rare_idx, carrier_offset, "rare idx");
            write_u32_le(out_rare_idx, record.n_carriers, "rare idx");

            ++rare_written;
        }

        ++total_written;
    }

    std::fclose(in_common_mks);
    std::fclose(in_rare_mks);
    std::fclose(in_common_bin);
    std::fclose(in_rare_bin);
    std::fclose(out_common_mks);
    std::fclose(out_rare_mks);
    std::fclose(out_common_idx);
    std::fclose(out_rare_idx);
    std::fclose(out_common_bin);
    std::fclose(out_rare_bin);
}

static void print_usage(const char* prog) {
    std::fprintf(
        stderr,
        "Usage:\n"
        "  %s in_prefix chr:start-end out_prefix\n"
        "  %s in_prefix chr start end out_prefix\n\n"
        "Coordinates are 1-based and inclusive.\n\n"
        "Example:\n"
        "  %s chr22 chr22:16000000-17000000 chr22.region\n",
        prog,
        prog,
        prog
    );
}

int main(int argc, char** argv) {
    if (argc != 4 && argc != 6) {
        print_usage(argv[0]);
        return 1;
    }

    std::string in_prefix = argv[1];
    int out_prefix_arg = 0;
    Region region = parse_region_args(argc, argv, out_prefix_arg);
    std::string out_prefix = argv[out_prefix_arg];

    Meta meta = read_meta(in_prefix + ".meta");
    copy_text_file(in_prefix + ".samples", out_prefix + ".samples");
    write_meta(meta, out_prefix + ".meta", in_prefix, region);

    std::vector<AncBlockRecord> selected_blocks;
    std::unordered_map<uint32_t, uint32_t> block_id_map;
    read_selected_anc_blocks(in_prefix, meta, region, selected_blocks, block_id_map);
    write_selected_anc_blocks(in_prefix, out_prefix, meta, selected_blocks);

    uint64_t common_written = 0;
    uint64_t rare_written = 0;
    uint64_t total_written = 0;
    write_selected_variants(
        in_prefix,
        out_prefix,
        meta,
        region,
        block_id_map,
        selected_blocks,
        common_written,
        rare_written,
        total_written
    );

    std::fprintf(stderr, "Finished.\n");
    std::fprintf(stderr, "Region:                %s\n", region.label.c_str());
    std::fprintf(stderr, "Global variants:       %llu\n", static_cast<unsigned long long>(total_written));
    std::fprintf(stderr, "Common variants:       %llu\n", static_cast<unsigned long long>(common_written));
    std::fprintf(stderr, "Rare variants:         %llu\n", static_cast<unsigned long long>(rare_written));
    std::fprintf(stderr, "Ancestry blocks:       %llu\n", static_cast<unsigned long long>(selected_blocks.size()));

    return 0;
}
