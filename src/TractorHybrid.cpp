// TractorHybrid.cpp
// designed by Kai, implemented by codex
//
// FELIX step2 reader for the ancestry-aware packed hybrid format.

// [[Rcpp::depends(RcppArmadillo)]]
#include "TractorHybrid.hpp"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#ifdef __linux__
#include <fcntl.h>
#include <unistd.h>
#endif

namespace TRACTOR_HYBRID {

struct TractorHybridClass::Meta {
    uint64_t n_samples = 0;
    uint64_t n_haps = 0;
    uint64_t n_words = 0;
    uint64_t n_ancestries = 0;
    uint32_t format_version = 0;
};

struct TractorHybridClass::VariantRecord {
    bool valid = false;
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
    uint32_t alt_count = 0;
};

static constexpr uint32_t kPackedHapMask = 0x07ffffffu;
static constexpr uint32_t kMaxStringLen = 64u * 1024u * 1024u;

class ScopedSeconds {
public:
    explicit ScopedSeconds(double& target)
        : target_(target), start_(std::chrono::steady_clock::now()) {}

    ~ScopedSeconds() {
        const auto end = std::chrono::steady_clock::now();
        target_ += std::chrono::duration<double>(end - start_).count();
    }

private:
    double& target_;
    std::chrono::steady_clock::time_point start_;
};

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
        Rcpp::stop("Invalid unsigned integer in hybrid meta field %s: %s", field, value.c_str());
    }
    return static_cast<uint64_t>(parsed);
}

static void set_buffer(FILE* fp, const std::string& path, size_t bytes) {
    if (std::setvbuf(fp, nullptr, _IOFBF, bytes) != 0) {
        Rcpp::warning("setvbuf failed for %s; continuing with libc default buffering", path.c_str());
    }
}

static void advise_sequential(FILE* fp, const std::string& path) {
#ifdef __linux__
    int fd = fileno(fp);
    if (fd >= 0) {
        int rc = posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
        if (rc != 0) {
            Rcpp::warning("posix_fadvise(SEQUENTIAL) failed for %s; continuing", path.c_str());
        }
    }
#else
    (void)fp;
    (void)path;
#endif
}

static FILE* open_or_stop(const std::string& path, const char* mode, size_t buffer_size) {
    FILE* fp = std::fopen(path.c_str(), mode);
    if (!fp) {
        Rcpp::stop("Cannot open %s: %s", path.c_str(), std::strerror(errno));
    }
    set_buffer(fp, path, buffer_size);
    advise_sequential(fp, path);
    return fp;
}

static void close_if_open(FILE*& fp) {
    if (fp) {
        std::fclose(fp);
        fp = nullptr;
    }
}

static void seek_or_stop(FILE* fp, uint64_t offset, const std::string& path) {
    if (offset > static_cast<uint64_t>(std::numeric_limits<off_t>::max())) {
        Rcpp::stop("Offset is too large for this platform while reading %s", path.c_str());
    }
    if (fseeko(fp, static_cast<off_t>(offset), SEEK_SET) != 0) {
        Rcpp::stop("Seek failed in %s: %s", path.c_str(), std::strerror(errno));
    }
}

static bool try_read_exact(FILE* fp, void* data, size_t bytes) {
    if (bytes == 0) return true;
    return std::fread(data, 1, bytes, fp) == bytes;
}

static void read_exact(FILE* fp, void* data, size_t bytes, const std::string& path) {
    if (!try_read_exact(fp, data, bytes)) {
        Rcpp::stop("Unexpected EOF while reading %s", path.c_str());
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
    if (len > kMaxStringLen) {
        Rcpp::stop("String length %u is implausibly large while reading %s", len, path.c_str());
    }
    std::string value(len, '\0');
    if (len > 0) read_exact(fp, &value[0], len, path);
    return value;
}

static void read_magic(FILE* fp, const std::string& path, const char expected[8]) {
    char got[8];
    read_exact(fp, got, sizeof(got), path);
    if (std::memcmp(got, expected, sizeof(got)) != 0) {
        Rcpp::stop("Bad magic header in %s", path.c_str());
    }
}

static bool bit_at(const std::vector<uint64_t>& words, uint64_t hap_id) {
    return (words[hap_id >> 6] >> (hap_id & 63)) & 1ULL;
}

static bool parse_one_based_suffix(
    const std::string& field,
    const char* prefix,
    uint32_t max_value,
    uint32_t& zero_based
) {
    size_t prefix_len = std::strlen(prefix);
    if (field.size() <= prefix_len || field.compare(0, prefix_len, prefix) != 0) {
        return false;
    }

    char* end = nullptr;
    unsigned long parsed = std::strtoul(field.c_str() + prefix_len, &end, 10);
    if (!end || *end != '\0' || parsed == 0 || parsed > max_value) {
        Rcpp::stop(
            "Unsupported hybrid admixed FORMAT field %s; expected %s1..%s%u",
            field.c_str(),
            prefix,
            prefix,
            max_value
        );
    }

    zero_based = static_cast<uint32_t>(parsed - 1);
    return true;
}

static TractorHybridClass::Meta read_meta_file(const std::string& path) {
    std::ifstream in(path);
    if (!in) Rcpp::stop("Cannot open %s", path.c_str());

    std::unordered_map<std::string, std::string> kv;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> fields = split_tabs(line);
        if (fields.size() >= 2) kv[fields[0]] = fields[1];
    }

    const char* required[] = {"format_version", "n_samples", "n_haps", "n_words", "n_ancestries"};
    for (const char* key : required) {
        if (!kv.count(key)) Rcpp::stop("Missing %s in %s", key, path.c_str());
    }

    TractorHybridClass::Meta meta;
    meta.format_version = static_cast<uint32_t>(parse_u64(kv["format_version"], "format_version"));
    meta.n_samples = parse_u64(kv["n_samples"], "n_samples");
    meta.n_haps = parse_u64(kv["n_haps"], "n_haps");
    meta.n_words = parse_u64(kv["n_words"], "n_words");
    meta.n_ancestries = parse_u64(kv["n_ancestries"], "n_ancestries");

    if (meta.format_version != 1) Rcpp::stop("Unsupported hybrid format_version: %u", meta.format_version);
    if (meta.n_samples == 0) Rcpp::stop("Hybrid n_samples must be positive");
    if (meta.n_haps != 2 * meta.n_samples) Rcpp::stop("Hybrid n_haps must equal 2 * n_samples");
    if (meta.n_words != (meta.n_haps + 63) / 64) Rcpp::stop("Hybrid n_words does not match n_haps");
    if (meta.n_ancestries == 0 || meta.n_ancestries > 32) {
        Rcpp::stop("Hybrid n_ancestries must be in [1, 32]");
    }
    if (meta.n_samples > std::numeric_limits<uint32_t>::max()) {
        Rcpp::stop("Hybrid n_samples exceeds uint32_t capacity expected by SAIGE");
    }

    return meta;
}

static std::vector<std::string> read_samples_file(const std::string& path, uint64_t expected_n) {
    std::ifstream in(path);
    if (!in) Rcpp::stop("Cannot open %s", path.c_str());

    std::vector<std::string> samples;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        samples.push_back(line);
    }

    if (samples.size() != expected_n) {
        Rcpp::stop(
            "Sample count mismatch in %s: got %llu, expected %llu",
            path.c_str(),
            static_cast<unsigned long long>(samples.size()),
            static_cast<unsigned long long>(expected_n)
        );
    }
    return samples;
}

class TractorHybridClass::MksStream {
public:
    MksStream(const std::string& path, bool common) : path_(path), common_(common) {
        fp_ = open_or_stop(path_, "rb", 8u * 1024u * 1024u);
        const char common_magic[8] = {'T', 'R', 'C', 'M', 'M', 'K', 'S', '1'};
        const char rare_magic[8] = {'T', 'R', 'R', 'A', 'M', 'K', 'S', '1'};
        read_magic(fp_, path_, common_ ? common_magic : rare_magic);
        read_next();
    }

    ~MksStream() {
        close_if_open(fp_);
    }

    const VariantRecord& current() const {
        return current_;
    }

    uint64_t count_read() const {
        return count_read_;
    }

    void advance() {
        read_next();
    }

private:
    void read_next() {
        current_ = VariantRecord{};

        int first = std::fgetc(fp_);
        if (first == EOF) return;
        if (std::ungetc(first, fp_) == EOF) {
            Rcpp::stop("ungetc failed while reading %s", path_.c_str());
        }

        current_.valid = true;
        current_.common = common_;
        current_.local_index = read_u64_le(fp_, path_);
        current_.global_variant_index = read_u32_le(fp_, path_);
        current_.chr = read_string(fp_, path_);
        current_.pos = read_i64_le(fp_, path_);
        current_.id = read_string(fp_, path_);
        current_.ref = read_string(fp_, path_);
        current_.alt = read_string(fp_, path_);
        current_.alt_index = read_u32_le(fp_, path_);

        if (common_) {
            current_.block_id = read_u32_le(fp_, path_);
            current_.payload_offset = read_u64_le(fp_, path_);
            current_.alt_count = read_u32_le(fp_, path_);
        } else {
            current_.payload_offset = read_u64_le(fp_, path_);
            current_.n_carriers = read_u32_le(fp_, path_);
            current_.alt_count = read_u32_le(fp_, path_);
        }

        ++count_read_;
    }

    std::string path_;
    bool common_ = false;
    FILE* fp_ = nullptr;
    VariantRecord current_;
    uint64_t count_read_ = 0;
};

TractorHybridClass::TractorHybridClass(
    std::string prefix,
    std::vector<std::string> sampleInModel
) : prefix_(std::move(prefix)) {
    meta_ = new Meta(read_meta_file(prefix_ + ".meta"));
    samples_ = read_samples_file(prefix_ + ".samples", meta_->n_samples);

    if (sampleInModel.empty()) {
        sampleInModel = samples_;
    }

    std::unordered_map<std::string, int32_t> raw_index_by_sample;
    raw_index_by_sample.reserve(samples_.size() * 2);
    for (size_t i = 0; i < samples_.size(); ++i) {
        if (!raw_index_by_sample.emplace(samples_[i], static_cast<int32_t>(i)).second) {
            Rcpp::stop("Duplicate sample ID in hybrid sample file: %s", samples_[i].c_str());
        }
    }

    model_to_raw_sample_.reserve(sampleInModel.size());
    raw_sample_to_model_.assign(samples_.size(), -1);
    std::unordered_map<std::string, bool> seen_model_samples;
    seen_model_samples.reserve(sampleInModel.size() * 2);

    for (size_t i = 0; i < sampleInModel.size(); ++i) {
        const std::string& sample = sampleInModel[i];
        if (!seen_model_samples.emplace(sample, true).second) {
            Rcpp::stop("Duplicate sample ID in sampleInModel: %s", sample.c_str());
        }
        auto it = raw_index_by_sample.find(sample);
        if (it == raw_index_by_sample.end()) {
            Rcpp::stop("Sample in model is absent from hybrid sample file: %s", sample.c_str());
        }
        model_to_raw_sample_.push_back(it->second);
        raw_sample_to_model_[it->second] = static_cast<int32_t>(i);
    }

    common_stream_ = new MksStream(prefix_ + ".common.variant.mks", true);
    rare_stream_ = new MksStream(prefix_ + ".rare.variant.mks", false);

    common_bin_ = open_or_stop(prefix_ + ".common.geno.bin", "rb", 64u * 1024u * 1024u);
    rare_bin_ = open_or_stop(prefix_ + ".rare.carrier.bin", "rb", 64u * 1024u * 1024u);
    anc_bin_ = open_or_stop(prefix_ + ".ancblock.bin", "rb", 64u * 1024u * 1024u);
    anc_idx_ = open_or_stop(prefix_ + ".ancblock.idx", "rb", 8u * 1024u * 1024u);

    const char anc_idx_magic[8] = {'T', 'R', 'A', 'N', 'I', 'D', 'X', '2'};
    read_magic(anc_idx_, prefix_ + ".ancblock.idx", anc_idx_magic);

    common_words_.reserve(static_cast<size_t>(meta_->n_words));
    read_anc_blocks();
}

TractorHybridClass::~TractorHybridClass() {
    closegenofile();
    delete common_stream_;
    delete rare_stream_;
    delete meta_;
    common_stream_ = nullptr;
    rare_stream_ = nullptr;
    meta_ = nullptr;
}

void TractorHybridClass::read_anc_blocks() {
    const std::string path = prefix_ + ".ancblock.mks";
    FILE* fp = open_or_stop(path, "rb", 8u * 1024u * 1024u);
    const char anc_mks_magic[8] = {'T', 'R', 'A', 'N', 'M', 'K', 'S', '1'};
    read_magic(fp, path, anc_mks_magic);

    anc_blocks_.clear();
    while (true) {
        int first = std::fgetc(fp);
        if (first == EOF) break;
        if (std::ungetc(first, fp) == EOF) {
            close_if_open(fp);
            Rcpp::stop("ungetc failed while reading %s", path.c_str());
        }

        AncBlockRecord block;
        block.block_id = read_u32_le(fp, path);
        block.chr = read_string(fp, path);
        block.start_pos = read_i64_le(fp, path);
        block.end_pos = read_i64_le(fp, path);
        block.anc_offset = read_u64_le(fp, path);

        if (block.end_pos < block.start_pos) {
            close_if_open(fp);
            Rcpp::stop(
                "Invalid ancestry interval in %s for block %u: %lld > %lld",
                path.c_str(),
                block.block_id,
                static_cast<long long>(block.start_pos),
                static_cast<long long>(block.end_pos)
            );
        }
        anc_blocks_.push_back(std::move(block));
    }

    close_if_open(fp);
    if (anc_blocks_.empty()) {
        Rcpp::stop("Hybrid ancestry block stream is empty: %s", path.c_str());
    }
}

void TractorHybridClass::closegenofile() {
    close_if_open(common_bin_);
    close_if_open(rare_bin_);
    close_if_open(anc_bin_);
    close_if_open(anc_idx_);
}

uint32_t TractorHybridClass::getN0() const {
    return static_cast<uint32_t>(meta_->n_samples);
}

uint32_t TractorHybridClass::getN() const {
    return static_cast<uint32_t>(model_to_raw_sample_.size());
}

uint32_t TractorHybridClass::getM0() const {
    return getM();
}

uint32_t TractorHybridClass::getM() const {
    uint64_t n = n_common_ + n_rare_;
    if (n > std::numeric_limits<uint32_t>::max()) {
        Rcpp::stop("Hybrid marker count exceeds uint32_t capacity");
    }
    return static_cast<uint32_t>(n);
}

const TractorHybridClass::VariantRecord* TractorHybridClass::current_record() const {
    const VariantRecord& c = common_stream_->current();
    const VariantRecord& r = rare_stream_->current();

    if (!c.valid && !r.valid) return nullptr;
    if (c.valid && !r.valid) return &c;
    if (!c.valid && r.valid) return &r;
    if (c.global_variant_index == r.global_variant_index) {
        Rcpp::stop("Hybrid global variant %u exists in both common and rare streams", c.global_variant_index);
    }
    return c.global_variant_index < r.global_variant_index ? &c : &r;
}

void TractorHybridClass::advance_current() {
    const VariantRecord* record = current_record();
    if (!record) return;

    if (record->common) {
        common_stream_->advance();
        ++n_common_;
    } else {
        rare_stream_->advance();
        ++n_rare_;
    }
}

bool TractorHybridClass::advance_to(uint32_t global_variant_index) {
    while (true) {
        const VariantRecord* record = current_record();
        if (!record) return false;
        if (record->global_variant_index == global_variant_index) return true;
        if (record->global_variant_index > global_variant_index) {
            Rcpp::stop(
                "Hybrid reader is past requested global_variant_index %u; current is %u. "
                "Step2 access must be in increasing marker order.",
                global_variant_index,
                record->global_variant_index
            );
        }
        advance_current();
    }
}

void TractorHybridClass::advance_to_next_in_filter() {
    if (!streaming_mode_ || !has_region_filter_) return;

    while (true) {
        const VariantRecord* record = current_record();
        if (!record) return;
        if (record->chr == region_chr_ &&
            record->pos >= region_start_ &&
            record->pos <= region_end_) {
            return;
        }
        if (record->chr == region_chr_ && record->pos > region_end_) {
            return;
        }
        advance_current();
    }
}

const TractorHybridClass::VariantRecord* TractorHybridClass::record_for_request(
    uint64_t t_gIndex,
    bool& t_isBoolRead
) {
    if (streaming_mode_) {
        advance_to_next_in_filter();
        t_isBoolRead = !check_iterator_end();
        if (!t_isBoolRead) return nullptr;
        return current_record();
    }

    if (t_gIndex > std::numeric_limits<uint32_t>::max()) {
        Rcpp::stop("Hybrid global variant index %llu exceeds uint32_t capacity",
                   static_cast<unsigned long long>(t_gIndex));
    }
    uint32_t target = static_cast<uint32_t>(t_gIndex);
    t_isBoolRead = advance_to(target);
    if (!t_isBoolRead) return nullptr;
    return current_record();
}

void TractorHybridClass::set_iterator(const std::string& chr, int64_t beg_pos, int64_t end_pos) {
    streaming_mode_ = true;
    has_region_filter_ = !chr.empty();
    region_chr_ = chr;
    region_start_ = beg_pos <= 0 ? 1 : beg_pos;
    region_end_ = end_pos <= 0 ? std::numeric_limits<int64_t>::max() : end_pos;
    advance_to_next_in_filter();
}

void TractorHybridClass::move_forward_iterator(int n_marker) {
    for (int i = 0; i < n_marker; ++i) {
        if (!current_record()) return;
        advance_current();
        advance_to_next_in_filter();
    }
}

bool TractorHybridClass::check_iterator_end() {
    advance_to_next_in_filter();
    const VariantRecord* record = current_record();
    if (!record) return true;
    if (!streaming_mode_ || !has_region_filter_) return false;
    if (record->chr != region_chr_) return true;
    if (record->pos < region_start_) return false;
    return record->pos > region_end_;
}

uint32_t TractorHybridClass::blockIdForMarker(uint64_t t_gIndex) {
    bool is_read = false;
    const VariantRecord* record = record_for_request(t_gIndex, is_read);
    if (!is_read || !record) return UINT32_MAX;
    return block_id_for_record(*record);
}

double TractorHybridClass::get_io_seconds() const {
    return io_seconds_;
}

double TractorHybridClass::get_decode_seconds() const {
    return decode_seconds_;
}

void TractorHybridClass::tracked_seek(
    FILE* fp,
    uint64_t offset,
    const std::string& path,
    uint64_t& current_offset
) {
    if (current_offset == offset) return;
    ScopedSeconds timer(io_seconds_);
    seek_or_stop(fp, offset, path);
    current_offset = offset;
}

void TractorHybridClass::tracked_read_exact(
    FILE* fp,
    void* data,
    size_t bytes,
    const std::string& path,
    uint64_t& current_offset
) {
    ScopedSeconds timer(io_seconds_);
    read_exact(fp, data, bytes, path);
    current_offset += bytes;
}

size_t TractorHybridClass::tracked_fread(
    FILE* fp,
    void* data,
    size_t size,
    size_t count,
    const std::string& path,
    uint64_t& current_offset
) {
    ScopedSeconds timer(io_seconds_);
    size_t got = std::fread(data, size, count, fp);
    current_offset += got * size;
    if (got != count) {
        Rcpp::stop("Failed reading %s", path.c_str());
    }
    return got;
}

uint32_t TractorHybridClass::tracked_read_u32_le(
    FILE* fp,
    const std::string& path,
    uint64_t& current_offset
) {
    uint8_t b[4];
    tracked_read_exact(fp, b, sizeof(b), path, current_offset);
    return static_cast<uint32_t>(b[0]) |
           (static_cast<uint32_t>(b[1]) << 8) |
           (static_cast<uint32_t>(b[2]) << 16) |
           (static_cast<uint32_t>(b[3]) << 24);
}

uint64_t TractorHybridClass::tracked_read_u64_le(
    FILE* fp,
    const std::string& path,
    uint64_t& current_offset
) {
    uint8_t b[8];
    tracked_read_exact(fp, b, sizeof(b), path, current_offset);
    return static_cast<uint64_t>(b[0]) |
           (static_cast<uint64_t>(b[1]) << 8) |
           (static_cast<uint64_t>(b[2]) << 16) |
           (static_cast<uint64_t>(b[3]) << 24) |
           (static_cast<uint64_t>(b[4]) << 32) |
           (static_cast<uint64_t>(b[5]) << 40) |
           (static_cast<uint64_t>(b[6]) << 48) |
           (static_cast<uint64_t>(b[7]) << 56);
}

uint32_t TractorHybridClass::block_id_for_record(const VariantRecord& record) {
    if (record.common) return record.block_id;

    auto covers = [&record](const AncBlockRecord& block) {
        return block.chr == record.chr &&
               block.start_pos <= record.pos &&
               record.pos <= block.end_pos;
    };

    if (!anc_blocks_.empty() && anc_block_cursor_ < anc_blocks_.size() &&
        covers(anc_blocks_[anc_block_cursor_])) {
        return anc_blocks_[anc_block_cursor_].block_id;
    }

    while (anc_block_cursor_ < anc_blocks_.size()) {
        const AncBlockRecord& block = anc_blocks_[anc_block_cursor_];
        if (covers(block)) return block.block_id;
        if (block.chr == record.chr && block.end_pos < record.pos) {
            ++anc_block_cursor_;
            continue;
        }
        break;
    }

    for (size_t i = 0; i < anc_blocks_.size(); ++i) {
        if (covers(anc_blocks_[i])) {
            anc_block_cursor_ = i;
            return anc_blocks_[i].block_id;
        }
    }

    Rcpp::stop(
        "No ancestry block covers hybrid variant %s:%lld (%s)",
        record.chr.c_str(),
        static_cast<long long>(record.pos),
        record.id.c_str()
    );
    return UINT32_MAX;
}

void TractorHybridClass::fill_marker_fields(
    const VariantRecord& record,
    std::string& t_ref,
    std::string& t_alt,
    std::string& t_marker,
    uint32_t& t_pd,
    std::string& t_chr,
    double& t_missingRate,
    double& t_imputeInfo,
    bool t_isImputation
) const {
    if (record.pos < 0 || record.pos > std::numeric_limits<uint32_t>::max()) {
        Rcpp::stop("Hybrid marker position exceeds uint32_t capacity at global variant %u", record.global_variant_index);
    }

    t_chr = record.chr;
    t_pd = static_cast<uint32_t>(record.pos);
    t_marker = record.id;
    t_ref = record.ref;
    t_alt = record.alt;
    t_missingRate = 0.0;
    t_imputeInfo = t_isImputation ? 1.0 : 1.0;
}

void TractorHybridClass::read_dosage(
    const VariantRecord& record,
    arma::vec& dosages,
    std::vector<uint>& indexForNonZero,
    double& altCounts
) {
    ensure_marker_dosage_cache(record);
    dosages = cached_total_dosage_;
    indexForNonZero = cached_total_nonzero_;
    altCounts = cached_total_alt_counts_;
}

void TractorHybridClass::load_common_words(const VariantRecord& record) {
    if (!record.common) {
        Rcpp::stop("Internal error: load_common_words called on a rare hybrid marker");
    }
    if (cached_common_global_variant_index_ == record.global_variant_index &&
        common_words_.size() == static_cast<size_t>(meta_->n_words)) {
        return;
    }

    common_words_.assign(static_cast<size_t>(meta_->n_words), 0);
    tracked_seek(common_bin_, record.payload_offset, prefix_ + ".common.geno.bin", common_bin_pos_);
    tracked_fread(
        common_bin_,
        common_words_.data(),
        sizeof(uint64_t),
        static_cast<size_t>(meta_->n_words),
        prefix_ + ".common.geno.bin",
        common_bin_pos_
    );
    cached_common_global_variant_index_ = record.global_variant_index;
}

void TractorHybridClass::load_rare_carriers(const VariantRecord& record) {
    if (record.common) {
        Rcpp::stop("Internal error: load_rare_carriers called on a common hybrid marker");
    }
    if (cached_rare_global_variant_index_ == record.global_variant_index) return;

    cached_rare_carriers_.assign(record.n_carriers, RareCarrierPacked{});
    tracked_seek(rare_bin_, record.payload_offset, prefix_ + ".rare.carrier.bin", rare_bin_pos_);
    if (!cached_rare_carriers_.empty()) {
        tracked_read_exact(
            rare_bin_,
            cached_rare_carriers_.data(),
            cached_rare_carriers_.size() * sizeof(RareCarrierPacked),
            prefix_ + ".rare.carrier.bin",
            rare_bin_pos_
        );
    }

    for (const RareCarrierPacked& carrier : cached_rare_carriers_) {
        if (carrier.pos_index != record.global_variant_index) {
            Rcpp::stop(
                "Rare carrier pos_index mismatch: got %u, expected %u",
                carrier.pos_index,
                record.global_variant_index
            );
        }
        uint32_t ancestry = carrier.anc_hap >> 27;
        uint32_t hap_id = carrier.anc_hap & kPackedHapMask;
        if (hap_id >= meta_->n_haps) Rcpp::stop("Rare carrier hap_id %u exceeds n_haps", hap_id);
        if (ancestry >= meta_->n_ancestries) Rcpp::stop("Rare carrier ancestry %u exceeds n_ancestries", ancestry);
    }

    cached_rare_global_variant_index_ = record.global_variant_index;
}

void TractorHybridClass::ensure_marker_dosage_cache(const VariantRecord& record) {
    if (cached_dosage_global_variant_index_ == record.global_variant_index) return;

    if (record.common) {
        load_ancestry_block(record.block_id);
        load_common_words(record);
    } else {
        load_rare_carriers(record);
    }

    ScopedSeconds timer(decode_seconds_);
    const uint32_t n_model = getN();
    const arma::uword n_anc = static_cast<arma::uword>(meta_->n_ancestries);

    cached_total_dosage_.zeros(n_model);
    cached_dosage_by_ancestry_.zeros(n_model, n_anc);
    cached_total_nonzero_.clear();
    cached_dosage_nonzero_by_ancestry_.assign(static_cast<size_t>(meta_->n_ancestries), std::vector<uint>{});
    cached_dosage_alt_counts_by_ancestry_.assign(static_cast<size_t>(meta_->n_ancestries), 0.0);
    cached_total_alt_counts_ = 0.0;

    auto add_alt_hap = [&](uint32_t model_i, uint32_t ancestry) {
        if (cached_total_dosage_[model_i] == 0.0) {
            cached_total_nonzero_.push_back(model_i);
        }
        if (cached_dosage_by_ancestry_(model_i, ancestry) == 0.0) {
            cached_dosage_nonzero_by_ancestry_[ancestry].push_back(model_i);
        }
        cached_total_dosage_[model_i] += 1.0;
        cached_dosage_by_ancestry_(model_i, ancestry) += 1.0;
        cached_total_alt_counts_ += 1.0;
        cached_dosage_alt_counts_by_ancestry_[ancestry] += 1.0;
    };

    if (record.common) {
        for (uint32_t model_i = 0; model_i < n_model; ++model_i) {
            uint64_t raw_i = static_cast<uint64_t>(model_to_raw_sample_[model_i]);
            uint64_t h0 = 2 * raw_i;
            uint64_t h1 = h0 + 1;

            if (bit_at(common_words_, h0)) {
                uint32_t ancestry = hap_ancestry_[static_cast<size_t>(h0)];
                add_alt_hap(model_i, ancestry);
            }
            if (bit_at(common_words_, h1)) {
                uint32_t ancestry = hap_ancestry_[static_cast<size_t>(h1)];
                add_alt_hap(model_i, ancestry);
            }
        }
    } else {
        for (const RareCarrierPacked& carrier : cached_rare_carriers_) {
            uint32_t ancestry = carrier.anc_hap >> 27;
            uint32_t hap_id = carrier.anc_hap & kPackedHapMask;
            uint32_t raw_sample = hap_id >> 1;
            int32_t model_i = raw_sample_to_model_[raw_sample];
            if (model_i < 0) continue;
            add_alt_hap(static_cast<uint32_t>(model_i), ancestry);
        }
    }

    cached_dosage_global_variant_index_ = record.global_variant_index;
}

void TractorHybridClass::load_ancestry_block(uint32_t block_id) {
    if (cached_anc_block_id_ == block_id) return;

    const std::string idx_path = prefix_ + ".ancblock.idx";
    const std::string bin_path = prefix_ + ".ancblock.bin";
    uint64_t idx_offset = 8ULL + static_cast<uint64_t>(block_id) * 20ULL;
    tracked_seek(anc_idx_, idx_offset, idx_path, anc_idx_pos_);

    uint32_t got_block_id = tracked_read_u32_le(anc_idx_, idx_path, anc_idx_pos_);
    (void)tracked_read_u64_le(anc_idx_, idx_path, anc_idx_pos_);
    uint64_t anc_offset = tracked_read_u64_le(anc_idx_, idx_path, anc_idx_pos_);
    if (got_block_id != block_id) {
        Rcpp::stop("Ancestry block index mismatch: got %u, expected %u", got_block_id, block_id);
    }

    std::vector<uint64_t> mask(static_cast<size_t>(meta_->n_words), 0);
    hap_ancestry_.assign(static_cast<size_t>(meta_->n_haps), UINT8_MAX);

    tracked_seek(anc_bin_, anc_offset, bin_path, anc_bin_pos_);
    for (uint32_t ancestry = 0; ancestry < meta_->n_ancestries; ++ancestry) {
        tracked_fread(
            anc_bin_,
            mask.data(),
            sizeof(uint64_t),
            static_cast<size_t>(meta_->n_words),
            bin_path,
            anc_bin_pos_
        );

        ScopedSeconds timer(decode_seconds_);
        for (uint64_t word_i = 0; word_i < meta_->n_words; ++word_i) {
            uint64_t word = mask[static_cast<size_t>(word_i)];
            while (word != 0) {
                unsigned bit = static_cast<unsigned>(__builtin_ctzll(word));
                uint64_t hap_id = (word_i << 6) + bit;
                if (hap_id < meta_->n_haps) {
                    hap_ancestry_[static_cast<size_t>(hap_id)] = static_cast<uint8_t>(ancestry);
                }
                word &= word - 1;
            }
        }
    }

    for (uint64_t hap_id = 0; hap_id < meta_->n_haps; ++hap_id) {
        if (hap_ancestry_[static_cast<size_t>(hap_id)] == UINT8_MAX) {
            Rcpp::stop("Ancestry block %u has no ancestry assignment for hap_id %llu",
                       block_id,
                       static_cast<unsigned long long>(hap_id));
        }
    }

    cached_anc_block_id_ = block_id;
}

void TractorHybridClass::read_dosage_for_ancestry(
    const VariantRecord& record,
    uint32_t ancestry,
    arma::vec& dosages,
    std::vector<uint>& indexForNonZero,
    double& altCounts
) {
    if (ancestry >= meta_->n_ancestries) {
        Rcpp::stop("Requested ancestry %u exceeds n_ancestries", ancestry);
    }

    ensure_marker_dosage_cache(record);
    dosages = cached_dosage_by_ancestry_.col(ancestry);
    indexForNonZero = cached_dosage_nonzero_by_ancestry_[ancestry];
    altCounts = cached_dosage_alt_counts_by_ancestry_[ancestry];
}

void TractorHybridClass::read_ancestry_count(
    const VariantRecord& record,
    uint32_t ancestry,
    arma::vec& ancestryCounts,
    std::vector<uint>& indexForNonZero,
    double& altCounts
) {
    if (ancestry >= meta_->n_ancestries) {
        Rcpp::stop("Requested ancestry %u exceeds n_ancestries", ancestry);
    }

    uint32_t block_id = block_id_for_record(record);
    ensure_ancestry_count_cache(block_id);
    ancestryCounts = cached_ancestry_counts_.col(ancestry);
    indexForNonZero = cached_ancestry_nonzero_by_ancestry_[ancestry];
    altCounts = cached_ancestry_counts_by_ancestry_[ancestry];
}

void TractorHybridClass::ensure_ancestry_count_cache(uint32_t block_id) {
    if (cached_ancestry_counts_block_id_ == block_id) return;

    load_ancestry_block(block_id);
    ScopedSeconds timer(decode_seconds_);
    const uint32_t n_model = getN();
    const arma::uword n_anc = static_cast<arma::uword>(meta_->n_ancestries);

    cached_ancestry_counts_.zeros(n_model, n_anc);
    cached_ancestry_nonzero_by_ancestry_.assign(static_cast<size_t>(meta_->n_ancestries), std::vector<uint>{});
    cached_ancestry_counts_by_ancestry_.assign(static_cast<size_t>(meta_->n_ancestries), 0.0);

    for (uint32_t model_i = 0; model_i < n_model; ++model_i) {
        uint64_t raw_i = static_cast<uint64_t>(model_to_raw_sample_[model_i]);
        uint64_t h0 = 2 * raw_i;
        uint64_t h1 = h0 + 1;

        uint32_t a0 = hap_ancestry_[static_cast<size_t>(h0)];
        uint32_t a1 = hap_ancestry_[static_cast<size_t>(h1)];

        if (cached_ancestry_counts_(model_i, a0) == 0.0) {
            cached_ancestry_nonzero_by_ancestry_[a0].push_back(model_i);
        }
        cached_ancestry_counts_(model_i, a0) += 1.0;
        cached_ancestry_counts_by_ancestry_[a0] += 1.0;

        if (cached_ancestry_counts_(model_i, a1) == 0.0) {
            cached_ancestry_nonzero_by_ancestry_[a1].push_back(model_i);
        }
        cached_ancestry_counts_(model_i, a1) += 1.0;
        cached_ancestry_counts_by_ancestry_[a1] += 1.0;
    }

    cached_ancestry_counts_block_id_ = block_id;
}

void TractorHybridClass::read_dosage_by_ancestry(
    const VariantRecord& record,
    arma::vec& dosages,
    arma::mat& dosageByAncestry,
    std::vector<uint>& indexForNonZero,
    double& altCounts
) {
    ensure_marker_dosage_cache(record);
    dosages = cached_total_dosage_;
    dosageByAncestry = cached_dosage_by_ancestry_;
    indexForNonZero = cached_total_nonzero_;
    altCounts = cached_total_alt_counts_;
}

void TractorHybridClass::getOneMarker(
    uint64_t& t_gIndex_prev,
    uint64_t& t_gIndex,
    std::string& t_ref,
    std::string& t_alt,
    std::string& t_marker,
    uint32_t& t_pd,
    std::string& t_chr,
    double& t_altFreq,
    double& t_altCounts,
    double& t_missingRate,
    double& t_imputeInfo,
    bool& t_isOutputIndexForMissing,
    std::vector<uint>& t_indexForMissing,
    bool& t_isOnlyOutputNonZero,
    std::vector<uint>& t_indexForNonZero,
    bool& t_isBoolRead,
    arma::vec& t_GVec,
    bool t_isImputation
) {
    (void)t_gIndex_prev;
    const VariantRecord* record_ptr = record_for_request(t_gIndex, t_isBoolRead);
    if (!record_ptr) return;
    const VariantRecord& record = *record_ptr;
    fill_marker_fields(record, t_ref, t_alt, t_marker, t_pd, t_chr, t_missingRate, t_imputeInfo, t_isImputation);

    t_indexForMissing.clear();
    read_dosage(record, t_GVec, t_indexForNonZero, t_altCounts);
    t_altFreq = getN() == 0 ? 0.0 : t_altCounts / (2.0 * static_cast<double>(getN()));
    t_missingRate = 0.0;

    if (!t_isOutputIndexForMissing) {
        t_indexForMissing.clear();
    }

    if (t_isOnlyOutputNonZero) {
        arma::vec nonzero(t_indexForNonZero.size());
        for (size_t i = 0; i < t_indexForNonZero.size(); ++i) {
            nonzero[i] = t_GVec[t_indexForNonZero[i]];
        }
        t_GVec = nonzero;
    }

    if (streaming_mode_) {
        move_forward_iterator(1);
    }
}

void TractorHybridClass::getOneMarkerAncestry(
    uint64_t& t_gIndex_prev,
    uint64_t& t_gIndex,
    std::string& t_ref,
    std::string& t_alt,
    std::string& t_marker,
    uint32_t& t_pd,
    std::string& t_chr,
    double& t_altFreq,
    double& t_altCounts,
    double& t_missingRate,
    double& t_imputeInfo,
    bool& t_isOutputIndexForMissing,
    std::vector<uint>& t_indexForMissing,
    bool& t_isOnlyOutputNonZero,
    std::vector<uint>& t_indexForNonZero,
    bool& t_isBoolRead,
    arma::vec& t_GVec,
    arma::mat& t_GByAncestry,
    bool t_isImputation
) {
    (void)t_gIndex_prev;
    const VariantRecord* record_ptr = record_for_request(t_gIndex, t_isBoolRead);
    if (!record_ptr) return;
    const VariantRecord& record = *record_ptr;
    fill_marker_fields(record, t_ref, t_alt, t_marker, t_pd, t_chr, t_missingRate, t_imputeInfo, t_isImputation);

    t_indexForMissing.clear();
    read_dosage_by_ancestry(record, t_GVec, t_GByAncestry, t_indexForNonZero, t_altCounts);
    t_altFreq = getN() == 0 ? 0.0 : t_altCounts / (2.0 * static_cast<double>(getN()));
    t_missingRate = 0.0;

    if (!t_isOutputIndexForMissing) {
        t_indexForMissing.clear();
    }

    if (t_isOnlyOutputNonZero) {
        arma::vec nonzero(t_indexForNonZero.size());
        arma::mat nonzeroByAncestry(t_indexForNonZero.size(), t_GByAncestry.n_cols);
        for (size_t i = 0; i < t_indexForNonZero.size(); ++i) {
            nonzero[i] = t_GVec[t_indexForNonZero[i]];
            nonzeroByAncestry.row(i) = t_GByAncestry.row(t_indexForNonZero[i]);
        }
        t_GVec = nonzero;
        t_GByAncestry = nonzeroByAncestry;
    }
}

void TractorHybridClass::getOneMarkerAdmixedField(
    uint64_t& t_gIndex_prev,
    uint64_t& t_gIndex,
    std::string& t_ref,
    std::string& t_alt,
    std::string& t_marker,
    uint32_t& t_pd,
    std::string& t_chr,
    double& t_altFreq,
    double& t_altCounts,
    double& t_missingRate,
    double& t_imputeInfo,
    bool& t_isOutputIndexForMissing,
    std::vector<uint>& t_indexForMissing,
    bool& t_isOnlyOutputNonZero,
    std::vector<uint>& t_indexForNonZero,
    bool& t_isBoolRead,
    arma::vec& t_GVec,
    bool t_isImputation,
    const std::string& t_vcfField
) {
    (void)t_gIndex_prev;
    const VariantRecord* record_ptr = record_for_request(t_gIndex, t_isBoolRead);
    if (!record_ptr) return;
    const VariantRecord& record = *record_ptr;
    fill_marker_fields(record, t_ref, t_alt, t_marker, t_pd, t_chr, t_missingRate, t_imputeInfo, t_isImputation);

    t_indexForMissing.clear();

    uint32_t ancestry = 0;
    if (t_vcfField == "DSALL" || t_vcfField == "DS" || t_vcfField == "GT") {
        read_dosage(record, t_GVec, t_indexForNonZero, t_altCounts);
    } else if (parse_one_based_suffix(t_vcfField, "DS", static_cast<uint32_t>(meta_->n_ancestries), ancestry)) {
        read_dosage_for_ancestry(record, ancestry, t_GVec, t_indexForNonZero, t_altCounts);
    } else if (parse_one_based_suffix(t_vcfField, "ANC", static_cast<uint32_t>(meta_->n_ancestries), ancestry)) {
        read_ancestry_count(record, ancestry, t_GVec, t_indexForNonZero, t_altCounts);
    } else {
        Rcpp::stop(
            "Unsupported hybrid admixed FORMAT field %s; expected ANC#, DS#, or DSALL",
            t_vcfField.c_str()
        );
    }

    t_altFreq = getN() == 0 ? 0.0 : t_altCounts / (2.0 * static_cast<double>(getN()));
    t_missingRate = 0.0;

    if (!t_isOutputIndexForMissing) {
        t_indexForMissing.clear();
    }

    if (t_isOnlyOutputNonZero) {
        arma::vec nonzero(t_indexForNonZero.size());
        for (size_t i = 0; i < t_indexForNonZero.size(); ++i) {
            nonzero[i] = t_GVec[t_indexForNonZero[i]];
        }
        t_GVec = nonzero;
    }
}

} // namespace TRACTOR_HYBRID
