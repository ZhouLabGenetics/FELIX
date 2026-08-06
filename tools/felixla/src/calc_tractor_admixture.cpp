// calc_tractor_admixture.cpp
//
// Compute per-sample global ancestry proportions from FLARE local ancestry
// VCFs with FORMAT/AN1 and FORMAT/AN2, or old SAIGE-TRACTOR dosage VCFs
// with scalar FORMAT/ANC1, FORMAT/ANC2, ... fields.

#include <htslib/hts.h>
#include <htslib/vcf.h>
#if __has_include(<htslib/hts_log.h>)
#include <htslib/hts_log.h>
#define HAVE_HTSLIB_LOG 1
#endif

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

struct Args {
    std::vector<std::string> vcfs;
    std::string pattern;
    std::string chroms = "1-22";
    std::vector<std::string> requested_fields;
    std::string output;
    uint64_t progress_every = 1000000;
    uint64_t thin = 1;
};

struct FieldInfo {
    std::string id;
    std::string description;
    int anc_number = 0;
};

struct Totals {
    std::vector<std::string> samples;
    std::unordered_map<std::string, size_t> sample_to_index;
    std::vector<std::vector<double>> counts;
    std::vector<double> denominators;
};

struct Stats {
    uint64_t records = 0;
    uint64_t used_records = 0;
    uint64_t skipped_format_rows = 0;
    uint64_t skipped_sample_values = 0;
    uint64_t skipped_short_rows = 0;
};

struct IndexRecordCount {
    bool available = false;
    uint64_t total = 0;
};

static void usage(std::ostream& os) {
    os << "Usage:\n"
       << "  calc_tractor_admixture [VCF ...] [--pattern TEMPLATE] [options]\n\n"
       << "Options:\n"
       << "  --pattern TEMPLATE       Input template, e.g. 'tractor/chr{chr}.vcf.gz'\n"
       << "  --chroms LIST            Chromosomes for --pattern, default: 1-22\n"
       << "  --anc-fields LIST        Comma-separated old dosage fields; default auto-detects FLARE AN1/AN2 or ANC1,ANC2,...\n"
       << "  -o, --output PATH        Output TSV, default: stdout\n"
       << "  --progress-every N       Report progress every N scanned records; percent/ETA if indexed; default 1000000; 0 disables\n"
       << "  --thin N                 Use every Nth variant; default 1 is exact\n"
       << "  -h, --help               Show this help\n";
}

[[noreturn]] static void die(const std::string& msg) {
    throw std::runtime_error(msg);
}

#ifdef HTS_IDX_SILENT_FAIL
static hts_idx_t* load_any_vcf_index(const std::string& path) {
    hts_idx_t* idx = hts_idx_load3(path.c_str(), nullptr, HTS_FMT_CSI, HTS_IDX_SILENT_FAIL);
    if (idx) return idx;
    return hts_idx_load3(path.c_str(), nullptr, HTS_FMT_TBI, HTS_IDX_SILENT_FAIL);
}
#endif

static IndexRecordCount get_index_record_count(const std::string& path) {
    IndexRecordCount result;
#ifdef HTS_IDX_SILENT_FAIL
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
            return IndexRecordCount{};
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
#endif
    return result;
}

static std::string format_duration(double seconds) {
    if (!std::isfinite(seconds) || seconds < 0.0) return "?";
    uint64_t rounded = static_cast<uint64_t>(seconds + 0.5);
    uint64_t h = rounded / 3600;
    uint64_t m = (rounded % 3600) / 60;
    uint64_t s = rounded % 60;

    std::ostringstream out;
    if (h > 0) {
        out << h << "h" << std::setw(2) << std::setfill('0') << m
            << "m" << std::setw(2) << std::setfill('0') << s << "s";
    } else if (m > 0) {
        out << m << "m" << std::setw(2) << std::setfill('0') << s << "s";
    } else {
        out << s << "s";
    }
    return out.str();
}

class ProgressReporter {
public:
    ProgressReporter(
        std::string path,
        size_t file_index,
        size_t file_count,
        uint64_t progress_every
    )
        : path_(std::move(path)),
          file_index_(file_index),
          file_count_(file_count),
          progress_every_(progress_every),
          total_(progress_every > 0 ? get_index_record_count(path_) : IndexRecordCount{}),
          start_(std::chrono::steady_clock::now()),
          last_report_(start_),
          next_record_report_(progress_every > 0 ? progress_every : 0) {
        if (progress_every_ == 0) return;
        std::cerr << "[progress] file " << file_index_ << "/" << file_count_ << " "
                  << path_ << ": ";
        if (total_.available) {
            std::cerr << "index reports " << total_.total << " records\n";
        } else {
            std::cerr << "index record count unavailable; percent/ETA disabled; "
                      << "reporting scanned records and speed\n";
        }
    }

    void maybe_report(uint64_t scanned_records, uint64_t used_records, bool force = false) {
        if (progress_every_ == 0 || scanned_records == 0) return;

        auto now = std::chrono::steady_clock::now();
        bool count_due = scanned_records >= next_record_report_;
        bool time_due = now - last_report_ >= std::chrono::seconds(5);
        if (!force && !count_due && !time_due) return;

        while (next_record_report_ <= scanned_records) {
            next_record_report_ += progress_every_;
        }
        last_report_ = now;

        double elapsed = std::chrono::duration<double>(now - start_).count();
        double rate = elapsed > 0.0 ? static_cast<double>(scanned_records) / elapsed : 0.0;

        std::cerr << "[progress] file " << file_index_ << "/" << file_count_ << " "
                  << path_ << ": " << scanned_records;
        if (total_.available && total_.total > 0) {
            double pct = 100.0 * static_cast<double>(scanned_records) /
                         static_cast<double>(total_.total);
            if (pct > 100.0) pct = 100.0;
            double eta = rate > 0.0
                ? (static_cast<double>(total_.total) - static_cast<double>(scanned_records)) / rate
                : std::numeric_limits<double>::quiet_NaN();
            if (eta < 0.0) eta = 0.0;
            std::cerr << "/" << total_.total << " ("
                      << std::fixed << std::setprecision(1) << pct << "%"
                      << std::defaultfloat << ")";
            std::cerr << ", ETA " << format_duration(eta);
        }
        std::cerr << ", used " << used_records
                  << ", " << std::fixed << std::setprecision(0) << rate
                  << std::defaultfloat << " records/s"
                  << ", elapsed " << format_duration(elapsed) << "\n";
    }

private:
    std::string path_;
    size_t file_index_ = 0;
    size_t file_count_ = 0;
    uint64_t progress_every_ = 0;
    IndexRecordCount total_;
    std::chrono::steady_clock::time_point start_;
    std::chrono::steady_clock::time_point last_report_;
    uint64_t next_record_report_ = 0;
};

static bool starts_with(const std::string& text, const std::string& prefix) {
    return text.rfind(prefix, 0) == 0;
}

static uint64_t parse_u64(const std::string& text, const std::string& opt) {
    if (text.empty()) die(opt + " requires a positive integer");
    size_t pos = 0;
    unsigned long long value = 0;
    try {
        value = std::stoull(text, &pos, 10);
    } catch (const std::exception&) {
        die(opt + " requires a positive integer, got '" + text + "'");
    }
    if (pos != text.size()) die(opt + " requires an integer, got '" + text + "'");
    return static_cast<uint64_t>(value);
}

static std::string require_value(int& i, int argc, char** argv, const std::string& opt) {
    if (i + 1 >= argc) die(opt + " requires a value");
    ++i;
    return argv[i];
}

static std::vector<std::string> split_csv(const std::string& text) {
    std::vector<std::string> out;
    size_t start = 0;
    while (start <= text.size()) {
        size_t comma = text.find(',', start);
        std::string part = text.substr(start, comma == std::string::npos ? comma : comma - start);
        size_t first = part.find_first_not_of(" \t\r\n");
        size_t last = part.find_last_not_of(" \t\r\n");
        if (first != std::string::npos) out.push_back(part.substr(first, last - first + 1));
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    return out;
}

static Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            usage(std::cout);
            std::exit(0);
        } else if (arg == "--pattern") {
            args.pattern = require_value(i, argc, argv, arg);
        } else if (starts_with(arg, "--pattern=")) {
            args.pattern = arg.substr(std::strlen("--pattern="));
        } else if (arg == "--chroms") {
            args.chroms = require_value(i, argc, argv, arg);
        } else if (starts_with(arg, "--chroms=")) {
            args.chroms = arg.substr(std::strlen("--chroms="));
        } else if (arg == "--anc-fields") {
            args.requested_fields = split_csv(require_value(i, argc, argv, arg));
        } else if (starts_with(arg, "--anc-fields=")) {
            args.requested_fields = split_csv(arg.substr(std::strlen("--anc-fields=")));
        } else if (arg == "-o" || arg == "--output") {
            args.output = require_value(i, argc, argv, arg);
        } else if (starts_with(arg, "--output=")) {
            args.output = arg.substr(std::strlen("--output="));
        } else if (arg == "--progress-every") {
            args.progress_every = parse_u64(require_value(i, argc, argv, arg), arg);
        } else if (starts_with(arg, "--progress-every=")) {
            args.progress_every = parse_u64(arg.substr(std::strlen("--progress-every=")), "--progress-every");
        } else if (arg == "--thin") {
            args.thin = parse_u64(require_value(i, argc, argv, arg), arg);
        } else if (starts_with(arg, "--thin=")) {
            args.thin = parse_u64(arg.substr(std::strlen("--thin=")), "--thin");
        } else if (!arg.empty() && arg[0] == '-') {
            die("unknown option: " + arg);
        } else {
            args.vcfs.push_back(arg);
        }
    }

    if (args.thin < 1) die("--thin must be >= 1");
    if (args.vcfs.empty() && args.pattern.empty()) {
        die("provide input VCFs or --pattern");
    }
    return args;
}

static bool is_integer(const std::string& text) {
    if (text.empty()) return false;
    size_t i = 0;
    if (text[0] == '-' || text[0] == '+') i = 1;
    if (i == text.size()) return false;
    for (; i < text.size(); ++i) {
        if (text[i] < '0' || text[i] > '9') return false;
    }
    return true;
}

static std::vector<std::string> parse_chroms(const std::string& chroms) {
    std::vector<std::string> out;
    for (const std::string& part_raw : split_csv(chroms)) {
        size_t dash = part_raw.find('-');
        if (dash != std::string::npos) {
            std::string a = part_raw.substr(0, dash);
            std::string b = part_raw.substr(dash + 1);
            if (!is_integer(a) || !is_integer(b)) {
                die("--chroms ranges must be numeric, got '" + part_raw + "'");
            }
            int start = std::stoi(a);
            int end = std::stoi(b);
            int step = start <= end ? 1 : -1;
            for (int chr = start;; chr += step) {
                out.push_back(std::to_string(chr));
                if (chr == end) break;
            }
        } else {
            out.push_back(part_raw);
        }
    }
    return out;
}

static void replace_all(std::string& text, const std::string& from, const std::string& to) {
    if (from.empty()) return;
    size_t pos = 0;
    while ((pos = text.find(from, pos)) != std::string::npos) {
        text.replace(pos, from.size(), to);
        pos += to.size();
    }
}

static std::vector<std::string> expand_inputs(const Args& args) {
    std::vector<std::string> inputs = args.vcfs;
    if (!args.pattern.empty()) {
        for (const std::string& chr : parse_chroms(args.chroms)) {
            std::string path = args.pattern;
            replace_all(path, "{chr}", chr);
            replace_all(path, "{chrom}", chr);
            replace_all(path, "{chromosome}", chr);
            replace_all(path, "{CHR}", chr);
            inputs.push_back(path);
        }
    }
    return inputs;
}

static bool file_exists(const std::string& path) {
    std::ifstream f(path);
    return f.good();
}

static bool parse_anc_number(const std::string& field, int& number) {
    if (field.size() <= 3 || field[0] != 'A' || field[1] != 'N' || field[2] != 'C') {
        return false;
    }
    int value = 0;
    for (size_t i = 3; i < field.size(); ++i) {
        if (field[i] < '0' || field[i] > '9') return false;
        value = value * 10 + (field[i] - '0');
    }
    number = value;
    return true;
}

static std::string lower_ascii(std::string text) {
    for (char& c : text) {
        if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
    }
    return text;
}

static bool label_char(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
           (c >= '0' && c <= '9') || c == '_' || c == '.' || c == '-';
}

static std::string extract_after_phrase(const std::string& desc, const std::string& phrase) {
    std::string lower = lower_ascii(desc);
    std::string lower_phrase = lower_ascii(phrase);
    size_t pos = lower.find(lower_phrase);
    if (pos == std::string::npos) return "";
    pos += phrase.size();
    while (pos < desc.size() && (desc[pos] == ' ' || desc[pos] == '\t')) ++pos;
    size_t end = pos;
    while (end < desc.size() && label_char(desc[end])) ++end;
    return desc.substr(pos, end - pos);
}

static std::string label_from_description(const std::string& field, const std::string& description) {
    std::string label = extract_after_phrase(description, "Ancestry of ");
    if (!label.empty()) return label;

    label = extract_after_phrase(description, "Dosage of ");
    int ignored = 0;
    if (!label.empty() && !parse_anc_number(label, ignored)) return label;
    return field;
}

static const char* hrec_value(bcf_hrec_t* hrec, const char* key) {
    int idx = bcf_hrec_find_key(hrec, key);
    if (idx < 0 || idx >= hrec->nkeys || hrec->vals[idx] == nullptr) return "";
    return hrec->vals[idx];
}

static bool has_format_field(bcf_hdr_t* hdr, const char* field) {
    for (int i = 0; i < hdr->nhrec; ++i) {
        bcf_hrec_t* hrec = hdr->hrec[i];
        if (!hrec || !hrec->key || std::strcmp(hrec->key, "FORMAT") != 0) continue;
        if (std::strcmp(hrec_value(hrec, "ID"), field) == 0) return true;
    }
    return false;
}

static std::vector<std::string> parse_flare_ancestry_labels(bcf_hdr_t* hdr) {
    std::vector<std::string> labels;
    for (int i = 0; i < hdr->nhrec; ++i) {
        bcf_hrec_t* hrec = hdr->hrec[i];
        if (!hrec || !hrec->key || std::strcmp(hrec->key, "ANCESTRY") != 0) continue;
        for (int j = 0; j < hrec->nkeys; ++j) {
            if (!hrec->keys[j] || !hrec->vals[j]) continue;
            char* end = nullptr;
            errno = 0;
            long value = std::strtol(hrec->vals[j], &end, 10);
            if (errno != 0 || end == hrec->vals[j] || value < 0) continue;
            size_t idx = static_cast<size_t>(value);
            if (labels.size() <= idx) labels.resize(idx + 1);
            labels[idx] = hrec->keys[j];
        }
    }
    for (size_t i = 0; i < labels.size(); ++i) {
        if (labels[i].empty()) labels[i] = "ANC" + std::to_string(i);
    }
    return labels;
}

static void ensure_ancestry_size(
    Totals& totals,
    std::vector<std::string>& labels,
    size_t n_anc
) {
    if (labels.size() < n_anc) {
        size_t old_size = labels.size();
        labels.resize(n_anc);
        for (size_t i = old_size; i < n_anc; ++i) {
            labels[i] = "ANC" + std::to_string(i);
        }
    }
    for (std::vector<double>& row : totals.counts) {
        if (row.size() < n_anc) row.resize(n_anc, 0.0);
    }
}

static bool is_missing_int32(int32_t value) {
    return value == bcf_int32_missing || value == bcf_int32_vector_end || value < 0;
}

static std::vector<FieldInfo> detect_anc_fields(bcf_hdr_t* hdr) {
    std::vector<FieldInfo> fields;
    for (int i = 0; i < hdr->nhrec; ++i) {
        bcf_hrec_t* hrec = hdr->hrec[i];
        if (!hrec || !hrec->key || std::strcmp(hrec->key, "FORMAT") != 0) continue;
        std::string id = hrec_value(hrec, "ID");
        int anc_number = 0;
        if (!parse_anc_number(id, anc_number)) continue;
        fields.push_back(FieldInfo{id, hrec_value(hrec, "Description"), anc_number});
    }

    std::sort(fields.begin(), fields.end(), [](const FieldInfo& a, const FieldInfo& b) {
        if (a.anc_number != b.anc_number) return a.anc_number < b.anc_number;
        return a.id < b.id;
    });
    return fields;
}

static std::vector<FieldInfo> requested_field_infos(
    bcf_hdr_t* hdr,
    const std::vector<std::string>& requested_fields
) {
    std::vector<FieldInfo> out;
    for (const std::string& field : requested_fields) {
        std::string description;
        for (int i = 0; i < hdr->nhrec; ++i) {
            bcf_hrec_t* hrec = hdr->hrec[i];
            if (!hrec || !hrec->key || std::strcmp(hrec->key, "FORMAT") != 0) continue;
            if (field == hrec_value(hrec, "ID")) {
                description = hrec_value(hrec, "Description");
                break;
            }
        }
        int anc_number = 0;
        parse_anc_number(field, anc_number);
        out.push_back(FieldInfo{field, description, anc_number});
    }
    return out;
}

static std::vector<size_t> register_samples(Totals& totals, bcf_hdr_t* hdr, size_t n_anc) {
    int n_samples = bcf_hdr_nsamples(hdr);
    std::vector<size_t> local_to_global(n_samples);
    for (int i = 0; i < n_samples; ++i) {
        std::string sample = hdr->samples[i];
        auto it = totals.sample_to_index.find(sample);
        if (it == totals.sample_to_index.end()) {
            size_t index = totals.samples.size();
            totals.sample_to_index.emplace(sample, index);
            totals.samples.push_back(sample);
            totals.counts.push_back(std::vector<double>(n_anc, 0.0));
            totals.denominators.push_back(0.0);
            local_to_global[i] = index;
        } else {
            local_to_global[i] = it->second;
            if (totals.counts[it->second].size() != n_anc) {
                die("internal error: sample count vector has wrong ancestry size");
            }
        }
    }
    return local_to_global;
}

static bool is_bad_float(float value) {
    return bcf_float_is_missing(value) ||
           bcf_float_is_vector_end(value) ||
           !std::isfinite(static_cast<double>(value)) ||
           value < 0.0f;
}

static std::vector<std::string> field_ids(const std::vector<FieldInfo>& infos) {
    std::vector<std::string> ids;
    ids.reserve(infos.size());
    for (const FieldInfo& info : infos) ids.push_back(info.id);
    return ids;
}

static void process_flare_file(
    const std::string& path,
    size_t file_index,
    size_t file_count,
    std::vector<std::string>& active_fields,
    std::vector<std::string>& labels,
    Totals& totals,
    Stats& stats,
    uint64_t progress_every,
    uint64_t thin
) {
    htsFile* fp = hts_open(path.c_str(), "r");
    if (!fp) die("failed to open " + path);

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        die(path + ": failed to read VCF header");
    }

    if (!has_format_field(hdr, "AN1") || !has_format_field(hdr, "AN2")) {
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        die(path + ": FLARE mode requires FORMAT/AN1 and FORMAT/AN2");
    }

    std::vector<std::string> flare_fields = {"AN1", "AN2"};
    if (active_fields.empty()) {
        active_fields = flare_fields;
        std::vector<std::string> parsed_labels = parse_flare_ancestry_labels(hdr);
        if (!parsed_labels.empty()) labels = parsed_labels;
    } else if (active_fields != flare_fields) {
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        die(path + ": input ancestry fields differ from previous files");
    } else if (labels.empty()) {
        std::vector<std::string> parsed_labels = parse_flare_ancestry_labels(hdr);
        if (!parsed_labels.empty()) labels = parsed_labels;
    }

    int n_samples = bcf_hdr_nsamples(hdr);
    if (n_samples <= 0) {
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        die(path + ": VCF header has no sample columns");
    }

    ensure_ancestry_size(totals, labels, labels.size());
    std::vector<size_t> local_to_global = register_samples(totals, hdr, labels.size());

    int32_t* an1 = nullptr;
    int32_t* an2 = nullptr;
    int n_an1 = 0;
    int n_an2 = 0;
    bcf1_t* rec = bcf_init();
    if (!rec) {
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        die("failed to allocate VCF record");
    }

    uint64_t records = 0;
    uint64_t used_records = 0;
    uint64_t skipped_format = 0;
    uint64_t skipped_sample_values = 0;
    ProgressReporter progress(path, file_index, file_count, progress_every);

    while (bcf_read(fp, hdr, rec) == 0) {
        ++records;
        if (thin > 1 && ((records - 1) % thin) != 0) {
            progress.maybe_report(records, used_records);
            continue;
        }

        int n1 = bcf_get_format_int32(hdr, rec, "AN1", &an1, &n_an1);
        int n2 = bcf_get_format_int32(hdr, rec, "AN2", &an2, &n_an2);
        if (n1 != n_samples || n2 != n_samples) {
            ++skipped_format;
            progress.maybe_report(records, used_records);
            continue;
        }

        ++used_records;
        for (int sample_i = 0; sample_i < n_samples; ++sample_i) {
            int32_t a1 = an1[sample_i];
            int32_t a2 = an2[sample_i];
            if (is_missing_int32(a1) || is_missing_int32(a2)) {
                ++skipped_sample_values;
                continue;
            }

            size_t idx1 = static_cast<size_t>(a1);
            size_t idx2 = static_cast<size_t>(a2);
            size_t needed = std::max(idx1, idx2) + 1;
            if (labels.size() < needed) {
                ensure_ancestry_size(totals, labels, needed);
            }

            size_t global_i = local_to_global[sample_i];
            totals.counts[global_i][idx1] += 1.0;
            totals.counts[global_i][idx2] += 1.0;
            totals.denominators[global_i] += 2.0;
        }

        progress.maybe_report(records, used_records);
    }

    std::free(an1);
    std::free(an2);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    if (hts_close(fp) != 0) die(path + ": error while closing input");

    stats.records += records;
    stats.used_records += used_records;
    stats.skipped_format_rows += skipped_format;
    stats.skipped_sample_values += skipped_sample_values;

    std::cerr << "[done] " << path << ": FLARE AN1/AN2, used " << used_records
              << "/" << records << " records, " << n_samples
              << " samples, skipped " << skipped_format
              << " rows without scalar AN1/AN2\n";
}

static void process_file(
    const std::string& path,
    size_t file_index,
    size_t file_count,
    const std::vector<std::string>& requested_fields,
    std::vector<std::string>& active_fields,
    std::vector<std::string>& labels,
    Totals& totals,
    Stats& stats,
    uint64_t progress_every,
    uint64_t thin
) {
    htsFile* fp = hts_open(path.c_str(), "r");
    if (!fp) die("failed to open " + path);

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        die(path + ": failed to read VCF header");
    }

    std::vector<FieldInfo> infos = requested_fields.empty()
        ? detect_anc_fields(hdr)
        : requested_field_infos(hdr, requested_fields);
    if (infos.empty()) {
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        die(path + ": no FORMAT fields named ANC1, ANC2, ... were found; use --anc-fields if needed");
    }

    std::vector<std::string> ids = field_ids(infos);
    if (active_fields.empty()) {
        active_fields = ids;
        labels.reserve(infos.size());
        for (const FieldInfo& info : infos) {
            labels.push_back(label_from_description(info.id, info.description));
        }
    } else if (active_fields != ids) {
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        std::string msg = path + ": ancestry fields differ from previous files";
        die(msg);
    }

    int n_samples = bcf_hdr_nsamples(hdr);
    if (n_samples <= 0) {
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        die(path + ": VCF header has no sample columns");
    }

    const size_t n_anc = active_fields.size();
    std::vector<size_t> local_to_global = register_samples(totals, hdr, n_anc);
    std::vector<float*> buffers(n_anc, nullptr);
    std::vector<int> capacities(n_anc, 0);
    std::vector<int> widths(n_anc, 0);

    bcf1_t* rec = bcf_init();
    if (!rec) {
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        die("failed to allocate VCF record");
    }

    uint64_t records = 0;
    uint64_t used_records = 0;
    uint64_t skipped_format = 0;
    uint64_t skipped_sample_values = 0;
    ProgressReporter progress(path, file_index, file_count, progress_every);

    while (bcf_read(fp, hdr, rec) == 0) {
        ++records;
        if (thin > 1 && ((records - 1) % thin) != 0) {
            progress.maybe_report(records, used_records);
            continue;
        }
        ++used_records;

        bool missing_format = false;
        for (size_t a = 0; a < n_anc; ++a) {
            int n = bcf_get_format_float(
                hdr,
                rec,
                active_fields[a].c_str(),
                &buffers[a],
                &capacities[a]
            );
            if (n < n_samples || n % n_samples != 0) {
                missing_format = true;
                break;
            }
            widths[a] = n / n_samples;
        }
        if (missing_format) {
            ++skipped_format;
            progress.maybe_report(records, used_records);
            continue;
        }

        for (int sample_i = 0; sample_i < n_samples; ++sample_i) {
            double denom = 0.0;
            bool ok = true;
            for (size_t a = 0; a < n_anc; ++a) {
                float value = buffers[a][sample_i * widths[a]];
                if (is_bad_float(value)) {
                    ok = false;
                    break;
                }
                denom += static_cast<double>(value);
            }
            if (!ok || denom <= 0.0) {
                ++skipped_sample_values;
                continue;
            }

            size_t global_i = local_to_global[sample_i];
            for (size_t a = 0; a < n_anc; ++a) {
                totals.counts[global_i][a] += static_cast<double>(buffers[a][sample_i * widths[a]]);
            }
            totals.denominators[global_i] += denom;
        }

        progress.maybe_report(records, used_records);
    }

    for (float* buffer : buffers) std::free(buffer);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    if (hts_close(fp) != 0) die(path + ": error while closing input");

    stats.records += records;
    stats.used_records += used_records;
    stats.skipped_format_rows += skipped_format;
    stats.skipped_sample_values += skipped_sample_values;

    std::cerr << "[done] " << path << ": used " << used_records << "/" << records
              << " records, " << n_samples << " samples, skipped " << skipped_format
              << " rows without all ANC fields\n";
}

static void write_output(
    const std::string& output,
    const std::vector<std::string>& labels,
    const Totals& totals
) {
    std::ofstream file_out;
    std::ostream* out = &std::cout;
    if (!output.empty()) {
        file_out.open(output);
        if (!file_out) die("failed to open output " + output + ": " + std::strerror(errno));
        out = &file_out;
    }

    *out << std::setprecision(10) << std::defaultfloat;
    *out << "sample_id\tn_called_alleles";
    for (const std::string& label : labels) *out << '\t' << label << "_prop";
    for (const std::string& label : labels) *out << '\t' << label << "_count";
    *out << '\n';

    std::vector<size_t> order(totals.samples.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return totals.samples[a] < totals.samples[b];
    });

    for (size_t idx : order) {
        double denom = totals.denominators[idx];
        *out << totals.samples[idx] << '\t' << denom;
        for (double count : totals.counts[idx]) {
            double prop = denom > 0.0
                ? count / denom
                : std::numeric_limits<double>::quiet_NaN();
            *out << '\t' << prop;
        }
        for (double count : totals.counts[idx]) *out << '\t' << count;
        *out << '\n';
    }
}

static bool should_read_as_flare(
    const std::string& path,
    const std::vector<std::string>& requested_fields
) {
    if (!requested_fields.empty()) return false;

    htsFile* fp = hts_open(path.c_str(), "r");
    if (!fp) die("failed to open " + path);
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        die(path + ": failed to read VCF header");
    }

    bool is_flare = has_format_field(hdr, "AN1") && has_format_field(hdr, "AN2");
    bcf_hdr_destroy(hdr);
    if (hts_close(fp) != 0) die(path + ": error while closing input");
    return is_flare;
}

int main(int argc, char** argv) {
    try {
#ifdef HAVE_HTSLIB_LOG
        hts_set_log_level(HTS_LOG_ERROR);
#endif
        Args args = parse_args(argc, argv);
        std::vector<std::string> inputs = expand_inputs(args);
        std::vector<std::string> missing;
        for (const std::string& input : inputs) {
            if (!file_exists(input)) missing.push_back(input);
        }
        if (!missing.empty()) {
            std::cerr << "ERROR: these input files do not exist:\n";
            for (const std::string& path : missing) std::cerr << "  " << path << '\n';
            return 2;
        }

        Totals totals;
        Stats stats;
        std::vector<std::string> active_fields;
        std::vector<std::string> labels;

        for (size_t i = 0; i < inputs.size(); ++i) {
            const std::string& input = inputs[i];
            if (should_read_as_flare(input, args.requested_fields)) {
                process_flare_file(
                    input,
                    i + 1,
                    inputs.size(),
                    active_fields,
                    labels,
                    totals,
                    stats,
                    args.progress_every,
                    args.thin
                );
            } else {
                process_file(
                    input,
                    i + 1,
                    inputs.size(),
                    args.requested_fields,
                    active_fields,
                    labels,
                    totals,
                    stats,
                    args.progress_every,
                    args.thin
                );
            }
        }

        write_output(args.output, labels, totals);
        std::cerr << "[summary] used " << stats.used_records << "/" << stats.records
                  << " records across " << inputs.size() << " files; skipped "
                  << stats.skipped_format_rows << " rows without ancestry fields and "
                  << stats.skipped_sample_values << " sample-site values\n";
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
