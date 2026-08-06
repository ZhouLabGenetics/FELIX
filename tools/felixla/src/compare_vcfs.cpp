// compare_vcfs.cpp
// designed by Kai, implemented by codex
//
// Compare two phased genotype VCF/BCF files, optionally after logical ALT split.

#include <htslib/hts.h>
#include <htslib/vcf.h>

#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <iostream>
#include <string>
#include <vector>

struct Options {
    const char* lhs_path = nullptr;
    const char* rhs_path = nullptr;
    bool split_multiallelic = false;
    uint64_t max_diffs = 20;
};

struct SplitVariant {
    std::string chr;
    int64_t pos = 0;
    std::string ref;
    std::string alt;
    std::vector<uint8_t> alt_haps;
};

struct VcfStream {
    const char* path = nullptr;
    htsFile* fp = nullptr;
    bcf_hdr_t* hdr = nullptr;
    bcf1_t* rec = nullptr;
    int32_t* gt_arr = nullptr;
    int ngt_arr = 0;
    int n_samples = 0;
    bool split_multiallelic = false;
    std::deque<SplitVariant> pending;
    bool eof = false;
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

static uint64_t parse_u64_arg(const char* value, const char* name) {
    char* end = nullptr;
    unsigned long long parsed = std::strtoull(value, &end, 10);
    if (!end || *end != '\0') die("invalid %s: %s", name, value);
    return static_cast<uint64_t>(parsed);
}

static void print_usage(const char* prog) {
    std::fprintf(
        stderr,
        "Usage:\n"
        "  %s lhs.vcf[.gz|.bcf] rhs.vcf[.gz|.bcf] [options]\n\n"
        "Options:\n"
        "  --split-multiallelic     Compare after logically splitting each ALT allele\n"
        "  --max-diffs N            Stop after reporting N differences (default: 20)\n",
        prog
    );
}

static Options parse_options(int argc, char** argv) {
    if (argc < 3) {
        print_usage(argv[0]);
        std::exit(1);
    }

    Options opt;
    opt.lhs_path = argv[1];
    opt.rhs_path = argv[2];

    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--split-multiallelic") {
            opt.split_multiallelic = true;
        } else if (arg == "--max-diffs") {
            if (i + 1 >= argc) die("--max-diffs requires a value");
            opt.max_diffs = parse_u64_arg(argv[++i], "--max-diffs");
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            die("unknown option: %s", arg.c_str());
        }
    }

    return opt;
}

static void open_stream(VcfStream& stream, const char* path, bool split_multiallelic) {
    stream.path = path;
    stream.split_multiallelic = split_multiallelic;
    stream.fp = bcf_open(path, "r");
    if (!stream.fp) die("cannot open %s", path);

    stream.hdr = bcf_hdr_read(stream.fp);
    if (!stream.hdr) die("cannot read header from %s", path);

    stream.n_samples = bcf_hdr_nsamples(stream.hdr);
    if (stream.n_samples <= 0) die("%s has no samples", path);

    stream.rec = bcf_init();
}

static void close_stream(VcfStream& stream) {
    if (stream.gt_arr) std::free(stream.gt_arr);
    if (stream.rec) bcf_destroy(stream.rec);
    if (stream.hdr) bcf_hdr_destroy(stream.hdr);
    if (stream.fp) bcf_close(stream.fp);
}

static bool check_sample_order(bcf_hdr_t* lhs, bcf_hdr_t* rhs) {
    int n_lhs = bcf_hdr_nsamples(lhs);
    int n_rhs = bcf_hdr_nsamples(rhs);
    if (n_lhs != n_rhs) return false;

    for (int i = 0; i < n_lhs; ++i) {
        if (std::strcmp(lhs->samples[i], rhs->samples[i]) != 0) {
            std::fprintf(stderr, "Sample mismatch at %d: %s vs %s\n", i, lhs->samples[i], rhs->samples[i]);
            return false;
        }
    }

    return true;
}

static void validate_extra_ploidy(const int32_t* sample_gt, int ploidy, const char* path, const char* chr, int64_t pos, int sample_index) {
    for (int p = 2; p < ploidy; ++p) {
        int32_t gt = sample_gt[p];
        if (gt == bcf_int32_vector_end) break;
        die("non-diploid genotype in %s at %s:%lld sample index %d", path, chr, static_cast<long long>(pos), sample_index);
    }
}

static void fill_pending_from_record(VcfStream& stream) {
    bcf_unpack(stream.rec, BCF_UN_STR);
    if (stream.rec->n_allele < 2) return;

    const char* chr = bcf_hdr_id2name(stream.hdr, stream.rec->rid);
    int64_t pos = static_cast<int64_t>(stream.rec->pos) + 1;
    const char* ref = stream.rec->d.allele[0];

    int ngt = bcf_get_genotypes(stream.hdr, stream.rec, &stream.gt_arr, &stream.ngt_arr);
    if (ngt <= 0) die("missing FORMAT/GT in %s at %s:%lld", stream.path, chr, static_cast<long long>(pos));
    if (ngt % stream.n_samples != 0) {
        die("GT field length not divisible by sample count in %s at %s:%lld", stream.path, chr, static_cast<long long>(pos));
    }

    int ploidy = ngt / stream.n_samples;
    if (ploidy < 2) die("expected diploid GT in %s at %s:%lld", stream.path, chr, static_cast<long long>(pos));

    int first_alt = 1;
    int last_alt = stream.split_multiallelic ? stream.rec->n_allele - 1 : 1;
    if (!stream.split_multiallelic && stream.rec->n_allele != 2) {
        die("%s has multi-allelic record at %s:%lld; rerun with --split-multiallelic", stream.path, chr, static_cast<long long>(pos));
    }

    for (int alt_idx = first_alt; alt_idx <= last_alt; ++alt_idx) {
        SplitVariant variant;
        variant.chr = chr;
        variant.pos = pos;
        variant.ref = ref;
        variant.alt = stream.rec->d.allele[alt_idx];
        variant.alt_haps.assign(static_cast<size_t>(2 * stream.n_samples), 0);

        for (int i = 0; i < stream.n_samples; ++i) {
            const int32_t* sample_gt = stream.gt_arr + static_cast<size_t>(i) * ploidy;
            validate_extra_ploidy(sample_gt, ploidy, stream.path, chr, pos, i);

            int32_t g0 = sample_gt[0];
            int32_t g1 = sample_gt[1];
            if (g0 == bcf_int32_vector_end || g1 == bcf_int32_vector_end ||
                bcf_gt_is_missing(g0) || bcf_gt_is_missing(g1)) {
                die("missing genotype in %s at %s:%lld sample index %d", stream.path, chr, static_cast<long long>(pos), i);
            }
            if (!bcf_gt_is_phased(g1)) {
                die("unphased genotype in %s at %s:%lld sample index %d", stream.path, chr, static_cast<long long>(pos), i);
            }

            if (bcf_gt_allele(g0) == alt_idx) variant.alt_haps[2 * i] = 1;
            if (bcf_gt_allele(g1) == alt_idx) variant.alt_haps[2 * i + 1] = 1;
        }

        stream.pending.push_back(std::move(variant));
    }
}

static bool next_variant(VcfStream& stream, SplitVariant& out) {
    while (stream.pending.empty() && !stream.eof) {
        int ret = bcf_read(stream.fp, stream.hdr, stream.rec);
        if (ret < 0) {
            stream.eof = true;
            break;
        }
        fill_pending_from_record(stream);
    }

    if (stream.pending.empty()) return false;

    out = std::move(stream.pending.front());
    stream.pending.pop_front();
    return true;
}

static std::string variant_key(const SplitVariant& v) {
    return v.chr + ":" + std::to_string(v.pos) + ":" + v.ref + ":" + v.alt;
}

static void report_diff(uint64_t& diffs, uint64_t max_diffs, const std::string& message) {
    ++diffs;
    if (diffs <= max_diffs) {
        std::cerr << "DIFF " << diffs << ": " << message << "\n";
    }
}

int main(int argc, char** argv) {
    Options opt = parse_options(argc, argv);

    VcfStream lhs;
    VcfStream rhs;
    open_stream(lhs, opt.lhs_path, opt.split_multiallelic);
    open_stream(rhs, opt.rhs_path, opt.split_multiallelic);

    if (!check_sample_order(lhs.hdr, rhs.hdr)) {
        die("sample IDs must be identical and in the same order");
    }

    uint64_t compared = 0;
    uint64_t diffs = 0;
    SplitVariant left;
    SplitVariant right;
    bool has_left = next_variant(lhs, left);
    bool has_right = next_variant(rhs, right);

    while (has_left || has_right) {
        if (!has_left) {
            report_diff(diffs, opt.max_diffs, std::string("left EOF before right variant ") + variant_key(right));
            has_right = next_variant(rhs, right);
            continue;
        }
        if (!has_right) {
            report_diff(diffs, opt.max_diffs, std::string("right EOF before left variant ") + variant_key(left));
            has_left = next_variant(lhs, left);
            continue;
        }

        ++compared;
        std::string left_key = variant_key(left);
        std::string right_key = variant_key(right);

        if (left_key != right_key) {
            report_diff(diffs, opt.max_diffs, "variant key mismatch: " + left_key + " vs " + right_key);
        } else if (left.alt_haps != right.alt_haps) {
            uint64_t first = 0;
            while (first < left.alt_haps.size() && left.alt_haps[first] == right.alt_haps[first]) ++first;
            report_diff(
                diffs,
                opt.max_diffs,
                "GT mismatch at " + left_key + " hap_id " + std::to_string(first)
            );
        }

        has_left = next_variant(lhs, left);
        has_right = next_variant(rhs, right);
    }

    close_stream(lhs);
    close_stream(rhs);

    std::cout << "Compared split variants: " << compared << "\n";
    std::cout << "Differences:             " << diffs << "\n";

    return diffs == 0 ? 0 : 2;
}
