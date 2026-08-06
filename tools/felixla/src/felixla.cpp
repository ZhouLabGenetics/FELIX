#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <unistd.h>

int felixla_pack_main(int argc, char** argv);
int felixla_rfmix_main(int argc, char** argv);
int felixla_dosage_main(int argc, char** argv);
int felixla_to_vcf_main(int argc, char** argv);
int felixla_extract_main(int argc, char** argv);
int felixla_query_main(int argc, char** argv);
int felixla_admixture_main(int argc, char** argv);
int felixla_compare_main(int argc, char** argv);

namespace {

using ToolMain = int (*)(int, char**);

[[noreturn]] void die(const std::string& message) {
    std::cerr << "ERROR: " << message << '\n';
    std::exit(1);
}

void usage(std::ostream& out) {
    out <<
R"(Usage:
  felixla [input flags] [action flags] [output/filter flags]

PLINK-style examples:
  felixla --phase-vcf genotype.phased.vcf.gz \
          --flare-vcf flare.anc.vcf.gz \
          --n-ancestries 5 \
          --make-felixla \
          --out hybrid/chr22

  felixla --felixla hybrid/chr22 \
          --query chr22:160500 \
          --ref A \
          --alt G \
          --nonzero-only

  felixla --felixla hybrid/chr22 --export vcf --out hybrid/chr22.roundtrip

Primary input flags:
  --phase-vcf PATH              Phased diploid genotype VCF/BCF.
  --flare-vcf PATH              FLARE local ancestry VCF/BCF with AN1/AN2; samples matched by ID.
  --rfmix-msp PATH              RFMix MSP file.
  --tractor-dosage-vcf PATH     TRACTOR dosage VCF/BCF with ANC#/DS# fields.
  --felixla PREFIX              Existing FELIXla/tractor_hybrid prefix.
  --vcf PATH                    VCF input for --admixture; may be repeated.

Action flags:
  --make-felixla                Create a FELIXla prefix from inputs, or extract from --felixla.
  --export vcf                  Export --felixla prefix to split-biallelic phased VCF.
  --query [CHR:POS]             Query DSALL and ancestry-specific DS1..DSk.
  --admixture                   Estimate per-sample global ancestry proportions.
  --compare-vcfs A B            Compare two VCFs.
  --shapeit-args                Build converter argument rows from SHAPEIT/GLIMPSE chunks.
  --recommend-mac-threshold     Print ceil(n_samples / 32).

Common output and parameter flags:
  --out PATH                    Output prefix or report path, depending on action.
  --n-ancestries INT            Number of ancestry labels.
  --mac-threshold INT           Sparse/dense MAC threshold. Default: auto.
  --n-samples INT               Sample count for threshold helpers and template-only jobs.
  --keep FILE                   Sample IDs to retain, one ID per line.
  --extract FILE                PVAR/VCF allele list to retain by CHROM/POS/REF/ALT.
  --region CHR:START-END        Region for conversion/extraction.
  --chr CHR --from-bp N --to-bp N
                               PLINK-like region form for extraction.
  --bp N                        Position for --query with --chr.
  --ref ALLELE --alt ALLELE     Disambiguate split-biallelic query target.
  --global-index N              Query by split variant global index.
  --nonzero-only                For --query, emit only samples with nonzero DSALL.

Compatibility mode:
  felixla from-flare ...        Equivalent to old flare_subset_to_tractor_hybrid arguments.
  felixla from-rfmix ...        Equivalent to old rfmix_msp_to_tractor_hybrid arguments.
  felixla from-tractor-dosage ...
  felixla to-vcf ...
  felixla extract ...
  felixla query ...
  felixla admixture ...
  felixla compare-vcfs ...
  felixla shapeit-args ...
)";
}

void print_version() {
    std::cout
        << "FELIXla CLI v0 for FELIX "
        << "(Full-cohort Efficient Local ancestry-Integrated miXed-model framework), "
        << "tractor_hybrid-compatible layout format_version=1\n";
}

bool starts_with_dash(const std::string& value) {
    return !value.empty() && value[0] == '-';
}

std::string require_value(int& i, int argc, char** argv, const std::string& flag) {
    if (i + 1 >= argc) die(flag + " requires a value");
    return argv[++i];
}

std::vector<char*> make_argv(std::vector<std::string>& args) {
    std::vector<char*> argv;
    argv.reserve(args.size() + 1);
    for (std::string& arg : args) argv.push_back(arg.data());
    argv.push_back(nullptr);
    return argv;
}

int run_tool(ToolMain tool, std::vector<std::string> args) {
    std::vector<char*> argv = make_argv(args);
    return tool(static_cast<int>(args.size()), argv.data());
}

int run_tool_stdout(const std::string& out_path, ToolMain tool, std::vector<std::string> args) {
    if (out_path.empty()) return run_tool(tool, std::move(args));

    std::fflush(stdout);
    int saved_stdout = dup(STDOUT_FILENO);
    if (saved_stdout < 0) {
        die(std::string("failed to save stdout: ") + std::strerror(errno));
    }
    if (!freopen(out_path.c_str(), "w", stdout)) {
        int open_errno = errno;
        close(saved_stdout);
        die("cannot open output file " + out_path + ": " + std::strerror(open_errno));
    }

    int rc = run_tool(tool, std::move(args));
    std::fflush(stdout);
    if (dup2(saved_stdout, STDOUT_FILENO) < 0) {
        close(saved_stdout);
        die(std::string("failed to restore stdout: ") + std::strerror(errno));
    }
    close(saved_stdout);
    clearerr(stdout);
    return rc;
}

std::string vcf_export_path(const std::string& out) {
    if (out.size() >= 4 && out.compare(out.size() - 4, 4, ".vcf") == 0) return out;
    if (out.size() >= 7 && out.compare(out.size() - 7, 7, ".vcf.gz") == 0) return out;
    if (out.size() >= 8 && out.compare(out.size() - 8, 8, ".vcf.bgz") == 0) return out;
    return out + ".vcf.gz";
}

bool is_nonnegative_integer(const std::string& text) {
    if (text.empty()) return false;
    for (char c : text) {
        if (c < '0' || c > '9') return false;
    }
    return true;
}

int parse_positive_int(const std::string& text, const std::string& label) {
    if (!is_nonnegative_integer(text)) die(label + " must be a non-negative integer");
    long long value = std::strtoll(text.c_str(), nullptr, 10);
    if (value <= 0) die(label + " must be positive");
    if (value > 2147483647LL) die(label + " is too large");
    return static_cast<int>(value);
}

std::string compute_mac_threshold(const std::string& n_samples) {
    int samples = parse_positive_int(n_samples, "n_samples");
    return std::to_string((samples + 31) / 32);
}

int run_mac_threshold(const std::vector<std::string>& args) {
    std::string n_samples;
    if (args.size() == 2 && is_nonnegative_integer(args[1])) {
        n_samples = args[1];
    } else {
        for (size_t i = 1; i < args.size(); ++i) {
            const std::string& arg = args[i];
            if (arg == "-h" || arg == "--help") {
                std::cout <<
R"(Usage:
  felixla --recommend-mac-threshold --n-samples N_SAMPLES
  felixla mac-threshold N_SAMPLES

Storage break-even threshold:
  mac_threshold = ceil(n_samples / 32)
)";
                return 0;
            }
            if (arg == "-n" || arg == "--n-samples") {
                if (i + 1 >= args.size()) die(arg + " requires a value");
                n_samples = args[++i];
            } else {
                die("unknown mac-threshold argument: " + arg);
            }
        }
    }
    std::string threshold = compute_mac_threshold(n_samples);
    std::cout << "n_samples\t" << n_samples << '\n';
    std::cout << "recommended_mac_threshold\t" << threshold << '\n';
    return 0;
}

std::string strip_chr_prefix(const std::string& chrom) {
    if (chrom.rfind("chr", 0) == 0) return chrom.substr(3);
    return chrom;
}

std::string style_chrom(const std::string& chrom, const std::string& chrom_style) {
    std::string base = strip_chr_prefix(chrom);
    if (chrom_style == "chr") return "chr" + base;
    if (chrom_style == "nochr") return base;
    return chrom;
}

std::string style_region(const std::string& region, const std::string& chrom_style) {
    size_t colon = region.find(':');
    if (colon == std::string::npos) die("region lacks chromosome prefix: " + region);
    return style_chrom(region.substr(0, colon), chrom_style) + ":" + region.substr(colon + 1);
}

void replace_all(std::string& value, const std::string& from, const std::string& to) {
    size_t pos = 0;
    while ((pos = value.find(from, pos)) != std::string::npos) {
        value.replace(pos, from.size(), to);
        pos += to.size();
    }
}

std::string zero_pad(int value, int width) {
    std::ostringstream out;
    out << std::setw(width) << std::setfill('0') << value;
    return out.str();
}

std::string render_template(
    std::string value,
    const std::string& chrom,
    const std::string& chunk,
    const std::string& chunk0,
    const std::string& global_chunk,
    const std::string& global_chunk0,
    const std::string& shapeit_chunk,
    const std::string& region,
    const std::string& core_region,
    const std::string& buffered_region
) {
    replace_all(value, "{chrom}", chrom);
    replace_all(value, "{chunk}", chunk);
    replace_all(value, "{chunk0}", chunk0);
    replace_all(value, "{global_chunk}", global_chunk);
    replace_all(value, "{global_chunk0}", global_chunk0);
    replace_all(value, "{shapeit_chunk}", shapeit_chunk);
    replace_all(value, "{region}", region);
    replace_all(value, "{core_region}", core_region);
    replace_all(value, "{buffered_region}", buffered_region);
    return value;
}

bool file_exists(const std::string& path) {
    std::ifstream in(path);
    return static_cast<bool>(in);
}

int run_shapeit_args(const std::vector<std::string>& raw_args, const std::string& plink_out_path = "") {
    std::string phase_template = "phase/{chrom}.phased.vcf.gz";
    std::string flare_template = "flare/{chrom}.flare.vcf.gz";
    std::string out_prefix_template = "hybrid/{chrom}.shapeit4cM.chunk{chunk0}";
    std::string region_column = "core";
    std::string chrom_style = "keep";
    std::string n_ancestries;
    std::string mac_threshold;
    std::string chunks_dir;
    std::string out_path = plink_out_path;
    int chunk_width = 4;
    bool header = false;
    std::vector<std::string> chunks;

    for (size_t i = 1; i < raw_args.size(); ++i) {
        const std::string& arg = raw_args[i];
        auto value = [&](const std::string& flag) -> std::string {
            if (i + 1 >= raw_args.size()) die(flag + " needs a value");
            return raw_args[++i];
        };

        if (arg == "--chunks") chunks.push_back(value(arg));
        else if (arg == "--chunks-dir") chunks_dir = value(arg);
        else if (arg == "--phase-template") phase_template = value(arg);
        else if (arg == "--flare-template") flare_template = value(arg);
        else if (arg == "--out-prefix-template") out_prefix_template = value(arg);
        else if (arg == "--n-ancestries") n_ancestries = value(arg);
        else if (arg == "--mac-threshold") mac_threshold = value(arg);
        else if (arg == "--region-column") region_column = value(arg);
        else if (arg == "--chrom-style") chrom_style = value(arg);
        else if (arg == "--chunk-width") chunk_width = parse_positive_int(value(arg), "--chunk-width");
        else if (arg == "--out") out_path = value(arg);
        else if (arg == "--header") header = true;
        else if (arg == "-h" || arg == "--help") {
            std::cout <<
R"(Usage:
  felixla shapeit-args --chunks FILE [--chunks FILE ...] [options]
  felixla shapeit-args --chunks-dir DIR [options]

Convert SHAPEIT5/GLIMPSE-style chunk files to FELIXla converter argument rows.

Options:
  --chunks FILE              Chunk file to read. May be repeated.
  --chunks-dir DIR           Directory with chunks_chr1.txt ... chunks_chr22.txt.
  --phase-template STR       Source phased VCF template. Default: phase/{chrom}.phased.vcf.gz
  --flare-template STR       Source FLARE VCF template. Default: flare/{chrom}.flare.vcf.gz
  --out-prefix-template STR  Output prefix template. Default: hybrid/{chrom}.shapeit4cM.chunk{chunk0}
  --n-ancestries INT         Required converter n_ancestries argument.
  --mac-threshold INT        Required converter mac_threshold argument.
  --region-column core|buffered
  --chrom-style keep|chr|nochr
  --chunk-width INT          Width for {chunk0} and {global_chunk0}. Default: 4.
  --header                   Emit a header row.
)";
            return 0;
        } else {
            die("unknown shapeit-args argument: " + arg);
        }
    }

    if (region_column != "core" && region_column != "buffered") {
        die("--region-column must be core or buffered");
    }
    if (chrom_style != "keep" && chrom_style != "chr" && chrom_style != "nochr") {
        die("--chrom-style must be keep, chr, or nochr");
    }
    if (n_ancestries.empty()) die("--n-ancestries is required");
    if (mac_threshold.empty()) die("--mac-threshold is required");

    if (!chunks_dir.empty()) {
        for (int chr = 1; chr <= 22; ++chr) {
            std::string chunk_file = chunks_dir;
            if (!chunk_file.empty() && chunk_file.back() != '/') chunk_file += "/";
            chunk_file += "chunks_chr" + std::to_string(chr) + ".txt";
            if (!file_exists(chunk_file)) die("missing chunk file: " + chunk_file);
            chunks.push_back(chunk_file);
        }
    }
    if (chunks.empty()) die("provide --chunks or --chunks-dir");

    std::ofstream out_file;
    std::ostream* out = &std::cout;
    if (!out_path.empty()) {
        out_file.open(out_path);
        if (!out_file) die("cannot open output file: " + out_path);
        out = &out_file;
    }

    if (header) {
        *out << "source_phase_vcf\tsource_flare_vcf\tn_ancestries\tmac_threshold\tout_prefix\tregion\n";
    }

    int global_chunk = 0;
    for (const std::string& chunk_file : chunks) {
        std::ifstream in(chunk_file);
        if (!in) die("chunk file not found: " + chunk_file);

        int per_chrom_chunk = 0;
        int line_no = 0;
        std::string line;
        while (std::getline(in, line)) {
            ++line_no;
            std::istringstream iss(line);
            std::string shapeit_chunk;
            std::string chrom;
            std::string buffered_region;
            std::string core_region;
            if (!(iss >> shapeit_chunk)) continue;
            if (!shapeit_chunk.empty() && shapeit_chunk[0] == '#') continue;
            if (!(iss >> chrom >> buffered_region >> core_region)) {
                die(chunk_file + ":" + std::to_string(line_no) + " needs at least four columns");
            }

            ++per_chrom_chunk;
            ++global_chunk;
            std::string chunk = std::to_string(per_chrom_chunk);
            std::string global = std::to_string(global_chunk);
            std::string chunk0 = zero_pad(per_chrom_chunk, chunk_width);
            std::string global0 = zero_pad(global_chunk, chunk_width);
            std::string out_chrom = style_chrom(chrom, chrom_style);
            std::string out_buffered_region = style_region(buffered_region, chrom_style);
            std::string out_core_region = style_region(core_region, chrom_style);
            std::string region = region_column == "buffered" ? out_buffered_region : out_core_region;

            std::string source_phase = render_template(
                phase_template, out_chrom, chunk, chunk0, global, global0, shapeit_chunk,
                region, out_core_region, out_buffered_region
            );
            std::string source_flare = render_template(
                flare_template, out_chrom, chunk, chunk0, global, global0, shapeit_chunk,
                region, out_core_region, out_buffered_region
            );
            std::string out_prefix = render_template(
                out_prefix_template, out_chrom, chunk, chunk0, global, global0, shapeit_chunk,
                region, out_core_region, out_buffered_region
            );

            *out << source_phase << '\t'
                 << source_flare << '\t'
                 << n_ancestries << '\t'
                 << mac_threshold << '\t'
                 << out_prefix << '\t'
                 << region << '\n';
        }
    }

    return 0;
}

struct PlinkArgs {
    std::string phase_vcf;
    std::string flare_vcf;
    std::string rfmix_msp;
    std::string tractor_dosage_vcf;
    std::string felixla_prefix;
    std::string out_path;
    std::string n_ancestries;
    std::string mac_threshold;
    std::string region;
    std::string chr;
    std::string from_bp;
    std::string to_bp;
    std::string bp;
    std::string ref;
    std::string alt;
    std::string global_index;
    std::string export_format;
    std::string query_target;
    std::string n_samples;
    std::string keep_path;
    std::string extract_path;
    std::string max_diffs;
    std::string pattern;
    std::string chroms;
    std::string anc_fields;
    std::string progress_every;
    std::string thin;
    std::string chunks_dir;
    std::string phase_template;
    std::string flare_template;
    std::string out_prefix_template;
    std::string region_column;
    std::string chrom_style;
    std::string chunk_width;
    bool make_felixla = false;
    bool query_action = false;
    bool admixture_action = false;
    bool shapeit_action = false;
    bool compare_action = false;
    bool recommend_action = false;
    bool nonzero_only = false;
    bool split_multiallelic = false;
    bool header = false;
    std::vector<std::string> vcfs;
    std::vector<std::string> compare_paths;
    std::vector<std::string> chunks;
};

PlinkArgs parse_plink_args(int argc, char** argv) {
    PlinkArgs args;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            usage(std::cout);
            std::exit(0);
        } else if (arg == "--version") {
            print_version();
            std::exit(0);
        } else if (arg == "--phase-vcf" || arg == "--genotype-vcf") {
            args.phase_vcf = require_value(i, argc, argv, arg);
        } else if (arg == "--flare-vcf" || arg == "--lai-vcf") {
            args.flare_vcf = require_value(i, argc, argv, arg);
        } else if (arg == "--rfmix-msp") {
            args.rfmix_msp = require_value(i, argc, argv, arg);
        } else if (arg == "--tractor-dosage-vcf") {
            args.tractor_dosage_vcf = require_value(i, argc, argv, arg);
        } else if (arg == "--felixla" || arg == "--tractor-hybrid") {
            args.felixla_prefix = require_value(i, argc, argv, arg);
        } else if (arg == "--vcf") {
            args.vcfs.push_back(require_value(i, argc, argv, arg));
        } else if (arg == "--make-felixla" || arg == "--make-tractor-hybrid") {
            args.make_felixla = true;
        } else if (arg == "--export") {
            args.export_format = require_value(i, argc, argv, arg);
        } else if (arg == "--query") {
            args.query_action = true;
            if (i + 1 < argc && !starts_with_dash(argv[i + 1])) {
                args.query_target = argv[++i];
            }
        } else if (arg == "--admixture") {
            args.admixture_action = true;
        } else if (arg == "--shapeit-args") {
            args.shapeit_action = true;
        } else if (arg == "--compare-vcfs") {
            if (i + 2 >= argc) die("--compare-vcfs requires two VCF paths");
            args.compare_action = true;
            args.compare_paths = {argv[++i], argv[++i]};
        } else if (arg == "--recommend-mac-threshold" || arg == "--choose-mac-threshold") {
            args.recommend_action = true;
        } else if (arg == "--out") {
            args.out_path = require_value(i, argc, argv, arg);
        } else if (arg == "--n-ancestries" || arg == "--n-ancestry" || arg == "--ancestries") {
            args.n_ancestries = require_value(i, argc, argv, arg);
        } else if (arg == "--mac-threshold") {
            args.mac_threshold = require_value(i, argc, argv, arg);
        } else if (arg == "--n-samples") {
            args.n_samples = require_value(i, argc, argv, arg);
        } else if (arg == "--keep") {
            args.keep_path = require_value(i, argc, argv, arg);
        } else if (arg == "--extract") {
            args.extract_path = require_value(i, argc, argv, arg);
        } else if (arg == "--region") {
            args.region = require_value(i, argc, argv, arg);
        } else if (arg == "--chr") {
            args.chr = require_value(i, argc, argv, arg);
        } else if (arg == "--from-bp") {
            args.from_bp = require_value(i, argc, argv, arg);
        } else if (arg == "--to-bp") {
            args.to_bp = require_value(i, argc, argv, arg);
        } else if (arg == "--bp") {
            args.bp = require_value(i, argc, argv, arg);
        } else if (arg == "--ref") {
            args.ref = require_value(i, argc, argv, arg);
        } else if (arg == "--alt") {
            args.alt = require_value(i, argc, argv, arg);
        } else if (arg == "--global-index") {
            args.global_index = require_value(i, argc, argv, arg);
            args.query_action = true;
        } else if (arg == "--nonzero-only") {
            args.nonzero_only = true;
        } else if (arg == "--pattern") {
            args.pattern = require_value(i, argc, argv, arg);
        } else if (arg == "--chroms") {
            args.chroms = require_value(i, argc, argv, arg);
        } else if (arg == "--anc-fields") {
            args.anc_fields = require_value(i, argc, argv, arg);
        } else if (arg == "--progress-every") {
            args.progress_every = require_value(i, argc, argv, arg);
        } else if (arg == "--thin") {
            args.thin = require_value(i, argc, argv, arg);
        } else if (arg == "--split-multiallelic") {
            args.split_multiallelic = true;
        } else if (arg == "--max-diffs") {
            args.max_diffs = require_value(i, argc, argv, arg);
        } else if (arg == "--chunks") {
            args.chunks.push_back(require_value(i, argc, argv, arg));
        } else if (arg == "--chunks-dir") {
            args.chunks_dir = require_value(i, argc, argv, arg);
        } else if (arg == "--phase-template") {
            args.phase_template = require_value(i, argc, argv, arg);
        } else if (arg == "--flare-template") {
            args.flare_template = require_value(i, argc, argv, arg);
        } else if (arg == "--out-prefix-template") {
            args.out_prefix_template = require_value(i, argc, argv, arg);
        } else if (arg == "--region-column") {
            args.region_column = require_value(i, argc, argv, arg);
        } else if (arg == "--chrom-style") {
            args.chrom_style = require_value(i, argc, argv, arg);
        } else if (arg == "--chunk-width") {
            args.chunk_width = require_value(i, argc, argv, arg);
        } else if (arg == "--header") {
            args.header = true;
        } else {
            die("unknown FELIXla flag: " + arg);
        }
    }
    return args;
}

int run_plink_style(int argc, char** argv) {
    PlinkArgs args = parse_plink_args(argc, argv);

    if (!args.keep_path.empty() || !args.extract_path.empty()) {
        if (!args.make_felixla ||
            args.query_action ||
            args.admixture_action ||
            args.shapeit_action ||
            args.compare_action ||
            args.recommend_action ||
            !args.export_format.empty() ||
            args.phase_vcf.empty() ||
            args.flare_vcf.empty()) {
            die("--keep/--extract are currently supported only for --phase-vcf + --flare-vcf --make-felixla");
        }
    }

    if (args.recommend_action) {
        if (args.n_samples.empty()) die("--recommend-mac-threshold requires --n-samples");
        return run_mac_threshold({"mac-threshold", "--n-samples", args.n_samples});
    }

    if (args.shapeit_action) {
        if (args.mac_threshold.empty() && !args.n_samples.empty()) {
            args.mac_threshold = compute_mac_threshold(args.n_samples);
            std::cerr << "FELIXla: --mac-threshold not supplied; using "
                      << args.mac_threshold << " = ceil(" << args.n_samples
                      << " / 32) from --n-samples.\n";
        }
        if (args.mac_threshold.empty()) die("--shapeit-args requires --mac-threshold or --n-samples");

        std::vector<std::string> shapeit_args = {"shapeit-args"};
        for (const std::string& chunk : args.chunks) {
            shapeit_args.push_back("--chunks");
            shapeit_args.push_back(chunk);
        }
        if (!args.chunks_dir.empty()) {
            shapeit_args.push_back("--chunks-dir");
            shapeit_args.push_back(args.chunks_dir);
        }
        if (!args.phase_template.empty()) {
            shapeit_args.push_back("--phase-template");
            shapeit_args.push_back(args.phase_template);
        }
        if (!args.flare_template.empty()) {
            shapeit_args.push_back("--flare-template");
            shapeit_args.push_back(args.flare_template);
        }
        if (!args.out_prefix_template.empty()) {
            shapeit_args.push_back("--out-prefix-template");
            shapeit_args.push_back(args.out_prefix_template);
        }
        if (!args.n_ancestries.empty()) {
            shapeit_args.push_back("--n-ancestries");
            shapeit_args.push_back(args.n_ancestries);
        }
        shapeit_args.push_back("--mac-threshold");
        shapeit_args.push_back(args.mac_threshold);
        if (!args.region_column.empty()) {
            shapeit_args.push_back("--region-column");
            shapeit_args.push_back(args.region_column);
        }
        if (!args.chrom_style.empty()) {
            shapeit_args.push_back("--chrom-style");
            shapeit_args.push_back(args.chrom_style);
        }
        if (!args.chunk_width.empty()) {
            shapeit_args.push_back("--chunk-width");
            shapeit_args.push_back(args.chunk_width);
        }
        if (args.header) shapeit_args.push_back("--header");
        return run_shapeit_args(shapeit_args, args.out_path);
    }

    if (args.compare_action) {
        std::vector<std::string> compare_args = {"compare_vcfs"};
        compare_args.insert(compare_args.end(), args.compare_paths.begin(), args.compare_paths.end());
        if (args.split_multiallelic) compare_args.push_back("--split-multiallelic");
        if (!args.max_diffs.empty()) {
            compare_args.push_back("--max-diffs");
            compare_args.push_back(args.max_diffs);
        }
        return run_tool_stdout(args.out_path, felixla_compare_main, compare_args);
    }

    if (args.admixture_action) {
        std::vector<std::string> admixture_args = {"calc_tractor_admixture"};
        admixture_args.insert(admixture_args.end(), args.vcfs.begin(), args.vcfs.end());
        if (!args.pattern.empty()) {
            admixture_args.push_back("--pattern");
            admixture_args.push_back(args.pattern);
        }
        if (!args.chroms.empty()) {
            admixture_args.push_back("--chroms");
            admixture_args.push_back(args.chroms);
        }
        if (!args.anc_fields.empty()) {
            admixture_args.push_back("--anc-fields");
            admixture_args.push_back(args.anc_fields);
        }
        if (!args.progress_every.empty()) {
            admixture_args.push_back("--progress-every");
            admixture_args.push_back(args.progress_every);
        }
        if (!args.thin.empty()) {
            admixture_args.push_back("--thin");
            admixture_args.push_back(args.thin);
        }
        if (!args.out_path.empty()) {
            admixture_args.push_back("--output");
            admixture_args.push_back(args.out_path);
        }
        return run_tool(felixla_admixture_main, admixture_args);
    }

    if (args.make_felixla) {
        if (args.out_path.empty()) die("--make-felixla requires --out");
        if (!args.phase_vcf.empty() && !args.flare_vcf.empty()) {
            if (args.n_ancestries.empty()) die("from-FLARE writing requires --n-ancestries");
            if (args.mac_threshold.empty()) args.mac_threshold = "auto";
            std::vector<std::string> pack_args = {
                "flare_subset_to_tractor_hybrid",
                args.phase_vcf,
                args.flare_vcf,
                args.n_ancestries,
                args.mac_threshold,
                args.out_path
            };
            if (!args.region.empty()) pack_args.push_back(args.region);
            if (!args.keep_path.empty()) {
                pack_args.push_back("--keep");
                pack_args.push_back(args.keep_path);
            }
            if (!args.extract_path.empty()) {
                pack_args.push_back("--extract");
                pack_args.push_back(args.extract_path);
            }
            return run_tool(felixla_pack_main, pack_args);
        }
        if (!args.phase_vcf.empty() && !args.rfmix_msp.empty()) {
            if (!args.region.empty()) die("--region is not supported for --rfmix-msp conversion");
            if (!args.keep_path.empty() || !args.extract_path.empty()) {
                die("--keep/--extract are currently supported for --phase-vcf + --flare-vcf writing");
            }
            if (args.n_ancestries.empty()) die("from-RFMix writing requires --n-ancestries");
            if (args.mac_threshold.empty()) args.mac_threshold = "auto";
            return run_tool(felixla_rfmix_main, {
                "rfmix_msp_to_tractor_hybrid",
                args.phase_vcf,
                args.rfmix_msp,
                args.n_ancestries,
                args.mac_threshold,
                args.out_path
            });
        }
        if (!args.tractor_dosage_vcf.empty()) {
            if (!args.keep_path.empty() || !args.extract_path.empty()) {
                die("--keep/--extract are currently supported for --phase-vcf + --flare-vcf writing");
            }
            if (args.n_ancestries.empty()) die("from-TRACTOR-dosage writing requires --n-ancestries");
            if (args.mac_threshold.empty()) args.mac_threshold = "auto";
            return run_tool(felixla_dosage_main, {
                "tractor_dosage_vcf_to_hybrid",
                args.tractor_dosage_vcf,
                args.n_ancestries,
                args.mac_threshold,
                args.out_path
            });
        }
        if (!args.felixla_prefix.empty()) {
            if (!args.keep_path.empty() || !args.extract_path.empty()) {
                die("--keep/--extract are currently supported for --phase-vcf + --flare-vcf writing");
            }
            if (!args.region.empty()) {
                return run_tool(felixla_extract_main, {
                    "tractor_hybrid_extract_region",
                    args.felixla_prefix,
                    args.region,
                    args.out_path
                });
            }
            if (args.chr.empty() || args.from_bp.empty() || args.to_bp.empty()) {
                die("extracting from --felixla requires --region or --chr/--from-bp/--to-bp");
            }
            return run_tool(felixla_extract_main, {
                "tractor_hybrid_extract_region",
                args.felixla_prefix,
                args.chr,
                args.from_bp,
                args.to_bp,
                args.out_path
            });
        }
        die("--make-felixla needs one input source: --phase-vcf + --flare-vcf, --phase-vcf + --rfmix-msp, --tractor-dosage-vcf, or --felixla");
    }

    if (!args.export_format.empty()) {
        if (args.export_format != "vcf") die("only --export vcf is currently supported");
        if (args.felixla_prefix.empty()) die("--export vcf requires --felixla PREFIX");
        if (args.out_path.empty()) die("--export vcf requires --out");
        return run_tool(felixla_to_vcf_main, {
            "tractor_hybrid_to_vcf",
            args.felixla_prefix,
            vcf_export_path(args.out_path)
        });
    }

    if (args.query_action) {
        if (args.felixla_prefix.empty()) die("--query requires --felixla PREFIX");
        std::vector<std::string> query_args = {"felixla_query", args.felixla_prefix};
        if (!args.global_index.empty()) {
            query_args.push_back("--global-index");
            query_args.push_back(args.global_index);
        } else if (!args.query_target.empty()) {
            query_args.push_back(args.query_target);
        } else if (!args.chr.empty() && !args.bp.empty()) {
            query_args.push_back(args.chr);
            query_args.push_back(args.bp);
        } else {
            die("--query requires CHR:POS, --chr/--bp, or --global-index");
        }
        if (!args.ref.empty() || !args.alt.empty()) {
            if (args.ref.empty() || args.alt.empty()) die("--ref and --alt must be supplied together");
            query_args.push_back(args.ref);
            query_args.push_back(args.alt);
        }
        if (args.nonzero_only) query_args.push_back("--nonzero-only");
        return run_tool_stdout(args.out_path, felixla_query_main, query_args);
    }

    usage(std::cerr);
    die("no action flag supplied");
}

std::vector<std::string> command_args(const std::string& command, int argc, char** argv) {
    std::vector<std::string> args;
    args.push_back(command);
    for (int i = 2; i < argc; ++i) args.push_back(argv[i]);
    return args;
}

int run_compat_command(int argc, char** argv) {
    std::string cmd = argv[1];
    if (cmd == "-h" || cmd == "--help" || cmd == "help") {
        usage(std::cout);
        return 0;
    }
    if (cmd == "--version" || cmd == "version") {
        print_version();
        return 0;
    }
    if (cmd == "from-flare" || cmd == "build-flare" || cmd == "pack-flare" ||
        cmd == "flare" || cmd == "flare_subset_to_tractor_hybrid") {
        return run_tool(felixla_pack_main, command_args("flare_subset_to_tractor_hybrid", argc, argv));
    }
    if (cmd == "from-rfmix" || cmd == "rfmix" || cmd == "rfmix-msp" ||
        cmd == "rfmix_msp_to_tractor_hybrid") {
        return run_tool(felixla_rfmix_main, command_args("rfmix_msp_to_tractor_hybrid", argc, argv));
    }
    if (cmd == "from-tractor-dosage" || cmd == "from-dosage" ||
        cmd == "tractor-dosage" || cmd == "dosage" ||
        cmd == "tractor_dosage_vcf_to_hybrid") {
        return run_tool(felixla_dosage_main, command_args("tractor_dosage_vcf_to_hybrid", argc, argv));
    }
    if (cmd == "to-vcf" || cmd == "export-vcf" || cmd == "vcf" ||
        cmd == "tractor_hybrid_to_vcf") {
        return run_tool(felixla_to_vcf_main, command_args("tractor_hybrid_to_vcf", argc, argv));
    }
    if (cmd == "extract" || cmd == "extract-region" || cmd == "region" ||
        cmd == "tractor_hybrid_extract_region") {
        return run_tool(felixla_extract_main, command_args("tractor_hybrid_extract_region", argc, argv));
    }
    if (cmd == "query" || cmd == "dosage-query" || cmd == "felixla_query") {
        return run_tool(felixla_query_main, command_args("felixla_query", argc, argv));
    }
    if (cmd == "admixture" || cmd == "calc-admixture" ||
        cmd == "global-ancestry" || cmd == "calc_tractor_admixture") {
        return run_tool(felixla_admixture_main, command_args("calc_tractor_admixture", argc, argv));
    }
    if (cmd == "compare" || cmd == "compare-vcfs" || cmd == "compare_vcfs") {
        return run_tool(felixla_compare_main, command_args("compare_vcfs", argc, argv));
    }
    if (cmd == "shapeit-args" || cmd == "shapeit_chunks_to_tractor_args" ||
        cmd == "shapeit_chunks_to_tractor_args.sh") {
        return run_shapeit_args(command_args("shapeit-args", argc, argv));
    }
    if (cmd == "mac-threshold" || cmd == "choose-mac-threshold") {
        return run_mac_threshold(command_args("mac-threshold", argc, argv));
    }

    usage(std::cerr);
    die("unknown FELIXla command: " + cmd);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 1) {
        usage(std::cerr);
        return 2;
    }
    std::string first = argv[1];
    if (starts_with_dash(first)) return run_plink_style(argc, argv);
    return run_compat_command(argc, argv);
}
