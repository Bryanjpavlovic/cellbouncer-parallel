#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <zlib.h>

#include "mt_ratio_model.h"
#include "mt_evidence.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <cctype>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr double kTiny = 1e-12;
constexpr int kReconciledPairPanelMinDepth = 10;
constexpr double kReconciledPairHomoplasmyAf = 0.95;

struct Options {
    std::string bam;
    std::string vcf;
    std::string assignments;
    std::string output_prefix;
    std::string empty_barcodes;
    std::string ambient_profile;
    std::string ambient_fraction_file;
    std::string mt_mask_bed;
    std::string library_id;
    std::string site_manifest;
    std::string site_calibration;
    std::string site_calibration_stratum;
    std::string calibration_reference;
    std::string rho_reference;
    std::string rho_mode = "free";
    std::string site_influence_mode = "full";
    std::string likelihood = "beta_binomial";
    std::string assay_mode = "RNA";
    std::string mito_chrom = "chrM";
    std::string barcode_tag = "CB";
    std::string umi_tag = "UB";
    int min_mapq = 30;
    int min_baseq = 20;
    int min_molecules = 20;
    int min_sites = 3;
    int ambient_min_molecules = 10;
    int min_ambient_only_sites = 1;
    int max_iterations = 500;
    int pileup_max_depth = 10000000;
    int threads = 4;
    double error_rate = 0.005;
    double tolerance = 1e-8;
    double overdispersion_initial = 0.02;
    double overdispersion_max = 0.25;
    double pooled_rho = std::numeric_limits<double>::quiet_NaN();
    double rho_prior_strength = 10.0;
    double ambient_qc_max = std::numeric_limits<double>::quiet_NaN();
    int rho_low_information_molecules = 50;
    double single_parent_epsilon = 0.02;
    double profile_grid_step = 0.005;
    bool ambient_none = false;
    bool allow_no_umi = false;
    bool keep_duplicates = false;
    bool allow_ambiguous_alignments = false;
    bool allow_unanchored_ambient = false;
    bool legacy_panel_gt = false;
    bool write_site_counts = true;
    bool write_profile_grid = false;
    bool atac_include_singletons = false;
};

struct Interval {
    int64_t start = 0;
    int64_t end = 0;
};

struct AlleleCount {
    uint64_t ref = 0;
    uint64_t alt = 0;
    uint64_t total() const { return ref + alt; }
};

struct RawFragmentSite {
    int site_index = -1;
    uint8_t allele = 0;
};

struct RawFragmentObservation {
    std::vector<RawFragmentSite> sites;
};

struct MtSite {
    int rid = -1;
    int64_t pos = -1;  // 0-based
    char ref = 'N';
    char alt = 'N';
    std::vector<int8_t> genotype_state;  // -1 unknown/heteroplasmic, 0 ref, 1 alt
};

struct AmbientSite {
    uint64_t ref = 0;
    uint64_t alt = 0;
    double alt_fraction = std::numeric_limits<double>::quiet_NaN();
    bool usable = false;
    uint64_t total() const { return ref + alt; }
};

struct CellInfo {
    std::string barcode;
    std::string identity;
    std::string parent1;
    std::string parent2;
    int parent1_index = -1;
    int parent2_index = -1;
    std::string calibration_true_parent;
    std::string original_identity;
    std::unordered_map<int, AlleleCount> counts;
    std::vector<RawFragmentObservation> fragments;
};

struct ReadStats {
    uint64_t seen = 0;
    uint64_t accepted_for_pileup = 0;
    uint64_t reject_unmapped = 0;
    uint64_t reject_secondary = 0;
    uint64_t reject_supplementary = 0;
    uint64_t reject_qcfail = 0;
    uint64_t reject_duplicate = 0;
    uint64_t reject_mapq = 0;
    uint64_t reject_multimapping = 0;
    uint64_t pileup_bases = 0;
    uint64_t reject_baseq = 0;
    uint64_t reject_missing_barcode = 0;
    uint64_t reject_irrelevant_barcode = 0;
    uint64_t reject_missing_umi = 0;
    uint64_t reject_nonpanel_allele = 0;
    uint64_t accepted_observations = 0;
    uint64_t conflicting_molecules = 0;
    uint64_t atac_reads_considered = 0;
    uint64_t atac_reject_unpaired = 0;
    uint64_t atac_reject_mate_off_mito = 0;
    uint64_t atac_orphan_fragments = 0;
    uint64_t atac_fragments_accepted = 0;
    uint64_t atac_fragments_multisite = 0;
    uint64_t atac_fragment_site_observations = 0;
    uint64_t atac_overlap_agreements = 0;
    uint64_t atac_overlap_conflicts = 0;
};

struct PileupReader {
    samFile* bam = nullptr;
    hts_itr_t* iterator = nullptr;
    int min_mapq = 30;
    bool keep_duplicates = false;
    bool allow_ambiguous_alignments = false;
    ReadStats* stats = nullptr;
};

enum ManifestSiteClass {
    MANIFEST_A_ALT_B_REF = 0,
    MANIFEST_A_REF_B_ALT = 1,
    MANIFEST_AMBIENT_ONLY = 2,
    MANIFEST_LEGACY_RATIO = 3
};

struct ManifestSite {
    int site_index = -1;
    int8_t canonical_parent1_state = -1;
    int8_t canonical_parent2_state = -1;
    ManifestSiteClass site_class = MANIFEST_LEGACY_RATIO;
};

struct PairManifest {
    std::string canonical_parent1;
    std::string canonical_parent2;
    std::unordered_map<int, ManifestSite> sites;
    int ratio_sites_available = 0;
    int ambient_only_sites_available = 0;
};

struct SiteCalibrationRecord {
    double parent1_alt_probability = std::numeric_limits<double>::quiet_NaN();
    double parent2_alt_probability = std::numeric_limits<double>::quiet_NaN();
    bool has_parent1 = false;
    bool has_parent2 = false;
    int parent1_priority = -1;
    int parent2_priority = -1;
    std::string calibration_id;
};

struct SiteCalibrationData {
    std::unordered_map<std::string, std::unordered_map<int, SiteCalibrationRecord>> by_pair;
    std::set<std::string> calibration_ids;
    uint64_t loaded_rows = 0;
};

struct RhoReferenceEntry {
    std::string assay_mode;
    std::string library_id;
    std::string canonical_parent1;
    std::string canonical_parent2;
    double rho = std::numeric_limits<double>::quiet_NaN();
};

struct LocalMolecule {
    int target_index = -1;  // -1 = empty-droplet ambient molecule
    uint8_t allele = 0;     // 0 ref, 1 alt, 2 conflict
};

void usage(FILE* out, int code) {
    std::fprintf(out,
        "mt_fusion_ratio --bam INPUT.bam --vcf MT_PANEL.bcf \\\n"
        "  --assignments DEMUX.assignments --library_id LIB \\\n"
        "  --site_manifest MT_PANEL.site_manifest.tsv --output_prefix PREFIX [OPTIONS]\n"
        "\n"
        "Estimates the two parental mitochondrial fractions in nuclear fusion calls.\n"
        "Production mode uses the library-specific three-class manifest emitted by\n"
        "downsample_vcf_parallel: both parent-discriminating directions plus sites\n"
        "where the two parents agree and another library donor differs (ambient-only).\n"
        "\n"
        "Required:\n"
        "  --bam FILE                    Coordinate-sorted, indexed BAM/CRAM\n"
        "  --vcf FILE                    Dedicated mitochondrial panel BCF\n"
        "  --assignments FILE            demux_parallel .assignments file\n"
        "  --library_id ID               Library key used in the site manifest\n"
        "  --site_manifest FILE          Per-library three-class site manifest\n"
        "  --output_prefix PREFIX        Output prefix\n"
        "\n"
        "MT ambient/background source (choose exactly one):\n"
        "  --empty_barcodes FILE         Same-library empty-droplet barcode list\n"
        "  --ambient_profile FILE        Site-level mt ambient-profile TSV\n"
        "  --ambient_none                No MT ambient/background source\n"
        "  --ambient_qc_max FLOAT        Exclude MT ratio inference above independent anchor-only c [disabled]\n"
        "\n"
        "Likelihood and identifiability:\n"
        "  --likelihood NAME             beta_binomial (default) or binomial\n"
        "  --overdispersion_initial X    Initial beta-binomial rho [0.02]\n"
        "  --overdispersion_max X        Upper bound for fitted rho [0.25]\n"
        "  --min_ambient_only_sites INT  Legacy compatibility option [1]\n"
        "  --allow_unanchored_ambient    Legacy compatibility option\n"
        "  --ambient_fraction_file FILE  Optional precomputed per-cell ambient QC fraction override\n"
        "  --legacy_panel_gt             Diagnostic only; ignore manifest and use GT-different sites\n"
        "  --single_parent_epsilon FLOAT Practical single-parent boundary [0.02]\n"
        "  --write_profile_grid          Write PREFIX.mt_profile.tsv.gz for population modeling\n"
        "  --profile_grid_step FLOAT     Population-profile grid spacing [0.005]\n"
        "  --assay_mode MODE             RNA (default), ATAC, or GENERIC\n"
        "  --site_calibration FILE       Optional site-specific parental ALT probabilities\n"
        "  --site_calibration_stratum S  Calibration stratum key [library_id]\n"
        "  --calibration_reference FILE  Score listed control barcodes against cross-parent pairs\n"
        "  --rho_mode MODE               free (default), fixed, low_information_fixed, or shrink\n"
        "  --pooled_rho FLOAT            Pooled rho used by fixed/shrink modes\n"
        "  --rho_reference FILE          Optional assay/library/pair-specific pooled-rho TSV\n"
        "  --rho_low_information_molecules INT  Fixed-rho cutoff [50]\n"
        "  --rho_prior_strength FLOAT    Shrinkage penalty strength [10]\n"
        "\n"
        "Counting and filters:\n"
        "  --mito_chrom NAME             Mitochondrial contig [chrM]\n"
        "  --barcode_tag TAG             Cell barcode BAM tag [CB]\n"
        "  --umi_tag TAG                 Corrected UMI BAM tag [UB]\n"
        "  --allow_no_umi                Use read name when UB is absent (RNA/GENERIC only)\n"
        "  --atac_include_singletons     In ATAC mode, retain unmatched/unpaired reads as fragments\n"
        "  --min_mapq INT                MAPQ floor [30]\n"
        "  --min_baseq INT               Base-quality floor [20]\n"
        "  --error_rate FLOAT            Parental-state error probability [0.005]\n"
        "  --min_molecules INT           Parent-discriminating molecules required [20]\n"
        "  --min_sites INT               Parent-discriminating sites observed [3]\n"
        "  --ambient_min_molecules INT   Empty-droplet molecules/site required [10]\n"
        "  --mt_mask_bed FILE            Optional chrM-coordinate ambiguity mask\n"
        "  --keep_duplicates             Retain duplicate-flagged reads\n"
        "  --allow_ambiguous_alignments  Retain NH>1, SA, or XA alignments\n"
        "  --max_iterations INT          Coordinate-ascent maximum [500]\n"
        "  --tolerance FLOAT             Convergence tolerance [1e-8]\n"
        "  --pileup_max_depth INT        Per-position pileup cap [10000000]\n"
        "  --threads INT                 HTSlib decode / post-count fit worker threads [4]\n"
        "  --no_site_counts              Suppress per-cell site-count output\n"
        "  --site_influence_mode MODE    full (default) or none; none keeps site rows but skips LOO refits\n"
        "  --help                        Show this help\n"
        "\n"
        "The autosomal numts.bed is not accepted here because it is in nuclear, not\n"
        "chrM, coordinate space. Default mapping filters and the optional chrM mask\n"
        "provide the read-level NUMT protection.\n");
    std::exit(code);
}

std::string trim(std::string value) {
    const char* ws = " \t\r\n";
    const auto first = value.find_first_not_of(ws);
    if (first == std::string::npos) return "";
    const auto last = value.find_last_not_of(ws);
    return value.substr(first, last - first + 1);
}

// Keep mitochondrial sidecar matching backward-compatible with legacy panel
// metadata that used bare numeric library IDs (for example, "19") while
// downstream code uses the canonical "lib19" form.  This intentionally
// mirrors mt_identity_score.cpp.
std::string canonical_library_key(const std::string& raw) {
    std::string s = trim(raw);
    if (s.size() > 3 &&
        (s[0] == 'l' || s[0] == 'L') &&
        (s[1] == 'i' || s[1] == 'I') &&
        (s[2] == 'b' || s[2] == 'B')) {
        bool numeric = true;
        for (size_t i = 3; i < s.size(); ++i) {
            if (s[i] < '0' || s[i] > '9') {
                numeric = false;
                break;
            }
        }
        if (numeric) s = s.substr(3);
    }
    return s;
}

bool same_library_key(const std::string& a, const std::string& b) {
    return a == b || canonical_library_key(a) == canonical_library_key(b);
}

std::vector<std::string> split_ws(const std::string& line) {
    std::stringstream ss(line);
    std::vector<std::string> fields;
    std::string field;
    while (ss >> field) fields.push_back(field);
    return fields;
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    size_t start = 0;
    while (true) {
        const size_t tab = line.find('\t', start);
        if (tab == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }
    return fields;
}

std::vector<std::string> read_text_lines(const std::string& filename) {
    std::vector<std::string> lines;
    if (filename.size() >= 3 && filename.substr(filename.size() - 3) == ".gz") {
        gzFile gz = gzopen(filename.c_str(), "rb");
        if (!gz) throw std::runtime_error("Could not open gzipped text file: " + filename);
        std::vector<char> buffer(1 << 20);
        while (gzgets(gz, buffer.data(), static_cast<int>(buffer.size())) != nullptr) {
            std::string line(buffer.data());
            while (!line.empty() && line.back() != '\n' && !gzeof(gz)) {
                if (gzgets(gz, buffer.data(), static_cast<int>(buffer.size())) == nullptr) break;
                line += buffer.data();
            }
            lines.push_back(trim(line));
        }
        gzclose(gz);
    } else {
        std::ifstream in(filename);
        if (!in) throw std::runtime_error("Could not open text file: " + filename);
        std::string line;
        while (std::getline(in, line)) lines.push_back(trim(line));
    }
    return lines;
}

std::vector<Interval> load_bed_for_contig(const std::string& filename,
                                          const std::string& contig) {
    std::vector<Interval> intervals;
    if (filename.empty()) return intervals;
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Could not open mask BED: " + filename);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::string chrom;
        int64_t start = 0, end = 0;
        if (!(ss >> chrom >> start >> end)) continue;
        if (chrom == contig && start >= 0 && end > start) intervals.push_back({start, end});
    }
    std::sort(intervals.begin(), intervals.end(), [](const Interval& a, const Interval& b) {
        return a.start < b.start || (a.start == b.start && a.end < b.end);
    });
    std::vector<Interval> merged;
    merged.reserve(intervals.size());
    for (const Interval& interval : intervals) {
        if (merged.empty() || interval.start > merged.back().end) {
            merged.push_back(interval);
        } else if (interval.end > merged.back().end) {
            merged.back().end = interval.end;
        }
    }
    return merged;
}

bool masked_position(int64_t pos, const std::vector<Interval>& intervals) {
    auto it = std::upper_bound(intervals.begin(), intervals.end(), pos,
        [](int64_t value, const Interval& interval) { return value < interval.start; });
    if (it == intervals.begin()) return false;
    --it;
    return it->start <= pos && pos < it->end;
}

bool valid_base(char base) {
    base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

int homoplasmic_state(const int32_t* gt, int ploidy) {
    int state = -1;
    bool saw = false;
    for (int k = 0; k < ploidy; ++k) {
        const int32_t value = gt[k];
        if (value == bcf_int32_vector_end) break;
        if (bcf_gt_is_missing(value)) return -1;
        const int allele = bcf_gt_allele(value);
        if (allele < 0 || allele > 1) return -1;
        if (!saw) {
            state = allele;
            saw = true;
        } else if (state != allele) {
            return -1;
        }
    }
    return saw ? state : -1;
}

Options parse_options(int argc, char** argv) {
    Options opt;
    static option long_options[] = {
        {"bam", required_argument, nullptr, 'b'},
        {"vcf", required_argument, nullptr, 'v'},
        {"assignments", required_argument, nullptr, 'a'},
        {"output_prefix", required_argument, nullptr, 'o'},
        {"empty_barcodes", required_argument, nullptr, 1000},
        {"ambient_profile", required_argument, nullptr, 1001},
        {"ambient_none", no_argument, nullptr, 1002},
        {"ambient_fraction_file", required_argument, nullptr, 1003},
        {"mito_chrom", required_argument, nullptr, 1004},
        {"barcode_tag", required_argument, nullptr, 1005},
        {"umi_tag", required_argument, nullptr, 1006},
        {"allow_no_umi", no_argument, nullptr, 1007},
        {"min_mapq", required_argument, nullptr, 1008},
        {"min_baseq", required_argument, nullptr, 1009},
        {"error_rate", required_argument, nullptr, 1010},
        {"min_molecules", required_argument, nullptr, 1011},
        {"min_sites", required_argument, nullptr, 1012},
        {"ambient_min_molecules", required_argument, nullptr, 1013},
        {"mt_mask_bed", required_argument, nullptr, 1014},
        {"keep_duplicates", no_argument, nullptr, 1015},
        {"allow_ambiguous_alignments", no_argument, nullptr, 1016},
        {"max_iterations", required_argument, nullptr, 1017},
        {"tolerance", required_argument, nullptr, 1018},
        {"pileup_max_depth", required_argument, nullptr, 1019},
        {"threads", required_argument, nullptr, 't'},
        {"no_site_counts", no_argument, nullptr, 1020},
        {"library_id", required_argument, nullptr, 1021},
        {"site_manifest", required_argument, nullptr, 1022},
        {"likelihood", required_argument, nullptr, 1023},
        {"overdispersion_initial", required_argument, nullptr, 1024},
        {"overdispersion_max", required_argument, nullptr, 1025},
        {"min_ambient_only_sites", required_argument, nullptr, 1026},
        {"allow_unanchored_ambient", no_argument, nullptr, 1027},
        {"legacy_panel_gt", no_argument, nullptr, 1028},
        {"single_parent_epsilon", required_argument, nullptr, 1029},
        {"write_profile_grid", no_argument, nullptr, 1030},
        {"profile_grid_step", required_argument, nullptr, 1031},
        {"assay_mode", required_argument, nullptr, 1032},
        {"site_calibration", required_argument, nullptr, 1033},
        {"site_calibration_stratum", required_argument, nullptr, 1034},
        {"calibration_reference", required_argument, nullptr, 1035},
        {"rho_mode", required_argument, nullptr, 1036},
        {"pooled_rho", required_argument, nullptr, 1037},
        {"rho_reference", required_argument, nullptr, 1038},
        {"rho_low_information_molecules", required_argument, nullptr, 1039},
        {"rho_prior_strength", required_argument, nullptr, 1040},
        {"atac_include_singletons", no_argument, nullptr, 1041},
        {"site_influence_mode", required_argument, nullptr, 1042},
        {"ambient_qc_max", required_argument, nullptr, 1043},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int ch = 0;
    while ((ch = getopt_long(argc, argv, "b:v:a:o:t:h", long_options, nullptr)) != -1) {
        switch (ch) {
            case 'b': opt.bam = optarg; break;
            case 'v': opt.vcf = optarg; break;
            case 'a': opt.assignments = optarg; break;
            case 'o': opt.output_prefix = optarg; break;
            case 't': opt.threads = std::stoi(optarg); break;
            case 'h': usage(stdout, 0); break;
            case 1000: opt.empty_barcodes = optarg; break;
            case 1001: opt.ambient_profile = optarg; break;
            case 1002: opt.ambient_none = true; break;
            case 1003: opt.ambient_fraction_file = optarg; break;
            case 1004: opt.mito_chrom = optarg; break;
            case 1005: opt.barcode_tag = optarg; break;
            case 1006: opt.umi_tag = optarg; break;
            case 1007: opt.allow_no_umi = true; break;
            case 1008: opt.min_mapq = std::stoi(optarg); break;
            case 1009: opt.min_baseq = std::stoi(optarg); break;
            case 1010: opt.error_rate = std::stod(optarg); break;
            case 1011: opt.min_molecules = std::stoi(optarg); break;
            case 1012: opt.min_sites = std::stoi(optarg); break;
            case 1013: opt.ambient_min_molecules = std::stoi(optarg); break;
            case 1014: opt.mt_mask_bed = optarg; break;
            case 1015: opt.keep_duplicates = true; break;
            case 1016: opt.allow_ambiguous_alignments = true; break;
            case 1017: opt.max_iterations = std::stoi(optarg); break;
            case 1018: opt.tolerance = std::stod(optarg); break;
            case 1019: opt.pileup_max_depth = std::stoi(optarg); break;
            case 1020: opt.write_site_counts = false; break;
            case 1021: opt.library_id = optarg; break;
            case 1022: opt.site_manifest = optarg; break;
            case 1023: opt.likelihood = optarg; break;
            case 1024: opt.overdispersion_initial = std::stod(optarg); break;
            case 1025: opt.overdispersion_max = std::stod(optarg); break;
            case 1026: opt.min_ambient_only_sites = std::stoi(optarg); break;
            case 1027: opt.allow_unanchored_ambient = true; break;
            case 1028: opt.legacy_panel_gt = true; break;
            case 1029: opt.single_parent_epsilon = std::stod(optarg); break;
            case 1030: opt.write_profile_grid = true; break;
            case 1031: opt.profile_grid_step = std::stod(optarg); break;
            case 1032: opt.assay_mode = optarg; break;
            case 1033: opt.site_calibration = optarg; break;
            case 1034: opt.site_calibration_stratum = optarg; break;
            case 1035: opt.calibration_reference = optarg; break;
            case 1036: opt.rho_mode = optarg; break;
            case 1037: opt.pooled_rho = std::stod(optarg); break;
            case 1038: opt.rho_reference = optarg; break;
            case 1039: opt.rho_low_information_molecules = std::stoi(optarg); break;
            case 1040: opt.rho_prior_strength = std::stod(optarg); break;
            case 1041: opt.atac_include_singletons = true; break;
            case 1042: opt.site_influence_mode = optarg; break;
            case 1043: opt.ambient_qc_max = std::stod(optarg); break;
            default: usage(stderr, 2);
        }
    }

    if (opt.bam.empty() || opt.vcf.empty() || opt.assignments.empty() ||
        opt.output_prefix.empty()) usage(stderr, 2);
    if (!opt.legacy_panel_gt && (opt.library_id.empty() || opt.site_manifest.empty())) {
        std::fprintf(stderr,
            "ERROR: production mode requires both --library_id and --site_manifest; "
            "--legacy_panel_gt is diagnostic only\n");
        std::exit(2);
    }
    if (opt.legacy_panel_gt && (!opt.library_id.empty() || !opt.site_manifest.empty())) {
        std::fprintf(stderr,
            "ERROR: do not combine --legacy_panel_gt with --library_id/--site_manifest\n");
        std::exit(2);
    }
    const int ambient_choices = (!opt.empty_barcodes.empty() ? 1 : 0) +
                                (!opt.ambient_profile.empty() ? 1 : 0) +
                                (opt.ambient_none ? 1 : 0);
    if (ambient_choices != 1) {
        std::fprintf(stderr,
            "ERROR: choose exactly one of --empty_barcodes, --ambient_profile, or --ambient_none\n");
        std::exit(2);
    }
    if (opt.ambient_none && !opt.ambient_fraction_file.empty()) {
        std::fprintf(stderr,
            "ERROR: --ambient_fraction_file cannot be used with --ambient_none\n");
        std::exit(2);
    }
    if (std::isfinite(opt.ambient_qc_max) && opt.ambient_none) {
        std::fprintf(stderr,
            "ERROR: --ambient_qc_max requires --empty_barcodes or --ambient_profile\n");
        std::exit(2);
    }
    if (opt.assay_mode == "ATAC" && opt.allow_no_umi) {
        std::fprintf(stderr,
            "ERROR: --allow_no_umi is an RNA/GENERIC fallback; ATAC mode uses paired fragments directly\n");
        std::exit(2);
    }
    if (opt.rho_mode != "free" && opt.rho_mode != "fixed" &&
        opt.rho_mode != "low_information_fixed" && opt.rho_mode != "shrink") {
        std::fprintf(stderr, "ERROR: --rho_mode must be free, fixed, low_information_fixed, or shrink\n");
        std::exit(2);
    }
    if (opt.assay_mode == "ATAC" && opt.rho_mode != "free") {
        std::fprintf(stderr,
            "ERROR: rho modes apply to RNA/GENERIC beta-binomial site counts; formal ATAC uses fragment likelihoods\n");
        std::exit(2);
    }
    if (opt.site_influence_mode != "full" && opt.site_influence_mode != "none") {
        std::fprintf(stderr, "ERROR: --site_influence_mode must be full or none\n");
        std::exit(2);
    }
    if (opt.barcode_tag.size() != 2 || opt.umi_tag.size() != 2 ||
        opt.min_mapq < 0 || opt.min_baseq < 0 || opt.min_molecules < 0 ||
        opt.min_sites < 1 || opt.ambient_min_molecules < 1 ||
        opt.min_ambient_only_sites < 0 ||
        opt.max_iterations < 1 || opt.pileup_max_depth < 1 || opt.threads < 1 ||
        opt.error_rate <= 0.0 || opt.error_rate >= 0.5 || opt.tolerance <= 0.0 ||
        opt.overdispersion_initial < 0.0 ||
        opt.overdispersion_max <= 0.0 || opt.overdispersion_max >= 1.0 ||
        opt.overdispersion_initial > opt.overdispersion_max ||
        (std::isfinite(opt.pooled_rho) && (opt.pooled_rho < 0.0 || opt.pooled_rho > opt.overdispersion_max)) ||
        opt.rho_low_information_molecules < 1 || opt.rho_prior_strength < 0.0 ||
        (std::isfinite(opt.ambient_qc_max) &&
         (opt.ambient_qc_max < 0.0 || opt.ambient_qc_max > 0.99)) ||
        opt.single_parent_epsilon <= 0.0 || opt.single_parent_epsilon >= 0.5 ||
        opt.profile_grid_step <= 0.0 || opt.profile_grid_step > 0.05 ||
        (opt.assay_mode != "RNA" && opt.assay_mode != "ATAC" && opt.assay_mode != "GENERIC") ||
        (opt.likelihood != "beta_binomial" && opt.likelihood != "binomial")) {
        std::fprintf(stderr, "ERROR: invalid numeric/tag option\n");
        std::exit(2);
    }
    if (opt.rho_mode != "free" && !std::isfinite(opt.pooled_rho) && opt.rho_reference.empty()) {
        std::fprintf(stderr,
            "ERROR: --rho_mode %s requires --pooled_rho or --rho_reference\n", opt.rho_mode.c_str());
        std::exit(2);
    }
    if (opt.site_calibration_stratum.empty()) opt.site_calibration_stratum = opt.library_id;
    return opt;
}

std::vector<MtSite> load_mt_sites(const Options& opt,
                                  std::vector<std::string>& sample_names,
                                  std::unordered_map<std::string, int>& sample_index,
                                  std::unordered_map<int64_t, int>& position_to_site) {
    htsFile* fp = bcf_open(opt.vcf.c_str(), "r");
    if (!fp) throw std::runtime_error("Could not open mt panel: " + opt.vcf);
    bcf_hdr_t* header = bcf_hdr_read(fp);
    if (!header) {
        hts_close(fp);
        throw std::runtime_error("Could not read mt-panel header");
    }
    const int n_samples = bcf_hdr_nsamples(header);
    for (int i = 0; i < n_samples; ++i) {
        sample_names.emplace_back(header->samples[i]);
        sample_index[header->samples[i]] = i;
    }
    const int mt_rid = bcf_hdr_name2id(header, opt.mito_chrom.c_str());
    if (mt_rid < 0) throw std::runtime_error(
        "Mitochondrial contig missing from mt panel: " + opt.mito_chrom);

    const std::vector<Interval> mask = load_bed_for_contig(opt.mt_mask_bed, opt.mito_chrom);
    std::vector<MtSite> sites;
    bcf1_t* record = bcf_init();
    int32_t* genotypes = nullptr;
    int genotype_capacity = 0;
    while (bcf_read(fp, header, record) >= 0) {
        if (record->rid != mt_rid || masked_position(record->pos, mask)) continue;
        bcf_unpack(record, BCF_UN_STR);
        if (record->n_allele != 2 || std::strlen(record->d.allele[0]) != 1 ||
            std::strlen(record->d.allele[1]) != 1) continue;
        const char ref = static_cast<char>(std::toupper(record->d.allele[0][0]));
        const char alt = static_cast<char>(std::toupper(record->d.allele[1][0]));
        if (!valid_base(ref) || !valid_base(alt) || ref == alt) continue;
        if (position_to_site.count(record->pos)) {
            std::fprintf(stderr, "WARNING: duplicate mt-panel position %s:%lld; keeping first\n",
                opt.mito_chrom.c_str(), static_cast<long long>(record->pos + 1));
            continue;
        }
        MtSite site;
        site.rid = record->rid;
        site.pos = record->pos;
        site.ref = ref;
        site.alt = alt;
        site.genotype_state.assign(n_samples, -1);
        // Production fits take parental states from the depth-qualified site
        // manifest.  Do not silently re-impose a GT requirement here.  GT is
        // parsed only for the explicitly requested legacy diagnostic mode.
        if (opt.legacy_panel_gt) {
            const int ngt = bcf_get_genotypes(
                header, record, &genotypes, &genotype_capacity);
            if (ngt > 0 && ngt % n_samples == 0) {
                const int ploidy = ngt / n_samples;
                for (int i = 0; i < n_samples; ++i) {
                    site.genotype_state[i] = static_cast<int8_t>(
                        homoplasmic_state(genotypes + i * ploidy, ploidy));
                }
            }
        }
        position_to_site[site.pos] = static_cast<int>(sites.size());
        sites.push_back(std::move(site));
    }
    std::free(genotypes);
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    hts_close(fp);
    if (sites.empty()) throw std::runtime_error("No usable chrM SNPs found in mt panel");
    return sites;
}


std::pair<std::string, std::string> canonical_parents(const std::string& a,
                                                       const std::string& b) {
    return a <= b ? std::make_pair(a, b) : std::make_pair(b, a);
}

std::string pair_manifest_key(const std::string& a, const std::string& b) {
    const std::pair<std::string, std::string> canonical = canonical_parents(a, b);
    return canonical.first + std::string(1, '\x1f') + canonical.second;
}

ManifestSiteClass parse_manifest_site_class(const std::string& value) {
    if (value == "A_ALT_B_REF") return MANIFEST_A_ALT_B_REF;
    if (value == "A_REF_B_ALT") return MANIFEST_A_REF_B_ALT;
    if (value == "AMBIENT_ONLY") return MANIFEST_AMBIENT_ONLY;
    throw std::runtime_error("Unknown mitochondrial site class in manifest: " + value);
}

const char* manifest_site_class_name(ManifestSiteClass value) {
    switch (value) {
        case MANIFEST_A_ALT_B_REF: return "A_ALT_B_REF";
        case MANIFEST_A_REF_B_ALT: return "A_REF_B_ALT";
        case MANIFEST_AMBIENT_ONLY: return "AMBIENT_ONLY";
        case MANIFEST_LEGACY_RATIO: return "LEGACY_GT_RATIO";
    }
    return "UNKNOWN";
}

int required_column(const std::unordered_map<std::string, int>& columns,
                    const std::string& name,
                    const std::string& filename) {
    const auto found = columns.find(name);
    if (found == columns.end()) {
        throw std::runtime_error("Site manifest is missing required column '" + name +
                                 "': " + filename);
    }
    return found->second;
}

std::unordered_map<std::string, PairManifest> load_site_manifest(
        const Options& opt,
        const std::vector<MtSite>& sites,
        const std::unordered_map<int64_t, int>& position_to_site) {
    std::unordered_map<std::string, PairManifest> result;
    if (opt.legacy_panel_gt) return result;

    const std::vector<std::string> lines = read_text_lines(opt.site_manifest);
    if (lines.empty()) {
        throw std::runtime_error("Mitochondrial site manifest is empty: " + opt.site_manifest);
    }
    std::vector<std::string> header = split_tab(lines[0]);
    if (header.size() < 2) header = split_ws(lines[0]);
    std::unordered_map<std::string, int> columns;
    for (size_t i = 0; i < header.size(); ++i) columns[header[i]] = static_cast<int>(i);

    const int library_col = required_column(columns, "library_id", opt.site_manifest);
    const int parent1_col = required_column(columns, "parent1", opt.site_manifest);
    const int parent2_col = required_column(columns, "parent2", opt.site_manifest);
    const int chrom_col = required_column(columns, "chrom", opt.site_manifest);
    const int pos_col = required_column(columns, "pos", opt.site_manifest);
    const int ref_col = required_column(columns, "ref", opt.site_manifest);
    const int alt_col = required_column(columns, "alt", opt.site_manifest);
    const int class_col = required_column(columns, "site_class", opt.site_manifest);
    const int state1_col = required_column(columns, "parent1_state", opt.site_manifest);
    const int state2_col = required_column(columns, "parent2_state", opt.site_manifest);
    const int max_col = std::max({library_col, parent1_col, parent2_col, chrom_col, pos_col,
                                  ref_col, alt_col, class_col, state1_col, state2_col});

    size_t matching_rows = 0;
    for (size_t line_index = 1; line_index < lines.size(); ++line_index) {
        if (lines[line_index].empty() || lines[line_index][0] == '#') continue;
        std::vector<std::string> fields = split_tab(lines[line_index]);
        if (fields.size() == 1) fields = split_ws(lines[line_index]);
        if (max_col >= static_cast<int>(fields.size())) {
            throw std::runtime_error("Malformed site-manifest row " +
                                     std::to_string(line_index + 1) + " in " +
                                     opt.site_manifest);
        }
        if (!same_library_key(fields[library_col], opt.library_id)) continue;
        ++matching_rows;
        if (fields[chrom_col] != opt.mito_chrom) {
            throw std::runtime_error("Manifest library " + opt.library_id +
                                     " contains non-mitochondrial contig: " +
                                     fields[chrom_col]);
        }

        int64_t pos0 = -1;
        int state1 = -1;
        int state2 = -1;
        try {
            pos0 = std::stoll(fields[pos_col]) - 1;
            state1 = std::stoi(fields[state1_col]);
            state2 = std::stoi(fields[state2_col]);
        } catch (...) {
            throw std::runtime_error("Invalid numeric field in site-manifest row " +
                                     std::to_string(line_index + 1));
        }
        if ((state1 != 0 && state1 != 1) || (state2 != 0 && state2 != 1)) {
            throw std::runtime_error("Manifest parental states must be 0 or 1 at row " +
                                     std::to_string(line_index + 1));
        }
        const auto site_found = position_to_site.find(pos0);
        if (site_found == position_to_site.end()) {
            throw std::runtime_error("Manifest site is absent from mitochondrial panel: " +
                                     fields[chrom_col] + ":" + fields[pos_col]);
        }
        const int site_index = site_found->second;
        const MtSite& panel_site = sites[site_index];
        if (fields[ref_col].size() != 1 || fields[alt_col].size() != 1 ||
            std::toupper(static_cast<unsigned char>(fields[ref_col][0])) != panel_site.ref ||
            std::toupper(static_cast<unsigned char>(fields[alt_col][0])) != panel_site.alt) {
            throw std::runtime_error("Manifest REF/ALT mismatch at " + fields[chrom_col] + ":" +
                                     fields[pos_col]);
        }

        const ManifestSiteClass declared_class = parse_manifest_site_class(fields[class_col]);
        if (fields[parent1_col].empty() || fields[parent2_col].empty() ||
            fields[parent1_col] == fields[parent2_col]) {
            throw std::runtime_error("Manifest parents must be distinct non-empty donor IDs at row " +
                                     std::to_string(line_index + 1));
        }
        if (state1 != state2) {
            const ManifestSiteClass expected_declared =
                state1 == 1 ? MANIFEST_A_ALT_B_REF : MANIFEST_A_REF_B_ALT;
            if (declared_class != expected_declared) {
                throw std::runtime_error("Manifest site class disagrees with parental states at " +
                                         fields[chrom_col] + ":" + fields[pos_col]);
            }
        } else if (declared_class != MANIFEST_AMBIENT_ONLY) {
            throw std::runtime_error("Ratio manifest class has equal parental states at " +
                                     fields[chrom_col] + ":" + fields[pos_col]);
        }

        std::string canonical1 = fields[parent1_col];
        std::string canonical2 = fields[parent2_col];
        if (canonical2 < canonical1) {
            std::swap(canonical1, canonical2);
            std::swap(state1, state2);
        }
        const ManifestSiteClass canonical_class = state1 == state2
            ? MANIFEST_AMBIENT_ONLY
            : (state1 == 1 ? MANIFEST_A_ALT_B_REF : MANIFEST_A_REF_B_ALT);

        const std::string key = pair_manifest_key(canonical1, canonical2);
        PairManifest& pair = result[key];
        if (pair.canonical_parent1.empty()) {
            pair.canonical_parent1 = canonical1;
            pair.canonical_parent2 = canonical2;
        }
        ManifestSite entry;
        entry.site_index = site_index;
        entry.canonical_parent1_state = static_cast<int8_t>(state1);
        entry.canonical_parent2_state = static_cast<int8_t>(state2);
        entry.site_class = canonical_class;
        const auto prior = pair.sites.find(site_index);
        if (prior != pair.sites.end()) {
            if (prior->second.canonical_parent1_state != entry.canonical_parent1_state ||
                prior->second.canonical_parent2_state != entry.canonical_parent2_state ||
                prior->second.site_class != entry.site_class) {
                throw std::runtime_error("Conflicting duplicate site-manifest row for pair " +
                                         canonical1 + "+" + canonical2 + " at " +
                                         opt.mito_chrom + ":" + std::to_string(pos0 + 1));
            }
            continue;
        }
        pair.sites[site_index] = entry;
        if (entry.site_class == MANIFEST_AMBIENT_ONLY) ++pair.ambient_only_sites_available;
        else ++pair.ratio_sites_available;
    }

    if (matching_rows == 0 || result.empty()) {
        throw std::runtime_error("No site-manifest rows found for library_id '" +
                                 opt.library_id + "' in " + opt.site_manifest);
    }
    return result;
}


PairManifest derive_pair_manifest_from_panel(
        const std::string& parent_a,
        const std::string& parent_b,
        const mt_evidence::Panel& panel,
        const std::vector<MtSite>& sites,
        const std::unordered_map<int64_t, int>& position_to_site) {
    const std::pair<std::string, std::string> canonical =
        canonical_parents(parent_a, parent_b);
    PairManifest pair;
    pair.canonical_parent1 = canonical.first;
    pair.canonical_parent2 = canonical.second;

    const auto p1_found = panel.sample_index.find(pair.canonical_parent1);
    const auto p2_found = panel.sample_index.find(pair.canonical_parent2);
    if (p1_found == panel.sample_index.end() || p2_found == panel.sample_index.end()) {
        throw std::runtime_error(
            "Reconciled pair parent is missing from the mitochondrial panel: " +
            pair.canonical_parent1 + "+" + pair.canonical_parent2);
    }
    const int p1_index = p1_found->second;
    const int p2_index = p2_found->second;

    for (const mt_evidence::Site& evidence_site : panel.sites) {
        const auto local_found = position_to_site.find(evidence_site.pos);
        if (local_found == position_to_site.end()) continue;
        const int site_index = local_found->second;
        const MtSite& local_site = sites[site_index];
        const char evidence_ref = static_cast<char>(
            std::toupper(static_cast<unsigned char>(evidence_site.ref)));
        const char evidence_alt = static_cast<char>(
            std::toupper(static_cast<unsigned char>(evidence_site.alt)));
        if (evidence_ref != local_site.ref || evidence_alt != local_site.alt) {
            throw std::runtime_error(
                "Mitochondrial panel REF/ALT mismatch while deriving reconciled pair " +
                pair.canonical_parent1 + "+" + pair.canonical_parent2 + " at " +
                panel.chrom + ":" + std::to_string(evidence_site.pos + 1));
        }

        const int state1 = evidence_site.state[p1_index];
        const int state2 = evidence_site.state[p2_index];
        if (state1 < 0 || state2 < 0) continue;

        ManifestSite entry;
        entry.site_index = site_index;
        entry.canonical_parent1_state = static_cast<int8_t>(state1);
        entry.canonical_parent2_state = static_cast<int8_t>(state2);

        if (state1 != state2) {
            entry.site_class = state1 == 1
                ? MANIFEST_A_ALT_B_REF : MANIFEST_A_REF_B_ALT;
            ++pair.ratio_sites_available;
        } else {
            bool other_donor_differs = false;
            for (int8_t state : evidence_site.state) {
                if (state >= 0 && state != state1) {
                    other_donor_differs = true;
                    break;
                }
            }
            if (!other_donor_differs) continue;
            entry.site_class = MANIFEST_AMBIENT_ONLY;
            ++pair.ambient_only_sites_available;
        }
        pair.sites[site_index] = entry;
    }
    return pair;
}

std::unordered_map<std::string, PairManifest> derive_missing_pair_manifests(
        const Options& opt,
        const std::vector<CellInfo>& cells,
        const std::unordered_map<std::string, PairManifest>& library_pair_manifests,
        const std::vector<MtSite>& sites,
        const std::unordered_map<int64_t, int>& position_to_site) {
    std::map<std::string, std::pair<std::string, std::string>> missing_pairs;
    for (const CellInfo& cell : cells) {
        if (cell.parent1 == cell.parent2) continue;
        const std::string key = pair_manifest_key(cell.parent1, cell.parent2);
        if (library_pair_manifests.count(key)) continue;
        missing_pairs.emplace(key, canonical_parents(cell.parent1, cell.parent2));
    }

    std::unordered_map<std::string, PairManifest> derived;
    if (missing_pairs.empty() || opt.legacy_panel_gt) return derived;

    // Match the donor-state interpretation used by mt_identity_score. Existing
    // library manifests remain authoritative for expected pairs; this fallback
    // is only for reconciliation-supported pairs absent from the historical
    // library design.
    const mt_evidence::Panel evidence_panel = mt_evidence::load_panel(
        opt.vcf, opt.mito_chrom, kReconciledPairPanelMinDepth,
        kReconciledPairHomoplasmyAf, nullptr);

    for (const auto& item : missing_pairs) {
        PairManifest pair = derive_pair_manifest_from_panel(
            item.second.first, item.second.second, evidence_panel, sites,
            position_to_site);
        std::fprintf(
            stderr,
            "Derived reconciled MT pair %s+%s from existing panel: %d ratio sites, %d ambient-only sites\n",
            pair.canonical_parent1.c_str(), pair.canonical_parent2.c_str(),
            pair.ratio_sites_available, pair.ambient_only_sites_available);
        derived.emplace(item.first, std::move(pair));
    }
    return derived;
}

std::vector<CellInfo> load_assignments(const std::string& filename,
                                       const std::string& calibration_reference,
                                       const std::string& library_id,
                                       const std::unordered_map<std::string, int>& sample_index,
                                       std::unordered_map<std::string, int>& barcode_to_cell) {
    std::unordered_map<std::string, std::string> assignment_identity;
    std::vector<std::pair<std::string, std::string>> assignment_rows;
    for (const std::string& line : read_text_lines(filename)) {
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> fields = split_ws(line);
        if (fields.size() < 2 || fields[0] == "barcode") continue;
        if (!assignment_identity.emplace(fields[0], fields[1]).second) {
            throw std::runtime_error("Duplicate barcode in assignments: " + fields[0]);
        }
        assignment_rows.emplace_back(fields[0], fields[1]);
    }

    std::vector<CellInfo> cells;
    if (!calibration_reference.empty()) {
        const std::vector<std::string> lines = read_text_lines(calibration_reference);
        if (lines.empty()) throw std::runtime_error(
            "Calibration-reference file is empty: " + calibration_reference);
        const std::vector<std::string> header = split_tab(lines[0]);
        auto find_col = [&](const std::string& name, bool required) {
            for (size_t i = 0; i < header.size(); ++i) if (header[i] == name) return static_cast<int>(i);
            if (required) throw std::runtime_error(
                "Calibration-reference file missing required column: " + name);
            return -1;
        };
        const int barcode_col = find_col("barcode", true);
        const int p1_col = find_col("canonical_parent1", true);
        const int p2_col = find_col("canonical_parent2", true);
        const int true_col = find_col("true_parent", true);
        const int lib_col = find_col("library_id", false);
        for (size_t line_index = 1; line_index < lines.size(); ++line_index) {
            if (lines[line_index].empty() || lines[line_index][0] == '#') continue;
            const std::vector<std::string> fields = split_tab(lines[line_index]);
            const int max_col = std::max({barcode_col, p1_col, p2_col, true_col, lib_col});
            if (static_cast<int>(fields.size()) <= max_col) {
                throw std::runtime_error("Short calibration-reference row " +
                                         std::to_string(line_index + 1));
            }
            if (lib_col >= 0 && !same_library_key(fields[lib_col], library_id) &&
                fields[lib_col] != "*" && fields[lib_col] != "ALL") continue;
            const std::string& barcode = fields[barcode_col];
            std::string parent1 = fields[p1_col];
            std::string parent2 = fields[p2_col];
            const std::string true_parent = fields[true_col];
            if (parent1.empty() || parent2.empty() || parent1 == parent2) {
                throw std::runtime_error("Calibration-reference parents must be distinct at row " +
                                         std::to_string(line_index + 1));
            }
            if (parent2 < parent1) std::swap(parent1, parent2);
            if (true_parent != parent1 && true_parent != parent2) {
                throw std::runtime_error("Calibration-reference true_parent must equal one reference parent at row " +
                                         std::to_string(line_index + 1));
            }
            const auto p1 = sample_index.find(parent1);
            const auto p2 = sample_index.find(parent2);
            if (p1 == sample_index.end() || p2 == sample_index.end()) {
                throw std::runtime_error("Calibration-reference parent missing from mt-panel VCF: " +
                                         parent1 + "+" + parent2);
            }
            const auto original = assignment_identity.find(barcode);
            if (original == assignment_identity.end()) {
                throw std::runtime_error("Calibration-reference barcode missing from assignments: " + barcode);
            }
            if (barcode_to_cell.count(barcode)) {
                throw std::runtime_error("Duplicate barcode in calibration-reference file: " + barcode);
            }
            CellInfo cell;
            cell.barcode = barcode;
            cell.identity = parent1 + "+" + parent2;
            cell.parent1 = parent1;
            cell.parent2 = parent2;
            cell.parent1_index = p1->second;
            cell.parent2_index = p2->second;
            cell.calibration_true_parent = true_parent;
            cell.original_identity = original->second;
            barcode_to_cell[barcode] = static_cast<int>(cells.size());
            cells.push_back(std::move(cell));
        }
        if (cells.empty()) throw std::runtime_error(
            "No calibration-reference rows matched library_id '" + library_id + "'");
        return cells;
    }

    for (const auto& item : assignment_rows) {
        const std::string& barcode = item.first;
        const std::string& identity = item.second;
        const size_t plus = identity.find('+');
        if (plus == std::string::npos || identity.find('+', plus + 1) != std::string::npos) continue;
        const std::string parent1 = identity.substr(0, plus);
        const std::string parent2 = identity.substr(plus + 1);
        const auto p1 = sample_index.find(parent1);
        const auto p2 = sample_index.find(parent2);
        if (p1 == sample_index.end() || p2 == sample_index.end()) {
            throw std::runtime_error("Assignment parent missing from mt-panel VCF: " + identity);
        }
        CellInfo cell;
        cell.barcode = barcode;
        cell.identity = identity;
        cell.parent1 = parent1;
        cell.parent2 = parent2;
        cell.parent1_index = p1->second;
        cell.parent2_index = p2->second;
        cell.original_identity = identity;
        barcode_to_cell[barcode] = static_cast<int>(cells.size());
        cells.push_back(std::move(cell));
    }
    if (cells.empty()) throw std::runtime_error(
        "No two-parent fusion assignments were found in: " + filename);
    return cells;
}


std::string file_fingerprint(const std::string& filename) {
    if (filename.empty()) return "NA";
    std::ifstream input(filename, std::ios::binary);
    if (!input) return "UNAVAILABLE";
    uint64_t hash = 1469598103934665603ULL;
    char buffer[1 << 16];
    while (input) {
        input.read(buffer, sizeof(buffer));
        const std::streamsize n = input.gcount();
        for (std::streamsize i = 0; i < n; ++i) {
            hash ^= static_cast<unsigned char>(buffer[i]);
            hash *= 1099511628211ULL;
        }
    }
    std::ostringstream out;
    out << "fnv1a64:" << std::hex << std::setfill('0') << std::setw(16) << hash;
    return out.str();
}

std::string path_basename(const std::string& path) {
    if (path.empty()) return "NA";
    const size_t pos = path.find_last_of("/\\");
    return pos == std::string::npos ? path : path.substr(pos + 1);
}

bool wildcard_match(const std::string& value, const std::string& requested) {
    return same_library_key(value, requested) || value == "*" || value == "ALL" || value.empty();
}

SiteCalibrationData load_site_calibration(
        const Options& opt,
        const std::unordered_map<int64_t, int>& position_to_site) {
    SiteCalibrationData data;
    if (opt.site_calibration.empty()) return data;
    const std::vector<std::string> lines = read_text_lines(opt.site_calibration);
    if (lines.empty()) throw std::runtime_error("Site-calibration file is empty: " + opt.site_calibration);
    const std::vector<std::string> header = split_tab(lines[0]);
    auto col = [&](const std::string& name, bool required) {
        for (size_t i = 0; i < header.size(); ++i) if (header[i] == name) return static_cast<int>(i);
        if (required) throw std::runtime_error("Site-calibration file missing required column: " + name);
        return -1;
    };
    const int id_col = col("calibration_id", true);
    const int assay_col = col("assay_mode", true);
    const int stratum_col = col("library_or_stratum", true);
    const int p1_col = col("canonical_parent1", true);
    const int p2_col = col("canonical_parent2", true);
    const int chrom_col = col("chrom", true);
    const int pos_col = col("pos", true);
    const int parent_col = col("parent", true);
    const int prob_col = col("alt_probability", true);
    for (size_t line_index = 1; line_index < lines.size(); ++line_index) {
        if (lines[line_index].empty() || lines[line_index][0] == '#') continue;
        const std::vector<std::string> fields = split_tab(lines[line_index]);
        const int max_col = std::max({id_col, assay_col, stratum_col, p1_col, p2_col,
                                      chrom_col, pos_col, parent_col, prob_col});
        if (static_cast<int>(fields.size()) <= max_col) {
            throw std::runtime_error("Short site-calibration row " + std::to_string(line_index + 1));
        }
        if (!wildcard_match(fields[assay_col], opt.assay_mode) ||
            !wildcard_match(fields[stratum_col], opt.site_calibration_stratum)) continue;
        if (fields[chrom_col] != opt.mito_chrom) continue;
        int64_t pos0 = -1;
        double probability = std::numeric_limits<double>::quiet_NaN();
        try {
            pos0 = std::stoll(fields[pos_col]) - 1;
            probability = std::stod(fields[prob_col]);
        } catch (...) {
            throw std::runtime_error("Invalid site-calibration numeric field at row " +
                                     std::to_string(line_index + 1));
        }
        if (probability <= 0.0 || probability >= 1.0) {
            throw std::runtime_error("Site-calibration alt_probability must be strictly between 0 and 1 at row " +
                                     std::to_string(line_index + 1));
        }
        const auto site_it = position_to_site.find(pos0);
        if (site_it == position_to_site.end()) continue;
        std::string p1 = fields[p1_col];
        std::string p2 = fields[p2_col];
        if (p1.empty() || p2.empty() || p1 == p2) {
            throw std::runtime_error("Invalid site-calibration parental pair at row " +
                                     std::to_string(line_index + 1));
        }
        if (p2 < p1) std::swap(p1, p2);
        const std::string parent = fields[parent_col];
        if (parent != p1 && parent != p2) {
            throw std::runtime_error("Site-calibration parent is not a member of its pair at row " +
                                     std::to_string(line_index + 1));
        }
        int priority = 0;
        if (fields[assay_col] == opt.assay_mode) priority += 2;
        if (same_library_key(fields[stratum_col], opt.site_calibration_stratum)) priority += 1;
        SiteCalibrationRecord& record = data.by_pair[pair_manifest_key(p1, p2)][site_it->second];
        bool* has = parent == p1 ? &record.has_parent1 : &record.has_parent2;
        double* dest = parent == p1 ? &record.parent1_alt_probability : &record.parent2_alt_probability;
        int* dest_priority = parent == p1 ? &record.parent1_priority : &record.parent2_priority;
        if (*has && priority == *dest_priority && std::fabs(*dest - probability) > 1e-12) {
            throw std::runtime_error("Conflicting site-calibration rows at equal priority for " +
                                     p1 + "+" + p2 + " " + opt.mito_chrom + ":" +
                                     std::to_string(pos0 + 1));
        }
        if (!*has || priority > *dest_priority) {
            *has = true;
            *dest = probability;
            *dest_priority = priority;
            if (record.calibration_id.empty()) record.calibration_id = fields[id_col];
            else if (record.calibration_id != fields[id_col]) record.calibration_id = "MIXED";
        }
        data.calibration_ids.insert(fields[id_col]);
        ++data.loaded_rows;
    }
    return data;
}

std::vector<RhoReferenceEntry> load_rho_reference(const Options& opt) {
    std::vector<RhoReferenceEntry> rows;
    if (opt.rho_reference.empty()) return rows;
    const std::vector<std::string> lines = read_text_lines(opt.rho_reference);
    if (lines.empty()) throw std::runtime_error("Rho-reference file is empty: " + opt.rho_reference);
    const std::vector<std::string> header = split_tab(lines[0]);
    auto col = [&](const std::string& name, bool required) {
        for (size_t i = 0; i < header.size(); ++i) if (header[i] == name) return static_cast<int>(i);
        if (required) throw std::runtime_error("Rho-reference file missing required column: " + name);
        return -1;
    };
    const int assay_col = col("assay_mode", true);
    const int lib_col = col("library_id", true);
    const int p1_col = col("canonical_parent1", true);
    const int p2_col = col("canonical_parent2", true);
    const int rho_col = col("rho", true);
    for (size_t i = 1; i < lines.size(); ++i) {
        if (lines[i].empty() || lines[i][0] == '#') continue;
        const auto fields = split_tab(lines[i]);
        const int max_col = std::max({assay_col, lib_col, p1_col, p2_col, rho_col});
        if (static_cast<int>(fields.size()) <= max_col) throw std::runtime_error(
            "Short rho-reference row " + std::to_string(i + 1));
        RhoReferenceEntry row;
        row.assay_mode = fields[assay_col];
        row.library_id = fields[lib_col];
        row.canonical_parent1 = fields[p1_col];
        row.canonical_parent2 = fields[p2_col];
        if (row.canonical_parent2 < row.canonical_parent1 && row.canonical_parent1 != "*" &&
            row.canonical_parent2 != "*") std::swap(row.canonical_parent1, row.canonical_parent2);
        try { row.rho = std::stod(fields[rho_col]); }
        catch (...) { throw std::runtime_error("Invalid rho-reference value at row " + std::to_string(i + 1)); }
        if (row.rho < 0.0 || row.rho > opt.overdispersion_max) throw std::runtime_error(
            "Rho-reference value outside configured bounds at row " + std::to_string(i + 1));
        rows.push_back(std::move(row));
    }
    return rows;
}

double resolve_pooled_rho(const Options& opt,
                          const std::vector<RhoReferenceEntry>& rows,
                          const std::string& canonical_parent1,
                          const std::string& canonical_parent2) {
    double best = opt.pooled_rho;
    int best_priority = std::isfinite(best) ? -1 : -1000;
    for (const RhoReferenceEntry& row : rows) {
        if (!wildcard_match(row.assay_mode, opt.assay_mode) ||
            !wildcard_match(row.library_id, opt.library_id)) continue;
        const bool pair_exact = row.canonical_parent1 == canonical_parent1 &&
                                row.canonical_parent2 == canonical_parent2;
        const bool pair_wild = (row.canonical_parent1 == "*" || row.canonical_parent1 == "ALL") &&
                               (row.canonical_parent2 == "*" || row.canonical_parent2 == "ALL");
        if (!pair_exact && !pair_wild) continue;
        int priority = 0;
        if (row.assay_mode == opt.assay_mode) priority += 4;
        if (same_library_key(row.library_id, opt.library_id)) priority += 2;
        if (pair_exact) priority += 1;
        if (priority > best_priority) {
            best_priority = priority;
            best = row.rho;
        }
    }
    return best;
}

std::unordered_set<std::string> load_empty_barcodes(
        const std::string& filename,
        const std::unordered_map<std::string, int>& target_barcodes) {
    std::unordered_set<std::string> result;
    if (filename.empty()) return result;
    uint64_t overlap = 0;
    for (const std::string& line : read_text_lines(filename)) {
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> fields = split_ws(line);
        if (fields.empty() || fields[0] == "barcode") continue;
        if (target_barcodes.count(fields[0])) {
            ++overlap;
            continue;
        }
        result.insert(fields[0]);
    }
    if (result.empty()) throw std::runtime_error(
        "Empty-barcode file yielded no non-cell barcodes: " + filename);
    if (overlap) std::fprintf(stderr,
        "WARNING: removed %llu target-cell barcodes from empty-droplet list\n",
        static_cast<unsigned long long>(overlap));
    return result;
}

std::unordered_map<std::string, double> load_fixed_ambient_fractions(
        const std::string& filename) {
    std::unordered_map<std::string, double> result;
    if (filename.empty()) return result;
    const std::vector<std::string> lines = read_text_lines(filename);
    if (lines.empty()) return result;

    int barcode_col = 0;
    int value_col = 1;
    size_t start = 0;
    const std::vector<std::string> header = split_ws(lines[0]);
    bool has_header = false;
    bool found_barcode = false;
    bool found_value = false;
    for (size_t i = 0; i < header.size(); ++i) {
        if (header[i] == "barcode") {
            barcode_col = static_cast<int>(i);
            has_header = true;
            found_barcode = true;
        }
        if (header[i] == "ambient_fraction" || header[i] == "contam_rate_c") {
            value_col = static_cast<int>(i);
            has_header = true;
            found_value = true;
        }
    }
    if (has_header) {
        if (!found_barcode || !found_value) {
            throw std::runtime_error(
                "Ambient-fraction header must contain barcode and either "
                "ambient_fraction or contam_rate_c: " + filename);
        }
        start = 1;
    }

    for (size_t i = start; i < lines.size(); ++i) {
        if (lines[i].empty() || lines[i][0] == '#') continue;
        const std::vector<std::string> fields = split_ws(lines[i]);
        if (barcode_col >= static_cast<int>(fields.size()) ||
            value_col >= static_cast<int>(fields.size())) continue;
        try {
            const double value = std::stod(fields[value_col]);
            if (value >= 0.0 && value <= 1.0 && std::isfinite(value)) {
                result[fields[barcode_col]] = value;
            }
        } catch (...) {
            continue;
        }
    }
    return result;
}

void load_ambient_profile(const std::string& filename,
                          const std::string& mito_chrom,
                          const std::vector<MtSite>& sites,
                          const std::unordered_map<int64_t, int>& position_to_site,
                          int min_molecules,
                          std::vector<AmbientSite>& ambient) {
    const std::vector<std::string> lines = read_text_lines(filename);
    if (lines.empty()) throw std::runtime_error("Ambient profile is empty: " + filename);
    std::vector<std::string> header = split_tab(lines[0]);
    if (header.size() < 2) header = split_ws(lines[0]);
    std::unordered_map<std::string, int> column;
    for (size_t i = 0; i < header.size(); ++i) column[header[i]] = static_cast<int>(i);
    if (!column.count("pos") || !column.count("alt_fraction")) {
        throw std::runtime_error(
            "Ambient profile must contain pos and alt_fraction columns");
    }
    const int chrom_col = column.count("chrom") ? column["chrom"] : -1;
    const int ref_col = column.count("ref") ? column["ref"] : -1;
    const int alt_col = column.count("alt") ? column["alt"] : -1;
    const int ref_count_col = column.count("ref_molecules") ? column["ref_molecules"] : -1;
    const int alt_count_col = column.count("alt_molecules") ? column["alt_molecules"] : -1;

    for (size_t line_index = 1; line_index < lines.size(); ++line_index) {
        if (lines[line_index].empty() || lines[line_index][0] == '#') continue;
        std::vector<std::string> fields = split_tab(lines[line_index]);
        if (fields.size() == 1) fields = split_ws(lines[line_index]);
        const int max_required = std::max(column["pos"], column["alt_fraction"]);
        if (max_required >= static_cast<int>(fields.size())) continue;
        if (chrom_col >= 0 && chrom_col < static_cast<int>(fields.size()) &&
            fields[chrom_col] != mito_chrom) continue;
        try {
            const int64_t pos0 = std::stoll(fields[column["pos"]]) - 1;
            const auto site_it = position_to_site.find(pos0);
            if (site_it == position_to_site.end()) continue;
            const int site_index = site_it->second;
            if (ref_col >= 0 && ref_col < static_cast<int>(fields.size()) &&
                !fields[ref_col].empty() &&
                std::toupper(fields[ref_col][0]) != sites[site_index].ref) continue;
            if (alt_col >= 0 && alt_col < static_cast<int>(fields.size()) &&
                !fields[alt_col].empty() &&
                std::toupper(fields[alt_col][0]) != sites[site_index].alt) continue;
            const double fraction = std::stod(fields[column["alt_fraction"]]);
            if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0) continue;
            AmbientSite& a = ambient[site_index];
            a.alt_fraction = fraction;
            if (ref_count_col >= 0 && alt_count_col >= 0 &&
                ref_count_col < static_cast<int>(fields.size()) &&
                alt_count_col < static_cast<int>(fields.size())) {
                a.ref = static_cast<uint64_t>(std::stoull(fields[ref_count_col]));
                a.alt = static_cast<uint64_t>(std::stoull(fields[alt_count_col]));
                a.usable = static_cast<int>(a.total()) >= min_molecules;
            } else {
                a.usable = true;
            }
        } catch (...) {
            continue;
        }
    }
}

int read_for_pileup(void* raw, bam1_t* record) {
    PileupReader* reader = static_cast<PileupReader*>(raw);
    while (sam_itr_next(reader->bam, reader->iterator, record) >= 0) {
        ReadStats& stats = *reader->stats;
        ++stats.seen;
        const uint16_t flag = record->core.flag;
        if (flag & BAM_FUNMAP) { ++stats.reject_unmapped; continue; }
        if (flag & BAM_FSECONDARY) { ++stats.reject_secondary; continue; }
        if (flag & BAM_FSUPPLEMENTARY) { ++stats.reject_supplementary; continue; }
        if (flag & BAM_FQCFAIL) { ++stats.reject_qcfail; continue; }
        if (!reader->keep_duplicates && (flag & BAM_FDUP)) {
            ++stats.reject_duplicate;
            continue;
        }
        if (record->core.qual < reader->min_mapq) {
            ++stats.reject_mapq;
            continue;
        }
        if (!reader->allow_ambiguous_alignments) {
            uint8_t* nh = bam_aux_get(record, "NH");
            if ((nh && bam_aux2i(nh) > 1) || bam_aux_get(record, "SA") ||
                bam_aux_get(record, "XA")) {
                ++stats.reject_multimapping;
                continue;
            }
        }
        ++stats.accepted_for_pileup;
        return 0;
    }
    return -1;
}

const char* aux_string(const bam1_t* record, const std::string& tag) {
    uint8_t* value = bam_aux_get(record, tag.c_str());
    if (!value) return nullptr;
    return bam_aux2Z(value);
}

void count_molecules(const Options& opt,
                     const std::vector<MtSite>& sites,
                     const std::unordered_map<int64_t, int>& position_to_site,
                     const std::unordered_map<std::string, int>& barcode_to_cell,
                     const std::unordered_set<std::string>& empty_barcodes,
                     std::vector<CellInfo>& cells,
                     std::vector<AmbientSite>& ambient,
                     ReadStats& stats) {
    samFile* bam = sam_open(opt.bam.c_str(), "r");
    if (!bam) throw std::runtime_error("Could not open BAM/CRAM: " + opt.bam);
    hts_set_threads(bam, opt.threads);
    sam_hdr_t* header = sam_hdr_read(bam);
    if (!header) {
        sam_close(bam);
        throw std::runtime_error("Could not read BAM header");
    }
    const int mt_tid = sam_hdr_name2tid(header, opt.mito_chrom.c_str());
    if (mt_tid < 0) throw std::runtime_error(
        "Mitochondrial contig missing from BAM header: " + opt.mito_chrom);
    hts_idx_t* index = sam_index_load(bam, opt.bam.c_str());
    if (!index) throw std::runtime_error("Could not load BAM/CRAM index: " + opt.bam);
    hts_itr_t* iterator = sam_itr_querys(index, header, opt.mito_chrom.c_str());
    if (!iterator) throw std::runtime_error("Could not create BAM chrM iterator");

    PileupReader reader;
    reader.bam = bam;
    reader.iterator = iterator;
    reader.min_mapq = opt.min_mapq;
    reader.keep_duplicates = opt.keep_duplicates;
    reader.allow_ambiguous_alignments = opt.allow_ambiguous_alignments;
    reader.stats = &stats;

    bam_plp_t pileup = bam_plp_init(read_for_pileup, &reader);
    if (!pileup) throw std::runtime_error("Could not initialize pileup");
    bam_plp_set_maxcnt(pileup, opt.pileup_max_depth);

    int tid = -1;
    int pos = -1;
    int depth = 0;
    const bam_pileup1_t* entries = nullptr;
    while ((entries = bam_plp_auto(pileup, &tid, &pos, &depth)) != nullptr) {
        if (tid != mt_tid) continue;
        const auto site_it = position_to_site.find(pos);
        if (site_it == position_to_site.end()) continue;
        const int site_index = site_it->second;
        const MtSite& site = sites[site_index];
        std::unordered_map<std::string, LocalMolecule> molecules;
        molecules.reserve(static_cast<size_t>(std::max(depth, 1)));

        for (int i = 0; i < depth; ++i) {
            const bam_pileup1_t& p = entries[i];
            if (p.is_del || p.is_refskip || p.qpos < 0) continue;
            ++stats.pileup_bases;
            const bam1_t* record = p.b;
            const uint8_t* qualities = bam_get_qual(record);
            if (qualities && qualities[p.qpos] < opt.min_baseq) {
                ++stats.reject_baseq;
                continue;
            }
            const char* barcode_ptr = aux_string(record, opt.barcode_tag);
            if (!barcode_ptr || !*barcode_ptr) {
                ++stats.reject_missing_barcode;
                continue;
            }
            const std::string barcode(barcode_ptr);
            const auto cell_it = barcode_to_cell.find(barcode);
            const bool is_empty = empty_barcodes.count(barcode) > 0;
            if (cell_it == barcode_to_cell.end() && !is_empty) {
                ++stats.reject_irrelevant_barcode;
                continue;
            }

            std::string molecule;
            const char* umi_ptr = aux_string(record, opt.umi_tag);
            if (umi_ptr && *umi_ptr) {
                molecule = umi_ptr;
            } else if (opt.allow_no_umi) {
                molecule = bam_get_qname(record);
            } else {
                ++stats.reject_missing_umi;
                continue;
            }

            const uint8_t* sequence = bam_get_seq(record);
            const char base = static_cast<char>(
                std::toupper(seq_nt16_str[bam_seqi(sequence, p.qpos)]));
            uint8_t allele = 2;
            if (base == site.ref) allele = 0;
            else if (base == site.alt) allele = 1;
            else {
                ++stats.reject_nonpanel_allele;
                continue;
            }

            std::string key;
            key.reserve(barcode.size() + molecule.size() + 1);
            key.append(barcode);
            key.push_back('\x1f');
            key.append(molecule);
            auto found = molecules.find(key);
            if (found == molecules.end()) {
                LocalMolecule observed;
                observed.target_index = is_empty ? -1 : cell_it->second;
                observed.allele = allele;
                molecules.emplace(std::move(key), observed);
            } else if (found->second.allele != allele) {
                found->second.allele = 2;
            }
        }

        for (const auto& entry : molecules) {
            const LocalMolecule& molecule = entry.second;
            if (molecule.allele > 1) {
                ++stats.conflicting_molecules;
                continue;
            }
            if (molecule.target_index >= 0) {
                AlleleCount& count = cells[molecule.target_index].counts[site_index];
                if (molecule.allele == 0) ++count.ref;
                else ++count.alt;
            } else {
                if (molecule.allele == 0) ++ambient[site_index].ref;
                else ++ambient[site_index].alt;
            }
            ++stats.accepted_observations;
        }
    }

    bam_plp_destroy(pileup);
    hts_itr_destroy(iterator);
    hts_idx_destroy(index);
    sam_hdr_destroy(header);
    sam_close(bam);
}


struct AtacFragmentAccumulator {
    int target_index = -2;  // >=0 cell, -1 empty, -2 unset
    bool seen_first = false;
    bool seen_second = false;
    int accepted_reads = 0;
    std::unordered_map<int, uint8_t> alleles;
};

bool accept_atac_read(const Options& opt, const bam1_t* record, int mt_tid, ReadStats& stats) {
    ++stats.seen;
    ++stats.atac_reads_considered;
    const uint16_t flag = record->core.flag;
    if (flag & BAM_FUNMAP) { ++stats.reject_unmapped; return false; }
    if (flag & BAM_FSECONDARY) { ++stats.reject_secondary; return false; }
    if (flag & BAM_FSUPPLEMENTARY) { ++stats.reject_supplementary; return false; }
    if (flag & BAM_FQCFAIL) { ++stats.reject_qcfail; return false; }
    if (!opt.keep_duplicates && (flag & BAM_FDUP)) { ++stats.reject_duplicate; return false; }
    if (record->core.qual < opt.min_mapq) { ++stats.reject_mapq; return false; }
    if (!opt.allow_ambiguous_alignments) {
        uint8_t* nh = bam_aux_get(record, "NH");
        if ((nh && bam_aux2i(nh) > 1) || bam_aux_get(record, "SA") || bam_aux_get(record, "XA")) {
            ++stats.reject_multimapping;
            return false;
        }
    }
    if (!(flag & BAM_FPAIRED)) {
        if (!opt.atac_include_singletons) { ++stats.atac_reject_unpaired; return false; }
    } else if ((flag & BAM_FMUNMAP) || record->core.mtid != mt_tid) {
        if (!opt.atac_include_singletons) { ++stats.atac_reject_mate_off_mito; return false; }
    }
    ++stats.accepted_for_pileup;
    return true;
}

std::unordered_map<int, uint8_t> extract_atac_site_alleles(
        const Options& opt,
        const bam1_t* record,
        const std::vector<MtSite>& sites,
        const std::unordered_map<int64_t, int>& position_to_site,
        ReadStats& stats) {
    std::unordered_map<int, uint8_t> result;
    int64_t ref_pos = record->core.pos;
    int query_pos = 0;
    const uint32_t* cigar = bam_get_cigar(record);
    const uint8_t* sequence = bam_get_seq(record);
    const uint8_t* qualities = bam_get_qual(record);
    for (uint32_t ci = 0; ci < record->core.n_cigar; ++ci) {
        const int op = bam_cigar_op(cigar[ci]);
        const int len = bam_cigar_oplen(cigar[ci]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (int k = 0; k < len; ++k) {
                const auto found = position_to_site.find(ref_pos + k);
                if (found == position_to_site.end()) continue;
                ++stats.pileup_bases;
                const int qpos = query_pos + k;
                if (qualities && qualities[qpos] < opt.min_baseq) {
                    ++stats.reject_baseq;
                    continue;
                }
                const MtSite& site = sites[found->second];
                const char base = static_cast<char>(
                    std::toupper(seq_nt16_str[bam_seqi(sequence, qpos)]));
                uint8_t allele = 2;
                if (base == site.ref) allele = 0;
                else if (base == site.alt) allele = 1;
                else {
                    ++stats.reject_nonpanel_allele;
                    continue;
                }
                result[found->second] = allele;
            }
            ref_pos += len;
            query_pos += len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            query_pos += len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            ref_pos += len;
        }
    }
    return result;
}

void count_atac_fragments(const Options& opt,
                          const std::vector<MtSite>& sites,
                          const std::unordered_map<int64_t, int>& position_to_site,
                          const std::unordered_map<std::string, int>& barcode_to_cell,
                          const std::unordered_set<std::string>& empty_barcodes,
                          std::vector<CellInfo>& cells,
                          std::vector<AmbientSite>& ambient,
                          ReadStats& stats) {
    samFile* bam = sam_open(opt.bam.c_str(), "r");
    if (!bam) throw std::runtime_error("Could not open BAM/CRAM: " + opt.bam);
    hts_set_threads(bam, opt.threads);
    sam_hdr_t* header = sam_hdr_read(bam);
    if (!header) { sam_close(bam); throw std::runtime_error("Could not read BAM header"); }
    const int mt_tid = sam_hdr_name2tid(header, opt.mito_chrom.c_str());
    if (mt_tid < 0) throw std::runtime_error(
        "Mitochondrial contig missing from BAM header: " + opt.mito_chrom);
    hts_idx_t* index = sam_index_load(bam, opt.bam.c_str());
    if (!index) throw std::runtime_error("Could not load BAM/CRAM index: " + opt.bam);
    hts_itr_t* iterator = sam_itr_querys(index, header, opt.mito_chrom.c_str());
    if (!iterator) throw std::runtime_error("Could not create BAM chrM iterator");

    std::unordered_map<std::string, AtacFragmentAccumulator> pending;
    pending.reserve(1 << 16);
    auto finalize = [&](AtacFragmentAccumulator& acc) {
        RawFragmentObservation fragment;
        fragment.sites.reserve(acc.alleles.size());
        for (const auto& item : acc.alleles) {
            if (item.second > 1) continue;
            fragment.sites.push_back({item.first, item.second});
        }
        if (fragment.sites.empty()) return;
        std::sort(fragment.sites.begin(), fragment.sites.end(),
                  [](const RawFragmentSite& a, const RawFragmentSite& b) {
                      return a.site_index < b.site_index;
                  });
        if (acc.target_index >= 0) {
            for (const RawFragmentSite& site : fragment.sites) {
                AlleleCount& count = cells[acc.target_index].counts[site.site_index];
                if (site.allele == 0) ++count.ref; else ++count.alt;
            }
            cells[acc.target_index].fragments.push_back(fragment);
        } else if (acc.target_index == -1) {
            for (const RawFragmentSite& site : fragment.sites) {
                if (site.allele == 0) ++ambient[site.site_index].ref;
                else ++ambient[site.site_index].alt;
            }
        }
        ++stats.atac_fragments_accepted;
        if (fragment.sites.size() > 1) ++stats.atac_fragments_multisite;
        stats.atac_fragment_site_observations += fragment.sites.size();
        stats.accepted_observations += fragment.sites.size();
    };

    bam1_t* record = bam_init1();
    while (sam_itr_next(bam, iterator, record) >= 0) {
        if (!accept_atac_read(opt, record, mt_tid, stats)) continue;
        const char* barcode_ptr = aux_string(record, opt.barcode_tag);
        if (!barcode_ptr || !*barcode_ptr) { ++stats.reject_missing_barcode; continue; }
        const std::string barcode(barcode_ptr);
        const auto cell_it = barcode_to_cell.find(barcode);
        const bool is_empty = empty_barcodes.count(barcode) > 0;
        if (cell_it == barcode_to_cell.end() && !is_empty) {
            ++stats.reject_irrelevant_barcode;
            continue;
        }
        std::unordered_map<int, uint8_t> alleles = extract_atac_site_alleles(
            opt, record, sites, position_to_site, stats);
        if (alleles.empty()) continue;
        const uint16_t flag = record->core.flag;
        const bool paired = flag & BAM_FPAIRED;
        if (!paired || ((flag & BAM_FMUNMAP) || record->core.mtid != mt_tid)) {
            AtacFragmentAccumulator singleton;
            singleton.target_index = is_empty ? -1 : cell_it->second;
            singleton.accepted_reads = 1;
            singleton.alleles = std::move(alleles);
            finalize(singleton);
            continue;
        }
        std::string key;
        const char* qname = bam_get_qname(record);
        key.reserve(barcode.size() + std::strlen(qname) + 1);
        key.append(barcode);
        key.push_back('\x1f');
        key.append(qname);
        AtacFragmentAccumulator& acc = pending[key];
        const int target_index = is_empty ? -1 : cell_it->second;
        if (acc.target_index == -2) acc.target_index = target_index;
        else if (acc.target_index != target_index) {
            throw std::runtime_error("ATAC fragment barcode target changed for QNAME " + std::string(qname));
        }
        ++acc.accepted_reads;
        if (flag & BAM_FREAD1) acc.seen_first = true;
        if (flag & BAM_FREAD2) acc.seen_second = true;
        for (const auto& item : alleles) {
            auto prior = acc.alleles.find(item.first);
            if (prior == acc.alleles.end()) acc.alleles[item.first] = item.second;
            else if (prior->second == item.second) ++stats.atac_overlap_agreements;
            else {
                prior->second = 2;
                ++stats.atac_overlap_conflicts;
                ++stats.conflicting_molecules;
            }
        }
        if ((acc.seen_first && acc.seen_second) || acc.accepted_reads >= 2) {
            finalize(acc);
            pending.erase(key);
        }
    }
    for (auto& item : pending) {
        ++stats.atac_orphan_fragments;
        if (opt.atac_include_singletons) finalize(item.second);
    }
    bam_destroy1(record);
    hts_itr_destroy(iterator);
    hts_idx_destroy(index);
    sam_hdr_destroy(header);
    sam_close(bam);
}

double clamp_probability_local(double value) {
    return std::max(kTiny, std::min(1.0 - kTiny, value));
}

struct AmbientQcEstimate {
    double fraction = std::numeric_limits<double>::quiet_NaN();
    int sites_observed = 0;
    int sites_used = 0;
    uint64_t molecules = 0;
    std::string status = "NOT_REQUESTED";
    std::string reason = "ambient_qc_not_requested";
};

struct AmbientAnchorObservation {
    uint64_t ref = 0;
    uint64_t alt = 0;
    double theta = 0.5;
    double ambient_alt = 0.5;
};

double ambient_anchor_log_likelihood(
        const std::vector<AmbientAnchorObservation>& observations,
        double ambient_fraction) {
    const double c = std::max(0.0, std::min(0.99, ambient_fraction));
    double value = 0.0;
    for (const AmbientAnchorObservation& obs : observations) {
        const double p = clamp_probability_local(
            (1.0 - c) * obs.theta + c * obs.ambient_alt);
        value += static_cast<double>(obs.alt) * std::log(p) +
                 static_cast<double>(obs.ref) * std::log(1.0 - p);
    }
    return value;
}

double ambient_anchor_score(
        const std::vector<AmbientAnchorObservation>& observations,
        double ambient_fraction) {
    const double c = std::max(0.0, std::min(0.99, ambient_fraction));
    double value = 0.0;
    for (const AmbientAnchorObservation& obs : observations) {
        const double delta = obs.ambient_alt - obs.theta;
        if (std::fabs(delta) <= kTiny) continue;
        const double p = clamp_probability_local(
            (1.0 - c) * obs.theta + c * obs.ambient_alt);
        value += delta * (
            static_cast<double>(obs.alt) / p -
            static_cast<double>(obs.ref) / (1.0 - p));
    }
    return value;
}

double maximize_ambient_anchor_fraction(
        const std::vector<AmbientAnchorObservation>& observations) {
    constexpr double lower = 0.0;
    constexpr double upper = 0.99;
    const double lower_score = ambient_anchor_score(observations, lower);
    const double upper_score = ambient_anchor_score(observations, upper);
    if (!std::isfinite(lower_score) || !std::isfinite(upper_score)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (lower_score <= 0.0) return lower;
    if (upper_score >= 0.0) return upper;

    double lo = lower;
    double hi = upper;
    for (int iteration = 0; iteration < 100 && hi - lo > 1e-10; ++iteration) {
        const double mid = 0.5 * (lo + hi);
        const double score = ambient_anchor_score(observations, mid);
        if (!std::isfinite(score)) return std::numeric_limits<double>::quiet_NaN();
        if (score > 0.0) lo = mid;
        else hi = mid;
    }
    const double candidate = 0.5 * (lo + hi);

    // The anchor binomial log likelihood is concave in c. Keep the explicit
    // boundary comparison so finite-precision score behavior reproduces a
    // bounded scalar maximum rather than forcing an interior value.
    double best_c = candidate;
    double best_ll = ambient_anchor_log_likelihood(observations, candidate);
    for (double boundary : {lower, upper}) {
        const double ll = ambient_anchor_log_likelihood(observations, boundary);
        if (ll > best_ll) {
            best_ll = ll;
            best_c = boundary;
        }
    }
    return best_c;
}

AmbientQcEstimate estimate_ambient_qc(
        const Options& opt,
        const CellInfo& cell,
        const PairManifest* pair_manifest,
        const std::vector<AmbientSite>& ambient,
        const double* precomputed_fraction) {
    AmbientQcEstimate out;
    if (!std::isfinite(opt.ambient_qc_max)) return out;

    out.status = "NOT_ESTIMABLE";
    out.reason = "no_usable_ambient_anchor_observations";
    if (!pair_manifest) return out;

    std::vector<AmbientAnchorObservation> observations;
    observations.reserve(static_cast<size_t>(
        std::max(0, pair_manifest->ambient_only_sites_available)));
    bool has_contrast = false;
    for (const auto& item : cell.counts) {
        const int site_index = item.first;
        const AlleleCount& count = item.second;
        if (count.total() == 0) continue;
        const auto manifest_found = pair_manifest->sites.find(site_index);
        if (manifest_found == pair_manifest->sites.end()) continue;
        const ManifestSite& manifest_site = manifest_found->second;
        if (manifest_site.site_class != MANIFEST_AMBIENT_ONLY ||
            manifest_site.canonical_parent1_state != manifest_site.canonical_parent2_state) {
            continue;
        }
        ++out.sites_observed;
        if (site_index < 0 || static_cast<size_t>(site_index) >= ambient.size() ||
            !ambient[site_index].usable ||
            !std::isfinite(ambient[site_index].alt_fraction)) {
            continue;
        }
        AmbientAnchorObservation obs;
        obs.ref = count.ref;
        obs.alt = count.alt;
        obs.theta = manifest_site.canonical_parent1_state == 1
            ? 1.0 - opt.error_rate : opt.error_rate;
        obs.ambient_alt = clamp_probability_local(ambient[site_index].alt_fraction);
        if (std::fabs(obs.ambient_alt - obs.theta) > kTiny) has_contrast = true;
        observations.push_back(obs);
        ++out.sites_used;
        out.molecules += count.total();
    }

    if (precomputed_fraction && std::isfinite(*precomputed_fraction)) {
        out.fraction = std::max(0.0, std::min(0.99, *precomputed_fraction));
        out.status = out.fraction > opt.ambient_qc_max ? "HIGH_MT_AMBIENT" : "PASS";
        out.reason = out.status == "PASS"
            ? "ambient_fraction_within_qc_max_precomputed"
            : "ambient_fraction_exceeds_qc_max";
        return out;
    }

    if (observations.empty() || out.molecules == 0 || !has_contrast) {
        if (!has_contrast && !observations.empty()) {
            out.reason = "ambient_anchor_profile_has_no_contrast";
        }
        return out;
    }

    out.fraction = maximize_ambient_anchor_fraction(observations);
    if (!std::isfinite(out.fraction)) {
        out.reason = "ambient_anchor_optimization_failed";
        return out;
    }
    out.status = out.fraction > opt.ambient_qc_max ? "HIGH_MT_AMBIENT" : "PASS";
    out.reason = out.status == "PASS"
        ? "ambient_fraction_within_qc_max"
        : "ambient_fraction_exceeds_qc_max";
    return out;
}

std::string format_double(double value, int precision = 8) {
    if (!std::isfinite(value)) return "NA";
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}

#ifndef CELLBOUNCER_SOURCE_REVISION
#define CELLBOUNCER_SOURCE_REVISION "unknown"
#endif

constexpr const char* kMtModelVersion = "mt_fusion_ratio_full_v5_anchor_ambient_qc";

std::string assay_estimand(const std::string& assay_mode) {
    if (assay_mode == "ATAC") return "mtDNA_parental_fragment_fraction";
    if (assay_mode == "GENERIC") return "parental_mitochondrial_alignment_fraction";
    return "mtRNA_parental_molecule_fraction";
}

struct CalibrationProbabilityResult {
    double parent1_alt_probability = std::numeric_limits<double>::quiet_NaN();
    double parent2_alt_probability = std::numeric_limits<double>::quiet_NaN();
    bool parent1_calibrated = false;
    bool parent2_calibrated = false;
    std::string status = "NOT_REQUESTED";
};

CalibrationProbabilityResult calibration_probabilities(
        const Options& opt,
        const SiteCalibrationData& calibration,
        const PairManifest& pair_manifest,
        int site_index,
        bool assignment_matches_canonical,
        int assignment_parent1_state,
        int assignment_parent2_state) {
    CalibrationProbabilityResult out;
    out.parent1_alt_probability = assignment_parent1_state == 1
        ? 1.0 - opt.error_rate : opt.error_rate;
    out.parent2_alt_probability = assignment_parent2_state == 1
        ? 1.0 - opt.error_rate : opt.error_rate;
    if (opt.site_calibration.empty()) return out;
    out.status = "FALLBACK";
    const auto pair_it = calibration.by_pair.find(
        pair_manifest_key(pair_manifest.canonical_parent1, pair_manifest.canonical_parent2));
    if (pair_it == calibration.by_pair.end()) return out;
    const auto site_it = pair_it->second.find(site_index);
    if (site_it == pair_it->second.end()) return out;
    const SiteCalibrationRecord& record = site_it->second;
    if (assignment_matches_canonical) {
        if (record.has_parent1) {
            out.parent1_alt_probability = record.parent1_alt_probability;
            out.parent1_calibrated = true;
        }
        if (record.has_parent2) {
            out.parent2_alt_probability = record.parent2_alt_probability;
            out.parent2_calibrated = true;
        }
    } else {
        if (record.has_parent2) {
            out.parent1_alt_probability = record.parent2_alt_probability;
            out.parent1_calibrated = true;
        }
        if (record.has_parent1) {
            out.parent2_alt_probability = record.parent1_alt_probability;
            out.parent2_calibrated = true;
        }
    }
    if (out.parent1_calibrated && out.parent2_calibrated) out.status = "FULL";
    else if (out.parent1_calibrated || out.parent2_calibrated) out.status = "PARTIAL";
    return out;
}

MtRatioFitConfig make_fit_config(const Options& opt,
                                 double pooled_rho,
                                 uint64_t ratio_molecules) {
    MtRatioFitConfig config;
    config.use_beta_binomial = opt.assay_mode != "ATAC" && opt.likelihood == "beta_binomial";
    config.ambient_enabled = false;
    config.ambient_fixed = false;
    config.fixed_ambient = 0.0;
    config.overdispersion_initial = opt.overdispersion_initial;
    config.overdispersion_max = opt.overdispersion_max;
    if (config.use_beta_binomial && opt.rho_mode == "fixed") {
        config.overdispersion_fixed = true;
        config.fixed_overdispersion = pooled_rho;
    } else if (config.use_beta_binomial && opt.rho_mode == "low_information_fixed" &&
               ratio_molecules < static_cast<uint64_t>(opt.rho_low_information_molecules)) {
        config.overdispersion_fixed = true;
        config.fixed_overdispersion = pooled_rho;
    } else if (config.use_beta_binomial && opt.rho_mode == "shrink") {
        config.overdispersion_penalized = true;
        config.overdispersion_prior_mean = pooled_rho;
        config.overdispersion_prior_strength = opt.rho_prior_strength;
        config.overdispersion_initial = pooled_rho;
    }
    config.max_iterations = opt.max_iterations;
    config.tolerance = opt.tolerance;
    return config;
}



std::unordered_map<int, double> compute_rna_site_influence(
        const std::vector<MtRatioObservation>& observations,
        const MtRatioFitConfig& config,
        const MtRatioFitResult& full_fit) {
    std::unordered_map<int, double> result;
    if (!full_fit.converged || !std::isfinite(full_fit.parent2_fraction)) return result;
    for (const MtRatioObservation& target : observations) {
        std::vector<MtRatioObservation> reduced;
        reduced.reserve(observations.size() > 0 ? observations.size() - 1 : 0);
        for (const MtRatioObservation& obs : observations) {
            if (obs.site_index != target.site_index) reduced.push_back(obs);
        }
        if (reduced.empty()) {
            result[target.site_index] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }
        MtRatioFitResult loo = fit_mt_ratio(reduced, config);
        result[target.site_index] = loo.converged && std::isfinite(loo.parent2_fraction)
            ? loo.parent2_fraction - full_fit.parent2_fraction
            : std::numeric_limits<double>::quiet_NaN();
    }
    return result;
}

std::unordered_map<int, double> compute_atac_site_influence(
        const std::vector<MtRatioFragmentObservation>& fragments,
        const MtRatioFitConfig& config,
        const MtRatioFitResult& full_fit,
        const std::unordered_set<int>& used_sites) {
    std::unordered_map<int, double> result;
    if (!full_fit.converged || !std::isfinite(full_fit.parent2_fraction)) return result;
    for (int site_index : used_sites) {
        std::vector<MtRatioFragmentObservation> reduced;
        reduced.reserve(fragments.size());
        for (const MtRatioFragmentObservation& fragment : fragments) {
            MtRatioFragmentObservation copy;
            copy.sites.reserve(fragment.sites.size());
            for (const MtRatioFragmentSite& site : fragment.sites) {
                if (site.site_index != site_index) copy.sites.push_back(site);
            }
            if (!copy.sites.empty()) reduced.push_back(std::move(copy));
        }
        if (reduced.empty()) {
            result[site_index] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }
        MtRatioFitResult loo = fit_mt_fragment_ratio(reduced, config);
        result[site_index] = loo.converged && std::isfinite(loo.parent2_fraction)
            ? loo.parent2_fraction - full_fit.parent2_fraction
            : std::numeric_limits<double>::quiet_NaN();
    }
    return result;
}

double site_deviance_residual(const MtRatioObservation& obs,
                              double predicted,
                              double rho,
                              bool use_beta) {
    const uint64_t total = obs.ref + obs.alt;
    if (total == 0) return std::numeric_limits<double>::quiet_NaN();
    const double observed = static_cast<double>(obs.alt) / static_cast<double>(total);
    const double fitted_ll = mt_ratio_count_log_likelihood(
        obs.ref, obs.alt, predicted, rho, use_beta);
    const double saturated_ll = mt_ratio_count_log_likelihood(
        obs.ref, obs.alt, std::max(kTiny, std::min(1.0 - kTiny, observed)), rho, use_beta);
    const double dev = std::max(0.0, 2.0 * (saturated_ll - fitted_ll));
    const double sign = observed >= predicted ? 1.0 : -1.0;
    return sign * std::sqrt(dev);
}

std::string calibration_ids_string(const SiteCalibrationData& data) {
    if (data.calibration_ids.empty()) return "NA";
    std::ostringstream out;
    bool first = true;
    for (const std::string& id : data.calibration_ids) {
        if (!first) out << ',';
        out << id;
        first = false;
    }
    return out.str();
}

MtRatioProfileResult canonicalize_profile(MtRatioProfileResult profile,
                                          bool assignment_matches_canonical) {
    if (assignment_matches_canonical) return profile;
    const double old_low = profile.parent2_ci_low;
    const double old_high = profile.parent2_ci_high;
    if (std::isfinite(old_low) && std::isfinite(old_high)) {
        profile.parent2_ci_low = 1.0 - old_high;
        profile.parent2_ci_high = 1.0 - old_low;
        profile.parent1_ci_low = 1.0 - profile.parent2_ci_high;
        profile.parent1_ci_high = 1.0 - profile.parent2_ci_low;
    }
    std::swap(profile.delta_ll_parent1_only, profile.delta_ll_parent2_only);
    if (profile.inheritance_class == "ONLY_PARENT1") {
        profile.inheritance_class = "ONLY_PARENT2";
    } else if (profile.inheritance_class == "ONLY_PARENT2") {
        profile.inheritance_class = "ONLY_PARENT1";
    }
    if (profile.inheritance_class_reason == "PROFILE_INTERVAL_WITHIN_PARENT1_REGION") {
        profile.inheritance_class_reason = "PROFILE_INTERVAL_WITHIN_PARENT2_REGION";
    } else if (profile.inheritance_class_reason == "PROFILE_INTERVAL_WITHIN_PARENT2_REGION") {
        profile.inheritance_class_reason = "PROFILE_INTERVAL_WITHIN_PARENT1_REGION";
    }
    if (!profile.grid_delta_log_likelihood.empty()) {
        std::reverse(profile.grid_delta_log_likelihood.begin(),
                     profile.grid_delta_log_likelihood.end());
    }
    return profile;
}

void gz_write_text(gzFile out, const std::string& text) {
    if (!out) return;
    const int written = gzwrite(out, text.data(), static_cast<unsigned int>(text.size()));
    if (written != static_cast<int>(text.size())) {
        throw std::runtime_error("Failed while writing compressed mt profile output");
    }
}

std::string profile_delta_csv(const std::vector<double>& values) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) out << ',';
        if (std::isfinite(values[i])) out << values[i];
        else out << "NA";
    }
    return out.str();
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options opt = parse_options(argc, argv);

        std::vector<std::string> sample_names;
        std::unordered_map<std::string, int> sample_index;
        std::unordered_map<int64_t, int> position_to_site;
        std::vector<MtSite> sites = load_mt_sites(
            opt, sample_names, sample_index, position_to_site);
        const std::unordered_map<std::string, PairManifest> pair_manifests =
            load_site_manifest(opt, sites, position_to_site);

        std::unordered_map<std::string, int> barcode_to_cell;
        std::vector<CellInfo> cells = load_assignments(
            opt.assignments, opt.calibration_reference, opt.library_id,
            sample_index, barcode_to_cell);
        const std::unordered_map<std::string, PairManifest> derived_pair_manifests =
            derive_missing_pair_manifests(
                opt, cells, pair_manifests, sites, position_to_site);
        std::unordered_set<std::string> empty_barcodes = load_empty_barcodes(
            opt.empty_barcodes, barcode_to_cell);
        std::unordered_map<std::string, double> fixed_ambient =
            load_fixed_ambient_fractions(opt.ambient_fraction_file);
        const SiteCalibrationData site_calibration =
            load_site_calibration(opt, position_to_site);
        const std::vector<RhoReferenceEntry> rho_reference = load_rho_reference(opt);
        const std::string panel_fingerprint = file_fingerprint(opt.vcf);
        const std::string manifest_fingerprint = opt.legacy_panel_gt
            ? "NA" : file_fingerprint(opt.site_manifest);
        const std::string calibration_fingerprint = opt.site_calibration.empty()
            ? "NA" : file_fingerprint(opt.site_calibration);
        const bool ambient_background_available = !opt.ambient_none;
        const bool ambient_qc_requested = std::isfinite(opt.ambient_qc_max);
        const std::string ambient_mode = opt.ambient_none ? "NONE" :
            (!opt.empty_barcodes.empty() ? "QC_EMPTY_BARCODES" : "QC_PROFILE_FILE");
        const std::string ambient_profile_id = !opt.ambient_profile.empty()
            ? path_basename(opt.ambient_profile)
            : (!opt.empty_barcodes.empty() ? "INLINE_EMPTY_BARCODES" : "NA");
        const std::string ambient_profile_fingerprint = !opt.ambient_profile.empty()
            ? file_fingerprint(opt.ambient_profile) : "NA";
        const std::string ambient_fraction_fingerprint = !opt.ambient_fraction_file.empty()
            ? file_fingerprint(opt.ambient_fraction_file) : "NA";
        const std::string experimental_mode =
            !opt.site_calibration.empty() && ambient_qc_requested
                ? "SITE_CALIBRATION_PLUS_AMBIENT_QC"
                : (!opt.site_calibration.empty() ? "SITE_CALIBRATION"
                   : (ambient_qc_requested ? "AMBIENT_QC" : "BASELINE"));

        std::vector<AmbientSite> ambient(sites.size());
        if (!opt.ambient_profile.empty()) {
            load_ambient_profile(opt.ambient_profile, opt.mito_chrom, sites,
                                 position_to_site, opt.ambient_min_molecules, ambient);
        }

        std::fprintf(stderr, "Mitochondrial panel sites: %zu\n", sites.size());
        std::fprintf(stderr, "Fusion cells in assignments: %zu\n", cells.size());
        if (!opt.legacy_panel_gt) {
            std::fprintf(stderr, "Library-specific manifest pairs for %s: %zu\n",
                         opt.library_id.c_str(), pair_manifests.size());
        } else {
            std::fprintf(stderr,
                "WARNING: legacy GT-only panel mode is diagnostic and has no ambient-only anchors\n");
        }
        if (!empty_barcodes.empty()) {
            std::fprintf(stderr, "Empty-droplet barcodes: %zu\n", empty_barcodes.size());
        }
        if (!fixed_ambient.empty()) {
            std::fprintf(stderr,
                "Precomputed ambient QC fractions loaded: %zu (optional per-cell QC override)\n",
                fixed_ambient.size());
        }

        ReadStats stats;
        if (opt.assay_mode == "ATAC") {
            count_atac_fragments(opt, sites, position_to_site, barcode_to_cell,
                                 empty_barcodes, cells, ambient, stats);
        } else {
            count_molecules(opt, sites, position_to_site, barcode_to_cell,
                            empty_barcodes, cells, ambient, stats);
        }
        if (stats.pileup_bases == 0) {
            throw std::runtime_error(
                "The BAM contains the mitochondrial contig but no reads overlap the selected "
                "mitochondrial panel. This commonly means a nuclear-only variant-site subset BAM "
                "was supplied; rebuild the subset BED with the emitted mt sites BED.");
        }
        if (stats.accepted_observations == 0) {
            throw std::runtime_error(
                opt.assay_mode == "ATAC"
                    ? "No mitochondrial ATAC fragment/site observations were accepted. Check paired-read flags, CB tags, mapping filters, and panel overlap."
                    : "No UMI-collapsed mitochondrial molecule/site observations were accepted. Check CB/UB tags, barcode agreement, mapping filters, and panel overlap.");
        }

        if (!empty_barcodes.empty()) {
            for (AmbientSite& a : ambient) {
                const uint64_t total = a.total();
                a.alt_fraction = (static_cast<double>(a.alt) + 0.5) /
                                 (static_cast<double>(total) + 1.0);
                a.usable = static_cast<int>(total) >= opt.ambient_min_molecules;
            }
        }

        // The selected experimental architecture uses MT ambient/background only
        // as an independent anchor-derived QC measurement. It is never enabled
        // in the parental-ratio likelihood.
        int usable_ambient_sites = 0;
        if (ambient_background_available) {
            for (const AmbientSite& a : ambient) if (a.usable) ++usable_ambient_sites;
            if (ambient_qc_requested && usable_ambient_sites == 0) {
                throw std::runtime_error(
                    "No mt-panel sites have a usable ambient profile; provide more empty-droplet "
                    "evidence or lower --ambient_min_molecules for --ambient_qc_max");
            }
        }

        const std::string ambient_path = opt.output_prefix + ".mt_ambient_profile.tsv";
        std::ofstream ambient_out(ambient_path);
        if (!ambient_out) throw std::runtime_error("Could not write: " + ambient_path);
        ambient_out << "chrom\tpos\tref\talt\tref_molecules\talt_molecules\talt_fraction\tstatus\n";
        for (size_t i = 0; i < sites.size(); ++i) {
            ambient_out << opt.mito_chrom << '\t' << sites[i].pos + 1 << '\t'
                        << sites[i].ref << '\t' << sites[i].alt << '\t'
                        << ambient[i].ref << '\t' << ambient[i].alt << '\t'
                        << format_double(ambient[i].alt_fraction) << '\t'
                        << (!ambient_background_available ? "NOT_USED" :
                            (ambient[i].usable ? "PASS" : "LOW_COUNTS")) << '\n';
        }

        const std::string ratio_path = opt.output_prefix + ".mt_ratio.tsv";
        std::ofstream ratio_out(ratio_path);
        if (!ratio_out) throw std::runtime_error("Could not write: " + ratio_path);
        ratio_out << "library_id\tbarcode\tidentity\tparent1\tparent2\tlikelihood"
                  << "\tparent1_fraction\tparent2_fraction\tratio_se\tratio_se_robust"
                  << "\tparent1_ci_low\tparent1_ci_high\tambient_fraction\tambient_se"
                  << "\tambient_se_robust\toverdispersion_rho\toverdispersion_se"
                  << "\tratio_sites_available\tambient_only_sites_available"
                  << "\tratio_sites_observed\tambient_only_sites_observed"
                  << "\tratio_sites_used\tambient_only_sites_used"
                  << "\tratio_molecules\ttotal_molecules_used\tsites_used"
                  << "\traw_parent1_support\traw_parent2_support\traw_parent1_fraction"
                  << "\tlog_likelihood\tinformation_condition\tmin_information_eigenvalue"
                  << "\titerations\tconverged\tfit_mode\tstatus"
                  << "\tsource_revision\tmodel_version\tassay_mode\testimand"
                  << "\tcanonical_parent1\tcanonical_parent2\tassignment_matches_canonical"
                  << "\tcanonical_parent1_fraction\tcanonical_parent2_fraction"
                  << "\tparent1_profile_ci_low\tparent1_profile_ci_high"
                  << "\tparent2_profile_ci_low\tparent2_profile_ci_high"
                  << "\tprofile_ci_level\tprofile_status\tprofile_evaluations"
                  << "\tprofile_failed_evaluations\tsingle_parent_epsilon"
                  << "\tinheritance_class\tinheritance_class_reason"
                  << "\tdelta_ll_parent1_only\tdelta_ll_both\tdelta_ll_parent2_only"
                  << "\tpanel_fingerprint\tmanifest_fingerprint\tsite_calibration_id"
                  << "\tsite_calibration_fingerprint\tsite_calibration_stratum"
                  << "\tsite_calibration_fraction\tsite_calibration_status"
                  << "\tuncalibrated_canonical_parent2_fraction\tuncalibrated_log_likelihood"
                  << "\tcalibration_shift_parent2\tmax_abs_site_influence"
                  << "\tmost_influential_site\tfraction_support_from_top_site"
                  << "\trho_mode\trho_mode_effective\tpooled_rho_used\trho_prior_strength"
                  << "\tfit_objective\tcalibration_reference_mode\tcalibration_true_parent"
                  << "\tcalibration_original_identity\tatac_fragments_used"
                  << "\tatac_multisite_fragments_used\tatac_fragment_site_observations_used"
                  << "\tpair_definition_source\texperimental_mode\tambient_mode"
                  << "\tambient_profile_id\tambient_profile_fingerprint"
                  << "\tambient_fraction_file_fingerprint"
                  << "\tambient_qc_max\tambient_qc_status\tambient_qc_reason"
                  << "\tambient_anchor_molecules\n";

        gzFile profile_out = nullptr;
        const std::string profile_path = opt.output_prefix + ".mt_profile.tsv.gz";
        if (opt.write_profile_grid) {
            profile_out = gzopen(profile_path.c_str(), "wb");
            if (!profile_out) throw std::runtime_error("Could not write: " + profile_path);
            gz_write_text(profile_out,
                "library_id\tbarcode\tidentity\tcanonical_parent1\tcanonical_parent2"
                "\tcanonical_parent2_fraction_mle\tgrid_start\tgrid_step\tgrid_points"
                "\tdelta_log_likelihood_csv\tprofile_status\tprofile_objective"
                "\tsite_calibration_status\trho_mode\n");
        }

        std::ofstream site_out;
        if (opt.write_site_counts) {
            const std::string site_path = opt.output_prefix + ".mt_site_counts.tsv";
            site_out.open(site_path);
            if (!site_out) throw std::runtime_error("Could not write: " + site_path);
            site_out << "library_id\tbarcode\tidentity\tchrom\tpos\tref\talt\tsite_class"
                     << "\tparent1_state\tparent2_state\tratio_informative\tambient_anchor"
                     << "\tref_molecules\talt_molecules\tambient_ref_molecules"
                     << "\tambient_alt_molecules\tambient_alt_fraction\tused_in_fit"
                     << "\tparent1_alt_probability\tparent2_alt_probability"
                     << "\tpredicted_alt_probability\tsite_log_likelihood"
                     << "\tsite_deviance_residual\tleave_one_site_out_delta_r"
                     << "\tsite_influence_status\tsite_calibration_status"
                     << "\tcalibration_true_parent\tassay_mode"
                     << "\tcanonical_parent1\tcanonical_parent2\n";
        }

        struct CellCounters {
            uint64_t pass_cells = 0;
            uint64_t weak_cells = 0;
            uint64_t low_cells = 0;
            uint64_t fixed_cells = 0;
            uint64_t pair_missing_cells = 0;
            uint64_t derived_pair_cells = 0;
            uint64_t collapsed_pair_cells = 0;
            uint64_t unanchored_cells = 0;
            uint64_t profile_attempted_cells = 0;
            uint64_t profile_pass_cells = 0;
            uint64_t profile_partial_cells = 0;
            uint64_t profile_failed_cells = 0;
            uint64_t profile_multimodal_cells = 0;
            uint64_t inheritance_parent1_cells = 0;
            uint64_t inheritance_both_cells = 0;
            uint64_t inheritance_parent2_cells = 0;
            uint64_t inheritance_ambiguous_cells = 0;
            uint64_t calibration_full_cells = 0;
            uint64_t calibration_partial_cells = 0;
            uint64_t calibration_fallback_cells = 0;
            uint64_t calibration_full_site_uses = 0;
            uint64_t calibration_partial_site_uses = 0;
            uint64_t calibration_fallback_site_uses = 0;
            uint64_t rho_fixed_cells = 0;
            uint64_t rho_shrunk_cells = 0;
            uint64_t high_mt_ambient_cells = 0;
            uint64_t ambient_qc_not_estimable_cells = 0;
        };

        struct CellResult {
            std::string ratio_row;
            std::string profile_row;
            std::string site_rows;
            CellCounters counters;
            std::exception_ptr error;
        };
        uint64_t pass_cells = 0;
        uint64_t weak_cells = 0;
        uint64_t low_cells = 0;
        uint64_t fixed_cells = 0;
        uint64_t pair_missing_cells = 0;
        uint64_t derived_pair_cells = 0;
        uint64_t collapsed_pair_cells = 0;
        uint64_t unanchored_cells = 0;
        uint64_t profile_attempted_cells = 0;
        uint64_t profile_pass_cells = 0;
        uint64_t profile_partial_cells = 0;
        uint64_t profile_failed_cells = 0;
        uint64_t profile_multimodal_cells = 0;
        uint64_t inheritance_parent1_cells = 0;
        uint64_t inheritance_both_cells = 0;
        uint64_t inheritance_parent2_cells = 0;
        uint64_t inheritance_ambiguous_cells = 0;
        uint64_t calibration_full_cells = 0;
        uint64_t calibration_partial_cells = 0;
        uint64_t calibration_fallback_cells = 0;
        uint64_t calibration_full_site_uses = 0;
        uint64_t calibration_partial_site_uses = 0;
        uint64_t calibration_fallback_site_uses = 0;
        uint64_t rho_fixed_cells = 0;
        uint64_t rho_shrunk_cells = 0;
        uint64_t high_mt_ambient_cells = 0;
        uint64_t ambient_qc_not_estimable_cells = 0;

        std::vector<CellResult> cell_results(cells.size());

        auto process_cell = [&](size_t cell_index) {
            CellResult& cell_result = cell_results[cell_index];
            try {
                const CellInfo& cell = cells[cell_index];
                std::ostringstream ratio_out;
                std::ostringstream site_out;
            auto& pass_cells = cell_result.counters.pass_cells;
            auto& weak_cells = cell_result.counters.weak_cells;
            auto& low_cells = cell_result.counters.low_cells;
            auto& fixed_cells = cell_result.counters.fixed_cells;
            auto& pair_missing_cells = cell_result.counters.pair_missing_cells;
            auto& derived_pair_cells = cell_result.counters.derived_pair_cells;
            auto& collapsed_pair_cells = cell_result.counters.collapsed_pair_cells;
            auto& profile_attempted_cells = cell_result.counters.profile_attempted_cells;
            auto& profile_pass_cells = cell_result.counters.profile_pass_cells;
            auto& profile_partial_cells = cell_result.counters.profile_partial_cells;
            auto& profile_failed_cells = cell_result.counters.profile_failed_cells;
            auto& profile_multimodal_cells = cell_result.counters.profile_multimodal_cells;
            auto& inheritance_parent1_cells = cell_result.counters.inheritance_parent1_cells;
            auto& inheritance_both_cells = cell_result.counters.inheritance_both_cells;
            auto& inheritance_parent2_cells = cell_result.counters.inheritance_parent2_cells;
            auto& inheritance_ambiguous_cells = cell_result.counters.inheritance_ambiguous_cells;
            auto& calibration_full_cells = cell_result.counters.calibration_full_cells;
            auto& calibration_partial_cells = cell_result.counters.calibration_partial_cells;
            auto& calibration_fallback_cells = cell_result.counters.calibration_fallback_cells;
            auto& calibration_full_site_uses = cell_result.counters.calibration_full_site_uses;
            auto& calibration_partial_site_uses = cell_result.counters.calibration_partial_site_uses;
            auto& calibration_fallback_site_uses = cell_result.counters.calibration_fallback_site_uses;
            auto& rho_fixed_cells = cell_result.counters.rho_fixed_cells;
            auto& rho_shrunk_cells = cell_result.counters.rho_shrunk_cells;
            auto& high_mt_ambient_cells = cell_result.counters.high_mt_ambient_cells;
            auto& ambient_qc_not_estimable_cells = cell_result.counters.ambient_qc_not_estimable_cells;
            std::string status = "PASS";
            std::string fit_mode;
            if (opt.legacy_panel_gt) fit_mode = "LEGACY_GT_";
            if (opt.assay_mode == "ATAC") {
                fit_mode += "ATAC_FRAGMENT_";
                fit_mode += "NO_AMBIENT";
            } else {
                fit_mode += opt.likelihood == "beta_binomial" ? "BETA_BINOMIAL_" : "BINOMIAL_";
                fit_mode += "NO_AMBIENT";
            }
            if (ambient_qc_requested) fit_mode += "_AMBIENT_QC";

            const bool same_parent = cell.parent1_index == cell.parent2_index;
            const PairManifest* pair_manifest = nullptr;
            std::string pair_definition_source = same_parent ? "NOT_APPLICABLE" : "UNAVAILABLE";
            PairManifest legacy_manifest;
            if (!same_parent) {
                if (opt.legacy_panel_gt) {
                    pair_definition_source = "LEGACY_GT_PANEL";
                    const std::pair<std::string, std::string> canonical =
                        canonical_parents(cell.parent1, cell.parent2);
                    legacy_manifest.canonical_parent1 = canonical.first;
                    legacy_manifest.canonical_parent2 = canonical.second;
                    for (size_t site_index = 0; site_index < sites.size(); ++site_index) {
                        int s1 = sites[site_index].genotype_state[cell.parent1_index];
                        int s2 = sites[site_index].genotype_state[cell.parent2_index];
                        if (s1 < 0 || s2 < 0 || s1 == s2) continue;
                        ManifestSite entry;
                        entry.site_index = static_cast<int>(site_index);
                        if (cell.parent1 == canonical.first) {
                            entry.canonical_parent1_state = static_cast<int8_t>(s1);
                            entry.canonical_parent2_state = static_cast<int8_t>(s2);
                        } else {
                            entry.canonical_parent1_state = static_cast<int8_t>(s2);
                            entry.canonical_parent2_state = static_cast<int8_t>(s1);
                        }
                        entry.site_class = MANIFEST_LEGACY_RATIO;
                        legacy_manifest.sites[static_cast<int>(site_index)] = entry;
                        ++legacy_manifest.ratio_sites_available;
                    }
                    pair_manifest = &legacy_manifest;
                } else {
                    const std::string pair_key = pair_manifest_key(cell.parent1, cell.parent2);
                    const auto found = pair_manifests.find(pair_key);
                    if (found != pair_manifests.end()) {
                        pair_manifest = &found->second;
                        pair_definition_source = "LIBRARY_SITE_MANIFEST";
                    } else {
                        const auto derived_found = derived_pair_manifests.find(pair_key);
                        if (derived_found != derived_pair_manifests.end()) {
                            pair_manifest = &derived_found->second;
                            pair_definition_source = "RECONCILED_PAIR_FROM_PANEL";
                            ++derived_pair_cells;
                        }
                    }
                }
            }

            const auto fixed_it = fixed_ambient.find(cell.barcode);
            const double* fixed_value = fixed_it == fixed_ambient.end() ? nullptr : &fixed_it->second;
            if (fixed_value && ambient_qc_requested) ++fixed_cells;

            int ratio_sites_available = pair_manifest ? pair_manifest->ratio_sites_available : 0;
            int ambient_sites_available = pair_manifest ? pair_manifest->ambient_only_sites_available : 0;
            int ratio_sites_observed = 0;
            int ambient_sites_observed = 0;
            int ratio_sites_used = 0;
            int ambient_sites_used = 0;
            uint64_t ratio_molecules = 0;
            uint64_t total_molecules_used = 0;
            uint64_t support_parent1 = 0;
            uint64_t support_parent2 = 0;
            std::vector<MtRatioObservation> observations;
            std::vector<MtRatioObservation> uncalibrated_observations;
            std::vector<MtRatioFragmentObservation> fragment_observations;
            std::vector<MtRatioFragmentObservation> uncalibrated_fragment_observations;
            std::unordered_set<int> used_site_indices;
            uint64_t fully_calibrated_used_sites = 0;
            uint64_t partially_calibrated_used_sites = 0;
            uint64_t atac_multisite_fragments_used = 0;
            uint64_t atac_fragment_site_observations_used = 0;

            if (same_parent) {
                status = "SAME_PARENT_UNIDENTIFIABLE";
            } else if (!pair_manifest) {
                status = "PAIR_NOT_IN_LIBRARY_MANIFEST";
                ++pair_missing_cells;
            } else if (ratio_sites_available == 0) {
                status = "MT_HAPLOTYPE_COLLAPSE";
                ++collapsed_pair_cells;
            }

            bool assignment_matches_canonical = false;
            if (pair_manifest) {
                assignment_matches_canonical =
                    cell.parent1 == pair_manifest->canonical_parent1 &&
                    cell.parent2 == pair_manifest->canonical_parent2;
                const bool assignment_reversed =
                    cell.parent2 == pair_manifest->canonical_parent1 &&
                    cell.parent1 == pair_manifest->canonical_parent2;
                if (!assignment_matches_canonical && !assignment_reversed) {
                    throw std::runtime_error("Internal pair-manifest orientation error for " +
                                             cell.identity);
                }
            }

            const std::pair<std::string, std::string> fallback_canonical =
                canonical_parents(cell.parent1, cell.parent2);
            const std::string canonical_parent1 = pair_manifest
                ? pair_manifest->canonical_parent1 : fallback_canonical.first;
            const std::string canonical_parent2 = pair_manifest
                ? pair_manifest->canonical_parent2 : fallback_canonical.second;

            if (pair_manifest) {
                for (const auto& item : cell.counts) {
                    const int site_index = item.first;
                    const AlleleCount& count = item.second;
                    const auto manifest_found = pair_manifest->sites.find(site_index);
                    if (manifest_found == pair_manifest->sites.end()) continue;
                    const ManifestSite& manifest_site = manifest_found->second;
                    int p1_state = manifest_site.canonical_parent1_state;
                    int p2_state = manifest_site.canonical_parent2_state;
                    if (!assignment_matches_canonical) std::swap(p1_state, p2_state);
                    const bool ratio_informative = p1_state != p2_state;
                    const bool ambient_anchor = manifest_site.site_class == MANIFEST_AMBIENT_ONLY;
                    if (count.total() == 0) continue;
                    if (ratio_informative) {
                        ++ratio_sites_observed;
                        if (opt.assay_mode != "ATAC") ratio_molecules += count.total();
                        if (p1_state == 0) {
                            support_parent1 += count.ref;
                            support_parent2 += count.alt;
                        } else {
                            support_parent1 += count.alt;
                            support_parent2 += count.ref;
                        }
                    } else if (ambient_anchor) {
                        ++ambient_sites_observed;
                    }

                    if (!ratio_informative) continue;
                    const CalibrationProbabilityResult probabilities = calibration_probabilities(
                        opt, site_calibration, *pair_manifest, site_index,
                        assignment_matches_canonical, p1_state, p2_state);
                    MtRatioObservation obs;
                    obs.site_index = site_index;
                    obs.ref = count.ref;
                    obs.alt = count.alt;
                    obs.parent1_alt_probability = probabilities.parent1_alt_probability;
                    obs.parent2_alt_probability = probabilities.parent2_alt_probability;
                    obs.ambient_alt_probability = 0.5;
                    obs.ratio_informative = ratio_informative;
                    obs.ambient_anchor = ambient_anchor;
                    observations.push_back(obs);

                    MtRatioObservation uncalibrated = obs;
                    uncalibrated.parent1_alt_probability = p1_state == 1
                        ? 1.0 - opt.error_rate : opt.error_rate;
                    uncalibrated.parent2_alt_probability = p2_state == 1
                        ? 1.0 - opt.error_rate : opt.error_rate;
                    uncalibrated_observations.push_back(uncalibrated);
                    used_site_indices.insert(site_index);
                    if (probabilities.parent1_calibrated && probabilities.parent2_calibrated) {
                        ++fully_calibrated_used_sites;
                    } else if (probabilities.parent1_calibrated || probabilities.parent2_calibrated) {
                        ++partially_calibrated_used_sites;
                    }
                    if (opt.assay_mode != "ATAC") total_molecules_used += count.total();
                    if (ratio_informative) ++ratio_sites_used;
                    if (ambient_anchor) ++ambient_sites_used;
                }

                if (opt.assay_mode == "ATAC") {
                    ratio_molecules = 0;
                    total_molecules_used = 0;
                    for (const RawFragmentObservation& raw_fragment : cell.fragments) {
                        MtRatioFragmentObservation fragment;
                        MtRatioFragmentObservation uncalibrated_fragment;
                        bool has_ratio_evidence = false;
                        for (const RawFragmentSite& raw_site : raw_fragment.sites) {
                            const int site_index = raw_site.site_index;
                            const auto manifest_found = pair_manifest->sites.find(site_index);
                            if (manifest_found == pair_manifest->sites.end()) continue;
                            const ManifestSite& manifest_site = manifest_found->second;
                            int p1_state = manifest_site.canonical_parent1_state;
                            int p2_state = manifest_site.canonical_parent2_state;
                            if (!assignment_matches_canonical) std::swap(p1_state, p2_state);
                            const bool ratio_informative = p1_state != p2_state;
                            const bool ambient_anchor = manifest_site.site_class == MANIFEST_AMBIENT_ONLY;
                            if (!ratio_informative) continue;
                            const CalibrationProbabilityResult probabilities = calibration_probabilities(
                                opt, site_calibration, *pair_manifest, site_index,
                                assignment_matches_canonical, p1_state, p2_state);
                            MtRatioFragmentSite fragment_site;
                            fragment_site.site_index = site_index;
                            fragment_site.allele = raw_site.allele;
                            fragment_site.parent1_alt_probability = probabilities.parent1_alt_probability;
                            fragment_site.parent2_alt_probability = probabilities.parent2_alt_probability;
                            fragment_site.ambient_alt_probability = 0.5;
                            fragment_site.ratio_informative = ratio_informative;
                            fragment_site.ambient_anchor = ambient_anchor;
                            fragment.sites.push_back(fragment_site);

                            MtRatioFragmentSite uncalibrated_site = fragment_site;
                            uncalibrated_site.parent1_alt_probability = p1_state == 1
                                ? 1.0 - opt.error_rate : opt.error_rate;
                            uncalibrated_site.parent2_alt_probability = p2_state == 1
                                ? 1.0 - opt.error_rate : opt.error_rate;
                            uncalibrated_fragment.sites.push_back(uncalibrated_site);
                            if (ratio_informative) has_ratio_evidence = true;
                        }
                        if (fragment.sites.empty()) continue;
                        if (has_ratio_evidence) ++ratio_molecules;
                        ++total_molecules_used;
                        atac_fragment_site_observations_used += fragment.sites.size();
                        if (fragment.sites.size() > 1) ++atac_multisite_fragments_used;
                        fragment_observations.push_back(std::move(fragment));
                        uncalibrated_fragment_observations.push_back(std::move(uncalibrated_fragment));
                    }
                }
            }

            if (status == "PASS" && ratio_sites_observed == 0) {
                status = "NO_RATIO_SITES_OBSERVED";
            } else if (status == "PASS" &&
                       ratio_molecules < static_cast<uint64_t>(opt.min_molecules)) {
                status = "LOW_RATIO_MOLECULES";
            } else if (status == "PASS" && ratio_sites_used < opt.min_sites) {
                status = "LOW_RATIO_SITES";
            }

            const AmbientQcEstimate ambient_qc = estimate_ambient_qc(
                opt, cell, pair_manifest, ambient, fixed_value);
            if (ambient_qc_requested) {
                ambient_sites_observed = ambient_qc.sites_observed;
                ambient_sites_used = ambient_qc.sites_used;
                if (ambient_qc.status == "HIGH_MT_AMBIENT") {
                    ++high_mt_ambient_cells;
                    if (status == "PASS") status = "HIGH_MT_AMBIENT";
                } else if (ambient_qc.status == "NOT_ESTIMABLE") {
                    ++ambient_qc_not_estimable_cells;
                    if (status == "PASS") status = "MT_AMBIENT_NOT_ESTIMABLE";
                }
            }

            const double pooled_rho_for_cell = pair_manifest
                ? resolve_pooled_rho(opt, rho_reference, canonical_parent1, canonical_parent2)
                : opt.pooled_rho;
            if (opt.assay_mode != "ATAC" && opt.rho_mode != "free" &&
                pair_manifest && !std::isfinite(pooled_rho_for_cell)) {
                throw std::runtime_error("No pooled rho reference matched pair " +
                    canonical_parent1 + "+" + canonical_parent2 + " for library " + opt.library_id);
            }
            const MtRatioFitConfig fit_config = make_fit_config(
                opt, pooled_rho_for_cell, ratio_molecules);
            std::string rho_mode_effective = "NOT_APPLICABLE";
            if (fit_config.use_beta_binomial) {
                if (fit_config.overdispersion_fixed) {
                    rho_mode_effective = opt.rho_mode == "low_information_fixed"
                        ? "FIXED_LOW_INFORMATION" : "FIXED";
                    ++rho_fixed_cells;
                    if (opt.rho_mode != "free") fit_mode += "_RHO_FIXED";
                } else if (fit_config.overdispersion_penalized) {
                    rho_mode_effective = "SHRINK";
                    ++rho_shrunk_cells;
                    fit_mode += "_RHO_SHRINK";
                } else {
                    rho_mode_effective = opt.rho_mode == "low_information_fixed"
                        ? "FREE_HIGH_INFORMATION" : "FREE";
                }
            }

            MtRatioFitResult fit;
            MtRatioFitResult uncalibrated_fit;
            const bool fit_attempted = status == "PASS";
            if (fit_attempted) {
                fit = opt.assay_mode == "ATAC"
                    ? fit_mt_fragment_ratio(fragment_observations, fit_config)
                    : fit_mt_ratio(observations, fit_config);
                if (!opt.site_calibration.empty()) {
                    uncalibrated_fit = opt.assay_mode == "ATAC"
                        ? fit_mt_fragment_ratio(uncalibrated_fragment_observations, fit_config)
                        : fit_mt_ratio(uncalibrated_observations, fit_config);
                }
                if (!fit.converged) {
                    status = "NOT_CONVERGED";
                } else if (fit_config.use_beta_binomial && !fit_config.overdispersion_fixed &&
                           fit.overdispersion_rho >= 0.999 * opt.overdispersion_max) {
                    status = "OVERDISPERSION_AT_BOUND";
                }
            }

            if (status == "PASS") {
                ++pass_cells;
            } else if (status == "WEAKLY_IDENTIFIABLE" ||
                       status == "OVERDISPERSION_AT_BOUND") {
                ++weak_cells;
            } else {
                ++low_cells;
            }

            const double parent1_fraction = std::isfinite(fit.parent2_fraction)
                ? 1.0 - fit.parent2_fraction : std::numeric_limits<double>::quiet_NaN();
            double se_for_ci = std::numeric_limits<double>::quiet_NaN();
            if (std::isfinite(fit.ratio_se)) se_for_ci = fit.ratio_se;
            if (std::isfinite(fit.ratio_se_robust)) {
                se_for_ci = std::isfinite(se_for_ci)
                    ? std::max(se_for_ci, fit.ratio_se_robust)
                    : fit.ratio_se_robust;
            }
            const double ci_low = std::isfinite(parent1_fraction) && std::isfinite(se_for_ci)
                ? std::max(0.0, parent1_fraction - 1.96 * se_for_ci)
                : std::numeric_limits<double>::quiet_NaN();
            const double ci_high = std::isfinite(parent1_fraction) && std::isfinite(se_for_ci)
                ? std::min(1.0, parent1_fraction + 1.96 * se_for_ci)
                : std::numeric_limits<double>::quiet_NaN();
            const uint64_t raw_total = support_parent1 + support_parent2;
            const double raw_fraction = raw_total > 0
                ? static_cast<double>(support_parent1) / raw_total
                : std::numeric_limits<double>::quiet_NaN();

            const double canonical_parent2_fraction =
                std::isfinite(fit.parent2_fraction) && pair_manifest
                    ? (assignment_matches_canonical
                        ? fit.parent2_fraction : 1.0 - fit.parent2_fraction)
                    : std::numeric_limits<double>::quiet_NaN();
            const double canonical_parent1_fraction =
                std::isfinite(canonical_parent2_fraction)
                    ? 1.0 - canonical_parent2_fraction
                    : std::numeric_limits<double>::quiet_NaN();
            const double uncalibrated_canonical_parent2_fraction =
                std::isfinite(uncalibrated_fit.parent2_fraction) && pair_manifest
                    ? (assignment_matches_canonical
                        ? uncalibrated_fit.parent2_fraction : 1.0 - uncalibrated_fit.parent2_fraction)
                    : std::numeric_limits<double>::quiet_NaN();
            const double calibration_shift_parent2 =
                std::isfinite(canonical_parent2_fraction) &&
                std::isfinite(uncalibrated_canonical_parent2_fraction)
                    ? canonical_parent2_fraction - uncalibrated_canonical_parent2_fraction
                    : std::numeric_limits<double>::quiet_NaN();

            std::string site_calibration_status = "NOT_REQUESTED";
            const uint64_t used_site_count = used_site_indices.size();
            if (!opt.site_calibration.empty()) {
                if (used_site_count == 0) site_calibration_status = "NO_USED_SITES";
                else if (fully_calibrated_used_sites == used_site_count) site_calibration_status = "FULL";
                else if (fully_calibrated_used_sites + partially_calibrated_used_sites > 0) {
                    site_calibration_status = "PARTIAL";
                } else site_calibration_status = "FALLBACK";
            }
            const double site_calibration_fraction = used_site_count > 0
                ? static_cast<double>(fully_calibrated_used_sites) /
                  static_cast<double>(used_site_count)
                : std::numeric_limits<double>::quiet_NaN();
            if (site_calibration_status == "FULL") ++calibration_full_cells;
            else if (site_calibration_status == "PARTIAL") ++calibration_partial_cells;
            else if (site_calibration_status == "FALLBACK") ++calibration_fallback_cells;
            calibration_full_site_uses += fully_calibrated_used_sites;
            calibration_partial_site_uses += partially_calibrated_used_sites;
            if (!opt.site_calibration.empty() && used_site_count >
                    fully_calibrated_used_sites + partially_calibrated_used_sites) {
                calibration_fallback_site_uses += used_site_count -
                    fully_calibrated_used_sites - partially_calibrated_used_sites;
            }

            MtRatioProfileResult profile;
            profile.profile_status = "NOT_ATTEMPTED";
            profile.inheritance_class = "AMBIGUOUS";
            profile.inheritance_class_reason = "PRIMARY_FIT_NOT_SUCCESSFUL";
            profile.single_parent_epsilon = opt.single_parent_epsilon;
            if (fit.converged && std::isfinite(fit.log_likelihood) && pair_manifest) {
                ++profile_attempted_cells;
                profile = opt.assay_mode == "ATAC"
                    ? profile_mt_fragment_ratio(fragment_observations, fit_config, fit,
                        opt.single_parent_epsilon, opt.profile_grid_step, opt.write_profile_grid)
                    : profile_mt_ratio(observations, fit_config, fit, opt.single_parent_epsilon,
                        opt.profile_grid_step, opt.write_profile_grid);
                profile = canonicalize_profile(profile, assignment_matches_canonical);
                if (profile.profile_status == "PASS") ++profile_pass_cells;
                else if (profile.profile_status == "PARTIAL") ++profile_partial_cells;
                else if (profile.profile_status == "MULTIMODAL") ++profile_multimodal_cells;
                else ++profile_failed_cells;
            }

            const std::unordered_map<int, double> site_influence =
                fit.converged && opt.write_site_counts && opt.site_influence_mode == "full"
                ? (opt.assay_mode == "ATAC"
                    ? compute_atac_site_influence(fragment_observations, fit_config, fit, used_site_indices)
                    : compute_rna_site_influence(observations, fit_config, fit))
                : std::unordered_map<int, double>();
            double max_abs_site_influence = std::numeric_limits<double>::quiet_NaN();
            double sum_abs_site_influence = 0.0;
            int most_influential_site_index = -1;
            for (const auto& item : site_influence) {
                if (!std::isfinite(item.second)) continue;
                const double value = std::fabs(item.second);
                sum_abs_site_influence += value;
                if (!std::isfinite(max_abs_site_influence) || value > max_abs_site_influence) {
                    max_abs_site_influence = value;
                    most_influential_site_index = item.first;
                }
            }
            const double fraction_support_from_top_site =
                std::isfinite(max_abs_site_influence) && sum_abs_site_influence > 0.0
                    ? max_abs_site_influence / sum_abs_site_influence
                    : (std::isfinite(max_abs_site_influence) ? 0.0
                       : std::numeric_limits<double>::quiet_NaN());
            const std::string most_influential_site = most_influential_site_index >= 0
                ? opt.mito_chrom + ":" + std::to_string(sites[most_influential_site_index].pos + 1)
                : "NA";

            if (profile.inheritance_class == "ONLY_PARENT1") ++inheritance_parent1_cells;
            else if (profile.inheritance_class == "BOTH") ++inheritance_both_cells;
            else if (profile.inheritance_class == "ONLY_PARENT2") ++inheritance_parent2_cells;
            else ++inheritance_ambiguous_cells;

            if (opt.write_profile_grid && profile_out && !profile.grid_delta_log_likelihood.empty()) {
                const double effective_grid_step = profile.grid_r.size() > 1
                    ? profile.grid_r[1] - profile.grid_r[0] : 1.0;
                const std::string profile_objective = opt.assay_mode == "ATAC"
                    ? "FRAGMENT_LOG_LIKELIHOOD"
                    : (fit_config.overdispersion_penalized
                        ? "PENALIZED_RHO_PROFILE" : "PROFILE_LOG_LIKELIHOOD");
                std::ostringstream row;
                row << (opt.legacy_panel_gt ? "NA" : opt.library_id) << '\t'
                    << cell.barcode << '\t' << cell.identity << '\t'
                    << canonical_parent1 << '\t' << canonical_parent2 << '\t'
                    << format_double(canonical_parent2_fraction) << '\t'
                    << "0.0" << '\t' << format_double(effective_grid_step, 8) << '\t'
                    << profile.grid_delta_log_likelihood.size() << '\t'
                    << profile_delta_csv(profile.grid_delta_log_likelihood) << '\t'
                    << profile.profile_status << '\t' << profile_objective << '\t'
                    << site_calibration_status << '\t' << opt.rho_mode << '\n';
                cell_result.profile_row = row.str();
            }

            const std::string likelihood_label = opt.assay_mode == "ATAC"
                ? "fragment_mixture" : opt.likelihood;
            const double reported_ambient_fraction = ambient_qc_requested
                ? ambient_qc.fraction : fit.ambient_fraction;
            ratio_out << (opt.legacy_panel_gt ? "NA" : opt.library_id) << '\t'
                      << cell.barcode << '\t' << cell.identity << '\t'
                      << cell.parent1 << '\t' << cell.parent2 << '\t'
                      << likelihood_label << '\t'
                      << format_double(parent1_fraction) << '\t'
                      << format_double(fit.parent2_fraction) << '\t'
                      << format_double(fit.ratio_se) << '\t'
                      << format_double(fit.ratio_se_robust) << '\t'
                      << format_double(ci_low) << '\t' << format_double(ci_high) << '\t'
                      << format_double(reported_ambient_fraction) << '\t'
                      << format_double(fit.ambient_se) << '\t'
                      << format_double(fit.ambient_se_robust) << '\t'
                      << format_double(fit.overdispersion_rho) << '\t'
                      << format_double(fit.overdispersion_se) << '\t'
                      << ratio_sites_available << '\t' << ambient_sites_available << '\t'
                      << ratio_sites_observed << '\t' << ambient_sites_observed << '\t'
                      << ratio_sites_used << '\t' << ambient_sites_used << '\t'
                      << ratio_molecules << '\t' << total_molecules_used << '\t'
                      << observations.size() << '\t' << support_parent1 << '\t'
                      << support_parent2 << '\t' << format_double(raw_fraction) << '\t'
                      << format_double(fit.log_likelihood) << '\t'
                      << format_double(fit.information_condition) << '\t'
                      << format_double(fit.min_information_eigenvalue) << '\t'
                      << fit.iterations << '\t' << (fit.converged ? 1 : 0) << '\t'
                      << fit_mode << '\t' << status << '\t'
                      << CELLBOUNCER_SOURCE_REVISION << '\t' << kMtModelVersion << '\t'
                      << opt.assay_mode << '\t' << assay_estimand(opt.assay_mode) << '\t'
                      << canonical_parent1 << '\t' << canonical_parent2 << '\t'
                      << (pair_manifest ? (assignment_matches_canonical ? "1" : "0") : "NA") << '\t'
                      << format_double(canonical_parent1_fraction) << '\t'
                      << format_double(canonical_parent2_fraction) << '\t'
                      << format_double(profile.parent1_ci_low) << '\t'
                      << format_double(profile.parent1_ci_high) << '\t'
                      << format_double(profile.parent2_ci_low) << '\t'
                      << format_double(profile.parent2_ci_high) << '\t'
                      << format_double(profile.profile_ci_level, 2) << '\t'
                      << profile.profile_status << '\t' << profile.evaluations << '\t'
                      << profile.failed_evaluations << '\t'
                      << format_double(opt.single_parent_epsilon, 6) << '\t'
                      << profile.inheritance_class << '\t'
                      << profile.inheritance_class_reason << '\t'
                      << format_double(profile.delta_ll_parent1_only) << '\t'
                      << format_double(profile.delta_ll_both) << '\t'
                      << format_double(profile.delta_ll_parent2_only) << '\t'
                      << panel_fingerprint << '\t' << manifest_fingerprint << '\t'
                      << calibration_ids_string(site_calibration) << '\t'
                      << calibration_fingerprint << '\t' << opt.site_calibration_stratum << '\t'
                      << format_double(site_calibration_fraction) << '\t'
                      << site_calibration_status << '\t'
                      << format_double(uncalibrated_canonical_parent2_fraction) << '\t'
                      << format_double(uncalibrated_fit.log_likelihood) << '\t'
                      << format_double(calibration_shift_parent2) << '\t'
                      << format_double(max_abs_site_influence) << '\t'
                      << most_influential_site << '\t'
                      << format_double(fraction_support_from_top_site) << '\t'
                      << opt.rho_mode << '\t' << rho_mode_effective << '\t'
                      << format_double((fit_config.overdispersion_fixed || fit_config.overdispersion_penalized)
                           ? pooled_rho_for_cell : std::numeric_limits<double>::quiet_NaN()) << '\t'
                      << format_double(fit_config.overdispersion_penalized
                           ? opt.rho_prior_strength : 0.0) << '\t'
                      << format_double(fit.objective_value) << '\t'
                      << (opt.calibration_reference.empty() ? 0 : 1) << '\t'
                      << (cell.calibration_true_parent.empty() ? "NA" : cell.calibration_true_parent) << '\t'
                      << (cell.original_identity.empty() ? "NA" : cell.original_identity) << '\t'
                      << (opt.assay_mode == "ATAC" ? fragment_observations.size() : 0) << '\t'
                      << (opt.assay_mode == "ATAC" ? atac_multisite_fragments_used : 0) << '\t'
                      << (opt.assay_mode == "ATAC" ? atac_fragment_site_observations_used : 0) << '\t'
                      << pair_definition_source << '\t'
                      << experimental_mode << '\t' << ambient_mode << '\t'
                      << ambient_profile_id << '\t' << ambient_profile_fingerprint << '\t'
                      << ambient_fraction_fingerprint << '\t'
                      << format_double(opt.ambient_qc_max) << '\t'
                      << ambient_qc.status << '\t' << ambient_qc.reason << '\t'
                      << ambient_qc.molecules << '\n';

            if (opt.write_site_counts && pair_manifest) {
                for (const auto& item : cell.counts) {
                    const int site_index = item.first;
                    const auto manifest_found = pair_manifest->sites.find(site_index);
                    if (manifest_found == pair_manifest->sites.end()) continue;
                    const ManifestSite& manifest_site = manifest_found->second;
                    int p1_state = manifest_site.canonical_parent1_state;
                    int p2_state = manifest_site.canonical_parent2_state;
                    if (!assignment_matches_canonical) std::swap(p1_state, p2_state);
                    const bool ratio_informative = p1_state != p2_state;
                    const bool ambient_anchor = manifest_site.site_class == MANIFEST_AMBIENT_ONLY;
                    const AlleleCount& count = item.second;
                    const MtSite& site = sites[site_index];
                    const CalibrationProbabilityResult probabilities = calibration_probabilities(
                        opt, site_calibration, *pair_manifest, site_index,
                        assignment_matches_canonical, p1_state, p2_state);
                    MtRatioObservation site_observation;
                    site_observation.site_index = site_index;
                    site_observation.ref = count.ref;
                    site_observation.alt = count.alt;
                    site_observation.parent1_alt_probability = probabilities.parent1_alt_probability;
                    site_observation.parent2_alt_probability = probabilities.parent2_alt_probability;
                    site_observation.ambient_alt_probability = 0.5;
                    site_observation.ratio_informative = ratio_informative;
                    site_observation.ambient_anchor = ambient_anchor;
                    const bool used = used_site_indices.count(site_index) > 0;
                    const double predicted = used && fit.converged
                        ? mt_ratio_predicted_alt_probability(site_observation,
                            fit.parent2_fraction, fit.ambient_fraction)
                        : std::numeric_limits<double>::quiet_NaN();
                    const bool site_beta = opt.assay_mode != "ATAC" && fit_config.use_beta_binomial;
                    const double site_rho = site_beta && std::isfinite(fit.overdispersion_rho)
                        ? fit.overdispersion_rho : 0.0;
                    const double site_ll = used && std::isfinite(predicted)
                        ? mt_ratio_count_log_likelihood(count.ref, count.alt, predicted, site_rho, site_beta)
                        : std::numeric_limits<double>::quiet_NaN();
                    const double deviance = used && std::isfinite(predicted)
                        ? site_deviance_residual(site_observation, predicted, site_rho, site_beta)
                        : std::numeric_limits<double>::quiet_NaN();
                    const auto influence_it = site_influence.find(site_index);
                    const double influence = influence_it == site_influence.end()
                        ? std::numeric_limits<double>::quiet_NaN() : influence_it->second;
                    const std::string influence_status = !used ? "NOT_USED" :
                        (opt.site_influence_mode == "none" ? "DISABLED" :
                         (!fit.converged ? "PRIMARY_FIT_UNAVAILABLE" :
                          (std::isfinite(influence) ? "PASS" : "LOO_FAILED")));
                    site_out << (opt.legacy_panel_gt ? "NA" : opt.library_id) << '\t'
                             << cell.barcode << '\t' << cell.identity << '\t'
                             << opt.mito_chrom << '\t' << site.pos + 1 << '\t'
                             << site.ref << '\t' << site.alt << '\t'
                             << manifest_site_class_name(manifest_site.site_class) << '\t'
                             << p1_state << '\t' << p2_state << '\t'
                             << (ratio_informative ? 1 : 0) << '\t'
                             << (ambient_anchor ? 1 : 0) << '\t'
                             << count.ref << '\t' << count.alt << '\t'
                             << ambient[site_index].ref << '\t' << ambient[site_index].alt << '\t'
                             << format_double(ambient[site_index].alt_fraction) << '\t'
                             << (used ? 1 : 0) << '\t'
                             << format_double(probabilities.parent1_alt_probability) << '\t'
                             << format_double(probabilities.parent2_alt_probability) << '\t'
                             << format_double(predicted) << '\t' << format_double(site_ll) << '\t'
                             << format_double(deviance) << '\t' << format_double(influence) << '\t'
                             << influence_status << '\t' << probabilities.status << '\t'
                             << (cell.calibration_true_parent.empty() ? "NA" : cell.calibration_true_parent) << '\t'
                             << opt.assay_mode << '\t'
                             << canonical_parent1 << '\t' << canonical_parent2
                             << '\n';
                }
            }
                cell_result.ratio_row = ratio_out.str();
                cell_result.site_rows = site_out.str();
            } catch (...) {
                cell_result.error = std::current_exception();
            }
        };

        const size_t worker_count = cells.empty() ? 0 : std::min(
            cells.size(), static_cast<size_t>(std::max(1, opt.threads)));
        if (worker_count <= 1) {
            for (size_t cell_index = 0; cell_index < cells.size(); ++cell_index) {
                process_cell(cell_index);
            }
        } else {
            std::atomic<size_t> next_cell{0};
            std::vector<std::thread> workers;
            workers.reserve(worker_count);
            for (size_t worker = 0; worker < worker_count; ++worker) {
                workers.emplace_back([&]() {
                    while (true) {
                        const size_t cell_index = next_cell.fetch_add(1);
                        if (cell_index >= cells.size()) break;
                        process_cell(cell_index);
                    }
                });
            }
            for (std::thread& worker : workers) worker.join();
        }

        // Preserve deterministic failure attribution and output ordering.
        for (const CellResult& cell_result : cell_results) {
            if (cell_result.error) std::rethrow_exception(cell_result.error);
        }
        for (const CellResult& cell_result : cell_results) {
            pass_cells += cell_result.counters.pass_cells;
            weak_cells += cell_result.counters.weak_cells;
            low_cells += cell_result.counters.low_cells;
            fixed_cells += cell_result.counters.fixed_cells;
            pair_missing_cells += cell_result.counters.pair_missing_cells;
            derived_pair_cells += cell_result.counters.derived_pair_cells;
            collapsed_pair_cells += cell_result.counters.collapsed_pair_cells;
            unanchored_cells += cell_result.counters.unanchored_cells;
            profile_attempted_cells += cell_result.counters.profile_attempted_cells;
            profile_pass_cells += cell_result.counters.profile_pass_cells;
            profile_partial_cells += cell_result.counters.profile_partial_cells;
            profile_failed_cells += cell_result.counters.profile_failed_cells;
            profile_multimodal_cells += cell_result.counters.profile_multimodal_cells;
            inheritance_parent1_cells += cell_result.counters.inheritance_parent1_cells;
            inheritance_both_cells += cell_result.counters.inheritance_both_cells;
            inheritance_parent2_cells += cell_result.counters.inheritance_parent2_cells;
            inheritance_ambiguous_cells += cell_result.counters.inheritance_ambiguous_cells;
            calibration_full_cells += cell_result.counters.calibration_full_cells;
            calibration_partial_cells += cell_result.counters.calibration_partial_cells;
            calibration_fallback_cells += cell_result.counters.calibration_fallback_cells;
            calibration_full_site_uses += cell_result.counters.calibration_full_site_uses;
            calibration_partial_site_uses += cell_result.counters.calibration_partial_site_uses;
            calibration_fallback_site_uses += cell_result.counters.calibration_fallback_site_uses;
            rho_fixed_cells += cell_result.counters.rho_fixed_cells;
            rho_shrunk_cells += cell_result.counters.rho_shrunk_cells;
            high_mt_ambient_cells += cell_result.counters.high_mt_ambient_cells;
            ambient_qc_not_estimable_cells += cell_result.counters.ambient_qc_not_estimable_cells;
            ratio_out << cell_result.ratio_row;
            if (profile_out && !cell_result.profile_row.empty()) {
                gz_write_text(profile_out, cell_result.profile_row);
            }
            if (opt.write_site_counts && !cell_result.site_rows.empty()) {
                site_out << cell_result.site_rows;
            }
        }
        if (profile_out) {
            if (gzclose(profile_out) != Z_OK) {
                throw std::runtime_error("Failed to close compressed mt profile output: " + profile_path);
            }
            profile_out = nullptr;
        }

        const std::string qc_path = opt.output_prefix + ".mt_qc.tsv";
        std::ofstream qc(qc_path);
        if (!qc) throw std::runtime_error("Could not write: " + qc_path);
        qc << "metric\tvalue\n";
        qc << "library_id\t" << (opt.legacy_panel_gt ? "NA" : opt.library_id) << '\n';
        qc << "panel_mode\t" << (opt.legacy_panel_gt ? "LEGACY_GT_DIAGNOSTIC" :
                                                "LIBRARY_SPECIFIC_THREE_CLASS_MANIFEST") << '\n';
        qc << "site_manifest\t" << (opt.legacy_panel_gt ? "NA" : opt.site_manifest) << '\n';
        qc << "likelihood\t" << opt.likelihood << '\n';
        qc << "overdispersion_initial\t" << opt.overdispersion_initial << '\n';
        qc << "overdispersion_max\t" << opt.overdispersion_max << '\n';
        qc << "min_ambient_only_sites\t" << opt.min_ambient_only_sites << '\n';
        qc << "allow_unanchored_ambient\t" << (opt.allow_unanchored_ambient ? 1 : 0) << '\n';
        qc << "panel_sites\t" << sites.size() << '\n';
        qc << "manifest_pairs\t" << pair_manifests.size() << '\n';
        qc << "derived_reconciled_pairs\t" << derived_pair_manifests.size() << '\n';
        qc << "derived_pair_panel_min_depth\t" << kReconciledPairPanelMinDepth << '\n';
        qc << "derived_pair_panel_homoplasmy_af\t" << kReconciledPairHomoplasmyAf << '\n';
        qc << "fusion_cells\t" << cells.size() << '\n';
        qc << "empty_barcodes\t" << empty_barcodes.size() << '\n';
        qc << "usable_ambient_sites\t" << usable_ambient_sites << '\n';
        qc << "fixed_ambient_cells\t" << fixed_cells << '\n';
        qc << "ambient_qc_enabled\t" << (ambient_qc_requested ? 1 : 0) << '\n';
        qc << "ambient_qc_max\t" << format_double(opt.ambient_qc_max) << '\n';
        qc << "high_mt_ambient_cells\t" << high_mt_ambient_cells << '\n';
        qc << "ambient_qc_not_estimable_cells\t" << ambient_qc_not_estimable_cells << '\n';
        qc << "ambient_background_used_in_ratio_likelihood\t0\n";
        qc << "pass_cells\t" << pass_cells << '\n';
        qc << "weak_cells\t" << weak_cells << '\n';
        qc << "low_or_unfit_cells\t" << low_cells << '\n';
        qc << "pair_missing_cells\t" << pair_missing_cells << '\n';
        qc << "derived_reconciled_pair_cells\t" << derived_pair_cells << '\n';
        qc << "haplotype_collapsed_pair_cells\t" << collapsed_pair_cells << '\n';
        qc << "unanchored_joint_fit_cells\t" << unanchored_cells << '\n';
        qc << "reads_seen_chrM\t" << stats.seen << '\n';
        qc << "reads_accepted_for_pileup\t" << stats.accepted_for_pileup << '\n';
        qc << "reject_unmapped\t" << stats.reject_unmapped << '\n';
        qc << "reject_secondary\t" << stats.reject_secondary << '\n';
        qc << "reject_supplementary\t" << stats.reject_supplementary << '\n';
        qc << "reject_qcfail\t" << stats.reject_qcfail << '\n';
        qc << "reject_duplicate\t" << stats.reject_duplicate << '\n';
        qc << "reject_mapq\t" << stats.reject_mapq << '\n';
        qc << "reject_multimapping\t" << stats.reject_multimapping << '\n';
        qc << "pileup_bases_at_panel_sites\t" << stats.pileup_bases << '\n';
        qc << "reject_baseq\t" << stats.reject_baseq << '\n';
        qc << "reject_missing_barcode\t" << stats.reject_missing_barcode << '\n';
        qc << "reject_irrelevant_barcode\t" << stats.reject_irrelevant_barcode << '\n';
        qc << "reject_missing_umi\t" << stats.reject_missing_umi << '\n';
        qc << "reject_nonpanel_allele\t" << stats.reject_nonpanel_allele << '\n';
        qc << "conflicting_molecules\t" << stats.conflicting_molecules << '\n';
        qc << "accepted_molecule_site_observations\t" << stats.accepted_observations << '\n';
        qc << "experimental_mode\t" << experimental_mode << '\n';
        qc << "ambient_mode\t" << ambient_mode << '\n';
        qc << "ambient_profile_id\t" << ambient_profile_id << '\n';
        qc << "ambient_profile_fingerprint\t" << ambient_profile_fingerprint << '\n';
        qc << "ambient_fraction_file\t"
           << (opt.ambient_fraction_file.empty() ? "NA" : opt.ambient_fraction_file) << '\n';
        qc << "ambient_fraction_file_fingerprint\t" << ambient_fraction_fingerprint << '\n';
        qc << "min_mapq\t" << opt.min_mapq << '\n';
        qc << "min_baseq\t" << opt.min_baseq << '\n';
        qc << "error_rate\t" << opt.error_rate << '\n';
        qc << "ambiguous_alignments_allowed\t" << (opt.allow_ambiguous_alignments ? 1 : 0) << '\n';
        qc << "duplicates_kept\t" << (opt.keep_duplicates ? 1 : 0) << '\n';
        qc << "source_revision\t" << CELLBOUNCER_SOURCE_REVISION << '\n';
        qc << "model_version\t" << kMtModelVersion << '\n';
        qc << "assay_mode\t" << opt.assay_mode << '\n';
        qc << "profile_summary_attempted_cells\t" << profile_attempted_cells << '\n';
        qc << "profile_summary_pass_cells\t" << profile_pass_cells << '\n';
        qc << "profile_summary_partial_cells\t" << profile_partial_cells << '\n';
        qc << "profile_summary_failed_cells\t" << profile_failed_cells << '\n';
        qc << "profile_multimodal_cells\t" << profile_multimodal_cells << '\n';
        qc << "profile_grid_written\t" << (opt.write_profile_grid ? 1 : 0) << '\n';
        qc << "profile_grid_step\t" << opt.profile_grid_step << '\n';
        qc << "single_parent_epsilon\t" << opt.single_parent_epsilon << '\n';
        qc << "inheritance_only_parent1_cells\t" << inheritance_parent1_cells << '\n';
        qc << "inheritance_both_cells\t" << inheritance_both_cells << '\n';
        qc << "inheritance_only_parent2_cells\t" << inheritance_parent2_cells << '\n';
        qc << "inheritance_ambiguous_cells\t" << inheritance_ambiguous_cells << '\n';
        qc << "estimand\t" << assay_estimand(opt.assay_mode) << '\n';
        qc << "panel_fingerprint\t" << panel_fingerprint << '\n';
        qc << "manifest_fingerprint\t" << manifest_fingerprint << '\n';
        qc << "site_calibration\t" << (opt.site_calibration.empty() ? "NA" : opt.site_calibration) << '\n';
        qc << "site_calibration_id\t" << calibration_ids_string(site_calibration) << '\n';
        qc << "site_calibration_fingerprint\t" << calibration_fingerprint << '\n';
        qc << "site_calibration_stratum\t" << (opt.site_calibration_stratum.empty() ? "NA" : opt.site_calibration_stratum) << '\n';
        qc << "site_calibration_full_cells\t" << calibration_full_cells << '\n';
        qc << "site_calibration_partial_cells\t" << calibration_partial_cells << '\n';
        qc << "site_calibration_fallback_cells\t" << calibration_fallback_cells << '\n';
        qc << "site_calibration_full_site_uses\t" << calibration_full_site_uses << '\n';
        qc << "site_calibration_partial_site_uses\t" << calibration_partial_site_uses << '\n';
        qc << "site_calibration_fallback_site_uses\t" << calibration_fallback_site_uses << '\n';
        qc << "calibration_reference_mode\t" << (opt.calibration_reference.empty() ? 0 : 1) << '\n';
        qc << "calibration_reference\t" << (opt.calibration_reference.empty() ? "NA" : opt.calibration_reference) << '\n';
        qc << "rho_mode\t" << opt.rho_mode << '\n';
        qc << "pooled_rho_option\t" << format_double(opt.pooled_rho) << '\n';
        qc << "rho_reference\t" << (opt.rho_reference.empty() ? "NA" : opt.rho_reference) << '\n';
        qc << "rho_low_information_molecules\t" << opt.rho_low_information_molecules << '\n';
        qc << "rho_prior_strength\t" << opt.rho_prior_strength << '\n';
        qc << "rho_fixed_cells\t" << rho_fixed_cells << '\n';
        qc << "rho_shrunk_cells\t" << rho_shrunk_cells << '\n';
        qc << "atac_formal_fragment_mode\t" << (opt.assay_mode == "ATAC" ? 1 : 0) << '\n';
        qc << "atac_include_singletons\t" << (opt.atac_include_singletons ? 1 : 0) << '\n';
        qc << "atac_reads_considered\t" << stats.atac_reads_considered << '\n';
        qc << "atac_reject_unpaired\t" << stats.atac_reject_unpaired << '\n';
        qc << "atac_reject_mate_off_mito\t" << stats.atac_reject_mate_off_mito << '\n';
        qc << "atac_orphan_fragments\t" << stats.atac_orphan_fragments << '\n';
        qc << "atac_fragments_accepted\t" << stats.atac_fragments_accepted << '\n';
        qc << "atac_fragments_multisite\t" << stats.atac_fragments_multisite << '\n';
        qc << "atac_fragment_site_observations\t" << stats.atac_fragment_site_observations << '\n';
        qc << "atac_overlap_agreements\t" << stats.atac_overlap_agreements << '\n';
        qc << "atac_overlap_conflicts\t" << stats.atac_overlap_conflicts << '\n';

        std::fprintf(stderr, "\nMitochondrial fusion-ratio estimation complete\n");
        std::fprintf(stderr, "  cells: pass=%llu weak=%llu other=%llu\n",
            static_cast<unsigned long long>(pass_cells),
            static_cast<unsigned long long>(weak_cells),
            static_cast<unsigned long long>(low_cells));
        if (opt.assay_mode == "ATAC") {
            std::fprintf(stderr, "  accepted ATAC fragments: %llu (%llu fragment-site observations)\n",
                static_cast<unsigned long long>(stats.atac_fragments_accepted),
                static_cast<unsigned long long>(stats.atac_fragment_site_observations));
        } else {
            std::fprintf(stderr, "  accepted molecule-site observations: %llu\n",
                static_cast<unsigned long long>(stats.accepted_observations));
        }
        std::fprintf(stderr, "  ratio table: %s\n", ratio_path.c_str());
        std::fprintf(stderr, "  ambient profile: %s\n", ambient_path.c_str());
        if (opt.write_profile_grid) {
            std::fprintf(stderr, "  profile grid: %s\n", profile_path.c_str());
        }
        std::fprintf(stderr, "  QC: %s\n", qc_path.c_str());
        return 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "ERROR: %s\n", error.what());
        return 1;
    }
}
