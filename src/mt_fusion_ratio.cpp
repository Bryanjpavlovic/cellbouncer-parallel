#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <zlib.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cctype>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
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
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr double kTiny = 1e-12;

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
    std::string likelihood = "beta_binomial";
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
    bool ambient_none = false;
    bool allow_no_umi = false;
    bool keep_duplicates = false;
    bool allow_ambiguous_alignments = false;
    bool allow_unanchored_ambient = false;
    bool legacy_panel_gt = false;
    bool write_site_counts = true;
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
    std::unordered_map<int, AlleleCount> counts;
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

struct SiteObservation {
    int site_index = -1;
    uint64_t ref = 0;
    uint64_t alt = 0;
    double parent1_alt_probability = 0.0;
    double parent2_alt_probability = 0.0;
    double ambient_alt_probability = 0.0;
    ManifestSiteClass site_class = MANIFEST_LEGACY_RATIO;
    bool ratio_informative = false;
    bool ambient_anchor = false;
};

struct FitResult {
    double parent1_weight = std::numeric_limits<double>::quiet_NaN();
    double parent2_weight = std::numeric_limits<double>::quiet_NaN();
    double ambient_weight = std::numeric_limits<double>::quiet_NaN();
    double parent2_fraction = std::numeric_limits<double>::quiet_NaN();
    double log_likelihood = std::numeric_limits<double>::quiet_NaN();
    double ratio_se = std::numeric_limits<double>::quiet_NaN();
    double ratio_se_robust = std::numeric_limits<double>::quiet_NaN();
    double ambient_se = std::numeric_limits<double>::quiet_NaN();
    double ambient_se_robust = std::numeric_limits<double>::quiet_NaN();
    double overdispersion_rho = std::numeric_limits<double>::quiet_NaN();
    double overdispersion_se = std::numeric_limits<double>::quiet_NaN();
    double information_condition = std::numeric_limits<double>::quiet_NaN();
    double min_information_eigenvalue = std::numeric_limits<double>::quiet_NaN();
    int iterations = 0;
    bool converged = false;
    bool joint_ambient = false;
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
        "Ambient profile (choose exactly one):\n"
        "  --empty_barcodes FILE         Same-library empty-droplet barcode list\n"
        "  --ambient_profile FILE        Site-level mt ambient-profile TSV\n"
        "  --ambient_none                Explicitly disable ambient correction\n"
        "\n"
        "Likelihood and identifiability:\n"
        "  --likelihood NAME             beta_binomial (default) or binomial\n"
        "  --overdispersion_initial X    Initial beta-binomial rho [0.02]\n"
        "  --overdispersion_max X        Upper bound for fitted rho [0.25]\n"
        "  --min_ambient_only_sites INT  Observed ambient anchors required for joint c [1]\n"
        "  --allow_unanchored_ambient    Diagnostic override; retains weak-ID labeling\n"
        "  --ambient_fraction_file FILE  Explicit fixed per-cell c_mito override\n"
        "  --legacy_panel_gt             Diagnostic only; ignore manifest and use GT-different sites\n"
        "\n"
        "Counting and filters:\n"
        "  --mito_chrom NAME             Mitochondrial contig [chrM]\n"
        "  --barcode_tag TAG             Cell barcode BAM tag [CB]\n"
        "  --umi_tag TAG                 Corrected UMI BAM tag [UB]\n"
        "  --allow_no_umi                Use read name when UB is absent\n"
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
        "  --threads INT                 HTSlib decompression threads [4]\n"
        "  --no_site_counts              Suppress per-cell site-count output\n"
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
    if (opt.barcode_tag.size() != 2 || opt.umi_tag.size() != 2 ||
        opt.min_mapq < 0 || opt.min_baseq < 0 || opt.min_molecules < 0 ||
        opt.min_sites < 1 || opt.ambient_min_molecules < 1 ||
        opt.min_ambient_only_sites < 0 ||
        opt.max_iterations < 1 || opt.pileup_max_depth < 1 || opt.threads < 1 ||
        opt.error_rate <= 0.0 || opt.error_rate >= 0.5 || opt.tolerance <= 0.0 ||
        opt.overdispersion_initial < 0.0 ||
        opt.overdispersion_max <= 0.0 || opt.overdispersion_max >= 1.0 ||
        opt.overdispersion_initial > opt.overdispersion_max ||
        (opt.likelihood != "beta_binomial" && opt.likelihood != "binomial")) {
        std::fprintf(stderr, "ERROR: invalid numeric/tag option\n");
        std::exit(2);
    }
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
        if (fields[library_col] != opt.library_id) continue;
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

std::vector<CellInfo> load_assignments(const std::string& filename,
                                       const std::unordered_map<std::string, int>& sample_index,
                                       std::unordered_map<std::string, int>& barcode_to_cell) {
    std::vector<CellInfo> cells;
    for (const std::string& line : read_text_lines(filename)) {
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> fields = split_ws(line);
        if (fields.size() < 2 || fields[0] == "barcode") continue;
        const std::string& barcode = fields[0];
        const std::string& identity = fields[1];
        const size_t plus = identity.find('+');
        if (plus == std::string::npos || identity.find('+', plus + 1) != std::string::npos) continue;
        const std::string parent1 = identity.substr(0, plus);
        const std::string parent2 = identity.substr(plus + 1);
        const auto p1 = sample_index.find(parent1);
        const auto p2 = sample_index.find(parent2);
        if (p1 == sample_index.end() || p2 == sample_index.end()) {
            throw std::runtime_error("Assignment parent missing from mt-panel VCF: " + identity);
        }
        if (barcode_to_cell.count(barcode)) {
            throw std::runtime_error("Duplicate barcode in assignments: " + barcode);
        }
        CellInfo cell;
        cell.barcode = barcode;
        cell.identity = identity;
        cell.parent1 = parent1;
        cell.parent2 = parent2;
        cell.parent1_index = p1->second;
        cell.parent2_index = p2->second;
        barcode_to_cell[barcode] = static_cast<int>(cells.size());
        cells.push_back(std::move(cell));
    }
    if (cells.empty()) throw std::runtime_error(
        "No two-parent fusion assignments were found in: " + filename);
    return cells;
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

double clamp_probability(double value) {
    return std::max(kTiny, std::min(1.0 - kTiny, value));
}

double predicted_alt_probability(const SiteObservation& obs,
                                 double parent2_fraction,
                                 double ambient_fraction) {
    const double r = std::max(0.0, std::min(1.0, parent2_fraction));
    const double c = std::max(0.0, std::min(0.99, ambient_fraction));
    const double parental = (1.0 - r) * obs.parent1_alt_probability +
                            r * obs.parent2_alt_probability;
    return clamp_probability((1.0 - c) * parental + c * obs.ambient_alt_probability);
}

double binomial_log_pmf(uint64_t alt, uint64_t ref, double q) {
    const double n = static_cast<double>(alt + ref);
    const double k = static_cast<double>(alt);
    q = clamp_probability(q);
    return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) -
           std::lgamma(n - k + 1.0) + k * std::log(q) +
           (n - k) * std::log1p(-q);
}

double beta_binomial_log_pmf(uint64_t alt,
                             uint64_t ref,
                             double q,
                             double rho) {
    if (rho <= 1e-9 || alt + ref <= 1) return binomial_log_pmf(alt, ref, q);
    q = clamp_probability(q);
    rho = std::max(1e-9, std::min(0.999999, rho));
    const double concentration = 1.0 / rho - 1.0;
    const double alpha = std::max(kTiny, q * concentration);
    const double beta = std::max(kTiny, (1.0 - q) * concentration);
    const double n = static_cast<double>(alt + ref);
    const double k = static_cast<double>(alt);
    return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) -
           std::lgamma(n - k + 1.0) +
           std::lgamma(k + alpha) + std::lgamma(n - k + beta) -
           std::lgamma(n + alpha + beta) -
           std::lgamma(alpha) - std::lgamma(beta) +
           std::lgamma(alpha + beta);
}

double site_log_likelihood(const SiteObservation& obs,
                           double parent2_fraction,
                           double ambient_fraction,
                           double rho,
                           bool use_beta_binomial) {
    const double q = predicted_alt_probability(obs, parent2_fraction, ambient_fraction);
    return use_beta_binomial
        ? beta_binomial_log_pmf(obs.alt, obs.ref, q, rho)
        : binomial_log_pmf(obs.alt, obs.ref, q);
}

double log_likelihood(const std::vector<SiteObservation>& observations,
                      double parent2_fraction,
                      double ambient_fraction,
                      double rho,
                      bool use_beta_binomial) {
    double value = 0.0;
    for (const SiteObservation& obs : observations) {
        value += site_log_likelihood(obs, parent2_fraction, ambient_fraction,
                                     rho, use_beta_binomial);
    }
    return value;
}

template <typename Function>
double golden_maximize(Function function,
                       double lower,
                       double upper,
                       int max_iterations,
                       double tolerance,
                       int& iterations) {
    if (!(upper > lower)) {
        iterations = 0;
        return lower;
    }
    const double phi = (std::sqrt(5.0) - 1.0) / 2.0;
    double left = lower;
    double right = upper;
    double x1 = right - phi * (right - left);
    double x2 = left + phi * (right - left);
    double f1 = function(x1);
    double f2 = function(x2);
    iterations = 0;
    while (iterations < max_iterations && right - left > tolerance) {
        if (f1 < f2) {
            left = x1;
            x1 = x2;
            f1 = f2;
            x2 = left + phi * (right - left);
            f2 = function(x2);
        } else {
            right = x2;
            x2 = x1;
            f2 = f1;
            x1 = right - phi * (right - left);
            f1 = function(x1);
        }
        ++iterations;
    }
    const double interior = 0.5 * (left + right);
    const double f_interior = function(interior);
    const double f_lower = function(lower);
    const double f_upper = function(upper);
    if (f_lower >= f_interior && f_lower >= f_upper) return lower;
    if (f_upper >= f_interior && f_upper >= f_lower) return upper;
    return interior;
}

struct FitStart {
    double ratio;
    double ambient;
    double rho;
};

FitResult fit_parameters(const std::vector<SiteObservation>& observations,
                         const Options& opt,
                         bool ambient_enabled,
                         const double* fixed_ambient) {
    const bool use_beta = opt.likelihood == "beta_binomial";
    const bool ambient_free = ambient_enabled && fixed_ambient == nullptr;
    const double fixed_c = !ambient_enabled ? 0.0 :
                           (fixed_ambient != nullptr ? *fixed_ambient : 0.0);
    const double rho_lower = 1e-8;
    std::vector<FitStart> starts;
    starts.push_back({0.50, ambient_free ? 0.05 : fixed_c, opt.overdispersion_initial});
    starts.push_back({0.20, ambient_free ? 0.10 : fixed_c, opt.overdispersion_initial});
    starts.push_back({0.80, ambient_free ? 0.10 : fixed_c, opt.overdispersion_initial});
    starts.push_back({0.50, ambient_free ? 0.30 : fixed_c, std::min(0.10, opt.overdispersion_max)});
    starts.push_back({0.10, ambient_free ? 0.50 : fixed_c, std::min(0.05, opt.overdispersion_max)});
    starts.push_back({0.90, ambient_free ? 0.50 : fixed_c, std::min(0.05, opt.overdispersion_max)});

    FitResult best;
    best.joint_ambient = ambient_free;
    double best_ll = -std::numeric_limits<double>::infinity();

    for (const FitStart& start : starts) {
        double r = std::max(0.0, std::min(1.0, start.ratio));
        double c = ambient_free ? std::max(0.0, std::min(0.99, start.ambient)) : fixed_c;
        double rho = use_beta
            ? std::max(rho_lower, std::min(opt.overdispersion_max, start.rho))
            : 0.0;
        bool converged = false;
        int iteration = 0;
        double previous_ll = log_likelihood(observations, r, c, rho, use_beta);

        for (; iteration < opt.max_iterations; ++iteration) {
            const double old_r = r;
            const double old_c = c;
            const double old_rho = rho;
            int inner = 0;
            r = golden_maximize(
                [&](double candidate) {
                    return log_likelihood(observations, candidate, c, rho, use_beta);
                }, 0.0, 1.0, 100, opt.tolerance, inner);

            if (ambient_free) {
                c = golden_maximize(
                    [&](double candidate) {
                        return log_likelihood(observations, r, candidate, rho, use_beta);
                    }, 0.0, 0.99, 100, opt.tolerance, inner);
            }

            if (use_beta) {
                rho = golden_maximize(
                    [&](double candidate) {
                        return log_likelihood(observations, r, c, candidate, true);
                    }, rho_lower, opt.overdispersion_max, 100, opt.tolerance, inner);
                const double boundary_snap = std::max(10.0 * opt.tolerance, 1e-7);
                if (rho <= rho_lower + boundary_snap) rho = rho_lower;
                if (rho >= opt.overdispersion_max - boundary_snap) {
                    rho = opt.overdispersion_max;
                }
            }

            const double current_ll = log_likelihood(observations, r, c, rho, use_beta);
            const double parameter_delta = std::max({std::fabs(r - old_r),
                                                     std::fabs(c - old_c),
                                                     std::fabs(rho - old_rho)});
            const double relative_ll_delta =
                std::fabs(current_ll - previous_ll) /
                (1.0 + std::fabs(previous_ll));
            previous_ll = current_ll;
            const double parameter_tolerance = std::max(10.0 * opt.tolerance, 1e-7);
            const double likelihood_tolerance = std::max(opt.tolerance, 1e-10);
            if (parameter_delta <= parameter_tolerance &&
                relative_ll_delta <= likelihood_tolerance) {
                converged = true;
                ++iteration;
                break;
            }
        }

        const double ll = log_likelihood(observations, r, c, rho, use_beta);
        if (ll > best_ll) {
            best_ll = ll;
            best.parent2_fraction = r;
            best.ambient_weight = c;
            best.parent1_weight = (1.0 - c) * (1.0 - r);
            best.parent2_weight = (1.0 - c) * r;
            best.overdispersion_rho = use_beta ? rho : 0.0;
            best.log_likelihood = ll;
            // Report coordinate-ascent sweeps. Golden-section evaluations are
            // an implementation detail and should not inflate this diagnostic.
            best.iterations = iteration;
            best.converged = converged;
            best.joint_ambient = ambient_free;
        }
    }
    return best;
}

double finite_step(double value, double lower, double upper) {
    const double distance = std::min(value - lower, upper - value);
    if (distance <= 1e-7) return 0.0;
    return std::min(1e-4, 0.25 * distance);
}

void compute_uncertainty(const std::vector<SiteObservation>& observations,
                         const Options& opt,
                         FitResult& fit,
                         bool ambient_parameter_free) {
    const bool use_beta = opt.likelihood == "beta_binomial";
    const double r = fit.parent2_fraction;
    const double c = fit.ambient_weight;
    const double rho = use_beta ? fit.overdispersion_rho : 0.0;
    const double hr = finite_step(r, 0.0, 1.0);
    const double hc = ambient_parameter_free ? finite_step(c, 0.0, 0.99) : 0.0;
    if (hr <= 0.0) return;

    const auto total_ll = [&](double rr, double cc) {
        return log_likelihood(observations, rr, cc, rho, use_beta);
    };
    const double center = total_ll(r, c);
    const double i_rr = -(total_ll(r + hr, c) - 2.0 * center +
                          total_ll(r - hr, c)) / (hr * hr);

    if (!ambient_parameter_free) {
        if (i_rr > kTiny && std::isfinite(i_rr)) {
            fit.ratio_se = std::sqrt(1.0 / i_rr);
            double b_rr = 0.0;
            for (const SiteObservation& obs : observations) {
                const double score = (site_log_likelihood(obs, r + hr, c, rho, use_beta) -
                                      site_log_likelihood(obs, r - hr, c, rho, use_beta)) /
                                     (2.0 * hr);
                b_rr += score * score;
            }
            fit.ratio_se_robust = std::sqrt(std::max(0.0, b_rr / (i_rr * i_rr)));
            fit.information_condition = 1.0;
            fit.min_information_eigenvalue = i_rr;
        }
    } else if (hc > 0.0) {
        const double i_cc = -(total_ll(r, c + hc) - 2.0 * center +
                              total_ll(r, c - hc)) / (hc * hc);
        const double i_rc = -(total_ll(r + hr, c + hc) - total_ll(r + hr, c - hc) -
                              total_ll(r - hr, c + hc) + total_ll(r - hr, c - hc)) /
                            (4.0 * hr * hc);
        const double trace = i_rr + i_cc;
        const double discriminant = std::sqrt(std::max(0.0,
            (i_rr - i_cc) * (i_rr - i_cc) + 4.0 * i_rc * i_rc));
        const double lambda_max = 0.5 * (trace + discriminant);
        const double lambda_min = 0.5 * (trace - discriminant);
        fit.min_information_eigenvalue = lambda_min;
        fit.information_condition = lambda_max / std::max(lambda_min, kTiny);
        const double determinant = i_rr * i_cc - i_rc * i_rc;
        if (determinant > kTiny && i_rr > 0.0 && i_cc > 0.0 &&
            std::isfinite(determinant)) {
            const double inv_rr = i_cc / determinant;
            const double inv_rc = -i_rc / determinant;
            const double inv_cc = i_rr / determinant;
            fit.ratio_se = std::sqrt(std::max(0.0, inv_rr));
            fit.ambient_se = std::sqrt(std::max(0.0, inv_cc));

            double b_rr = 0.0;
            double b_rc = 0.0;
            double b_cc = 0.0;
            for (const SiteObservation& obs : observations) {
                const double score_r =
                    (site_log_likelihood(obs, r + hr, c, rho, use_beta) -
                     site_log_likelihood(obs, r - hr, c, rho, use_beta)) / (2.0 * hr);
                const double score_c =
                    (site_log_likelihood(obs, r, c + hc, rho, use_beta) -
                     site_log_likelihood(obs, r, c - hc, rho, use_beta)) / (2.0 * hc);
                b_rr += score_r * score_r;
                b_rc += score_r * score_c;
                b_cc += score_c * score_c;
            }
            const double v_rr = inv_rr * (b_rr * inv_rr + b_rc * inv_rc) +
                                inv_rc * (b_rc * inv_rr + b_cc * inv_rc);
            const double v_cc = inv_rc * (b_rr * inv_rc + b_rc * inv_cc) +
                                inv_cc * (b_rc * inv_rc + b_cc * inv_cc);
            fit.ratio_se_robust = std::sqrt(std::max(0.0, v_rr));
            fit.ambient_se_robust = std::sqrt(std::max(0.0, v_cc));
        }
    }

    if (use_beta) {
        const double h_rho = finite_step(rho, 1e-8, opt.overdispersion_max);
        if (h_rho > 0.0) {
            const double ll_plus = log_likelihood(observations, r, c, rho + h_rho, true);
            const double ll_minus = log_likelihood(observations, r, c, rho - h_rho, true);
            const double info = -(ll_plus - 2.0 * center + ll_minus) /
                                (h_rho * h_rho);
            if (info > kTiny && std::isfinite(info)) fit.overdispersion_se = std::sqrt(1.0 / info);
        }
    }
}

FitResult fit_cell(const std::vector<SiteObservation>& observations,
                   const Options& opt,
                   bool ambient_enabled,
                   const double* fixed_ambient) {
    FitResult fit = fit_parameters(observations, opt, ambient_enabled, fixed_ambient);
    compute_uncertainty(observations, opt, fit,
                        ambient_enabled && fixed_ambient == nullptr);
    return fit;
}

std::string format_double(double value, int precision = 8) {
    if (!std::isfinite(value)) return "NA";
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
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
            opt.assignments, sample_index, barcode_to_cell);
        std::unordered_set<std::string> empty_barcodes = load_empty_barcodes(
            opt.empty_barcodes, barcode_to_cell);
        std::unordered_map<std::string, double> fixed_ambient =
            load_fixed_ambient_fractions(opt.ambient_fraction_file);

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
                "Fixed ambient fractions loaded: %zu (explicit per-cell override)\n",
                fixed_ambient.size());
        }

        ReadStats stats;
        count_molecules(opt, sites, position_to_site, barcode_to_cell,
                        empty_barcodes, cells, ambient, stats);
        if (stats.pileup_bases == 0) {
            throw std::runtime_error(
                "The BAM contains the mitochondrial contig but no reads overlap the selected "
                "mitochondrial panel. This commonly means a nuclear-only variant-site subset BAM "
                "was supplied; rebuild the subset BED with the emitted mt sites BED.");
        }
        if (stats.accepted_observations == 0) {
            throw std::runtime_error(
                "No UMI-collapsed mitochondrial molecule/site observations were accepted. "
                "Check CB/UB tags, barcode agreement, mapping filters, and panel overlap.");
        }

        if (!empty_barcodes.empty()) {
            for (AmbientSite& a : ambient) {
                const uint64_t total = a.total();
                a.alt_fraction = (static_cast<double>(a.alt) + 0.5) /
                                 (static_cast<double>(total) + 1.0);
                a.usable = static_cast<int>(total) >= opt.ambient_min_molecules;
            }
        }

        const bool ambient_enabled = !opt.ambient_none;
        int usable_ambient_sites = 0;
        if (ambient_enabled) {
            for (const AmbientSite& a : ambient) if (a.usable) ++usable_ambient_sites;
            if (usable_ambient_sites == 0) {
                throw std::runtime_error(
                    "No mt-panel sites have a usable ambient profile; provide more empty-droplet "
                    "evidence, lower --ambient_min_molecules, or explicitly use --ambient_none");
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
                        << (!ambient_enabled ? "NOT_USED" :
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
                  << "\titerations\tconverged\tfit_mode\tstatus\n";

        std::ofstream site_out;
        if (opt.write_site_counts) {
            const std::string site_path = opt.output_prefix + ".mt_site_counts.tsv";
            site_out.open(site_path);
            if (!site_out) throw std::runtime_error("Could not write: " + site_path);
            site_out << "library_id\tbarcode\tidentity\tchrom\tpos\tref\talt\tsite_class"
                     << "\tparent1_state\tparent2_state\tratio_informative\tambient_anchor"
                     << "\tref_molecules\talt_molecules\tambient_ref_molecules"
                     << "\tambient_alt_molecules\tambient_alt_fraction\tused_in_fit\n";
        }

        uint64_t pass_cells = 0;
        uint64_t weak_cells = 0;
        uint64_t low_cells = 0;
        uint64_t fixed_cells = 0;
        uint64_t pair_missing_cells = 0;
        uint64_t collapsed_pair_cells = 0;
        uint64_t unanchored_cells = 0;

        for (CellInfo& cell : cells) {
            std::string status = "PASS";
            std::string fit_mode;
            if (opt.legacy_panel_gt) fit_mode = "LEGACY_GT_";
            fit_mode += opt.likelihood == "beta_binomial" ? "BETA_BINOMIAL_" : "BINOMIAL_";
            fit_mode += ambient_enabled ? "JOINT_MT_AMBIENT" : "NO_AMBIENT";

            const bool same_parent = cell.parent1_index == cell.parent2_index;
            const PairManifest* pair_manifest = nullptr;
            PairManifest legacy_manifest;
            if (!same_parent) {
                if (opt.legacy_panel_gt) {
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
                    const auto found = pair_manifests.find(
                        pair_manifest_key(cell.parent1, cell.parent2));
                    if (found != pair_manifests.end()) pair_manifest = &found->second;
                }
            }

            const auto fixed_it = fixed_ambient.find(cell.barcode);
            const double* fixed_value = fixed_it == fixed_ambient.end() ? nullptr : &fixed_it->second;
            if (fixed_value) {
                fit_mode = (opt.legacy_panel_gt ? "LEGACY_GT_" : "") +
                           std::string(opt.likelihood == "beta_binomial" ?
                                       "BETA_BINOMIAL_FIXED_AMBIENT" : "BINOMIAL_FIXED_AMBIENT");
                ++fixed_cells;
            }
            const bool joint_ambient = ambient_enabled && fixed_value == nullptr;

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
            std::vector<SiteObservation> observations;
            std::unordered_set<int> used_site_indices;

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
                    const bool ambient_anchor =
                        manifest_site.site_class == MANIFEST_AMBIENT_ONLY;
                    if (count.total() == 0) continue;
                    if (ratio_informative) {
                        ++ratio_sites_observed;
                        ratio_molecules += count.total();
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

                    const bool ambient_ok = !ambient_enabled || ambient[site_index].usable;
                    const bool class_used = ratio_informative || (ambient_enabled && ambient_anchor);
                    if (!ambient_ok || !class_used) continue;
                    SiteObservation obs;
                    obs.site_index = site_index;
                    obs.ref = count.ref;
                    obs.alt = count.alt;
                    obs.parent1_alt_probability = p1_state == 1
                        ? 1.0 - opt.error_rate : opt.error_rate;
                    obs.parent2_alt_probability = p2_state == 1
                        ? 1.0 - opt.error_rate : opt.error_rate;
                    obs.ambient_alt_probability = ambient_enabled
                        ? clamp_probability(ambient[site_index].alt_fraction) : 0.5;
                    obs.site_class = manifest_site.site_class;
                    obs.ratio_informative = ratio_informative;
                    obs.ambient_anchor = ambient_anchor;
                    observations.push_back(obs);
                    used_site_indices.insert(site_index);
                    total_molecules_used += count.total();
                    if (ratio_informative) ++ratio_sites_used;
                    if (ambient_anchor) ++ambient_sites_used;
                }
            }

            if (status == "PASS" && ratio_sites_observed == 0) {
                status = "NO_RATIO_SITES_OBSERVED";
            } else if (status == "PASS" &&
                       ratio_molecules < static_cast<uint64_t>(opt.min_molecules)) {
                status = "LOW_RATIO_MOLECULES";
            } else if (status == "PASS" && ratio_sites_used < opt.min_sites) {
                status = ambient_enabled ? "LOW_AMBIENT_COVERED_RATIO_SITES" :
                                           "LOW_RATIO_SITES";
            }

            const bool missing_anchor = joint_ambient &&
                ambient_sites_used < opt.min_ambient_only_sites;
            if (status == "PASS" && missing_anchor && !opt.allow_unanchored_ambient) {
                status = "NO_AMBIENT_ANCHOR_OBSERVED";
                ++unanchored_cells;
            }

            FitResult fit;
            if (status == "PASS") {
                fit = fit_cell(observations, opt, ambient_enabled, fixed_value);
                // Missing ambient-only evidence is a structural limitation, not
                // merely an optimizer outcome.  When the user explicitly asks
                // for this diagnostic fit, preserve that scientific status even
                // if the optimizer drifts along the expected shallow ridge.
                if (missing_anchor) {
                    status = "WEAKLY_IDENTIFIABLE_UNANCHORED";
                    ++unanchored_cells;
                } else if (!fit.converged) {
                    status = "NOT_CONVERGED";
                } else if (fit.joint_ambient &&
                           (fit.information_condition > 1e8 ||
                            fit.min_information_eigenvalue < 1e-4 ||
                            !std::isfinite(fit.ratio_se))) {
                    status = "WEAKLY_IDENTIFIABLE";
                } else if (fit.ambient_weight > 0.98) {
                    status = "PARENTAL_SIGNAL_ABSENT";
                } else if (opt.likelihood == "beta_binomial" &&
                           fit.overdispersion_rho >= 0.999 * opt.overdispersion_max) {
                    status = "OVERDISPERSION_AT_BOUND";
                }
            }

            if (status == "PASS") {
                ++pass_cells;
            } else if (status == "WEAKLY_IDENTIFIABLE" ||
                       status == "WEAKLY_IDENTIFIABLE_UNANCHORED" ||
                       status == "PARENTAL_SIGNAL_ABSENT" ||
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

            ratio_out << (opt.legacy_panel_gt ? "NA" : opt.library_id) << '\t'
                      << cell.barcode << '\t' << cell.identity << '\t'
                      << cell.parent1 << '\t' << cell.parent2 << '\t'
                      << opt.likelihood << '\t'
                      << format_double(parent1_fraction) << '\t'
                      << format_double(fit.parent2_fraction) << '\t'
                      << format_double(fit.ratio_se) << '\t'
                      << format_double(fit.ratio_se_robust) << '\t'
                      << format_double(ci_low) << '\t' << format_double(ci_high) << '\t'
                      << format_double(fit.ambient_weight) << '\t'
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
                      << fit_mode << '\t' << status << '\n';

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
                    const bool ambient_anchor =
                        manifest_site.site_class == MANIFEST_AMBIENT_ONLY;
                    const AlleleCount& count = item.second;
                    const MtSite& site = sites[site_index];
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
                             << (used_site_indices.count(site_index) ? 1 : 0) << '\n';
                }
            }
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
        qc << "fusion_cells\t" << cells.size() << '\n';
        qc << "empty_barcodes\t" << empty_barcodes.size() << '\n';
        qc << "usable_ambient_sites\t" << usable_ambient_sites << '\n';
        qc << "fixed_ambient_cells\t" << fixed_cells << '\n';
        qc << "pass_cells\t" << pass_cells << '\n';
        qc << "weak_cells\t" << weak_cells << '\n';
        qc << "low_or_unfit_cells\t" << low_cells << '\n';
        qc << "pair_missing_cells\t" << pair_missing_cells << '\n';
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
        qc << "ambient_mode\t" << (opt.ambient_none ? "NONE" :
            (!opt.empty_barcodes.empty() ? "EMPTY_BARCODES" : "PROFILE_FILE")) << '\n';
        qc << "min_mapq\t" << opt.min_mapq << '\n';
        qc << "min_baseq\t" << opt.min_baseq << '\n';
        qc << "error_rate\t" << opt.error_rate << '\n';
        qc << "ambiguous_alignments_allowed\t" << (opt.allow_ambiguous_alignments ? 1 : 0) << '\n';
        qc << "duplicates_kept\t" << (opt.keep_duplicates ? 1 : 0) << '\n';

        std::fprintf(stderr, "\nMitochondrial fusion-ratio estimation complete\n");
        std::fprintf(stderr, "  cells: pass=%llu weak=%llu other=%llu\n",
            static_cast<unsigned long long>(pass_cells),
            static_cast<unsigned long long>(weak_cells),
            static_cast<unsigned long long>(low_cells));
        std::fprintf(stderr, "  accepted molecule-site observations: %llu\n",
            static_cast<unsigned long long>(stats.accepted_observations));
        std::fprintf(stderr, "  ratio table: %s\n", ratio_path.c_str());
        std::fprintf(stderr, "  ambient profile: %s\n", ambient_path.c_str());
        std::fprintf(stderr, "  QC: %s\n", qc_path.c_str());
        return 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "ERROR: %s\n", error.what());
        return 1;
    }
}
