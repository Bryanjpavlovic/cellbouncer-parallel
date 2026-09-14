#include <getopt.h>
#include <argp.h>
#include <string>
#include <algorithm>
#include <array>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <bitset>
#include <deque>
#include <cstdlib>
#include <random>
#include <functional>
#include <utility>
#include <cmath>
#include <float.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include <htslib/kseq.h>
#include <zlib.h>
#include <htswrapper/bc.h>
#include <htswrapper/bam.h>
#include <htswrapper/gzreader.h>
#include <mixtureDist/functions.h>
#include <optimML/multivar_ml.h>
#include <optimML/brent.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include <iomanip>
#include <numeric>
#include "common.h"
#include "downsample_vcf_parallel.h"
// ============================================================================
// MITOCHONDRIAL PANEL SELECTION (integrated)
//
// This implementation is intentionally compiled into downsample_vcf_parallel.
// There is no standalone panel-selection executable or second VCF scan.
// Supplying --mt_output enables the optional mitochondrial panel during the
// same source-VCF load used for the nuclear panels.
// ============================================================================
#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>

namespace mtpanel {

struct Config {
    std::string output;
    std::string pool_combinations;
    std::string mito_chrom = "chrM";
    std::string coverage;
    std::string blacklist_bed;
    std::string mt_mask_bed;
    std::string pair_audit;
    std::string site_manifest;
    std::string haplotype_groups;
    std::string haplotype_pairwise;
    std::string sample_audit;
    std::string sites_bed;
    int threads = 4;
    int min_depth = 10;
    int min_pair_sites = 20;
    int min_ambient_sites = 1;
    double homoplasmy_af = 0.95;
    double min_coverage = 0.0;
    double min_qual = 0.0;
    bool require_pair_targets = false;
    bool require_coverage = true;
    bool all_pairs = false;
};

struct Summary {
    int64_t records_scanned = 0;
    int64_t selected_sites = 0;
    int64_t rejected_blacklist = 0;
    int64_t rejected_mt_mask = 0;
    int64_t rejected_non_snv = 0;
    int64_t rejected_quality = 0;
    int64_t rejected_coverage = 0;
    int64_t rejected_no_use = 0;
    int libraries = 0;
    int expected_donors = 0;
    int expected_pairs = 0;
    int pairs_below_ratio_target = 0;
    int pairs_without_ambient_anchor = 0;
};

// Build the production mitochondrial panel from records already read during the
// downsample_vcf_parallel input pass. Records are borrowed and remain owned by
// the caller. Only records on Config::mito_chrom should normally be supplied.
Summary build_panel(bcf_hdr_t* header,
                    const std::vector<bcf1_t*>& records,
                    const Config& config);

}  // namespace mtpanel

namespace mtpanel {
namespace {

struct Interval {
    int64_t start;
    int64_t end;
    Interval() : start(0), end(0) {}
    Interval(int64_t start_in, int64_t end_in) : start(start_in), end(end_in) {}
};

struct CoverageInterval {
    int64_t start;
    int64_t end;
    double value;
    CoverageInterval() : start(0), end(0), value(0.0) {}
    CoverageInterval(int64_t start_in, int64_t end_in, double value_in)
        : start(start_in), end(end_in), value(value_in) {}
};

struct PairDef {
    std::string first;
    std::string second;
    int first_index = -1;
    int second_index = -1;
};

struct LibraryDef {
    std::string id;
    std::vector<std::string> donors;
    std::vector<int> donor_indices;
    std::vector<PairDef> pairs;
};

struct Evidence {
    int8_t state = -1;  // 0 REF, 1 ALT, -1 not a clean homoplasmic call
    int ref_depth = -1;
    int alt_depth = -1;
    int total_depth = -1;
    double alt_fraction = std::numeric_limits<double>::quiet_NaN();
    std::string depth_source;
    std::string failure;
};

enum SiteClass {
    CLASS_A_ALT_B_REF = 0,
    CLASS_A_REF_B_ALT = 1,
    CLASS_AMBIENT_ONLY = 2
};

const char* class_name(SiteClass value) {
    switch (value) {
        case CLASS_A_ALT_B_REF: return "A_ALT_B_REF";
        case CLASS_A_REF_B_ALT: return "A_REF_B_ALT";
        case CLASS_AMBIENT_ONLY: return "AMBIENT_ONLY";
    }
    return "UNKNOWN";
}

struct SiteUse {
    int library_index = -1;
    int pair_index = -1;
    SiteClass site_class = CLASS_AMBIENT_ONLY;
    std::vector<std::string> opposite_donors;
};

struct Candidate {
    bcf1_t* record = nullptr;  // borrowed
    int64_t pos = -1;
    char ref = 'N';
    char alt = 'N';
    double qual = 0.0;
    double coverage = 0.0;
    std::vector<Evidence> evidence;
    std::vector<SiteUse> uses;
};

struct PairCounts {
    int a_alt_b_ref = 0;
    int a_ref_b_alt = 0;
    int ambient_only = 0;
};

std::string trim(std::string value) {
    const char* ws = " \t\r\n";
    const std::string::size_type first = value.find_first_not_of(ws);
    if (first == std::string::npos) return "";
    const std::string::size_type last = value.find_last_not_of(ws);
    return value.substr(first, last - first + 1);
}


std::string canonical_library_id(std::string value) {
    value = trim(value);
    if (value.empty()) return value;

    size_t digit_start = 0;
    if (value.size() >= 3 &&
        std::tolower(static_cast<unsigned char>(value[0])) == 'l' &&
        std::tolower(static_cast<unsigned char>(value[1])) == 'i' &&
        std::tolower(static_cast<unsigned char>(value[2])) == 'b') {
        digit_start = 3;
    }
    if (digit_start >= value.size()) return value;

    for (size_t i = digit_start; i < value.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(value[i]))) return value;
    }

    std::string digits = value.substr(digit_start);
    const size_t first_nonzero = digits.find_first_not_of('0');
    if (first_nonzero == std::string::npos) digits = "0";
    else if (first_nonzero > 0) digits.erase(0, first_nonzero);
    return "lib" + digits;
}

std::vector<std::string> split(const std::string& value, char delim) {
    std::vector<std::string> out;
    std::stringstream ss(value);
    std::string item;
    while (std::getline(ss, item, delim)) out.push_back(item);
    return out;
}

std::string join(const std::vector<std::string>& values, const char* delim) {
    std::ostringstream out;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) out << delim;
        out << values[i];
    }
    return out.str();
}

std::string default_path(const std::string& configured,
                         const std::string& output,
                         const std::string& suffix) {
    return configured.empty() ? output + suffix : configured;
}

bool valid_base_string(const char* allele) {
    if (!allele || allele[0] == '\0' || allele[1] != '\0') return false;
    const char b = static_cast<char>(std::toupper(static_cast<unsigned char>(allele[0])));
    return b == 'A' || b == 'C' || b == 'G' || b == 'T';
}

std::vector<Interval> load_bed_for_contig(const std::string& filename,
                                          const std::string& contig) {
    std::vector<Interval> result;
    if (filename.empty()) return result;
    std::ifstream in(filename.c_str());
    if (!in) throw std::runtime_error("Could not open chrM mask BED: " + filename);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::string chrom;
        int64_t start = 0, end = 0;
        if (!(ss >> chrom >> start >> end)) continue;
        if (chrom == contig && start >= 0 && end > start) result.push_back(Interval{start, end});
    }
    std::sort(result.begin(), result.end(), [](const Interval& a, const Interval& b) {
        return a.start < b.start || (a.start == b.start && a.end < b.end);
    });
    std::vector<Interval> merged;
    for (size_t i = 0; i < result.size(); ++i) {
        if (merged.empty() || result[i].start > merged.back().end) {
            merged.push_back(result[i]);
        } else if (result[i].end > merged.back().end) {
            merged.back().end = result[i].end;
        }
    }
    return merged;
}

bool in_intervals(int64_t pos, const std::vector<Interval>& intervals) {
    std::vector<Interval>::const_iterator it = std::upper_bound(
        intervals.begin(), intervals.end(), pos,
        [](int64_t value, const Interval& interval) { return value < interval.start; });
    if (it == intervals.begin()) return false;
    --it;
    return it->start <= pos && pos < it->end;
}

std::vector<CoverageInterval> load_coverage_for_contig(const std::string& filename,
                                                       const std::string& contig) {
    std::vector<CoverageInterval> result;
    if (filename.empty()) return result;
    const bool compressed =
        (filename.size() >= 3 && filename.substr(filename.size() - 3) == ".gz") ||
        (filename.size() >= 4 && filename.substr(filename.size() - 4) == ".bgz");
    if (compressed) {
        tbx_t* tbx = tbx_index_load(filename.c_str());
        if (!tbx) throw std::runtime_error(
            "Compressed coverage track requires a tabix index: " + filename);
        htsFile* fp = hts_open(filename.c_str(), "r");
        if (!fp) {
            tbx_destroy(tbx);
            throw std::runtime_error("Could not open coverage track: " + filename);
        }
        hts_itr_t* itr = tbx_itr_querys(tbx, contig.c_str());
        if (!itr) {
            hts_close(fp);
            tbx_destroy(tbx);
            return result;
        }
        kstring_t line = {0, 0, NULL};
        while (tbx_itr_next(fp, tbx, itr, &line) >= 0) {
            std::stringstream ss(std::string(line.s, line.l));
            std::string chrom;
            int64_t start = 0, end = 0;
            double value = 0.0;
            if (ss >> chrom >> start >> end >> value) {
                if (chrom == contig && start >= 0 && end > start) {
                    result.push_back(CoverageInterval{start, end, value});
                }
            }
        }
        std::free(line.s);
        hts_itr_destroy(itr);
        hts_close(fp);
        tbx_destroy(tbx);
    } else {
        std::ifstream in(filename.c_str());
        if (!in) throw std::runtime_error("Could not open coverage track: " + filename);
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::stringstream ss(line);
            std::string chrom;
            int64_t start = 0, end = 0;
            double value = 0.0;
            if (ss >> chrom >> start >> end >> value) {
                if (chrom == contig && start >= 0 && end > start) {
                    result.push_back(CoverageInterval{start, end, value});
                }
            }
        }
    }
    std::sort(result.begin(), result.end(), [](const CoverageInterval& a,
                                               const CoverageInterval& b) {
        return a.start < b.start || (a.start == b.start && a.end < b.end);
    });
    return result;
}

double coverage_at(int64_t pos, const std::vector<CoverageInterval>& coverage) {
    if (coverage.empty()) return 0.0;
    std::vector<CoverageInterval>::const_iterator it = std::upper_bound(
        coverage.begin(), coverage.end(), pos,
        [](int64_t value, const CoverageInterval& interval) {
            return value < interval.start;
        });
    if (it == coverage.begin()) return 0.0;
    --it;
    return (it->start <= pos && pos < it->end) ? it->value : 0.0;
}

void add_unique(std::vector<std::string>& values,
                std::unordered_set<std::string>& seen,
                const std::string& value) {
    if (!value.empty() && seen.insert(value).second) values.push_back(value);
}

std::vector<LibraryDef> load_libraries(const Config& config,
                                       const std::vector<std::string>& samples,
                                       std::set<std::string>& expected_donors) {
    std::unordered_map<std::string, int> sample_index;
    for (size_t i = 0; i < samples.size(); ++i) sample_index[samples[i]] = static_cast<int>(i);

    std::vector<LibraryDef> libraries;
    if (config.all_pairs) {
        LibraryDef lib;
        lib.id = "ALL_VCF_SAMPLES";
        lib.donors = samples;
        for (size_t i = 0; i < samples.size(); ++i) {
            lib.donor_indices.push_back(static_cast<int>(i));
            expected_donors.insert(samples[i]);
        }
        for (size_t i = 0; i < samples.size(); ++i) {
            for (size_t j = i + 1; j < samples.size(); ++j) {
                PairDef pair;
                pair.first = samples[i];
                pair.second = samples[j];
                pair.first_index = static_cast<int>(i);
                pair.second_index = static_cast<int>(j);
                lib.pairs.push_back(pair);
            }
        }
        libraries.push_back(lib);
        return libraries;
    }

    std::ifstream in(config.pool_combinations.c_str());
    if (!in) throw std::runtime_error(
        "Could not open metadata-derived pool combinations: " + config.pool_combinations);
    std::string line;
    bool first_data_line = true;
    std::set<std::string> seen_library_ids;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const std::string::size_type tab = line.find('\t');
        if (tab == std::string::npos) continue;
        const std::string raw_library_id = trim(line.substr(0, tab));
        const std::string identities_field = trim(line.substr(tab + 1));
        if (first_data_line && raw_library_id == "library_id") {
            first_data_line = false;
            continue;
        }
        first_data_line = false;
        const std::string library_id = canonical_library_id(raw_library_id);
        if (library_id.empty() || identities_field.empty()) continue;
        if (!seen_library_ids.insert(library_id).second) {
            throw std::runtime_error("Duplicate library_id in pool combinations: " + library_id);
        }

        LibraryDef lib;
        lib.id = library_id;
        std::unordered_set<std::string> donor_seen;
        std::set<std::pair<std::string, std::string> > pair_seen;
        const std::vector<std::string> identities = split(identities_field, ',');
        for (size_t i = 0; i < identities.size(); ++i) {
            const std::string identity = trim(identities[i]);
            if (identity.empty()) continue;
            const std::vector<std::string> parts_raw = split(identity, '+');
            std::vector<std::string> parts;
            for (size_t j = 0; j < parts_raw.size(); ++j) {
                const std::string donor = trim(parts_raw[j]);
                if (!donor.empty()) parts.push_back(donor);
            }
            if (parts.empty() || parts.size() > 2) {
                throw std::runtime_error("Unsupported identity in library " + library_id + ": " + identity);
            }
            for (size_t j = 0; j < parts.size(); ++j) add_unique(lib.donors, donor_seen, parts[j]);
            if (parts.size() == 2 && parts[0] != parts[1]) {
                std::string a = parts[0];
                std::string b = parts[1];
                if (a > b) std::swap(a, b);
                if (pair_seen.insert(std::make_pair(a, b)).second) {
                    PairDef pair;
                    pair.first = a;
                    pair.second = b;
                    lib.pairs.push_back(pair);
                }
            }
        }
        if (lib.donors.empty()) {
            throw std::runtime_error("No donors resolved for library: " + library_id);
        }
        for (size_t i = 0; i < lib.donors.size(); ++i) {
            const std::unordered_map<std::string, int>::const_iterator it = sample_index.find(lib.donors[i]);
            if (it == sample_index.end()) {
                throw std::runtime_error("Metadata donor is absent from VCF samples: " + lib.donors[i] +
                                         " (library " + library_id + ")");
            }
            lib.donor_indices.push_back(it->second);
            expected_donors.insert(lib.donors[i]);
        }
        for (size_t i = 0; i < lib.pairs.size(); ++i) {
            lib.pairs[i].first_index = sample_index[lib.pairs[i].first];
            lib.pairs[i].second_index = sample_index[lib.pairs[i].second];
        }
        libraries.push_back(lib);
    }
    if (libraries.empty()) throw std::runtime_error("No library definitions were loaded");
    return libraries;
}

bool valid_format_value(int32_t value) {
    return value != bcf_int32_missing && value != bcf_int32_vector_end && value >= 0;
}

int format_value(const int32_t* data, int total_values, int n_samples,
                 int sample, int offset) {
    if (!data || total_values <= 0 || n_samples <= 0 || total_values % n_samples != 0) return -1;
    const int stride = total_values / n_samples;
    if (offset < 0 || offset >= stride) return -1;
    const int32_t value = data[sample * stride + offset];
    return valid_format_value(value) ? static_cast<int>(value) : -1;
}

int gt_homoplasmic_state(const int32_t* gt, int total_values, int n_samples, int sample) {
    if (!gt || total_values <= 0 || total_values % n_samples != 0) return -1;
    const int stride = total_values / n_samples;
    int state = -1;
    bool saw = false;
    for (int k = 0; k < stride; ++k) {
        const int32_t value = gt[sample * stride + k];
        if (value == bcf_int32_vector_end) break;
        if (bcf_gt_is_missing(value)) return -1;
        const int allele = bcf_gt_allele(value);
        if (allele < 0 || allele > 1) return -3;
        if (!saw) {
            state = allele;
            saw = true;
        } else if (state != allele) {
            return -2;
        }
    }
    return saw ? state : -1;
}

std::vector<Evidence> extract_evidence(bcf_hdr_t* header,
                                       bcf1_t* record,
                                       const Config& config) {
    const int n_samples = bcf_hdr_nsamples(header);
    int32_t* ro = NULL;
    int32_t* ao = NULL;
    int32_t* ad = NULL;
    int32_t* gt = NULL;
    int n_ro = 0, n_ao = 0, n_ad = 0, n_gt = 0;
    const int ret_ro = bcf_get_format_int32(header, record, "RO", &ro, &n_ro);
    const int ret_ao = bcf_get_format_int32(header, record, "AO", &ao, &n_ao);
    const int ret_ad = bcf_get_format_int32(header, record, "AD", &ad, &n_ad);
    const int ret_gt = bcf_get_genotypes(header, record, &gt, &n_gt);

    std::vector<Evidence> result(n_samples);
    for (int i = 0; i < n_samples; ++i) {
        Evidence evidence;
        int ref_depth = -1;
        int alt_depth = -1;
        if (ret_ro > 0 && ret_ao > 0) {
            ref_depth = format_value(ro, ret_ro, n_samples, i, 0);
            alt_depth = format_value(ao, ret_ao, n_samples, i, 0);
            if (ref_depth >= 0 && alt_depth >= 0) evidence.depth_source = "RO_AO";
        }
        if ((ref_depth < 0 || alt_depth < 0) && ret_ad > 0) {
            ref_depth = format_value(ad, ret_ad, n_samples, i, 0);
            alt_depth = format_value(ad, ret_ad, n_samples, i, 1);
            if (ref_depth >= 0 && alt_depth >= 0) evidence.depth_source = "AD";
        }
        evidence.ref_depth = ref_depth;
        evidence.alt_depth = alt_depth;
        if (ref_depth < 0 || alt_depth < 0) {
            evidence.failure = "NO_ALLELE_DEPTH";
            result[i] = evidence;
            continue;
        }
        evidence.total_depth = ref_depth + alt_depth;
        if (evidence.total_depth < config.min_depth) {
            evidence.failure = "LOW_DEPTH";
            result[i] = evidence;
            continue;
        }
        evidence.alt_fraction = static_cast<double>(alt_depth) /
                                static_cast<double>(evidence.total_depth);
        int depth_state = -1;
        if (evidence.alt_fraction <= 1.0 - config.homoplasmy_af) depth_state = 0;
        if (evidence.alt_fraction >= config.homoplasmy_af) depth_state = 1;
        if (depth_state < 0) {
            evidence.failure = "MIXED_ALLELE_FRACTION";
            result[i] = evidence;
            continue;
        }
        const int gt_state = gt_homoplasmic_state(gt, ret_gt, n_samples, i);
        if (gt_state == -2) {
            evidence.failure = "HETEROZYGOUS_GT";
            result[i] = evidence;
            continue;
        }
        if (gt_state == -3) {
            evidence.failure = "UNSUPPORTED_GT_ALLELE";
            result[i] = evidence;
            continue;
        }
        if (gt_state >= 0 && gt_state != depth_state) {
            evidence.failure = "GT_DEPTH_CONFLICT";
            result[i] = evidence;
            continue;
        }
        evidence.state = static_cast<int8_t>(depth_state);
        evidence.failure = "PASS";
        result[i] = evidence;
    }

    std::free(ro);
    std::free(ao);
    std::free(ad);
    std::free(gt);
    return result;
}


uint64_t fnv1a64(const std::string& text) {
    uint64_t hash = 1469598103934665603ULL;
    for (size_t i = 0; i < text.size(); ++i) {
        hash ^= static_cast<unsigned char>(text[i]);
        hash *= 1099511628211ULL;
    }
    return hash;
}

std::string hex64(uint64_t value) {
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << value;
    return out.str();
}

std::string output_mode(const std::string& output) {
    if (output.size() >= 7 && output.substr(output.size() - 7) == ".vcf.gz") return "wz";
    if (output.size() >= 4 && output.substr(output.size() - 4) == ".vcf") return "w";
    return "wb";
}

void write_sample_audit(const std::string& path,
                        const std::vector<std::string>& samples,
                        const std::vector<LibraryDef>& libraries,
                        const std::set<std::string>& expected_donors) {
    std::map<std::string, std::vector<std::string> > donor_libraries;
    for (size_t l = 0; l < libraries.size(); ++l) {
        for (size_t d = 0; d < libraries[l].donors.size(); ++d) {
            donor_libraries[libraries[l].donors[d]].push_back(libraries[l].id);
        }
    }
    std::ofstream out(path.c_str());
    if (!out) throw std::runtime_error("Could not write mt sample audit: " + path);
    out << "vcf_sample\texpected_in_any_library\tlibrary_count\tlibraries\tselection_role\n";
    for (size_t i = 0; i < samples.size(); ++i) {
        const bool expected = expected_donors.count(samples[i]) != 0;
        const std::vector<std::string>& libs = donor_libraries[samples[i]];
        out << samples[i] << '\t' << (expected ? 1 : 0) << '\t' << libs.size() << '\t'
            << join(libs, ",") << '\t'
            << (expected ? "METADATA_EXPECTED_DONOR" : "EXCLUDED_UNRELATED_VCF_SAMPLE") << '\n';
    }
}

}  // namespace

Summary build_panel(bcf_hdr_t* header,
                    const std::vector<bcf1_t*>& records,
                    const Config& config) {
    if (!header) throw std::runtime_error("Mito panel builder received a null VCF header");
    if (config.output.empty()) throw std::runtime_error("Mito panel output path is empty");
    if (!config.all_pairs && config.pool_combinations.empty()) {
        throw std::runtime_error("Mito panel requires metadata-derived --pool_combinations");
    }
    if (config.min_depth < 1 || config.homoplasmy_af <= 0.5 ||
        config.homoplasmy_af > 1.0 || config.min_pair_sites < 0 ||
        config.min_ambient_sites < 0 || config.min_coverage < 0.0 ||
        config.threads < 1) {
        throw std::runtime_error("Invalid mitochondrial panel threshold");
    }

    Summary summary;
    const int n_samples = bcf_hdr_nsamples(header);
    if (n_samples < 2) throw std::runtime_error("Mito panel requires at least two VCF samples");
    std::vector<std::string> samples;
    for (int i = 0; i < n_samples; ++i) samples.push_back(header->samples[i]);

    std::set<std::string> expected_donors;
    const std::vector<LibraryDef> libraries = load_libraries(config, samples, expected_donors);
    summary.libraries = static_cast<int>(libraries.size());
    summary.expected_donors = static_cast<int>(expected_donors.size());
    for (size_t l = 0; l < libraries.size(); ++l) {
        summary.expected_pairs += static_cast<int>(libraries[l].pairs.size());
    }
    if (summary.expected_pairs == 0) {
        throw std::runtime_error("No heterotypic fusion pairs were found in library metadata");
    }

    const int mt_rid = bcf_hdr_name2id(header, config.mito_chrom.c_str());
    if (mt_rid < 0) throw std::runtime_error(
        "Mitochondrial contig is absent from VCF header: " + config.mito_chrom);
    const std::vector<Interval> blacklist =
        load_bed_for_contig(config.blacklist_bed, config.mito_chrom);
    const std::vector<Interval> mt_mask =
        load_bed_for_contig(config.mt_mask_bed, config.mito_chrom);
    const std::vector<CoverageInterval> coverage =
        load_coverage_for_contig(config.coverage, config.mito_chrom);
    if (config.require_coverage && coverage.empty()) {
        throw std::runtime_error(
            "Mito panel selection requires a coverage track containing " + config.mito_chrom +
            "; use require_coverage=false only for an explicit diagnostic run");
    }

    const std::string pair_audit_path = default_path(
        config.pair_audit, config.output, ".pair_audit.tsv");
    const std::string site_manifest_path = default_path(
        config.site_manifest, config.output, ".site_manifest.tsv");
    const std::string haplotype_groups_path = default_path(
        config.haplotype_groups, config.output, ".haplotype_groups.tsv");
    const std::string haplotype_pairwise_path = default_path(
        config.haplotype_pairwise, config.output, ".haplotype_pairwise.tsv");
    const std::string sample_audit_path = default_path(
        config.sample_audit, config.output, ".sample_audit.tsv");
    const std::string sites_bed_path = default_path(
        config.sites_bed, config.output, ".sites.bed");

    write_sample_audit(sample_audit_path, samples, libraries, expected_donors);

    std::vector<Candidate> candidates;
    for (size_t r = 0; r < records.size(); ++r) {
        bcf1_t* record = records[r];
        if (!record) continue;
        ++summary.records_scanned;
        if (record->rid != mt_rid) continue;
        if (in_intervals(record->pos, blacklist)) {
            ++summary.rejected_blacklist;
            continue;
        }
        if (in_intervals(record->pos, mt_mask)) {
            ++summary.rejected_mt_mask;
            continue;
        }
        bcf_unpack(record, BCF_UN_STR | BCF_UN_FMT);
        if (record->n_allele != 2 || !valid_base_string(record->d.allele[0]) ||
            !valid_base_string(record->d.allele[1]) ||
            record->d.allele[0][0] == record->d.allele[1][0]) {
            ++summary.rejected_non_snv;
            continue;
        }
        const double qual = bcf_float_is_missing(record->qual) ? 0.0 : record->qual;
        if (qual < config.min_qual) {
            ++summary.rejected_quality;
            continue;
        }
        const double cov = coverage.empty() ? 0.0 : coverage_at(record->pos, coverage);
        if (!coverage.empty() && cov < config.min_coverage) {
            ++summary.rejected_coverage;
            continue;
        }

        Candidate candidate;
        candidate.record = record;
        candidate.pos = record->pos;
        candidate.ref = static_cast<char>(std::toupper(record->d.allele[0][0]));
        candidate.alt = static_cast<char>(std::toupper(record->d.allele[1][0]));
        candidate.qual = qual;
        candidate.coverage = cov;
        candidate.evidence = extract_evidence(header, record, config);

        for (size_t l = 0; l < libraries.size(); ++l) {
            const LibraryDef& lib = libraries[l];
            for (size_t p = 0; p < lib.pairs.size(); ++p) {
                const PairDef& pair = lib.pairs[p];
                const int a = candidate.evidence[pair.first_index].state;
                const int b = candidate.evidence[pair.second_index].state;
                if (a < 0 || b < 0) continue;
                SiteUse use;
                use.library_index = static_cast<int>(l);
                use.pair_index = static_cast<int>(p);
                if (a != b) {
                    use.site_class = (a == 1) ? CLASS_A_ALT_B_REF : CLASS_A_REF_B_ALT;
                    candidate.uses.push_back(use);
                    continue;
                }
                for (size_t d = 0; d < lib.donors.size(); ++d) {
                    const std::string& donor = lib.donors[d];
                    if (donor == pair.first || donor == pair.second) continue;
                    const int donor_state = candidate.evidence[lib.donor_indices[d]].state;
                    if (donor_state >= 0 && donor_state != a) use.opposite_donors.push_back(donor);
                }
                if (!use.opposite_donors.empty()) {
                    use.site_class = CLASS_AMBIENT_ONLY;
                    candidate.uses.push_back(use);
                }
            }
        }
        if (candidate.uses.empty()) {
            ++summary.rejected_no_use;
            continue;
        }
        candidates.push_back(candidate);
    }

    std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
        if (a.record->rid != b.record->rid) return a.record->rid < b.record->rid;
        return a.pos < b.pos;
    });
    summary.selected_sites = static_cast<int64_t>(candidates.size());
    if (candidates.empty()) {
        throw std::runtime_error("No chrM sites passed depth, homoplasmy, expression, and library-use filters");
    }

    std::vector<std::vector<PairCounts> > pair_counts(libraries.size());
    for (size_t l = 0; l < libraries.size(); ++l) pair_counts[l].resize(libraries[l].pairs.size());
    for (size_t c = 0; c < candidates.size(); ++c) {
        for (size_t u = 0; u < candidates[c].uses.size(); ++u) {
            const SiteUse& use = candidates[c].uses[u];
            PairCounts& counts = pair_counts[use.library_index][use.pair_index];
            if (use.site_class == CLASS_A_ALT_B_REF) ++counts.a_alt_b_ref;
            else if (use.site_class == CLASS_A_REF_B_ALT) ++counts.a_ref_b_alt;
            else ++counts.ambient_only;
        }
    }

    {
        std::ofstream out(pair_audit_path.c_str());
        if (!out) throw std::runtime_error("Could not write mt pair audit: " + pair_audit_path);
        out << "library_id\tparent1\tparent2\tpair\ta_alt_b_ref_sites\ta_ref_b_alt_sites"
            << "\tratio_sites\tambient_only_sites\trequested_min_ratio_sites"
            << "\trequested_min_ambient_sites\tstatus\n";
        for (size_t l = 0; l < libraries.size(); ++l) {
            for (size_t p = 0; p < libraries[l].pairs.size(); ++p) {
                const PairDef& pair = libraries[l].pairs[p];
                const PairCounts& counts = pair_counts[l][p];
                const int ratio_total = counts.a_alt_b_ref + counts.a_ref_b_alt;
                std::vector<std::string> failures;
                if (ratio_total < config.min_pair_sites) {
                    failures.push_back(ratio_total == 0 ? "NO_RATIO_MARKERS" : "BELOW_RATIO_TARGET");
                    ++summary.pairs_below_ratio_target;
                }
                if (counts.ambient_only < config.min_ambient_sites) {
                    failures.push_back("NO_AMBIENT_ANCHOR");
                    ++summary.pairs_without_ambient_anchor;
                }
                out << libraries[l].id << '\t' << pair.first << '\t' << pair.second << '\t'
                    << pair.first << '+' << pair.second << '\t'
                    << counts.a_alt_b_ref << '\t' << counts.a_ref_b_alt << '\t'
                    << ratio_total << '\t' << counts.ambient_only << '\t'
                    << config.min_pair_sites << '\t' << config.min_ambient_sites << '\t'
                    << (failures.empty() ? "PASS" : join(failures, ";")) << '\n';
            }
        }
    }

    {
        std::ofstream out(site_manifest_path.c_str());
        if (!out) throw std::runtime_error("Could not write mt site manifest: " + site_manifest_path);
        out << "library_id\tparent1\tparent2\tchrom\tpos\tref\talt\tsite_class"
            << "\tparent1_state\tparent2_state\topposite_donors"
            << "\tparent1_depth\tparent1_alt_fraction\tparent1_depth_source"
            << "\tparent2_depth\tparent2_alt_fraction\tparent2_depth_source"
            << "\tcoverage\tqual\n";
        out << std::setprecision(10);
        for (size_t c = 0; c < candidates.size(); ++c) {
            const Candidate& candidate = candidates[c];
            for (size_t u = 0; u < candidate.uses.size(); ++u) {
                const SiteUse& use = candidate.uses[u];
                const LibraryDef& lib = libraries[use.library_index];
                const PairDef& pair = lib.pairs[use.pair_index];
                const Evidence& a = candidate.evidence[pair.first_index];
                const Evidence& b = candidate.evidence[pair.second_index];
                out << lib.id << '\t' << pair.first << '\t' << pair.second << '\t'
                    << config.mito_chrom << '\t' << candidate.pos + 1 << '\t'
                    << candidate.ref << '\t' << candidate.alt << '\t'
                    << class_name(use.site_class) << '\t' << static_cast<int>(a.state) << '\t'
                    << static_cast<int>(b.state) << '\t' << join(use.opposite_donors, ",") << '\t'
                    << a.total_depth << '\t' << a.alt_fraction << '\t' << a.depth_source << '\t'
                    << b.total_depth << '\t' << b.alt_fraction << '\t' << b.depth_source << '\t'
                    << candidate.coverage << '\t' << candidate.qual << '\n';
            }
        }
    }

    {
        std::ofstream groups(haplotype_groups_path.c_str());
        std::ofstream pairwise(haplotype_pairwise_path.c_str());
        if (!groups) throw std::runtime_error("Could not write mt haplotype groups: " + haplotype_groups_path);
        if (!pairwise) throw std::runtime_error("Could not write mt haplotype pairwise report: " + haplotype_pairwise_path);
        groups << "library_id\thaplotype_group\tdonors\tdonor_count\tcommon_callable_sites"
               << "\tsignature_hash\tresolution_status\n";
        pairwise << "library_id\tdonor1\tdonor2\tjoint_callable_sites\tdistinguishing_sites"
                 << "\tidentical_on_retained_panel\tstatus\n";
        for (size_t l = 0; l < libraries.size(); ++l) {
            const LibraryDef& lib = libraries[l];
            std::vector<size_t> common_sites;
            for (size_t c = 0; c < candidates.size(); ++c) {
                bool all_called = true;
                for (size_t d = 0; d < lib.donor_indices.size(); ++d) {
                    if (candidates[c].evidence[lib.donor_indices[d]].state < 0) {
                        all_called = false;
                        break;
                    }
                }
                if (all_called) common_sites.push_back(c);
            }
            std::map<std::string, std::vector<std::string> > signature_groups;
            for (size_t d = 0; d < lib.donors.size(); ++d) {
                std::string signature;
                signature.reserve(common_sites.size());
                for (size_t k = 0; k < common_sites.size(); ++k) {
                    const int state = candidates[common_sites[k]].evidence[lib.donor_indices[d]].state;
                    signature.push_back(state == 1 ? '1' : '0');
                }
                signature_groups[signature].push_back(lib.donors[d]);
            }
            int group_number = 0;
            for (std::map<std::string, std::vector<std::string> >::const_iterator it =
                     signature_groups.begin(); it != signature_groups.end(); ++it) {
                ++group_number;
                std::ostringstream group_id;
                group_id << "MTG" << std::setw(3) << std::setfill('0') << group_number;
                std::string status;
                if (common_sites.empty()) status = "NO_COMMON_CALLABLE_SITES";
                else if (it->second.size() > 1) status = "COLLAPSED_SHARED_HAPLOTYPE";
                else status = "UNIQUE_ON_RETAINED_PANEL";
                groups << lib.id << '\t' << group_id.str() << '\t' << join(it->second, ",") << '\t'
                       << it->second.size() << '\t' << common_sites.size() << '\t'
                       << hex64(fnv1a64(it->first)) << '\t' << status << '\n';
            }
            for (size_t i = 0; i < lib.donors.size(); ++i) {
                for (size_t j = i + 1; j < lib.donors.size(); ++j) {
                    int callable = 0;
                    int differences = 0;
                    for (size_t c = 0; c < candidates.size(); ++c) {
                        const int a = candidates[c].evidence[lib.donor_indices[i]].state;
                        const int b = candidates[c].evidence[lib.donor_indices[j]].state;
                        if (a < 0 || b < 0) continue;
                        ++callable;
                        if (a != b) ++differences;
                    }
                    const bool identical = callable > 0 && differences == 0;
                    const char* status = callable == 0 ? "NO_JOINT_CALLABLE_SITES" :
                        (identical ? "COLLAPSED_SHARED_HAPLOTYPE" : "DISTINGUISHABLE");
                    pairwise << lib.id << '\t' << lib.donors[i] << '\t' << lib.donors[j] << '\t'
                             << callable << '\t' << differences << '\t' << (identical ? 1 : 0)
                             << '\t' << status << '\n';
                }
            }
        }
    }

    if (config.require_pair_targets &&
        (summary.pairs_below_ratio_target > 0 || summary.pairs_without_ambient_anchor > 0)) {
        throw std::runtime_error(
            "At least one library/pair failed mitochondrial ratio or ambient-anchor targets; see " +
            pair_audit_path);
    }

    bcf_hdr_t* out_header = bcf_hdr_dup(header);
    if (!out_header) throw std::runtime_error("Could not duplicate mt panel header");
    bcf_hdr_append(out_header,
        "##INFO=<ID=MT_LIBRARY_COUNT,Number=1,Type=Integer,Description=\"Libraries for which this site is useful\">");
    bcf_hdr_append(out_header,
        "##INFO=<ID=MT_PAIR_USE_COUNT,Number=1,Type=Integer,Description=\"Library-parent-pair uses for this site\">");
    bcf_hdr_append(out_header,
        "##INFO=<ID=MT_RATIO_USE_COUNT,Number=1,Type=Integer,Description=\"Parent-discriminating library-pair uses\">");
    bcf_hdr_append(out_header,
        "##INFO=<ID=MT_AMBIENT_USE_COUNT,Number=1,Type=Integer,Description=\"Ambient-only library-pair uses\">");
    bcf_hdr_append(out_header,
        "##INFO=<ID=MT_COV_SCORE,Number=1,Type=Float,Description=\"Mitochondrial expression coverage score\">");
    bcf_hdr_append(out_header,
        "##INFO=<ID=MT_DEPTH_FILTER,Number=1,Type=Integer,Description=\"Minimum parental allele depth required by selector\">");
    bcf_hdr_append(out_header,
        "##INFO=<ID=MT_HOMOPLASMY_AF_FILTER,Number=1,Type=Float,Description=\"Near-fixed parental allele-fraction threshold\">");
    bcf_hdr_append(out_header,
        "##mt_panel_selection=depth_based_homoplasmy;library_specific_three_site_classes;metadata_expected_donors_only");
    if (bcf_hdr_sync(out_header) < 0) {
        bcf_hdr_destroy(out_header);
        throw std::runtime_error("Could not sync mt panel output header");
    }

    htsFile* output = hts_open(config.output.c_str(), output_mode(config.output).c_str());
    if (!output) {
        bcf_hdr_destroy(out_header);
        throw std::runtime_error("Could not open mt panel output: " + config.output);
    }
    hts_set_threads(output, config.threads);
    if (bcf_hdr_write(output, out_header) < 0) {
        hts_close(output);
        bcf_hdr_destroy(out_header);
        throw std::runtime_error("Could not write mt panel header");
    }
    std::ofstream bed(sites_bed_path.c_str());
    if (!bed) {
        hts_close(output);
        bcf_hdr_destroy(out_header);
        throw std::runtime_error("Could not write mt site BED: " + sites_bed_path);
    }
    for (size_t c = 0; c < candidates.size(); ++c) {
        const Candidate& candidate = candidates[c];
        std::set<int> library_ids;
        int ratio_uses = 0;
        int ambient_uses = 0;
        for (size_t u = 0; u < candidate.uses.size(); ++u) {
            library_ids.insert(candidate.uses[u].library_index);
            if (candidate.uses[u].site_class == CLASS_AMBIENT_ONLY) ++ambient_uses;
            else ++ratio_uses;
        }
        bcf1_t* record = bcf_dup(candidate.record);
        if (!record) {
            hts_close(output);
            bcf_hdr_destroy(out_header);
            throw std::runtime_error("Could not duplicate mt panel record");
        }
        int32_t library_count = static_cast<int32_t>(library_ids.size());
        int32_t pair_use_count = static_cast<int32_t>(candidate.uses.size());
        int32_t ratio_use_count = ratio_uses;
        int32_t ambient_use_count = ambient_uses;
        int32_t depth_filter = config.min_depth;
        float cov = static_cast<float>(candidate.coverage);
        float hom_af = static_cast<float>(config.homoplasmy_af);
        bcf_update_info_int32(out_header, record, "MT_LIBRARY_COUNT", &library_count, 1);
        bcf_update_info_int32(out_header, record, "MT_PAIR_USE_COUNT", &pair_use_count, 1);
        bcf_update_info_int32(out_header, record, "MT_RATIO_USE_COUNT", &ratio_use_count, 1);
        bcf_update_info_int32(out_header, record, "MT_AMBIENT_USE_COUNT", &ambient_use_count, 1);
        bcf_update_info_float(out_header, record, "MT_COV_SCORE", &cov, 1);
        bcf_update_info_int32(out_header, record, "MT_DEPTH_FILTER", &depth_filter, 1);
        bcf_update_info_float(out_header, record, "MT_HOMOPLASMY_AF_FILTER", &hom_af, 1);
        if (bcf_write(output, out_header, record) < 0) {
            bcf_destroy(record);
            hts_close(output);
            bcf_hdr_destroy(out_header);
            throw std::runtime_error("Failed while writing mt panel record");
        }
        bcf_destroy(record);
        bed << config.mito_chrom << '\t' << candidate.pos << '\t' << candidate.pos + 1 << '\n';
    }
    hts_close(output);
    bcf_hdr_destroy(out_header);
    if (bcf_index_build(config.output.c_str(), 14) < 0) {
        std::fprintf(stderr, "WARNING: could not build CSI index for %s\n", config.output.c_str());
    }

    std::fprintf(stderr, "\nMitochondrial panel selection complete\n");
    std::fprintf(stderr, "  metadata libraries: %d\n", summary.libraries);
    std::fprintf(stderr, "  expected donors: %d of %d VCF samples\n", summary.expected_donors, n_samples);
    std::fprintf(stderr, "  expected heterotypic pairs: %d\n", summary.expected_pairs);
    std::fprintf(stderr, "  chrM records supplied: %lld\n", static_cast<long long>(summary.records_scanned));
    std::fprintf(stderr, "  selected union sites: %lld\n", static_cast<long long>(summary.selected_sites));
    std::fprintf(stderr,
        "  rejected: blacklist=%lld mt_mask=%lld non_snv=%lld quality=%lld coverage=%lld no_library_use=%lld\n",
        static_cast<long long>(summary.rejected_blacklist),
        static_cast<long long>(summary.rejected_mt_mask),
        static_cast<long long>(summary.rejected_non_snv),
        static_cast<long long>(summary.rejected_quality),
        static_cast<long long>(summary.rejected_coverage),
        static_cast<long long>(summary.rejected_no_use));
    std::fprintf(stderr, "  pair target failures: ratio=%d ambient_anchor=%d\n",
        summary.pairs_below_ratio_target, summary.pairs_without_ambient_anchor);
    std::fprintf(stderr, "  output: %s\n", config.output.c_str());
    std::fprintf(stderr, "  site manifest: %s\n", site_manifest_path.c_str());
    std::fprintf(stderr, "  haplotype groups: %s\n", haplotype_groups_path.c_str());
    return summary;
}

}  // namespace mtpanel

#include <chrono>
#include <omp.h>
#include <mutex>
#include <atomic>

using std::cout;
using std::endl;
using namespace std;

// ============================================================================
// PACKED GENOTYPE ACCESS
// ============================================================================
// Each sample has 2 alleles; each allele takes 2 bits.
// Byte b in the packed vector holds 4 allele values: positions 0..3 within byte
// b correspond to global allele index 4b..4b+3.
//
// For sample s, ploidy 2:
//   allele 0 = global index 2s
//   allele 1 = global index 2s + 1
//   byte index = (2s + a) / 4
//   bit shift  = ((2s + a) % 4) * 2
//   value      = (packed[byte] >> shift) & 0x3
//
// Encoding per allele slot (2 bits):
//   00 = ref (allele index 0)
//   01 = alt (allele index 1)
//   10 = missing (bcf_gt_is_missing was true)
//   11 = reserved (unused)

inline uint8_t gt_packed_get(const vector<uint8_t>& packed, int sample, int allele) {
    int gi = 2 * sample + allele;
    return (packed[gi / 4] >> ((gi % 4) * 2)) & 0x3;
}
inline bool gt_packed_is_missing(const vector<uint8_t>& packed, int sample, int allele) {
    return gt_packed_get(packed, sample, allele) == 0x2;
}
inline int gt_packed_allele(const vector<uint8_t>& packed, int sample, int allele) {
    // Returns 0 (ref), 1 (alt), or -1 (missing). For biallelic only.
    uint8_t v = gt_packed_get(packed, sample, allele);
    if (v == 0x2) return -1;
    return (int)v;
}
inline void gt_packed_set(vector<uint8_t>& packed, int sample, int allele, uint8_t val) {
    int gi = 2 * sample + allele;
    int byte_idx = gi / 4;
    int shift = (gi % 4) * 2;
    packed[byte_idx] = (packed[byte_idx] & ~(0x3 << shift)) | ((val & 0x3) << shift);
}

// Pack htslib int32_t genotypes into 2-bit-per-allele form.
// For multi-allelic sites (alt index > 1), stores 0x1 (alt) for any non-ref
// non-missing allele. This is acceptable because multi-allelic sites are
// filtered before genotypes are consumed for selection.
inline void pack_genotypes(const int32_t* raw_gts, int n_samples, vector<uint8_t>& out) {
    int n_alleles_total = n_samples * 2;
    int n_bytes = (n_alleles_total + 3) / 4;
    out.assign(n_bytes, 0);
    for (int s = 0; s < n_samples; s++) {
        for (int a = 0; a < 2; a++) {
            int32_t gt_raw = raw_gts[s * 2 + a];
            uint8_t val;
            if (bcf_gt_is_missing(gt_raw)) {
                val = 0x2; // missing
            } else {
                int allele_idx = bcf_gt_allele(gt_raw);
                if (allele_idx == 0) {
                    val = 0x0; // ref
                } else {
                    val = 0x1; // alt (any non-ref)
                }
            }
            gt_packed_set(out, s, a, val);
        }
    }
}

// Unpack 2-bit genotypes back to htslib int32_t format for BCF output.
// Produces phased-looking output matching bcf_gt_phased/bcf_gt_unphased conventions.
inline void unpack_genotypes(const vector<uint8_t>& packed, int n_samples, vector<int32_t>& out) {
    out.resize(n_samples * 2);
    for (int s = 0; s < n_samples; s++) {
        for (int a = 0; a < 2; a++) {
            uint8_t v = gt_packed_get(packed, s, a);
            if (v == 0x2) {
                out[s * 2 + a] = bcf_gt_missing;
            } else {
                out[s * 2 + a] = bcf_gt_unphased((int)v);
            }
        }
    }
}

// ============================================================================
// HET STATS HELPER (C1)
// ============================================================================
struct HetStats {
    int n_het;
    int n_called;
    int n_missing;
    float het_freq;
    float missing_freq;
};

HetStats compute_het_stats_packed(const vector<uint8_t>& packed, int n_samples) {
    HetStats hs = {0, 0, 0, 0.0f, 0.0f};
    for (int i = 0; i < n_samples; ++i) {
        int a1 = gt_packed_allele(packed, i, 0);
        int a2 = gt_packed_allele(packed, i, 1);
        if (a1 < 0 || a2 < 0) {
            hs.n_missing++;
        } else {
            hs.n_called++;
            if ((a1 == 0 && a2 != 0) || (a1 != 0 && a2 == 0)) hs.n_het++;
        }
    }
    hs.het_freq = (hs.n_called > 0) ? (float)hs.n_het / hs.n_called : 0.0f;
    hs.missing_freq = (float)hs.n_missing / n_samples;
    return hs;
}

// ============================================================================
// DOSAGE FLIP HELPER (A3)
// ============================================================================
// Flips REF/ALT polarity under dosage-aware encoding.
// 0/0 (0,0) -> 1/1 (1,1)
// 0/1 (1,0) -> 0/1 (1,0)  [het unchanged]
// 1/1 (1,1) -> 0/0 (0,0)
// Missing stays (0,0) since it's guarded by present[i].
void dosage_flip(const bitset<NBITS>& src, const bitset<NBITS>& present,
                 int num_samples, bitset<NBITS>& dst) {
    dst.reset();
    for (int i = 0; i < num_samples; i++) {
        if (!present[i]) continue;
        bool alt_present = src[2*i];
        bool homo_alt   = src[2*i + 1];
        bool new_alt_present = !homo_alt;
        bool new_homo_alt   = !alt_present;
        if (new_alt_present) dst.set(2*i);
        if (new_homo_alt)    dst.set(2*i + 1);
    }
}

// Compute total_alt from a dosage-aware bitset (sum of bit-pair values)
inline int compute_total_alt(const bitset<NBITS>& bs, int num_samples) {
    int total = 0;
    for (int i = 0; i < num_samples; i++) {
        if (bs[2*i])     total++;
        if (bs[2*i + 1]) total++;
    }
    return total;
}

// ============================================================================
// BJPs ADDITIONS: GTF and Coverage support structures
// ============================================================================

struct Interval {
    int64_t start;
    int64_t end;
    bool operator<(const Interval& other) const {
        return start < other.start ||
               (start == other.start && end < other.end);
    }
};

struct BlacklistExclusions {
    set<string> whole_scaffolds;
    map<string, vector<Interval>> intervals;

    bool empty() const {
        return whole_scaffolds.empty() && intervals.empty();
    }
};

// hts_open/hts_getline transparently handle both ordinary text and BGZF/gzip
// input.  The production annotation is .gtf.gz, so using std::ifstream here
// would parse compressed bytes as text and silently produce an empty map.
template <typename LineConsumer>
void for_each_annotation_line(const string& filename,
                              const char* description,
                              LineConsumer consume) {
    htsFile* input = hts_open(filename.c_str(), "r");
    if (input == NULL) {
        fprintf(stderr, "ERROR: Could not open %s: %s\n",
            description, filename.c_str());
        exit(1);
    }
    kstring_t buffer = {0, 0, NULL};
    int read_status = -1;
    long line_number = 0;
    while ((read_status = hts_getline(input, KS_SEP_LINE, &buffer)) >= 0) {
        line_number++;
        consume(string(buffer.s == NULL ? "" : buffer.s, buffer.l),
                line_number);
    }
    free(buffer.s);
    const int close_status = hts_close(input);
    if (read_status < -1) {
        fprintf(stderr, "ERROR: Failed while reading %s at line %ld: %s\n",
            description, line_number + 1, filename.c_str());
        exit(1);
    }
    if (close_status != 0) {
        fprintf(stderr, "ERROR: Failed to close %s cleanly: %s\n",
            description, filename.c_str());
        exit(1);
    }
}

void load_gtf(const string& gtf_file, map<string, vector<Interval>>& annotations) {
    long annotation_count = 0;
    for_each_annotation_line(gtf_file, "GTF annotation",
        [&](const string& line, long /*line_number*/) {
        if (line.empty() || line[0] == '#') return;
        stringstream ss(line);
        string seqname, source, feature, start_s, end_s;
        ss >> seqname >> source >> feature >> start_s >> end_s;
        
        if (feature == "gene" || feature == "exon") {
            try {
                Interval iv = {stoi(start_s), stoi(end_s)};
                annotations[seqname].push_back(iv);
                annotation_count++;
            } catch (...) { return; }
        }
    });
    if (annotation_count == 0 || annotations.empty()) {
        fprintf(stderr,
            "ERROR: GTF contains no parseable gene/exon annotations: %s\n",
            gtf_file.c_str());
        exit(1);
    }
    for (auto& pair : annotations) {
        sort(pair.second.begin(), pair.second.end());

        // Gene and exon records are commonly nested.  The interval lookup below
        // checks only the lower-bound record and its immediate predecessor, so
        // retaining nested records can hide an earlier, longer gene interval.
        // Collapse the union of overlapping 1-based inclusive gene-body
        // intervals so that lookup remains correct without changing the intended
        // gene-body-versus-intergenic scoring semantics.
        vector<Interval> merged;
        merged.reserve(pair.second.size());
        for (const Interval& iv : pair.second) {
            if (merged.empty() || iv.start > merged.back().end) {
                merged.push_back(iv);
            } else if (iv.end > merged.back().end) {
                merged.back().end = iv.end;
            }
        }
        pair.second.swap(merged);
    }
    fprintf(stderr,
        "Loaded %ld gene/exon annotations for %lu chromosomes from %s.\n",
        annotation_count, annotations.size(), gtf_file.c_str());
}

// Build the ATAC default regulatory prior from strand-aware core promoters.
// Coordinates are stored as 0-based half-open intervals, matching BED and the
// 0-based BCF positions passed to is_blacklisted().  This is a ranking prior,
// never an ATAC eligibility filter.
void load_gtf_promoters(const string& gtf_file,
                        map<string, vector<Interval>>& promoters,
                        int upstream_bp = 2000) {
    long count = 0;
    for_each_annotation_line(gtf_file, "GTF promoter annotation",
        [&](const string& line, long /*line_number*/) {
        if (line.empty() || line[0] == '#') return;
        vector<string> fields;
        string field;
        stringstream ss(line);
        while (getline(ss, field, '\t')) fields.push_back(field);
        if (fields.size() < 7 || fields[2] != "gene") return;
        try {
            const int start1 = stoi(fields[3]);
            const int end1 = stoi(fields[4]);
            const char strand = fields[6].empty() ? '+' : fields[6][0];
            Interval iv;
            if (strand == '-') {
                const int tss0 = end1 - 1;
                iv.start = tss0 + 1;
                iv.end = tss0 + 1 + upstream_bp;
            } else {
                const int tss0 = start1 - 1;
                iv.start = max(0, tss0 - upstream_bp);
                iv.end = tss0;
            }
            if (iv.end > iv.start) {
                promoters[fields[0]].push_back(iv);
                count++;
            }
        } catch (...) { return; }
    });
    if (count == 0 || promoters.empty()) {
        fprintf(stderr,
            "ERROR: GTF contains no parseable gene records for promoter "
            "derivation: %s\n", gtf_file.c_str());
        exit(1);
    }
    for (auto& pair : promoters) {
        sort(pair.second.begin(), pair.second.end());
        vector<Interval> merged;
        merged.reserve(pair.second.size());
        for (const Interval& iv : pair.second) {
            if (merged.empty() || iv.start > merged.back().end) {
                merged.push_back(iv);
            } else if (iv.end > merged.back().end) {
                merged.back().end = iv.end;
            }
        }
        pair.second.swap(merged);
    }
    fprintf(stderr,
        "Derived %ld strand-aware %d-bp upstream promoter intervals for %lu "
        "chromosomes from %s.\n",
        count, upstream_bp, promoters.size(), gtf_file.c_str());
}

// B11: Load BED file (0-based half-open intervals). This permissive loader is
// retained for NUMT and ATAC-regulatory inputs; their behavior is independent
// of the user blacklist.
void load_bed(const string& bed_file, map<string, vector<Interval>>& blacklist,
              const char* label = "blacklist") {
    ifstream file(bed_file);
    if (!file.is_open()) {
        fprintf(stderr, "ERROR: Could not open %s BED file: %s\n", label, bed_file.c_str());
        exit(1);
    }
    string line;
    long count = 0;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        string seqname, start_s, end_s;
        ss >> seqname >> start_s >> end_s;
        try {
            Interval iv = {stoi(start_s), stoi(end_s)};
            blacklist[seqname].push_back(iv);
            count++;
        } catch (...) { continue; }
    }
    for (auto& pair : blacklist) {
        sort(pair.second.begin(), pair.second.end());
        // Merge overlapping/adjacent records so interval lookup remains correct
        // even when a BED contains nested annotations.
        vector<Interval> merged;
        merged.reserve(pair.second.size());
        for (const Interval& iv : pair.second) {
            if (merged.empty() || iv.start > merged.back().end) {
                merged.push_back(iv);
            } else if (iv.end > merged.back().end) {
                merged.back().end = iv.end;
            }
        }
        pair.second.swap(merged);
    }
    fprintf(stderr, "Loaded %ld %s intervals for %lu chromosomes.\n",
        count, label, blacklist.size());
}

bool parse_blacklist_coordinate(const string& text, int64_t& value) {
    if (text.empty()) return false;
    errno = 0;
    char* end = NULL;
    const long long parsed = strtoll(text.c_str(), &end, 10);
    if (errno == ERANGE || end == text.c_str() || *end != '\0') return false;
    value = static_cast<int64_t>(parsed);
    return true;
}

BlacklistExclusions load_blacklist(const string& blacklist_file) {
    ifstream file(blacklist_file.c_str());
    if (!file.is_open()) {
        fprintf(stderr, "ERROR: Could not open blacklist file: %s\n",
            blacklist_file.c_str());
        exit(1);
    }

    BlacklistExclusions result;
    string line;
    long line_number = 0;
    long interval_records = 0;
    while (getline(file, line)) {
        line_number++;
        const string::size_type first = line.find_first_not_of(" \t\r\n");
        if (first == string::npos || line[first] == '#') continue;

        vector<string> fields;
        string field;
        stringstream ss(line);
        while (ss >> field) fields.push_back(field);
        if (fields.size() == 1) {
            result.whole_scaffolds.insert(fields[0]);
            continue;
        }
        if (fields.size() == 2) {
            fprintf(stderr,
                "ERROR: Malformed blacklist line %ld in %s: two-column "
                "records are not allowed\n",
                line_number, blacklist_file.c_str());
            exit(1);
        }
        if (fields.size() < 3) continue;

        int64_t start = 0;
        int64_t end = 0;
        if (!parse_blacklist_coordinate(fields[1], start) ||
            !parse_blacklist_coordinate(fields[2], end) ||
            start < 0 || end <= start) {
            fprintf(stderr,
                "ERROR: Invalid 0-based half-open blacklist coordinates at "
                "line %ld in %s: %s %s %s\n",
                line_number, blacklist_file.c_str(), fields[0].c_str(),
                fields[1].c_str(), fields[2].c_str());
            exit(1);
        }
        result.intervals[fields[0]].push_back(Interval{start, end});
        interval_records++;
    }
    if (file.bad()) {
        fprintf(stderr, "ERROR: Failed while reading blacklist file: %s\n",
            blacklist_file.c_str());
        exit(1);
    }

    // A whole-scaffold record supersedes every interval for that scaffold.
    for (const string& scaffold : result.whole_scaffolds) {
        result.intervals.erase(scaffold);
    }

    size_t merged_intervals = 0;
    for (auto& pair : result.intervals) {
        sort(pair.second.begin(), pair.second.end());
        vector<Interval> merged;
        merged.reserve(pair.second.size());
        for (const Interval& iv : pair.second) {
            if (merged.empty() || iv.start > merged.back().end) {
                merged.push_back(iv);
            } else if (iv.end > merged.back().end) {
                merged.back().end = iv.end;
            }
        }
        pair.second.swap(merged);
        merged_intervals += pair.second.size();
    }

    fprintf(stderr, "Blacklist whole scaffolds loaded: %zu\n",
        result.whole_scaffolds.size());
    fprintf(stderr, "Blacklist interval records loaded: %ld\n",
        interval_records);
    fprintf(stderr, "Blacklist merged intervals: %zu\n", merged_intervals);
    return result;
}

void validate_blacklist_scaffolds(
    const BlacklistExclusions& blacklist,
    const map<string, int>& vcf_scaffolds,
    const string& blacklist_file) {
    set<string> unknown;
    for (const string& scaffold : blacklist.whole_scaffolds) {
        if (vcf_scaffolds.count(scaffold) == 0) unknown.insert(scaffold);
    }
    for (const auto& pair : blacklist.intervals) {
        if (vcf_scaffolds.count(pair.first) == 0) unknown.insert(pair.first);
    }
    if (!unknown.empty()) {
        string names;
        for (const string& scaffold : unknown) {
            if (!names.empty()) names += ",";
            names += scaffold;
        }
        fprintf(stderr,
            "ERROR: Blacklist scaffold(s) absent from input VCF header in %s: %s\n",
            blacklist_file.c_str(), names.c_str());
        exit(1);
    }
}

bool is_blacklisted(const string& chrom, int pos,
                    const map<string, vector<Interval>>& bl) {
    if (bl.count(chrom) == 0) return false;
    const auto& ivs = bl.at(chrom);
    auto it = lower_bound(ivs.begin(), ivs.end(), Interval{pos, pos});
    if (it != ivs.end() && it->start <= pos && it->end > pos) return true;
    if (it != ivs.begin()) {
        auto prev = std::prev(it);
        if (prev->start <= pos && prev->end > pos) return true;
    }
    return false;
}

bool is_blacklisted(const string& chrom, int pos,
                    const BlacklistExclusions& blacklist) {
    if (blacklist.whole_scaffolds.count(chrom) != 0) return true;
    return is_blacklisted(chrom, pos, blacklist.intervals);
}

// Nuclear-panel exclusion bits. The annotated-all output retains excluded
// records and reports the reason, while every selected nuclear panel omits them.
enum NuclearExclusion : uint8_t {
    NUCLEAR_EXCLUDE_NONE      = 0,
    NUCLEAR_EXCLUDE_MITO      = 1u << 0,
    NUCLEAR_EXCLUDE_NUMT      = 1u << 1,
    NUCLEAR_EXCLUDE_BLACKLIST = 1u << 2
};

uint8_t nuclear_exclusion_mask(const string& chrom, int pos,
                               const string& mito_chrom,
                               bool include_mito,
                               const map<string, vector<Interval>>& numts,
                               bool include_numts,
                               const BlacklistExclusions& blacklist) {
    uint8_t mask = NUCLEAR_EXCLUDE_NONE;
    if (!include_mito && chrom == mito_chrom) {
        mask |= NUCLEAR_EXCLUDE_MITO;
    }
    if (!include_numts && !numts.empty() && is_blacklisted(chrom, pos, numts)) {
        mask |= NUCLEAR_EXCLUDE_NUMT;
    }
    if (!blacklist.empty() && is_blacklisted(chrom, pos, blacklist)) {
        mask |= NUCLEAR_EXCLUDE_BLACKLIST;
    }
    return mask;
}

string nuclear_exclusion_reason(uint8_t mask) {
    if (mask == NUCLEAR_EXCLUDE_NONE) return "PASS";
    string reason;
    if (mask & NUCLEAR_EXCLUDE_MITO) reason += "MITO";
    if (mask & NUCLEAR_EXCLUDE_NUMT) {
        if (!reason.empty()) reason += "+";
        reason += "NUMT";
    }
    if (mask & NUCLEAR_EXCLUDE_BLACKLIST) {
        if (!reason.empty()) reason += "+";
        reason += "BLACKLIST";
    }
    return reason;
}

float get_annot_score(const string& chrom, int pos, const map<string, vector<Interval>>& annotations) {
    if (annotations.count(chrom) == 0) return 0.1;

    const vector<Interval>& ivs = annotations.at(chrom);
    auto it = lower_bound(ivs.begin(), ivs.end(), Interval{pos, pos});

    if (it != ivs.end() && it->start <= pos && it->end >= pos) return 1.0;
    if (it != ivs.begin()) {
        auto prev = std::prev(it);
        if (prev->start <= pos && prev->end >= pos) return 1.0;
    }

    return 0.1;
}

// Coverage interval with score
struct CoverageInterval {
    int start;
    int end;
    float score;
    bool operator<(const CoverageInterval& other) const {
        return start < other.start;
    }
};

// Load coverage bedgraph into memory
map<string, vector<CoverageInterval>> load_coverage(const string& filename) {
    map<string, vector<CoverageInterval>> cov_map;

    htsFile* fp = hts_open(filename.c_str(), "r");
    if (!fp) return cov_map;
    
    kstring_t str = {0, 0, 0};
    long count = 0;
    
    while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
        if (str.l == 0 || str.s[0] == '#') continue;
        
        char* chrom = str.s;
        char* p = str.s;
        int col = 0;
        int start = 0, end = 0;
        float score = 0;
        
        while (*p) {
            if (*p == '\t') {
                *p = '\0';
                col++;
                if (col == 1) start = atoi(p + 1);
                else if (col == 2) end = atoi(p + 1);
                else if (col == 3) {
                    score = atof(p + 1);
                    break;
                }
            }
            p++;
        }
        
        if (col >= 3) {
            cov_map[string(chrom)].push_back({start, end, score});
            count++;
            if (count % 10000000 == 0) {
                fprintf(stderr, "  Loaded %ld coverage intervals...\r", count);
            }
        }
    }
    
    free(str.s);
    hts_close(fp);
    
    // Sort each chromosome's intervals
    for (auto& kv : cov_map) {
        sort(kv.second.begin(), kv.second.end());
    }
    
    fprintf(stderr, "  Loaded %ld coverage intervals from %zu chromosomes\n", count, cov_map.size());
    return cov_map;
}

// Fast in-memory coverage lookup
float get_coverage_score_fast(const string& chrom, int pos, const map<string, vector<CoverageInterval>>& coverage) {
    auto it = coverage.find(chrom);
    if (it == coverage.end()) return 0.0;
    
    const vector<CoverageInterval>& ivs = it->second;
    auto cit = lower_bound(ivs.begin(), ivs.end(), CoverageInterval{pos, 0, 0});
    
    // Check current interval
    if (cit != ivs.end() && cit->start <= pos && cit->end > pos) {
        return cit->score;
    }
    // Check previous interval
    if (cit != ivs.begin()) {
        auto prev = std::prev(cit);
        if (prev->start <= pos && prev->end > pos) {
            return prev->score;
        }
    }
    
    return 0.0;
}

float get_coverage_score(const string& chrom, int pos, tbx_t* tbx, htsFile* fp) {
    if (!tbx || !fp) return 0.0;
    
    int tid = tbx_name2id(tbx, chrom.c_str());
    if (tid < 0) return 0.0;

    hts_itr_t* itr = tbx_itr_queryi(tbx, tid, pos, pos + 1);
    if (!itr) return 0.0;

    kstring_t str = {0, 0, 0};
    float cov = 0.0;
    
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
        char* p = str.s;
        int col = 0;
        while (*p) {
            if (col == 3) {
                cov = atof(p);
                break;
            }
            if (*p == '\t') col++;
            p++;
        }
    }
    
    tbx_itr_destroy(itr);
    free(str.s);
    return cov;
}

// ============================================================================
// PANEL METADATA AND POOL COMBINATIONS STRUCTURES
// ============================================================================

// PanelMetadata struct and load_panel_metadata() are in common.h / common.cpp

// Pool combinations: library -> list of cell identities
struct PoolCombinations {
    // Each identity is either a single indiv or a pair (heterotypic fusion)
    struct CellIdentity {
        vector<string> components; // 1 element = singlet/homotypic, 2 = heterotypic
    };
    
    map<int, vector<CellIdentity>> lib_identities; // lib_num -> identities
};

static string join_identity_components(const vector<string>& components) {
    string label;
    for (size_t i = 0; i < components.size(); ++i) {
        if (i) label += "+";
        label += components[i];
    }
    return label;
}

// Compact donor pairs that coexist in at least one supplied library. These are
// retained solely for post-selection diagnostics; they never create ATAC
// selection objectives, floors, eligibility gates, or rescue targets.
struct PoolDonorPair {
    int first_index;
    int second_index;
    string first_name;
    string second_name;
    vector<int> library_ids;
};

vector<PoolDonorPair> build_pool_donor_pairs(
    const PoolCombinations& pool_combos,
    const map<string, int>& sample_name_to_idx,
    const vector<string>& vcf_samples,
    vector<int>& roster_sample_indices) {

    map<pair<int, int>, set<int>> pair_libraries;
    set<int> roster_samples;
    for (const auto& lib_entry : pool_combos.lib_identities) {
        set<int> lib_samples;
        for (const auto& identity : lib_entry.second) {
            for (const string& component : identity.components) {
                auto found = sample_name_to_idx.find(component);
                if (found == sample_name_to_idx.end()) continue;
                lib_samples.insert(found->second);
                roster_samples.insert(found->second);
            }
        }
        vector<int> donors(lib_samples.begin(), lib_samples.end());
        for (size_t i = 0; i < donors.size(); ++i) {
            for (size_t j = i + 1; j < donors.size(); ++j) {
                pair_libraries[make_pair(donors[i], donors[j])].insert(
                    lib_entry.first);
            }
        }
    }

    roster_sample_indices.assign(roster_samples.begin(), roster_samples.end());
    vector<PoolDonorPair> result;
    result.reserve(pair_libraries.size());
    for (const auto& entry : pair_libraries) {
        PoolDonorPair pair;
        pair.first_index = entry.first.first;
        pair.second_index = entry.first.second;
        pair.first_name = vcf_samples[pair.first_index];
        pair.second_name = vcf_samples[pair.second_index];
        pair.library_ids.assign(entry.second.begin(), entry.second.end());
        result.push_back(std::move(pair));
    }
    return result;
}

PoolCombinations load_pool_combinations(const string& filename,
                                         const vector<string>& vcf_samples) {
    PoolCombinations pc;
    
    set<string> sample_set(vcf_samples.begin(), vcf_samples.end());
    
    ifstream f(filename);
    if (!f.is_open()) {
        fprintf(stderr, "ERROR: Could not open pool combinations: %s\n", filename.c_str());
        exit(1);
    }
    
    string line;
    getline(f, line); // skip header
    auto trim_field = [](string value) {
        const size_t first = value.find_first_not_of(" \t\r\n");
        if (first == string::npos) return string();
        const size_t last = value.find_last_not_of(" \t\r\n");
        return value.substr(first, last - first + 1);
    };
    long line_number = 1;
    
    while (getline(f, line)) {
        line_number++;
        if (trim_field(line).empty()) continue;
        size_t tab = line.find('\t');
        if (tab == string::npos) {
            fprintf(stderr,
                "ERROR: Pool-combinations line %ld has no tab delimiter\n",
                line_number);
            exit(1);
        }
        
        const string library_field = trim_field(line.substr(0, tab));
        char* library_end = NULL;
        errno = 0;
        const long parsed_library = strtol(
            library_field.c_str(), &library_end, 10);
        if (errno != 0 || library_end == library_field.c_str() ||
            *library_end != '\0' || parsed_library < 0 ||
            parsed_library > numeric_limits<int>::max()) {
            fprintf(stderr,
                "ERROR: Invalid library ID '%s' on pool-combinations line %ld\n",
                library_field.c_str(), line_number);
            exit(1);
        }
        const int lib_num = (int)parsed_library;
        if (pc.lib_identities.count(lib_num)) {
            fprintf(stderr,
                "ERROR: Duplicate library ID %d in pool combinations "
                "(line %ld)\n", lib_num, line_number);
            exit(1);
        }
        const string identities_str = trim_field(line.substr(tab + 1));
        if (identities_str.empty()) {
            fprintf(stderr,
                "ERROR: Library %d has an empty identity roster\n", lib_num);
            exit(1);
        }
        if (identities_str.front() == ',' || identities_str.back() == ',') {
            fprintf(stderr,
                "ERROR: Library %d contains an empty identity token\n",
                lib_num);
            exit(1);
        }
        
        vector<PoolCombinations::CellIdentity> identities;
        set<string> seen_identity_labels;
        stringstream ss(identities_str);
        string token;
        while (getline(ss, token, ',')) {
            token = trim_field(token);
            if (token.empty()) {
                fprintf(stderr,
                    "ERROR: Library %d contains an empty identity token\n",
                    lib_num);
                exit(1);
            }
            PoolCombinations::CellIdentity ci;
            size_t plus = token.find('+');
            if (plus != string::npos) {
                if (token.find('+', plus + 1) != string::npos) {
                    fprintf(stderr,
                        "ERROR: Identity '%s' in library %d has more than "
                        "two components\n", token.c_str(), lib_num);
                    exit(1);
                }
                ci.components.push_back(trim_field(token.substr(0, plus)));
                ci.components.push_back(trim_field(token.substr(plus + 1)));
            } else {
                ci.components.push_back(trim_field(token));
            }
            for (const string& component : ci.components) {
                if (component.empty()) {
                    fprintf(stderr,
                        "ERROR: Identity '%s' in library %d contains an "
                        "empty component\n", token.c_str(), lib_num);
                    exit(1);
                }
            }
            if (ci.components.size() == 2 &&
                ci.components[1] < ci.components[0]) {
                swap(ci.components[0], ci.components[1]);
            }
            const string canonical_label =
                join_identity_components(ci.components);
            if (!seen_identity_labels.insert(canonical_label).second) {
                continue;
            }
            
            // Validate all components exist in VCF
            for (const auto& comp : ci.components) {
                if (sample_set.count(comp) == 0) {
                    fprintf(stderr, "ERROR: Identity component '%s' in library %d not found in VCF samples\n",
                        comp.c_str(), lib_num);
                    exit(1);
                }
            }
            identities.push_back(ci);
        }
        if (identities.empty()) {
            fprintf(stderr,
                "ERROR: Library %d has no valid unique identities\n", lib_num);
            exit(1);
        }
        
        pc.lib_identities.emplace(lib_num, std::move(identities));
    }
    
    fprintf(stderr, "Pool combinations: %zu libraries\n", pc.lib_identities.size());
    return pc;
}

// ============================================================================
// SPECIES CANDIDATE STRUCTURE
// ============================================================================

struct SpeciesCandidate {
    string chrom;
    int64_t pos;
    int chrom_idx;
    int bin_idx;
    float cov_score;             // coverage used for selection (RNA or ATAC)
    vector<float> pair_score;    // legacy multiplicative score per species pair (recorded only)
    vector<float> pair_discrim;  // PRIMARY ranking score per pair
    float max_pair_score;        // max across all pairs (legacy)
    float max_pair_discrim;      // max pair_discrim across all pairs (primary ranking)
    float selection_score;       // max_pair_score * min(log2(cov+1), cov_cap) (legacy, retained)
    bool atac_regulatory;        // overlaps supplied regulatory BED/default 2-kb promoter prior
    int assigned_pair;           // index into species_pairs for best-fit assignment
    string ref_allele;
    string alt_allele;
    size_t variant_idx;          // index into all_variants for deferred genotype fetch
    int n_samples;
    
    // Per-pair dosage separation (closed-form): |sp_mean_alt[A] - sp_mean_alt[B]| / 2
    // Used as the multi-pair contribution weight in global utility scoring.
    vector<float> pair_dosage_sep;
    
    // Per-species stats for output/audit
    vector<float> sp_mean_alt;            // weighted mean alt dosage per species
    vector<float> sp_major_dosage_frac;   // within-species consistency (weighted)
    vector<float> sp_call_rate;           // call rate per species (weighted)
    vector<int>   sp_n_called;            // unweighted count of called individuals per species
    
    // For sorting by selection_score descending
    bool operator<(const SpeciesCandidate& other) const {
        return selection_score > other.selection_score;
    }
};

// ============================================================================
// HET SITE CANDIDATE STRUCTURE
// ============================================================================

struct HetCandidate {
    string chrom;
    int64_t pos;
    int chrom_idx;          // Index into chroms vector for ordering
    int bin_idx;            // Bin index within chromosome
    float het_score;        // Combined het score
    float het_freq;         // Fraction of called individuals that are het
    float missing_freq;     // Fraction of individuals with missing genotype
    float annot_score;      // Annotation score (if GTF provided)
    float cov_score;        // Coverage score (if coverage provided)
    float demux_score;      // Clade rarity score (computed after Pass 1)
    float atac_cov_score;   // ATAC coverage score (if --atac_cov set)
    float pool_discrim_norm;// Pool discrimination score (if --enable_pool_scoring set)
    float het_sel_score;    // Selection ranking score (for HET_SEL_SCORE INFO field)
    float atac_het_sel_score;// ATAC-specific ranking score used by ATAC-het selection
    bool atac_regulatory;    // overlaps supplied regulatory BED/default 2-kb promoter prior
    bool passes_rna_het_hwe; // legacy pan-panel HWE gate; ATAC selection does not use it
    string ref_allele;
    string alt_allele;
    size_t variant_idx;     // index into all_variants for deferred genotype fetch
    int n_samples;
    
    // For sorting by het_score descending
    bool operator<(const HetCandidate& other) const {
        return het_score > other.het_score;  // Higher score = better
    }
};

// ============================================================================
// STORED VARIANT STRUCTURE - for in-memory VCF processing
// ============================================================================

struct VariantPattern {
    bitset<NBITS> alt;
    bitset<NBITS> alt_flip;
    bitset<NBITS> present;
};

struct StoredVariant {
    int chrom_idx;              // bcf_record->rid
    int64_t pos;                // bcf_record->pos
    int n_alleles;              // bcf_record->n_allele
    float qual;                 // bcf_record->qual (preserved for output)
    string alleles;             // All alleles comma-separated: "A,T" or "A,T,C"
    vector<uint8_t> genotypes_packed;  // Packed 2-bit-per-allele genotypes
    vector<int32_t> genotypes_raw;     // Raw genotypes for non-biallelic (annotate_all only)
    
    // Precomputed by proc_bcf_record logic (only valid for biallelic)
    bool is_biallelic;          // n_alleles == 2
    bool passes_qc;             // What proc_bcf_record would return
    uint8_t nuclear_exclusion;  // NuclearExclusion bit mask (chrM/NUMT/blacklist)
    // Exact genotype patterns are interned once instead of retaining three
    // fixed 2048-bit values in every source record.
    const VariantPattern* pattern;
    int16_t total_alt;          // Total alt allele count across called samples
    int16_t total_alt_flip;     // Total alt allele count for flipped polarity

    // Exact, immutable lookup results reused by every selector/output pass.
    // Coverage and interval lookup functions already return float/bool, so
    // caching them changes neither precision nor selection semantics.
    float rna_cov_score;
    float atac_cov_score;
    float annot_score;
    bool atac_regulatory;
    bool atac_soft_overlap;
    
    // For stratified selection (Option C)
    size_t variant_idx;         // Index into all_variants for back-reference
};

bool operator<(const std::bitset<NBITS>& a, const std::bitset<NBITS>& b);

inline int pattern_dosage(const VariantPattern* pattern, int sample_index) {
    if (pattern == NULL || !pattern->present.test(sample_index)) return -1;
    if (pattern->alt.test(2 * sample_index + 1)) return 2;
    if (pattern->alt.test(2 * sample_index)) return 1;
    return 0;
}

inline float operational_roster_call_rate(
    const VariantPattern* pattern,
    const vector<int>& roster_sample_indices) {
    if (pattern == NULL || roster_sample_indices.empty()) return 0.0f;
    int called = 0;
    for (int sample_index : roster_sample_indices) {
        if (pattern->present.test(sample_index)) called++;
    }
    return (float)called / roster_sample_indices.size();
}

inline bool donor_pair_is_distinguished(
    const VariantPattern* pattern,
    const PoolDonorPair& pair,
    float min_pair_score) {
    const int first = pattern_dosage(pattern, pair.first_index);
    const int second = pattern_dosage(pattern, pair.second_index);
    if (first < 0 || second < 0) return false;
    return fabsf((float)first / 2.0f - (float)second / 2.0f) >=
        min_pair_score;
}

inline pair<bitset<NBITS>, bitset<NBITS>> native_atac_clade_key(
    const StoredVariant& variant) {
    bitset<NBITS> oriented;
    if (variant.total_alt < variant.total_alt_flip ||
        (variant.total_alt == variant.total_alt_flip &&
         variant.pattern->alt < variant.pattern->alt_flip)) {
        oriented = variant.pattern->alt;
    } else {
        oriented = variant.pattern->alt_flip;
    }
    // The present mask is part of the key: once partially called sites are
    // eligible, treating missing as homozygous reference would collapse
    // biologically different patterns and corrupt allocation/capacity counts.
    return make_pair(variant.pattern->present, oriented);
}

inline void append_atac_soft_overlap_headers(
    bcf_hdr_t* header, float multiplier) {
    bcf_hdr_append(header,
        "##INFO=<ID=ATAC_RNA_SOFT_OVERLAP,Number=1,Type=Integer,Description=\"1 when this locus occurs in the supplied selected RNA-demux panel\">");
    bcf_hdr_append(header,
        "##cellbouncer_atac_rna_overlap_policy=soft_penalty");
    bcf_hdr_append(header,
        "##cellbouncer_atac_pool_objective_policy=rna_style_scalar_pool_discrimination_no_identity_pair_floors");
    char multiplier_header[128];
    snprintf(multiplier_header, sizeof(multiplier_header),
        "##cellbouncer_atac_soft_overlap_multiplier=%.6f", multiplier);
    bcf_hdr_append(header, multiplier_header);
}

inline float apply_atac_soft_overlap_penalty(
    float score, float multiplier) {
    if (score >= 0.0f) return score * multiplier;
    if (multiplier > 0.0f) return score / multiplier;
    return -numeric_limits<float>::infinity();
}

// ============================================================================
// PARALLEL CODE
// ============================================================================

bool operator<(const std::bitset<NBITS>& a, const std::bitset<NBITS>& b) {
    for (int i = NBITS-1; i >= 0; --i) {
        if (a[i] != b[i]) return b[i];
    }
    return false;
}

bool parse_nonnegative_int_option(const char* text, int& value) {
    if (text == NULL || *text == '\0') return false;
    char* end = NULL;
    errno = 0;
    const long parsed = std::strtol(text, &end, 10);
    if (errno == ERANGE || end == text || *end != '\0' || parsed < 0 ||
        parsed > std::numeric_limits<int>::max()) {
        return false;
    }
    value = static_cast<int>(parsed);
    return true;
}

void help(int code){
    fprintf(stderr, "downsample_vcf_parallel [OPTIONS]\n");
    fprintf(stderr, "Parallel version of downsample_vcf - uses OpenMP for chromosome-level parallelization.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Given a VCF file containing variants segregating within a panel of\n");
    fprintf(stderr, "deeply divergent genomes, intelligently downsamples the variant panel.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --vcf -v     A VCF/BCF file listing variants (must be indexed).\n");
    fprintf(stderr, "    --output -o  Output VCF file name. Will be gzipped.\n");
    fprintf(stderr, "    --num -n     Desired final number of SNPs (not required with --annotate_only).\n");
    fprintf(stderr, "                 Recommend > 1M for good clade coverage.\n");
    fprintf(stderr, "===== OPTIONAL =====\n");
    fprintf(stderr, "    --threads -t Number of threads (default: all available).\n");
    fprintf(stderr, "    --gtf -g     GTF annotation file. Adds ANNOT_SCORE INFO field.\n");
    fprintf(stderr, "                 Required if --output or --het_output is set.\n");
    fprintf(stderr, "    --cov -c     Tabix-indexed bedgraph coverage file (RNA). Adds COV_SCORE INFO field.\n");
    fprintf(stderr, "    --min_cov -m Minimum coverage threshold (default: 1.0). Sites below this are excluded.\n");
    fprintf(stderr, "    --seed -s    Random seed (currently unused; selection is deterministic).\n");
    fprintf(stderr, "    --blacklist  Additional BED regions to exclude (0-based half-open).\n");
    fprintf(stderr, "    --mito_chrom Contig name for the mitochondrial genome (default: chrM).\n");
    fprintf(stderr, "    --include_mito  Allow mitochondrial variants in nuclear panels (default: excluded).\n");
    fprintf(stderr, "    --numts_bed BED of nuclear mitochondrial insertions to exclude (0-based half-open).\n");
    fprintf(stderr, "                May also be supplied via CELLBOUNCER_NUMTS_BED.\n");
    fprintf(stderr, "    --include_numts  Do not apply --numts_bed/CELLBOUNCER_NUMTS_BED.\n");
    fprintf(stderr, "===== MITOCHONDRIAL FUSION-RATIO PANEL (optional) =====\n");
    fprintf(stderr, "    --mt_output FILE              Emit the dedicated chrM panel in this same VCF load.\n");
    fprintf(stderr, "    --mt_min_depth INT            Minimum parental RO+AO/AD depth (default: 10).\n");
    fprintf(stderr, "    --mt_homoplasmy_af FLOAT      Near-fixed parental allele fraction (default: 0.95).\n");
    fprintf(stderr, "    --mt_min_coverage FLOAT       chrM expression floor (default: inherit --min_cov).\n");
    fprintf(stderr, "    --mt_min_pair_sites INT       Requested ratio markers per library/pair (default: 20).\n");
    fprintf(stderr, "    --mt_min_ambient_sites INT    Requested ambient-only anchors per pair (default: 1).\n");
    fprintf(stderr, "    --mt_require_pair_targets     Fail after writing audits if either target is missed.\n");
    fprintf(stderr, "    --mt_allow_missing_coverage   Diagnostic only; permit no chrM coverage track.\n");
    fprintf(stderr, "    --mt_mask_bed FILE            Optional chrM-coordinate ambiguity mask.\n");
    fprintf(stderr, "    --mt_site_manifest FILE       Three-class per-library site manifest.\n");
    fprintf(stderr, "    --mt_pair_audit FILE          Per-library/pair marker audit.\n");
    fprintf(stderr, "    --mt_haplotype_groups FILE    Per-library mt haplotype equivalence groups.\n");
    fprintf(stderr, "    --mt_haplotype_pairwise FILE  Pairwise retained-panel distinction report.\n");
    fprintf(stderr, "    --mt_sample_audit FILE        VCF sample versus metadata role report.\n");
    fprintf(stderr, "    --mt_sites_bed FILE           BED for retaining chrM sites in subset BAMs.\n");
    fprintf(stderr, "===== HET SITE SELECTION =====\n");
    fprintf(stderr, "    --het_output -H  Output file for het-enriched SNP set (optional).\n");
    fprintf(stderr, "    --het_num -N     Target number of het sites (default: half of --num).\n");
    fprintf(stderr, "    --bin_size -b    Bin size (bp) for evenness scoring (default: 100000).\n");
    fprintf(stderr, "    --het_max_excess_z  Excess-het z-score threshold (default: 3.0). Sites above are excluded.\n");
    fprintf(stderr, "===== ATAC OUTPUTS (optional) =====\n");
    fprintf(stderr, "    --atac_cov -C        ATAC coverage bedGraph (tabix-indexed).\n");
    fprintf(stderr, "    --atac_output -O     ATAC-demux output BCF.\n");
    fprintf(stderr, "    --atac_num -M        Target N for ATAC-demux. Required if --atac_output set.\n");
    fprintf(stderr, "    --atac_het_output -T ATAC-het output BCF.\n");
    fprintf(stderr, "    --atac_het_num -U    Target N for ATAC-het (default: atac_num/2).\n");
    fprintf(stderr, "    --atac_candidates_output  Indexed BCF containing the reusable union of all\n");
    fprintf(stderr, "                              ATAC demux/het/species input candidates.\n");
    fprintf(stderr, "    --atac_candidates_only    Write that reusable candidate BCF and exit.\n");
    fprintf(stderr, "    --atac_min_cov       Minimum ATAC coverage threshold (default: same as --min_cov).\n");
    fprintf(stderr, "    --atac_exclude_vcf   Existing RNA panel BCF/VCF whose loci are forbidden in native ATAC outputs.\n");
    fprintf(stderr, "                         Repeatable legacy hard-exclusion input.\n");
    fprintf(stderr, "    --atac_soft_overlap_vcf Exact selected RNA-demux BCF/VCF whose loci remain eligible but receive a ranking penalty. Repeatable.\n");
    fprintf(stderr, "    --atac_soft_overlap_multiplier Multiplicative ranking factor for those loci (default: 0.95).\n");
    fprintf(stderr, "    --atac_audit_prefix  Prefix for compact candidate-SFS, pool-pair capacity, and species-tier TSV audits.\n");
    fprintf(stderr, "    --atac_regulatory_bed Optional BED used as an ATAC ranking prior.\n");
    fprintf(stderr, "    --atac_regulatory_required Restrict ATAC nuclear selection to the regulatory prior (default: off).\n");
    fprintf(stderr, "                          Without it, strand-aware 2-kb upstream promoters are derived from --gtf.\n");
    fprintf(stderr, "    --atac_numt_output   Separate ATAC-covered NUMT diagnostic BCF; NUMTs remain excluded from nuclear panels.\n");
    fprintf(stderr, "===== POOL-AWARE SCORING (optional) =====\n");
    fprintf(stderr, "    --enable_pool_scoring -E  Enable pool-aware het scoring.\n");
    fprintf(stderr, "    --pool_combinations -L    TSV: library_id<tab>comma-separated identities.\n");
    fprintf(stderr, "    --min_pair_score -Q       Threshold for pair discrimination (default: 0.3).\n");
    fprintf(stderr, "    --het_pool_weight         Weight on pool_discrim_norm when --enable_pool_scoring (default: 0.3).\n");
    fprintf(stderr, "                              Final het component = (1-w)*het_freq + w*pool_discrim_norm.\n");
    fprintf(stderr, "===== PAIRWISE FLOOR (optional) =====\n");
    fprintf(stderr, "    --min_pairwise -W    Generic/RNA: legacy individual-pair rescue may add sites beyond --num.\n");
    fprintf(stderr, "                         Native ATAC: not used as a selection floor; legal pool identities\n");
    fprintf(stderr, "                         contribute only to the RNA-style normalized scalar ranking score.\n");
    fprintf(stderr, "    --max_pairwise_rescue  Hard cap on total rescue variants (default: num/10).\n");
    fprintf(stderr, "    --pairwise_rescue_max_per_bin  Bin-spacing for rescue (default: species_max_per_bin).\n");
    fprintf(stderr, "    --pairwise_rescue_min_call_rate  Call-rate floor for rescue candidates (default: 0.5).\n");
    fprintf(stderr, "===== SPECIES OUTPUT (optional) =====\n");
    fprintf(stderr, "    --panel_metadata -P   TSV: indiv_id<tab>species.\n");
    fprintf(stderr, "    --species_output -S   Species discrimination output BCF.\n");
    fprintf(stderr, "    --species_num -R      Target N for species set. Required if --species_output set.\n");
    fprintf(stderr, "    --species_coverage -K Coverage source for species: 'rna' or 'atac' (default: rna).\n");
    fprintf(stderr, "    --species_outgroup Species label whose outgroup-only sites are capped during backfill.\n");
    fprintf(stderr, "    --species_outgroup_max_fraction Maximum fraction of the species budget assigned only to that outgroup.\n");
    fprintf(stderr, "    --species_min_pair_score  DEPRECATED: alias for --species_min_pair_discrim.\n");
    fprintf(stderr, "    --species_min_pair_discrim  Minimum pair_discrim for species selection (default: 0.5). Primary ranking score.\n");
    fprintf(stderr, "    --species_major_dosage_frac_floor  Per-species major-dosage-fraction floor (default: 0.7).\n");
    fprintf(stderr, "    --species_call_rate_floor          Per-species call-rate floor (default: 0.7).\n");
    fprintf(stderr, "    --species_min_called      Per-species floor on called individuals (default: 1; capped at n_in_species so singleton species can pass).\n");
    fprintf(stderr, "    --species_cov_cap         Coverage score cap (log2) for species ranking (default: 5.0).\n");
    fprintf(stderr, "    --species_max_per_bin     Max species sites per genomic bin (default: 5, 0 = unlimited).\n");
    fprintf(stderr, "    --max_het_per_bin         Max het sites per genomic bin (default: 10, 0 = unlimited).\n");
    fprintf(stderr, "===== ANNOTATION MODES =====\n");
    fprintf(stderr, "    --annotate_only  Annotate ALL sites with scores and exit (no downsampling).\n");
    fprintf(stderr, "    --annotate_all FILE  Annotate ALL sites first, save to FILE, then downsample.\n");
    fprintf(stderr, "    --help -h    Display this message and exit.\n");
    exit(code);
}

// proc_bcf_record() and proc_bcf_record_het() removed: all call sites now use
// precomputed dosage-aware bitsets from StoredVariant (loaded at VCF ingest time).

void count_branch(unordered_map<bitset<NBITS>, double>& branchcounts,
    bitset<NBITS>& clade,
    bitset<NBITS>& clade_flip,
    int total_alt,
    int total_alt_flip,
    double count){
   
    if (total_alt < total_alt_flip || (total_alt == total_alt_flip && clade < clade_flip)){
        if (branchcounts.count(clade) == 0){
            branchcounts.insert(make_pair(clade, count));
        }
        else{
            branchcounts[clade] += count;
        }
    }
    else{
        if (branchcounts.count(clade_flip) == 0){
            branchcounts.insert(make_pair(clade_flip, count));
        }
        else{
            branchcounts[clade_flip] += count;
        }
    }
}

void count_branch_missing(unordered_map<pair<bitset<NBITS>, bitset<NBITS> >, double>& branchcounts,
    unordered_map<bitset<NBITS>, bitset<NBITS> >& miss2flip,
    bitset<NBITS>& present,
    bitset<NBITS>& clade,
    bitset<NBITS>& clade_flip,
    int total_alt,
    int total_alt_flip,
    double count){
    
    pair<bitset<NBITS>, bitset<NBITS> > key;
    
    if (total_alt < total_alt_flip || (total_alt == total_alt_flip && clade < clade_flip)){
        key = make_pair(present, clade);
        if (miss2flip.count(clade) == 0){
            miss2flip.insert(make_pair(clade, clade_flip));
        }
    }
    else{
        key = make_pair(present, clade_flip);
        if (miss2flip.count(clade_flip) == 0){
            miss2flip.insert(make_pair(clade_flip, clade));
        }
    }
    if (branchcounts.count(key) == 0){
        branchcounts.insert(make_pair(key, count));
    }
    else{
        branchcounts[key] += count;
    }
}

// Merge thread-local branchcounts into global
void merge_branchcounts(
    unordered_map<bitset<NBITS>, double>& global,
    unordered_map<bitset<NBITS>, double>& local,
    mutex& mtx){
    
    lock_guard<mutex> lock(mtx);
    for (auto& kv : local){
        if (global.count(kv.first) == 0){
            global.insert(kv);
        }
        else{
            global[kv.first] += kv.second;
        }
    }
}

void merge_branchcounts_missing(
    unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, double>& global,
    unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, double>& local,
    mutex& mtx){
    
    lock_guard<mutex> lock(mtx);
    for (auto& kv : local){
        if (global.count(kv.first) == 0){
            global.insert(kv);
        }
        else{
            global[kv.first] += kv.second;
        }
    }
}

void merge_miss2flip(
    unordered_map<bitset<NBITS>, bitset<NBITS>>& global,
    unordered_map<bitset<NBITS>, bitset<NBITS>>& local,
    mutex& mtx){
    
    lock_guard<mutex> lock(mtx);
    for (auto& kv : local){
        if (global.count(kv.first) == 0){
            global.insert(kv);
        }
    }
}

// Merge thread-local het candidates into global
void merge_het_candidates(
    vector<HetCandidate>& global,
    vector<HetCandidate>& local,
    mutex& mtx){
    
    lock_guard<mutex> lock(mtx);
    global.insert(global.end(), 
        make_move_iterator(local.begin()), 
        make_move_iterator(local.end()));
}

// Merge thread-local species candidates into global
void merge_species_candidates(
    vector<SpeciesCandidate>& global,
    vector<SpeciesCandidate>& local,
    mutex& mtx){
    
    lock_guard<mutex> lock(mtx);
    global.insert(global.end(),
        make_move_iterator(local.begin()),
        make_move_iterator(local.end()));
}

// RNA-production pool score: mean expected-dosage separation across the legal
// identities listed in each real library.  The comparisons contribute to one
// scalar site-ranking score; they are never materialized as independent floor
// or rescue objectives.
float compute_pool_discrim(
    const vector<uint8_t>& genotypes_packed,
    int num_samples,
    const PoolCombinations& pool_combos,
    const map<string, int>& sample_name_to_idx,
    float min_pair_score) {
    
    // Compute alt_count per sample: 0/0=0, 0/1=1, 1/1=2, missing=-1
    vector<int> alt_count(num_samples, -1);
    for (int s = 0; s < num_samples; s++) {
        int a1 = gt_packed_allele(genotypes_packed, s, 0);
        int a2 = gt_packed_allele(genotypes_packed, s, 1);
        if (a1 < 0 || a2 < 0) {
            alt_count[s] = -1;
        } else {
            alt_count[s] = (a1 > 0 ? 1 : 0) + (a2 > 0 ? 1 : 0);
        }
    }
    
    double total_pair_sum = 0.0;
    long total_pair_count = 0;
    for (const auto& library : pool_combos.lib_identities) {
        const vector<PoolCombinations::CellIdentity>& identities =
            library.second;
        vector<float> expected;
        expected.reserve(identities.size());
        for (const auto& identity : identities) {
            if (identity.components.size() == 2) {
                auto first_found = sample_name_to_idx.find(identity.components[0]);
                auto second_found = sample_name_to_idx.find(identity.components[1]);
                if (first_found == sample_name_to_idx.end() ||
                    second_found == sample_name_to_idx.end()) continue;
                const int first = alt_count[first_found->second];
                const int second = alt_count[second_found->second];
                if (first < 0 || second < 0) continue;
                expected.push_back((float)(first + second) / 4.0f);
            } else {
                auto found = sample_name_to_idx.find(identity.components[0]);
                if (found == sample_name_to_idx.end()) continue;
                const int count = alt_count[found->second];
                if (count < 0) continue;
                expected.push_back((float)count / 2.0f);
            }
        }
        const long n_valid = (long)expected.size();
        if (n_valid < 2) {
            const long n_all = (long)identities.size();
            total_pair_count += (n_all * (n_all - 1)) / 2;
            continue;
        }
        total_pair_count += (n_valid * (n_valid - 1)) / 2;
        for (long i = 0; i < n_valid; ++i) {
            for (long j = i + 1; j < n_valid; ++j) {
                const float difference = fabsf(expected[i] - expected[j]);
                if (difference >= min_pair_score)
                    total_pair_sum += difference;
            }
        }
    }
    return total_pair_count > 0
        ? (float)(total_pair_sum / total_pair_count) : 0.0f;
}

// Compute species pair discrimination for a site.
// New scoring: pair_score = abs(mean_alt_A - mean_alt_B)/2
//              * call_rate_A * call_rate_B
//              * major_dosage_frac_A * major_dosage_frac_B
//
// Also outputs per-species stats for audit and the legacy raw pairwise fraction.
void compute_species_discrim(
    const vector<uint8_t>& genotypes_packed,
    int num_samples,
    const PanelMetadata& panel,
    const vector<int>& min_called_per_species,
    float major_dosage_frac_floor,
    float call_rate_floor,
    vector<float>& pair_score_out,
    vector<float>& pair_discrim_out,
    vector<float>& sp_mean_alt_out,
    vector<float>& sp_major_dosage_frac_out,
    vector<float>& sp_call_rate_out,
    vector<int>&   sp_n_called_out) {
    
    int n_sp = (int)panel.species_list.size();
    pair_score_out.assign(panel.n_pairs, 0.0f);
    pair_discrim_out.assign(panel.n_pairs, 0.0f);
    sp_mean_alt_out.assign(n_sp, 0.0f);
    sp_major_dosage_frac_out.assign(n_sp, 0.0f);
    sp_call_rate_out.assign(n_sp, 0.0f);
    sp_n_called_out.assign(n_sp, 0);

    // Per-sample alt count (-1 = missing). HALF-MISSING TREATED AS MISSING.
    vector<int> alt_count(num_samples, -1);
    for (int s = 0; s < num_samples; s++) {
        int a1 = gt_packed_allele(genotypes_packed, s, 0);
        int a2 = gt_packed_allele(genotypes_packed, s, 1);
        if (a1 < 0 || a2 < 0) {
            alt_count[s] = -1;
        } else {
            alt_count[s] = (a1 > 0 ? 1 : 0) + (a2 > 0 ? 1 : 0);
        }
    }

    // Per-species WEIGHTED statistics using PanelMetadata::get_weight().
    // Per-species weighted dosage counts are retained so the pair_discrim
    // computation below can use a closed-form O(3x3) calculation instead of
    // an O(N_a * N_b) double loop over individuals.
    vector<float> sp_total_weight(n_sp, 0.0f);
    vector<array<double, 3>> sp_wdc(n_sp, {0.0, 0.0, 0.0});  // weighted dosage counts per species
    for (int si = 0; si < n_sp; si++) {
        const string& sp = panel.species_list[si];
        auto it = panel.species_to_sample_indices.find(sp);
        if (it == panel.species_to_sample_indices.end()) continue;

        double weight_total = 0.0;
        double weight_called = 0.0;
        double weighted_alt_sum = 0.0;
        double weighted_dosage_counts[3] = {0.0, 0.0, 0.0};
        int    n_called_int = 0;

        for (int idx : it->second) {
            double w = panel.get_weight(sp, idx);
            weight_total += w;
            if (alt_count[idx] >= 0) {
                weight_called   += w;
                weighted_alt_sum += w * alt_count[idx];
                weighted_dosage_counts[alt_count[idx]] += w;
                n_called_int++;
            }
        }

        sp_total_weight[si] = (float)weight_total;
        sp_n_called_out[si] = n_called_int;
        sp_wdc[si][0] = weighted_dosage_counts[0];
        sp_wdc[si][1] = weighted_dosage_counts[1];
        sp_wdc[si][2] = weighted_dosage_counts[2];

        if (weight_called > 0.0) {
            sp_mean_alt_out[si] = (float)(weighted_alt_sum / weight_called);
            double mx = max({weighted_dosage_counts[0],
                             weighted_dosage_counts[1],
                             weighted_dosage_counts[2]});
            sp_major_dosage_frac_out[si] = (float)(mx / weight_called);
        }
        sp_call_rate_out[si] = (weight_total > 0.0)
                               ? (float)(weight_called / weight_total) : 0.0f;
    }

    // For each species pair: compute pair_discrim (primary score) and
    // legacy multiplicative pair_score (recorded only). Apply per-species
    // gates: n_called >= min_called_per_species[si], call_rate >= floor,
    // major_dosage_frac >= floor.
    for (int pi = 0; pi < panel.n_pairs; pi++) {
        const string& sp_a = panel.species_pairs[pi].first;
        const string& sp_b = panel.species_pairs[pi].second;

        int si_a = -1, si_b = -1;
        for (int si = 0; si < n_sp; si++) {
            if (panel.species_list[si] == sp_a) si_a = si;
            if (panel.species_list[si] == sp_b) si_b = si;
        }
        if (si_a < 0 || si_b < 0) continue;

        // Hard gates
        if (sp_n_called_out[si_a] < min_called_per_species[si_a]) continue;
        if (sp_n_called_out[si_b] < min_called_per_species[si_b]) continue;
        if (sp_call_rate_out[si_a] < call_rate_floor) continue;
        if (sp_call_rate_out[si_b] < call_rate_floor) continue;
        if (sp_major_dosage_frac_out[si_a] < major_dosage_frac_floor) continue;
        if (sp_major_dosage_frac_out[si_b] < major_dosage_frac_floor) continue;

        // pair_discrim: weighted fraction of cross-species (A,B) individual
        // pairs that are distinguishable (different alt counts).
        // Closed form: sum_{i,j in {0,1,2}, i!=j} wdc_A[i] * wdc_B[j]
        //              divided by (sum wdc_A) * (sum wdc_B).
        // This is mathematically identical to the O(N_a * N_b) double loop
        // (each individual contributes its weight to exactly one dosage bin)
        // but runs in O(3 * 3) per pair instead of O(N_a * N_b).
        auto it_a = panel.species_to_sample_indices.find(sp_a);
        auto it_b = panel.species_to_sample_indices.find(sp_b);
        if (it_a == panel.species_to_sample_indices.end() ||
            it_b == panel.species_to_sample_indices.end()) continue;

        double num = 0.0;
        double denom = 0.0;
        for (int i = 0; i < 3; i++) {
            double wa_i = sp_wdc[si_a][i];
            if (wa_i <= 0.0) continue;
            for (int j = 0; j < 3; j++) {
                double wb_j = sp_wdc[si_b][j];
                if (wb_j <= 0.0) continue;
                double w = wa_i * wb_j;
                denom += w;
                if (i != j) num += w;
            }
        }
        pair_discrim_out[pi] = (denom > 0.0) ? (float)(num / denom) : 0.0f;

        // Legacy multiplicative score (recorded for backward-compat output)
        float between = fabsf(sp_mean_alt_out[si_a] - sp_mean_alt_out[si_b]) / 2.0f;
        float consistency = sp_major_dosage_frac_out[si_a] * sp_major_dosage_frac_out[si_b];
        float callrate = sp_call_rate_out[si_a] * sp_call_rate_out[si_b];
        pair_score_out[pi] = between * consistency * callrate;
    }
}

double brentfun(double param, const map<string, double>& data_d, const map<string, int>& data_i){
    int num = data_i.at("num");
    char buf[50];
    string bufstr;
    double sum = 0.0;
    for (int i = 0; i < num; ++i){
        snprintf(buf, sizeof(buf), "x%d", i);
        bufstr = buf;
        double count = data_d.at(bufstr);
        sum += pow(count, param);
    }
    return sum - (double)data_i.at("target");
}

double dbrentfun(double param, const map<string, double>& data_d, const map<string, int>& data_i){
    int num = data_i.at("num");
    char buf[50];
    string bufstr;
    double d_df = 0.0;
    for (int i = 0; i < num; ++i){
        snprintf(buf, sizeof(buf), "x%d", i);
        bufstr = buf;
        double count = data_d.at(bufstr);
        d_df += pow(count, param) * log(count);
    }
    return d_df;
}

// Helper to get bin index from position
inline int get_bin_idx(int64_t pos, int bin_size) {
    return (int)(pos / bin_size);
}

// Helper to create BIN_ID string
inline string make_bin_id(const string& chrom, int bin_idx) {
    return chrom + ":" + to_string(bin_idx);
}

// Compute within-clade selection score for Option C stratified selection
// Coverage dominates, annotation is a small tiebreaker bonus
inline float compute_selection_score(float cov_score, float annot_score) {
    // annot_score is 1.0 for genic, 0.1 for intergenic
    // Convert to annot_boost: 1.0 if genic, 0.0 if intergenic
    float annot_boost = (annot_score >= 1.0f) ? 1.0f : 0.0f;
    return log2f(cov_score + 1.0f) + (0.1f * annot_boost);
}

// Structure for tracking SNPs within a clade for stratified selection
struct CladeCandidate {
    size_t variant_idx;      // Index into all_variants
    float selection_score;   // log2(cov+1) + 0.1*annot_boost
    
    bool operator<(const CladeCandidate& other) const {
        return selection_score > other.selection_score;  // Descending order
    }
};

// Load a compact position-only exclusion set from existing panel BCF/VCFs.
// Keys use the source VCF contig index in the high 32 bits and the 0-based
// position in the low 32 bits. Position-only matching is intentional: a new
// ATAC marker must not reuse an RNA locus under a different allele encoding.
vector<uint64_t> load_panel_exclusion_positions(
    const vector<string>& filenames,
    const map<string, int>& source_chrom_to_idx,
    int n_threads) {
    vector<uint64_t> positions;
    for (const string& filename : filenames) {
        fprintf(stderr, "Loading RNA-panel exclusion loci from %s...\n",
            filename.c_str());
        htsFile* input = bcf_open(filename.c_str(), "r");
        if (!input) {
            fprintf(stderr, "ERROR: Could not open ATAC exclusion panel: %s\n",
                filename.c_str());
            exit(1);
        }
        hts_set_threads(input, max(1, min(n_threads, 8)));
        bcf_hdr_t* header = bcf_hdr_read(input);
        if (!header) {
            fprintf(stderr, "ERROR: Could not read ATAC exclusion panel header: %s\n",
                filename.c_str());
            hts_close(input);
            exit(1);
        }
        bcf1_t* record = bcf_init();
        uint64_t file_count = 0;
        while (bcf_read(input, header, record) >= 0) {
            const char* chrom = bcf_hdr_id2name(header, record->rid);
            map<string, int>::const_iterator found =
                source_chrom_to_idx.find(chrom ? chrom : "");
            if (found == source_chrom_to_idx.end()) continue;
            if (record->pos < 0 ||
                static_cast<uint64_t>(record->pos) > 0xffffffffULL) {
                fprintf(stderr,
                    "ERROR: exclusion-panel position is outside the supported "
                    "32-bit coordinate range: %s:%lld\n",
                    chrom, static_cast<long long>(record->pos + 1));
                bcf_destroy(record);
                bcf_hdr_destroy(header);
                hts_close(input);
                exit(1);
            }
            const uint64_t key =
                (static_cast<uint64_t>(static_cast<uint32_t>(found->second)) << 32) |
                static_cast<uint32_t>(record->pos);
            positions.push_back(key);
            ++file_count;
        }
        bcf_destroy(record);
        bcf_hdr_destroy(header);
        hts_close(input);
        fprintf(stderr, "  Read %llu exclusion records\n",
            static_cast<unsigned long long>(file_count));
    }
    sort(positions.begin(), positions.end());
    positions.erase(unique(positions.begin(), positions.end()), positions.end());
    fprintf(stderr, "ATAC/RNA exclusion union: %zu unique loci from %zu panel(s)\n",
        positions.size(), filenames.size());
    return positions;
}

inline bool panel_position_is_excluded(
    const vector<uint64_t>& sorted_positions,
    int chrom_idx,
    int64_t pos) {
    if (sorted_positions.empty() || chrom_idx < 0 || pos < 0) return false;
    const uint64_t key =
        (static_cast<uint64_t>(static_cast<uint32_t>(chrom_idx)) << 32) |
        static_cast<uint32_t>(pos);
    return binary_search(sorted_positions.begin(), sorted_positions.end(), key);
}


int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"vcf", required_argument, 0, 'v'},
       {"output", required_argument, 0, 'o'},
       {"num", required_argument, 0, 'n'},
       {"threads", required_argument, 0, 't'},
       {"gtf", required_argument, 0, 'g'},
       {"cov", required_argument, 0, 'c'},
       {"seed", required_argument, 0, 's'},
       {"het_output", required_argument, 0, 'H'},
       {"het_num", required_argument, 0, 'N'},
       {"bin_size", required_argument, 0, 'b'},
       {"min_cov", required_argument, 0, 'm'},
       {"annotate_only", no_argument, 0, 'A'},
       {"annotate_all", required_argument, 0, 'a'},
       {"atac_cov", required_argument, 0, 'C'},
       {"panel_metadata", required_argument, 0, 'P'},
       {"pool_combinations", required_argument, 0, 'L'},
       {"enable_pool_scoring", no_argument, 0, 'E'},
       {"atac_output", required_argument, 0, 'O'},
       {"atac_num", required_argument, 0, 'M'},
       {"atac_het_output", required_argument, 0, 'T'},
       {"atac_het_num", required_argument, 0, 'U'},
       {"atac_min_cov", required_argument, 0, 1204},
       {"atac_exclude_vcf", required_argument, 0, 1500},
       {"atac_candidates_output", required_argument, 0, 1501},
       {"atac_candidates_only", no_argument, 0, 1502},
       {"atac_regulatory_bed", required_argument, 0, 1503},
       {"atac_numt_output", required_argument, 0, 1504},
       {"species_outgroup", required_argument, 0, 1505},
       {"species_outgroup_max_fraction", required_argument, 0, 1506},
       {"atac_regulatory_required", no_argument, 0, 1507},
       {"atac_soft_overlap_vcf", required_argument, 0, 1508},
       {"atac_soft_overlap_multiplier", required_argument, 0, 1509},
       {"atac_audit_prefix", required_argument, 0, 1510},
       {"species_output", required_argument, 0, 'S'},
       {"species_num", required_argument, 0, 'R'},
       {"species_coverage", required_argument, 0, 'K'},
       {"species_min_pair_score", required_argument, 0, 1100},
       {"species_min_pair_discrim", required_argument, 0, 1106},
       {"species_major_dosage_frac_floor", required_argument, 0, 1104},
       {"species_call_rate_floor", required_argument, 0, 1105},
       {"species_min_called", required_argument, 0, 1101},
       {"species_cov_cap", required_argument, 0, 1102},
       {"species_max_per_bin", required_argument, 0, 1103},
       {"max_het_per_bin", required_argument, 0, 1200},
       {"blacklist", required_argument, 0, 1203},
       {"numts_bed", required_argument, 0, 1205},
       {"include_mito", no_argument, 0, 1206},
       {"include_numts", no_argument, 0, 1207},
       {"mito_chrom", required_argument, 0, 1208},
       {"mt_output", required_argument, 0, 1400},
       {"mt_site_manifest", required_argument, 0, 1401},
       {"mt_pair_audit", required_argument, 0, 1402},
       {"mt_haplotype_groups", required_argument, 0, 1403},
       {"mt_haplotype_pairwise", required_argument, 0, 1404},
       {"mt_sample_audit", required_argument, 0, 1405},
       {"mt_sites_bed", required_argument, 0, 1406},
       {"mt_mask_bed", required_argument, 0, 1407},
       {"mt_min_depth", required_argument, 0, 1408},
       {"mt_homoplasmy_af", required_argument, 0, 1409},
       {"mt_min_coverage", required_argument, 0, 1410},
       {"mt_min_pair_sites", required_argument, 0, 1411},
       {"mt_min_ambient_sites", required_argument, 0, 1412},
       {"mt_require_pair_targets", no_argument, 0, 1413},
       {"mt_allow_missing_coverage", no_argument, 0, 1414},
       {"het_pool_weight", required_argument, 0, 1201},
       {"het_max_excess_z", required_argument, 0, 1202},
       {"min_pair_score", required_argument, 0, 'Q'},
       {"min_pairwise", required_argument, 0, 'W'},
       {"max_pairwise_rescue", required_argument, 0, 1300},
       {"pairwise_rescue_max_per_bin", required_argument, 0, 1301},
       {"pairwise_rescue_min_call_rate", required_argument, 0, 1302},
       {0, 0, 0, 0} 
    };
    
    string vcf_file = "";
    string outfile = "";
    string gtf_file = "";
    string cov_file = "";
    string het_outfile = "";
    string annotate_all_file = "";
    string atac_cov_file = "";
    string panel_metadata_file = "";
    string pool_combinations_file = "";
    string atac_outfile = "";
    string atac_het_outfile = "";
    string atac_candidates_outfile = "";
    string atac_regulatory_bed_file = "";
    string atac_numt_outfile = "";
    vector<string> atac_exclude_vcfs;
    vector<string> atac_soft_overlap_vcfs;
    string atac_audit_prefix = "";
    string species_outgroup = "";
    string species_outfile = "";
    string blacklist_file = "";  // Additional user blacklist
    string numts_bed_file = "";
    string mito_chrom = "chrM";
    string mt_outfile = "";
    string mt_site_manifest_file = "";
    string mt_pair_audit_file = "";
    string mt_haplotype_groups_file = "";
    string mt_haplotype_pairwise_file = "";
    string mt_sample_audit_file = "";
    string mt_sites_bed_file = "";
    string mt_mask_bed_file = "";
    string species_coverage_mode = "rna";
    int num = -1;
    int het_num = -1;
    int atac_num = -1;
    int atac_het_num = -1;
    int species_num = -1;
    int bin_size = 100000;
    float min_cov = 1.0f;
    float atac_min_cov = -1.0f;  // B14: -1 = default to min_cov
    float min_pair_score = 0.3f;
    float species_min_pair_score = 0.05f;
    float species_min_pair_discrim = 0.5f;
    float species_major_dosage_frac_floor = 0.7f;
    float species_call_rate_floor = 0.7f;
    int species_min_called = 1;
    float species_cov_cap = 5.0f;
    float species_outgroup_max_fraction = 1.0f;
    float atac_soft_overlap_multiplier = 0.95f;
    int species_max_per_bin = 5;  // nonzero default avoids dense regional clusters
    int max_het_per_bin = 10;     // het bin cap (0 = unlimited)
    float het_pool_weight = 0.3f; // weight on pool_discrim_norm in het score
    float het_max_excess_z = 3.0f; // B7: excess-het z-score filter threshold
    long min_pairwise = 0;  // 0 = disabled
    long max_pairwise_rescue = -1;  // B8: -1 = default to num/10
    int pairwise_rescue_max_per_bin = -1;  // B8: -1 = use species_max_per_bin
    float pairwise_rescue_min_call_rate = 0.5f;  // B8
    int mt_min_depth = 10;
    int mt_min_pair_sites = 20;
    int mt_min_ambient_sites = 1;
    float mt_homoplasmy_af = 0.95f;
    float mt_min_coverage = -1.0f;  // -1 = inherit --min_cov
    int n_threads = omp_get_max_threads();
    int seed = -1;
    bool annotate_only = false;
    bool enable_pool_scoring = false;
    bool include_mito = false;
    bool include_numts = false;
    bool mt_require_pair_targets = false;
    bool mt_require_coverage = true;
    bool atac_candidates_only = false;
    bool atac_regulatory_required = false;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:n:o:t:g:c:s:H:N:b:m:Aa:hC:P:L:EO:M:T:U:S:R:K:Q:W:", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                break;
            case 'h':
                help(0);
                break;
            case 'v':
                vcf_file = optarg;
                break;
            case 'n':
                num = atoi(optarg);
                break;
            case 'o':
                outfile = optarg;
                break;
            case 't':
                n_threads = atoi(optarg);
                break;
            case 'g':
                gtf_file = optarg;
                break;
            case 'c':
                cov_file = optarg;
                break;
            case 's':
                seed = atoi(optarg);
                break;
            case 'H':
                het_outfile = optarg;
                break;
            case 'N':
                het_num = atoi(optarg);
                break;
            case 'b':
                bin_size = atoi(optarg);
                break;
            case 'm':
                min_cov = atof(optarg);
                break;
            case 'A':
                annotate_only = true;
                break;
            case 'a':
                annotate_all_file = optarg;
                break;
            case 'C':
                atac_cov_file = optarg;
                break;
            case 'P':
                panel_metadata_file = optarg;
                break;
            case 'L':
                pool_combinations_file = optarg;
                break;
            case 'E':
                enable_pool_scoring = true;
                break;
            case 'O':
                atac_outfile = optarg;
                break;
            case 'M':
                atac_num = atoi(optarg);
                break;
            case 'T':
                atac_het_outfile = optarg;
                break;
            case 'U':
                atac_het_num = atoi(optarg);
                break;
            case 1204:
                atac_min_cov = atof(optarg);
                break;
            case 1500:
                atac_exclude_vcfs.push_back(optarg);
                break;
            case 1501:
                atac_candidates_outfile = optarg;
                break;
            case 1502:
                atac_candidates_only = true;
                break;
            case 1503:
                atac_regulatory_bed_file = optarg;
                break;
            case 1504:
                atac_numt_outfile = optarg;
                break;
            case 1505:
                species_outgroup = optarg;
                break;
            case 1506:
                species_outgroup_max_fraction = atof(optarg);
                break;
            case 1507:
                atac_regulatory_required = true;
                break;
            case 1508:
                atac_soft_overlap_vcfs.push_back(optarg);
                break;
            case 1509:
                atac_soft_overlap_multiplier = atof(optarg);
                break;
            case 1510:
                atac_audit_prefix = optarg;
                break;
            case 'S':
                species_outfile = optarg;
                break;
            case 'R':
                species_num = atoi(optarg);
                break;
            case 'K':
                species_coverage_mode = optarg;
                break;
            case 1100:
                species_min_pair_score = atof(optarg);
                species_min_pair_discrim = atof(optarg);
                fprintf(stderr,
                    "WARNING: --species_min_pair_score is deprecated; use "
                    "--species_min_pair_discrim. Both set the same threshold.\n");
                break;
            case 1101:
                species_min_called = atoi(optarg);
                break;
            case 1102:
                species_cov_cap = atof(optarg);
                break;
            case 1103:
                if (!parse_nonnegative_int_option(optarg, species_max_per_bin)) {
                    fprintf(stderr,
                        "ERROR: --species_max_per_bin requires a nonnegative integer "
                        "(0 = unlimited); got '%s'\n", optarg);
                    exit(1);
                }
                break;
            case 1104:
                species_major_dosage_frac_floor = atof(optarg);
                break;
            case 1105:
                species_call_rate_floor = atof(optarg);
                break;
            case 1106:
                species_min_pair_discrim = atof(optarg);
                break;
            case 1200:
                if (!parse_nonnegative_int_option(optarg, max_het_per_bin)) {
                    fprintf(stderr,
                        "ERROR: --max_het_per_bin requires a nonnegative integer "
                        "(0 = unlimited); got '%s'\n", optarg);
                    exit(1);
                }
                break;
            case 1203:
                blacklist_file = optarg;
                break;
            case 1205:
                numts_bed_file = optarg;
                break;
            case 1206:
                include_mito = true;
                break;
            case 1207:
                include_numts = true;
                break;
            case 1208:
                mito_chrom = optarg;
                break;
            case 1400:
                mt_outfile = optarg;
                break;
            case 1401:
                mt_site_manifest_file = optarg;
                break;
            case 1402:
                mt_pair_audit_file = optarg;
                break;
            case 1403:
                mt_haplotype_groups_file = optarg;
                break;
            case 1404:
                mt_haplotype_pairwise_file = optarg;
                break;
            case 1405:
                mt_sample_audit_file = optarg;
                break;
            case 1406:
                mt_sites_bed_file = optarg;
                break;
            case 1407:
                mt_mask_bed_file = optarg;
                break;
            case 1408:
                mt_min_depth = atoi(optarg);
                break;
            case 1409:
                mt_homoplasmy_af = atof(optarg);
                break;
            case 1410:
                mt_min_coverage = atof(optarg);
                break;
            case 1411:
                mt_min_pair_sites = atoi(optarg);
                break;
            case 1412:
                mt_min_ambient_sites = atoi(optarg);
                break;
            case 1413:
                mt_require_pair_targets = true;
                break;
            case 1414:
                mt_require_coverage = false;
                break;
            case 1201:
                het_pool_weight = atof(optarg);
                // C10: het_pool_weight is clamped to [0, 1]. pool_discrim_norm is bounded
                // in [0, 1] because each contributing |expected_alts[i] - expected_alts[j]|
                // is bounded by 1 for biallelic dosage (max diff = |0 - 1| = 1).
                if (het_pool_weight < 0.0f) het_pool_weight = 0.0f;
                if (het_pool_weight > 1.0f) het_pool_weight = 1.0f;
                break;
            case 1202:
                het_max_excess_z = atof(optarg);
                break;
            case 'Q':
                min_pair_score = atof(optarg);
                break;
            case 'W':
                min_pairwise = atol(optarg);
                break;
            case 1300:
                max_pairwise_rescue = atol(optarg);
                break;
            case 1301:
                pairwise_rescue_max_per_bin = atoi(optarg);
                break;
            case 1302:
                pairwise_rescue_min_call_rate = atof(optarg);
                break;
            default:
                help(0);
                break;
        }    
    }
    
    if (numts_bed_file.empty()) {
        const char* env_numts = getenv("CELLBOUNCER_NUMTS_BED");
        if (env_numts && *env_numts) numts_bed_file = env_numts;
    }

    if (vcf_file == ""){
        fprintf(stderr, "ERROR: VCF file required\n");
        exit(1);
    }
    const bool has_native_atac_output = !atac_outfile.empty()
                                     || !atac_het_outfile.empty()
                                     || !atac_candidates_outfile.empty()
                                     || !atac_numt_outfile.empty();
    const bool wants_generic_demux = !outfile.empty();
    const bool wants_generic_allocation = wants_generic_demux
                                        || !het_outfile.empty()
                                        || !annotate_all_file.empty();

    if (num <= 0 && !annotate_only && !atac_candidates_only &&
        wants_generic_allocation){
        fprintf(stderr, "ERROR: --num/-n must be positive\n");
        exit(1);
    }
    if (num > 0 && num < 100000 && !annotate_only) {
        fprintf(stderr, "WARNING: --num=%d is quite small; recommend > 100K for adequate "
            "clade coverage. Continue at your own risk.\n", num);
    }
    if (outfile == "" && !atac_candidates_only && !has_native_atac_output){
        fprintf(stderr, "ERROR: at least one output is required.\n");
        exit(1);
    }
    if (!outfile.empty() && outfile.rfind(".bcf") == string::npos && outfile.rfind(".vcf") == string::npos){
        outfile += ".bcf";
    }

    if (atac_candidates_only && atac_candidates_outfile.empty()) {
        fprintf(stderr,
            "ERROR: --atac_candidates_only requires --atac_candidates_output\n");
        exit(1);
    }
    
    // Handle annotate_all output file extension
    if (!annotate_all_file.empty()) {
        if (annotate_all_file.rfind(".bcf") == string::npos && annotate_all_file.rfind(".vcf") == string::npos){
            annotate_all_file += ".bcf";
        }
    }
    
    // Set het_num default to half of num if not specified
    if (het_num < 0) {
        het_num = num / 2;
    }
    
    // Handle het output file extension
    if (!het_outfile.empty()) {
        if (het_outfile.rfind(".bcf") == string::npos && het_outfile.rfind(".vcf") == string::npos){
            het_outfile += ".het.bcf";
        }
    }
    
    // Validate new flag dependencies
    if ((wants_generic_demux || !het_outfile.empty()) &&
        gtf_file.empty() && !annotate_only) {
        fprintf(stderr, "ERROR: --gtf/-g is required when --output or --het_output is set\n");
        exit(1);
    }
    
    if (!atac_outfile.empty() || !atac_het_outfile.empty() ||
        !atac_candidates_outfile.empty() || !atac_numt_outfile.empty()) {
        if (atac_cov_file.empty()) {
            fprintf(stderr, "ERROR: --atac_cov/-C is required for ATAC outputs\n");
            exit(1);
        }
    }
    if (!atac_outfile.empty() && atac_num <= 0) {
        fprintf(stderr, "ERROR: --atac_num/-M is required when --atac_output is set\n");
        exit(1);
    }
    if (!atac_het_outfile.empty() && atac_het_num < 0) {
        if (atac_num > 0) {
            atac_het_num = atac_num / 2;
        } else {
            fprintf(stderr, "ERROR: --atac_het_num/-U is required when --atac_het_output is set\n");
            exit(1);
        }
    }
    
    // D8: Require --cov when RNA outputs are requested and min_cov > 0
    // annotate_only mode can annotate without coverage (fields will be 0/absent)
    bool wants_rna_output = !annotate_only
                           && (wants_generic_demux
                               || !het_outfile.empty()
                               || (!species_outfile.empty() && species_coverage_mode == "rna"));
    if (wants_rna_output && min_cov > 0.0f && cov_file.empty()) {
        fprintf(stderr,
            "ERROR: --cov/-c is required when --output, --het_output, or RNA-mode "
            "--species_output is set and --min_cov > 0.\n"
            "       To run without RNA coverage filtering, pass --min_cov 0.\n");
        exit(1);
    }
    
    if (enable_pool_scoring && pool_combinations_file.empty()) {
        fprintf(stderr, "ERROR: --pool_combinations/-L is required when --enable_pool_scoring is set\n");
        exit(1);
    }
    if (!atac_outfile.empty() && pool_combinations_file.empty()) {
        fprintf(stderr,
            "ERROR: native --atac_output requires --pool_combinations so "
            "fusion-aware eligibility and rescue have a defined operational "
            "objective universe\n");
        exit(1);
    }
    if (atac_soft_overlap_multiplier < 0.0f ||
        atac_soft_overlap_multiplier > 1.0f) {
        fprintf(stderr,
            "ERROR: --atac_soft_overlap_multiplier must be in [0,1]\n");
        exit(1);
    }
    if (!atac_audit_prefix.empty()) {
        if (atac_cov_file.empty() || pool_combinations_file.empty() ||
            panel_metadata_file.empty()) {
            fprintf(stderr,
                "ERROR: --atac_audit_prefix requires --atac_cov, "
                "--pool_combinations, and --panel_metadata\n");
            exit(1);
        }
    }
    if (!mt_outfile.empty()) {
        if (annotate_only) {
            fprintf(stderr, "ERROR: --mt_output is not supported with --annotate_only; run the normal multi-output path\n");
            exit(1);
        }
        if (pool_combinations_file.empty()) {
            fprintf(stderr, "ERROR: --pool_combinations/-L is required when --mt_output is set\n");
            exit(1);
        }
        if (mt_min_depth < 1 || mt_homoplasmy_af <= 0.5f || mt_homoplasmy_af > 1.0f ||
            mt_min_pair_sites < 0 || mt_min_ambient_sites < 0) {
            fprintf(stderr, "ERROR: invalid mitochondrial panel threshold\n");
            exit(1);
        }
        if (mt_min_coverage < 0.0f) mt_min_coverage = min_cov;
        if (mt_min_coverage < 0.0f) {
            fprintf(stderr, "ERROR: --mt_min_coverage must be nonnegative\n");
            exit(1);
        }
        if (mt_require_coverage && cov_file.empty() && atac_cov_file.empty()) {
            fprintf(stderr,
                "ERROR: --mt_output requires --cov or --atac_cov with chrM coverage. "
                "Use --mt_allow_missing_coverage only for an explicit diagnostic run.\n");
            exit(1);
        }
        if (mt_outfile.rfind(".bcf") == string::npos && mt_outfile.rfind(".vcf") == string::npos) {
            mt_outfile += ".bcf";
        }
    }
    
    if (!species_outfile.empty()) {
        if (panel_metadata_file.empty()) {
            fprintf(stderr, "ERROR: --panel_metadata/-P is required when --species_output is set\n");
            exit(1);
        }
        if (species_num <= 0) {
            fprintf(stderr, "ERROR: --species_num/-R is required when --species_output is set\n");
            exit(1);
        }
        if (species_coverage_mode != "rna" && species_coverage_mode != "atac") {
            fprintf(stderr, "ERROR: --species_coverage must be 'rna' or 'atac'\n");
            exit(1);
        }
        if (species_coverage_mode == "atac" && atac_cov_file.empty()) {
            fprintf(stderr, "ERROR: --atac_cov/-C is required when --species_coverage is 'atac'\n");
            exit(1);
        }
        if (species_outgroup_max_fraction < 0.0f ||
            species_outgroup_max_fraction > 1.0f) {
            fprintf(stderr,
                "ERROR: --species_outgroup_max_fraction must be in [0,1]\n");
            exit(1);
        }
        if (!species_outgroup.empty() && species_outgroup_max_fraction >= 1.0f) {
            fprintf(stderr,
                "WARNING: --species_outgroup is set but its maximum fraction is 1.0; "
                "outgroup-only backfill is effectively uncapped.\n");
        }
    }
    
    // Handle ATAC output file extensions
    if (!atac_outfile.empty()) {
        if (atac_outfile.rfind(".bcf") == string::npos && atac_outfile.rfind(".vcf") == string::npos)
            atac_outfile += ".bcf";
    }
    if (!atac_het_outfile.empty()) {
        if (atac_het_outfile.rfind(".bcf") == string::npos && atac_het_outfile.rfind(".vcf") == string::npos)
            atac_het_outfile += ".bcf";
    }
    if (!atac_candidates_outfile.empty()) {
        if (atac_candidates_outfile.rfind(".bcf") == string::npos &&
            atac_candidates_outfile.rfind(".vcf") == string::npos)
            atac_candidates_outfile += ".bcf";
    }
    if (!atac_numt_outfile.empty()) {
        if (atac_numt_outfile.rfind(".bcf") == string::npos &&
            atac_numt_outfile.rfind(".vcf") == string::npos)
            atac_numt_outfile += ".bcf";
    }
    if (!species_outfile.empty()) {
        if (species_outfile.rfind(".bcf") == string::npos && species_outfile.rfind(".vcf") == string::npos)
            species_outfile += ".bcf";
    }

    fprintf(stderr, "Using %d threads\n", n_threads);
    fprintf(stderr, "Bin size for evenness scoring: %d bp\n", bin_size);
    fprintf(stderr, "Minimum coverage threshold: %.1f\n", min_cov);
    fprintf(stderr, "Mitochondrial contig: %s (%s in nuclear panels)\n",
        mito_chrom.c_str(), include_mito ? "included" : "excluded by default");
    if (!mt_outfile.empty()) {
        fprintf(stderr,
            "Dedicated mt panel: %s (same VCF load; min_depth=%d, homoplasmy_af=%.3f, "
            "min_coverage=%.3f, ratio_target=%d, ambient_anchor_target=%d)\n",
            mt_outfile.c_str(), mt_min_depth, mt_homoplasmy_af, mt_min_coverage,
            mt_min_pair_sites, mt_min_ambient_sites);
    }
    if (include_numts) {
        fprintf(stderr, "NUMT exclusion: disabled by --include_numts\n");
    } else if (!numts_bed_file.empty()) {
        fprintf(stderr, "NUMT exclusion BED: %s\n", numts_bed_file.c_str());
    } else {
        fprintf(stderr,
            "WARNING: no NUMT BED was supplied. chrM is still excluded, but nuclear "
            "NUMT loci cannot be removed. Use --numts_bed or CELLBOUNCER_NUMTS_BED.\n");
    }
    if (atac_min_cov < 0.0f) atac_min_cov = min_cov;  // B14: default to min_cov
    if (seed != -1) {
        fprintf(stderr, "WARNING: --seed is currently unused; selection is fully "
            "deterministic regardless of seed value.\n");
    }
    if (min_pairwise > 0) {
        if (!atac_outfile.empty() && !wants_generic_allocation) {
            fprintf(stderr,
                "Native ATAC identity-pair floors are disabled; --min_pairwise=%ld "
                "is retained only in diagnostic metadata\n", min_pairwise);
        } else {
            fprintf(stderr,
                "Generic/RNA pairwise floor: %ld distinguishing SNPs per "
                "individual pair\n", min_pairwise);
        }
        // B8: Set defaults for rescue caps
        if (max_pairwise_rescue < 0) {
            const int rescue_budget_basis = wants_generic_allocation ? num : atac_num;
            max_pairwise_rescue = rescue_budget_basis / 10;
        }
        if (pairwise_rescue_max_per_bin < 0) pairwise_rescue_max_per_bin = species_max_per_bin;
        fprintf(stderr, "  Rescue caps: max_total=%ld, max_per_bin=%d, min_call_rate=%.2f\n",
            max_pairwise_rescue, pairwise_rescue_max_per_bin, pairwise_rescue_min_call_rate);
    }
    if (annotate_only) {
        fprintf(stderr, "Mode: ANNOTATE ONLY (all sites will be annotated, no downsampling)\n");
    } else if (!annotate_all_file.empty()) {
        fprintf(stderr, "Mode: ANNOTATE ALL first -> %s, then downsample\n", annotate_all_file.c_str());
    }
    if (!het_outfile.empty() && !annotate_only) {
        fprintf(stderr, "Het output enabled: %s (target: %d sites)\n", het_outfile.c_str(), het_num);
    }
    omp_set_num_threads(n_threads);

    // ========================================================================
    // Load optional annotation and coverage files
    // ========================================================================
    map<string, vector<Interval>> annotations;
    if (!gtf_file.empty() && (annotate_only || wants_generic_allocation)) {
        load_gtf(gtf_file, annotations);
    }

    // Load optional interval exclusions. NUMTs are kept separate from the
    // generic blacklist so annotate_all can report the exact exclusion reason.
    BlacklistExclusions blacklist;
    if (!blacklist_file.empty()) {
        blacklist = load_blacklist(blacklist_file);
    }
    map<string, vector<Interval>> numts;
    if (!include_numts && !numts_bed_file.empty()) {
        load_bed(numts_bed_file, numts, "NUMT");
    }
    map<string, vector<Interval>> atac_regulatory_regions;
    string atac_regulatory_source = "NONE";
    if (!atac_cov_file.empty()) {
        if (!atac_regulatory_bed_file.empty()) {
            load_bed(atac_regulatory_bed_file, atac_regulatory_regions,
                     "ATAC regulatory");
            atac_regulatory_source = atac_regulatory_bed_file;
        } else if (!gtf_file.empty()) {
            load_gtf_promoters(gtf_file, atac_regulatory_regions, 2000);
            atac_regulatory_source = "GTF_2KB_UPSTREAM_PROMOTERS";
        }
        fprintf(stderr,
            "ATAC regulatory prior: %s (%s; required=%s)\n",
            atac_regulatory_source.c_str(),
            atac_regulatory_regions.empty() ? "no intervals" : "ranking bonus",
            atac_regulatory_required ? "yes" : "no");
        if (atac_regulatory_required && atac_regulatory_regions.empty()) {
            fprintf(stderr,
                "ERROR: --atac_regulatory_required was set but the regulatory prior is empty\n");
            exit(1);
        }
    }

    // Load coverage into memory for fast lookup
    map<string, vector<CoverageInterval>> coverage_map;
    if (!cov_file.empty()) {
        fprintf(stderr, "Loading RNA coverage from %s...\n", cov_file.c_str());
        coverage_map = load_coverage(cov_file);
    }
    
    // Load ATAC coverage if requested
    map<string, vector<CoverageInterval>> atac_coverage_map;
    if (!atac_cov_file.empty()) {
        fprintf(stderr, "Loading ATAC coverage from %s...\n", atac_cov_file.c_str());
        atac_coverage_map = load_coverage(atac_cov_file);
    }

    // ========================================================================
    // ANNOTATE_ONLY MODE: Stream through VCF, add scores to all sites, exit
    // ========================================================================
    if (annotate_only) {
        fprintf(stderr, "\n========================================\n");
        fprintf(stderr, "ANNOTATE ONLY MODE\n");
        fprintf(stderr, "========================================\n");
        
        auto start_time = chrono::high_resolution_clock::now();
        
        // Open input VCF
        htsFile* in_fp = bcf_open(vcf_file.c_str(), "r");
        if (!in_fp) {
            fprintf(stderr, "ERROR: Could not open %s\n", vcf_file.c_str());
            exit(1);
        }
        hts_set_threads(in_fp, n_threads);
        
        bcf_hdr_t* hdr = bcf_hdr_read(in_fp);
        if (!hdr) {
            fprintf(stderr, "ERROR: Could not read VCF header\n");
            exit(1);
        }

        if (!blacklist_file.empty()) {
            map<string, int> header_scaffolds;
            int n_header_scaffolds = 0;
            const char** header_names = bcf_hdr_seqnames(
                hdr, &n_header_scaffolds);
            for (int i = 0; i < n_header_scaffolds; ++i) {
                header_scaffolds[header_names[i]] = i;
            }
            free(header_names);
            validate_blacklist_scaffolds(
                blacklist, header_scaffolds, blacklist_file);
        }
        
        int n_samples = bcf_hdr_nsamples(hdr);
        fprintf(stderr, "Samples: %d\n", n_samples);
        
        // Add INFO headers
        bcf_hdr_append(hdr, "##INFO=<ID=HET_FREQ,Number=1,Type=Float,Description=\"Fraction of called individuals that are heterozygous\">");
        bcf_hdr_append(hdr, "##INFO=<ID=MISSING_FREQ,Number=1,Type=Float,Description=\"Fraction of individuals with missing genotype\">");
        bcf_hdr_append(hdr, "##INFO=<ID=HET_SCORE,Number=1,Type=Float,Description=\"Combined het score: het_freq * (1 - missing_freq * 0.5)\">");
        bcf_hdr_append(hdr, "##INFO=<ID=N_HET,Number=1,Type=Integer,Description=\"Number of heterozygous individuals\">");
        bcf_hdr_append(hdr, "##INFO=<ID=N_CALLED,Number=1,Type=Integer,Description=\"Number of individuals with called genotype\">");
        bcf_hdr_append(hdr, "##INFO=<ID=HET_SEL_SCORE,Number=1,Type=Float,Description=\"Selection ranking score (annotation-only estimate)\">");
        if (!annotations.empty()) {
            bcf_hdr_append(hdr, "##INFO=<ID=ANNOT_SCORE,Number=1,Type=Float,Description=\"Annotation weight (1.0=genic, 0.1=intergenic)\">");
        }
        if (!coverage_map.empty()) {
            bcf_hdr_append(hdr, "##INFO=<ID=COV_SCORE,Number=1,Type=Float,Description=\"Coverage score from bedgraph\">");
        }
        bcf_hdr_append(hdr, "##INFO=<ID=BIN_ID,Number=1,Type=String,Description=\"Chromosome:bin_index for evenness tracking\">");
        bcf_hdr_append(hdr, "##INFO=<ID=NUCLEAR_PANEL_STATUS,Number=1,Type=String,Description=\"PASS or exclusion reason(s): MITO, NUMT, BLACKLIST\">");
        
        if (bcf_hdr_sync(hdr) < 0) {
            fprintf(stderr, "ERROR: Failed to sync header\n");
            exit(1);
        }
        
        // Open output BCF
        htsFile* out_fp = hts_open(outfile.c_str(), "wb");
        if (!out_fp) {
            fprintf(stderr, "ERROR: Could not open output file %s\n", outfile.c_str());
            exit(1);
        }
        hts_set_threads(out_fp, n_threads);
        
        if (bcf_hdr_write(out_fp, hdr) < 0) {
            fprintf(stderr, "ERROR: Failed to write header\n");
            exit(1);
        }
        
        // Process records
        bcf1_t* rec = bcf_init();
        int32_t* gts = NULL;
        int n_gts = 0;
        
        long n_records = 0;
        long n_with_het = 0;
        
        fprintf(stderr, "Processing variants...\n");
        
        while (bcf_read(in_fp, hdr, rec) >= 0) {
            
            // Get genotypes
            int ret = bcf_get_genotypes(hdr, rec, &gts, &n_gts);
            
            int n_het = 0;
            int n_called = 0;
            int n_missing = 0;
            
            if (ret > 0) {
                int ploidy = ret / n_samples;
                
                for (int i = 0; i < n_samples; i++) {
                    int32_t* gt_ptr = gts + i * ploidy;
                    
                    // Check for missing
                    if (bcf_gt_is_missing(gt_ptr[0]) || 
                        (ploidy > 1 && bcf_gt_is_missing(gt_ptr[1]))) {
                        n_missing++;
                        continue;
                    }
                    
                    n_called++;
                    
                    // Check for het (different alleles)
                    if (ploidy >= 2) {
                        int a1 = bcf_gt_allele(gt_ptr[0]);
                        int a2 = bcf_gt_allele(gt_ptr[1]);
                        if (a1 != a2) {
                            n_het++;
                        }
                    }
                }
            } else {
                // No genotype data - all missing
                n_missing = n_samples;
            }
            
            // Compute metrics
            float het_freq = (n_called > 0) ? (float)n_het / n_called : 0.0f;
            float missing_freq = (float)n_missing / n_samples;
            float missing_penalty = 1.0f - (missing_freq * 0.5f);
            float het_score = het_freq * missing_penalty;
            
            // Add INFO fields
            bcf_update_info_float(hdr, rec, "HET_FREQ", &het_freq, 1);
            bcf_update_info_float(hdr, rec, "MISSING_FREQ", &missing_freq, 1);
            bcf_update_info_float(hdr, rec, "HET_SCORE", &het_score, 1);
            bcf_update_info_int32(hdr, rec, "N_HET", &n_het, 1);
            bcf_update_info_int32(hdr, rec, "N_CALLED", &n_called, 1);
            
            // C5: HET_SEL_SCORE parity with main pipeline
            {
                float call_rate_ao = 1.0f - missing_freq;
                float cov_ao = 0.0f;
                if (!coverage_map.empty()) {
                    const char* chrom_ao = bcf_hdr_id2name(hdr, rec->rid);
                    cov_ao = get_coverage_score_fast(chrom_ao, rec->pos, coverage_map);
                }
                float annot_ao = 0.1f;
                if (!annotations.empty()) {
                    const char* chrom_ao2 = bcf_hdr_id2name(hdr, rec->rid);
                    annot_ao = get_annot_score(chrom_ao2, rec->pos + 1, annotations);
                }
                float annot_boost_ao = (annot_ao >= 1.0f) ? 1.0f : 0.0f;
                float sel_score_ao = het_freq * call_rate_ao * (log2f(cov_ao + 1.0f) + 0.1f * annot_boost_ao);
                bcf_update_info_float(hdr, rec, "HET_SEL_SCORE", &sel_score_ao, 1);
            }
            
            // Add ANNOT_SCORE if GTF provided
            if (!annotations.empty()) {
                const char* chrom = bcf_hdr_id2name(hdr, rec->rid);
                float annot_score = get_annot_score(chrom, rec->pos + 1, annotations);
                bcf_update_info_float(hdr, rec, "ANNOT_SCORE", &annot_score, 1);
            }
            
            // Add COV_SCORE if coverage provided
            if (!coverage_map.empty()) {
                const char* chrom = bcf_hdr_id2name(hdr, rec->rid);
                float cov_score = get_coverage_score_fast(chrom, rec->pos, coverage_map);
                bcf_update_info_float(hdr, rec, "COV_SCORE", &cov_score, 1);
            }
            
            // Add BIN_ID and nuclear-panel exclusion status.
            const char* chrom = bcf_hdr_id2name(hdr, rec->rid);
            int bin_idx_val = get_bin_idx(rec->pos, bin_size);
            string bin_id = make_bin_id(chrom, bin_idx_val);
            bcf_update_info_string(hdr, rec, "BIN_ID", bin_id.c_str());
            uint8_t exclusion = nuclear_exclusion_mask(
                chrom, rec->pos, mito_chrom, include_mito,
                numts, include_numts, blacklist);
            string exclusion_reason = nuclear_exclusion_reason(exclusion);
            bcf_update_info_string(hdr, rec, "NUCLEAR_PANEL_STATUS", exclusion_reason.c_str());
            
            // Write record
            if (bcf_write(out_fp, hdr, rec) < 0) {
                fprintf(stderr, "ERROR: Failed to write record\n");
                exit(1);
            }
            
            n_records++;
            if (n_het > 0) n_with_het++;
            
            if (n_records % 1000000 == 0) {
                fprintf(stderr, "  Processed %ld million variants...\r", n_records / 1000000);
                fflush(stderr);
            }
        }
        
        fprintf(stderr, "\n");
        
        // Cleanup
        free(gts);
        bcf_destroy(rec);
        bcf_hdr_destroy(hdr);
        hts_close(in_fp);
        hts_close(out_fp);
        
        // Index output BCF
        fprintf(stderr, "Indexing output BCF...\n");
        if (bcf_index_build(outfile.c_str(), 14) < 0) {
            fprintf(stderr, "WARNING: Failed to create index for %s\n", outfile.c_str());
        }
        
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
        
        fprintf(stderr, "\n========================================\n");
        fprintf(stderr, "ANNOTATION COMPLETE\n");
        fprintf(stderr, "========================================\n");
        fprintf(stderr, "  Total variants: %ld\n", n_records);
        fprintf(stderr, "  With het individuals: %ld (%.1f%%)\n", 
                n_with_het, 100.0 * n_with_het / n_records);
        fprintf(stderr, "  Time: %ld seconds\n", duration.count());
        fprintf(stderr, "  Output: %s\n", outfile.c_str());
        
        return 0;
    }

    // Note: --annotate_all is now handled in Pass 2 (parallel, with all scores including DEMUX_SCORE)

    // ========================================================================
    // Get chromosome list from VCF
    // ========================================================================
    htsFile* bcf_reader = bcf_open(vcf_file.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR interpreting %s as BCF format.\n", vcf_file.c_str());
        exit(1);
    }
    bcf_hdr_t* bcf_header = bcf_hdr_read(bcf_reader);
    int num_samples = bcf_hdr_nsamples(bcf_header);
    
    vector<string> samples;
    for (int i = 0; i < num_samples; ++i){
        samples.push_back(bcf_header->samples[i]);
    }
    
    if (samples.size() * 4 > NBITS){
        fprintf(stderr, "ERROR: too many samples in VCF. Please recompile with NBITS = %ld or higher.\n",
            samples.size() * 4);
        exit(1); 
    }

    // Get chromosome names
    vector<string> chroms;
    map<string, int> chrom_to_idx;  // For ordering chromosomes
    int n_seqs = 0;
    const char** seqnames = bcf_hdr_seqnames(bcf_header, &n_seqs);
    for (int i = 0; i < n_seqs; ++i){
        chrom_to_idx[seqnames[i]] = chroms.size();
        chroms.push_back(seqnames[i]);
    }
    free(seqnames);

    if (!blacklist_file.empty()) {
        validate_blacklist_scaffolds(
            blacklist, chrom_to_idx, blacklist_file);
    }

    vector<uint64_t> atac_excluded_positions;
    if (!atac_exclude_vcfs.empty()) {
        atac_excluded_positions = load_panel_exclusion_positions(
            atac_exclude_vcfs, chrom_to_idx, n_threads);
        if (atac_excluded_positions.empty()) {
            fprintf(stderr,
                "ERROR: --atac_exclude_vcf was supplied but its union has no loci "
                "on source-VCF contigs\n");
            exit(1);
        }
    }
    vector<uint64_t> atac_soft_overlap_positions;
    if (!atac_soft_overlap_vcfs.empty()) {
        atac_soft_overlap_positions = load_panel_exclusion_positions(
            atac_soft_overlap_vcfs, chrom_to_idx, n_threads);
        if (atac_soft_overlap_positions.empty()) {
            fprintf(stderr,
                "ERROR: --atac_soft_overlap_vcf was supplied but has no loci "
                "on source-VCF contigs\n");
            exit(1);
        }
        fprintf(stderr,
            "ATAC soft-overlap policy: %zu RNA-demux loci remain eligible "
            "with ranking multiplier %.6f\n",
            atac_soft_overlap_positions.size(),
            atac_soft_overlap_multiplier);
    }
    
    bcf_hdr_destroy(bcf_header);
    hts_close(bcf_reader);

    fprintf(stderr, "Found %d chromosomes, %d samples\n", (int)chroms.size(), num_samples);
    
    // Build sample name -> index map for pool/species lookups
    map<string, int> sample_name_to_idx;
    for (int i = 0; i < num_samples; i++) {
        sample_name_to_idx[samples[i]] = i;
    }
    
    // Load panel metadata if needed
    PanelMetadata panel_metadata;
    bool has_panel_metadata = false;
    vector<int> species_min_called_per_species;
    if (!panel_metadata_file.empty()) {
        // The species selector treats the Chinobo hybrid as its own explicit
        // class.  Folding it into both parental species is appropriate for
        // some mixture models, but conflicts with the species-consistency
        // gate and destroys bonobo/chimp discrimination here.
        panel_metadata = load_panel_metadata(panel_metadata_file, samples, false);
        has_panel_metadata = true;
        
        // CHANGE 1: Compute per-species minimum-called thresholds.
        // For singleton species (n_in_sp == 1), require min_called=1 instead of
        // the previous max(2, ceil(0.5 * n_in_sp)) which made singletons impossible.
        int n_sp_local = (int)panel_metadata.species_list.size();
        species_min_called_per_species.resize(n_sp_local, 1);
        for (int si = 0; si < n_sp_local; si++) {
            const string& sp = panel_metadata.species_list[si];
            auto it = panel_metadata.species_to_sample_indices.find(sp);
            int n_in_sp = (it != panel_metadata.species_to_sample_indices.end())
                          ? (int)it->second.size() : 0;
            int per_species_min;
            if (n_in_sp <= 1) {
                per_species_min = 1;
                fprintf(stderr, "  WARNING: species %s has only %d individual(s); "
                    "--species_min_called is forced to 1 for this species "
                    "(low-confidence diagnostic markers).\n",
                    sp.c_str(), n_in_sp);
            } else {
                per_species_min = max(1, (int)ceil(0.5 * n_in_sp));
                // The user-supplied --species_min_called acts as an additional floor,
                // but is capped at n_in_sp so it cannot make the species impossible.
                int user_floor = min(n_in_sp, species_min_called);
                if (user_floor > per_species_min) per_species_min = user_floor;
            }
            species_min_called_per_species[si] = per_species_min;
            fprintf(stderr, "  Species %s: n_in_panel=%d, min_called_required=%d\n",
                sp.c_str(), n_in_sp, per_species_min);
        }
    }
    
    // Load pool combinations if needed
    PoolCombinations pool_combos;
    bool has_pool_combos = false;
    if (!pool_combinations_file.empty()) {
        pool_combos = load_pool_combinations(pool_combinations_file, samples);
        has_pool_combos = true;
    }
    vector<int> atac_roster_sample_indices;
    vector<PoolDonorPair> atac_pool_donor_pairs;
    if (has_pool_combos) {
        atac_pool_donor_pairs = build_pool_donor_pairs(
            pool_combos, sample_name_to_idx, samples,
            atac_roster_sample_indices);
        fprintf(stderr,
            "ATAC operational roster: %zu donors, %zu unique co-pool donor pairs\n",
            atac_roster_sample_indices.size(), atac_pool_donor_pairs.size());
        fprintf(stderr,
            "ATAC pool policy: RNA-production scalar pool-discrimination "
            "ranking; no legal-identity pair floors or rescue objectives\n");
    } else {
        atac_roster_sample_indices.resize(num_samples);
        for (int i = 0; i < num_samples; ++i) atac_roster_sample_indices[i] = i;
    }
    if (has_panel_metadata && has_pool_combos) {
        const set<int> roster_set(atac_roster_sample_indices.begin(),
                                  atac_roster_sample_indices.end());
        vector<string> retained_species;
        for (const string& species : panel_metadata.species_list) {
            vector<int>& indices =
                panel_metadata.species_to_sample_indices[species];
            indices.erase(remove_if(indices.begin(), indices.end(),
                [&](int sample_index) {
                    return roster_set.count(sample_index) == 0;
                }), indices.end());
            if (!indices.empty()) retained_species.push_back(species);
        }
        panel_metadata.species_list.swap(retained_species);
        panel_metadata.species_pairs.clear();
        panel_metadata.pair_to_index.clear();
        panel_metadata.n_pairs = 0;
        for (size_t i = 0; i < panel_metadata.species_list.size(); ++i) {
            for (size_t j = i + 1;
                 j < panel_metadata.species_list.size(); ++j) {
                pair<string, string> species_pair = make_pair(
                    panel_metadata.species_list[i],
                    panel_metadata.species_list[j]);
                panel_metadata.pair_to_index[species_pair] =
                    panel_metadata.n_pairs++;
                panel_metadata.species_pairs.push_back(species_pair);
            }
        }
        species_min_called_per_species.assign(
            panel_metadata.species_list.size(), 1);
        for (size_t si = 0; si < panel_metadata.species_list.size(); ++si) {
            const string& species = panel_metadata.species_list[si];
            const int n_in_species = (int)
                panel_metadata.species_to_sample_indices[species].size();
            if (n_in_species <= 1) {
                species_min_called_per_species[si] = 1;
            } else {
                species_min_called_per_species[si] = max(
                    max(1, (int)ceil(0.5 * n_in_species)),
                    min(n_in_species, species_min_called));
            }
            fprintf(stderr,
                "  Operational species %s: n_in_pool_roster=%d, "
                "min_called_required=%d\n",
                species.c_str(), n_in_species,
                species_min_called_per_species[si]);
        }
        fprintf(stderr,
            "Operational species selector: %zu classes, %d pairs "
            "(pool-roster intersection; hybrid not folded)\n",
            panel_metadata.species_list.size(), panel_metadata.n_pairs);
    }

    // ========================================================================
    // Load VCF index (CSI required for parallel region queries)
    // ========================================================================

    fprintf(stderr, "Loading VCF index...\n");
    hts_idx_t* vcf_idx = bcf_index_load(vcf_file.c_str());
    if (!vcf_idx) {
        fprintf(stderr, "ERROR: Could not load VCF index (.csi required for parallel queries)\n");
        fprintf(stderr, "  Create with: bcftools index -c %s\n", vcf_file.c_str());
        exit(1);
    }
    fprintf(stderr, "  Index loaded successfully\n");

    // ========================================================================
    // LOAD ENTIRE VCF INTO MEMORY
    // ========================================================================
    
    fprintf(stderr, "\n========================================\n");
    fprintf(stderr, "Loading VCF into memory...\n");
    fprintf(stderr, "========================================\n");
    
    auto load_start = chrono::high_resolution_clock::now();
    
    htsFile* load_fp = bcf_open(vcf_file.c_str(), "r");
    if (!load_fp) {
        fprintf(stderr, "ERROR: Could not open %s for loading\n", vcf_file.c_str());
        exit(1);
    }
    hts_set_threads(load_fp, n_threads);
    
    bcf_hdr_t* load_hdr = bcf_hdr_read(load_fp);
    bcf1_t* load_rec = bcf_init();
    int32_t* load_gts = NULL;
    int load_n_gts = 0;
    
    vector<StoredVariant> all_variants;
    deque<VariantPattern> variant_patterns;
    unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, const VariantPattern*>
        interned_patterns;
    // chrM is tiny (~2K source records). Preserve only those raw records so the
    // dedicated mitochondrial panel can be emitted from this same full VCF load
    // with RO/AO/AD intact, rather than reopening and rescanning the 73-GB panel.
    vector<bcf1_t*> mt_source_records;
    if (!mt_outfile.empty()) mt_source_records.reserve(4096);
    // Reserve from CSI/BCF index statistics instead of unconditionally reserving
    // space for 150M records. This preserves the no-reallocation fast path on
    // production panels while avoiding an enormous virtual allocation for small
    // tests or targeted VCFs.
    uint64_t indexed_variant_count = 0;
    for (int rid = 0; rid < (int)chroms.size(); ++rid) {
        uint64_t mapped = 0, unmapped = 0;
        if (hts_idx_get_stat(vcf_idx, rid, &mapped, &unmapped) == 0) {
            indexed_variant_count += mapped + unmapped;
        }
    }
    const size_t reserve_count = indexed_variant_count > 0
        ? (size_t)min<uint64_t>(indexed_variant_count + indexed_variant_count / 100 + 1024,
                                150000000ULL)
        : 5000000ULL;
    all_variants.reserve(reserve_count);
    interned_patterns.reserve(min<size_t>(reserve_count, 2000000ULL));
    fprintf(stderr, "  Reserving in-memory capacity for %zu variants (index estimate: %llu)\n",
        reserve_count, (unsigned long long)indexed_variant_count);
    
    long n_loaded = 0;
    long n_biallelic = 0;
    long n_pass_qc = 0;
    long n_excluded_mito = 0;
    long n_excluded_numt = 0;
    long n_excluded_blacklist = 0;
    
    while (bcf_read(load_fp, load_hdr, load_rec) >= 0) {
        StoredVariant sv;
        sv.chrom_idx = load_rec->rid;
        sv.pos = load_rec->pos;
        sv.n_alleles = load_rec->n_allele;
        sv.is_biallelic = (load_rec->n_allele == 2);
        sv.passes_qc = false;
        sv.pattern = NULL;
        const string& load_chrom = chroms[sv.chrom_idx];
        sv.nuclear_exclusion = nuclear_exclusion_mask(
            load_chrom, sv.pos, mito_chrom, include_mito,
            numts, include_numts, blacklist);
        sv.qual = load_rec->qual;
        sv.rna_cov_score = coverage_map.empty()
            ? 0.0f : get_coverage_score_fast(load_chrom, sv.pos, coverage_map);
        sv.atac_cov_score = atac_coverage_map.empty()
            ? 0.0f : get_coverage_score_fast(load_chrom, sv.pos, atac_coverage_map);
        sv.annot_score = annotations.empty()
            ? 0.1f : get_annot_score(load_chrom, sv.pos + 1, annotations);
        sv.atac_regulatory = !atac_regulatory_regions.empty() &&
            is_blacklisted(load_chrom, sv.pos, atac_regulatory_regions);
        sv.atac_soft_overlap = panel_position_is_excluded(
            atac_soft_overlap_positions, sv.chrom_idx, sv.pos);
        
        // Get alleles - store ALL alleles for ALL variants
        bcf_unpack(load_rec, BCF_UN_STR | BCF_UN_FMT);
        if (!mt_outfile.empty() && load_chrom == mito_chrom) {
            bcf1_t* mt_copy = bcf_dup(load_rec);
            if (!mt_copy) {
                fprintf(stderr, "ERROR: failed to preserve chrM record for mitochondrial panel\n");
                exit(1);
            }
            mt_source_records.push_back(mt_copy);
        }
        if (load_rec->n_allele > 0) {
            sv.alleles = load_rec->d.allele[0];
            for (int i = 1; i < load_rec->n_allele; i++) {
                sv.alleles += ",";
                sv.alleles += load_rec->d.allele[i];
            }
        }
        if (sv.is_biallelic) {
            n_biallelic++;
        }
        
        // Get genotypes
        int num_loaded_gt = bcf_get_genotypes(load_hdr, load_rec, &load_gts, &load_n_gts);
        if (num_loaded_gt > 0) {
            // Pack raw genotypes into 2-bit-per-allele storage
            pack_genotypes(load_gts, num_samples, sv.genotypes_packed);
            
            // For non-biallelic records, also keep raw genotypes so annotate_all
            // can write them without the biallelic-collapsing data loss (Fix 3).
            if (!sv.is_biallelic && !annotate_all_file.empty()) {
                sv.genotypes_raw.assign(load_gts, load_gts + num_loaded_gt);
            }
            
            // For biallelic variants, precompute dosage-aware bitsets
            if (sv.is_biallelic) {
                bool pass = true;
                
                // Check alleles are valid A/C/G/T
                for (int i = 0; i < load_rec->n_allele; ++i) {
                    if (strcmp(load_rec->d.allele[i], "A") != 0 &&
                        strcmp(load_rec->d.allele[i], "C") != 0 &&
                        strcmp(load_rec->d.allele[i], "G") != 0 && 
                        strcmp(load_rec->d.allele[i], "T") != 0) {
                        pass = false;
                        break;
                    }
                }
                // Check ref != alt
                if (load_rec->d.allele[0][0] == load_rec->d.allele[1][0]) {
                    pass = false;
                }
                
                if (pass) {
                    // Dosage-aware bitset encoding (A1/A5):
                    // bit 2i     = 1 iff sample i is called AND carries >= 1 alt allele
                    // bit 2i + 1 = 1 iff sample i is called AND is homozygous alt
                    VariantPattern pattern;
                    pattern.alt.reset();
                    pattern.alt_flip.reset();
                    pattern.present.reset();
                    
                    bool dosage_seen[3] = {false, false, false};
                    int total_alt_val = 0;
                    
                    for (int i = 0; i < num_samples; ++i) {
                        bool miss_a = gt_packed_is_missing(sv.genotypes_packed, i, 0);
                        bool miss_b = gt_packed_is_missing(sv.genotypes_packed, i, 1);
                        
                        if (!miss_a && !miss_b) {
                            pattern.present.set(i);
                            int a1 = gt_packed_allele(sv.genotypes_packed, i, 0);
                            int a2 = gt_packed_allele(sv.genotypes_packed, i, 1);
                            int dosage = (a1 != 0 ? 1 : 0) + (a2 != 0 ? 1 : 0);
                            
                            if (dosage >= 1) pattern.alt.set(2 * i);
                            if (dosage == 2) pattern.alt.set(2 * i + 1);
                            
                            total_alt_val += dosage;
                            dosage_seen[dosage] = true;
                        }
                    }
                    
                    // Build alt_flip using dosage_flip helper
                    dosage_flip(pattern.alt, pattern.present, num_samples,
                                pattern.alt_flip);

                    const pair<bitset<NBITS>, bitset<NBITS>> pattern_key =
                        make_pair(pattern.present, pattern.alt);
                    auto pattern_it = interned_patterns.find(pattern_key);
                    if (pattern_it == interned_patterns.end()) {
                        variant_patterns.push_back(std::move(pattern));
                        sv.pattern = &variant_patterns.back();
                        interned_patterns.emplace(pattern_key, sv.pattern);
                    } else {
                        sv.pattern = pattern_it->second;
                    }
                    
                    // Store total_alt and total_alt_flip (A6)
                    int n_called = (int)sv.pattern->present.count();
                    sv.total_alt = (int16_t)total_alt_val;
                    sv.total_alt_flip = (int16_t)(n_called * 2 - total_alt_val);
                    
                    const int distinct_dosages =
                        (dosage_seen[0] ? 1 : 0) +
                        (dosage_seen[1] ? 1 : 0) +
                        (dosage_seen[2] ? 1 : 0);
                    if (distinct_dosages >= 2) {
                        sv.passes_qc = true;
                        n_pass_qc++;
                    }
                }
            }
        }
        
        if (sv.nuclear_exclusion & NUCLEAR_EXCLUDE_MITO) n_excluded_mito++;
        if (sv.nuclear_exclusion & NUCLEAR_EXCLUDE_NUMT) n_excluded_numt++;
        if (sv.nuclear_exclusion & NUCLEAR_EXCLUDE_BLACKLIST) n_excluded_blacklist++;

        sv.variant_idx = n_loaded;  // Store index for back-reference
        
        all_variants.push_back(std::move(sv));
        n_loaded++;
        
        if (n_loaded % 5000000 == 0) {
            fprintf(stderr, "  Loaded %ld variants (%ld biallelic, %ld pass QC)...\r", 
                n_loaded, n_biallelic, n_pass_qc);
            fflush(stderr);
        }
    }
    
    if (!mt_outfile.empty()) {
        mtpanel::Config mt_config;
        mt_config.output = mt_outfile;
        mt_config.pool_combinations = pool_combinations_file;
        mt_config.mito_chrom = mito_chrom;
        // ATAC-only builds intentionally do not load/pass a duplicate generic
        // coverage track.  The mitochondrial builder consumes the identical
        // ATAC aggregate directly in that mode.
        mt_config.coverage = !cov_file.empty() ? cov_file : atac_cov_file;
        mt_config.blacklist_bed = blacklist_file;
        mt_config.mt_mask_bed = mt_mask_bed_file;
        mt_config.pair_audit = mt_pair_audit_file;
        mt_config.site_manifest = mt_site_manifest_file;
        mt_config.haplotype_groups = mt_haplotype_groups_file;
        mt_config.haplotype_pairwise = mt_haplotype_pairwise_file;
        mt_config.sample_audit = mt_sample_audit_file;
        mt_config.sites_bed = mt_sites_bed_file;
        mt_config.threads = n_threads;
        mt_config.min_depth = mt_min_depth;
        mt_config.min_pair_sites = mt_min_pair_sites;
        mt_config.min_ambient_sites = mt_min_ambient_sites;
        mt_config.homoplasmy_af = mt_homoplasmy_af;
        mt_config.min_coverage = mt_min_coverage;
        mt_config.require_pair_targets = mt_require_pair_targets;
        mt_config.require_coverage = mt_require_coverage;
        try {
            mtpanel::build_panel(load_hdr, mt_source_records, mt_config);
        } catch (const std::exception& error) {
            for (size_t i = 0; i < mt_source_records.size(); ++i) {
                bcf_destroy(mt_source_records[i]);
            }
            free(load_gts);
            bcf_destroy(load_rec);
            bcf_hdr_destroy(load_hdr);
            hts_close(load_fp);
            fprintf(stderr, "ERROR: mitochondrial panel build failed: %s\n", error.what());
            exit(1);
        }
    }
    for (size_t i = 0; i < mt_source_records.size(); ++i) {
        bcf_destroy(mt_source_records[i]);
    }
    mt_source_records.clear();

    free(load_gts);
    bcf_destroy(load_rec);
    bcf_hdr_destroy(load_hdr);
    hts_close(load_fp);
    
    // all_variants.shrink_to_fit() removed: at 150M variants the reallocation
    // copy spike can exceed available memory. The reserve(150M) overshoot is
    // a smaller cost than a transient 2x peak.
    
    auto load_end = chrono::high_resolution_clock::now();
    auto load_duration = chrono::duration_cast<chrono::seconds>(load_end - load_start);
    
    fprintf(stderr, "\n  Loaded %ld variants in %ld seconds\n", n_loaded, load_duration.count());
    fprintf(stderr, "  Biallelic: %ld (%.1f%%)\n", n_biallelic, 100.0 * n_biallelic / n_loaded);
    fprintf(stderr, "  Pass QC: %ld (%.1f%%)\n", n_pass_qc, 100.0 * n_pass_qc / n_loaded);
    fprintf(stderr, "  Interned genotype patterns: %zu (shared across %ld records)\n",
        variant_patterns.size(), n_loaded);
    fprintf(stderr, "  Excluded from nuclear panels: chrM=%ld, NUMT=%ld, blacklist=%ld\n",
        n_excluded_mito, n_excluded_numt, n_excluded_blacklist);
    if (!blacklist_file.empty()) {
        long blacklist_rna_demux_candidates = 0;
        long blacklist_rna_het_candidates = 0;
        long blacklist_species_candidates = 0;
        long blacklist_atac_demux_candidates = 0;
        long blacklist_atac_het_candidates = 0;
        for (const StoredVariant& sv : all_variants) {
            if (!(sv.nuclear_exclusion & NUCLEAR_EXCLUDE_BLACKLIST) ||
                !sv.is_biallelic || sv.pattern == NULL) {
                continue;
            }
            if (!outfile.empty() && sv.passes_qc &&
                sv.rna_cov_score >= min_cov) {
                blacklist_rna_demux_candidates++;
            }
            if (!het_outfile.empty() && !sv.genotypes_packed.empty() &&
                sv.rna_cov_score >= min_cov) {
                blacklist_rna_het_candidates++;
            }
            if (!species_outfile.empty() && !sv.genotypes_packed.empty()) {
                const float species_cov = species_coverage_mode == "atac"
                    ? sv.atac_cov_score : sv.rna_cov_score;
                const float species_min_cov = species_coverage_mode == "atac"
                    ? atac_min_cov : min_cov;
                if (species_cov >= species_min_cov) {
                    blacklist_species_candidates++;
                }
            }
            if (!atac_outfile.empty() && sv.passes_qc &&
                sv.atac_cov_score >= atac_min_cov) {
                blacklist_atac_demux_candidates++;
            }
            if (!atac_het_outfile.empty() && !sv.genotypes_packed.empty() &&
                sv.atac_cov_score >= atac_min_cov) {
                blacklist_atac_het_candidates++;
            }
        }
        if (!outfile.empty()) {
            fprintf(stderr,
                "BLACKLIST_EXCLUDED_CANDIDATES panel=RNA_DEMUX count=%ld\n",
                blacklist_rna_demux_candidates);
        }
        if (!het_outfile.empty()) {
            fprintf(stderr,
                "BLACKLIST_EXCLUDED_CANDIDATES panel=RNA_HET count=%ld\n",
                blacklist_rna_het_candidates);
        }
        if (!species_outfile.empty()) {
            fprintf(stderr,
                "BLACKLIST_EXCLUDED_CANDIDATES panel=SPECIES count=%ld\n",
                blacklist_species_candidates);
        }
        if (!atac_outfile.empty()) {
            fprintf(stderr,
                "BLACKLIST_EXCLUDED_CANDIDATES panel=ATAC_DEMUX count=%ld\n",
                blacklist_atac_demux_candidates);
        }
        if (!atac_het_outfile.empty()) {
            fprintf(stderr,
                "BLACKLIST_EXCLUDED_CANDIDATES panel=ATAC_HET count=%ld\n",
                blacklist_atac_het_candidates);
        }
    }
    
    // Build chromosome -> variant range index
    vector<pair<size_t, size_t>> chrom_ranges(chroms.size(), {0, 0});
    {
        size_t start_idx = 0;
        int current_chrom = -1;
        for (size_t i = 0; i < all_variants.size(); i++) {
            if (all_variants[i].chrom_idx != current_chrom) {
                if (current_chrom >= 0 && current_chrom < (int)chroms.size()) {
                    chrom_ranges[current_chrom] = {start_idx, i};
                }
                current_chrom = all_variants[i].chrom_idx;
                start_idx = i;
            }
        }
        if (current_chrom >= 0 && current_chrom < (int)chroms.size()) {
            chrom_ranges[current_chrom] = {start_idx, all_variants.size()};
        }
    }

    // ========================================================================
    // REUSABLE ATAC CANDIDATE POOL
    //
    // This is deliberately a union input pool, not merely the stricter
    // demux subset.  Keeping every ATAC-covered,
    // nuclear-safe biallelic SNV with usable genotypes preserves the inputs
    // needed to re-run demux, het, and species selection without reloading the
    // much larger source VCF.
    // ========================================================================

    long atac_candidates_written = 0;
    long atac_demux_eligible_written = 0;
    struct CandidateAuditPatternCount {
        const VariantPattern* pattern;
        size_t representative_variant_idx;
        long count;
        long soft_overlap_count;
    };
    vector<CandidateAuditPatternCount> candidate_audit_union_patterns;
    vector<CandidateAuditPatternCount> candidate_audit_demux_patterns;
    unordered_map<const VariantPattern*, size_t>
        candidate_audit_union_pattern_index;
    unordered_map<const VariantPattern*, size_t>
        candidate_audit_demux_pattern_index;
    const bool collect_candidate_only_audits = atac_candidates_only &&
        !atac_audit_prefix.empty();
    const int candidate_audit_n_gates = 7;
    long candidate_audit_gate_counts[candidate_audit_n_gates] =
        {0,0,0,0,0,0,0};
    long candidate_audit_gate_overlap[candidate_audit_n_gates] =
        {0,0,0,0,0,0,0};
    if (!atac_candidates_outfile.empty()) {
        fprintf(stderr, "\n========================================\n");
        fprintf(stderr, "Writing reusable ATAC candidate pool...\n");
        fprintf(stderr, "========================================\n");

        htsFile* candidate_header_reader = hts_open(vcf_file.c_str(), "r");
        if (!candidate_header_reader) {
            fprintf(stderr, "ERROR: Failed to reopen source header for ATAC candidates\n");
            exit(1);
        }
        bcf_hdr_t* candidate_header = bcf_hdr_read(candidate_header_reader);
        hts_close(candidate_header_reader);
        if (!candidate_header) {
            fprintf(stderr, "ERROR: Failed to read source header for ATAC candidates\n");
            exit(1);
        }

        bcf_hdr_append(candidate_header,
            "##INFO=<ID=ATAC_COV_SCORE,Number=1,Type=Float,Description=\"ATAC coverage score used for candidate eligibility\">");
        bcf_hdr_append(candidate_header,
            "##INFO=<ID=ATAC_BIN_ID,Number=1,Type=String,Description=\"Chromosome:bin_index for the ATAC candidate pool\">");
        bcf_hdr_append(candidate_header,
            "##INFO=<ID=ATAC_N_CALLED,Number=1,Type=Integer,Description=\"Number of panel individuals with called diploid genotypes\">");
        bcf_hdr_append(candidate_header,
            "##INFO=<ID=ATAC_N_HET,Number=1,Type=Integer,Description=\"Number of called panel individuals with heterozygous genotypes\">");
        bcf_hdr_append(candidate_header,
            "##INFO=<ID=ATAC_HET_FREQ,Number=1,Type=Float,Description=\"Fraction of called panel individuals that are heterozygous\">");
        bcf_hdr_append(candidate_header,
            "##INFO=<ID=ATAC_MISSING_FREQ,Number=1,Type=Float,Description=\"Fraction of panel individuals without a complete diploid genotype\">");
        bcf_hdr_append(candidate_header,
            "##INFO=<ID=ATAC_DEMUX_ELIGIBLE,Number=1,Type=Integer,Description=\"1 when the site is polymorphic and clears the operational-roster call-rate gate for ATAC-demux\">");
        bcf_hdr_append(candidate_header,
            "##INFO=<ID=ATAC_ROSTER_CALL_RATE,Number=1,Type=Float,Description=\"Complete diploid-genotype call rate among donors used by supplied pool combinations\">");
        bcf_hdr_append(candidate_header,
            "##INFO=<ID=ATAC_REGULATORY,Number=1,Type=Integer,Description=\"1 when the site overlaps the supplied regulatory BED or default strand-aware 2-kb upstream promoter prior; ranking annotation, not an eligibility filter\">");
        bcf_hdr_append(candidate_header,
            "##INFO=<ID=ATAC_RNA_SOFT_OVERLAP,Number=1,Type=Integer,Description=\"1 when this locus occurs in the supplied selected RNA-demux panel\">");
        bcf_hdr_append(candidate_header,
            "##cellbouncer_atac_candidate_pool_schema=atac_candidate_union_v2");
        bcf_hdr_append(candidate_header,
            "##cellbouncer_atac_candidate_scope=demux_het_species_union");
        bcf_hdr_append(candidate_header,
            "##cellbouncer_atac_candidate_filters=biallelic_snv,distinct_called_dosages_or_heterozygous,usable_gt,atac_coverage,nuclear_safe,numt_safe");
        bcf_hdr_append(candidate_header,
            "##cellbouncer_atac_demux_annotation_gate=operational_roster_call_rate");
        bcf_hdr_append(candidate_header,
            "##cellbouncer_atac_demux_selection_gate=atac_coverage_qc_nuclear_safety_with_scalar_pool_ranking");
        bcf_hdr_append(candidate_header,
            "##cellbouncer_atac_pool_objective_policy=rna_style_scalar_pool_discrimination_no_identity_pair_floors");
        bcf_hdr_append(candidate_header,
            "##cellbouncer_atac_rna_overlap_policy=soft_penalty");
        {
            char overlap_multiplier_header[128];
            snprintf(overlap_multiplier_header,
                sizeof(overlap_multiplier_header),
                "##cellbouncer_atac_soft_overlap_multiplier=%.6f",
                atac_soft_overlap_multiplier);
            bcf_hdr_append(candidate_header, overlap_multiplier_header);
        }
        {
            string regulatory_header =
                "##cellbouncer_atac_regulatory_source=" + atac_regulatory_source;
            bcf_hdr_append(candidate_header, regulatory_header.c_str());
        }
        if (!atac_excluded_positions.empty()) {
            bcf_hdr_append(candidate_header,
                "##cellbouncer_atac_rna_locus_exclusion=position_only");
        }
        if (bcf_hdr_sync(candidate_header) < 0) {
            fprintf(stderr, "ERROR: Failed to sync ATAC candidate header\n");
            exit(1);
        }

        htsFile* candidate_output = hts_open(atac_candidates_outfile.c_str(), "wb");
        if (!candidate_output) {
            fprintf(stderr, "ERROR: Failed to open ATAC candidate output: %s\n",
                atac_candidates_outfile.c_str());
            exit(1);
        }
        hts_set_threads(candidate_output, n_threads);
        if (bcf_hdr_write(candidate_output, candidate_header) < 0) {
            fprintf(stderr, "ERROR: Failed to write ATAC candidate header\n");
            exit(1);
        }

        bcf1_t* candidate_record = bcf_init();
        vector<int32_t> candidate_gt;
        for (size_t c = 0; c < chroms.size(); ++c) {
            size_t range_start = chrom_ranges[c].first;
            size_t range_end = chrom_ranges[c].second;
            for (size_t vi = range_start; vi < range_end; ++vi) {
                const StoredVariant& sv = all_variants[vi];
                if (collect_candidate_only_audits) {
                    const bool overlap = sv.atac_soft_overlap;
                    candidate_audit_gate_counts[0]++;
                    if (overlap) candidate_audit_gate_overlap[0]++;
                    if (sv.is_biallelic && sv.pattern != NULL) {
                        candidate_audit_gate_counts[1]++;
                        if (overlap) candidate_audit_gate_overlap[1]++;
                        if (sv.nuclear_exclusion == NUCLEAR_EXCLUDE_NONE) {
                            candidate_audit_gate_counts[2]++;
                            if (overlap) candidate_audit_gate_overlap[2]++;
                            if (sv.atac_cov_score >= atac_min_cov) {
                                candidate_audit_gate_counts[3]++;
                                if (overlap) candidate_audit_gate_overlap[3]++;
                                if (sv.passes_qc) {
                                    candidate_audit_gate_counts[4]++;
                                    if (overlap)
                                        candidate_audit_gate_overlap[4]++;
                                    if (!panel_position_is_excluded(
                                            atac_excluded_positions,
                                            sv.chrom_idx, sv.pos)) {
                                        candidate_audit_gate_counts[5]++;
                                        if (overlap)
                                            candidate_audit_gate_overlap[5]++;
                                        if (operational_roster_call_rate(
                                                sv.pattern,
                                                atac_roster_sample_indices) >=
                                            pairwise_rescue_min_call_rate) {
                                            candidate_audit_gate_counts[6]++;
                                            if (overlap)
                                                candidate_audit_gate_overlap[6]++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if (sv.nuclear_exclusion != NUCLEAR_EXCLUDE_NONE) continue;
                if (!sv.is_biallelic || sv.genotypes_packed.empty()) continue;
                if (panel_position_is_excluded(
                        atac_excluded_positions, sv.chrom_idx, sv.pos)) continue;

                size_t comma = sv.alleles.find(',');
                if (comma != 1 || sv.alleles.size() != 3) continue;
                const char ref = sv.alleles[0];
                const char alt = sv.alleles[2];
                const bool valid_ref = ref == 'A' || ref == 'C' || ref == 'G' || ref == 'T';
                const bool valid_alt = alt == 'A' || alt == 'C' || alt == 'G' || alt == 'T';
                if (!valid_ref || !valid_alt || ref == alt) continue;

                float atac_cov = sv.atac_cov_score;
                if (atac_cov < atac_min_cov) continue;

                HetStats hs = compute_het_stats_packed(
                    sv.genotypes_packed, num_samples);
                if (!sv.passes_qc && hs.n_het <= 0) continue;

                bcf_clear(candidate_record);
                candidate_record->rid = sv.chrom_idx;
                candidate_record->pos = sv.pos;
                candidate_record->qual = sv.qual;
                bcf_update_alleles_str(
                    candidate_header, candidate_record, sv.alleles.c_str());
                unpack_genotypes(
                    sv.genotypes_packed, num_samples, candidate_gt);
                bcf_update_genotypes(
                    candidate_header, candidate_record,
                    candidate_gt.data(), candidate_gt.size());

                int32_t n_called = hs.n_called;
                int32_t n_het = hs.n_het;
                float het_freq = hs.het_freq;
                float missing_freq = hs.missing_freq;
                int roster_called = 0;
                for (int sample_index : atac_roster_sample_indices) {
                    if (sv.pattern && sv.pattern->present.test(sample_index)) {
                        roster_called++;
                    }
                }
                float roster_call_rate = atac_roster_sample_indices.empty()
                    ? 0.0f
                    : (float)roster_called / atac_roster_sample_indices.size();
                int32_t demux_eligible = sv.passes_qc &&
                    roster_call_rate >= pairwise_rescue_min_call_rate ? 1 : 0;
                int32_t regulatory = sv.atac_regulatory ? 1 : 0;
                int32_t soft_overlap = sv.atac_soft_overlap ? 1 : 0;
                int bin_idx = get_bin_idx(sv.pos, bin_size);
                string bin_id = make_bin_id(chroms[c], bin_idx);

                if (collect_candidate_only_audits &&
                    (!atac_regulatory_required || sv.atac_regulatory)) {
                    auto union_found =
                        candidate_audit_union_pattern_index.find(sv.pattern);
                    if (union_found ==
                        candidate_audit_union_pattern_index.end()) {
                        CandidateAuditPatternCount summary;
                        summary.pattern = sv.pattern;
                        summary.representative_variant_idx = vi;
                        summary.count = 1;
                        summary.soft_overlap_count = soft_overlap;
                        candidate_audit_union_pattern_index[sv.pattern] =
                            candidate_audit_union_patterns.size();
                        candidate_audit_union_patterns.push_back(summary);
                    } else {
                        CandidateAuditPatternCount& summary =
                            candidate_audit_union_patterns[union_found->second];
                        summary.count++;
                        summary.soft_overlap_count += soft_overlap;
                    }
                    if (demux_eligible) {
                        auto demux_found =
                            candidate_audit_demux_pattern_index.find(sv.pattern);
                        if (demux_found ==
                            candidate_audit_demux_pattern_index.end()) {
                            CandidateAuditPatternCount summary;
                            summary.pattern = sv.pattern;
                            summary.representative_variant_idx = vi;
                            summary.count = 1;
                            summary.soft_overlap_count = soft_overlap;
                            candidate_audit_demux_pattern_index[sv.pattern] =
                                candidate_audit_demux_patterns.size();
                            candidate_audit_demux_patterns.push_back(summary);
                        } else {
                            CandidateAuditPatternCount& summary =
                                candidate_audit_demux_patterns[
                                    demux_found->second];
                            summary.count++;
                            summary.soft_overlap_count += soft_overlap;
                        }
                    }
                }

                bcf_update_info_float(
                    candidate_header, candidate_record, "ATAC_COV_SCORE", &atac_cov, 1);
                bcf_update_info_string(
                    candidate_header, candidate_record, "ATAC_BIN_ID", bin_id.c_str());
                bcf_update_info_int32(
                    candidate_header, candidate_record, "ATAC_N_CALLED", &n_called, 1);
                bcf_update_info_int32(
                    candidate_header, candidate_record, "ATAC_N_HET", &n_het, 1);
                bcf_update_info_float(
                    candidate_header, candidate_record, "ATAC_HET_FREQ", &het_freq, 1);
                bcf_update_info_float(
                    candidate_header, candidate_record, "ATAC_MISSING_FREQ", &missing_freq, 1);
                bcf_update_info_int32(
                    candidate_header, candidate_record,
                    "ATAC_DEMUX_ELIGIBLE", &demux_eligible, 1);
                bcf_update_info_float(
                    candidate_header, candidate_record,
                    "ATAC_ROSTER_CALL_RATE", &roster_call_rate, 1);
                bcf_update_info_int32(
                    candidate_header, candidate_record,
                    "ATAC_REGULATORY", &regulatory, 1);
                bcf_update_info_int32(
                    candidate_header, candidate_record,
                    "ATAC_RNA_SOFT_OVERLAP", &soft_overlap, 1);

                if (bcf_write(candidate_output, candidate_header, candidate_record) < 0) {
                    fprintf(stderr, "ERROR: Failed to write ATAC candidate record\n");
                    exit(1);
                }
                atac_candidates_written++;
                atac_demux_eligible_written += demux_eligible;
            }
        }

        bcf_destroy(candidate_record);
        hts_close(candidate_output);
        bcf_hdr_destroy(candidate_header);
        if (atac_candidates_written <= 0) {
            fprintf(stderr, "ERROR: Reusable ATAC candidate pool is empty\n");
            exit(1);
        }
        fprintf(stderr, "Indexing reusable ATAC candidate BCF...\n");
        if (bcf_index_build(atac_candidates_outfile.c_str(), 14) < 0) {
            fprintf(stderr, "ERROR: Failed to index ATAC candidate BCF: %s\n",
                atac_candidates_outfile.c_str());
            exit(1);
        }
        fprintf(stderr,
            "ATAC candidate pool: %ld total union candidates; %ld pass the "
            "operational-roster call-rate demux pre-gate\n",
            atac_candidates_written, atac_demux_eligible_written);

        if (atac_candidates_only) {
            if (!atac_audit_prefix.empty()) {
                // Candidate-only audits describe the complete ATAC-eligible
                // pool. Pool discrimination is a scalar ranking annotation,
                // never a gate or a collection of optimization objectives.
                vector<uint8_t> pattern_informative(
                    candidate_audit_demux_patterns.size(), 1);
                const string candidate_audit_path =
                    atac_audit_prefix + ".candidate_sfs.tsv";
                ofstream candidate_audit(candidate_audit_path.c_str());
                if (!candidate_audit) {
                    fprintf(stderr,
                        "ERROR: Could not write candidate-only ATAC SFS audit: %s\n",
                        candidate_audit_path.c_str());
                    exit(1);
                }
                candidate_audit <<
                    "section\tmetric\tn_sites\tsoft_overlap_sites\n";
                const char* gate_names[candidate_audit_n_gates] = {
                    "source_loaded", "biallelic_with_genotypes",
                    "nuclear_safe", "atac_covered",
                    "distinct_called_dosages", "legacy_hard_exclusion_safe",
                    "roster_call_rate"
                };
                for (int gate = 0; gate < candidate_audit_n_gates; ++gate) {
                    candidate_audit << "gate\t" << gate_names[gate] << '\t'
                        << candidate_audit_gate_counts[gate] << '\t'
                        << candidate_audit_gate_overlap[gate] << '\n';
                }
                long informative_sites = 0;
                long informative_overlap = 0;
                const char* mac_names[6] = {
                    "MAC_1", "MAC_2", "MAC_3_5", "MAC_6_10",
                    "MAC_11_20", "MAC_21_PLUS"
                };
                long mac_sites[6] = {0,0,0,0,0,0};
                long mac_overlap[6] = {0,0,0,0,0,0};
                const char* call_names[4] = {
                    "CALL_RATE_0.50_0.75", "CALL_RATE_0.75_0.90",
                    "CALL_RATE_0.90_LT1", "CALL_RATE_1.00"
                };
                long call_sites[4] = {0,0,0,0};
                long call_overlap[4] = {0,0,0,0};
                for (size_t pattern_index = 0;
                     pattern_index < candidate_audit_demux_patterns.size();
                     ++pattern_index) {
                    if (!pattern_informative[pattern_index]) continue;
                    const CandidateAuditPatternCount& summary =
                        candidate_audit_demux_patterns[pattern_index];
                    informative_sites += summary.count;
                    informative_overlap += summary.soft_overlap_count;
                    int called_count = 0;
                    int alt_count = 0;
                    for (int sample_index : atac_roster_sample_indices) {
                        const int dosage = pattern_dosage(
                            summary.pattern, sample_index);
                        if (dosage < 0) continue;
                        called_count++;
                        alt_count += dosage;
                    }
                    const int mac = min(
                        alt_count, 2 * called_count - alt_count);
                    const int mac_bin = mac <= 1 ? 0 : mac == 2 ? 1 :
                        mac <= 5 ? 2 : mac <= 10 ? 3 : mac <= 20 ? 4 : 5;
                    const float call_rate = atac_roster_sample_indices.empty()
                        ? 0.0f : (float)called_count /
                          atac_roster_sample_indices.size();
                    const int call_bin = call_rate < 0.75f ? 0 :
                        call_rate < 0.90f ? 1 :
                        call_rate < 1.0f ? 2 : 3;
                    mac_sites[mac_bin] += summary.count;
                    mac_overlap[mac_bin] += summary.soft_overlap_count;
                    call_sites[call_bin] += summary.count;
                    call_overlap[call_bin] += summary.soft_overlap_count;
                }
                candidate_audit << "gate\tatac_eligible_scalar_pool_ranked\t"
                    << informative_sites << '\t' << informative_overlap << '\n';
                for (int bin = 0; bin < 6; ++bin) {
                    candidate_audit << "mac\tcandidate_" << mac_names[bin]
                        << '\t' << mac_sites[bin] << '\t'
                        << mac_overlap[bin] << '\n';
                    candidate_audit << "mac\tselected_" << mac_names[bin]
                        << "\t0\t0\n";
                }
                for (int bin = 0; bin < 4; ++bin) {
                    candidate_audit << "call_rate\tcandidate_"
                        << call_names[bin] << '\t' << call_sites[bin] << '\t'
                        << call_overlap[bin] << '\n';
                    candidate_audit << "call_rate\tselected_"
                        << call_names[bin] << "\t0\t0\n";
                }
                candidate_audit <<
                    "selection\tatac_demux_selected\t0\t0\n";
                candidate_audit.close();

                const string pair_audit_path =
                    atac_audit_prefix + ".pair_capacity.tsv";
                ofstream pair_audit(pair_audit_path.c_str());
                if (!pair_audit) {
                    fprintf(stderr,
                        "ERROR: Could not write candidate-only pair audit: %s\n",
                        pair_audit_path.c_str());
                    exit(1);
                }
                pair_audit <<
                    "scope\tdataset\tobjective_id\tsignature_a\tsignature_b\t"
                    "library_ids\tcandidate_joint_callable\t"
                    "candidate_distinguishing\tcandidate_rna_demux_overlap\t"
                    "candidate_rna_disjoint\tcandidate_soft_distinguishing\t"
                    "candidate_soft_rna_demux_overlap\t"
                    "selected_distinguishing\tselected_soft_distinguishing\t"
                    "requested_floor\teffective_floor\tstatus\n";
                for (const PoolDonorPair& donor_pair : atac_pool_donor_pairs) {
                    long joint = 0, strong = 0, overlap = 0;
                    long soft = 0, soft_overlap = 0;
                    for (size_t pattern_index = 0;
                         pattern_index < candidate_audit_demux_patterns.size();
                         ++pattern_index) {
                        if (!pattern_informative[pattern_index]) continue;
                        const CandidateAuditPatternCount& summary =
                            candidate_audit_demux_patterns[pattern_index];
                        const int first = pattern_dosage(
                            summary.pattern, donor_pair.first_index);
                        const int second = pattern_dosage(
                            summary.pattern, donor_pair.second_index);
                        if (first < 0 || second < 0) continue;
                        joint += summary.count;
                        const float difference = fabsf(
                            (float)first / 2.0f - (float)second / 2.0f);
                        if (difference >= min_pair_score) {
                            strong += summary.count;
                            overlap += summary.soft_overlap_count;
                        } else if (difference >= 0.25f) {
                            soft += summary.count;
                            soft_overlap += summary.soft_overlap_count;
                        }
                    }
                    string library_list;
                    for (size_t li = 0;
                         li < donor_pair.library_ids.size(); ++li) {
                        if (li) library_list += ",";
                        library_list += to_string(donor_pair.library_ids[li]);
                    }
                    string first_signature = donor_pair.first_name + ":1";
                    string second_signature = donor_pair.second_name + ":1";
                    if (second_signature < first_signature)
                        swap(first_signature, second_signature);
                    pair_audit << "donor_pair\tatac_demux\t"
                        << first_signature << "||" << second_signature << '\t'
                        << first_signature << '\t' << second_signature << '\t'
                        << library_list << '\t' << joint << '\t' << strong
                        << '\t' << overlap << '\t' << strong - overlap << '\t'
                        << soft << '\t' << soft_overlap
                        << "\t0\t0\t0\t0\t"
                        << "DIAGNOSTIC_ONLY_NOT_A_SELECTION_OBJECTIVE\n";
                }
                pair_audit.close();

                const string species_audit_path =
                    atac_audit_prefix + ".species_tiers.tsv";
                ofstream species_audit(species_audit_path.c_str());
                if (!species_audit) {
                    fprintf(stderr,
                        "ERROR: Could not write candidate-only species audit: %s\n",
                        species_audit_path.c_str());
                    exit(1);
                }
                species_audit <<
                    "species_a\tspecies_b\toutgroup_pair\ttier\t"
                    "min_pair_discrim\tmax_pair_discrim\tcandidate_sites\t"
                    "soft_overlap_sites\tselected_sites\n";
                const int n_species_pairs = panel_metadata.n_pairs;
                vector<array<long, 2>> tier_candidate(
                    n_species_pairs, {0, 0});
                vector<array<long, 2>> tier_overlap(
                    n_species_pairs, {0, 0});
                vector<uint8_t> pair_has_outgroup(n_species_pairs, 0);
                for (int pi = 0; pi < n_species_pairs; ++pi) {
                    const pair<string, string>& species_pair =
                        panel_metadata.species_pairs[pi];
                    if (!species_outgroup.empty() &&
                        (species_pair.first == species_outgroup ||
                         species_pair.second == species_outgroup)) {
                        pair_has_outgroup[pi] = 1;
                    }
                }
                const float species_soft_threshold =
                    species_min_pair_discrim * 0.5f;
                long outgroup_only_candidate = 0;
                long outgroup_only_overlap = 0;
                long union_candidate_sites = 0;
                for (const CandidateAuditPatternCount& summary :
                     candidate_audit_union_patterns) {
                    union_candidate_sites += summary.count;
                    const StoredVariant& representative = all_variants[
                        summary.representative_variant_idx];
                    vector<float> pair_score, pair_discrim, mean_alt,
                        major_dosage_fraction, call_rate;
                    vector<int> n_called;
                    compute_species_discrim(
                        representative.genotypes_packed, num_samples,
                        panel_metadata, species_min_called_per_species,
                        species_major_dosage_frac_floor,
                        species_call_rate_floor, pair_score, pair_discrim,
                        mean_alt, major_dosage_fraction, call_rate, n_called);
                    bool any_soft = false;
                    bool any_non_outgroup_soft = false;
                    for (int pi = 0; pi < n_species_pairs; ++pi) {
                        int tier = -1;
                        if (pair_discrim[pi] >= species_min_pair_discrim)
                            tier = 0;
                        else if (pair_discrim[pi] >= species_soft_threshold)
                            tier = 1;
                        if (tier < 0) continue;
                        tier_candidate[pi][tier] += summary.count;
                        tier_overlap[pi][tier] +=
                            summary.soft_overlap_count;
                        any_soft = true;
                        if (!pair_has_outgroup[pi])
                            any_non_outgroup_soft = true;
                    }
                    if (any_soft && !any_non_outgroup_soft) {
                        outgroup_only_candidate += summary.count;
                        outgroup_only_overlap += summary.soft_overlap_count;
                    }
                }
                for (int pi = 0; pi < n_species_pairs; ++pi) {
                    for (int tier = 0; tier < 2; ++tier) {
                        species_audit
                            << panel_metadata.species_pairs[pi].first << '\t'
                            << panel_metadata.species_pairs[pi].second << '\t'
                            << (pair_has_outgroup[pi] ? 1 : 0) << '\t'
                            << (tier == 0 ? "strict" : "soft") << '\t'
                            << (tier == 0 ? species_min_pair_discrim
                                          : species_soft_threshold) << '\t'
                            << (tier == 0 ? 1.0f
                                          : species_min_pair_discrim) << '\t'
                            << tier_candidate[pi][tier] << '\t'
                            << tier_overlap[pi][tier] << "\t0\n";
                    }
                }
                species_audit
                    << "__SUMMARY__\t__SUMMARY__\t0\tACTUAL_FINAL_TOTAL\t"
                    << "0\t1\t" << union_candidate_sites << "\t0\t0\n";
                species_audit
                    << (species_outgroup.empty() ? "O" : species_outgroup)
                    << "\t__ALL__\t1\tOUTGROUP_ONLY_ACTUAL\t"
                    << species_soft_threshold << "\t1\t"
                    << outgroup_only_candidate << '\t'
                    << outgroup_only_overlap << "\t0\n";
                species_audit.close();

                fprintf(stderr,
                    "Candidate-only audits: %s, %s, %s\n",
                    candidate_audit_path.c_str(), pair_audit_path.c_str(),
                    species_audit_path.c_str());
            }
            hts_idx_destroy(vcf_idx);
            fprintf(stderr, "ATAC candidate-only mode complete: %s\n",
                atac_candidates_outfile.c_str());
            return 0;
        }
    }

    // ========================================================================
    // PASS 1: Count clades AND collect het candidates (PARALLEL BY CHROMOSOME)
    // Now using in-memory variants instead of disk I/O
    // ========================================================================

    unordered_map<bitset<NBITS>, double> branchcounts;
    unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, double> branchcounts_missing;
    unordered_map<bitset<NBITS>, bitset<NBITS>> miss2flip;
    vector<HetCandidate> all_het_candidates;  // Global het candidates
    vector<SpeciesCandidate> all_species_candidates; // Global species candidates

    mutex bc_mutex, bcm_mutex, m2f_mutex, het_mutex, species_mutex;
    atomic<long> total_snps(0);
    atomic<long> total_het_candidates(0);
    atomic<long> total_species_candidates(0);
    atomic<int> chroms_done(0);

    auto start_time = chrono::high_resolution_clock::now();

    bool collect_het = !het_outfile.empty() || !atac_het_outfile.empty();
    bool collect_species = !species_outfile.empty();

    fprintf(stderr, "Pass 1: Counting mutations on branches%s%s (parallel, %d threads)...\n",
        collect_het ? " and collecting het candidates" : "",
        collect_species ? " and collecting species candidates" : "",
        n_threads);

    #pragma omp parallel num_threads(n_threads)
    {
        // Thread-local data structures (SAME as before)
        unordered_map<bitset<NBITS>, double> local_branchcounts;
        unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, double> local_branchcounts_missing;
        unordered_map<bitset<NBITS>, bitset<NBITS>> local_miss2flip;
        vector<HetCandidate> local_het_candidates;
        vector<SpeciesCandidate> local_species_candidates;
        
        #pragma omp for schedule(dynamic, 1)
        for (size_t c = 0; c < chroms.size(); c++){
            const string& chrom = chroms[c];
            
            // Get range of variants for this chromosome from memory
            size_t range_start = chrom_ranges[c].first;
            size_t range_end = chrom_ranges[c].second;
            
            if (range_start >= range_end) {
                chroms_done++;
                continue;
            }
            
            long local_snps = 0;
            long local_het = 0;
            
            // Iterate over stored variants for this chromosome
            for (size_t vi = range_start; vi < range_end; vi++) {
                const StoredVariant& sv = all_variants[vi];
                if (sv.nuclear_exclusion != NUCLEAR_EXCLUDE_NONE) continue;
                
                // Skip non-biallelic (same as original: if n_allele != 2 continue)
                if (!sv.is_biallelic) continue;
                
                // Use precomputed proc_bcf_record result
                if (sv.passes_qc) {
                    if (wants_generic_allocation) {
                    // Make copies since count_branch/count_branch_missing take non-const refs
                    bitset<NBITS> alt = sv.pattern->alt;
                    bitset<NBITS> alt_flip = sv.pattern->alt_flip;
                    bitset<NBITS> present = sv.pattern->present;
                    
                    int ta = sv.total_alt;
                    int ta_flip = sv.total_alt_flip;
                
                    if (present.count() == (size_t)num_samples){
                        count_branch(local_branchcounts, alt, alt_flip, ta, ta_flip, 1.0);
                    }
                    else{
                        count_branch_missing(local_branchcounts_missing, local_miss2flip, 
                            present, alt, alt_flip, ta, ta_flip, 1.0);
                    }
                    }
                    local_snps++;
                }
                
                // Collect het candidates if enabled
                // Using packed genotype accessors
                if (collect_het && !sv.genotypes_packed.empty() && sv.is_biallelic) {
                    // Parse ref and alt from alleles string
                    size_t comma_pos = sv.alleles.find(',');
                    string ref_allele = sv.alleles.substr(0, comma_pos);
                    string alt_allele = sv.alleles.substr(comma_pos + 1);
                    
                    // Check alleles are valid
                    bool het_pass = true;
                    if (ref_allele.length() != 1 || alt_allele.length() != 1) {
                        het_pass = false;
                    } else {
                        char r = ref_allele[0];
                        char a = alt_allele[0];
                        if ((r != 'A' && r != 'C' && r != 'G' && r != 'T') ||
                            (a != 'A' && a != 'C' && a != 'G' && a != 'T') ||
                            r == a) {
                            het_pass = false;
                        }
                    }
                    
                    if (het_pass) {
                        if (sv.nuclear_exclusion != NUCLEAR_EXCLUDE_NONE) {
                            het_pass = false;
                        }
                    }
                    
                    if (het_pass) {
                        HetStats hs = compute_het_stats_packed(sv.genotypes_packed, num_samples);
                        int n_het = hs.n_het;
                        int n_called = hs.n_called;
                        int n_missing = hs.n_missing;
                        
                        // Only proceed if at least one het
                        if (n_het > 0) {
                            // Preserve the historical RNA-het excess-
                            // heterozygosity gate, but do not apply a
                            // pan-species HWE assumption to ATAC-het.  The
                            // flag is carried on the shared candidate and is
                            // enforced only by the RNA selector below.
                            bool passes_rna_het_hwe = true;
                            if (n_called >= 5) {
                                float n_alt_alleles_b7 = 0;
                                for (int ii = 0; ii < num_samples; ++ii) {
                                    int aa1 = gt_packed_allele(sv.genotypes_packed, ii, 0);
                                    int aa2 = gt_packed_allele(sv.genotypes_packed, ii, 1);
                                    if (aa1 < 0 || aa2 < 0) continue;
                                    if (aa1 > 0) n_alt_alleles_b7 += 1;
                                    if (aa2 > 0) n_alt_alleles_b7 += 1;
                                }
                                float p_af = n_alt_alleles_b7 / (2.0f * n_called);
                                float expected_het_hw = 2.0f * p_af * (1.0f - p_af);
                                float observed_het_hw = (float)n_het / n_called;
                                float var_het_hw = expected_het_hw * (1.0f - expected_het_hw) / n_called;
                                if (var_het_hw > 0) {
                                    float z_hw = (observed_het_hw - expected_het_hw) / sqrtf(var_het_hw);
                                    if (z_hw > het_max_excess_z) {
                                        passes_rna_het_hwe = false;
                                    }
                                }
                            }
                            if (!passes_rna_het_hwe && atac_het_outfile.empty()) {
                                continue;
                            }
                            
                            float het_freq = hs.het_freq;
                            float missing_freq = hs.missing_freq;
                            
                            const float annot_score = sv.annot_score;
                            const float cov_score = sv.rna_cov_score;
                            
                            // D9: Do NOT hard-filter here. A candidate is eligible if RNA cov OR ATAC cov
                            // passes its respective threshold; downstream selectors apply their own gate.
                            bool rna_pass = (cov_score >= min_cov);
                            bool atac_pass = !atac_coverage_map.empty() &&
                                             sv.atac_cov_score >= atac_min_cov;
                            if (!rna_pass && !atac_pass) {
                                continue;  // Skip only if neither modality clears the threshold
                            }
                            
                            // Get ATAC coverage if available
                            float atac_cov = sv.atac_cov_score;
                            
                            // Compute pool discrimination if enabled
                            float pool_dn = 0.0f;
                            if (enable_pool_scoring && has_pool_combos) {
                                pool_dn = compute_pool_discrim(sv.genotypes_packed, num_samples,
                                    pool_combos, sample_name_to_idx,
                                    min_pair_score);
                            }
                            
                            // Compute final_het_score with weighted pool combination
                            // (Item F): final = (1 - w) * het_freq + w * pool_discrim_norm
                            // bounded in [0, 1+] but each summand is approximately [0,1].
                            float final_het_score = het_freq;
                            if (enable_pool_scoring) {
                                final_het_score = (1.0f - het_pool_weight) * het_freq
                                                + het_pool_weight * pool_dn;
                            }
                            
                            // Het scoring with call-rate penalty (Item F):
                            // het_score = final_het_score * call_rate * (log2(cov+1) + 0.1 * annot_boost)
                            // Penalizes sites called in few individuals so they cannot dominate
                            // simply by having high het_freq among a tiny called set.
                            float call_rate = 1.0f - missing_freq;
                            float annot_boost = (annot_score >= 1.0f) ? 1.0f : 0.0f;
                            float selection_component = log2f(cov_score + 1.0f) + (0.1f * annot_boost);
                            float het_score = final_het_score * call_rate * selection_component;
                            
                            HetCandidate hc;
                            hc.chrom = chrom;
                            hc.pos = sv.pos;
                            hc.chrom_idx = c;
                            hc.bin_idx = get_bin_idx(sv.pos, bin_size);
                            hc.het_score = het_score;
                            hc.het_freq = het_freq;
                            hc.missing_freq = missing_freq;
                            hc.annot_score = annot_score;
                            hc.cov_score = cov_score;  // Store RAW coverage, same as demux
                            hc.atac_cov_score = atac_cov;
                            hc.atac_regulatory = sv.atac_regulatory;
                            hc.passes_rna_het_hwe = passes_rna_het_hwe;
                            hc.pool_discrim_norm = pool_dn;
                            hc.het_sel_score = het_score; // selection ranking score
                            hc.atac_het_sel_score = 0.0f;
                            hc.ref_allele = ref_allele;
                            hc.alt_allele = alt_allele;
                            hc.variant_idx = vi;  // deferred genotype fetch from all_variants
                            hc.n_samples = num_samples;
                            
                            local_het_candidates.push_back(std::move(hc));
                            local_het++;
                        }
                    }
                }
            
            // Collect species candidates if enabled
            if (collect_species && !sv.genotypes_packed.empty() && sv.is_biallelic) {
                // Parse ref/alt
                size_t comma_pos2 = sv.alleles.find(',');
                string ref_a2 = sv.alleles.substr(0, comma_pos2);
                string alt_a2 = sv.alleles.substr(comma_pos2 + 1);
                
                // Check alleles are valid single-base
                bool sp_pass = true;
                if (ref_a2.length() != 1 || alt_a2.length() != 1) {
                    sp_pass = false;
                } else {
                    char r = ref_a2[0], a = alt_a2[0];
                    if ((r != 'A' && r != 'C' && r != 'G' && r != 'T') ||
                        (a != 'A' && a != 'C' && a != 'G' && a != 'T') || r == a) {
                        sp_pass = false;
                    }
                }
                
                if (sp_pass) {
                    // B11: Blacklist filter
                    if (!blacklist.empty() && is_blacklisted(chrom, sv.pos, blacklist)) {
                        sp_pass = false;
                    }
                }

                if (sp_pass && species_coverage_mode == "atac" &&
                    panel_position_is_excluded(
                        atac_excluded_positions, sv.chrom_idx, sv.pos)) {
                    sp_pass = false;
                }
                if (sp_pass && species_coverage_mode == "atac" &&
                    atac_regulatory_required &&
                    !sv.atac_regulatory) {
                    sp_pass = false;
                }
                
                if (sp_pass) {
                    // Get coverage for species eligibility
                    float sp_cov = 0.0f;
                    if (species_coverage_mode == "rna" && !coverage_map.empty()) {
                        sp_cov = sv.rna_cov_score;
                    } else if (species_coverage_mode == "atac" && !atac_coverage_map.empty()) {
                        sp_cov = sv.atac_cov_score;
                    }
                    
                    float species_cov_min = (species_coverage_mode == "atac")
                                           ? atac_min_cov : min_cov;
                    if (sp_cov >= species_cov_min) {
                        // Compute pair discrimination (primary) and recorded legacy fields
                        vector<float> ps, pd, sp_ma, sp_mdf, sp_cr;
                        vector<int> sp_nc;
                        compute_species_discrim(
                            sv.genotypes_packed, num_samples, panel_metadata,
                            species_min_called_per_species,
                            species_major_dosage_frac_floor,
                            species_call_rate_floor,
                            ps, pd, sp_ma, sp_mdf, sp_cr, sp_nc
                        );

                        // Find max pair_discrim (primary ranking) and max pair_score (legacy)
                        float max_pd = 0.0f;
                        float max_ps = 0.0f;
                        for (float v : pd) { if (v > max_pd) max_pd = v; }
                        for (float v : ps) { if (v > max_ps) max_ps = v; }

                        // Pre-filter using pair_discrim (primary). Permissive gate at half threshold.
                        // Note: this pre-filter is intentionally permissive so backfill candidates
                        // remain available for the global utility pass (see species selection block).
                        if (max_pd >= species_min_pair_discrim * 0.5f) {
                            float cov_term = min(log2f(sp_cov + 1.0f), species_cov_cap);
                            
                            // Compute per-pair dosage separation while sp_ma is still local.
                            // dosage_sep[pi] = |mean_alt[A] - mean_alt[B]| / 2 in [0, 1].
                            // Used as multi-pair contribution weight in global utility (Item A).
                            int n_pairs_local = (int)pd.size();
                            vector<float> pds(n_pairs_local, 0.0f);
                            for (int pi = 0; pi < n_pairs_local; pi++) {
                                const string& spa = panel_metadata.species_pairs[pi].first;
                                const string& spb = panel_metadata.species_pairs[pi].second;
                                int sai = -1, sbi = -1;
                                for (int si = 0; si < (int)panel_metadata.species_list.size(); si++) {
                                    if (panel_metadata.species_list[si] == spa) sai = si;
                                    if (panel_metadata.species_list[si] == spb) sbi = si;
                                }
                                if (sai >= 0 && sbi >= 0) {
                                    pds[pi] = fabsf(sp_ma[sai] - sp_ma[sbi]) / 2.0f;
                                }
                            }
                            
                            SpeciesCandidate sc;
                            sc.chrom = chrom;
                            sc.pos = sv.pos;
                            sc.chrom_idx = c;
                            sc.bin_idx = get_bin_idx(sv.pos, bin_size);
                            sc.cov_score = sp_cov;
                            sc.pair_score = std::move(ps);
                            sc.pair_discrim = std::move(pd);
                            sc.pair_dosage_sep = std::move(pds);
                            sc.max_pair_score = max_ps;
                            sc.max_pair_discrim = max_pd;
                            sc.selection_score = max_ps * cov_term;  // legacy field, retained
                            sc.atac_regulatory =
                                species_coverage_mode == "atac" && sv.atac_regulatory;
                            sc.assigned_pair = -1;
                            sc.ref_allele = ref_a2;
                            sc.alt_allele = alt_a2;
                            sc.variant_idx = vi;  // deferred genotype fetch from all_variants
                            sc.n_samples = num_samples;
                            sc.sp_mean_alt = std::move(sp_ma);
                            sc.sp_major_dosage_frac = std::move(sp_mdf);
                            sc.sp_call_rate = std::move(sp_cr);
                            sc.sp_n_called = std::move(sp_nc);
                            local_species_candidates.push_back(std::move(sc));
                        }
                    }
                }
            }
            } // end per-variant for loop
            
            total_snps += local_snps;
            total_het_candidates += local_het;
            int done = ++chroms_done;
            
            if (local_snps > 0 || done % 100 == 0) {
                #pragma omp critical(stderr_log)
                {
                fprintf(stderr, "  Chromosome %s: %ld SNPs, %ld het, %ld sp [%d/%zu]\n", 
                    chrom.c_str(), local_snps, local_het, (long)local_species_candidates.size(), done, chroms.size());
                }
            }
        }
        
        // Merge thread-local counts into global (SAME as original)
        if (wants_generic_allocation) {
            merge_branchcounts(branchcounts, local_branchcounts, bc_mutex);
            merge_branchcounts_missing(branchcounts_missing, local_branchcounts_missing, bcm_mutex);
            merge_miss2flip(miss2flip, local_miss2flip, m2f_mutex);
        }
        if (collect_het) {
            merge_het_candidates(all_het_candidates, local_het_candidates, het_mutex);
        }
        if (collect_species) {
            merge_species_candidates(all_species_candidates, local_species_candidates, species_mutex);
        }
    }
    
    // Stable sort by genomic position for reproducibility across runs
    if (collect_het) {
        sort(all_het_candidates.begin(), all_het_candidates.end(),
            [](const HetCandidate& a, const HetCandidate& b) {
                if (a.chrom_idx != b.chrom_idx) return a.chrom_idx < b.chrom_idx;
                return a.pos < b.pos;
            });
    }
    if (collect_species) {
        sort(all_species_candidates.begin(), all_species_candidates.end(),
            [](const SpeciesCandidate& a, const SpeciesCandidate& b) {
                if (a.chrom_idx != b.chrom_idx) return a.chrom_idx < b.chrom_idx;
                return a.pos < b.pos;
            });
    }
    
    // vcf_idx no longer needed after Pass 1 (C4: free early)
    hts_idx_destroy(vcf_idx);
    vcf_idx = NULL;
    
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    
    fprintf(stderr, "\nPass 1 complete: %ld SNPs in %ld seconds\n", total_snps.load(), duration.count());
    if (wants_generic_allocation) {
        fprintf(stderr, "Found %lu unique clade patterns\n", branchcounts.size());
    } else {
        fprintf(stderr, "Skipped generic clade accounting (native ATAC-only outputs)\n");
    }
    if (collect_het) {
        fprintf(stderr, "Found %ld het candidates\n", total_het_candidates.load());
    }
    if (collect_species) {
        fprintf(stderr, "Found %zu species candidates\n", all_species_candidates.size());
    }
    
    // B3: too-many-clades handling is implemented in the slot allocation section below,
    // where it allocates 1 slot per clade for the top num most-common clades.

    // ========================================================================
    // HET SITE SELECTION WITH BIN BALANCING
    // (Output writing happens after DEMUX_SCORE computation)
    // ========================================================================
    
    set<pair<string, int64_t>> selected_het_sites;  // (chrom, pos) of selected het sites
    vector<size_t> selected_het_indices;  // Indices into all_het_candidates
    
    if (!het_outfile.empty() && !all_het_candidates.empty()) {
        fprintf(stderr, "\n========================================\n");
        fprintf(stderr, "Selecting het sites with bin balancing...\n");
        start_time = chrono::high_resolution_clock::now();
        
        // Sort candidates by het_score descending with stable (chrom_idx, pos) tiebreaker
        fprintf(stderr, "  Sorting %zu het candidates by score...\n", all_het_candidates.size());
        sort(all_het_candidates.begin(), all_het_candidates.end(),
            [](const HetCandidate& a, const HetCandidate& b) {
                if (a.het_score != b.het_score) return a.het_score > b.het_score;
                if (a.chrom_idx != b.chrom_idx) return a.chrom_idx < b.chrom_idx;
                return a.pos < b.pos;
            });
        
        // Bin tracking: map from (chrom, bin_idx) -> count of selected SNPs
        map<pair<string, int>, int> bins_selected;
        
        // Greedy selection with bin balancing
        selected_het_indices.reserve(het_num);
        
        fprintf(stderr, "  Selecting up to %d het sites (max %d per bin)...\n",
            het_num, max_het_per_bin);

        for (size_t i = 0; i < all_het_candidates.size() && (int)selected_het_indices.size() < het_num; ++i) {
            const HetCandidate& hc = all_het_candidates[i];

            if (!hc.passes_rna_het_hwe) continue;

            // D9: apply RNA cov_score filter here (was previously applied at candidate
            // creation, which broke ATAC-het output if --cov was not provided).
            if (hc.cov_score < min_cov) continue;

            pair<string, int> bin_key = make_pair(hc.chrom, hc.bin_idx);

            // D7: hard bin cap (replaces dead evenness_score code)
            if (max_het_per_bin > 0 && bins_selected[bin_key] >= max_het_per_bin) continue;

            selected_het_indices.push_back(i);
            bins_selected[bin_key]++;
            selected_het_sites.insert(make_pair(hc.chrom, hc.pos));

            if ((selected_het_indices.size() % 100000) == 0) {
                fprintf(stderr, "    Selected %zu het sites...\r", selected_het_indices.size());
                fflush(stderr);
            }
        }
        
        end_time = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
        fprintf(stderr, "  Selected %zu het sites from %zu bins in %ld seconds\n", 
            selected_het_indices.size(), bins_selected.size(), duration.count());
        fprintf(stderr, "  (Het output will be written after DEMUX_SCORE computation)\n");
    }

    // ========================================================================
    // ATAC-HET SITE SELECTION (if --atac_het_output set)
    // Uses same het candidate pool but scored with ATAC coverage, no annot boost
    // ========================================================================
    
    vector<size_t> selected_atac_het_indices;
    set<pair<string, int64_t>> selected_atac_het_sites;
    
    if (!atac_het_outfile.empty() && !all_het_candidates.empty()) {
        fprintf(stderr, "\n========================================\n");
        fprintf(stderr, "Selecting ATAC-het sites...\n");
        start_time = chrono::high_resolution_clock::now();
        
        // Create a score vector for ATAC-het ranking with call-rate penalty
        // and weighted pool combination, consistent with RNA-het scoring (Item F).
        // ATAC-het score = final_hs * call_rate * log2(atac_cov + 1)
        vector<pair<float, size_t>> atac_het_scored; // (score, index)
        atac_het_scored.reserve(all_het_candidates.size());
        
        for (size_t i = 0; i < all_het_candidates.size(); i++) {
            HetCandidate& hc = all_het_candidates[i];
            if (panel_position_is_excluded(
                    atac_excluded_positions, hc.chrom_idx, hc.pos)) continue;
            if (hc.atac_cov_score < atac_min_cov) continue;
            if (atac_regulatory_required && !hc.atac_regulatory) continue;

            const StoredVariant& atac_het_variant =
                all_variants[hc.variant_idx];
            int roster_called = 0;
            int roster_het = 0;
            for (int sample_index : atac_roster_sample_indices) {
                const int dosage = pattern_dosage(
                    atac_het_variant.pattern, sample_index);
                if (dosage < 0) continue;
                roster_called++;
                if (dosage == 1) roster_het++;
            }
            if (roster_called <= 0 || roster_het <= 0) continue;
            const float roster_het_freq =
                (float)roster_het / roster_called;
            const float roster_call_rate =
                (float)roster_called / atac_roster_sample_indices.size();
            if (roster_call_rate < pairwise_rescue_min_call_rate) continue;
            
            float final_hs = roster_het_freq;
            if (enable_pool_scoring) {
                final_hs = (1.0f - het_pool_weight) * roster_het_freq
                         + het_pool_weight * hc.pool_discrim_norm;
            }
            float atac_cov_term = log2f(hc.atac_cov_score + 1.0f)
                                + (hc.atac_regulatory ? 0.1f : 0.0f);
            float atac_sel = final_hs * roster_call_rate * atac_cov_term;
            if (all_variants[hc.variant_idx].atac_soft_overlap) {
                atac_sel = apply_atac_soft_overlap_penalty(
                    atac_sel, atac_soft_overlap_multiplier);
            }
            hc.atac_het_sel_score = atac_sel;
            atac_het_scored.push_back(make_pair(atac_sel, i));
        }
        
        // Sort descending by score
        sort(atac_het_scored.begin(), atac_het_scored.end(),
            [](const pair<float, size_t>& a, const pair<float, size_t>& b) {
                if (a.first != b.first) return a.first > b.first;
                return a.second < b.second;
            });
        
        // Greedy selection with bin tracking and bin cap enforcement (Item G)
        map<pair<string, int>, int> atac_het_bins;
        selected_atac_het_indices.reserve(atac_het_num);
        
        for (auto& scored : atac_het_scored) {
            if ((int)selected_atac_het_indices.size() >= atac_het_num) break;
            size_t idx = scored.second;
            const HetCandidate& hc = all_het_candidates[idx];
            pair<string,int> bin_key = make_pair(hc.chrom, hc.bin_idx);
            if (max_het_per_bin > 0 && atac_het_bins[bin_key] >= max_het_per_bin) continue;
            
            selected_atac_het_indices.push_back(idx);
            selected_atac_het_sites.insert(make_pair(hc.chrom, hc.pos));
            atac_het_bins[bin_key]++;
        }
        
        end_time = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
        fprintf(stderr, "  Selected %zu ATAC-het sites from %zu bins in %ld seconds\n",
            selected_atac_het_indices.size(), atac_het_bins.size(), duration.count());
    }

    // ========================================================================
    // SPECIES SITE SELECTION (Items A, B, C in redesign)
    //
    // Replaces the previous per-pair-queue approach with a deficit-weighted
    // greedy global optimization. At each step the candidate maximizing
    //
    //   utility(c) = [ sum_pi max(0, deficit[pi]) * pair_discrim[c][pi]
    //                                            * pair_dosage_sep[c][pi] ]
    //              * coverage_term(c)
    //              * spacing_penalty(c)
    //
    // is selected, deficits are decremented for every pair the candidate
    // satisfies above species_min_pair_discrim, and the loop continues until
    // species_num candidates have been chosen or no more candidates exist.
    //
    // After per-pair deficits are exhausted, the utility falls back to a
    // pure coverage * dosage * max_pair_discrim "global utility" term so the
    // loop continues backfilling to species_num (Item A backfill).
    //
    // Pair-target allocation (Item B): when species_num >= n_pairs every
    // pair gets a Hamilton-apportioned share; when species_num < n_pairs
    // only the pairs with the most strong candidates receive nonzero targets.
    //
    // strong_threshold (Item C): min(1.0f, species_min_pair_discrim + 0.25f).
    // ========================================================================
    
    vector<size_t> selected_species_indices;
    unordered_set<size_t> selected_species_sites;  // keyed by variant_idx
    
    if (collect_species && !all_species_candidates.empty()) {
        fprintf(stderr, "\n========================================\n");
        fprintf(stderr, "Selecting species sites (two-pass sort-scan deficit/backfill)...\n");
        start_time = chrono::high_resolution_clock::now();
        
        int n_sp_pairs = panel_metadata.n_pairs;
        int n_sp = (int)panel_metadata.species_list.size();
        
        // ----- Item C: adaptive strong_threshold capped at 1.0 -----
        float strong_threshold = std::min(1.0f, species_min_pair_discrim + 0.25f);
        const float soft_species_threshold = species_min_pair_discrim * 0.5f;
        
        // Count strong candidates per pair for proportional target allocation
        vector<long> pair_strong_count(n_sp_pairs, 0);
        for (const auto& sc : all_species_candidates) {
            for (int pi = 0; pi < n_sp_pairs; pi++) {
                if (sc.pair_discrim[pi] >= strong_threshold) pair_strong_count[pi]++;
            }
        }
        long total_strong = 0;
        for (long v : pair_strong_count) total_strong += v;
        
        // If no pair has any candidate above strong_threshold, fall back to
        // the 90th-percentile pair_discrim per pair as a softer strong cutoff.
        if (total_strong == 0) {
            fprintf(stderr, "  No candidates above strong_threshold=%.4f; "
                "falling back to per-pair 90th-percentile pair_discrim.\n",
                strong_threshold);
            for (int pi = 0; pi < n_sp_pairs; pi++) {
                vector<float> vals;
                vals.reserve(all_species_candidates.size());
                for (const auto& sc : all_species_candidates) {
                    if (sc.pair_discrim[pi] > 0.0f) vals.push_back(sc.pair_discrim[pi]);
                }
                if (vals.empty()) continue;
                sort(vals.begin(), vals.end());
                size_t p90 = (size_t)(vals.size() * 0.9);
                if (p90 >= vals.size()) p90 = vals.size() - 1;
                float thr = vals[p90];
                for (const auto& sc : all_species_candidates) {
                    if (sc.pair_discrim[pi] >= thr) pair_strong_count[pi]++;
                }
                total_strong += pair_strong_count[pi];
            }
        }
        
        // ----- Item B: budget-feasible pair-target allocation -----
        vector<long> pair_target(n_sp_pairs, 0);
        if (species_num >= n_sp_pairs) {
            // Hamilton apportionment: every pair gets floor(species_num / n_pairs)
            // and the remainder is distributed by strong-count rank.
            long base = species_num / n_sp_pairs;
            long remainder = species_num - base * n_sp_pairs;
            for (int pi = 0; pi < n_sp_pairs; pi++) pair_target[pi] = base;
            
            // Distribute remainder to pairs with the most strong candidates
            vector<int> order(n_sp_pairs);
            for (int i = 0; i < n_sp_pairs; i++) order[i] = i;
            sort(order.begin(), order.end(), [&](int a, int b){
                if (pair_strong_count[a] != pair_strong_count[b])
                    return pair_strong_count[a] > pair_strong_count[b];
                return a < b;  // deterministic
            });
            for (long k = 0; k < remainder && k < n_sp_pairs; k++) {
                pair_target[order[(size_t)k]]++;
            }
        } else {
            // species_num < n_sp_pairs: targets cannot cover every pair.
            // Only the top species_num pairs (by strong count) get target=1;
            // the rest get target=0 and participate in backfill only.
            vector<int> order(n_sp_pairs);
            for (int i = 0; i < n_sp_pairs; i++) order[i] = i;
            sort(order.begin(), order.end(), [&](int a, int b){
                if (pair_strong_count[a] != pair_strong_count[b])
                    return pair_strong_count[a] > pair_strong_count[b];
                return a < b;
            });
            for (int k = 0; k < species_num; k++) {
                pair_target[order[(size_t)k]] = 1;
            }
            fprintf(stderr, "  NOTE: species_num (%d) < n_sp_pairs (%d); "
                "%d pair(s) start with target=0 and rely on multi-pair credit "
                "or backfill.\n", species_num, n_sp_pairs,
                n_sp_pairs - species_num);
        }
        
        // Sanity: targets must sum to <= species_num
        long target_sum = 0;
        for (long v : pair_target) target_sum += v;
        if (target_sum > species_num) {
            // Should not happen with the above allocation, but trim defensively.
            long excess = target_sum - species_num;
            for (int pi = n_sp_pairs - 1; pi >= 0 && excess > 0; pi--) {
                long take = std::min(excess, pair_target[pi]);
                pair_target[pi] -= take;
                excess -= take;
            }
        }
        
        // B9: Precompute pair confidence (singleton species downweighted)
        vector<float> species_confidence(n_sp, 1.0f);
        for (int si = 0; si < n_sp; si++) {
            const string& sp = panel_metadata.species_list[si];
            auto it = panel_metadata.species_to_sample_indices.find(sp);
            int n_in_sp = (it != panel_metadata.species_to_sample_indices.end())
                          ? (int)it->second.size() : 0;
            if (n_in_sp <= 1) {
                species_confidence[si] = 0.5f;
            } else {
                species_confidence[si] = std::min(1.0f, sqrtf((float)n_in_sp / 3.0f));
            }
        }
        vector<float> pair_confidence(n_sp_pairs, 1.0f);
        vector<char> pair_has_outgroup(n_sp_pairs, 0);
        for (int pi = 0; pi < n_sp_pairs; pi++) {
            const string& spa = panel_metadata.species_pairs[pi].first;
            const string& spb = panel_metadata.species_pairs[pi].second;
            if (!species_outgroup.empty() &&
                (spa == species_outgroup || spb == species_outgroup)) {
                pair_has_outgroup[pi] = 1;
            }
            int si_a_conf = -1, si_b_conf = -1;
            for (int si = 0; si < n_sp; si++) {
                if (panel_metadata.species_list[si] == spa) si_a_conf = si;
                if (panel_metadata.species_list[si] == spb) si_b_conf = si;
            }
            if (si_a_conf >= 0 && si_b_conf >= 0) {
                pair_confidence[pi] = species_confidence[si_a_conf] * species_confidence[si_b_conf];
            }
        }
        bool species_outgroup_cap_enabled = !species_outgroup.empty() &&
            species_outgroup_max_fraction < 1.0f;
        long species_outgroup_max_sites = species_outgroup_cap_enabled
            ? (long)floor((double)species_num * species_outgroup_max_fraction)
            : species_num;
        long species_outgroup_only_selected = 0;
        long species_non_outgroup_or_mixed_selected = 0;
        auto actual_outgroup_limit = [&]() -> long {
            if (!species_outgroup_cap_enabled) return species_num;
            if (species_outgroup_max_fraction <= 0.0f) return 0;
            return (long)floor(
                ((double)species_outgroup_max_fraction /
                 (1.0 - (double)species_outgroup_max_fraction)) *
                species_non_outgroup_or_mixed_selected);
        };
        auto is_outgroup_only_candidate = [&](
            const SpeciesCandidate& sc, float threshold) -> bool {
            bool helps_outgroup = false;
            bool helps_non_outgroup = false;
            for (int pi = 0; pi < n_sp_pairs; pi++) {
                if (sc.pair_discrim[pi] < threshold) continue;
                if (pair_has_outgroup[pi]) helps_outgroup = true;
                else helps_non_outgroup = true;
            }
            return helps_outgroup && !helps_non_outgroup;
        };
        
        fprintf(stderr, "  %zu candidates, %d species, %d pairs, total budget %d\n",
            all_species_candidates.size(), n_sp, n_sp_pairs, species_num);
        fprintf(stderr, "  Strong threshold: %.4f (total strong candidates: %ld)\n",
            strong_threshold, total_strong);
        fprintf(stderr, "  Pair targets (budget-feasible, sum=%ld):\n", target_sum);
        for (int pi = 0; pi < n_sp_pairs; pi++) {
            fprintf(stderr, "    %s-%s: target=%ld (strong=%ld)\n",
                panel_metadata.species_pairs[pi].first.c_str(),
                panel_metadata.species_pairs[pi].second.c_str(),
                pair_target[pi], pair_strong_count[pi]);
        }
        fprintf(stderr, "  species_min_pair_discrim=%.4f, species_cov_cap=%.1f, "
            "species_max_per_bin=%d\n",
            species_min_pair_discrim, species_cov_cap, species_max_per_bin);
        if (species_outgroup_cap_enabled) {
            fprintf(stderr,
                "  Outgroup-only backfill cap: %s <= %.1f%% of the actual final set "
                "(%ld is only the requested-budget upper bound); sites that also "
                "distinguish a non-outgroup pair do not consume this cap\n",
                species_outgroup.c_str(),
                100.0f * species_outgroup_max_fraction,
                species_outgroup_max_sites);
        }
        
        // Tracking
        vector<long> pair_counts(n_sp_pairs, 0);              // satisfaction count per pair
        vector<long> pair_credited(n_sp_pairs, 0);            // total credited (audit)
        unordered_map<uint64_t, int> sp_bins_selected;        // key = (chrom_idx << 32) | bin_idx
        
        selected_species_indices.reserve(species_num);
        
        // ----- Sort-then-scan with two-pass deficit/backfill -----
        //
        // Replaces the O(species_num * n_candidates) full-rescan greedy loop
        // that became intractable at 54M candidates x 20M target.
        //
        // Structure:
        //   1. Compute a stable priority score per candidate.
        //   2. Sort by priority descending with deterministic tiebreaks.
        //   3. Pass A (deficit fill): scan sorted list once, accept only
        //      candidates that help at least one pair with positive deficit.
        //   4. Pass B (backfill): scan sorted list again, accept remaining
        //      best candidates regardless of deficit, up to species_num.
        //
        // Complexity: O(n log n + 2 * n * n_pairs), vs previous O(k * n * n_pairs)
        // where k = species_num.
        
        // Helper: pack (chrom_idx, bin_idx) into uint64_t key
        auto bin_key = [](int chrom_idx, int bin_idx) -> uint64_t {
            return ((uint64_t)(unsigned)chrom_idx << 32) | (uint64_t)(unsigned)bin_idx;
        };
        
        // Step 1: compute priority score for sorting.
        // priority = max_pi(pair_discrim[pi] * dosage_sep[pi] * pair_confidence[pi])
        //          * min(log2(cov_score + 1), species_cov_cap)
        // Non-outgroup contrasts sort ahead of outgroup-only contrasts when a
        // cap is active, preventing a singleton outgroup from dominating the
        // backfill before within-ingroups contrasts have been considered.
        fprintf(stderr, "  Computing priority scores for %zu candidates...\n",
            all_species_candidates.size());
        
        vector<float> priority_scores(all_species_candidates.size());
        vector<float> soft_priority_scores(all_species_candidates.size());
        for (size_t i = 0; i < all_species_candidates.size(); i++) {
            const SpeciesCandidate& sc = all_species_candidates[i];
            float best_contrib = 0.0f;
            float best_non_outgroup_contrib = 0.0f;
            float best_outgroup_contrib = 0.0f;
            float best_soft_contrib = 0.0f;
            float best_soft_non_outgroup_contrib = 0.0f;
            float best_soft_outgroup_contrib = 0.0f;
            for (int pi = 0; pi < n_sp_pairs; pi++) {
                float pd = sc.pair_discrim[pi];
                float ds = (pi < (int)sc.pair_dosage_sep.size())
                           ? sc.pair_dosage_sep[pi] : 0.0f;
                float v = pd * ds * pair_confidence[pi];
                if (pd >= species_min_pair_discrim) {
                    if (v > best_contrib) best_contrib = v;
                    if (pair_has_outgroup[pi]) {
                        if (v > best_outgroup_contrib) best_outgroup_contrib = v;
                    } else if (v > best_non_outgroup_contrib) {
                        best_non_outgroup_contrib = v;
                    }
                } else if (pd >= soft_species_threshold) {
                    if (v > best_soft_contrib) best_soft_contrib = v;
                    if (pair_has_outgroup[pi]) {
                        if (v > best_soft_outgroup_contrib)
                            best_soft_outgroup_contrib = v;
                    } else if (v > best_soft_non_outgroup_contrib) {
                        best_soft_non_outgroup_contrib = v;
                    }
                }
            }
            if (species_outgroup_cap_enabled) {
                best_contrib = best_non_outgroup_contrib > 0.0f
                    ? best_non_outgroup_contrib
                    : best_outgroup_contrib * 0.25f;
                best_soft_contrib = best_soft_non_outgroup_contrib > 0.0f
                    ? best_soft_non_outgroup_contrib
                    : best_soft_outgroup_contrib * 0.25f;
            }
            float cov_term = std::min(log2f(sc.cov_score + 1.0f), species_cov_cap);
            if (species_coverage_mode == "atac" && sc.atac_regulatory) {
                cov_term += 0.1f;
            }
            priority_scores[i] = best_contrib * cov_term;
            soft_priority_scores[i] = best_soft_contrib * cov_term;
            if (species_coverage_mode == "atac" &&
                all_variants[sc.variant_idx].atac_soft_overlap) {
                priority_scores[i] = apply_atac_soft_overlap_penalty(
                    priority_scores[i], atac_soft_overlap_multiplier);
                soft_priority_scores[i] = apply_atac_soft_overlap_penalty(
                    soft_priority_scores[i], atac_soft_overlap_multiplier);
            }
        }
        
        // Step 2: build index array and sort by priority descending.
        // Deterministic tiebreaks: priority desc, then chrom_idx asc, pos asc, variant_idx asc.
        fprintf(stderr, "  Sorting %zu species candidates by priority score...\n",
            all_species_candidates.size());
        
        vector<size_t> sorted_idx(all_species_candidates.size());
        for (size_t i = 0; i < sorted_idx.size(); i++) sorted_idx[i] = i;
        sort(sorted_idx.begin(), sorted_idx.end(), [&](size_t a, size_t b) {
            if (priority_scores[a] != priority_scores[b])
                return priority_scores[a] > priority_scores[b];
            const SpeciesCandidate& sa = all_species_candidates[a];
            const SpeciesCandidate& sb = all_species_candidates[b];
            if (sa.chrom_idx != sb.chrom_idx) return sa.chrom_idx < sb.chrom_idx;
            if (sa.pos != sb.pos) return sa.pos < sb.pos;
            return sa.variant_idx < sb.variant_idx;
        });
        
        // Pre-selection feasibility check: estimate max selectable sites
        // given species_max_per_bin constraint.
        if (species_max_per_bin > 0) {
            unordered_set<uint64_t> eligible_bins;
            for (size_t i = 0; i < all_species_candidates.size(); i++) {
                if (priority_scores[i] <= 0.0f) continue;
                const SpeciesCandidate& sc = all_species_candidates[i];
                eligible_bins.insert(bin_key(sc.chrom_idx, sc.bin_idx));
            }
            long max_possible = (long)eligible_bins.size() * species_max_per_bin;
            fprintf(stderr, "  Feasibility: %zu eligible bins x %d max_per_bin = %ld max possible\n",
                eligible_bins.size(), species_max_per_bin, max_possible);
            if (species_num > (int)max_possible) {
                fprintf(stderr, "  WARNING: species_num=%d exceeds bin-cap maximum %ld. "
                    "Selection will stop at %ld. Consider raising --species_max_per_bin.\n",
                    species_num, max_possible, max_possible);
            }
        }
        
        // Step 3: two-pass selection.
        fprintf(stderr, "  Pass A: deficit fill...\n");
        
        long blocked_by_bin = 0;
        long skipped_no_eligible = 0;
        vector<char> used_candidates(all_species_candidates.size(), 0);
        
        // Pass A: accept only candidates that help at least one pair with
        // positive remaining deficit. This ensures rare/weak pairs get
        // filled before the budget is consumed by already-satisfied pairs.
        for (size_t rank = 0; rank < sorted_idx.size() &&
             (int)selected_species_indices.size() < species_num; rank++) {
            
            size_t i = sorted_idx[rank];
            SpeciesCandidate& sc = all_species_candidates[i];
            
            if (priority_scores[i] <= 0.0f) break;
            
            // Duplicate-site guard (keyed by variant_idx)
            if (selected_species_sites.count(sc.variant_idx) > 0) {
                used_candidates[i] = 1;
                continue;
            }

            bool outgroup_only = is_outgroup_only_candidate(
                sc, soft_species_threshold);
            if (species_outgroup_cap_enabled && outgroup_only &&
                species_outgroup_only_selected >= actual_outgroup_limit()) {
                continue;
            }
            
            // Bin cap enforcement
            if (species_max_per_bin > 0) {
                uint64_t bk = bin_key(sc.chrom_idx, sc.bin_idx);
                auto it_bin = sp_bins_selected.find(bk);
                if (it_bin != sp_bins_selected.end() &&
                    it_bin->second >= species_max_per_bin) {
                    blocked_by_bin++;
                    continue;
                }
            }
            
            // Check if this candidate helps any pair with positive deficit.
            // Find the eligible pair with the largest remaining deficit.
            int best_pi = -1;
            long best_deficit = 0;
            float best_score_for_tie = -1.0f;
            bool helps_any_deficit = false;
            
            for (int pi = 0; pi < n_sp_pairs; pi++) {
                float pd = sc.pair_discrim[pi];
                if (pd < species_min_pair_discrim) continue;
                
                long d = pair_target[pi] - pair_counts[pi];
                if (d <= 0) continue;  // Pass A: skip satisfied pairs entirely
                
                helps_any_deficit = true;
                float tie_score = pd * pair_confidence[pi];
                
                if (d > best_deficit || (d == best_deficit && tie_score > best_score_for_tie)) {
                    best_deficit = d;
                    best_pi = pi;
                    best_score_for_tie = tie_score;
                }
            }
            
            if (!helps_any_deficit) continue;  // Pass A: defer to backfill
            
            // Accept candidate
            sc.assigned_pair = best_pi;
            used_candidates[i] = 1;
            selected_species_indices.push_back(i);
            selected_species_sites.insert(sc.variant_idx);
            sp_bins_selected[bin_key(sc.chrom_idx, sc.bin_idx)]++;
            if (outgroup_only) species_outgroup_only_selected++;
            else species_non_outgroup_or_mixed_selected++;
            
            // Multi-pair credit: increment for every pair above threshold
            for (int pj = 0; pj < n_sp_pairs; pj++) {
                if (sc.pair_discrim[pj] >= species_min_pair_discrim) {
                    pair_counts[pj]++;
                    pair_credited[pj]++;
                }
            }
            
            if ((selected_species_indices.size() % 1000000) == 0) {
                fprintf(stderr, "    [deficit] Selected %zu / %d species sites...\r",
                    selected_species_indices.size(), species_num);
                fflush(stderr);
            }
            
            // Check if all deficits are now satisfied
            bool all_satisfied = true;
            for (int pi = 0; pi < n_sp_pairs; pi++) {
                if (pair_counts[pi] < pair_target[pi]) { all_satisfied = false; break; }
            }
            if (all_satisfied) {
                fprintf(stderr, "    All pair deficits satisfied at %zu selected sites. "
                    "Switching to backfill.\n", selected_species_indices.size());
                break;
            }
        }
        
        long deficit_fill_count = (long)selected_species_indices.size();
        fprintf(stderr, "  Pass A complete: %ld sites selected (deficit fill)\n", deficit_fill_count);
        
        // Per-pair deficit status after pass A
        for (int pi = 0; pi < n_sp_pairs; pi++) {
            long remaining = pair_target[pi] - pair_counts[pi];
            if (remaining > 0) {
                fprintf(stderr, "    %s-%s: deficit remaining %ld / %ld\n",
                    panel_metadata.species_pairs[pi].first.c_str(),
                    panel_metadata.species_pairs[pi].second.c_str(),
                    remaining, pair_target[pi]);
            }
        }
        
        // Pass B: backfill. Scan the same sorted list again and accept
        // remaining best candidates regardless of deficit status.
        if ((int)selected_species_indices.size() < species_num) {
            fprintf(stderr, "  Pass B: backfill (%d remaining)...\n",
                species_num - (int)selected_species_indices.size());
            
            long backfill_blocked_by_bin = 0;
            
            for (size_t rank = 0; rank < sorted_idx.size() &&
                 (int)selected_species_indices.size() < species_num; rank++) {
                
                size_t i = sorted_idx[rank];
                if (used_candidates[i]) continue;
                
                SpeciesCandidate& sc = all_species_candidates[i];
                
                if (priority_scores[i] <= 0.0f) break;
                
                // Duplicate-site guard
                if (selected_species_sites.count(sc.variant_idx) > 0) {
                    used_candidates[i] = 1;
                    continue;
                }

                bool outgroup_only = is_outgroup_only_candidate(
                    sc, soft_species_threshold);
                if (species_outgroup_cap_enabled && outgroup_only &&
                    species_outgroup_only_selected >= actual_outgroup_limit()) {
                    continue;
                }
                
                // Bin cap enforcement
                if (species_max_per_bin > 0) {
                    uint64_t bk = bin_key(sc.chrom_idx, sc.bin_idx);
                    auto it_bin = sp_bins_selected.find(bk);
                    if (it_bin != sp_bins_selected.end() &&
                        it_bin->second >= species_max_per_bin) {
                        backfill_blocked_by_bin++;
                        continue;
                    }
                }
                
                // Find best eligible pair by discrim * confidence (no deficit filter)
                int best_pi = -1;
                float best_score = -1.0f;
                
                for (int pi = 0; pi < n_sp_pairs; pi++) {
                    float pd = sc.pair_discrim[pi];
                    if (pd < species_min_pair_discrim) continue;
                    float tie_score = pd * pair_confidence[pi];
                    if (tie_score > best_score) {
                        best_score = tie_score;
                        best_pi = pi;
                    }
                }
                
                if (best_pi < 0) {
                    skipped_no_eligible++;
                    continue;
                }
                
                // Accept candidate
                sc.assigned_pair = best_pi;
                used_candidates[i] = 1;
                selected_species_indices.push_back(i);
                selected_species_sites.insert(sc.variant_idx);
                sp_bins_selected[bin_key(sc.chrom_idx, sc.bin_idx)]++;
                if (outgroup_only) species_outgroup_only_selected++;
                else species_non_outgroup_or_mixed_selected++;
                
                for (int pj = 0; pj < n_sp_pairs; pj++) {
                    if (sc.pair_discrim[pj] >= species_min_pair_discrim) {
                        pair_counts[pj]++;
                        pair_credited[pj]++;
                    }
                }
                
                if ((selected_species_indices.size() % 1000000) == 0) {
                    fprintf(stderr, "    [backfill] Selected %zu / %d species sites...\r",
                        selected_species_indices.size(), species_num);
                    fflush(stderr);
                }
            }
            
            blocked_by_bin += backfill_blocked_by_bin;
            fprintf(stderr, "  Pass B complete: %zu sites selected (%zu backfill)\n",
                selected_species_indices.size(),
                selected_species_indices.size() - (size_t)deficit_fill_count);
        }

        // Pass C: only after every selectable strict candidate has been
        // considered, use the deliberately retained half-threshold candidates
        // as a transparent soft backfill.  These sites are assigned and
        // audited separately; they never masquerade as >=threshold credits.
        vector<long> pair_soft_credited(n_sp_pairs, 0);
        long soft_backfill_count = 0;
        if ((int)selected_species_indices.size() < species_num) {
            vector<size_t> soft_sorted_idx(all_species_candidates.size());
            for (size_t i = 0; i < soft_sorted_idx.size(); ++i) {
                soft_sorted_idx[i] = i;
            }
            sort(soft_sorted_idx.begin(), soft_sorted_idx.end(),
                [&](size_t a, size_t b) {
                    if (soft_priority_scores[a] != soft_priority_scores[b])
                        return soft_priority_scores[a] > soft_priority_scores[b];
                    const SpeciesCandidate& sa = all_species_candidates[a];
                    const SpeciesCandidate& sb = all_species_candidates[b];
                    if (sa.chrom_idx != sb.chrom_idx)
                        return sa.chrom_idx < sb.chrom_idx;
                    if (sa.pos != sb.pos) return sa.pos < sb.pos;
                    return sa.variant_idx < sb.variant_idx;
                });
            fprintf(stderr,
                "  Pass C: soft species backfill at pair_discrim >= %.4f "
                "(%d remaining)...\n",
                soft_species_threshold,
                species_num - (int)selected_species_indices.size());

            for (size_t rank = 0; rank < soft_sorted_idx.size() &&
                 (int)selected_species_indices.size() < species_num; ++rank) {
                size_t i = soft_sorted_idx[rank];
                if (used_candidates[i]) continue;
                if (soft_priority_scores[i] <= 0.0f) break;
                SpeciesCandidate& sc = all_species_candidates[i];
                if (selected_species_sites.count(sc.variant_idx) > 0) {
                    used_candidates[i] = 1;
                    continue;
                }
                bool outgroup_only = is_outgroup_only_candidate(
                    sc, soft_species_threshold);
                if (species_outgroup_cap_enabled && outgroup_only &&
                    species_outgroup_only_selected >=
                        actual_outgroup_limit()) {
                    continue;
                }
                if (species_max_per_bin > 0) {
                    uint64_t bk = bin_key(sc.chrom_idx, sc.bin_idx);
                    auto found_bin = sp_bins_selected.find(bk);
                    if (found_bin != sp_bins_selected.end() &&
                        found_bin->second >= species_max_per_bin) {
                        blocked_by_bin++;
                        continue;
                    }
                }

                int best_pi = -1;
                float best_score = -1.0f;
                for (int pi = 0; pi < n_sp_pairs; ++pi) {
                    float pd = sc.pair_discrim[pi];
                    if (pd < soft_species_threshold ||
                        pd >= species_min_pair_discrim) continue;
                    float tie_score = pd * pair_confidence[pi];
                    if (tie_score > best_score) {
                        best_score = tie_score;
                        best_pi = pi;
                    }
                }
                if (best_pi < 0) continue;

                sc.assigned_pair = best_pi;
                used_candidates[i] = 1;
                selected_species_indices.push_back(i);
                selected_species_sites.insert(sc.variant_idx);
                sp_bins_selected[bin_key(sc.chrom_idx, sc.bin_idx)]++;
                if (outgroup_only) species_outgroup_only_selected++;
                else species_non_outgroup_or_mixed_selected++;
                for (int pj = 0; pj < n_sp_pairs; ++pj) {
                    if (sc.pair_discrim[pj] >= soft_species_threshold &&
                        sc.pair_discrim[pj] < species_min_pair_discrim) {
                        pair_soft_credited[pj]++;
                    }
                }
                soft_backfill_count++;
            }
            fprintf(stderr,
                "  Pass C complete: %ld soft sites; total %zu / %d\n",
                soft_backfill_count, selected_species_indices.size(),
                species_num);
            for (int pi = 0; pi < n_sp_pairs; ++pi) {
                if (pair_soft_credited[pi] <= 0) continue;
                fprintf(stderr, "    %s-%s: soft credits=%ld\n",
                    panel_metadata.species_pairs[pi].first.c_str(),
                    panel_metadata.species_pairs[pi].second.c_str(),
                    pair_soft_credited[pi]);
            }
        }
        
        if (blocked_by_bin > 0) {
            fprintf(stderr, "  %ld candidates blocked by species_max_per_bin=%d\n",
                blocked_by_bin, species_max_per_bin);
        }
        if (skipped_no_eligible > 0) {
            fprintf(stderr, "  %ld candidates skipped (no pair above species_min_pair_discrim=%.4f)\n",
                skipped_no_eligible, species_min_pair_discrim);
        }
        if (species_outgroup_cap_enabled) {
            fprintf(stderr,
                "  Outgroup-only sites selected for %s: %ld / %zu final "
                "(%.3f%%; actual-set limit=%ld)\n",
                species_outgroup.c_str(), species_outgroup_only_selected,
                selected_species_indices.size(),
                selected_species_indices.empty() ? 0.0 :
                    100.0 * species_outgroup_only_selected /
                    selected_species_indices.size(),
                actual_outgroup_limit());
        }
        
        end_time = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
        fprintf(stderr, "  Selected %zu species sites from %zu bins in %ld seconds\n",
            selected_species_indices.size(), sp_bins_selected.size(), duration.count());
        
        if ((int)selected_species_indices.size() < species_num) {
            fprintf(stderr, "  WARNING: selected only %zu / %d species SNPs.\n",
                selected_species_indices.size(), species_num);
        }
        
        // Per-pair summary (both assigned and credited counts; Item from review section H)
        fprintf(stderr, "\n  Per-pair counts (assigned / credited / target):\n");
        // Recompute assigned counts (since we may have changed assignment logic)
        vector<long> pair_assigned(n_sp_pairs, 0);
        for (size_t idx : selected_species_indices) {
            int ap = all_species_candidates[idx].assigned_pair;
            if (ap >= 0 && ap < n_sp_pairs) pair_assigned[ap]++;
        }
        for (int pi = 0; pi < n_sp_pairs; pi++) {
            float fill_pct = pair_target[pi] > 0
                ? 100.0f * pair_credited[pi] / pair_target[pi] : 0.0f;
            fprintf(stderr, "    %s-%s: assigned=%ld credited=%ld target=%ld (%.1f%%)\n",
                panel_metadata.species_pairs[pi].first.c_str(),
                panel_metadata.species_pairs[pi].second.c_str(),
                pair_assigned[pi], pair_credited[pi], pair_target[pi], fill_pct);
        }
        
        // Post-selection species-balance audit (counts BOTH assigned and credited)
        fprintf(stderr, "\n  Species-balance audit:\n");
        map<string, long> species_site_counts_assigned;
        map<string, long> species_site_counts_credited;
        for (const string& sp : panel_metadata.species_list) {
            species_site_counts_assigned[sp] = 0;
            species_site_counts_credited[sp] = 0;
        }
        for (size_t idx : selected_species_indices) {
            const SpeciesCandidate& sc = all_species_candidates[idx];
            if (sc.assigned_pair >= 0) {
                species_site_counts_assigned[panel_metadata.species_pairs[sc.assigned_pair].first]++;
                species_site_counts_assigned[panel_metadata.species_pairs[sc.assigned_pair].second]++;
            }
            for (int pi = 0; pi < n_sp_pairs; pi++) {
                if (sc.pair_discrim[pi] >= species_min_pair_discrim) {
                    species_site_counts_credited[panel_metadata.species_pairs[pi].first]++;
                    species_site_counts_credited[panel_metadata.species_pairs[pi].second]++;
                }
            }
        }
        for (const string& sp : panel_metadata.species_list) {
            float pct_a = selected_species_indices.empty() ? 0.0f :
                100.0f * species_site_counts_assigned[sp] / selected_species_indices.size();
            float pct_c = selected_species_indices.empty() ? 0.0f :
                100.0f * species_site_counts_credited[sp] / selected_species_indices.size();
            fprintf(stderr, "    %s: %ld assigned (%.1f%%) / %ld credited (%.1f%%)\n",
                sp.c_str(),
                species_site_counts_assigned[sp], pct_a,
                species_site_counts_credited[sp], pct_c);
        }
        
        // Warn if any pair (with nonzero target) is severely underfilled
        for (int pi = 0; pi < n_sp_pairs; pi++) {
            if (pair_target[pi] > 0 && pair_credited[pi] < (long)(pair_target[pi] * 0.5)) {
                fprintf(stderr, "  WARNING: %s-%s is below 50%% of target (credited %ld / %ld)\n",
                    panel_metadata.species_pairs[pi].first.c_str(),
                    panel_metadata.species_pairs[pi].second.c_str(),
                    pair_credited[pi], pair_target[pi]);
            }
        }
        
        // Top 10 most-dense bins
        if (!sp_bins_selected.empty()) {
            vector<pair<int, uint64_t>> bin_ranked;
            bin_ranked.reserve(sp_bins_selected.size());
            for (auto& kv : sp_bins_selected) {
                bin_ranked.push_back(make_pair(kv.second, kv.first));
            }
            sort(bin_ranked.rbegin(), bin_ranked.rend());
            fprintf(stderr, "\n  Top 10 densest bins:\n");
            for (int i = 0; i < min(10, (int)bin_ranked.size()); i++) {
                int ci = (int)(bin_ranked[i].second >> 32);
                int bi = (int)(bin_ranked[i].second & 0xFFFFFFFF);
                const string& cname = (ci >= 0 && ci < (int)chroms.size())
                                      ? chroms[ci] : string("?");
                fprintf(stderr, "    %s:%d -> %d sites\n",
                    cname.c_str(), bi, bin_ranked[i].first);
            }
        }

        if (!atac_audit_prefix.empty()) {
            const string species_audit_path =
                atac_audit_prefix + ".species_tiers.tsv";
            ofstream species_audit(species_audit_path.c_str());
            if (!species_audit) {
                fprintf(stderr,
                    "ERROR: Could not write ATAC species-tier audit: %s\n",
                    species_audit_path.c_str());
                exit(1);
            }
            species_audit <<
                "species_a\tspecies_b\toutgroup_pair\ttier\t"
                "min_pair_discrim\tmax_pair_discrim\tcandidate_sites\t"
                "soft_overlap_sites\tselected_sites\n";
            vector<array<long, 2>> tier_candidates(n_sp_pairs, {0, 0});
            vector<array<long, 2>> tier_overlap(n_sp_pairs, {0, 0});
            vector<array<long, 2>> tier_selected(n_sp_pairs, {0, 0});
            long outgroup_only_candidates = 0;
            long outgroup_only_overlap = 0;
            for (size_t candidate_index = 0;
                 candidate_index < all_species_candidates.size();
                 ++candidate_index) {
                const SpeciesCandidate& candidate =
                    all_species_candidates[candidate_index];
                const bool overlap =
                    all_variants[candidate.variant_idx].atac_soft_overlap;
                const bool selected =
                    selected_species_sites.count(candidate.variant_idx) > 0;
                for (int pi = 0; pi < n_sp_pairs; ++pi) {
                    int tier = -1;
                    if (candidate.pair_discrim[pi] >=
                        species_min_pair_discrim) tier = 0;
                    else if (candidate.pair_discrim[pi] >=
                             soft_species_threshold) tier = 1;
                    if (tier < 0) continue;
                    tier_candidates[pi][tier]++;
                    if (overlap) tier_overlap[pi][tier]++;
                    if (selected) tier_selected[pi][tier]++;
                }
                if (is_outgroup_only_candidate(
                        candidate, soft_species_threshold)) {
                    outgroup_only_candidates++;
                    if (overlap) outgroup_only_overlap++;
                }
            }
            for (int pi = 0; pi < n_sp_pairs; ++pi) {
                for (int tier = 0; tier < 2; ++tier) {
                    species_audit
                        << panel_metadata.species_pairs[pi].first << '\t'
                        << panel_metadata.species_pairs[pi].second << '\t'
                        << (pair_has_outgroup[pi] ? 1 : 0) << '\t'
                        << (tier == 0 ? "strict" : "soft") << '\t'
                        << (tier == 0 ? species_min_pair_discrim
                                      : soft_species_threshold) << '\t'
                        << (tier == 0 ? 1.0f
                                      : species_min_pair_discrim) << '\t'
                        << tier_candidates[pi][tier] << '\t'
                        << tier_overlap[pi][tier] << '\t'
                        << tier_selected[pi][tier] << '\n';
                }
            }
            species_audit
                << "__SUMMARY__\t__SUMMARY__\t0\tACTUAL_FINAL_TOTAL\t0\t1\t"
                << all_species_candidates.size() << "\t0\t"
                << selected_species_indices.size() << '\n';
            species_audit
                << (species_outgroup.empty() ? "O" : species_outgroup)
                << "\t__ALL__\t1\tOUTGROUP_ONLY_ACTUAL\t"
                << soft_species_threshold << "\t1\t"
                << outgroup_only_candidates << '\t'
                << outgroup_only_overlap << '\t'
                << species_outgroup_only_selected << '\n';
            species_audit.close();
            fprintf(stderr, "  Species-tier audit: %s\n",
                species_audit_path.c_str());
        }
    }

    // ========================================================================
    // Process missing data and compute downsampling probabilities
    // (Same as original - this is fast)
    // ========================================================================

    // Generic/RNA allocation state remains available to the RNA writers below,
    // but the entire block is bypassed for native ATAC-only builds.
    unordered_map<bitset<NBITS>, long> clade_slot_allocation;
    unordered_map<bitset<NBITS>, double> downsample_prob;
    unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, double> downsample_miss_prob;
    double brent_exponent = 1.0;
    vector<size_t> selected_variant_indices;
    vector<uint8_t> selected_variant_mask;

    if (wants_generic_allocation) {
    selected_variant_indices.reserve(
        min((size_t)max(0, num), all_variants.size()));
    selected_variant_mask.assign(all_variants.size(), 0);
    auto select_generic_variant = [&](size_t variant_idx) -> bool {
        if (selected_variant_mask[variant_idx]) return false;
        selected_variant_mask[variant_idx] = 1;
        selected_variant_indices.push_back(variant_idx);
        return true;
    };
    
    // Under dosage encoding, index by total_alt instead of popcount.
    // A missing-data clade with total_alt = T and m missing samples can be
    // a parent of any fully-genotyped clade with total_alt in [T, T + 2m].
    map<int, vector<bitset<NBITS>>> total_alt2clades;
    for (auto& bc : branchcounts){
        int ta = compute_total_alt(bc.first, num_samples);
        total_alt2clades[ta].push_back(bc.first);
    }
    
    // NOTE: Missing-genotype sites contribute to clade SLOT ALLOCATION through
    // branchcounts_missing (the count_branch_missing pass distributes fractional
    // weight across compatible parent clades), but are skipped during the per-clade
    // candidate SELECTION pass below (only fully-genotyped sites enter
    // clade_to_variants). This is intentional: allocation uses fractional weights
    // to avoid undercounting common patterns, while selection requires unambiguous
    // assignment. The asymmetry produces shortfalls when missing-data sites are
    // abundant, which the redistribution step below corrects by shifting unused
    // slots to clades with surplus candidates.
    fprintf(stderr, "Adjusting counts using data from sites with missing genotypes (%zu entries)...\n", branchcounts_missing.size());
    fflush(stderr);
    
    unordered_map<bitset<NBITS>, float> to_add;
    unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, unordered_map<bitset<NBITS>, float>> miss_clade_probs;
    
    long adj_count = 0;
    for (auto& bcm : branchcounts_missing){
        adj_count++;
        if (adj_count % 100000 == 0) {
            fprintf(stderr, "  Adjusting: %ld/%zu\r", adj_count, branchcounts_missing.size());
            fflush(stderr);
        }
        unordered_map<bitset<NBITS>, float> mx;
        miss_clade_probs.insert(make_pair(bcm.first, mx));
        double parent_tot = 0.0;

        // bcm.first.first = present bitset (1-bit-per-sample)
        // bcm.first.second = canonical clade (dosage-aware, 2-bits-per-sample)
        // Build observed_mask: for each present sample i, set bits 2i and 2i+1
        bitset<NBITS> observed_mask;
        int n_missing_samples = 0;
        for (int i = 0; i < num_samples; i++) {
            if (bcm.first.first[i]) {
                observed_mask.set(2*i);
                observed_mask.set(2*i + 1);
            } else {
                n_missing_samples++;
            }
        }
        bitset<NBITS> observed_clade = bcm.first.second & observed_mask;
        int ta_observed = compute_total_alt(bcm.first.second, num_samples);

        // A compatible parent P must match on observed positions.
        // Its total_alt can be in [ta_observed, ta_observed + 2*n_missing_samples].
        // Search total_alt2clades for each possible total_alt value.
        auto check_parent = [&](const bitset<NBITS>& v) {
            if ((v & observed_mask) == observed_clade) {
                double cladecount = branchcounts[v];
                parent_tot += cladecount;
                miss_clade_probs[bcm.first].insert(make_pair(v, (float)cladecount));
            }
        };

        for (int ta_try = ta_observed; ta_try <= ta_observed + 2 * n_missing_samples; ta_try++) {
            auto it = total_alt2clades.find(ta_try);
            if (it == total_alt2clades.end()) continue;
            for (auto& v : it->second) {
                check_parent(v);
            }
        }

        // Also check flipped polarity parents
        if (miss2flip.count(bcm.first.second) > 0){
            bitset<NBITS> flipped_clade = miss2flip[bcm.first.second];
            bitset<NBITS> observed_flipped = flipped_clade & observed_mask;
            int ta_flip_observed = compute_total_alt(flipped_clade, num_samples);

            for (int ta_try = ta_flip_observed; ta_try <= ta_flip_observed + 2 * n_missing_samples; ta_try++) {
                auto it = total_alt2clades.find(ta_try);
                if (it == total_alt2clades.end()) continue;
                for (auto& v : it->second) {
                    if ((v & observed_mask) == observed_flipped) {
                        double cladecount = branchcounts[v];
                        parent_tot += cladecount;
                        miss_clade_probs[bcm.first].insert(make_pair(v, (float)cladecount));
                    }
                }
            }
        }
        
        if (parent_tot > 0){
            for (auto& mcp2 : miss_clade_probs[bcm.first]){
                float frac = mcp2.second / parent_tot;
                if (to_add.count(mcp2.first) == 0){
                    to_add.insert(make_pair(mcp2.first, frac * bcm.second));
                }
                else{
                    to_add[mcp2.first] += frac * bcm.second;
                }
                miss_clade_probs[bcm.first][mcp2.first] = frac;
            }
        }
        else{
            miss_clade_probs.erase(bcm.first);
        }
    }
    
    for (auto& ta : to_add){
        branchcounts[ta.first] += ta.second;
    }
    to_add.clear();
    
    // ========================================================================
    // Compute slot allocations using Brent's method (Option C: Stratified Deterministic)
    // ========================================================================
    
    fprintf(stderr, "Determine number of sites to downsample...\n");
    fflush(stderr);

    double allbc2 = 0.0;
    vector<bitset<NBITS>> bs_idx;
    vector<pair<double, int>> clcountsort;  // (negative_count, index)
    
    // Build sorted list
    fprintf(stderr, "Building clade count list from %zu entries...\n", branchcounts.size());
    fflush(stderr);
    for (auto bc = branchcounts.begin(); bc != branchcounts.end(); ){
        allbc2 += bc->second;
        clcountsort.push_back(make_pair(-bc->second, bs_idx.size()));
        bs_idx.push_back(bc->first);
        branchcounts.erase(bc++);
    }
    fprintf(stderr, "  Built list with %zu entries\n", clcountsort.size());
    fflush(stderr);
    
    fprintf(stderr, "%f total clade counts\n", allbc2);
    fflush(stderr);
    
    // Track whether we're in keep-all mode (total variants <= target)
    bool keep_all = false;
    
    if (allbc2 > (double)num){
        fprintf(stderr, "Sorting %zu clade counts...\n", clcountsort.size());
        fflush(stderr);
        sort(clcountsort.begin(), clcountsort.end());
        fprintf(stderr, "Sort complete.\n");
        fflush(stderr);
        
        // Extract counts into a simple vector for fast access
        vector<double> counts(clcountsort.size());
        for (size_t i = 0; i < clcountsort.size(); i++) {
            counts[i] = -clcountsort[i].first;
        }
        
        // B3: Too many distinct clades for requested target.
        // Allocate 1 slot to the top num most-common clades and skip Brent optimization.
        if (clcountsort.size() >= (size_t)num) {
            fprintf(stderr,
                "WARNING: %zu distinct allele sharing patterns found, but %d SNPs requested.\n",
                clcountsort.size(), num);
            fprintf(stderr,
                "  Allocating 1 slot per clade for the top %d most-common clades; rare clades dropped.\n",
                num);

            for (int i = 0; i < num && i < (int)clcountsort.size(); ++i) {
                const bitset<NBITS>& key = bs_idx[clcountsort[i].second];
                clade_slot_allocation[key] = 1;
                downsample_prob[key] = 1.0 / (-clcountsort[i].first);
            }

            brent_exponent = 0.0;
        } else {
        // Normal Brent/bisection optimization path
        
        // Pre-compute cumulative sums
        fprintf(stderr, "Computing cumulative sums...\n");
        fflush(stderr);
        vector<double> cumsum(counts.size() + 1, 0.0);
        for (size_t i = 0; i < counts.size(); i++) {
            cumsum[i + 1] = cumsum[i] + counts[i];
        }
        
        // Find starting point where sum_lo < num
        double threshold = allbc2 - (double)num;
        int start_n_hi = 1;
        for (size_t i = 1; i <= counts.size(); i++) {
            if (cumsum[i] > threshold) {
                start_n_hi = i;
                break;
            }
        }
        fprintf(stderr, "Starting optimization at n_hi=%d (sum_hi=%.0f, need > %.0f)\n", 
            start_n_hi, cumsum[start_n_hi], threshold);
        fprintf(stderr, "Top 10 clade counts: ");
        for (int i = 0; i < min(10, (int)counts.size()); i++) {
            fprintf(stderr, "%.0f ", counts[i]);
        }
        fprintf(stderr, "\n");
        fflush(stderr);
        
        auto brent_start = chrono::high_resolution_clock::now();
        
        // Lambda to compute sum of counts[0..n-1]^x (parallelized)
        auto power_sum = [&counts](int n, double x) -> double {
            double sum = 0.0;
            #pragma omp parallel for reduction(+:sum)
            for (int i = 0; i < n; i++) {
                sum += pow(counts[i], x);
            }
            return sum;
        };
        
        struct RootResult {
            double result;
            double upper;
        };
        auto solve_root = [&](int n_hi, double target) -> RootResult {
            double x_lo = 0.0;
            double x_hi = 1.0;
            double result = 0.5;
            // Deliberately retain the legacy 20 bisection steps and returned
            // midpoint. Candidate pruning below only avoids roots that cannot
            // satisfy the legacy strict boundary test.
            for (int iter = 0; iter < 20; ++iter) {
                result = (x_lo + x_hi) / 2.0;
                const double ps = power_sum(n_hi, result);
                if (ps < target) x_lo = result;
                else x_hi = result;
            }
            RootResult rr = {result, x_hi};
            return rr;
        };

        // The root is monotone nondecreasing as the prefix grows. Solve the
        // full prefix once to obtain a conservative upper bound, then retain
        // only boundaries that could possibly pass pow(last,x) > next. This
        // replaces ~one million full-prefix roots with one cheap boundary scan
        // while preserving the first accepted legacy prefix exactly.
        const int full_n = (int)counts.size();
        const RootResult full_root = solve_root(full_n, (double)num);
        const double root_guard = ldexp(1.0, -20);
        const double exponent_upper = min(
            1.0, nextafter(full_root.upper + root_guard,
                           numeric_limits<double>::infinity()));

        vector<int> boundary_candidates;
        boundary_candidates.reserve(4096);
        for (int n_hi = start_n_hi; n_hi <= full_n; ++n_hi) {
            const double sum_hi = cumsum[n_hi];
            const double sum_lo = allbc2 - sum_hi;
            const double target = (double)num - sum_lo;
            if (target < n_hi || target > sum_hi) continue;
            if (n_hi == full_n ||
                pow(counts[n_hi - 1], exponent_upper) > counts[n_hi]) {
                boundary_candidates.push_back(n_hi);
            }
        }
        fprintf(stderr,
            "  Exact boundary pruning: %d eligible prefixes -> %zu root solves\n",
            full_n - start_n_hi + 1, boundary_candidates.size());

        int roots_solved = 0;
        int final_n_hi = -1;
        for (int n_hi : boundary_candidates) {
            const double sum_hi = cumsum[n_hi];
            const double sum_lo = allbc2 - sum_hi;
            const double target = (double)num - sum_lo;
            const RootResult root = (n_hi == full_n)
                ? full_root : solve_root(n_hi, target);
            roots_solved++;
            const double res = root.result;
            const double lastcount = counts[n_hi - 1];
            const double nextcount = n_hi < full_n ? counts[n_hi] : -1.0;

            if (res >= 1.0 || pow(lastcount, res) <= nextcount) continue;

            auto brent_end = chrono::high_resolution_clock::now();
            auto brent_dur = chrono::duration_cast<chrono::seconds>(
                brent_end - brent_start);
            fprintf(stderr,
                "Found solution: x = %f at n_hi=%d (after %ld seconds, "
                "%d exact roots solved)\n",
                res, n_hi, brent_dur.count(), roots_solved);

            brent_exponent = res;
            final_n_hi = n_hi;
            fprintf(stderr, "Building slot allocation map for %d clades...\n", n_hi);

            vector<pair<double, int>> remainders;
            remainders.reserve(n_hi);
            vector<long> slots_vec(n_hi, 0);
            long sum_floor = 0;
            for (int i = 0; i < n_hi; ++i) {
                const double slots_float = pow(counts[i], res);
                const double slots_floor = floor(slots_float);
                long slots = (long)slots_floor;
                if (slots < 1) slots = 1;
                slots_vec[i] = slots;
                sum_floor += slots;
                remainders.push_back(make_pair(slots_float - slots_floor, i));
                downsample_prob.insert(make_pair(
                    bs_idx[clcountsort[i].second], slots_float / counts[i]));
            }

            const long target_long = (long)round(target);
            const long leftover = target_long - sum_floor;
            if (leftover > 0) {
                sort(remainders.begin(), remainders.end(),
                    [](const pair<double,int>& a, const pair<double,int>& b) {
                        return a.first > b.first;
                    });
                long added = 0;
                for (const auto& remainder : remainders) {
                    if (added >= leftover) break;
                    slots_vec[remainder.second] += 1;
                    added++;
                }
            }

            long total_slots_allocated = 0;
            for (int i = 0; i < n_hi; ++i) {
                clade_slot_allocation.insert(make_pair(
                    bs_idx[clcountsort[i].second], slots_vec[i]));
                total_slots_allocated += slots_vec[i];
            }
            fprintf(stderr,
                "  Allocated %ld initial slots across %zu clades "
                "(Hamilton remainder added %ld)\n",
                total_slots_allocated, clade_slot_allocation.size(),
                max(0L, leftover));
            break;
        }
        if (final_n_hi < 0) {
            fprintf(stderr, "ERROR: no valid clade-allocation boundary was found\n");
            exit(1);
        }
        } // end B3 else (normal Brent/bisection path)
    }
    else{
        fprintf(stderr, "%f total clade counts; allocating all sites (no downsampling).\n", allbc2);
        keep_all = true;
        for (size_t i = 0; i < bs_idx.size(); i++) {
            double count = -clcountsort[i].first;
            clade_slot_allocation[bs_idx[clcountsort[i].second]] = (long)count;
            downsample_prob[bs_idx[clcountsort[i].second]] = 1.0;
        }
        brent_exponent = 1.0;
    }
    
    // Compute downsample probabilities for missing genotype sites
    fprintf(stderr, "Computing downsample probabilities for missing genotype sites (%zu entries)...\n", miss_clade_probs.size());
    fflush(stderr);
    
    long miss_count = 0;
    for (auto& mcp_entry : miss_clade_probs){
        double p = 0.0;
        for (auto& mcp2 : mcp_entry.second){
            if (downsample_prob.count(mcp2.first) > 0){
                p += mcp2.second * downsample_prob[mcp2.first];
            }
            else{
                p += mcp2.second;
            }
        }
        if (p < 1.0){
            downsample_miss_prob.insert(make_pair(mcp_entry.first, p));
        }
        miss_count++;
        if (miss_count % 100000 == 0) {
            fprintf(stderr, "  Processed %ld missing entries...\r", miss_count);
            fflush(stderr);
        }
    }
    fprintf(stderr, "  Computed %zu missing genotype probabilities\n", downsample_miss_prob.size());
    fflush(stderr);
    
    if (clade_slot_allocation.size() == 0){
        fprintf(stderr, "No clades were found with more than %d occurrences.\n", num);
        fprintf(stderr, "There is nothing to downsample.\n");
        return 0;
    }

    // ========================================================================
    // STRATIFIED SELECTION (Option C): Select top SNPs within each clade
    // ========================================================================
    
    fprintf(stderr, "\n========================================\n");
    fprintf(stderr, "Stratified Selection: Selecting top SNPs within each clade...\n");
    start_time = chrono::high_resolution_clock::now();
    
    // Map from clade bitset -> vector of (variant_idx, selection_score)
    unordered_map<bitset<NBITS>, vector<CladeCandidate>> clade_to_variants;
    
    if (keep_all) {
        // Fix 4: True keep-all path. Total eligible variants <= target, so bypass
        // clade-based stratified selection entirely. Insert every QC-passing,
        // coverage-passing, non-blacklisted biallelic variant.
        fprintf(stderr, "  keep_all mode: inserting all eligible variants directly...\n");
        long keep_all_count = 0;
        for (size_t vi = 0; vi < all_variants.size(); vi++) {
            const StoredVariant& sv = all_variants[vi];
            if (!sv.is_biallelic || !sv.passes_qc) continue;
            
            float cov_score = sv.rna_cov_score;
            if (cov_score < min_cov) continue;
            if (sv.nuclear_exclusion != NUCLEAR_EXCLUDE_NONE) continue;
            
            select_generic_variant(vi);
            keep_all_count++;
        }
        fprintf(stderr, "  Inserted %ld variants in keep_all mode\n", keep_all_count);
    } else {
    // Normal stratified selection path
    
    fprintf(stderr, "  Grouping %zu variants by clade...\n", all_variants.size());
    fflush(stderr);
    
    long grouped_count = 0;
    long skipped_low_cov = 0;
    
    for (size_t vi = 0; vi < all_variants.size(); vi++) {
        const StoredVariant& sv = all_variants[vi];
        
        // Only consider biallelic SNPs that pass QC
        if (!sv.is_biallelic || !sv.passes_qc) continue;
        
        const string& chrom = chroms[sv.chrom_idx];
        const float cov_score = sv.rna_cov_score;
        const float annot_score = sv.annot_score;
        
        // Apply min_cov hard filter
        if (cov_score < min_cov) {
            skipped_low_cov++;
            continue;
        }
        
        // B11: Blacklist filter
        if (sv.nuclear_exclusion != NUCLEAR_EXCLUDE_NONE) continue;
        
        // Determine the canonical clade bitset (dosage-aware: A8)
        int ta = sv.total_alt;
        int ta_flip = sv.total_alt_flip;
        bitset<NBITS> clade_key;
        
        if (sv.pattern->present.count() == (size_t)num_samples) {
            // No missing data - use direct clade
            if (ta < ta_flip || (ta == ta_flip && sv.pattern->alt < sv.pattern->alt_flip)) {
                clade_key = sv.pattern->alt;
            } else {
                clade_key = sv.pattern->alt_flip;
            }
        } else {
            // Has missing data: assign to most-likely parent clade (B5)
            bitset<NBITS> observed_mask_sel;
            for (int ii = 0; ii < num_samples; ii++) {
                if (sv.pattern->present[ii]) {
                    observed_mask_sel.set(2*ii);
                    observed_mask_sel.set(2*ii + 1);
                }
            }
            int ta_miss = sv.total_alt;
            int ta_miss_flip = sv.total_alt_flip;
            bitset<NBITS> miss_clade_canonical;
            if (ta_miss < ta_miss_flip ||
                (ta_miss == ta_miss_flip && sv.pattern->alt < sv.pattern->alt_flip)) {
                miss_clade_canonical = sv.pattern->alt;
            } else {
                miss_clade_canonical = sv.pattern->alt_flip;
            }
            auto miss_key = make_pair(sv.pattern->present, miss_clade_canonical);
            auto mcp_it = miss_clade_probs.find(miss_key);
            if (mcp_it == miss_clade_probs.end()) continue;  // no compatible parent
            // Find max-weight parent
            bitset<NBITS> best_parent;
            float best_weight = -1.0f;
            for (auto& pp : mcp_it->second) {
                if (pp.second > best_weight) {
                    best_weight = pp.second;
                    best_parent = pp.first;
                }
            }
            if (clade_slot_allocation.count(best_parent) > 0) {
                CladeCandidate cc;
                cc.variant_idx = vi;
                float present_fraction = (float)sv.pattern->present.count() / num_samples;
                cc.selection_score = compute_selection_score(cov_score, annot_score)
                                    * present_fraction;  // missingness penalty
                clade_to_variants[best_parent].push_back(cc);
                grouped_count++;
            }
            continue;  // Fix 6: do not fall through to fully-genotyped clade block
        }
        
        // Only add if this clade has a slot allocation
        if (clade_slot_allocation.count(clade_key) > 0) {
            CladeCandidate cc;
            cc.variant_idx = vi;
            cc.selection_score = compute_selection_score(cov_score, annot_score);
            clade_to_variants[clade_key].push_back(cc);
            grouped_count++;
        }
        
        if (grouped_count % 5000000 == 0 && grouped_count > 0) {
            fprintf(stderr, "    Grouped %ld variants...\r", grouped_count);
            fflush(stderr);
        }
    }
    
    fprintf(stderr, "  Grouped %ld variants into %zu clades (skipped %ld below min_cov=%.1f)\n", 
        grouped_count, clade_to_variants.size(), skipped_low_cov, min_cov);
    fflush(stderr);
    
    // Step 2: Sort each clade's variants by selection score and select top N
    fprintf(stderr, "  Selecting top variants within each clade...\n");
    fflush(stderr);
    
    long total_shortfall = 0;
    long clades_with_shortfall = 0;
    vector<pair<bitset<NBITS>, long>> clades_with_surplus;  // (clade, extra_capacity)
    
    for (auto& kv : clade_to_variants) {
        const bitset<NBITS>& clade = kv.first;
        vector<CladeCandidate>& candidates = kv.second;
        
        long slots = clade_slot_allocation[clade];
        long available = candidates.size();
        
        // Sort by selection score (descending)
        sort(candidates.begin(), candidates.end());
        
        // Select top min(slots, available)
        long to_select = min(slots, available);
        for (long i = 0; i < to_select; i++) {
            select_generic_variant(candidates[i].variant_idx);
        }
        
        if (available < slots) {
            // Shortfall: couldn't fill all allocated slots
            long shortfall = slots - available;
            total_shortfall += shortfall;
            clades_with_shortfall++;
        } else if (available > slots) {
            // This clade has extra capacity for redistribution
            clades_with_surplus.push_back(make_pair(clade, available - slots));
        }
    }
    
    fprintf(stderr, "  Initial selection: %zu variants\n", selected_variant_indices.size());
    fprintf(stderr, "  Shortfall: %ld slots from %ld clades\n", total_shortfall, clades_with_shortfall);
    fflush(stderr);
    
    // Step 3: Redistribute shortfall proportionally to clades with extra capacity
    if (total_shortfall > 0 && !clades_with_surplus.empty()) {
        fprintf(stderr, "  Redistributing %ld shortfall slots...\n", total_shortfall);
        fflush(stderr);
        
        // Calculate total extra capacity
        long total_extra_capacity = 0;
        for (auto& cn : clades_with_surplus) {
            total_extra_capacity += cn.second;
        }
        
        // Redistribute proportionally
        long redistributed = 0;
        for (auto& cn : clades_with_surplus) {
            const bitset<NBITS>& clade = cn.first;
            long extra_capacity = cn.second;
            
            // Proportional share of shortfall
            long extra_slots = (long)round((double)total_shortfall * extra_capacity / total_extra_capacity);
            if (extra_slots > extra_capacity) extra_slots = extra_capacity;
            
            // Get this clade's candidates (already sorted)
            vector<CladeCandidate>& candidates = clade_to_variants[clade];
            long already_selected = clade_slot_allocation[clade];
            
            // Add more from this clade
            for (long i = already_selected; i < already_selected + extra_slots && i < (long)candidates.size(); i++) {
                if (select_generic_variant(candidates[i].variant_idx)) {
                    redistributed++;
                }
            }
        }
        
        fprintf(stderr, "  Redistributed %ld additional slots\n", redistributed);
        fprintf(stderr, "  Final selection: %zu variants\n", selected_variant_indices.size());
    }
    
    end_time = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    fprintf(stderr, "Stratified selection complete in %ld seconds\n", duration.count());
    fprintf(stderr, "  Target: %d, Selected: %zu, Shortfall redistributed: %ld\n", 
        num, selected_variant_indices.size(), total_shortfall);
    fflush(stderr);

    // A11: Write clade audit table
    {
        string audit_path = outfile + ".clade_audit.tsv";
        FILE* audit_fp = fopen(audit_path.c_str(), "w");
        if (audit_fp) {
            fprintf(audit_fp, "clade_idx\ttotal_alt\tn_called_samples\tslots_allocated\t"
                "candidates_available\tselected\tmean_selection_score\n");
            int clade_idx = 0;
            for (auto& kv : clade_to_variants) {
                const bitset<NBITS>& ck = kv.first;
                const vector<CladeCandidate>& cands = kv.second;
                long slots = 0;
                if (clade_slot_allocation.count(ck) > 0) slots = clade_slot_allocation[ck];
                long sel = min(slots, (long)cands.size());
                int ta_audit = compute_total_alt(ck, num_samples);
                // Count called samples from the clade key: a sample is "called" if
                // any of its two bits is set or it appears in a present mask. For
                // fully-genotyped clades, n_called = num_samples.
                int n_called_audit = num_samples; // approximation for fully-genotyped
                double mean_score = 0.0;
                for (long i = 0; i < sel; i++) mean_score += cands[i].selection_score;
                if (sel > 0) mean_score /= sel;
                fprintf(audit_fp, "%d\t%d\t%d\t%ld\t%zu\t%ld\t%.4f\n",
                    clade_idx, ta_audit, n_called_audit, slots, cands.size(), sel, mean_score);
                clade_idx++;
            }
            fclose(audit_fp);
            fprintf(stderr, "  Clade audit: %d clades written to %s\n", clade_idx, audit_path.c_str());
        } else {
            fprintf(stderr, "  WARNING: Could not open %s for clade audit\n", audit_path.c_str());
        }
    }

    } // end else (normal stratified selection, not keep_all)

    // ========================================================================
    // PAIRWISE FLOOR ENFORCEMENT (if --min_pairwise set)
    // After clade-balanced selection, check each individual pair's
    // distinguishing SNP count. For deficit pairs, add unused variants
    // that distinguish that pair, ranked by coverage score.
    // ========================================================================
    
    if (min_pairwise > 0) {
        fprintf(stderr, "\n========================================\n");
        fprintf(stderr, "Pairwise Floor Enforcement (target: %ld per pair)...\n", min_pairwise);
        start_time = chrono::high_resolution_clock::now();
        
        int ploidy = 2;
        int n_pairs_total = (num_samples * (num_samples - 1)) / 2;
        
        // Step 1: Compute pairwise distinguishing counts for selected variants
        fprintf(stderr, "  Computing pairwise counts for %zu selected variants across %d pairs...\n",
            selected_variant_indices.size(), n_pairs_total);
        
        // Use a flat array for pair counts: pair(i,j) where i<j -> index = i*num_samples - i*(i+1)/2 + j - i - 1
        vector<long> pair_counts(n_pairs_total, 0);
        
        auto pair_idx = [&](int i, int j) -> int {
            // i < j required
            return i * num_samples - i * (i + 1) / 2 + j - i - 1;
        };
        
        vector<vector<long>> thread_pair_counts(
            n_threads, vector<long>(n_pairs_total, 0));
        vector<long> thread_variants_scanned(n_threads, 0);
        #pragma omp parallel num_threads(n_threads)
        {
            const int thread_id = omp_get_thread_num();
            vector<int> ac(num_samples, -1);
            vector<long>& local_counts = thread_pair_counts[thread_id];
            #pragma omp for schedule(static)
            for (size_t selected_idx = 0;
                 selected_idx < selected_variant_indices.size(); ++selected_idx) {
                const StoredVariant& sv =
                    all_variants[selected_variant_indices[selected_idx]];
                if (sv.genotypes_packed.empty()) continue;
                fill(ac.begin(), ac.end(), -1);
                for (int s = 0; s < num_samples; ++s) {
                    const int a1 = gt_packed_allele(sv.genotypes_packed, s, 0);
                    const int a2 = gt_packed_allele(sv.genotypes_packed, s, 1);
                    if (a1 >= 0 && a2 >= 0) {
                        ac[s] = (a1 > 0 ? 1 : 0) + (a2 > 0 ? 1 : 0);
                    }
                }
                for (int i = 0; i < num_samples; ++i) {
                    if (ac[i] < 0) continue;
                    for (int j = i + 1; j < num_samples; ++j) {
                        if (ac[j] >= 0 && ac[i] != ac[j]) {
                            local_counts[pair_idx(i, j)]++;
                        }
                    }
                }
                thread_variants_scanned[thread_id]++;
            }
        }
        long variants_scanned = 0;
        for (int thread_id = 0; thread_id < n_threads; ++thread_id) {
            variants_scanned += thread_variants_scanned[thread_id];
            for (int pair = 0; pair < n_pairs_total; ++pair) {
                pair_counts[pair] += thread_pair_counts[thread_id][pair];
            }
        }
        fprintf(stderr, "    Scanned %ld selected variants for pairwise counts\n", variants_scanned);
        
        // Step 2: Identify deficit pairs
        vector<pair<int, int>> deficit_pairs; // (i, j) pairs below floor
        long worst_count = min_pairwise;
        int worst_i = -1, worst_j = -1;
        
        for (int i = 0; i < num_samples; i++) {
            for (int j = i + 1; j < num_samples; j++) {
                long cnt = pair_counts[pair_idx(i, j)];
                if (cnt < min_pairwise) {
                    deficit_pairs.push_back(make_pair(i, j));
                }
                if (cnt < worst_count) {
                    worst_count = cnt;
                    worst_i = i;
                    worst_j = j;
                }
            }
        }
        
        fprintf(stderr, "  %zu / %d pairs below floor of %ld\n", 
            deficit_pairs.size(), n_pairs_total, min_pairwise);
        if (worst_i >= 0) {
            fprintf(stderr, "  Worst pair: %s vs %s (%ld distinguishing SNPs)\n",
                samples[worst_i].c_str(), samples[worst_j].c_str(), worst_count);
        }
        
        if (!deficit_pairs.empty()) {
            // Step 3: For each deficit pair, find unused variants that distinguish them
            // Sort deficit pairs by count ascending (worst first)
            sort(deficit_pairs.begin(), deficit_pairs.end(),
                [&](const pair<int,int>& a, const pair<int,int>& b) {
                    return pair_counts[pair_idx(a.first, a.second)] < 
                           pair_counts[pair_idx(b.first, b.second)];
                });
            
            // Build a set of candidate variants not yet selected.
            // For each variant, store (rescue_score, variant_idx) per deficit pair where
            // rescue_score = coverage * num_deficit_pairs_helped * dosage_difference_at_this_pair.
            // (Item J): adds multi-pair utility and dosage strength to the scoring,
            // replacing pure-coverage ranking.
            fprintf(stderr, "  Scanning variant pool for pairwise rescue candidates...\n");
            
            // First pass: for each variant, determine which deficit pairs it distinguishes
            // and the dosage difference at each. Cache results so we don't recompute below.
            struct RescueVariantInfo {
                size_t vi;
                float cov;
                vector<pair<int, int>> distinguished_pairs;  // (pidx, dosage_diff)
            };
            vector<RescueVariantInfo> variant_info_cache;
            variant_info_cache.reserve(all_variants.size() / 10);  // estimate
            
            long candidates_found = 0;
            vector<int> rescue_ac(num_samples, -1);
            for (size_t vi = 0; vi < all_variants.size(); vi++) {
                if (selected_variant_mask[vi]) continue; // already selected
                
                const StoredVariant& sv = all_variants[vi];
                if (sv.nuclear_exclusion != NUCLEAR_EXCLUDE_NONE) continue;
                if (!sv.is_biallelic || !sv.passes_qc) continue;
                if (sv.genotypes_packed.empty()) continue;
                
                // Get coverage
                const float cov = sv.rna_cov_score;
                if (cov < min_cov) continue;
                
                // B8: call-rate floor for rescue candidates
                int n_called_rescue = 0;
                for (int s = 0; s < num_samples; s++) {
                    if (gt_packed_allele(sv.genotypes_packed, s, 0) >= 0 &&
                        gt_packed_allele(sv.genotypes_packed, s, 1) >= 0) {
                        n_called_rescue++;
                    }
                }
                float call_rate_rescue = (float)n_called_rescue / num_samples;
                if (call_rate_rescue < pairwise_rescue_min_call_rate) continue;
                
                // Compute alt_count
                fill(rescue_ac.begin(), rescue_ac.end(), -1);
                for (int s = 0; s < num_samples; s++) {
                    int a1 = gt_packed_allele(sv.genotypes_packed, s, 0);
                    int a2 = gt_packed_allele(sv.genotypes_packed, s, 1);
                    if (a1 >= 0 && a2 >= 0) {
                        rescue_ac[s] = (a1 > 0 ? 1 : 0) + (a2 > 0 ? 1 : 0);
                    }
                }
                
                // Check which deficit pairs this variant distinguishes
                RescueVariantInfo info;
                info.vi = vi;
                info.cov = cov;
                for (auto& dp : deficit_pairs) {
                    int a = dp.first, b = dp.second;
                    if (rescue_ac[a] >= 0 && rescue_ac[b] >= 0 &&
                        rescue_ac[a] != rescue_ac[b]) {
                        int pidx = pair_idx(a, b);
                        int dosage_diff = std::abs(
                            rescue_ac[a] - rescue_ac[b]);  // 1 or 2
                        info.distinguished_pairs.push_back(make_pair(pidx, dosage_diff));
                        candidates_found++;
                    }
                }
                if (!info.distinguished_pairs.empty()) {
                    variant_info_cache.push_back(std::move(info));
                }
                
                if (vi % 5000000 == 0 && vi > 0) {
                    fprintf(stderr, "    Scanned %zu / %zu variants, %ld candidates...\r",
                        vi, all_variants.size(), candidates_found);
                    fflush(stderr);
                }
            }
            fprintf(stderr, "    Found %ld rescue candidates across %zu deficit pairs\n",
                candidates_found, deficit_pairs.size());
            
            // For each deficit pair, build a scored candidate list:
            //   score = cov * num_deficit_pairs_helped * dosage_diff_at_this_pair
            map<int, vector<pair<float, size_t>>> pair_candidates;
            for (auto& dp : deficit_pairs) {
                pair_candidates[pair_idx(dp.first, dp.second)] = vector<pair<float, size_t>>();
            }
            for (const auto& info : variant_info_cache) {
                int n_helped = (int)info.distinguished_pairs.size();
                for (const auto& kv : info.distinguished_pairs) {
                    int pidx = kv.first;
                    int dosage_diff = kv.second;
                    float score = info.cov * (float)n_helped * (float)dosage_diff;
                    pair_candidates[pidx].push_back(make_pair(score, info.vi));
                }
            }
            
            // Step 4: For each deficit pair, sort candidates by rescue score descending,
            // add top-N needed to close the gap
            long total_added = 0;
            long pairs_rescued = 0;
            long pairs_partial = 0;
            map<pair<string, int>, int> rescue_bins;  // B8: bin tracking for rescue
            
            for (auto& dp : deficit_pairs) {
                // B8: global rescue cap
                if (total_added >= max_pairwise_rescue) break;
                
                int pidx = pair_idx(dp.first, dp.second);
                long current = pair_counts[pidx];
                long needed = min_pairwise - current;
                if (needed <= 0) continue; // already met (from additions for other pairs)
                
                auto& candidates = pair_candidates[pidx];
                
                // Sort by rescue score descending
                sort(candidates.begin(), candidates.end(),
                    [](const pair<float, size_t>& a, const pair<float, size_t>& b) {
                        return a.first > b.first;
                    });
                
                long added_this_pair = 0;
                for (auto& cand : candidates) {
                    if (added_this_pair >= needed) break;
                    if (total_added >= max_pairwise_rescue) break;
                    size_t vi = cand.second;
                    
                    // Check if already added by a previous deficit pair's rescue
                    if (selected_variant_mask[vi]) {
                        continue;
                    }
                    
                    // B8: bin-spacing constraint for rescue
                    if (pairwise_rescue_max_per_bin > 0) {
                        const StoredVariant& sv_check = all_variants[vi];
                        string rc = chroms[sv_check.chrom_idx];
                        int rb = get_bin_idx(sv_check.pos, bin_size);
                        auto rbk = make_pair(rc, rb);
                        if (rescue_bins[rbk] >= pairwise_rescue_max_per_bin) continue;
                    }
                    
                    select_generic_variant(vi);
                    added_this_pair++;
                    total_added++;
                    
                    // B8: update rescue bin tracking
                    {
                        const StoredVariant& sv_bin = all_variants[vi];
                        string rbc = chroms[sv_bin.chrom_idx];
                        int rbb = get_bin_idx(sv_bin.pos, bin_size);
                        rescue_bins[make_pair(rbc, rbb)]++;
                    }
                    
                    // Update pairwise counts for ALL pairs this new variant distinguishes
                    const StoredVariant& sv = all_variants[vi];
                    vector<int> ac(num_samples, -1);
                    for (int s = 0; s < num_samples; s++) {
                        int a1 = gt_packed_allele(sv.genotypes_packed, s, 0);
                        int a2 = gt_packed_allele(sv.genotypes_packed, s, 1);
                        if (a1 >= 0 && a2 >= 0) {
                            ac[s] = (a1 > 0 ? 1 : 0) + (a2 > 0 ? 1 : 0);
                        }
                    }
                    for (int i = 0; i < num_samples; i++) {
                        if (ac[i] < 0) continue;
                        for (int j = i + 1; j < num_samples; j++) {
                            if (ac[j] < 0) continue;
                            if (ac[i] != ac[j]) {
                                pair_counts[pair_idx(i, j)]++;
                            }
                        }
                    }
                }
                
                long final_count = pair_counts[pidx];
                if (final_count >= min_pairwise) {
                    pairs_rescued++;
                } else {
                    pairs_partial++;
                    fprintf(stderr, "  WARNING: %s vs %s: only reached %ld / %ld (pool exhausted)\n",
                        samples[dp.first].c_str(), samples[dp.second].c_str(), 
                        final_count, min_pairwise);
                }
            }
            
            fprintf(stderr, "  Added %ld rescue variants\n", total_added);
            fprintf(stderr, "  Pairs fully rescued: %ld, partially rescued: %ld\n", 
                pairs_rescued, pairs_partial);
            fprintf(stderr, "  New total selected: %zu (was %d target)\n", 
                selected_variant_indices.size(), num);
        }
        
        // Print final pairwise summary (bottom 10)
        vector<tuple<long, int, int>> all_pair_counts;
        for (int i = 0; i < num_samples; i++) {
            for (int j = i + 1; j < num_samples; j++) {
                all_pair_counts.push_back(make_tuple(pair_counts[pair_idx(i, j)], i, j));
            }
        }
        sort(all_pair_counts.begin(), all_pair_counts.end());
        
        fprintf(stderr, "\n  Bottom 10 pairwise counts (after rescue):\n");
        for (int k = 0; k < min(10, (int)all_pair_counts.size()); k++) {
            fprintf(stderr, "    %s vs %s: %ld\n",
                samples[get<1>(all_pair_counts[k])].c_str(),
                samples[get<2>(all_pair_counts[k])].c_str(),
                get<0>(all_pair_counts[k]));
        }
        
        end_time = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
        fprintf(stderr, "\nPairwise floor enforcement complete in %ld seconds\n", duration.count());
    }

    } // end generic/RNA allocation and selection

    // ========================================================================
    // ATAC-DEMUX STRATIFIED SELECTION (if --atac_output set)
    // Same clade allocation, ranked by ATAC coverage (no annotation boost)
    // ========================================================================
    
    vector<size_t> selected_atac_demux_indices;
    vector<uint8_t> selected_atac_demux_mask;
    
    if (!atac_outfile.empty() && !atac_coverage_map.empty()) {
        fprintf(stderr, "\n========================================\n");
        fprintf(stderr, "ATAC-Demux Stratified Selection...\n");
        start_time = chrono::high_resolution_clock::now();
        selected_atac_demux_indices.reserve(
            min((size_t)max(0, atac_num), all_variants.size()));
        selected_atac_demux_mask.assign(all_variants.size(), 0);
        auto select_atac_variant = [&](size_t variant_idx) -> bool {
            if (selected_atac_demux_mask[variant_idx]) return false;
            selected_atac_demux_mask[variant_idx] = 1;
            selected_atac_demux_indices.push_back(variant_idx);
            return true;
        };
        
        // Build the clade allocation from the actual ATAC-eligible pool. RNA
        // demux overlaps remain eligible with the configured soft rank penalty.
        // Reusing the generic/full-source allocation and then redistributing
        // its empty clades by raw surplus capacity caused singleton/outgroup clades
        // to absorb most of the 20M budget.
        typedef pair<bitset<NBITS>, bitset<NBITS>> NativeAtacCladeKey;
        unordered_map<NativeAtacCladeKey, long> atac_clade_slots;
        const long atac_base_target = min((long)atac_num,
                                          (long)all_variants.size());
        
        // Group variants by clade with ATAC coverage scoring
        unordered_map<NativeAtacCladeKey, vector<CladeCandidate>> atac_clade_variants;
        long atac_grouped = 0;
        
        for (size_t vi = 0; vi < all_variants.size(); vi++) {
            const StoredVariant& sv = all_variants[vi];
            if (sv.nuclear_exclusion != NUCLEAR_EXCLUDE_NONE) continue;
            if (!sv.is_biallelic || !sv.passes_qc) continue;
            if (panel_position_is_excluded(
                    atac_excluded_positions, sv.chrom_idx, sv.pos)) continue;
            const float roster_call_rate = operational_roster_call_rate(
                sv.pattern, atac_roster_sample_indices);
            if (roster_call_rate < pairwise_rescue_min_call_rate) continue;
            
            const string& chrom = chroms[sv.chrom_idx];
            const float atac_cov = sv.atac_cov_score;
            if (atac_cov < atac_min_cov) continue;
            if (atac_regulatory_required && !sv.atac_regulatory) continue;
            
            NativeAtacCladeKey clade_key = native_atac_clade_key(sv);
            
            CladeCandidate cc;
            cc.variant_idx = vi;
            cc.selection_score = log2f(atac_cov + 1.0f)
                               + (sv.atac_regulatory ? 0.1f : 0.0f);
            cc.selection_score *= roster_call_rate;
            if (sv.atac_soft_overlap) {
                cc.selection_score = apply_atac_soft_overlap_penalty(
                    cc.selection_score, atac_soft_overlap_multiplier);
            }
            atac_clade_variants[clade_key].push_back(cc);
            atac_grouped++;
        }
        
        fprintf(stderr, "  Grouped %ld variants into %zu clades for ATAC-demux\n",
            atac_grouped, atac_clade_variants.size());
        
        // Derive a power-balanced allocation from the actual ATAC-eligible
        // pool.  This preserves clade diversity without allowing the largest
        // singleton/outgroup clades to absorb every unfilled slot.
        struct AtacAllocationEntry {
            NativeAtacCladeKey clade_key;
            long capacity;
            long slots;
            double remainder;
            float pool_relevance;
            size_t representative_variant_idx;
        };
        vector<AtacAllocationEntry> atac_allocation;
        atac_allocation.reserve(atac_clade_variants.size());
        for (auto& kv : atac_clade_variants) {
            sort(kv.second.begin(), kv.second.end());
            AtacAllocationEntry entry;
            entry.clade_key = kv.first;
            entry.capacity = (long)kv.second.size();
            entry.slots = 0;
            entry.remainder = 0.0;
            entry.pool_relevance = 0.0f;
            entry.representative_variant_idx = kv.second.front().variant_idx;
            atac_allocation.push_back(entry);
        }

        // Match the working RNA semantics: pool information is one normalized
        // site-ranking score.  It is never an ATAC eligibility gate and never
        // creates per-identity or per-pair floor/rescue targets.
        #pragma omp parallel for schedule(static) num_threads(n_threads)
        for (size_t ai = 0; ai < atac_allocation.size(); ++ai) {
            AtacAllocationEntry& entry = atac_allocation[ai];
            const StoredVariant& representative = all_variants[
                entry.representative_variant_idx];
            entry.pool_relevance = compute_pool_discrim(
                representative.genotypes_packed, num_samples, pool_combos,
                sample_name_to_idx, min_pair_score);
        }
        fprintf(stderr,
            "  ATAC eligibility remains coverage/QC/nuclear-safety based: "
            "%ld variants in %zu clades; pool discrimination is ranking-only\n",
            atac_grouped, atac_allocation.size());

        long atac_expected = min((long)atac_num, atac_grouped);
        long atac_allocation_target = min(atac_expected, atac_base_target);
        double atac_power = 1.0;
        if (atac_grouped > atac_allocation_target && !atac_allocation.empty()) {
            double lo = 0.0, hi = 1.0;
            for (int iter = 0; iter < 80; iter++) {
                double mid = (lo + hi) / 2.0;
                double total = 0.0;
                for (const auto& entry : atac_allocation) {
                    total += min((double)entry.capacity,
                                 max(1.0, pow((double)entry.capacity, mid) *
                                          (1.0 + entry.pool_relevance)));
                }
                if (total < (double)atac_allocation_target) lo = mid;
                else hi = mid;
            }
            atac_power = (lo + hi) / 2.0;
        }

        long allocated = 0;
        for (auto& entry : atac_allocation) {
            double raw = (atac_grouped <= atac_allocation_target)
                ? (double)entry.capacity
                : min((double)entry.capacity,
                      max(1.0, pow((double)entry.capacity, atac_power) *
                               (1.0 + entry.pool_relevance)));
            entry.slots = min(entry.capacity, (long)floor(raw));
            entry.remainder = raw - (double)entry.slots;
            allocated += entry.slots;
        }
        sort(atac_allocation.begin(), atac_allocation.end(),
            [](const AtacAllocationEntry& a, const AtacAllocationEntry& b) {
                if (a.remainder != b.remainder) return a.remainder > b.remainder;
                if (a.pool_relevance != b.pool_relevance)
                    return a.pool_relevance > b.pool_relevance;
                if (a.capacity != b.capacity) return a.capacity < b.capacity;
                if (a.clade_key.first != b.clade_key.first)
                    return a.clade_key.first < b.clade_key.first;
                return a.clade_key.second < b.clade_key.second;
            });
        for (auto& entry : atac_allocation) {
            if (allocated >= atac_allocation_target) break;
            if (entry.slots < entry.capacity) {
                entry.slots++;
                allocated++;
            }
        }
        if (allocated != atac_allocation_target) {
            fprintf(stderr,
                "ERROR: ATAC-native allocation produced %ld / %ld slots\n",
                allocated, atac_allocation_target);
            exit(1);
        }
        for (const auto& entry : atac_allocation) {
            atac_clade_slots[entry.clade_key] = entry.slots;
            const vector<CladeCandidate>& candidates =
                atac_clade_variants.at(entry.clade_key);
            for (long i = 0; i < entry.slots; i++) {
                select_atac_variant(candidates[i].variant_idx);
            }
        }
        fprintf(stderr,
            "  ATAC-native clade allocation: %zu variants (power=%.6f; "
            "no identity-pair rescue reserve)\n",
            selected_atac_demux_indices.size(), atac_power);

        if ((long)selected_atac_demux_indices.size() != atac_expected) {
            fprintf(stderr,
                "ERROR: ATAC-demux target invariant failed: selected %zu, "
                "expected %ld from %ld eligible candidates\n",
                selected_atac_demux_indices.size(), atac_expected,
                atac_grouped);
            exit(1);
        }
        fprintf(stderr,
            "  Native ATAC identity-pair floors: disabled; pool information "
            "used only as the RNA-style scalar ranking score\n");

        if (!atac_audit_prefix.empty()) {
            vector<long> entry_soft_counts(atac_allocation.size(), 0);
            vector<long> entry_selected_counts(atac_allocation.size(), 0);
            vector<long> entry_selected_soft_counts(atac_allocation.size(), 0);
            for (size_t ai = 0; ai < atac_allocation.size(); ++ai) {
                const vector<CladeCandidate>& candidates =
                    atac_clade_variants.at(atac_allocation[ai].clade_key);
                for (const CladeCandidate& candidate : candidates) {
                    const bool overlap =
                        all_variants[candidate.variant_idx].atac_soft_overlap;
                    const bool selected =
                        selected_atac_demux_mask[candidate.variant_idx] != 0;
                    if (overlap) entry_soft_counts[ai]++;
                    if (selected) entry_selected_counts[ai]++;
                    if (selected && overlap)
                        entry_selected_soft_counts[ai]++;
                }
            }

            const string candidate_audit_path =
                atac_audit_prefix + ".candidate_sfs.tsv";
            ofstream candidate_audit(candidate_audit_path.c_str());
            if (!candidate_audit) {
                fprintf(stderr,
                    "ERROR: Could not write ATAC candidate/SFS audit: %s\n",
                    candidate_audit_path.c_str());
                exit(1);
            }
            candidate_audit <<
                "section\tmetric\tn_sites\tsoft_overlap_sites\n";
            const int n_gate_metrics = 7;
            const char* gate_names[n_gate_metrics] = {
                "source_loaded", "biallelic_with_genotypes",
                "nuclear_safe", "atac_covered", "distinct_called_dosages",
                "legacy_hard_exclusion_safe", "roster_call_rate"
            };
            long gate_counts[n_gate_metrics] = {0,0,0,0,0,0,0};
            long gate_overlap[n_gate_metrics] = {0,0,0,0,0,0,0};
            for (const StoredVariant& variant : all_variants) {
                const bool overlap = variant.atac_soft_overlap;
                gate_counts[0]++;
                if (overlap) gate_overlap[0]++;
                if (!variant.is_biallelic || variant.pattern == NULL) continue;
                gate_counts[1]++;
                if (overlap) gate_overlap[1]++;
                if (variant.nuclear_exclusion != NUCLEAR_EXCLUDE_NONE) continue;
                gate_counts[2]++;
                if (overlap) gate_overlap[2]++;
                if (variant.atac_cov_score < atac_min_cov) continue;
                gate_counts[3]++;
                if (overlap) gate_overlap[3]++;
                if (!variant.passes_qc) continue;
                gate_counts[4]++;
                if (overlap) gate_overlap[4]++;
                if (panel_position_is_excluded(
                        atac_excluded_positions, variant.chrom_idx,
                        variant.pos)) continue;
                gate_counts[5]++;
                if (overlap) gate_overlap[5]++;
                if (operational_roster_call_rate(
                        variant.pattern, atac_roster_sample_indices) <
                    pairwise_rescue_min_call_rate) continue;
                gate_counts[6]++;
                if (overlap) gate_overlap[6]++;
            }
            for (int gate = 0; gate < n_gate_metrics; ++gate) {
                candidate_audit << "gate\t" << gate_names[gate] << '\t'
                    << gate_counts[gate] << '\t'
                    << gate_overlap[gate] << '\n';
            }
            long eligible_overlap = 0;
            for (long value : entry_soft_counts) eligible_overlap += value;
            candidate_audit << "gate\tatac_eligible_scalar_pool_ranked\t"
                << atac_grouped << '\t' << eligible_overlap << '\n';

            const char* mac_names[6] = {
                "MAC_1", "MAC_2", "MAC_3_5", "MAC_6_10",
                "MAC_11_20", "MAC_21_PLUS"
            };
            long mac_candidate[6] = {0,0,0,0,0,0};
            long mac_candidate_overlap[6] = {0,0,0,0,0,0};
            long mac_selected[6] = {0,0,0,0,0,0};
            long mac_selected_overlap[6] = {0,0,0,0,0,0};
            const char* call_names[4] = {
                "CALL_RATE_0.50_0.75", "CALL_RATE_0.75_0.90",
                "CALL_RATE_0.90_LT1", "CALL_RATE_1.00"
            };
            long call_candidate[4] = {0,0,0,0};
            long call_candidate_overlap[4] = {0,0,0,0};
            long call_selected[4] = {0,0,0,0};
            long call_selected_overlap[4] = {0,0,0,0};
            for (size_t ai = 0; ai < atac_allocation.size(); ++ai) {
                const AtacAllocationEntry& entry = atac_allocation[ai];
                const VariantPattern* pattern = all_variants[
                    entry.representative_variant_idx].pattern;
                int called = 0;
                int alt_count = 0;
                for (int sample_index : atac_roster_sample_indices) {
                    const int dosage = pattern_dosage(pattern, sample_index);
                    if (dosage < 0) continue;
                    called++;
                    alt_count += dosage;
                }
                const int mac = min(alt_count, 2 * called - alt_count);
                int mac_bin = mac <= 1 ? 0 : mac == 2 ? 1 :
                    mac <= 5 ? 2 : mac <= 10 ? 3 : mac <= 20 ? 4 : 5;
                const float call_rate = atac_roster_sample_indices.empty()
                    ? 0.0f : (float)called /
                               atac_roster_sample_indices.size();
                int call_bin = call_rate < 0.75f ? 0 :
                    call_rate < 0.90f ? 1 : call_rate < 1.0f ? 2 : 3;
                mac_candidate[mac_bin] += entry.capacity;
                mac_candidate_overlap[mac_bin] += entry_soft_counts[ai];
                mac_selected[mac_bin] += entry_selected_counts[ai];
                mac_selected_overlap[mac_bin] +=
                    entry_selected_soft_counts[ai];
                call_candidate[call_bin] += entry.capacity;
                call_candidate_overlap[call_bin] += entry_soft_counts[ai];
                call_selected[call_bin] += entry_selected_counts[ai];
                call_selected_overlap[call_bin] +=
                    entry_selected_soft_counts[ai];
            }
            for (int bin = 0; bin < 6; ++bin) {
                candidate_audit << "mac\tcandidate_" << mac_names[bin]
                    << '\t' << mac_candidate[bin] << '\t'
                    << mac_candidate_overlap[bin] << '\n';
                candidate_audit << "mac\tselected_" << mac_names[bin]
                    << '\t' << mac_selected[bin] << '\t'
                    << mac_selected_overlap[bin] << '\n';
            }
            for (int bin = 0; bin < 4; ++bin) {
                candidate_audit << "call_rate\tcandidate_"
                    << call_names[bin] << '\t' << call_candidate[bin]
                    << '\t' << call_candidate_overlap[bin] << '\n';
                candidate_audit << "call_rate\tselected_"
                    << call_names[bin] << '\t' << call_selected[bin]
                    << '\t' << call_selected_overlap[bin] << '\n';
            }
            candidate_audit << "selection\tatac_demux_selected\t"
                << selected_atac_demux_indices.size() << '\t'
                << accumulate(entry_selected_soft_counts.begin(),
                              entry_selected_soft_counts.end(), 0L) << '\n';
            candidate_audit.close();

            const string pair_audit_path =
                atac_audit_prefix + ".pair_capacity.tsv";
            ofstream pair_audit(pair_audit_path.c_str());
            if (!pair_audit) {
                fprintf(stderr,
                    "ERROR: Could not write ATAC pair-capacity audit: %s\n",
                    pair_audit_path.c_str());
                exit(1);
            }
            pair_audit <<
                "scope\tdataset\tobjective_id\tsignature_a\tsignature_b\t"
                "library_ids\tcandidate_joint_callable\t"
                "candidate_distinguishing\tcandidate_rna_demux_overlap\t"
                "candidate_rna_disjoint\tcandidate_soft_distinguishing\t"
                "candidate_soft_rna_demux_overlap\t"
                "selected_distinguishing\tselected_soft_distinguishing\t"
                "requested_floor\teffective_floor\tstatus\n";
            // Donor-pair counts are diagnostics only and never affect ATAC
            // eligibility, clade allocation, ranking, or target filling.
            for (const PoolDonorPair& donor_pair : atac_pool_donor_pairs) {
                long joint = 0, capacity = 0, overlap = 0, selected = 0;
                long soft_capacity = 0, soft_overlap = 0,
                     soft_selected = 0;
                for (size_t ai = 0; ai < atac_allocation.size(); ++ai) {
                    const AtacAllocationEntry& entry = atac_allocation[ai];
                    const VariantPattern* pattern = all_variants[
                        entry.representative_variant_idx].pattern;
                    const int first = pattern_dosage(
                        pattern, donor_pair.first_index);
                    const int second = pattern_dosage(
                        pattern, donor_pair.second_index);
                    if (first < 0 || second < 0) continue;
                    joint += entry.capacity;
                    const float difference = fabsf((float)first / 2.0f -
                                                   (float)second / 2.0f);
                    if (difference >= min_pair_score) {
                        capacity += entry.capacity;
                        overlap += entry_soft_counts[ai];
                        selected += entry_selected_counts[ai];
                    } else if (difference >= 0.25f) {
                        soft_capacity += entry.capacity;
                        soft_overlap += entry_soft_counts[ai];
                        soft_selected += entry_selected_counts[ai];
                    }
                }
                string library_list;
                for (size_t li = 0; li < donor_pair.library_ids.size(); ++li) {
                    if (li) library_list += ",";
                    library_list += to_string(donor_pair.library_ids[li]);
                }
                string first_signature = donor_pair.first_name + ":1";
                string second_signature = donor_pair.second_name + ":1";
                if (second_signature < first_signature)
                    swap(first_signature, second_signature);
                pair_audit << "donor_pair\tatac_demux\t"
                    << first_signature << "||" << second_signature << '\t'
                    << first_signature << '\t' << second_signature << '\t'
                    << library_list << '\t' << joint << '\t' << capacity
                    << '\t' << overlap << '\t' << capacity - overlap << '\t'
                    << soft_capacity << '\t' << soft_overlap << '\t'
                    << selected << '\t' << soft_selected
                    << "\t0\t0\tDIAGNOSTIC_ONLY_NOT_A_SELECTION_OBJECTIVE\n";
            }
            pair_audit.close();
            fprintf(stderr, "  ATAC candidate/SFS audit: %s\n",
                candidate_audit_path.c_str());
            fprintf(stderr, "  ATAC pair-capacity audit: %s\n",
                pair_audit_path.c_str());
        }
        
        end_time = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
        fprintf(stderr, "  ATAC-demux: selected %zu / %d target variants in %ld seconds\n",
            selected_atac_demux_indices.size(), atac_num, duration.count());
    }

    // ========================================================================
    // COMPUTE DEMUX_SCORE FOR HET CANDIDATES AND WRITE HET OUTPUT
    // ========================================================================
    
    if (!het_outfile.empty() && !all_het_candidates.empty() && !selected_het_indices.empty()) {
        fprintf(stderr, "\n========================================\n");
        fprintf(stderr, "Computing DEMUX_SCORE for het candidates and writing output...\n");
        start_time = chrono::high_resolution_clock::now();
        
        // Compute DEMUX_SCORE for each het candidate from stored genotypes
        fprintf(stderr, "  Computing DEMUX_SCORE for %zu selected het candidates...\n", selected_het_indices.size());
        
        #pragma omp parallel for schedule(dynamic, 1000)
        for (size_t i = 0; i < selected_het_indices.size(); ++i) {
            HetCandidate& hc = all_het_candidates[selected_het_indices[i]];
            const StoredVariant& sv_ref = all_variants[hc.variant_idx];
            
            // Use precomputed dosage-aware bitsets from StoredVariant
            // Look up downsample probability
            double samp_prob = 1.0;
            int ta = sv_ref.total_alt;
            int ta_flip = sv_ref.total_alt_flip;
            
            if (sv_ref.pattern->present.count() == (size_t)hc.n_samples) {
                // No missing data
                if (ta < ta_flip || (ta == ta_flip && sv_ref.pattern->alt < sv_ref.pattern->alt_flip)) {
                    if (downsample_prob.count(sv_ref.pattern->alt) > 0) {
                        samp_prob = downsample_prob.at(sv_ref.pattern->alt);
                    }
                } else {
                    if (downsample_prob.count(sv_ref.pattern->alt_flip) > 0) {
                        samp_prob = downsample_prob.at(sv_ref.pattern->alt_flip);
                    }
                }
            } else {
                // Has missing data
                pair<bitset<NBITS>, bitset<NBITS>> key;
                if (ta < ta_flip || (ta == ta_flip && sv_ref.pattern->alt < sv_ref.pattern->alt_flip)) {
                    key = make_pair(sv_ref.pattern->present, sv_ref.pattern->alt);
                } else {
                    key = make_pair(sv_ref.pattern->present, sv_ref.pattern->alt_flip);
                }
                if (downsample_miss_prob.count(key) > 0) {
                    samp_prob = downsample_miss_prob.at(key);
                }
            }
            
            hc.demux_score = (samp_prob < 1.0) ? (float)samp_prob : 1.0f;
        }
        
        fprintf(stderr, "  Writing het output to %s...\n", het_outfile.c_str());
        
        // Create output header
        htsFile* hdr_reader = hts_open(vcf_file.c_str(), "r");
        bcf_hdr_t* het_header = bcf_hdr_read(hdr_reader);
        hts_close(hdr_reader);
        
        // Add INFO headers - ALL fields
        bcf_hdr_append(het_header, "##INFO=<ID=DEMUX_SCORE,Number=1,Type=Float,Description=\"Clade rarity score\">");
        bcf_hdr_append(het_header, "##INFO=<ID=HET_SCORE,Number=1,Type=Float,Description=\"het_freq * (1 - missing_freq * 0.5)\">");
        bcf_hdr_append(het_header, "##INFO=<ID=HET_FREQ,Number=1,Type=Float,Description=\"Fraction of called individuals that are het\">");
        bcf_hdr_append(het_header, "##INFO=<ID=MISSING_FREQ,Number=1,Type=Float,Description=\"Fraction of individuals with missing genotype\">");
        bcf_hdr_append(het_header, "##INFO=<ID=ANNOT_SCORE,Number=1,Type=Float,Description=\"Annotation weight (1.0=genic, 0.1=intergenic)\">");
        bcf_hdr_append(het_header, "##INFO=<ID=COV_SCORE,Number=1,Type=Float,Description=\"RNA coverage score\">");
        bcf_hdr_append(het_header, "##INFO=<ID=BIN_ID,Number=1,Type=String,Description=\"Chromosome:bin_index for evenness tracking\">");
        bcf_hdr_append(het_header, "##INFO=<ID=HET_SEL_SCORE,Number=1,Type=Float,Description=\"Selection ranking score: final_het_score * call_rate * (log2(cov+1) + 0.1*annot_boost)\">");
        bcf_hdr_append(het_header, "##INFO=<ID=HET_RAW_SCORE,Number=1,Type=Float,Description=\"Same as HET_SCORE; explicit name.\">");
        if (bcf_hdr_sync(het_header) < 0) {
            fprintf(stderr, "ERROR: Failed to sync het header\n");
            exit(1);
        }
        
        // Open output BCF file
        htsFile* het_out = hts_open(het_outfile.c_str(), "wb");
        if (het_out == NULL) {
            fprintf(stderr, "ERROR: Failed to open het output: %s\n", het_outfile.c_str());
            bcf_hdr_destroy(het_header);
            exit(1);
        }
        hts_set_threads(het_out, n_threads);
        if (bcf_hdr_write(het_out, het_header) < 0) {
            fprintf(stderr, "ERROR: Failed to write het header\n");
            exit(1);
        }
        
        // Sort selected candidates by chromosome then position for ordered output
        vector<HetCandidate*> selected_sorted;
        for (size_t idx : selected_het_indices) {
            selected_sorted.push_back(&all_het_candidates[idx]);
        }
        sort(selected_sorted.begin(), selected_sorted.end(), 
            [](const HetCandidate* a, const HetCandidate* b) {
                if (a->chrom_idx != b->chrom_idx) return a->chrom_idx < b->chrom_idx;
                return a->pos < b->pos;
            });
        
        // Write records
        bcf1_t* het_rec = bcf_init();
        for (const HetCandidate* hc : selected_sorted) {
            bcf_clear(het_rec);
            
            // Set chromosome and position
            het_rec->rid = bcf_hdr_name2id(het_header, hc->chrom.c_str());
            het_rec->pos = hc->pos;
            het_rec->qual = all_variants[hc->variant_idx].qual;
            
            // Set alleles
            string alleles = hc->ref_allele + "," + hc->alt_allele;
            bcf_update_alleles_str(het_header, het_rec, alleles.c_str());
            
            // Set genotypes (unpack from packed storage)
            vector<int32_t> het_gt;
            unpack_genotypes(all_variants[hc->variant_idx].genotypes_packed, hc->n_samples, het_gt);
            bcf_update_genotypes(het_header, het_rec, het_gt.data(), het_gt.size());
            
            // Set INFO fields - ALL scores
            float demux_score = hc->demux_score;
            // HET_SCORE uses the INFO formula: het_freq * (1 - missing_freq * 0.5)
            float het_score_info = hc->het_freq * (1.0f - hc->missing_freq * 0.5f);
            float het_freq = hc->het_freq;
            float missing_freq = hc->missing_freq;
            float annot = hc->annot_score;
            float cov = hc->cov_score;
            
            bcf_update_info_float(het_header, het_rec, "DEMUX_SCORE", &demux_score, 1);
            bcf_update_info_float(het_header, het_rec, "HET_SCORE", &het_score_info, 1);
            bcf_update_info_float(het_header, het_rec, "HET_FREQ", &het_freq, 1);
            bcf_update_info_float(het_header, het_rec, "MISSING_FREQ", &missing_freq, 1);
            bcf_update_info_float(het_header, het_rec, "ANNOT_SCORE", &annot, 1);
            bcf_update_info_float(het_header, het_rec, "COV_SCORE", &cov, 1);
            
            string bin_id = make_bin_id(hc->chrom, hc->bin_idx);
            bcf_update_info_string(het_header, het_rec, "BIN_ID", bin_id.c_str());
            
            if (enable_pool_scoring) {
                float sel_score = hc->het_sel_score;
                bcf_update_info_float(het_header, het_rec, "HET_SEL_SCORE", &sel_score, 1);
            } else {
                float sel_score = hc->het_sel_score;
                bcf_update_info_float(het_header, het_rec, "HET_SEL_SCORE", &sel_score, 1);
            }
            {
                float raw_score = hc->het_freq * (1.0f - hc->missing_freq * 0.5f);
                bcf_update_info_float(het_header, het_rec, "HET_RAW_SCORE", &raw_score, 1);
            }
            
            if (bcf_write(het_out, het_header, het_rec) < 0) {
                fprintf(stderr, "ERROR: Failed to write het record\n");
                exit(1);
            }
        }
        
        bcf_destroy(het_rec);
        hts_close(het_out);
        bcf_hdr_destroy(het_header);
        
        // Index het output BCF
        fprintf(stderr, "  Indexing het output BCF...\n");
        if (bcf_index_build(het_outfile.c_str(), 14) < 0) {
            fprintf(stderr, "WARNING: Failed to create index for %s\n", het_outfile.c_str());
        }
        
        // Retain het candidates until the optional ATAC-het writer has finished.
        // selected_atac_het_indices indexes this same vector; freeing it here caused
        // a use-after-free/segmentation fault whenever RNA-het and ATAC-het outputs
        // were requested in the same run.
        
        end_time = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
        fprintf(stderr, "Het output complete: %zu sites in %ld seconds\n", 
            selected_het_sites.size(), duration.count());
        fprintf(stderr, "Het output: %s\n", het_outfile.c_str());
    }

    // ========================================================================
    // WRITE ATAC-HET OUTPUT (if --atac_het_output set)
    // ========================================================================
    
    if (!atac_het_outfile.empty() && !selected_atac_het_indices.empty()) {
        fprintf(stderr, "\n========================================\n");
        fprintf(stderr, "Writing ATAC-het output to %s...\n", atac_het_outfile.c_str());
        start_time = chrono::high_resolution_clock::now();
        
        htsFile* ahdr_reader = hts_open(vcf_file.c_str(), "r");
        bcf_hdr_t* atac_het_header = bcf_hdr_read(ahdr_reader);
        hts_close(ahdr_reader);
        
        bcf_hdr_append(atac_het_header, "##INFO=<ID=HET_SCORE,Number=1,Type=Float,Description=\"het_freq * (1 - missing_freq * 0.5)\">");
        bcf_hdr_append(atac_het_header, "##INFO=<ID=HET_FREQ,Number=1,Type=Float,Description=\"Fraction of called operational pool donors that are heterozygous\">");
        bcf_hdr_append(atac_het_header, "##INFO=<ID=MISSING_FREQ,Number=1,Type=Float,Description=\"Fraction of operational pool donors without a complete genotype\">");
        bcf_hdr_append(atac_het_header, "##INFO=<ID=COV_SCORE,Number=1,Type=Float,Description=\"ATAC coverage score\">");
        bcf_hdr_append(atac_het_header, "##INFO=<ID=BIN_ID,Number=1,Type=String,Description=\"Chromosome:bin_index\">");
        bcf_hdr_append(atac_het_header, "##INFO=<ID=HET_SEL_SCORE,Number=1,Type=Float,Description=\"ATAC-het selection ranking score: final_het_score * call_rate * log2(atac_cov+1)\">");
        bcf_hdr_append(atac_het_header, "##INFO=<ID=HET_RAW_SCORE,Number=1,Type=Float,Description=\"Same as HET_SCORE; explicit name.\">");
        bcf_hdr_append(atac_het_header, "##INFO=<ID=ATAC_REGULATORY,Number=1,Type=Integer,Description=\"1 when the site overlaps the supplied regulatory BED or default strand-aware 2-kb upstream promoter prior; ranking annotation, not an eligibility filter\">");
        append_atac_soft_overlap_headers(
            atac_het_header, atac_soft_overlap_multiplier);
        {
            string regulatory_header =
                "##cellbouncer_atac_regulatory_source=" + atac_regulatory_source;
            bcf_hdr_append(atac_het_header, regulatory_header.c_str());
        }
        if (!atac_excluded_positions.empty()) {
            bcf_hdr_append(atac_het_header,
                "##cellbouncer_atac_rna_locus_exclusion=position_only");
        }
        if (bcf_hdr_sync(atac_het_header) < 0) {
            fprintf(stderr, "ERROR: Failed to sync ATAC-het header\n");
            exit(1);
        }
        
        htsFile* atac_het_out = hts_open(atac_het_outfile.c_str(), "wb");
        hts_set_threads(atac_het_out, n_threads);
        if (bcf_hdr_write(atac_het_out, atac_het_header) < 0) {
            fprintf(stderr, "ERROR: Failed to write ATAC-het header\n");
            exit(1);
        }
        
        // Sort selected ATAC-het candidates by position
        vector<HetCandidate*> atac_het_sorted;
        for (size_t idx : selected_atac_het_indices) {
            atac_het_sorted.push_back(&all_het_candidates[idx]);
        }
        sort(atac_het_sorted.begin(), atac_het_sorted.end(),
            [](const HetCandidate* a, const HetCandidate* b) {
                if (a->chrom_idx != b->chrom_idx) return a->chrom_idx < b->chrom_idx;
                return a->pos < b->pos;
            });
        
        bcf1_t* ah_rec = bcf_init();
        for (const HetCandidate* hc : atac_het_sorted) {
            bcf_clear(ah_rec);
            ah_rec->rid = bcf_hdr_name2id(atac_het_header, hc->chrom.c_str());
            ah_rec->pos = hc->pos;
            ah_rec->qual = all_variants[hc->variant_idx].qual;
            
            string alleles_str = hc->ref_allele + "," + hc->alt_allele;
            bcf_update_alleles_str(atac_het_header, ah_rec, alleles_str.c_str());
            const vector<uint8_t>& ah_gt_packed = all_variants[hc->variant_idx].genotypes_packed;
            vector<int32_t> ah_gt;
            unpack_genotypes(ah_gt_packed, num_samples, ah_gt);
            bcf_update_genotypes(atac_het_header, ah_rec, ah_gt.data(), ah_gt.size());
            
            const StoredVariant& ah_variant =
                all_variants[hc->variant_idx];
            int ah_roster_called = 0;
            int ah_roster_het = 0;
            for (int sample_index : atac_roster_sample_indices) {
                const int dosage = pattern_dosage(
                    ah_variant.pattern, sample_index);
                if (dosage < 0) continue;
                ah_roster_called++;
                if (dosage == 1) ah_roster_het++;
            }
            float hf = ah_roster_called > 0
                ? (float)ah_roster_het / ah_roster_called : 0.0f;
            float mf = atac_roster_sample_indices.empty() ? 1.0f :
                1.0f - (float)ah_roster_called /
                       atac_roster_sample_indices.size();
            float hs_info = hf * (1.0f - mf * 0.5f);
            float ac = hc->atac_cov_score;
            bcf_update_info_float(atac_het_header, ah_rec, "HET_SCORE", &hs_info, 1);
            bcf_update_info_float(atac_het_header, ah_rec, "HET_FREQ", &hf, 1);
            bcf_update_info_float(atac_het_header, ah_rec, "MISSING_FREQ", &mf, 1);
            bcf_update_info_float(atac_het_header, ah_rec, "COV_SCORE", &ac, 1);
            
            string bin_id = make_bin_id(hc->chrom, hc->bin_idx);
            bcf_update_info_string(atac_het_header, ah_rec, "BIN_ID", bin_id.c_str());
            
            // Fix 5: HET_SEL_SCORE and HET_RAW_SCORE parity with RNA-het output
            float ah_sel_score = hc->atac_het_sel_score;
            bcf_update_info_float(atac_het_header, ah_rec, "HET_SEL_SCORE", &ah_sel_score, 1);
            float ah_raw_score = hs_info;
            bcf_update_info_float(atac_het_header, ah_rec, "HET_RAW_SCORE", &ah_raw_score, 1);
            int32_t ah_regulatory = hc->atac_regulatory ? 1 : 0;
            bcf_update_info_int32(atac_het_header, ah_rec,
                "ATAC_REGULATORY", &ah_regulatory, 1);
            int32_t ah_soft_overlap = ah_variant.atac_soft_overlap ? 1 : 0;
            bcf_update_info_int32(atac_het_header, ah_rec,
                "ATAC_RNA_SOFT_OVERLAP", &ah_soft_overlap, 1);
            
            if (bcf_write(atac_het_out, atac_het_header, ah_rec) < 0) {
                fprintf(stderr, "ERROR: Failed to write ATAC-het record\n");
                exit(1);
            }
        }
        
        bcf_destroy(ah_rec);
        hts_close(atac_het_out);
        bcf_hdr_destroy(atac_het_header);
        
        fprintf(stderr, "  Indexing ATAC-het output BCF...\n");
        bcf_index_build(atac_het_outfile.c_str(), 14);
        
        end_time = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
        fprintf(stderr, "ATAC-het output complete: %zu sites in %ld seconds\n",
            selected_atac_het_sites.size(), duration.count());
    }

    // Het and ATAC-het writers are both complete; release their shared
    // candidate storage before the remaining large output passes.
    if (collect_het) {
        all_het_candidates.clear();
        all_het_candidates.shrink_to_fit();
        selected_het_indices.clear();
        selected_het_indices.shrink_to_fit();
        selected_atac_het_indices.clear();
        selected_atac_het_indices.shrink_to_fit();
    }

    // ========================================================================
    // WRITE SPECIES OUTPUT
    // ========================================================================
    
    if (collect_species && !selected_species_indices.empty()) {
        fprintf(stderr, "\n========================================\n");
        fprintf(stderr, "Writing species output to %s...\n", species_outfile.c_str());
        start_time = chrono::high_resolution_clock::now();
        
        // Create output header
        htsFile* sp_hdr_reader = hts_open(vcf_file.c_str(), "r");
        bcf_hdr_t* sp_header = bcf_hdr_read(sp_hdr_reader);
        hts_close(sp_hdr_reader);
        
        // Add INFO headers
        bcf_hdr_append(sp_header, "##INFO=<ID=HET_FREQ,Number=1,Type=Float,Description=\"Fraction of called individuals that are het\">");
        bcf_hdr_append(sp_header, "##INFO=<ID=MISSING_FREQ,Number=1,Type=Float,Description=\"Fraction of individuals with missing genotype\">");
        bcf_hdr_append(sp_header, "##INFO=<ID=COV_SCORE,Number=1,Type=Float,Description=\"Coverage score used for selection\">");
        bcf_hdr_append(sp_header, "##INFO=<ID=BIN_ID,Number=1,Type=String,Description=\"Chromosome:bin_index\">");
        bcf_hdr_append(sp_header, "##INFO=<ID=PAIR_ASSIGNED,Number=1,Type=String,Description=\"Species pair this site was assigned to during selection\">");
        bcf_hdr_append(sp_header, "##INFO=<ID=PAIR_CREDITED,Number=1,Type=String,Description=\"Comma-separated species pairs this site discriminates above species_min_pair_discrim\">");
        bcf_hdr_append(sp_header, "##INFO=<ID=SP_MAX_PAIR_SCORE,Number=1,Type=Float,Description=\"Maximum species pair score across all pairs\">");
        bcf_hdr_append(sp_header, "##INFO=<ID=SP_SELECTION_SCORE,Number=1,Type=Float,Description=\"Legacy aggregate score (max_pair_score * capped_coverage); recorded only, NOT used for ranking. Primary ranking uses pair_discrim.\">");
        bcf_hdr_append(sp_header, "##INFO=<ID=SP_MAX_PAIR_DISCRIM,Number=1,Type=Float,Description=\"Maximum pair_discrim across all pairs (primary ranking score)\">");
        bcf_hdr_append(sp_header, "##INFO=<ID=ASSIGNED_PAIR_DISCRIM,Number=1,Type=Float,Description=\"pair_discrim for the pair that pulled this site\">");
        bcf_hdr_append(sp_header, "##INFO=<ID=ASSIGNED_PAIR_SCORE,Number=1,Type=Float,Description=\"Legacy multiplicative pair_score for the assigned pair\">");
        bcf_hdr_append(sp_header, "##INFO=<ID=SPECIES_SELECTION_TIER,Number=1,Type=String,Description=\"strict when assigned pair_discrim clears the requested threshold; soft_backfill when it clears only half-threshold\">");
        if (species_coverage_mode == "atac") {
            bcf_hdr_append(sp_header, "##INFO=<ID=ATAC_REGULATORY,Number=1,Type=Integer,Description=\"1 when the site overlaps the supplied regulatory BED or default strand-aware 2-kb upstream promoter prior; ranking annotation, not an eligibility filter\">");
            append_atac_soft_overlap_headers(
                sp_header, atac_soft_overlap_multiplier);
            string regulatory_header =
                "##cellbouncer_atac_regulatory_source=" + atac_regulatory_source;
            bcf_hdr_append(sp_header, regulatory_header.c_str());
            if (!species_outgroup.empty() &&
                species_outgroup_max_fraction < 1.0f) {
                string outgroup_header =
                    "##cellbouncer_species_outgroup=" + species_outgroup;
                string cap_header =
                    "##cellbouncer_species_outgroup_only_max_fraction=" +
                    to_string(species_outgroup_max_fraction);
                bcf_hdr_append(sp_header, outgroup_header.c_str());
                bcf_hdr_append(sp_header, cap_header.c_str());
            }
        }
        if (species_coverage_mode == "atac" && !atac_excluded_positions.empty()) {
            bcf_hdr_append(sp_header,
                "##cellbouncer_atac_rna_locus_exclusion=position_only");
        }
        
        // Per-pair PAIR_DISCRIM fields
        for (int pi = 0; pi < panel_metadata.n_pairs; pi++) {
            string field_name = "PAIR_DISCRIM_" + panel_metadata.species_pairs[pi].first + 
                                "_" + panel_metadata.species_pairs[pi].second;
            string hdr_line = "##INFO=<ID=" + field_name + 
                ",Number=1,Type=Float,Description=\"Pair discrimination score for species pair (" +
                panel_metadata.species_pairs[pi].first + "," + 
                panel_metadata.species_pairs[pi].second + ")\">";
            bcf_hdr_append(sp_header, hdr_line.c_str());
        }

        // Per-pair PAIR_SCORE (legacy companion to PAIR_DISCRIM)
        for (int pi = 0; pi < panel_metadata.n_pairs; pi++) {
            string field_name = "PAIR_SCORE_" + panel_metadata.species_pairs[pi].first +
                                "_" + panel_metadata.species_pairs[pi].second;
            string hdr_line = "##INFO=<ID=" + field_name +
                ",Number=1,Type=Float,Description=\"Legacy multiplicative pair_score (" +
                panel_metadata.species_pairs[pi].first + "," +
                panel_metadata.species_pairs[pi].second + ")\">";
            bcf_hdr_append(sp_header, hdr_line.c_str());
        }

        // Per-species stats
        for (const string& sp : panel_metadata.species_list) {
            string mf = "##INFO=<ID=SP_MEAN_ALT_" + sp + ",Number=1,Type=Float,Description=\"Weighted mean alt dosage in species " + sp + "\">";
            string cf = "##INFO=<ID=SP_CALL_RATE_" + sp + ",Number=1,Type=Float,Description=\"Weighted call rate in species " + sp + "\">";
            string df = "##INFO=<ID=SP_MAJOR_DOSAGE_FRAC_" + sp + ",Number=1,Type=Float,Description=\"Weighted major dosage fraction in species " + sp + "\">";
            string nf = "##INFO=<ID=SP_N_CALLED_" + sp + ",Number=1,Type=Integer,Description=\"Unweighted count of called individuals in species " + sp + "\">";
            bcf_hdr_append(sp_header, mf.c_str());
            bcf_hdr_append(sp_header, cf.c_str());
            bcf_hdr_append(sp_header, df.c_str());
            bcf_hdr_append(sp_header, nf.c_str());
        }
        
        if (bcf_hdr_sync(sp_header) < 0) {
            fprintf(stderr, "ERROR: Failed to sync species header\n");
            exit(1);
        }
        
        // Open output
        htsFile* sp_out = hts_open(species_outfile.c_str(), "wb");
        hts_set_threads(sp_out, n_threads);
        if (bcf_hdr_write(sp_out, sp_header) < 0) {
            fprintf(stderr, "ERROR: Failed to write species header\n");
            exit(1);
        }
        
        // Sort by position
        vector<SpeciesCandidate*> sp_sorted;
        for (size_t idx : selected_species_indices) {
            sp_sorted.push_back(&all_species_candidates[idx]);
        }
        sort(sp_sorted.begin(), sp_sorted.end(),
            [](const SpeciesCandidate* a, const SpeciesCandidate* b) {
                if (a->chrom_idx != b->chrom_idx) return a->chrom_idx < b->chrom_idx;
                return a->pos < b->pos;
            });
        
        // Write records
        bcf1_t* sp_rec = bcf_init();
        for (const SpeciesCandidate* sc : sp_sorted) {
            bcf_clear(sp_rec);
            
            sp_rec->rid = bcf_hdr_name2id(sp_header, sc->chrom.c_str());
            sp_rec->pos = sc->pos;
            sp_rec->qual = all_variants[sc->variant_idx].qual;
            
            string alleles_str = sc->ref_allele + "," + sc->alt_allele;
            bcf_update_alleles_str(sp_header, sp_rec, alleles_str.c_str());
            const vector<uint8_t>& sp_gt_packed = all_variants[sc->variant_idx].genotypes_packed;
            vector<int32_t> sp_gt;
            unpack_genotypes(sp_gt_packed, num_samples, sp_gt);
            bcf_update_genotypes(sp_header, sp_rec, sp_gt.data(), sp_gt.size());
            
            // Compute het_freq and missing_freq for the record
            HetStats sp_hs = compute_het_stats_packed(sp_gt_packed, sc->n_samples);
            float hf = sp_hs.het_freq;
            float mf = sp_hs.missing_freq;
            float cv = sc->cov_score;
            
            bcf_update_info_float(sp_header, sp_rec, "HET_FREQ", &hf, 1);
            bcf_update_info_float(sp_header, sp_rec, "MISSING_FREQ", &mf, 1);
            bcf_update_info_float(sp_header, sp_rec, "COV_SCORE", &cv, 1);
            if (species_coverage_mode == "atac") {
                int32_t sp_regulatory = sc->atac_regulatory ? 1 : 0;
                bcf_update_info_int32(sp_header, sp_rec,
                    "ATAC_REGULATORY", &sp_regulatory, 1);
                int32_t sp_soft_overlap =
                    all_variants[sc->variant_idx].atac_soft_overlap ? 1 : 0;
                bcf_update_info_int32(sp_header, sp_rec,
                    "ATAC_RNA_SOFT_OVERLAP", &sp_soft_overlap, 1);
            }
            
            string bin_id = make_bin_id(sc->chrom, sc->bin_idx);
            bcf_update_info_string(sp_header, sp_rec, "BIN_ID", bin_id.c_str());
            
            // PAIR_ASSIGNED
            if (sc->assigned_pair >= 0) {
                string pair_str = panel_metadata.species_pairs[sc->assigned_pair].first + 
                                  "-" + panel_metadata.species_pairs[sc->assigned_pair].second;
                bcf_update_info_string(sp_header, sp_rec, "PAIR_ASSIGNED", pair_str.c_str());
                const string selection_tier =
                    sc->pair_discrim[sc->assigned_pair] >=
                        species_min_pair_discrim
                    ? "strict" : "soft_backfill";
                bcf_update_info_string(sp_header, sp_rec,
                    "SPECIES_SELECTION_TIER", selection_tier.c_str());
            }
            
            // PAIR_CREDITED: comma-separated list of all pairs this site discriminates
            // above species_min_pair_discrim. Used by the audit summary so multi-pair
            // contribution is not lost when only PAIR_ASSIGNED is inspected.
            {
                string credited_list;
                for (int pi = 0; pi < panel_metadata.n_pairs; pi++) {
                    if (sc->pair_discrim[pi] >= species_min_pair_discrim) {
                        if (!credited_list.empty()) credited_list += ",";
                        credited_list += panel_metadata.species_pairs[pi].first +
                                         "-" + panel_metadata.species_pairs[pi].second;
                    }
                }
                if (!credited_list.empty()) {
                    bcf_update_info_string(sp_header, sp_rec, "PAIR_CREDITED",
                        credited_list.c_str());
                }
            }
            
            // Per-pair PAIR_DISCRIM fields (legacy) and new PAIR_SCORE fields
            for (int pi = 0; pi < panel_metadata.n_pairs; pi++) {
                string field_name = "PAIR_DISCRIM_" + panel_metadata.species_pairs[pi].first +
                                    "_" + panel_metadata.species_pairs[pi].second;
                float pd = sc->pair_discrim[pi];
                bcf_update_info_float(sp_header, sp_rec, field_name.c_str(), &pd, 1);
            }
            
            // New aggregate score fields
            float mps = sc->max_pair_score;
            float ss = sc->selection_score;
            bcf_update_info_float(sp_header, sp_rec, "SP_MAX_PAIR_SCORE", &mps, 1);
            bcf_update_info_float(sp_header, sp_rec, "SP_SELECTION_SCORE", &ss, 1);
            
            // New audit fields
            float mpd = sc->max_pair_discrim;
            bcf_update_info_float(sp_header, sp_rec, "SP_MAX_PAIR_DISCRIM", &mpd, 1);

            if (sc->assigned_pair >= 0) {
                float apd = sc->pair_discrim[sc->assigned_pair];
                float aps = sc->pair_score[sc->assigned_pair];
                bcf_update_info_float(sp_header, sp_rec, "ASSIGNED_PAIR_DISCRIM", &apd, 1);
                bcf_update_info_float(sp_header, sp_rec, "ASSIGNED_PAIR_SCORE", &aps, 1);
            }

            // Per-pair legacy pair_score (companion to existing PAIR_DISCRIM_<A>_<B>)
            for (int pi = 0; pi < panel_metadata.n_pairs; pi++) {
                string fld = "PAIR_SCORE_" + panel_metadata.species_pairs[pi].first +
                             "_" + panel_metadata.species_pairs[pi].second;
                float v = sc->pair_score[pi];
                bcf_update_info_float(sp_header, sp_rec, fld.c_str(), &v, 1);
            }

            // Per-species stats
            for (int si = 0; si < (int)panel_metadata.species_list.size(); si++) {
                const string& sp = panel_metadata.species_list[si];
                float ma = sc->sp_mean_alt[si];
                float cr = sc->sp_call_rate[si];
                float md = sc->sp_major_dosage_frac[si];
                int   nc = sc->sp_n_called[si];
                string fm = "SP_MEAN_ALT_" + sp;
                string fc = "SP_CALL_RATE_" + sp;
                string fd = "SP_MAJOR_DOSAGE_FRAC_" + sp;
                string fn = "SP_N_CALLED_" + sp;
                bcf_update_info_float(sp_header, sp_rec, fm.c_str(), &ma, 1);
                bcf_update_info_float(sp_header, sp_rec, fc.c_str(), &cr, 1);
                bcf_update_info_float(sp_header, sp_rec, fd.c_str(), &md, 1);
                bcf_update_info_int32(sp_header, sp_rec, fn.c_str(), &nc, 1);
            }
            
            if (bcf_write(sp_out, sp_header, sp_rec) < 0) {
                fprintf(stderr, "ERROR: Failed to write species record\n");
                exit(1);
            }
        }
        
        bcf_destroy(sp_rec);
        hts_close(sp_out);
        bcf_hdr_destroy(sp_header);
        
        // Index species output
        fprintf(stderr, "  Indexing species output BCF...\n");
        if (bcf_index_build(species_outfile.c_str(), 14) < 0) {
            fprintf(stderr, "WARNING: Failed to create index for %s\n", species_outfile.c_str());
        }
        
        // Free species candidates
        all_species_candidates.clear();
        all_species_candidates.shrink_to_fit();
        selected_species_indices.clear();
        
        end_time = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
        fprintf(stderr, "Species output complete: %zu sites in %ld seconds\n",
            selected_species_sites.size(), duration.count());
        fprintf(stderr, "Species output: %s\n", species_outfile.c_str());
    }

    // ========================================================================
    // PASS 2: Write selected variants (Option C: Deterministic from stratified selection)
    // Also writes annotate_all output if requested
    // ========================================================================

    long final_count = 0;
    bcf_hdr_t* out_header = NULL;
    if (wants_generic_demux || !annotate_all_file.empty()) {
    
    fprintf(stderr, "\n========================================\n");
    fprintf(stderr, "Pass 2: Writing selected variants (deterministic, parallel, %d threads)...\n", n_threads);
    if (!annotate_all_file.empty()) {
        fprintf(stderr, "  Also writing full annotated output to: %s\n", annotate_all_file.c_str());
    }
    fflush(stderr);
    start_time = chrono::high_resolution_clock::now();
    
    // vcf_idx is reused from Pass 1
    
    // Get header for INFO field additions
    htsFile* hdr_reader = hts_open(vcf_file.c_str(), "r");
    out_header = bcf_hdr_read(hdr_reader);
    hts_close(hdr_reader);
    
    // Add ALL INFO headers (used for both demux and annotate_all outputs)
    bcf_hdr_append(out_header, "##INFO=<ID=DEMUX_SCORE,Number=1,Type=Float,Description=\"Clade rarity score\">");
    bcf_hdr_append(out_header, "##INFO=<ID=HET_FREQ,Number=1,Type=Float,Description=\"Fraction of called individuals that are heterozygous\">");
    bcf_hdr_append(out_header, "##INFO=<ID=MISSING_FREQ,Number=1,Type=Float,Description=\"Fraction of individuals with missing genotype\">");
    bcf_hdr_append(out_header, "##INFO=<ID=HET_SCORE,Number=1,Type=Float,Description=\"Combined het score: het_freq * (1 - missing_freq * 0.5)\">");
    bcf_hdr_append(out_header, "##INFO=<ID=ANNOT_SCORE,Number=1,Type=Float,Description=\"Annotation weight (1.0=genic, 0.1=intergenic)\">");
    bcf_hdr_append(out_header, "##INFO=<ID=COV_SCORE,Number=1,Type=Float,Description=\"Coverage score from bedgraph\">");
    bcf_hdr_append(out_header, "##INFO=<ID=BIN_ID,Number=1,Type=String,Description=\"Chromosome:bin_index for evenness tracking\">");
    bcf_hdr_append(out_header, "##INFO=<ID=SELECTION_SCORE,Number=1,Type=Float,Description=\"Within-clade selection score: log2(cov+1) + 0.1*annot_boost\">");
    bcf_hdr_append(out_header, "##INFO=<ID=NUCLEAR_PANEL_STATUS,Number=1,Type=String,Description=\"PASS or exclusion reason(s): MITO, NUMT, BLACKLIST\">");
    if (bcf_hdr_sync(out_header) < 0) {
        fprintf(stderr, "ERROR: Failed to sync header\n");
        exit(1);
    }
    
    // Create temp directory for per-chromosome demux outputs
    string temp_dir = outfile + ".tmp";
    mkdir(temp_dir.c_str(), 0755);
    fprintf(stderr, "  Demux temp directory: %s\n", temp_dir.c_str());
    
    // Create temp directory for per-chromosome annotate_all outputs (if needed)
    string annot_temp_dir;
    if (!annotate_all_file.empty()) {
        annot_temp_dir = annotate_all_file + ".tmp";
        mkdir(annot_temp_dir.c_str(), 0755);
        fprintf(stderr, "  Annotate temp directory: %s\n", annot_temp_dir.c_str());
    }
    fflush(stderr);
    
    atomic<long> total_kept(0);
    atomic<long> total_elim(0);
    atomic<long> total_annotated(0);  // For annotate_all tracking
    atomic<int> chroms_written(0);
    
    // Vector to track which chromosomes have output (for concatenation order)
    vector<string> chrom_files(chroms.size());
    vector<string> annot_chrom_files(chroms.size());  // For annotate_all
    vector<uint8_t> chrom_has_data(chroms.size(), 0);
    vector<uint8_t> annot_chrom_has_data(chroms.size(), 0);
    
    #pragma omp parallel num_threads(n_threads)
    {
        // Thread-local copy of output header
        bcf_hdr_t* thread_out_header = bcf_hdr_dup(out_header);
        
        // Thread-local output record
        bcf1_t* out_rec = bcf_init();
        
        #pragma omp for schedule(dynamic, 1)
        for (size_t c = 0; c < chroms.size(); c++) {
            const string& chrom = chroms[c];
            
            // Get range of variants for this chromosome from memory
            size_t range_start = chrom_ranges[c].first;
            size_t range_end = chrom_ranges[c].second;
            
            if (range_start >= range_end) {
                chroms_written++;
                continue;
            }
            
            // Create temp output file for demux (filtered) output
            string temp_file = temp_dir + "/chr_" + to_string(c) + ".bcf";
            chrom_files[c] = temp_file;
            htsFile* temp_out = hts_open(temp_file.c_str(), "wb");
            (void)bcf_hdr_write(temp_out, thread_out_header);
            
            // Create temp output file for annotate_all (all sites) if requested
            htsFile* annot_out = NULL;
            if (!annotate_all_file.empty()) {
                string annot_file = annot_temp_dir + "/chr_" + to_string(c) + ".bcf";
                annot_chrom_files[c] = annot_file;
                annot_out = hts_open(annot_file.c_str(), "wb");
                (void)bcf_hdr_write(annot_out, thread_out_header);
            }
            
            long local_kept = 0;
            long local_elim = 0;
            long local_annotated = 0;
            
            // Iterate over stored variants for this chromosome
            for (size_t vi = range_start; vi < range_end; vi++) {
                const StoredVariant& sv = all_variants[vi];
                
                // For annotate_all, we process ALL records (including multi-allelic)
                // For demux, we only process biallelic
                bool is_biallelic = sv.is_biallelic;
                
                // Compute DEMUX_SCORE (for reporting) and filtering decision
                // Option C: Use pre-computed selected_variant_indices set (deterministic)
                double samp_prob = 1.0;
                bool passes_filter = false;
                
                if (sv.passes_qc) {
                    // Use precomputed bitsets (same as what proc_bcf_record would return)
                    int p2_ta = sv.total_alt;
                    int p2_ta_flip = sv.total_alt_flip;
                    
                    if (sv.pattern->present.count() == (size_t)num_samples) {
                        if (p2_ta < p2_ta_flip || (p2_ta == p2_ta_flip && sv.pattern->alt < sv.pattern->alt_flip)) {
                            if (downsample_prob.count(sv.pattern->alt) > 0) {
                                samp_prob = downsample_prob[sv.pattern->alt];
                            }
                        } else {
                            if (downsample_prob.count(sv.pattern->alt_flip) > 0) {
                                samp_prob = downsample_prob[sv.pattern->alt_flip];
                            }
                        }
                    } else {
                        pair<bitset<NBITS>, bitset<NBITS>> key;
                        if (p2_ta < p2_ta_flip || (p2_ta == p2_ta_flip && sv.pattern->alt < sv.pattern->alt_flip)) {
                            key = make_pair(sv.pattern->present, sv.pattern->alt);
                        } else {
                            key = make_pair(sv.pattern->present, sv.pattern->alt_flip);
                        }
                        if (downsample_miss_prob.count(key) > 0) {
                            samp_prob = downsample_miss_prob[key];
                        }
                    }
                    
                    // Option C: Deterministic selection using pre-computed set
                    // (replaces random sampling)
                    passes_filter = selected_variant_mask[vi] != 0;
                }
                
                // Compute ALL scores for this record
                float demux_score = (samp_prob < 1.0) ? (float)samp_prob : 1.0f;
                const float annot_score = sv.annot_score;
                const float cov_score = sv.rna_cov_score;
                float selection_score = compute_selection_score(cov_score, annot_score);
                int bin_idx = get_bin_idx(sv.pos, bin_size);
                string bin_id = make_bin_id(chrom, bin_idx);
                
                // Compute HET scores from packed genotypes using helper (C1)
                float het_freq = 0.0f;
                float missing_freq = 0.0f;
                float het_score = 0.0f;
                
                if (!sv.genotypes_packed.empty()) {
                    HetStats hs_p2 = compute_het_stats_packed(sv.genotypes_packed, num_samples);
                    het_freq = hs_p2.het_freq;
                    missing_freq = hs_p2.missing_freq;
                    float missing_penalty = 1.0f - (missing_freq * 0.5f);
                    het_score = het_freq * missing_penalty;
                }
                
                // Build output record from stored data
                bcf_clear(out_rec);
                out_rec->rid = sv.chrom_idx;
                out_rec->pos = sv.pos;
                out_rec->qual = sv.qual;
                
                // Set alleles
                if (!sv.alleles.empty()) {
                    bcf_update_alleles_str(thread_out_header, out_rec, sv.alleles.c_str());
                }
                
                // Set genotypes: use raw for non-biallelic (lossless), packed for biallelic
                if (!sv.genotypes_raw.empty()) {
                    // Non-biallelic: use preserved raw genotypes (avoids multi-allelic collapse)
                    bcf_update_genotypes(thread_out_header, out_rec,
                        sv.genotypes_raw.data(), sv.genotypes_raw.size());
                } else if (!sv.genotypes_packed.empty()) {
                    vector<int32_t> unpacked_gt;
                    unpack_genotypes(sv.genotypes_packed, num_samples, unpacked_gt);
                    bcf_update_genotypes(thread_out_header, out_rec, unpacked_gt.data(), unpacked_gt.size());
                }
                
                // Add ALL INFO fields
                bcf_update_info_float(thread_out_header, out_rec, "DEMUX_SCORE", &demux_score, 1);
                bcf_update_info_float(thread_out_header, out_rec, "HET_FREQ", &het_freq, 1);
                bcf_update_info_float(thread_out_header, out_rec, "MISSING_FREQ", &missing_freq, 1);
                bcf_update_info_float(thread_out_header, out_rec, "HET_SCORE", &het_score, 1);
                bcf_update_info_float(thread_out_header, out_rec, "ANNOT_SCORE", &annot_score, 1);
                bcf_update_info_float(thread_out_header, out_rec, "COV_SCORE", &cov_score, 1);
                bcf_update_info_float(thread_out_header, out_rec, "SELECTION_SCORE", &selection_score, 1);
                bcf_update_info_string(thread_out_header, out_rec, "BIN_ID", bin_id.c_str());
                string exclusion_reason = nuclear_exclusion_reason(sv.nuclear_exclusion);
                bcf_update_info_string(thread_out_header, out_rec, "NUCLEAR_PANEL_STATUS", exclusion_reason.c_str());
                
                // Write to annotate_all output (ALL records)
                if (annot_out) {
                    bcf_write1(annot_out, thread_out_header, out_rec);
                    local_annotated++;
                }
                
                // Write to demux output (filtered biallelic records only)
                if (is_biallelic && passes_filter) {
                    bcf_write1(temp_out, thread_out_header, out_rec);
                    local_kept++;
                } else if (is_biallelic) {
                    local_elim++;
                }
            }
            
            hts_close(temp_out);
            if (annot_out) {
                hts_close(annot_out);
            }
            
            if (local_kept > 0) {
                chrom_has_data[c] = 1;
            } else {
                // Remove empty file
                remove(temp_file.c_str());
            }
            
            if (local_annotated > 0) {
                annot_chrom_has_data[c] = 1;
            } else if (annot_out) {
                remove(annot_chrom_files[c].c_str());
            }
            
            total_kept += local_kept;
            total_elim += local_elim;
            total_annotated += local_annotated;
            int done = ++chroms_written;
            
            if (local_kept > 0 || done % 100 == 0) {
                #pragma omp critical(stderr_log)
                {
                fprintf(stderr, "  Chromosome %s: kept %ld [%d/%zu]\n", 
                    chrom.c_str(), local_kept, done, chroms.size());
                fflush(stderr);
                }
            }
        }
        
        bcf_destroy(out_rec);
        bcf_hdr_destroy(thread_out_header);
    }
    
    // vcf_idx already freed after Pass 1 (C4)
    
    end_time = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    
    fprintf(stderr, "\nPass 2 complete: %ld demux kept, %ld eliminated", 
        total_kept.load(), total_elim.load());
    if (!annotate_all_file.empty()) {
        fprintf(stderr, ", %ld annotated", total_annotated.load());
    }
    fprintf(stderr, " in %ld seconds\n", duration.count());
    
    // ========================================================================
    // Concatenate chromosome files in order
    // ========================================================================
    
    fprintf(stderr, "Concatenating chromosome files...\n");
    start_time = chrono::high_resolution_clock::now();
    
    htsFile* outf = hts_open(outfile.c_str(), "wb");
    hts_set_threads(outf, n_threads);
    if (bcf_hdr_write(outf, out_header) < 0) {
        fprintf(stderr, "ERROR: Failed to write output header\n");
        exit(1);
    }
    
    bcf1_t* rec = bcf_init();
    for (size_t c = 0; c < chroms.size(); c++) {
        if (!chrom_has_data[c]) continue;
        
        htsFile* temp_in = hts_open(chrom_files[c].c_str(), "r");
        bcf_hdr_t* temp_hdr = bcf_hdr_read(temp_in);
        
        while (bcf_read(temp_in, temp_hdr, rec) >= 0) {
            bcf_write1(outf, out_header, rec);
            final_count++;
        }
        
        bcf_hdr_destroy(temp_hdr);
        hts_close(temp_in);
        
        // Remove temp file
        remove(chrom_files[c].c_str());
    }
    
    bcf_destroy(rec);
    hts_close(outf);
    
    // Index demux output BCF
    fprintf(stderr, "Indexing demux output BCF...\n");
    if (bcf_index_build(outfile.c_str(), 14) < 0) {
        fprintf(stderr, "WARNING: Failed to create index for %s\n", outfile.c_str());
    }
    
    // Remove demux temp directory
    rmdir(temp_dir.c_str());
    
    fprintf(stderr, "Demux output: %ld SNPs written to %s\n", final_count, outfile.c_str());
    
    // ========================================================================
    // Concatenate annotate_all chromosome files (if requested)
    // ========================================================================
    
    if (!annotate_all_file.empty()) {
        fprintf(stderr, "\nConcatenating annotate_all chromosome files...\n");
        
        htsFile* annot_outf = hts_open(annotate_all_file.c_str(), "wb");
        hts_set_threads(annot_outf, n_threads);
        if (bcf_hdr_write(annot_outf, out_header) < 0) {
            fprintf(stderr, "ERROR: Failed to write annotate_all header\n");
            exit(1);
        }
        
        bcf1_t* annot_rec = bcf_init();
        long annot_final_count = 0;
        
        for (size_t c = 0; c < chroms.size(); c++) {
            if (!annot_chrom_has_data[c]) continue;
            
            htsFile* temp_in = hts_open(annot_chrom_files[c].c_str(), "r");
            bcf_hdr_t* temp_hdr = bcf_hdr_read(temp_in);
            
            while (bcf_read(temp_in, temp_hdr, annot_rec) >= 0) {
                bcf_write1(annot_outf, out_header, annot_rec);
                annot_final_count++;
            }
            
            bcf_hdr_destroy(temp_hdr);
            hts_close(temp_in);
            
            // Remove temp file
            remove(annot_chrom_files[c].c_str());
        }
        
        bcf_destroy(annot_rec);
        hts_close(annot_outf);
        
        // Index annotate_all output BCF
        fprintf(stderr, "Indexing annotate_all output BCF...\n");
        if (bcf_index_build(annotate_all_file.c_str(), 14) < 0) {
            fprintf(stderr, "WARNING: Failed to create index for %s\n", annotate_all_file.c_str());
        }
        
        // Remove annotate_all temp directory
        rmdir(annot_temp_dir.c_str());
        
        fprintf(stderr, "Annotate_all output: %ld variants written to %s\n", annot_final_count, annotate_all_file.c_str());
    }

    } // end generic demux/annotate output pass
    
    // ========================================================================
    // ATAC-DEMUX OUTPUT (if --atac_output set)
    // Write selected ATAC-demux variants directly from in-memory data
    // ========================================================================
    
    if (!atac_outfile.empty() && !selected_atac_demux_indices.empty()) {
        fprintf(stderr, "\nWriting ATAC-demux output...\n");
        
        htsFile* atac_hdr_reader = hts_open(vcf_file.c_str(), "r");
        bcf_hdr_t* atac_out_header = bcf_hdr_read(atac_hdr_reader);
        hts_close(atac_hdr_reader);
        
        bcf_hdr_append(atac_out_header, "##INFO=<ID=COV_SCORE,Number=1,Type=Float,Description=\"ATAC coverage score\">");
        bcf_hdr_append(atac_out_header, "##INFO=<ID=BIN_ID,Number=1,Type=String,Description=\"Chromosome:bin_index\">");
        bcf_hdr_append(atac_out_header, "##INFO=<ID=ATAC_REGULATORY,Number=1,Type=Integer,Description=\"1 when the site overlaps the supplied regulatory BED or default strand-aware 2-kb upstream promoter prior; ranking annotation, not an eligibility filter\">");
        append_atac_soft_overlap_headers(
            atac_out_header, atac_soft_overlap_multiplier);
        {
            string regulatory_header =
                "##cellbouncer_atac_regulatory_source=" + atac_regulatory_source;
            bcf_hdr_append(atac_out_header, regulatory_header.c_str());
        }
        if (!atac_excluded_positions.empty()) {
            bcf_hdr_append(atac_out_header,
                "##cellbouncer_atac_rna_locus_exclusion=position_only");
        }
        if (bcf_hdr_sync(atac_out_header) < 0) {
            fprintf(stderr, "ERROR: Failed to sync ATAC-demux header\n");
            exit(1);
        }
        
        htsFile* atac_outf = hts_open(atac_outfile.c_str(), "wb");
        hts_set_threads(atac_outf, n_threads);
        if (bcf_hdr_write(atac_outf, atac_out_header) < 0) {
            fprintf(stderr, "ERROR: Failed to write ATAC-demux header\n");
            exit(1);
        }
        
        bcf1_t* atac_rec = bcf_init();
        long atac_count = 0;
        
        // Iterate in chromosome order
        for (size_t c = 0; c < chroms.size(); c++) {
            size_t range_start = chrom_ranges[c].first;
            size_t range_end = chrom_ranges[c].second;
            
            for (size_t vi = range_start; vi < range_end; vi++) {
                if (!selected_atac_demux_mask[vi]) continue;
                const StoredVariant& sv = all_variants[vi];
                
                bcf_clear(atac_rec);
                atac_rec->rid = sv.chrom_idx;
                atac_rec->pos = sv.pos;
                atac_rec->qual = sv.qual;
                
                if (!sv.alleles.empty()) {
                    bcf_update_alleles_str(atac_out_header, atac_rec, sv.alleles.c_str());
                }
                if (!sv.genotypes_packed.empty()) {
                    vector<int32_t> atac_unpacked_gt;
                    unpack_genotypes(sv.genotypes_packed, num_samples, atac_unpacked_gt);
                    bcf_update_genotypes(atac_out_header, atac_rec, atac_unpacked_gt.data(), atac_unpacked_gt.size());
                }
                
                float atac_cov = sv.atac_cov_score;
                bcf_update_info_float(atac_out_header, atac_rec, "COV_SCORE", &atac_cov, 1);
                
                int bin_idx = get_bin_idx(sv.pos, bin_size);
                string bid = make_bin_id(chroms[c], bin_idx);
                bcf_update_info_string(atac_out_header, atac_rec, "BIN_ID", bid.c_str());
                int32_t atac_regulatory = sv.atac_regulatory ? 1 : 0;
                bcf_update_info_int32(atac_out_header, atac_rec,
                    "ATAC_REGULATORY", &atac_regulatory, 1);
                int32_t atac_soft_overlap = sv.atac_soft_overlap ? 1 : 0;
                bcf_update_info_int32(atac_out_header, atac_rec,
                    "ATAC_RNA_SOFT_OVERLAP", &atac_soft_overlap, 1);
                
                if (bcf_write(atac_outf, atac_out_header, atac_rec) < 0) {
                    fprintf(stderr, "ERROR: Failed to write ATAC-demux record\n");
                    exit(1);
                }
                atac_count++;
            }
        }
        
        bcf_destroy(atac_rec);
        hts_close(atac_outf);
        bcf_hdr_destroy(atac_out_header);
        
        fprintf(stderr, "Indexing ATAC-demux output BCF...\n");
        bcf_index_build(atac_outfile.c_str(), 14);
        fprintf(stderr, "ATAC-demux output: %ld SNPs written to %s\n", atac_count, atac_outfile.c_str());
    }

    // NUMT-overlapping sites remain excluded from every nuclear panel.  Keep
    // the covered, callable sites in a separately labelled diagnostic BCF.
    long atac_numt_count = 0;
    if (!atac_numt_outfile.empty() && !atac_coverage_map.empty()) {
        fprintf(stderr, "\nWriting separate ATAC NUMT diagnostic output...\n");
        htsFile* numt_hdr_reader = hts_open(vcf_file.c_str(), "r");
        bcf_hdr_t* numt_header = bcf_hdr_read(numt_hdr_reader);
        hts_close(numt_hdr_reader);
        bcf_hdr_append(numt_header,
            "##INFO=<ID=COV_SCORE,Number=1,Type=Float,Description=\"ATAC coverage score\">");
        bcf_hdr_append(numt_header,
            "##INFO=<ID=NUMT_PANEL,Number=1,Type=Integer,Description=\"1 for an explicitly NUMT-overlapping diagnostic site excluded from nuclear panels\">");
        bcf_hdr_append(numt_header,
            "##INFO=<ID=ATAC_ROSTER_CALL_RATE,Number=1,Type=Float,Description=\"Complete genotype call rate among operational pool donors\">");
        bcf_hdr_append(numt_header,
            "##cellbouncer_panel_role=atac_numt_diagnostic");
        if (bcf_hdr_sync(numt_header) < 0) {
            fprintf(stderr, "ERROR: Failed to sync ATAC NUMT header\n");
            exit(1);
        }
        htsFile* numt_out = hts_open(atac_numt_outfile.c_str(), "wb");
        hts_set_threads(numt_out, n_threads);
        if (bcf_hdr_write(numt_out, numt_header) < 0) {
            fprintf(stderr, "ERROR: Failed to write ATAC NUMT header\n");
            exit(1);
        }
        bcf1_t* numt_rec = bcf_init();
        for (size_t vi = 0; vi < all_variants.size(); vi++) {
            const StoredVariant& sv = all_variants[vi];
            if (!(sv.nuclear_exclusion & NUCLEAR_EXCLUDE_NUMT)) continue;
            if (sv.nuclear_exclusion & NUCLEAR_EXCLUDE_MITO) continue;
            if (!sv.is_biallelic || !sv.passes_qc) continue;
            const float roster_call_rate = operational_roster_call_rate(
                sv.pattern, atac_roster_sample_indices);
            if (roster_call_rate < pairwise_rescue_min_call_rate) continue;
            float atac_cov = sv.atac_cov_score;
            if (atac_cov < atac_min_cov) continue;

            bcf_clear(numt_rec);
            numt_rec->rid = sv.chrom_idx;
            numt_rec->pos = sv.pos;
            numt_rec->qual = sv.qual;
            bcf_update_alleles_str(numt_header, numt_rec, sv.alleles.c_str());
            vector<int32_t> numt_gt;
            unpack_genotypes(sv.genotypes_packed, num_samples, numt_gt);
            bcf_update_genotypes(numt_header, numt_rec,
                numt_gt.data(), numt_gt.size());
            bcf_update_info_float(numt_header, numt_rec,
                "COV_SCORE", &atac_cov, 1);
            bcf_update_info_float(numt_header, numt_rec,
                "ATAC_ROSTER_CALL_RATE", &roster_call_rate, 1);
            int32_t one = 1;
            bcf_update_info_int32(numt_header, numt_rec,
                "NUMT_PANEL", &one, 1);
            if (bcf_write(numt_out, numt_header, numt_rec) < 0) {
                fprintf(stderr, "ERROR: Failed to write ATAC NUMT record\n");
                exit(1);
            }
            atac_numt_count++;
        }
        bcf_destroy(numt_rec);
        hts_close(numt_out);
        bcf_hdr_destroy(numt_header);
        if (bcf_index_build(atac_numt_outfile.c_str(), 14) < 0) {
            fprintf(stderr, "WARNING: Failed to index %s\n",
                atac_numt_outfile.c_str());
        }
        fprintf(stderr, "ATAC NUMT diagnostic output: %ld SNPs written to %s\n",
            atac_numt_count, atac_numt_outfile.c_str());
    }
    
    if (out_header) bcf_hdr_destroy(out_header);
    
    end_time = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    
    fprintf(stderr, "\n========================================\n");
    fprintf(stderr, "Pipeline complete in %ld seconds\n", duration.count());
    fprintf(stderr, "========================================\n");
    if (!outfile.empty()) {
        fprintf(stderr, "Demux output: %s (%ld SNPs)\n", outfile.c_str(), final_count);
    }
    if (!annotate_all_file.empty()) {
        fprintf(stderr, "Annotated output: %s (all variants)\n", annotate_all_file.c_str());
    }
    if (!het_outfile.empty()) {
        fprintf(stderr, "Het output: %s (%zu sites)\n", het_outfile.c_str(), selected_het_sites.size());
    }
    if (!species_outfile.empty()) {
        fprintf(stderr, "Species output: %s (%zu sites)\n", species_outfile.c_str(), selected_species_sites.size());
    }
    if (!atac_outfile.empty()) {
        fprintf(stderr, "ATAC-demux output: %s (%zu sites)\n", atac_outfile.c_str(), selected_atac_demux_indices.size());
    }
    if (!atac_het_outfile.empty()) {
        fprintf(stderr, "ATAC-het output: %s (%zu sites)\n", atac_het_outfile.c_str(), selected_atac_het_sites.size());
    }
    if (!atac_candidates_outfile.empty()) {
        fprintf(stderr,
            "ATAC candidate pool: %s (%ld union sites; %ld strict demux-eligible)\n",
            atac_candidates_outfile.c_str(), atac_candidates_written,
            atac_demux_eligible_written);
    }
    if (!atac_numt_outfile.empty()) {
        fprintf(stderr, "ATAC NUMT diagnostic output: %s (%ld sites)\n",
            atac_numt_outfile.c_str(), atac_numt_count);
    }

    return 0;
}
