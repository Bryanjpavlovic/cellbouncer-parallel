// demux_parallel.cpp
//
// Provenance: CellBouncer / Tet2025. Revision history at the bottom of this
// header block. Stable pipeline source name; version is tracked in-file.
//
// Revision history
//   V2.13.3 Targeted reliability repair: skip BAM targets with no active panel
//          SNPs, fall back to synchronous I/O when optional HTSlib helper-thread
//          setup fails, and preserve staged-output cleanup even when a legacy
//          deep helper exits instead of returning through C++ stack unwinding.
//   V2.13.2 BAM-index scheduling-stat compatibility: valid BAI/CSI indexes may
//          lack mapped/unmapped metadata for empty header targets. Both single-
//          and dual-panel counting now use a zero read-count estimate only for
//          work-unit sizing/ordering instead of aborting; target iterators and
//          allele counting remain unchanged.
//   V2.13.1 Compatibility rollback: demux VCF loading no longer reads,
//          validates, or filters on FORMAT/GQ. Historical Float GQ, canonical
//          Integer GQ, absent GQ, and other GQ encodings are all ignored. A
//          record with unavailable or unsupported GT structure is skipped rather
//          than aborting the complete panel load.
//   V2.13  --dump_pileup dual-panel support. The pileup sidecars
//          (.pileup_sites.tsv.gz / .pileup_obs.tsv.gz, consumed by the ambient
//          benchmark variant-consistency producer) were previously emitted only
//          by the single-panel count_alleles_parallel path; the production
//          launch uses --species_counts_output, which routes counting through
//          count_alleles_parallel_dual, so the flag was silently ignored.
//          count_alleles_parallel_dual now accepts dump_pileup/pileup_prefix
//          (appended trailing parameters, defaulted off, so every existing
//          caller is unchanged), writes the interindividual-only sites sidecar
//          from the combined panel (panel_id == 0), records per-thread
//          per-(cell,SNP) evidence at panel-0 sites in the same fixed-point
//          units as the single-panel path, and flushes observations after the
//          count merge. Both sidecars stage through the output transaction as
//          before. Two fail-closed guards added in main: --dump_pileup with
//          loaded counts (reuse path) and --dump_pileup with --threads <= 1
//          are refused, because both paths would otherwise publish a bundle
//          without the promised pileup files.
//   V2.12  Harden the production demux_parallel path: strict half-open alignment
//          mapping, fail-closed parallel counting,
//          authoritative pairwise diagnostics, complete-edge maximin acceptance,
//          transactional output publication, and strict preservation of the
//          declared -I identity set as the only legal final-answer set.

#include <getopt.h>
#include <argp.h>
#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <errno.h>
#include <map>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <float.h>
#include <chrono>
#include <limits>
#include <stdexcept>
#include <omp.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <zlib.h>
#include <htswrapper/bc.h>
#include <htswrapper/bam.h>
#include <htswrapper/gzreader.h>
#include <mixtureDist/mixtureDist.h>
#include <mixtureDist/mixtureModel.h>
#include <mixtureDist/functions.h>
#include <optimML/multivar_ml.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "io.h"
#include "vcf_hts.h"
#include "genotype_llr.h"

using std::cout;
using std::endl;
using namespace std;

// Version/build information
#ifndef CELLBOUNCER_SOURCE_REVISION
#define CELLBOUNCER_SOURCE_REVISION "unknown"
#endif

const string VERSION = "2.13.3";
const string VERSION_MESSAGE = "parallel production demultiplexer with fail-closed counting and explicit comparison semantics";
const string VERSION_NEW = "v2.13.3: SNP-bearing BAM targets only, synchronous fallback for optional HTSlib helper threads, and exit-safe staged-output cleanup; includes the v2.13.2 index-statistics fix";

// Global verbose flag (defined in demux_parallel_llr.cpp)
extern bool g_verbose;

// Global debug flag (defined in demux_parallel_llr.cpp)
extern bool g_debug;

// Species panel mode (used in main and potentially in helper functions)
enum class SpeciesPanelMode { NONE, COUNT_ONLY, FILTER, AUGMENT, BOTH };


class OutputTransaction {
    struct PrefixPair {
        string final_prefix;
        string staged_prefix;
    };
    struct PublishedFile {
        string staged;
        string final_name;
        string backup;
        bool had_backup;
        bool installed;
    };

    vector<PrefixPair> prefixes_;
    string token_;
    bool committed_;
    OutputTransaction* previous_active_;

    static OutputTransaction* active_;
    static bool exit_handler_registered_;

    static void cleanup_active_at_exit() {
        if (active_ != NULL && !active_->committed_) {
            active_->cleanup_staged();
        }
    }

    static void split_prefix(const string& prefix, string& directory, string& base) {
        const size_t slash = prefix.find_last_of('/');
        if (slash == string::npos) {
            directory = ".";
            base = prefix;
        } else {
            directory = slash == 0 ? "/" : prefix.substr(0, slash);
            base = prefix.substr(slash + 1);
        }
    }

    static string join_path(const string& directory, const string& name) {
        if (directory == "/") return "/" + name;
        if (directory == ".") return name;
        return directory + "/" + name;
    }

    vector<PublishedFile> discover_files() const {
        vector<PublishedFile> files;
        set<string> staged_seen;
        for (const PrefixPair& pair : prefixes_) {
            string directory;
            string staged_base;
            string final_directory;
            string final_base;
            split_prefix(pair.staged_prefix, directory, staged_base);
            split_prefix(pair.final_prefix, final_directory, final_base);
            DIR* dirp = opendir(directory.c_str());
            if (dirp == NULL) continue;
            while (dirent* entry = readdir(dirp)) {
                const string name = entry->d_name;
                if (name.compare(0, staged_base.size(), staged_base) != 0) continue;
                const string staged_path = join_path(directory, name);
                if (!staged_seen.insert(staged_path).second) continue;
                const string suffix = name.substr(staged_base.size());
                PublishedFile file;
                file.staged = staged_path;
                file.final_name = join_path(final_directory, final_base + suffix);
                file.backup = file.final_name + ".backup.demux_parallel." + token_;
                file.had_backup = false;
                file.installed = false;
                files.push_back(file);
            }
            closedir(dirp);
        }
        sort(files.begin(), files.end(), [](const PublishedFile& a, const PublishedFile& b) {
            return a.final_name < b.final_name;
        });
        return files;
    }

    void cleanup_staged() const {
        const vector<PublishedFile> files = discover_files();
        for (const PublishedFile& file : files) unlink(file.staged.c_str());
    }

public:
    OutputTransaction() : committed_(false), previous_active_(active_) {
        const long long stamp = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now().time_since_epoch()).count();
        token_ = std::to_string((long long)getpid()) + "." + std::to_string(stamp);
        active_ = this;
        if (!exit_handler_registered_) {
            if (std::atexit(&OutputTransaction::cleanup_active_at_exit) == 0) {
                exit_handler_registered_ = true;
            } else {
                fprintf(stderr,
                    "WARNING: could not register emergency staged-output cleanup; "
                    "normal return-path cleanup remains active\n");
            }
        }
    }

    ~OutputTransaction() {
        if (!committed_) cleanup_staged();
        if (active_ == this) active_ = previous_active_;
    }

    string stage_prefix(const string& final_prefix) {
        for (const PrefixPair& pair : prefixes_) {
            if (pair.final_prefix == final_prefix) return pair.staged_prefix;
        }
        PrefixPair pair;
        pair.final_prefix = final_prefix;
        pair.staged_prefix = final_prefix + ".tmp.demux_parallel." + token_;
        prefixes_.push_back(pair);
        return pair.staged_prefix;
    }

    bool publish(string* error_message = NULL) {
        vector<PublishedFile> files = discover_files();
        if (files.empty()) {
            if (error_message) *error_message = "no staged output files were produced";
            return false;
        }
        for (PublishedFile& file : files) {
            struct stat st;
            if (lstat(file.final_name.c_str(), &st) == 0) {
                unlink(file.backup.c_str());
                if (rename(file.final_name.c_str(), file.backup.c_str()) != 0) {
                    if (error_message) {
                        *error_message = "could not back up " + file.final_name + ": " + strerror(errno);
                    }
                    for (PublishedFile& rollback : files) {
                        if (rollback.had_backup) rename(rollback.backup.c_str(), rollback.final_name.c_str());
                    }
                    return false;
                }
                file.had_backup = true;
            }
        }

        for (PublishedFile& file : files) {
            if (rename(file.staged.c_str(), file.final_name.c_str()) != 0) {
                if (error_message) {
                    *error_message = "could not publish " + file.final_name + ": " + strerror(errno);
                }
                for (PublishedFile& rollback : files) {
                    if (rollback.installed) unlink(rollback.final_name.c_str());
                }
                for (PublishedFile& rollback : files) {
                    if (rollback.had_backup) rename(rollback.backup.c_str(), rollback.final_name.c_str());
                }
                return false;
            }
            file.installed = true;
        }

        for (const PublishedFile& file : files) {
            if (file.had_backup) unlink(file.backup.c_str());
        }
        committed_ = true;
        return true;
    }
};

OutputTransaction* OutputTransaction::active_ = NULL;
bool OutputTransaction::exit_handler_registered_ = false;

static string optional_number(double value) {
    if (!std::isfinite(value)) return "NA";
    char buffer[128];
    snprintf(buffer, sizeof(buffer), "%.17g", value);
    return string(buffer);
}

static string optional_identity(int identity, vector<string>& samples) {
    return identity >= 0 ? idx2name(identity, samples) : "NA";
}

static string join_identity_names(const vector<int>& identities, vector<string>& samples) {
    if (identities.empty()) return "none";
    string result;
    for (size_t i = 0; i < identities.size(); ++i) {
        if (i > 0) result += ",";
        result += idx2name(identities[i], samples);
    }
    return result;
}

static bool publish_outputs(OutputTransaction& transaction, const char* context) {
    string error;
    if (!transaction.publish(&error)) {
        fprintf(stderr, "ERROR: could not publish %s outputs atomically: %s\n",
            context, error.c_str());
        return false;
    }
    return true;
}

struct PanelUsabilitySummary {
    size_t nonempty_panel_contigs;
    size_t shared_bam_contigs;
    size_t usable_snps;
    size_t removed_nonshared_contigs;

    PanelUsabilitySummary()
        : nonempty_panel_contigs(0), shared_bam_contigs(0), usable_snps(0),
          removed_nonshared_contigs(0) {}
};

static bool panel_snp_is_usable(const SNPData& snp) {
    // A loaded record whose every sample has unavailable GT is not usable
    // evidence. Keep partially covered sites; downstream comparison
    // completeness decides whether a particular identity edge is evaluable.
    return snp.data.haps_covered.any();
}

static size_t count_panel_snps(
    const robin_hood::unordered_map<int, ChromSNPs>& panel,
    size_t* nonempty_contigs = NULL) {
    size_t contigs = 0;
    size_t snps = 0;
    for (const auto& kv : panel) {
        size_t usable_on_contig = 0;
        for (const SNPData& snp : kv.second.snps) {
            if (panel_snp_is_usable(snp)) ++usable_on_contig;
        }
        if (usable_on_contig == 0) continue;
        ++contigs;
        snps += usable_on_contig;
    }
    if (nonempty_contigs != NULL) *nonempty_contigs = contigs;
    return snps;
}

static bool validate_panel_without_bam(
    const robin_hood::unordered_map<int, ChromSNPs>& panel,
    const string& panel_label,
    bool allow_empty) {
    size_t nonempty_contigs = 0;
    const size_t usable_snps = count_panel_snps(panel, &nonempty_contigs);
    if (nonempty_contigs == 0 || usable_snps == 0) {
        if (allow_empty) {
            fprintf(stderr,
                "WARNING: optional %s contains no usable SNPs after loading/filtering; "
                "diagnostic panel will be unavailable\n",
                panel_label.c_str());
            return true;
        }
        fprintf(stderr,
            "ERROR: required %s contains zero usable SNPs after loading/filtering\n",
            panel_label.c_str());
        return false;
    }
    fprintf(stderr, "%s usability: %lu nonempty contigs, %lu usable SNPs\n",
        panel_label.c_str(), (unsigned long)nonempty_contigs,
        (unsigned long)usable_snps);
    return true;
}

static bool validate_and_restrict_panel_to_bam(
    robin_hood::unordered_map<int, ChromSNPs>& panel,
    const map<string, int>& bam_seq2tid,
    const string& panel_label,
    bool allow_empty) {
    set<int> bam_tids;
    for (const auto& kv : bam_seq2tid) bam_tids.insert(kv.second);

    PanelUsabilitySummary summary;
    for (auto it = panel.begin(); it != panel.end();) {
        vector<SNPData>& snps = it->second.snps;
        const bool panel_contig_has_records = !snps.empty();
        if (panel_contig_has_records) {
            ++summary.nonempty_panel_contigs;
            if (bam_tids.count(it->first) > 0) {
                // Track physical BAM/panel overlap independently from whether
                // any records retain usable GT data. This keeps the two
                // release-blocker states distinguishable.
                ++summary.shared_bam_contigs;
            }
        }

        snps.erase(
            std::remove_if(snps.begin(), snps.end(),
                [](const SNPData& snp) { return !panel_snp_is_usable(snp); }),
            snps.end());
        if (bam_tids.count(it->first) == 0) {
            if (panel_contig_has_records) ++summary.removed_nonshared_contigs;
            panel.erase(it++);
            continue;
        }
        if (snps.empty()) {
            panel.erase(it++);
            continue;
        }
        summary.usable_snps += snps.size();
        ++it;
    }

    if (summary.shared_bam_contigs == 0) {
        if (allow_empty) {
            fprintf(stderr,
                "WARNING: optional %s has zero shared BAM/panel contigs; "
                "diagnostic panel will be unavailable\n",
                panel_label.c_str());
            return true;
        }
        fprintf(stderr,
            "ERROR: required %s has zero shared BAM/panel contigs\n",
            panel_label.c_str());
        return false;
    }
    if (summary.usable_snps == 0) {
        if (allow_empty) {
            fprintf(stderr,
                "WARNING: optional %s has zero usable SNPs after loading/filtering; "
                "diagnostic panel will be unavailable\n",
                panel_label.c_str());
            return true;
        }
        fprintf(stderr,
            "ERROR: required %s contains zero usable SNPs after loading/filtering\n",
            panel_label.c_str());
        return false;
    }

    fprintf(stderr,
        "%s usability: %lu shared BAM/panel contigs, %lu usable SNPs",
        panel_label.c_str(), (unsigned long)summary.shared_bam_contigs,
        (unsigned long)summary.usable_snps);
    if (summary.removed_nonshared_contigs > 0) {
        fprintf(stderr, " (%lu nonshared contigs ignored)",
            (unsigned long)summary.removed_nonshared_contigs);
    }
    fprintf(stderr, "\n");
    return true;
}

static set<string> shared_contig_names(
    const set<string>& panel_contigs,
    const map<string, int>& bam_seq2tid) {
    set<string> shared;
    for (const string& contig : panel_contigs) {
        if (bam_seq2tid.count(contig) > 0) shared.insert(contig);
    }
    return shared;
}

static bool write_samples_checked(const string& filename, const vector<string>& samples) {
    FILE* outf = fopen(filename.c_str(), "w");
    if (outf == NULL) {
        fprintf(stderr, "ERROR: could not open %s for writing: %s\n",
            filename.c_str(), strerror(errno));
        return false;
    }
    bool ok = true;
    for (const string& sample : samples) {
        if (fprintf(outf, "%s\n", sample.c_str()) < 0) {
            ok = false;
            break;
        }
    }
    if (ferror(outf)) ok = false;
    if (fclose(outf) != 0) ok = false;
    if (!ok) {
        fprintf(stderr, "ERROR: failed writing %s: %s\n",
            filename.c_str(), strerror(errno));
    }
    return ok;
}

/**
 * Log likelihood function for computing error rates
 */
double ll_err(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p0 = data_d.at("exp");
    double e_r = params[0];
    double e_a = params[1];
    double p = p0 - p0*e_a + (1.0 - p0)*e_r;
    if (p <= 0){
        p = DBL_MIN*1e6;
    }
    else if (p >= 1){
        p = 1.0-DBL_MIN*1e6;
    }
    double ll = dbinom(n, k, p)/log2(exp(1));
    return ll;
}

void dll_err(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i, vector<double>& results){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double e_r = params[0];
    double e_a = params[1];
    double p0 = data_d.at("exp");
    double p = p0 - p0*e_a + (1.0 - p0)*e_r;
    if (p <= 0){
        p = DBL_MIN*1e6;
    }
    else if (p >= 1.0){
        p = 1.0-DBL_MIN*1e6;
    }
    
    double dy_dp = (k-n*p)/(p-p*p);
    double dp_de_a = -p0;
    double dp_de_r = 1.0 - p0;
    results[0] = dy_dp * dp_de_r;
    results[1] = dy_dp * dp_de_a;
}

/**
 * Re-infer error rates from assignments using optimized structure
 */
pair<double, double> infer_error_rates_optimized(
    robin_hood::unordered_map<unsigned long, CellCounts>& cell_counts,
    int n_samples,
    robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    double error_ref,
    double error_alt,
    double error_sigma){
   
    vector<double> n;
    vector<double> k;
    vector<double> expected;
    vector<double> weights_llr;

    int doublet_points = 0, singlet_points = 0;
    for (auto& a : assn){
        double weight = assn_llr[a.first];
        const CellCounts& counts = cell_counts[a.first];
        
        bool is_combo = a.second >= n_samples;
        pair<int, int> combo;
        if (is_combo){
            combo = idx_to_hap_comb(a.second, n_samples);
        }
        
        // Debug: show assignment details for each cell
        if (g_debug){
            fprintf(stderr, "DEBUG cell: assn=%d n_samples=%d is_combo=%d", 
                    a.second, n_samples, is_combo);
            if (is_combo) {
                fprintf(stderr, " combo=(%d,%d)", combo.first, combo.second);
            }
            fprintf(stderr, "\n");
        }
        
        if (is_combo){
            // Doublet - use pairwise counts for the two component individuals
            int i = combo.first;
            int j = combo.second;
            
            for (int nalt_i = 0; nalt_i < 3; ++nalt_i){
                for (int nalt_j = 0; nalt_j < 3; ++nalt_j){
                    auto pair_counts = counts.get(i, nalt_i, j, nalt_j);
                    float ref_count = pair_counts.first;
                    float alt_count = pair_counts.second;
                    
                    if (ref_count + alt_count > 0){
                        double this_expected = (double)(nalt_i + nalt_j) / 4.0;
                        expected.push_back(this_expected);
                        n.push_back(ref_count + alt_count);
                        k.push_back(alt_count);
                        weights_llr.push_back(weight);
                        doublet_points++;
                    }
                }
            }
        }
        else{
            // Singlet - use total counts for this individual
            for (int nalt = 0; nalt < 3; ++nalt){
                auto total = counts.get_total(a.second, nalt);
                float ref_count = total.first;
                float alt_count = total.second;
                
                if (ref_count + alt_count > 0){
                    double this_expected = (double)nalt / 2.0;
                    expected.push_back(this_expected);
                    n.push_back(ref_count + alt_count);
                    k.push_back(alt_count);
                    weights_llr.push_back(weight);
                    singlet_points++;
                }
            }
        }
    }
    if (g_debug){
        fprintf(stderr, "DEBUG: doublet_points=%d singlet_points=%d\n", doublet_points, singlet_points);
    }

    if (n.empty()){
        return make_pair(error_ref, error_alt);
    }

    // Diagnostic: print totals going into optimizer
    double total_n = 0, total_k = 0, total_weight = 0;
    std::map<double, std::pair<double, double>> by_expected;  // expected -> (n, k)
    for (size_t i = 0; i < n.size(); i++) {
        total_n += n[i];
        total_k += k[i];
        total_weight += weights_llr[i];
        by_expected[expected[i]].first += n[i];
        by_expected[expected[i]].second += k[i];
    }
    if (g_debug){
        fprintf(stderr, "DEBUG infer_error_rates_optimized:\n");
        fprintf(stderr, "  data_points=%lu total_n=%.2f total_k=%.2f total_weight=%.2f\n", 
                n.size(), total_n, total_k, total_weight);
        fprintf(stderr, "  overall_alt_frac=%.6f\n", total_k / total_n);
        for (auto& e : by_expected) {
            fprintf(stderr, "  expected=%.2f: n=%.2f k=%.2f alt_frac=%.6f\n",
                    e.first, e.second.first, e.second.second, 
                    e.second.second / e.second.first);
        }
    }

    optimML::multivar_ml_solver solver({error_ref, error_alt}, ll_err, dll_err);
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("exp", expected);
    solver.add_weights(weights_llr);
    solver.constrain_01(0);
    solver.constrain_01(1);
    solver.add_normal_prior(0, error_ref, error_sigma, 0.0, 1.0);
    solver.add_normal_prior(1, error_alt, error_sigma, 0.0, 1.0);
    
    solver.set_silent(true);

    double sigma_curr = error_sigma;
    bool success = false;
    while (!success){
        success = true;
        try{ 
            solver.solve();
        } 
        catch (const int& err){
            if (err == optimML::OPTIMML_MATH_ERR){
                sigma_curr *= 0.5;
                fprintf(stderr, "Decreasing prior sd to %f...\n", sigma_curr);
                solver.set_prior_param(0, "sigma", sigma_curr);
                solver.set_prior_param(1, "sigma", sigma_curr);
                solver.set_param(0, error_ref);
                solver.set_param(1, error_alt);        
                success = false;
            }
            else{
                fprintf(stderr, "Unknown error encountered while inferring error rates\n");
                return make_pair(error_ref, error_alt);
            }
        }
    }
    return make_pair(solver.results[0], solver.results[1]);
}

/**
 * Mann-Whitney U test for comparing LLR distributions between identities.
 * Returns the p-value for testing whether id1's LLRs are lower than id2's (or all others if id2 == -1).
 */
double mannwhitney_llr(int id1, 
    int id2, 
    robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr){
    
    double n1 = 0.0;
    double n2 = 0.0;

    vector<pair<double, int> > llrs;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){
        if (assn_llr[a->first] > 0){
            if (a->second == id1){
                llrs.push_back(make_pair(assn_llr[a->first], id1));
                n1++;
            }
            else if (id2 == -1 || a->second == id2){
                llrs.push_back(make_pair(assn_llr[a->first], id2));
                n2++;
            }
        }
    }
    
    if (n1 == 0){
        return 0.0;
    }
    else if (n2 == 0){
        return 1.0;
    }

    sort(llrs.begin(), llrs.end());
    
    // Assign ranks
    vector<double> ranks;
    int rank = 1;
    for (int i = 0; i < llrs.size(); ++i){
        ranks.push_back((double)rank);
        rank++;
    }

    // Deal with ties in ranks
    double prevllr = 0.0;
    double ranksum = 0.0;
    int ranktot = 0;
    bool ties = false;

    vector<double> nties;

    for (int i = 0; i < llrs.size(); ++i){
        if (llrs[i].first == prevllr){
            ranksum += ranks[i];
            ranktot++;
        }
        else{
            if (ranktot > 1){
                nties.push_back((double)ranktot);
                ties = true;
                double rankmean = ranksum / (double)ranktot;
                for (int j = i - 1; j >= i - 1 - (ranktot-1); --j){
                    ranks[j] = rankmean;
                }
                ranksum = 0.0;
                ranktot = 0;
            }
            else{
                if (i > 0){
                    nties.push_back(0);
                    ranks[i-1] = ranksum;
                }
                ranksum = 0;
                ranktot = 0;
            }
            ranksum += ranks[i];
            ranktot++;

        }
        prevllr = llrs[i].first;    
    }
    // Handle last one.
    if (ranktot > 1){
        ties = true;
        nties.push_back(ranktot);
        double rankmean = ranksum / (double)ranktot;
        for (int j = llrs.size()-1; j >= (int)llrs.size()-1 - (ranktot-1); --j){
            ranks[j] = rankmean;
        }
    }
    else{
        nties.push_back(0);
        ranks[llrs.size()-1] = ranksum;
    }
    
    double sum_id1 = 0;
    double sum_id2 = 0;
    for (int i = 0; i < ranks.size(); ++i){
        if (llrs[i].second == id1){
            sum_id1 += ranks[i];
        }
        else{
            sum_id2 += ranks[i];
        }
    } 
    
    double U1 = sum_id1 - (n1*(n1+1))/2.0;
    double U2 = sum_id2 - (n2*(n2+1))/2.0;
    double m_u = (n1*n2)/2.0;
    double sigma_u = sqrt((n1*n2*(n1+n2 +1))/(12.0));
    if (ties){
        double tsum = 0;
        for (int i = 0; i < nties.size(); ++i){
            if (nties[i] > 0){
                tsum += (pow(nties[i], 3) - nties[i]);
            }
        }
        double term1 = (n1*n2*(n1+n2+1))/12;
        double term2 = (n1*n2*tsum)/(12*(n1+n2)*(n1+n2-1));
        sigma_u = sqrt(term1-term2);
    }
    
    if (n1 < 3 || n2 < 3){
        // Can't get a reliable sigma here.
        if (U1 < m_u){
            return 0.0;
        }
        else{
            return 1.0;
        }
    }

    // Test whether id1 < id2
    double p = pnorm(U1, m_u, sigma_u);
    return p;
}

/**
 * Does some QC on assignments - for each ID in the assignment file,
 * checks for significantly lower LLR distribution than the rest of
 * the data set. Also checks for significantly lower numbers of cells.
 */
void id_qc(robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    map<int, double>& pois_p,
    map<int, double>& mannwhitney_p){
    
    // Get total num cells for each ID 
    map<int, int> idsizes;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){
        if (assn_llr[a->first] > 0){
            if (idsizes.count(a->second) == 0){
                idsizes.insert(make_pair(a->second, 0));
            }
            idsizes[a->second]++;
        }
    }

    for (map<int, int>::iterator ids = idsizes.begin(); ids != idsizes.end(); ++ids){
        double mean_othersize = 0.0;
        double mean_othersize_tot = 0.0;
        for (map<int, int>::iterator ids2 = idsizes.begin(); ids2 != idsizes.end(); ++ids2){
            if (ids2->first != ids->first){
                mean_othersize += (double)ids2->second;
                mean_othersize_tot++;
            }
            
        }
        double p1 = ppois(ids->second, mean_othersize_tot);
        pois_p.insert(make_pair(ids->first, p1));
        double p2 = mannwhitney_llr(ids->first, -1, assn, assn_llr);
        mannwhitney_p.insert(make_pair(ids->first, p2));
    } 
}

/**
 * The -I file is an exact final-answer contract. Component singlets may remain
 * available internally so declared fusion likelihoods can be constructed, but
 * they are never promoted into allowed_ids2 and can never become final calls.
 */

/**
 * Dump cell counts from optimized CellCounts structure to file.
 * Format matches original: cell_barcode indv1 nalt1 indv2 nalt2 ref_count alt_count
 */
void dump_cellcounts_optimized(gzFile& out_cell,
    robin_hood::unordered_map<unsigned long, CellCounts>& cell_counts,
    int n_samples){
    
    char linebuf[1024];

    for (auto& cell : cell_counts){
        unsigned long bc = cell.first;
        const CellCounts& counts = cell.second;
        
        // Output totals (indv2 = -1, nalt2 = -1)
        for (int indv = 0; indv < n_samples; ++indv){
            for (int nalt = 0; nalt < 3; ++nalt){
                auto total = counts.get_total(indv, nalt);
                if (total.first > 0 || total.second > 0){
                    sprintf(&linebuf[0], "%lu\t%d\t%d\t%d\t%d\t%f\t%f\n", 
                        bc, indv, nalt, -1, -1, total.first, total.second);
                    gzwrite(out_cell, &linebuf[0], strlen(linebuf));
                }
            }
        }
        
        // Output pairwise counts
        for (int indv1 = 0; indv1 < n_samples; ++indv1){
            for (int nalt1 = 0; nalt1 < 3; ++nalt1){
                for (int indv2 = indv1 + 1; indv2 < n_samples; ++indv2){
                    for (int nalt2 = 0; nalt2 < 3; ++nalt2){
                        auto pair_counts = counts.get(indv1, nalt1, indv2, nalt2);
                        if (pair_counts.first > 0 || pair_counts.second > 0){
                            sprintf(&linebuf[0], "%lu\t%d\t%d\t%d\t%d\t%f\t%f\n", 
                                bc, indv1, nalt1, indv2, nalt2, 
                                pair_counts.first, pair_counts.second);
                            gzwrite(out_cell, &linebuf[0], strlen(linebuf));
                        }
                    }
                }
            }
        }
    } 
}


static bool text_file_contains(const string& filename, const string& token){
    std::ifstream in(filename.c_str());
    if (!in.good()) return false;
    string line;
    while (std::getline(in, line)){
        if (line.find(token) != string::npos) return true;
    }
    return false;
}

static bool validate_gzip_readable(const string& filename){
    gzFile test_gz = gzopen(filename.c_str(), "r");
    if (!test_gz){
        return false;
    }
    char buf[1024 * 64];
    bool gz_ok = true;
    while (true){
        int n_read = gzread(test_gz, buf, sizeof(buf));
        if (n_read < 0){
            gz_ok = false;
            break;
        }
        if (n_read == 0) break;
    }
    int close_code = gzclose(test_gz);
    return gz_ok && close_code == Z_OK;
}

/**
 * Load cell counts from gzipped file into CellCounts structures.
 * Reverses dump_cellcounts_optimized: reads 7-column format, scales
 * float values back to int64 fixed-point representation.
 *
 * Format per line: barcode indv1 nalt1 indv2 nalt2 ref_count alt_count
 * When indv2==-1 && nalt2==-1: total row -> add_total(indv1, nalt1, ...)
 * Otherwise: pairwise row -> add(indv1, nalt1, indv2, nalt2, ...)
 */
int load_cellcounts_optimized(const string& filename,
    robin_hood::unordered_map<unsigned long, CellCounts>& cell_counts,
    int n_samples){
    
    gzreader reader(filename);
    int n_lines = 0;
    
    while(reader.next()){
        istringstream splitter(reader.line);
        string field;
        int idx = 0;
        
        unsigned long bc = 0;
        int indv1 = 0, nalt1 = 0, indv2 = 0, nalt2 = 0;
        float ref_val = 0.0f, alt_val = 0.0f;
        
        while(getline(splitter, field, '\t')){
            switch(idx){
                case 0: bc = strtoul(field.c_str(), NULL, 10); break;
                case 1: indv1 = atoi(field.c_str()); break;
                case 2: nalt1 = atoi(field.c_str()); break;
                case 3: indv2 = atoi(field.c_str()); break;
                case 4: nalt2 = atoi(field.c_str()); break;
                case 5: ref_val = atof(field.c_str()); break;
                case 6: alt_val = atof(field.c_str()); break;
            }
            idx++;
        }
        
        if (idx < 7) continue;
        
        // Initialize CellCounts for this barcode if needed
        if (cell_counts.count(bc) == 0){
            cell_counts.emplace(bc, CellCounts(n_samples));
        }
        
        // Convert float values back to fixed-point int64
        int64_t ref_scaled = (int64_t)round((double)ref_val * FIXED_POINT_SCALE);
        int64_t alt_scaled = (int64_t)round((double)alt_val * FIXED_POINT_SCALE);
        
        if (indv2 == -1 && nalt2 == -1){
            // Total row
            cell_counts[bc].add_total(indv1, nalt1, ref_scaled, alt_scaled);
        }
        else{
            // Pairwise row
            cell_counts[bc].add(indv1, nalt1, indv2, nalt2, ref_scaled, alt_scaled);
        }
        
        n_lines++;
    }
    
    return n_lines;
}


// ============================================================================
// Native species artifact helpers (Path Separation Contract V1_R3)
// ============================================================================

static int64_t scaled_from_float(double x){
    return (int64_t)llround(x * (double)FIXED_POINT_SCALE);
}

static vector<vector<pair<int, double>>> build_indiv_to_species_weights(
    const PanelMetadata& pm,
    int n_indiv){

    vector<vector<pair<int, double>>> out(n_indiv);
    for (int sp_idx = 0; sp_idx < (int)pm.species_list.size(); ++sp_idx){
        const string& sp = pm.species_list[sp_idx];
        auto it = pm.species_to_sample_indices.find(sp);
        if (it == pm.species_to_sample_indices.end()) continue;
        for (int indiv_idx : it->second){
            if (indiv_idx < 0 || indiv_idx >= n_indiv) continue;
            double w = pm.get_weight(sp, indiv_idx);
            if (w <= 0.0) continue;
            out[indiv_idx].push_back(make_pair(sp_idx, w));
        }
    }
    return out;
}


struct NativeSpeciesNormalization {
    vector<double> singlet_weight;
    vector<vector<double>> pair_weight;
};

/**
 * Compute the total panel-membership weight represented by each native species
 * singlet row and each heterotypic species-pair row.
 *
 * These are the exact multiplicities introduced when individual-native count
 * rows are folded into species-native rows.  Hybrid samples can contribute to
 * multiple species with fractional weights.  For pair rows, only distinct
 * individual pairs are represented by CellCounts, so the normalization is
 * accumulated over i<j rather than using a simple product of species totals.
 */
static NativeSpeciesNormalization build_native_species_normalization(
    const vector<vector<pair<int, double>>>& i2sp,
    int n_species){

    NativeSpeciesNormalization norm;
    norm.singlet_weight.assign(n_species, 0.0);
    norm.pair_weight.assign(n_species, vector<double>(n_species, 0.0));

    for (int i = 0; i < (int)i2sp.size(); ++i){
        for (const auto& mi : i2sp[i]){
            if (mi.first < 0 || mi.first >= n_species || mi.second <= 0.0) continue;
            norm.singlet_weight[mi.first] += mi.second;
        }
    }

    for (int i = 0; i < (int)i2sp.size(); ++i){
        for (int j = i + 1; j < (int)i2sp.size(); ++j){
            for (const auto& mi : i2sp[i]){
                for (const auto& mj : i2sp[j]){
                    if (mi.first == mj.first || mi.second <= 0.0 || mj.second <= 0.0) continue;
                    int a = std::min(mi.first, mj.first);
                    int b = std::max(mi.first, mj.first);
                    norm.pair_weight[a][b] += mi.second * mj.second;
                }
            }
        }
    }

    for (int sp = 0; sp < n_species; ++sp){
        if (norm.singlet_weight[sp] <= 0.0){
            fprintf(stderr,
                "ERROR: native species %d has zero panel-membership weight; cannot normalize counts\n",
                sp);
            exit(1);
        }
    }
    for (int a = 0; a < n_species; ++a){
        for (int b = a + 1; b < n_species; ++b){
            if (norm.pair_weight[a][b] <= 0.0){
                fprintf(stderr,
                    "WARNING: native species pair %d+%d has zero distinct-individual panel weight; "
                    "its pair rows will remain empty\n",
                    a, b);
            }
        }
    }
    return norm;
}

static void report_native_species_normalization(
    const NativeSpeciesNormalization& norm,
    const vector<string>& species_names,
    const char* context){

    fprintf(stderr, "Native species normalization (%s):\n", context);
    for (int sp = 0; sp < (int)species_names.size(); ++sp){
        fprintf(stderr, "  singlet %s weight %.9g\n",
            species_names[sp].c_str(), norm.singlet_weight[sp]);
    }
    for (int a = 0; a < (int)species_names.size(); ++a){
        for (int b = a + 1; b < (int)species_names.size(); ++b){
            fprintf(stderr, "  pair %s+%s weight %.9g\n",
                species_names[a].c_str(), species_names[b].c_str(),
                norm.pair_weight[a][b]);
        }
    }
}

/**
 * Convert a legacy, panel-multiplicity-inflated native species count bundle to
 * the normalized v2.09 representation in memory.  This is used only by the
 * no-BAM reassignment diagnostic when --panel_metadata is supplied.
 */
struct NativeSpeciesTargetStats {
    uint64_t n_sites = 0;
    uint64_t n_sites_with_any_species = 0;
    vector<uint64_t> singlet_available_sites;
    vector<vector<uint64_t>> pair_available_sites;
};

/**
 * Precompute native-species accumulation targets independently at every SNP.
 *
 * Each species singlet and heterotypic species-pair hypothesis receives one
 * normalized unit of evidence per observed SNP, distributed across genotype
 * states according to the panel members that have valid genotypes at that SNP.
 * This makes count evidence invariant to unequal panel representation and to
 * site-specific missing genotypes.  It is the count-side analogue of the
 * target-species normalization used by compute_species_condf_native().
 */
static NativeSpeciesTargetStats precompute_native_species_targets(
    robin_hood::unordered_map<int, ChromSNPs>& species_snpdat,
    const PanelMetadata& pm,
    int n_indiv,
    NativeSpeciesTargetTable& target_table,
    int required_panel_id = -1){

    const int n_species = (int)pm.species_list.size();
    const int n_pairs = n_species * (n_species - 1) / 2;
    const uint64_t singlet_bins = (uint64_t)n_species * GENOTYPE_STATES;
    const uint64_t pair_bins =
        (uint64_t)n_pairs * GENOTYPE_STATES * GENOTYPE_STATES;
    const uint64_t weights_per_site = singlet_bins + pair_bins;
    const vector<vector<pair<int, double>>> i2sp =
        build_indiv_to_species_weights(pm, n_indiv);

    if (n_species <= 0 || weights_per_site == 0 ||
        weights_per_site > (uint64_t)std::numeric_limits<uint32_t>::max()){
        throw std::runtime_error(
            "invalid native-species target dimensions for per-site normalization");
    }

    auto pair_order = [n_species](int a, int b) -> int {
        return a * (2 * n_species - a - 1) / 2 + (b - a - 1);
    };

    target_table.clear();

    NativeSpeciesTargetStats stats;
    stats.singlet_available_sites.assign(n_species, 0);
    stats.pair_available_sites.assign(
        n_species, vector<uint64_t>(n_species, 0));

    for (auto& chrom_kv : species_snpdat){
        const int tid = chrom_kv.first;
        ChromSNPs& chrom = chrom_kv.second;
        NativeSpeciesChromTargets chrom_targets;
        chrom_targets.n_species = n_species;
        chrom_targets.n_pairs = n_pairs;
        chrom_targets.weights_per_site = (uint32_t)weights_per_site;
        chrom_targets.site_offsets.assign(chrom.snps.size(), UINT64_MAX);

        size_t qualifying_sites = 0;
        for (const auto& snp : chrom.snps){
            if (required_panel_id < 0 || (int)snp.panel_id == required_panel_id){
                ++qualifying_sites;
            }
        }
        if (qualifying_sites > 0){
            const uint64_t reserve_elems =
                (uint64_t)qualifying_sites * weights_per_site;
            if (reserve_elems > (uint64_t)chrom_targets.weights.max_size()){
                throw std::runtime_error(
                    "native-species target table exceeds vector capacity");
            }
            chrom_targets.weights.reserve((size_t)reserve_elems);
        }

        // Reused per-SNP work buffers. species_geno_weight[sp,g] is the
        // available panel-membership mass in one genotype bin. The pair table
        // is then obtained analytically as the outer product of two species'
        // genotype masses, minus same-individual hybrid membership. This is
        // exactly the distinct-individual sum used by the historical pair loop
        // without its O(n_individuals^2) cost.
        vector<double> species_geno_weight((size_t)singlet_bins, 0.0);
        vector<double> same_indiv_cross(
            (size_t)n_pairs * GENOTYPE_STATES, 0.0);
        vector<double> pair_raw(
            GENOTYPE_STATES * GENOTYPE_STATES, 0.0);
        vector<int32_t> block((size_t)weights_per_site, 0);

        auto scale_hypothesis = [&](const double* raw, int n_bins,
                                    double denom, size_t block_offset){
            if (denom <= 0.0) return;
            int best_bin = -1;
            double best_weight = -1.0;
            int64_t scaled_sum = 0;
            for (int bin = 0; bin < n_bins; ++bin){
                double weight = raw[bin] / denom;
                if (weight <= 0.0) continue;
                const int64_t scaled = scaled_from_float(weight);
                block[block_offset + (size_t)bin] = (int32_t)scaled;
                scaled_sum += scaled;
                if (weight > best_weight){
                    best_weight = weight;
                    best_bin = bin;
                }
            }
            if (best_bin >= 0){
                const int64_t corrected =
                    (int64_t)block[block_offset + (size_t)best_bin] +
                    (FIXED_POINT_SCALE - scaled_sum);
                block[block_offset + (size_t)best_bin] = (int32_t)corrected;
            }
        };

        for (size_t snp_index = 0; snp_index < chrom.snps.size(); ++snp_index){
            SNPData& snp = chrom.snps[snp_index];
            if (required_panel_id >= 0 && (int)snp.panel_id != required_panel_id){
                continue;
            }
            if ((int)snp.geno.size() != n_indiv){
                snp.precompute_genotypes(n_indiv);
            }
            ++stats.n_sites;

            std::fill(
                species_geno_weight.begin(), species_geno_weight.end(), 0.0);
            std::fill(same_indiv_cross.begin(), same_indiv_cross.end(), 0.0);
            std::fill(block.begin(), block.end(), 0);

            for (int i = 0; i < n_indiv; ++i){
                const int gi = (int)snp.geno[i];
                if (gi < 0 || gi >= GENOTYPE_STATES) continue;

                for (const auto& membership : i2sp[i]){
                    const int sp = membership.first;
                    const double weight = membership.second;
                    if (sp < 0 || sp >= n_species || weight <= 0.0) continue;
                    const size_t bin =
                        (size_t)sp * GENOTYPE_STATES + (size_t)gi;
                    species_geno_weight[bin] += weight;
                }

                for (size_t m = 0; m < i2sp[i].size(); ++m){
                    for (size_t n = m + 1; n < i2sp[i].size(); ++n){
                        const auto& left = i2sp[i][m];
                        const auto& right = i2sp[i][n];
                        if (left.first == right.first || left.second <= 0.0 ||
                            right.second <= 0.0){
                            continue;
                        }
                        const int a = std::min(left.first, right.first);
                        const int b = std::max(left.first, right.first);
                        const int ord = pair_order(a, b);
                        same_indiv_cross[
                            (size_t)ord * GENOTYPE_STATES + (size_t)gi] +=
                            left.second * right.second;
                    }
                }
            }

            bool any_hypothesis = false;
            bool any_species = false;

            for (int sp = 0; sp < n_species; ++sp){
                const double* raw = species_geno_weight.data() +
                    (size_t)sp * GENOTYPE_STATES;
                double denom = 0.0;
                for (int g = 0; g < GENOTYPE_STATES; ++g) denom += raw[g];
                if (denom <= 0.0) continue;
                ++stats.singlet_available_sites[sp];
                any_species = true;
                any_hypothesis = true;
                scale_hypothesis(
                    raw, GENOTYPE_STATES, denom,
                    (size_t)sp * GENOTYPE_STATES);
            }
            if (any_species) ++stats.n_sites_with_any_species;

            int ord = 0;
            for (int a = 0; a < n_species; ++a){
                for (int b = a + 1; b < n_species; ++b, ++ord){
                    std::fill(pair_raw.begin(), pair_raw.end(), 0.0);
                    double denom = 0.0;
                    for (int ga = 0; ga < GENOTYPE_STATES; ++ga){
                        const double wa = species_geno_weight[
                            (size_t)a * GENOTYPE_STATES + (size_t)ga];
                        for (int gb = 0; gb < GENOTYPE_STATES; ++gb){
                            const double wb = species_geno_weight[
                                (size_t)b * GENOTYPE_STATES + (size_t)gb];
                            double raw = wa * wb;
                            if (ga == gb){
                                raw -= same_indiv_cross[
                                    (size_t)ord * GENOTYPE_STATES + (size_t)ga];
                            }
                            if (raw < 0.0 && raw > -1e-12) raw = 0.0;
                            if (raw <= 0.0) continue;
                            const int bin = ga * GENOTYPE_STATES + gb;
                            pair_raw[(size_t)bin] = raw;
                            denom += raw;
                        }
                    }
                    if (denom <= 0.0) continue;
                    ++stats.pair_available_sites[a][b];
                    any_hypothesis = true;
                    scale_hypothesis(
                        pair_raw.data(),
                        GENOTYPE_STATES * GENOTYPE_STATES,
                        denom,
                        (size_t)singlet_bins +
                            (size_t)ord * GENOTYPE_STATES * GENOTYPE_STATES);
                }
            }

            if (!any_hypothesis) continue;
            chrom_targets.site_offsets[snp_index] =
                (uint64_t)chrom_targets.weights.size();
            chrom_targets.weights.insert(
                chrom_targets.weights.end(), block.begin(), block.end());
        }

        target_table[tid] = std::move(chrom_targets);
    }

    fprintf(stderr,
        "Native species per-site targets: %llu sites (%llu with at least one species hypothesis), %d species, %llu fixed weights/site\n",
        (unsigned long long)stats.n_sites,
        (unsigned long long)stats.n_sites_with_any_species,
        n_species,
        (unsigned long long)weights_per_site);
    for (int sp = 0; sp < n_species; ++sp){
        fprintf(stderr, "  singlet %s available at %llu sites\n",
            pm.species_list[sp].c_str(),
            (unsigned long long)stats.singlet_available_sites[sp]);
    }
    for (int a = 0; a < n_species; ++a){
        for (int b = a + 1; b < n_species; ++b){
            fprintf(stderr, "  pair %s+%s available at %llu sites\n",
                pm.species_list[a].c_str(), pm.species_list[b].c_str(),
                (unsigned long long)stats.pair_available_sites[a][b]);
        }
    }
    return stats;
}

static void normalize_existing_species_counts_in_place(
    robin_hood::unordered_map<unsigned long, CellCounts>& species_counts,
    const NativeSpeciesNormalization& norm,
    int n_species){

    for (auto& kv : species_counts){
        CellCounts& counts = kv.second;
        if (counts.n_samples != n_species){
            fprintf(stderr,
                "ERROR: native species count dimension mismatch while normalizing: expected %d, got %d\n",
                n_species, counts.n_samples);
            exit(1);
        }

        for (int sp = 0; sp < n_species; ++sp){
            const double denom = norm.singlet_weight[sp];
            for (int g = 0; g < GENOTYPE_STATES; ++g){
                const int idx = sp * GENOTYPE_STATES + g;
                counts.total_ref[idx] = (int64_t)llround((double)counts.total_ref[idx] / denom);
                counts.total_alt[idx] = (int64_t)llround((double)counts.total_alt[idx] / denom);
            }
        }

        for (int a = 0; a < n_species; ++a){
            for (int b = a + 1; b < n_species; ++b){
                const double denom = norm.pair_weight[a][b];
                if (denom <= 0.0) continue;
                for (int ga = 0; ga < GENOTYPE_STATES; ++ga){
                    const int idx_a = a * GENOTYPE_STATES + ga;
                    for (int gb = 0; gb < GENOTYPE_STATES; ++gb){
                        const int idx_b = b * GENOTYPE_STATES + gb;
                        const size_t flat = (size_t)idx_a * counts.state_count + idx_b;
                        counts.ref_counts[flat] =
                            (int64_t)llround((double)counts.ref_counts[flat] / denom);
                        counts.alt_counts[flat] =
                            (int64_t)llround((double)counts.alt_counts[flat] / denom);
                    }
                }
            }
        }
    }
}

static ConditionalWeightStats compute_species_condf_native(
    robin_hood::unordered_map<int, ChromSNPs>& species_snpdat,
    map<pair<int, int>, map<int, float>>& species_condf,
    const PanelMetadata& pm,
    int n_indiv,
    const AcceptedSiteWeightMap* accepted_site_weights = nullptr){

    ConditionalWeightStats stats;
    const int n_species = (int)pm.species_list.size();
    vector<vector<pair<int, double>>> i2sp = build_indiv_to_species_weights(pm, n_indiv);

    for (auto& kv : species_snpdat){
        for (auto& snp : kv.second.snps){
            if ((int)snp.geno.size() != n_indiv){
                snp.precompute_genotypes(n_indiv);
            }
        }
    }

    map<pair<int, int>, map<int, double>> sums;
    map<pair<int, int>, map<int, double>> tots;

    for (auto& chrom_kv : species_snpdat){
        const int tid = chrom_kv.first;
        for (auto& snp : chrom_kv.second.snps){
            double site_weight = 1.0;
            if (accepted_site_weights != nullptr){
                auto wit = accepted_site_weights->find(accepted_site_weight_key(tid, snp.pos));
                if (wit == accepted_site_weights->end() || wit->second <= 0) continue;
                site_weight = (double)wit->second;
            }
            ++stats.observed_sites;
            stats.accepted_weight += accepted_site_weights != nullptr
                ? site_weight / (double)FIXED_POINT_SCALE
                : 0.0;

            const int8_t* geno = snp.geno.data();

            vector<double> target_avg(n_species, 0.0);
            vector<double> target_wsum(n_species, 0.0);
            for (int j = 0; j < n_indiv; ++j){
                int8_t gj = geno[j];
                if (gj < 0 || gj >= GENOTYPE_STATES) continue;
                for (const auto& mj : i2sp[j]){
                    target_avg[mj.first] += mj.second * ((double)gj / 2.0);
                    target_wsum[mj.first] += mj.second;
                }
            }
            for (int sp = 0; sp < n_species; ++sp){
                if (target_wsum[sp] > 0.0) target_avg[sp] /= target_wsum[sp];
            }

            for (int i = 0; i < n_indiv; ++i){
                int8_t gi = geno[i];
                if (gi < 0 || gi >= GENOTYPE_STATES) continue;
                if (i2sp[i].empty()) continue;
                for (const auto& mi : i2sp[i]){
                    const double source_denom = target_wsum[mi.first];
                    if (source_denom <= 0.0 || mi.second <= 0.0) continue;
                    const double source_weight = mi.second / source_denom;
                    pair<int, int> row = make_pair(mi.first, (int)gi);
                    for (int sp_t = 0; sp_t < n_species; ++sp_t){
                        if (target_wsum[sp_t] <= 0.0) continue;
                        sums[row][sp_t] += site_weight * source_weight * target_avg[sp_t];
                        tots[row][sp_t] += site_weight * source_weight;
                    }
                }
            }
        }
    }

    species_condf.clear();
    for (auto& row_kv : sums){
        for (auto& col_kv : row_kv.second){
            double denom = tots[row_kv.first][col_kv.first];
            if (denom > 0.0){
                species_condf[row_kv.first][col_kv.first] = (float)(col_kv.second / denom);
            }
        }
    }

    fprintf(stderr, "Native species condf: %lu row entries, %d species, %llu observed sites, %.6f accepted weight\n",
        species_condf.size(), n_species,
        (unsigned long long)stats.observed_sites, stats.accepted_weight);
    return stats;
}

static void load_panel_metadata_if_needed(
    const string& panel_metadata_file,
    const vector<string>& samples,
    PanelMetadata& panel_meta,
    bool& panel_meta_loaded){

    if (!panel_meta_loaded){
        if (panel_metadata_file.empty()){
            fprintf(stderr, "ERROR: panel metadata required for native species artifacts\n");
            exit(1);
        }
        panel_meta = load_panel_metadata(panel_metadata_file, samples);
        panel_meta_loaded = true;
    }
}

static bool write_native_species_assignments(
    const string& output_prefix,
    robin_hood::unordered_map<unsigned long, CellCounts>& species_counts_native,
    const vector<string>& species_samples,
    const string& species_idfile,
    double doublet_rate,
    double error_ref,
    double error_alt,
    int n_threads,
    int n_target,
    const string& barcode_group,
    bool cellranger,
    bool seurat,
    bool underscore){

    if (species_counts_native.empty()){
        fprintf(stderr, "ERROR: native species assignment requested but no species counts are available\n");
        return false;
    }

    robin_hood::unordered_map<unsigned long, int> sp_assn;
    robin_hood::unordered_map<unsigned long, double> sp_assn_llr;
    set<int> allowed_ids;
    set<int> allowed_ids2;
    map<int, double> prior_weights;

    // Species assignment remains a separate native output. Candidate restriction
    // is the panel-agnostic analogue of -i/-I and does not alter individual calls.
    if (species_idfile.length() > 0){
        string idfile_copy = species_idfile;
        vector<string> samples_copy = species_samples;
        parse_idfile(idfile_copy, samples_copy, allowed_ids, allowed_ids2, false);
        if (allowed_ids.empty()){
            fprintf(stderr, "ERROR: no valid species labels found in %s; refusing to ignore --species_ids\n",
                species_idfile.c_str());
            return false;
        }
        fprintf(stderr, "Species candidate restriction active: %lu singlet label(s), "
            "%lu pair label(s) from %s\n",
            allowed_ids.size(), allowed_ids2.size(), species_idfile.c_str());
    }

    vector<string> tmp_species_samples = species_samples;
    if (!assign_ids_parallel(species_counts_native, tmp_species_samples,
            sp_assn, sp_assn_llr, allowed_ids, allowed_ids2, doublet_rate,
            error_ref, error_alt, false, prior_weights, n_threads, n_target)) {
        fprintf(stderr, "ERROR: native species identity scoring failed\n");
        return false;
    }

    const string fname = output_prefix + ".species_assignments";
    FILE* outf = fopen(fname.c_str(), "w");
    if (!outf){
        fprintf(stderr, "ERROR: could not open %s for writing: %s\n",
            fname.c_str(), strerror(errno));
        return false;
    }
    string tmp_group = barcode_group;
    dump_assignments(outf, sp_assn, sp_assn_llr, tmp_species_samples,
        tmp_group, cellranger, seurat, underscore);
    if (fclose(outf) != 0){
        fprintf(stderr, "ERROR: failed closing %s: %s\n", fname.c_str(), strerror(errno));
        return false;
    }
    fprintf(stderr, "Wrote native species assignments for %lu cells to %s\n",
        sp_assn.size(), fname.c_str());
    return true;
}


/**
 * Re-run native species assignment from an existing species-native count bundle.
 *
 * This diagnostic path deliberately reads only:
 *   INPUT_PREFIX.species_samples
 *   INPUT_PREFIX.species_counts
 * and writes a distinct output prefix. It never opens a BAM or VCF and never
 * modifies the source count bundle.
 */
static int reassign_native_species_from_existing(
    const string& input_prefix,
    const string& output_prefix,
    double species_doublet_rate,
    double error_ref,
    double error_alt,
    int n_threads,
    int n_target,
    int n_runner_ups,
    double close_threshold,
    const string& barcode_group,
    bool cellranger,
    bool seurat,
    bool underscore,
    const string& panel_metadata_file){

    if (input_prefix.empty()){
        fprintf(stderr, "ERROR: --species_reassign_from_prefix requires a non-empty prefix\n");
        return 1;
    }
    if (input_prefix == output_prefix){
        fprintf(stderr,
            "ERROR: diagnostic output prefix must differ from the source count-bundle prefix\n"
            "  source/output: %s\n",
            input_prefix.c_str());
        return 1;
    }
    if (species_doublet_rate < 0.0 || species_doublet_rate > 1.0){
        fprintf(stderr, "ERROR: --species_doublet_rate must be between 0 and 1\n");
        return 1;
    }

    string samples_file = input_prefix + ".species_samples";
    string counts_file = input_prefix + ".species_counts";
    if (!file_exists(samples_file)){
        fprintf(stderr, "ERROR: missing native species sample file: %s\n", samples_file.c_str());
        return 1;
    }
    if (!file_exists(counts_file)){
        fprintf(stderr, "ERROR: missing native species counts file: %s\n", counts_file.c_str());
        return 1;
    }
    if (!validate_gzip_readable(counts_file)){
        fprintf(stderr, "ERROR: native species counts file is not gzip-readable: %s\n", counts_file.c_str());
        return 1;
    }

    vector<string> species_samples;
    load_samples(samples_file, species_samples);
    if (species_samples.size() < 2){
        fprintf(stderr, "ERROR: expected at least two native species in %s; found %lu\n",
            samples_file.c_str(), species_samples.size());
        return 1;
    }
    string allocation_error;
    if (!validate_identity_and_allocation_request(
            (int)species_samples.size(), NULL, NULL, &allocation_error)){
        fprintf(stderr, "ERROR: native species identity universe is unsupported: %s\n",
            allocation_error.c_str());
        return 1;
    }

    robin_hood::unordered_map<unsigned long, CellCounts> species_counts;
    const int n_lines = load_cellcounts_optimized(
        counts_file, species_counts, species_samples.size());
    if (n_lines <= 0 || species_counts.empty()){
        fprintf(stderr, "ERROR: no native species counts loaded from %s\n", counts_file.c_str());
        return 1;
    }

    bool panel_normalized = false;
    NativeSpeciesNormalization reassign_norm;
    if (!panel_metadata_file.empty()){
        string indiv_samples_file = input_prefix + ".samples";
        if (!file_exists(indiv_samples_file)){
            fprintf(stderr,
                "ERROR: --panel_metadata normalization requires the original individual sample file: %s\n",
                indiv_samples_file.c_str());
            return 1;
        }
        vector<string> indiv_samples;
        load_samples(indiv_samples_file, indiv_samples);
        PanelMetadata pm = load_panel_metadata(panel_metadata_file, indiv_samples, true);
        if (pm.species_list != species_samples){
            fprintf(stderr,
                "ERROR: panel metadata native species order does not match %s\n",
                samples_file.c_str());
            return 1;
        }
        vector<vector<pair<int, double>>> i2sp =
            build_indiv_to_species_weights(pm, indiv_samples.size());
        reassign_norm = build_native_species_normalization(i2sp, species_samples.size());
        report_native_species_normalization(
            reassign_norm, species_samples, "legacy count-bundle diagnostic correction");
        normalize_existing_species_counts_in_place(
            species_counts, reassign_norm, species_samples.size());
        panel_normalized = true;
    }

    fprintf(stderr, "Native species reassignment diagnostic\n");
    fprintf(stderr, "  input prefix: %s\n", input_prefix.c_str());
    fprintf(stderr, "  output prefix: %s\n", output_prefix.c_str());
    fprintf(stderr, "  species: %lu\n", species_samples.size());
    fprintf(stderr, "  cells: %lu\n", species_counts.size());
    fprintf(stderr, "  count rows: %d\n", n_lines);
    fprintf(stderr, "  species doublet rate: %.9g\n", species_doublet_rate);
    fprintf(stderr, "  panel-normalized legacy counts: %s\n", panel_normalized ? "true" : "false");
    fprintf(stderr, "  error ref/alt: %.9g / %.9g\n", error_ref, error_alt);

    vector<unsigned long> barcodes;
    barcodes.reserve(species_counts.size());
    for (const auto& kv : species_counts) barcodes.push_back(kv.first);
    sort(barcodes.begin(), barcodes.end());

    enum ReassignStatus {
        REASSIGN_ASSIGNED = 0,
        REASSIGN_EMPTY_COUNTS = 1,
        REASSIGN_LLR_FAILED = 2,
        REASSIGN_NONPOSITIVE_MARGIN = 3,
        REASSIGN_NO_WINNER = 4,
        REASSIGN_MISSING_COMPARISONS = 5
    };

    vector<int> statuses(barcodes.size(), REASSIGN_LLR_FAILED);
    vector<int> winners(barcodes.size(), -1);
    vector<double> winner_scores(
        barcodes.size(), std::numeric_limits<double>::quiet_NaN());
    vector<CellDiagnostics> diagnostics(barcodes.size());
    vector<vector<RunnerUp> > runners(barcodes.size());

    const int n_species = (int)species_samples.size();
    set<int> allowed_ids;
    set<int> allowed_ids2;

    fprintf(stderr, "Assigning native species identities to %lu cells using %d threads...\n",
        barcodes.size(), n_threads);

    #pragma omp parallel for num_threads(n_threads) schedule(dynamic, 100)
    for (size_t idx = 0; idx < barcodes.size(); ++idx){
        const unsigned long bc = barcodes[idx];
        const CellCounts& counts = species_counts.at(bc);
        if (counts.is_empty()){
            statuses[idx] = REASSIGN_EMPTY_COUNTS;
            continue;
        }

        map<int, map<int, double> > llrs;
        llr_table tab(n_species);
        const bool success = populate_llr_table_optimized(
            counts, llrs, tab, n_species,
            allowed_ids, allowed_ids2,
            species_doublet_rate, error_ref, error_alt,
            NULL, false, 0.0, 0.0, NULL, n_target);
        if (!success){
            statuses[idx] = REASSIGN_LLR_FAILED;
            continue;
        }

        int winner = -1;
        double score = -std::numeric_limits<double>::infinity();
        tab.get_max(winner, score);
        winners[idx] = winner;
        winner_scores[idx] = score;
        if (winner < 0){
            statuses[idx] = REASSIGN_NO_WINNER;
            continue;
        }

        CellDiagnostics diag;
        vector<RunnerUp> cell_runners;
        get_diagnostics_from_llrs(
            llrs, tab, winner, score, n_species,
            n_runner_ups, close_threshold, diag, cell_runners);
        compute_margin_softmax_scores(llrs, tab, winner, n_species, diag);
        tab.get_max_by_max_llr_comparator(
            diag.max_llr_comparator_winner,
            diag.max_llr_comparator_score);
        diag.total_depth = compute_total_depth(counts, n_species);
        diagnostics[idx] = diag;
        runners[idx] = cell_runners;

        if (diag.selection_resolved){
            statuses[idx] = REASSIGN_ASSIGNED;
        }
        else if (!(score > 0.0) || !std::isfinite(score)){
            statuses[idx] = REASSIGN_NONPOSITIVE_MARGIN;
        }
        else if (!diag.missing_comparison_alternatives.empty()){
            statuses[idx] = REASSIGN_MISSING_COMPARISONS;
        }
        else{
            statuses[idx] = REASSIGN_NO_WINNER;
        }
    }

    auto status_name = [](int status) -> const char* {
        switch(status){
            case REASSIGN_ASSIGNED: return "assigned";
            case REASSIGN_EMPTY_COUNTS: return "empty_counts";
            case REASSIGN_LLR_FAILED: return "llr_construction_failed";
            case REASSIGN_NONPOSITIVE_MARGIN: return "nonpositive_winner_margin";
            case REASSIGN_NO_WINNER: return "no_winner";
            case REASSIGN_MISSING_COMPARISONS: return "missing_required_comparisons";
            default: return "unknown";
        }
    };

    robin_hood::unordered_map<unsigned long, int> assignments;
    robin_hood::unordered_map<unsigned long, double> assignments_llr;
    map<string, long long> status_counts;
    map<string, long long> identity_counts;
    long long n_singlet = 0;
    long long n_pair = 0;

    for (size_t idx = 0; idx < barcodes.size(); ++idx){
        status_counts[status_name(statuses[idx])]++;
        if (statuses[idx] != REASSIGN_ASSIGNED) continue;
        assignments.emplace(barcodes[idx], winners[idx]);
        assignments_llr.emplace(barcodes[idx], winner_scores[idx]);
        const string identity = idx2name(winners[idx], species_samples);
        identity_counts[identity]++;
        if (winners[idx] < n_species) ++n_singlet;
        else ++n_pair;
    }

    OutputTransaction transaction;
    const string staged_prefix = transaction.stage_prefix(output_prefix);
    const string assignment_file_staged = staged_prefix + ".species_assignments";
    const string diag_file_staged = staged_prefix + ".species_assignment_diagnostics.tsv";
    const string summary_file_staged = staged_prefix + ".species_assignment_summary.tsv";

    FILE* assignment_out = fopen(assignment_file_staged.c_str(), "w");
    if (!assignment_out){
        fprintf(stderr, "ERROR: could not open staged species assignments: %s\n", strerror(errno));
        return 1;
    }
    vector<string> tmp_samples = species_samples;
    string tmp_group = barcode_group;
    dump_assignments(assignment_out, assignments, assignments_llr, tmp_samples,
        tmp_group, cellranger, seurat, underscore);
    if (fclose(assignment_out) != 0){
        fprintf(stderr, "ERROR: failed closing staged species assignments: %s\n", strerror(errno));
        return 1;
    }

    FILE* diag_out = fopen(diag_file_staged.c_str(), "w");
    if (!diag_out){
        fprintf(stderr, "ERROR: could not open staged species diagnostics: %s\n", strerror(errno));
        return 1;
    }
    fprintf(diag_out,
        "barcode\tstatus\twinning_identity\tsinglet_doublet\twinner_min_margin\t"
        "llr_vs_runner_up\trunner_up_comparison_state\tworst_competitor\t"
        "worst_comparison_state\tn_close\ttotal_depth\tmargin_softmax_score\t"
        "margin_entropy\tselection_resolved\tmissing_comparison_alternatives\t"
        "runner_up_1\trunner_up_1_llr_vs_winner\trunner_up_1_comparison_state\t"
        "runner_up_1_min_margin\tmax_llr_comparator_winner\t"
        "max_llr_comparator_score\tschema_version\n");

    for (size_t idx = 0; idx < barcodes.size(); ++idx){
        string bc_str = bc2str(barcodes[idx]);
        mod_bc_libname(bc_str, barcode_group, cellranger, seurat, underscore);
        const CellDiagnostics& diag = diagnostics[idx];
        const string winner_name = optional_identity(winners[idx], species_samples);
        const string sd = winners[idx] < 0 ? "NA" :
            (winners[idx] < n_species ? "S" : "D");
        const string worst_name = optional_identity(diag.worst_competitor, species_samples);
        const string missing_names = join_identity_names(
            diag.missing_comparison_alternatives, species_samples);
        const RunnerUp* runner = runners[idx].empty() ? NULL : &runners[idx][0];
        const string runner_name = runner == NULL
            ? "NA" : optional_identity(runner->identity, species_samples);
        const string comparator_name = optional_identity(
            diag.max_llr_comparator_winner, species_samples);

        fprintf(diag_out,
            "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%s\t"
            "%s\t%s\t%s\t%s\t%s\t%s\tdemux_parallel_species_assignment_v3\n",
            bc_str.c_str(), status_name(statuses[idx]), winner_name.c_str(), sd.c_str(),
            optional_number(winner_scores[idx]).c_str(),
            optional_number(diag.llr_vs_runnerup).c_str(),
            comparison_state_name(diag.runnerup_comparison_state),
            worst_name.c_str(), comparison_state_name(diag.worst_comparison_state),
            diag.n_close, optional_number(diag.total_depth).c_str(),
            optional_number(diag.margin_softmax_score).c_str(),
            optional_number(diag.margin_entropy).c_str(),
            diag.selection_resolved ? 1 : 0, missing_names.c_str(), runner_name.c_str(),
            runner == NULL ? "NA" : optional_number(runner->llr_vs_winner).c_str(),
            runner == NULL ? comparison_state_name(ComparisonState::NOT_APPLICABLE)
                           : comparison_state_name(runner->comparison_state),
            runner == NULL ? "NA" : optional_number(runner->min_margin).c_str(),
            comparator_name.c_str(), optional_number(diag.max_llr_comparator_score).c_str());
    }
    if (fclose(diag_out) != 0){
        fprintf(stderr, "ERROR: failed closing staged species diagnostics: %s\n", strerror(errno));
        return 1;
    }

    FILE* summary_out = fopen(summary_file_staged.c_str(), "w");
    if (!summary_out){
        fprintf(stderr, "ERROR: could not open staged species summary: %s\n", strerror(errno));
        return 1;
    }
    fprintf(summary_out, "section\tkey\tvalue\n");
    fprintf(summary_out, "setting\tschema_version\tdemux_parallel_species_assignment_v3\n");
    fprintf(summary_out, "setting\tinput_prefix\t%s\n", input_prefix.c_str());
    fprintf(summary_out, "setting\toutput_prefix\t%s\n", output_prefix.c_str());
    fprintf(summary_out, "setting\tspecies_doublet_rate\t%.17g\n", species_doublet_rate);
    fprintf(summary_out, "setting\tpanel_normalized_legacy_counts\t%s\n",
        panel_normalized ? "true" : "false");
    if (panel_normalized){
        fprintf(summary_out, "setting\tpanel_metadata\t%s\n", panel_metadata_file.c_str());
        for (int sp = 0; sp < (int)species_samples.size(); ++sp){
            fprintf(summary_out, "normalization_singlet\t%s\t%.17g\n",
                species_samples[sp].c_str(), reassign_norm.singlet_weight[sp]);
        }
        for (int a = 0; a < (int)species_samples.size(); ++a){
            for (int b = a + 1; b < (int)species_samples.size(); ++b){
                fprintf(summary_out, "normalization_pair\t%s+%s\t%.17g\n",
                    species_samples[a].c_str(), species_samples[b].c_str(),
                    reassign_norm.pair_weight[a][b]);
            }
        }
    }
    fprintf(summary_out, "setting\terror_ref\t%.17g\n", error_ref);
    fprintf(summary_out, "setting\terror_alt\t%.17g\n", error_alt);
    fprintf(summary_out, "count\tinput_cells\t%lu\n", barcodes.size());
    fprintf(summary_out, "count\tassigned_cells\t%lu\n", assignments.size());
    fprintf(summary_out, "count\tassigned_singlets\t%lld\n", n_singlet);
    fprintf(summary_out, "count\tassigned_pairs\t%lld\n", n_pair);
    for (const auto& kv : status_counts){
        fprintf(summary_out, "status\t%s\t%lld\n", kv.first.c_str(), kv.second);
    }
    for (const auto& kv : identity_counts){
        fprintf(summary_out, "identity\t%s\t%lld\n", kv.first.c_str(), kv.second);
    }
    if (fclose(summary_out) != 0){
        fprintf(stderr, "ERROR: failed closing staged species summary: %s\n", strerror(errno));
        return 1;
    }

    if (!publish_outputs(transaction, "native species reassignment")) return 1;

    fprintf(stderr, "Native species reassignment complete\n");
    fprintf(stderr, "  assigned: %lu / %lu (%lld singlets, %lld pairs)\n",
        assignments.size(), barcodes.size(), n_singlet, n_pair);
    fprintf(stderr, "  assignments: %s.species_assignments\n", output_prefix.c_str());
    fprintf(stderr, "  diagnostics: %s.species_assignment_diagnostics.tsv\n", output_prefix.c_str());
    fprintf(stderr, "  summary: %s.species_assignment_summary.tsv\n", output_prefix.c_str());
    return 0;
}


/**
 * Print help message
 */
void help(int code){
    fprintf(stderr, "demux_parallel version %s\n", VERSION.c_str());
    fprintf(stderr, "%s\n\n", VERSION_MESSAGE.c_str());
    fprintf(stderr, "demux_parallel [OPTIONS]\n");
    fprintf(stderr, "Supported production parallel demultiplexer for genotype-based cell assignment.\n");
    fprintf(stderr, "\n[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --bam -b The BAM file of interest\n");
    fprintf(stderr, "    --vcf -v A VCF/BCF file containing genotypes\n");
    fprintf(stderr, "    --output_prefix -o Base name for output files\n");
    fprintf(stderr, "\n===== OPTIONAL =====\n");
    fprintf(stderr, "    --threads -t Number of threads for parallel processing [auto]\n");
    fprintf(stderr, "    --barcodes -B A file listing cell barcodes (one per line)\n");
    fprintf(stderr, "    --ids -i A file listing allowed individual IDs (singlets)\n");
    fprintf(stderr, "    --ids_doublet -I A file listing allowed doublet combinations\n");
    fprintf(stderr, "    --qual -q Minimum QUAL score for variants [50]\n");
    fprintf(stderr, "    --doublet_rate -D Prior probability of doublet [0.5]\n");
    fprintf(stderr, "    --error_ref -e Prior error rate for ref allele [0.005]\n");
    fprintf(stderr, "    --error_alt -E Prior error rate for alt allele [0.005]\n");
    fprintf(stderr, "    --error_sigma -s Sigma for error rate prior [0.1]\n");
    fprintf(stderr, "    --no_preload -P Disable VCF preloading (use with --shared_vcf for low memory)\n");
    fprintf(stderr, "    --shared_vcf -S Name of shared memory VCF to attach to\n");
    fprintf(stderr, "    --shared_het_vcf NAME Shared memory het VCF for ploidy stats (from vcf_loader_daemon)\n");
    fprintf(stderr, "    --vcf_chroms -c File listing chromosomes to use from VCF\n");
    fprintf(stderr, "    --libname -n Library name to append to barcodes\n");
    fprintf(stderr, "    --cellranger -C Format barcodes for CellRanger\n");
    fprintf(stderr, "    --seurat -R Format barcodes for Seurat\n");
    fprintf(stderr, "    --underscore -U Use underscore instead of hyphen for libname\n");
    fprintf(stderr, "    --disable_conditional -f Disable computing conditional match fractions\n");
    fprintf(stderr, "    --dump_conditional -F Load VCF, compute conditional match fractions, write\n");
    fprintf(stderr, "    --accepted_weighted_conditional Build .condf using accepted BAM-observation site weights (parallel recount only)\n");
    fprintf(stderr, "       .condf file, and exit. No BAM required. Use this to generate .condf after\n");
    fprintf(stderr, "       a run that used -f.\n");
    fprintf(stderr, "    --n_target -N Max singlets to keep before doublet eval [-1=auto, 0=no limit]\n");
    fprintf(stderr, "    --verbose -V Enable verbose output\n");
    fprintf(stderr, "\n===== DIAGNOSTIC OUTPUT (NEW) =====\n");
    fprintf(stderr, "    --debug          Enable DEBUG spam to stderr\n");
    fprintf(stderr, "    --diagnostics    Write diagnostic files (default: ON)\n");
    fprintf(stderr, "    --no-diagnostics Disable diagnostic file output\n");
    fprintf(stderr, "    --het_vcf FILE   Het VCF for ploidy stats (from downsample_vcf_parallel)\n");
    fprintf(stderr, "    --het_method M   Het balance method: welford (default) or persite\n");
    fprintf(stderr, "                     welford: Online variance, low memory\n");
    fprintf(stderr, "                     persite: Store per-site counts, more memory\n");
    fprintf(stderr, "    --min_het_sites N  Minimum het sites for variance (default: 100)\n");
    fprintf(stderr, "    --min_het_depth D  Minimum depth per het site for persite method (default: 5.0)\n");
    fprintf(stderr, "    --n_runner_ups N Number of runner-ups to report [8]\n");
    fprintf(stderr, "    --close_threshold F LLR threshold for n_close [20.0]\n");
    fprintf(stderr, "\n===== SKIP ASSIGNMENT =====\n");
    fprintf(stderr, "    --skip_assignment -K Write .counts and exit (no assignment)\n");
    fprintf(stderr, "\n===== PILEUP (variant-consistency benchmark) =====\n");
    fprintf(stderr, "    --dump_pileup PREFIX Emit PREFIX.pileup_sites.tsv.gz and PREFIX.pileup_obs.tsv.gz\n");
    fprintf(stderr, "                         (requires an actual counting pass: --force_recount and --threads > 1;\n");
    fprintf(stderr, "                          supported on the single-panel and dual-panel parallel paths)\n");
    fprintf(stderr, "    --dump_selection_audit  Append maximin-vs-max_llr_comparator columns to .diagnostics.gz\n");
    fprintf(stderr, "    --dump_source_observations PREFIX  Emit PREFIX.source_observations.tsv.gz from accepted YI-tagged read-SNP observations\n");
    fprintf(stderr, "    --source_provenance_tag TAG  Two-character BAM tag carrying injected source identity [YI]\n");
    fprintf(stderr, "    --synthetic_id_tag TAG  Two-character BAM tag carrying synthetic unit ID [YS]\n");
    fprintf(stderr, "    --expected_synthetic_id ID  Required unit ID; non-native YI reads must carry matching YS; YI=__NATIVE__ is retained and exempt\n");
    fprintf(stderr, "    --source_receiver_map FILE  Optional two-column barcode/receiver_identity map; adds receiver and (nalt_A,nalt_B) to individual source rows\n");
    fprintf(stderr, "    --source_reconciliation_mode  Account every accepted individual-panel observation into explicit provenance buckets and emit PREFIX.source_reconciliation.tsv.gz\n");
    fprintf(stderr, "    --source_donor_site_audit  Audit every accepted injected observation against its YI donor genotype; emit PREFIX.donor_genotype_audit.tsv.gz and a deterministic exact-site sample\n");
    fprintf(stderr, "    --source_donor_site_sample_mod N  Retain exact-site evidence when stable_hash(site) %% N == 0 [256; 1=all sites]\n");
    fprintf(stderr, "                            (maximin_winner, maximin_score, max_llr_comparator_winner, max_llr_comparator_score, selection_agree).\n");
    fprintf(stderr, "                            The comparator is non-mutating argmax(maxllr), not the historical destructive selector.\n");
    fprintf(stderr, "                         for the interindividual panel. Use with -B (candidate\n");
    fprintf(stderr, "                         barcodes) and -K. Routes through the single-panel path,\n");
    fprintf(stderr, "                         so do not pass a species VCF on the pileup run.\n");
    fprintf(stderr, "\n===== COUNTS FILE SAFETY =====\n");
    fprintf(stderr, "    --reuse_counts       Load existing .counts file after validating integrity\n");
    fprintf(stderr, "    --force_recount      Recount to staged files; replace existing finals only after success\n");
    fprintf(stderr, "    (default: error if .counts exists, to prevent stale/truncated reuse)\n");
    fprintf(stderr, "\n===== ATAC DUAL-MODALITY (2A) =====\n");
    fprintf(stderr, "    --atac_bam FILE      ATAC BAM for diagnostic/count-only evidence; does not alter individual assignment\n");
    fprintf(stderr, "    --atac_vcf FILE      ATAC-demux VCF panel (required if --atac_bam set)\n");
    fprintf(stderr, "    --atac_het_vcf FILE  ATAC-side het VCF for ATAC het balance diagnostics\n");
    fprintf(stderr, "    --atac_shared_vcf NAME   Shared-memory ATAC VCF (mutually exclusive with --atac_vcf)\n");
    fprintf(stderr, "    --atac_shared_het_vcf NAME   Shared-memory ATAC het VCF\n");
    fprintf(stderr, "\n===== IDENTITY PRIOR (2B) =====\n");
    fprintf(stderr, "    --identity_prior FILE  NOT IMPLEMENTED / NOT USED IN IDENTITY SCORING; passing it is a fatal error\n");
    fprintf(stderr, "\n===== SPECIES PANEL (2C) =====\n");
    fprintf(stderr, "    --species_vcf FILE         Species-discrimination VCF panel\n");
    fprintf(stderr, "    --species_shared_vcf NAME  Shared-memory species VCF\n");
    fprintf(stderr, "    --species_panel_mode MODE  count_only (native V1_R3 mode; required if --species_vcf set)\n");
    fprintf(stderr, "       count_only: count/write native species panel artifacts without using species SNPs in individual assignment\n");
    fprintf(stderr, "       filter|augment|both are legacy mixed modes and require --allow_legacy_mixed_species_panel\n");
    fprintf(stderr, "    --allow_legacy_mixed_species_panel  Permit old species-panel filter/augment/both behavior (not V1_R3 native)\n");
    fprintf(stderr, "    --species_assignment_output Produce .species_assignments (species-only scoring pass)\n");
    fprintf(stderr, "    --species_counts_output     Write .species_counts, .species_condf alongside .counts\n");
    fprintf(stderr, "    --panel_metadata FILE  Tab-separated indiv_id/species mapping (required for species output;\n");
    fprintf(stderr, "                           in --species_reassign_from_prefix mode, normalize legacy counts in memory)\n");
    fprintf(stderr, "    --species_ids FILE  A file listing allowed species labels (singlets) and\n");
    fprintf(stderr, "                        species pairs (A+B). Species analogue of -i/-I; without it\n");
    fprintf(stderr, "                        the species pass considers every label and every pair.\n");
    fprintf(stderr, "    --species_doublet_rate F  Species-assignment pair prior; defaults to --doublet_rate\n");
    fprintf(stderr, "    --species_reassign_from_prefix PREFIX\n");
    fprintf(stderr, "                         Re-run species assignment from PREFIX.species_counts and\n");
    fprintf(stderr, "                         PREFIX.species_samples only; no BAM/VCF access or recount.\n");
    fprintf(stderr, "                         -o must be a distinct diagnostic output prefix.\n");
    fprintf(stderr, "\n    --help -h Display this message and exit\n");
    exit(code);
}

// Timing helper
void print_elapsed(const std::chrono::steady_clock::time_point& start, const char* step) {
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
    int hours = elapsed / 3600;
    int mins = (elapsed % 3600) / 60;
    int secs = elapsed % 60;
    fprintf(stderr, "[%02d:%02d:%02d] %s\n", hours, mins, secs, step);
}

// ============================================================================
// CORRECTED DIAGNOSTIC OUTPUT FUNCTIONS
// ============================================================================

/** Write per-cell diagnostics with explicit comparison/missing-value semantics. */
bool write_diagnostics_gz(
    const string& filename,
    vector<string>& samples,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, double>& assignments_llr,
    robin_hood::unordered_map<unsigned long, CellDiagnostics>& diagnostics,
    int n_samples,
    const string& libname,
    bool cellranger,
    bool seurat,
    bool underscore,
    robin_hood::unordered_map<unsigned long, CellHetData>* atac_het_data = NULL,
    robin_hood::unordered_map<int, ChromSNPs>* atac_het_snpdat_ptr = NULL,
    vector<pair<int, int>>* atac_idx_to_site_ptr = NULL,
    HetBalanceMethod atac_het_method = HetBalanceMethod::WELFORD,
    int atac_min_het_sites = 100,
    double atac_min_het_depth = 5.0,
    bool dump_selection_audit = false){

    (void)assignments_llr;
    gzFile outf = gzopen(filename.c_str(), "w");
    if (!outf){
        fprintf(stderr, "ERROR: could not open %s for writing: %s\n",
            filename.c_str(), strerror(errno));
        return false;
    }

    const bool write_atac_het = atac_het_data != NULL;
    gzprintf(outf,
        "barcode\tassignment\tsinglet_doublet\tllr_vs_runner_up\tmin_margin\tworst_competitor\t"
        "n_close\ttotal_depth\thet_balance_var\tn_het_sites\thet_total_depth\t"
        "margin_softmax_score\tmargin_entropy");
    if (write_atac_het){
        gzprintf(outf,
            "\tatac_het_balance_var\tatac_n_het_sites\tatac_het_total_depth\t"
            "atac_het_diagnostic_state");
    }
    gzprintf(outf,
        "\trunnerup_comparison_state\tworst_comparison_state\thet_diagnostic_state\t"
        "selection_resolved\tcandidate_identity\tmissing_comparison_alternatives\t"
        "schema_version");
    if (dump_selection_audit){
        gzprintf(outf,
            "\tmaximin_winner\tmaximin_score\tmax_llr_comparator_winner\t"
            "max_llr_comparator_score\tselection_agree");
    }
    gzprintf(outf, "\n");

    vector<unsigned long> barcodes;
    barcodes.reserve(diagnostics.size());
    for (const auto& kv : diagnostics) barcodes.push_back(kv.first);
    sort(barcodes.begin(), barcodes.end());

    for (unsigned long bc : barcodes){
        const CellDiagnostics& diag = diagnostics.at(bc);
        const auto assignment_it = assignments.find(bc);
        const bool assigned = assignment_it != assignments.end();
        const int assignment = assigned ? assignment_it->second : -1;

        string bc_str = bc2str(bc);
        mod_bc_libname(bc_str, libname, cellranger, seurat, underscore);
        const string assignment_name = assigned
            ? idx2name(assignment, samples) : "unresolved";
        const string sd = assigned ? (assignment < n_samples ? "S" : "D") : "U";
        const string worst_name = optional_identity(diag.worst_competitor, samples);
        const string candidate_name = optional_identity(diag.maximin_candidate, samples);
        const string missing_names = join_identity_names(
            diag.missing_comparison_alternatives, samples);

        const string het_var = diag.het_diagnostic_available
            ? optional_number(diag.het_balance_var) : "NA";
        const string het_sites = diag.het_diagnostic_available
            ? std::to_string(diag.n_het_sites) : "NA";
        const string het_depth = diag.het_diagnostic_available
            ? optional_number(diag.het_total_depth) : "NA";

        gzprintf(outf,
            "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s",
            bc_str.c_str(), assignment_name.c_str(), sd.c_str(),
            optional_number(diag.llr_vs_runnerup).c_str(),
            optional_number(diag.min_margin).c_str(), worst_name.c_str(),
            diag.n_close, optional_number(diag.total_depth).c_str(),
            het_var.c_str(), het_sites.c_str(), het_depth.c_str(),
            optional_number(diag.margin_softmax_score).c_str(),
            optional_number(diag.margin_entropy).c_str());

        if (write_atac_het){
            CellDiagnostics atac_diag;
            if (assigned){
                const auto atac_it = atac_het_data->find(bc);
                if (atac_it != atac_het_data->end()){
                    const CellHetData& cell_atac_het = atac_it->second;
                    if (atac_het_method == HetBalanceMethod::PERSITE &&
                        atac_het_snpdat_ptr != NULL && atac_idx_to_site_ptr != NULL){
                        compute_het_balance_persite(cell_atac_het.persite_data,
                            *atac_het_snpdat_ptr, *atac_idx_to_site_ptr,
                            assignment, n_samples, atac_min_het_depth,
                            atac_min_het_sites, atac_diag);
                    }
                    else{
                        compute_het_balance_welford(cell_atac_het.welford_stats,
                            assignment, n_samples, atac_min_het_sites, atac_diag);
                    }
                }
            }
            const string atac_var = atac_diag.het_diagnostic_available
                ? optional_number(atac_diag.het_balance_var) : "NA";
            const string atac_sites = atac_diag.het_diagnostic_available
                ? std::to_string(atac_diag.n_het_sites) : "NA";
            const string atac_depth = atac_diag.het_diagnostic_available
                ? optional_number(atac_diag.het_total_depth) : "NA";
            gzprintf(outf, "\t%s\t%s\t%s\t%s",
                atac_var.c_str(), atac_sites.c_str(), atac_depth.c_str(),
                atac_diag.het_diagnostic_available ? "available" : "unavailable");
        }

        gzprintf(outf, "\t%s\t%s\t%s\t%d\t%s\t%s\t%s",
            comparison_state_name(diag.runnerup_comparison_state),
            comparison_state_name(diag.worst_comparison_state),
            diag.het_diagnostic_available ? "available" : "unavailable",
            diag.selection_resolved ? 1 : 0,
            candidate_name.c_str(), missing_names.c_str(),
            "demux_parallel_diagnostics_v3");

        if (dump_selection_audit){
            const string comparator_name = optional_identity(
                diag.max_llr_comparator_winner, samples);
            string agree = "NA";
            if (diag.maximin_candidate >= 0 && diag.max_llr_comparator_winner >= 0){
                agree = diag.maximin_candidate == diag.max_llr_comparator_winner ? "1" : "0";
            }
            gzprintf(outf, "\t%s\t%s\t%s\t%s\t%s",
                candidate_name.c_str(), optional_number(diag.maximin_score).c_str(),
                comparator_name.c_str(),
                optional_number(diag.max_llr_comparator_score).c_str(),
                agree.c_str());
        }
        gzprintf(outf, "\n");
    }

    if (gzclose(outf) != Z_OK){
        fprintf(stderr, "ERROR: failed closing %s\n", filename.c_str());
        return false;
    }
    fprintf(stderr, "Wrote diagnostics for %lu evaluated cells to %s\n",
        diagnostics.size(), filename.c_str());
    return true;
}

/** Write ranked retained alternatives with authoritative direct comparisons. */
bool write_runner_ups_gz(
    const string& filename,
    vector<string>& samples,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, CellDiagnostics>& diagnostics,
    robin_hood::unordered_map<unsigned long, vector<RunnerUp> >& runner_ups,
    int n_samples,
    const string& libname,
    bool cellranger,
    bool seurat,
    bool underscore){

    (void)assignments;
    (void)n_samples;
    gzFile outf = gzopen(filename.c_str(), "w");
    if (!outf){
        fprintf(stderr, "ERROR: could not open %s for writing: %s\n",
            filename.c_str(), strerror(errno));
        return false;
    }

    gzprintf(outf,
        "barcode\trank\tidentity\tllr_vs_winner\tmin_margin\tcomparison_state\t"
        "winner_candidate\tselection_resolved\tschema_version\n");

    vector<unsigned long> barcodes;
    barcodes.reserve(runner_ups.size());
    for (const auto& kv : runner_ups) barcodes.push_back(kv.first);
    sort(barcodes.begin(), barcodes.end());

    size_t rows_written = 0;
    for (unsigned long bc : barcodes){
        const vector<RunnerUp>& runners = runner_ups.at(bc);
        const auto diag_it = diagnostics.find(bc);
        const CellDiagnostics* diag = diag_it == diagnostics.end() ? NULL : &diag_it->second;
        string bc_str = bc2str(bc);
        mod_bc_libname(bc_str, libname, cellranger, seurat, underscore);
        const string candidate_name = diag == NULL
            ? "NA" : optional_identity(diag->maximin_candidate, samples);
        const int resolved = diag != NULL && diag->selection_resolved ? 1 : 0;

        for (size_t i = 0; i < runners.size(); ++i){
            const RunnerUp& runner = runners[i];
            const string identity = optional_identity(runner.identity, samples);
            gzprintf(outf, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%s\n",
                bc_str.c_str(), (int)(i + 1), identity.c_str(),
                optional_number(runner.llr_vs_winner).c_str(),
                optional_number(runner.min_margin).c_str(),
                comparison_state_name(runner.comparison_state),
                candidate_name.c_str(), resolved,
                "demux_parallel_runner_ups_v3");
            ++rows_written;
        }
    }

    if (gzclose(outf) != Z_OK){
        fprintf(stderr, "ERROR: failed closing %s\n", filename.c_str());
        return false;
    }
    fprintf(stderr, "Wrote %lu runner-up rows for %lu evaluated cells to %s\n",
        rows_written, runner_ups.size(), filename.c_str());
    return true;
}


int main(int argc, char *argv[]) {    
    
    // Start timing
    auto start_time = std::chrono::steady_clock::now();
    
    // Print version info
    fprintf(stderr, "demux_parallel version %s\n", VERSION.c_str());
    fprintf(stderr, "%s\n", VERSION_MESSAGE.c_str());
    fprintf(stderr, "New: %s\n", VERSION_NEW.c_str());
#ifdef _OPENMP
    fprintf(stderr, "Build: compiler=%s; OpenMP=%d; HTSlib=%s; source_revision=%s\n",
        __VERSION__, _OPENMP, hts_version(), CELLBOUNCER_SOURCE_REVISION);
#else
    fprintf(stderr, "Build: compiler=%s; OpenMP=disabled; HTSlib=%s; source_revision=%s\n",
        __VERSION__, hts_version(), CELLBOUNCER_SOURCE_REVISION);
#endif
   
    static struct option long_options[] = {
       {"bam", required_argument, 0, 'b'},
       {"vcf", required_argument, 0, 'v'},
       {"output_prefix", required_argument, 0, 'o'},
       {"barcodes", required_argument, 0, 'B'},
       {"doublet_rate", required_argument, 0, 'D'},
       {"ids", required_argument, 0, 'i'},
       {"ids_doublet", required_argument, 0, 'I'},
       {"qual", required_argument, 0, 'q'},
       {"libname", required_argument, 0, 'n'},
       {"cellranger", no_argument, 0, 'C'},
       {"seurat", no_argument, 0, 'R'},
       {"underscore", no_argument, 0, 'U'},
       {"error_ref", required_argument, 0, 'e'},
       {"error_alt", required_argument, 0, 'E'},
       {"error_sigma", required_argument, 0, 's'},
       {"disable_conditional", no_argument, 0, 'f'},
       {"dump_conditional", no_argument, 0, 'F'},
       {"no_preload", no_argument, 0, 'P'},
       {"vcf_chroms", required_argument, 0, 'c'},
       {"threads", required_argument, 0, 't'},
       {"shared_vcf", required_argument, 0, 'S'},
       {"shared_het_vcf", required_argument, 0, 1007},
       {"n_target", required_argument, 0, 'N'},
       {"verbose", no_argument, 0, 'V'},
       // NEW diagnostic options
       {"debug", no_argument, 0, 1001},
       {"diagnostics", no_argument, 0, 1002},
       {"no-diagnostics", no_argument, 0, 1003},
       {"het_vcf", required_argument, 0, 1004},
       {"n_runner_ups", required_argument, 0, 1005},
       {"close_threshold", required_argument, 0, 1006},
       {"het_method", required_argument, 0, 1010},
       {"min_het_sites", required_argument, 0, 1011},
       {"min_het_depth", required_argument, 0, 1012},
       // Step 0a: skip assignment
       {"skip_assignment", no_argument, 0, 'K'},
       // Variant-consistency benchmark: per-SNP pileup sidecars
       {"dump_pileup", required_argument, 0, 1099},
       {"dump_selection_audit", no_argument, 0, 1100},
       {"dump_source_observations", required_argument, 0, 1101},
       {"source_provenance_tag", required_argument, 0, 1102},
       {"synthetic_id_tag", required_argument, 0, 1103},
       {"expected_synthetic_id", required_argument, 0, 1104},
       {"source_receiver_map", required_argument, 0, 1105},
       {"source_reconciliation_mode", no_argument, 0, 1106},
       {"source_donor_site_audit", no_argument, 0, 1107},
       {"source_donor_site_sample_mod", required_argument, 0, 1108},
       {"accepted_weighted_conditional", no_argument, 0, 1109},
       // 2A: ATAC dual-modality
       {"atac_bam", required_argument, 0, 1020},
       {"atac_vcf", required_argument, 0, 1021},
       {"atac_het_vcf", required_argument, 0, 1022},
       {"atac_shared_vcf", required_argument, 0, 1023},
       {"atac_shared_het_vcf", required_argument, 0, 1024},
       // 2B: identity prior
       {"identity_prior", required_argument, 0, 1030},
       // 2C: species panel
       {"species_vcf", required_argument, 0, 1040},
       {"species_shared_vcf", required_argument, 0, 1041},
       {"species_panel_mode", required_argument, 0, 1042},
       {"allow_legacy_mixed_species_panel", no_argument, 0, 1059},
       {"species_assignment_output", no_argument, 0, 1043},
       {"species_counts_output", no_argument, 0, 1045},
       {"panel_metadata", required_argument, 0, 1044},
       {"species_doublet_rate", required_argument, 0, 1060},
       {"species_ids", required_argument, 0, 1062},
       {"species_reassign_from_prefix", required_argument, 0, 1061},
       {"reuse_counts", no_argument, 0, 1050},
       {"force_recount", no_argument, 0, 1051},
       {"version", no_argument, 0, 1063},
       {"help", no_argument, 0, 'h'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string bamfile = "";
    string vcf_file = "";
    bool cell_barcode = false;
    string cell_barcode_file = "";
    string output_prefix = "";
    int vq = 50;
    string idfile;
    string idfile_doublet;
    bool idfile_given = false;
    bool idfile_doublet_given = false;
    double doublet_rate = 0.5;
    string barcode_group = "";
    double error_ref = 0.005;
    double error_alt = 0.005;
    double error_sigma = 0.1;
    
    bool cellranger = false;
    bool seurat = false;
    bool underscore = false;

    bool disable_conditional = false;
    bool dump_conditional = false;

    bool no_preload = false;
    string vcf_chroms_file = "";
    bool vcf_chroms_given = false;
    
    // New parallel processing options
    int n_threads = omp_get_num_procs();
    int htslib_threads = 0;  // 0 = auto-calculate after parsing args
    string shared_vcf_name = "";
    string shared_het_vcf_name = "";
    int n_target = -1;    // -1=auto, 0=no prune, >0=use value
    bool verbose = false;
    
    // NEW diagnostic options
    bool debug_mode = false;
    bool write_diagnostics = true;  // Default ON
    string het_vcf_file = "";
    HetBalanceMethod het_method = HetBalanceMethod::WELFORD;  // Default to Welford
    int min_het_sites = 100;   // Minimum sites for reliable variance
    double min_het_depth = 5.0;  // Minimum depth per site (persite method only)
    int n_runner_ups = 8;
    double close_threshold = 20.0;

    // Step 0a: skip assignment
    bool skip_assignment = false;
    bool dump_selection_audit = false;

    // Variant-consistency benchmark: emit per-SNP pileup sidecars
    bool dump_pileup = false;
    string dump_pileup_prefix = "";
    bool dump_source_observations = false;
    string dump_source_observations_prefix = "";
    string source_provenance_tag = "YI";
    string synthetic_id_tag = "YS";
    string expected_synthetic_id = "";
    string source_receiver_map = "";
    bool source_reconciliation_mode = false;
    bool source_donor_site_audit = false;
    int source_donor_site_sample_mod = 256;
    bool accepted_weighted_conditional = false;

    // 2A: ATAC dual-modality
    string atac_bamfile = "";
    string atac_vcf_file = "";
    string atac_het_vcf_file = "";
    string atac_shared_vcf_name = "";
    string atac_shared_het_vcf_name = "";

    // 2B: identity prior
    string identity_prior_file = "";

    // 2C: species panel
    string species_vcf_file = "";
    string species_shared_vcf_name = "";
    string species_panel_mode_str = "";
    bool allow_legacy_mixed_species_panel = false;
    bool species_assignment_output = false;
    bool species_counts_output = false;
    string panel_metadata_file = "";
    double species_doublet_rate = -1.0;  // <0 means inherit --doublet_rate
    string species_idfile;               // empty means unrestricted (legacy)
    string species_reassign_from_prefix = "";

    // Counts file reuse safety
    bool reuse_counts = false;
    bool force_recount = false;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:v:o:B:i:I:q:D:n:e:E:s:c:t:S:N:fFPCRUVhK", 
        long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                break;
            case 'h':
                help(0);
                break;
            case 1063:  // --version
                return 0;
            case 'b':
                bamfile = optarg;
                break;
            case 'v':
                vcf_file = optarg;
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 'n':
                barcode_group = optarg;
                break;
            case 'C':
                cellranger = true;
                break;
            case 'R':
                seurat = true;
                break;
            case 'U':
                underscore = true;
                break;
            case 'B':
                cell_barcode = true;
                cell_barcode_file = optarg;
                break;
            case 'D':
                doublet_rate = atof(optarg);
                break;
            case 'i':
                idfile_given = true;
                idfile = optarg;
                break;
            case 'I':
                idfile_doublet_given = true;
                idfile_doublet = optarg;
                break;
            case 'q':
                vq = atoi(optarg);
                break;
            case 'e':
                error_ref = atof(optarg);
                break;
            case 'E':
                error_alt = atof(optarg);
                break;
            case 's':
                error_sigma = atof(optarg);
                break;
            case 'f':
                disable_conditional = true;
                break;
            case 'F':
                dump_conditional = true;
                break;
            case 'P':
                no_preload = true;
                break;
            case 'c':
                vcf_chroms_given = true;
                vcf_chroms_file = optarg;
                break;
            case 't':
                n_threads = atoi(optarg);
                break;
            case 'S':
                shared_vcf_name = optarg;
                break;
            case 'N':
                n_target = atoi(optarg);
                break;
            case 'V':
                verbose = true;
                break;
            // NEW diagnostic options
            case 1001:  // --debug
                debug_mode = true;
                break;
            case 1002:  // --diagnostics
                write_diagnostics = true;
                break;
            case 1003:  // --no-diagnostics
                write_diagnostics = false;
                break;
            case 1004:  // --het_vcf
                het_vcf_file = optarg;
                break;
            case 1005:  // --n_runner_ups
                n_runner_ups = atoi(optarg);
                break;
            case 1006:  // --close_threshold
                close_threshold = atof(optarg);
                break;
            case 1007:  // --shared_het_vcf
                shared_het_vcf_name = optarg;
                break;
            case 1010:  // --het_method
                if (strcmp(optarg, "persite") == 0) {
                    het_method = HetBalanceMethod::PERSITE;
                } else if (strcmp(optarg, "welford") == 0) {
                    het_method = HetBalanceMethod::WELFORD;
                } else {
                    fprintf(stderr, "ERROR: Unknown het_method: %s\n", optarg);
                    fprintf(stderr, "Valid options: welford, persite\n");
                    return 1;
                }
                break;
            case 1011:  // --min_het_sites
                min_het_sites = atoi(optarg);
                break;
            case 1012:  // --min_het_depth
                min_het_depth = atof(optarg);
                break;
            case 'K':  // --skip_assignment
                skip_assignment = true;
                break;
            case 1099:  // --dump_pileup
                dump_pileup = true;
                dump_pileup_prefix = optarg;
                break;
            case 1100:  // --dump_selection_audit
                dump_selection_audit = true;
                break;
            case 1101:  // --dump_source_observations
                dump_source_observations = true;
                dump_source_observations_prefix = optarg;
                break;
            case 1102:  // --source_provenance_tag
                source_provenance_tag = optarg;
                break;
            case 1103:  // --synthetic_id_tag
                synthetic_id_tag = optarg;
                break;
            case 1104:  // --expected_synthetic_id
                expected_synthetic_id = optarg;
                break;
            case 1105:  // --source_receiver_map
                source_receiver_map = optarg;
                break;
            case 1106:  // --source_reconciliation_mode
                source_reconciliation_mode = true;
                break;
            case 1107:  // --source_donor_site_audit
                source_donor_site_audit = true;
                break;
            case 1108:  // --source_donor_site_sample_mod
                source_donor_site_sample_mod = atoi(optarg);
                break;
            case 1109:  // --accepted_weighted_conditional
                accepted_weighted_conditional = true;
                break;
            case 1020:  // --atac_bam
                atac_bamfile = optarg;
                break;
            case 1021:  // --atac_vcf
                atac_vcf_file = optarg;
                break;
            case 1022:  // --atac_het_vcf
                atac_het_vcf_file = optarg;
                break;
            case 1023:  // --atac_shared_vcf
                atac_shared_vcf_name = optarg;
                break;
            case 1024:  // --atac_shared_het_vcf
                atac_shared_het_vcf_name = optarg;
                break;
            case 1030:  // --identity_prior
                identity_prior_file = optarg;
                break;
            case 1040:  // --species_vcf
                species_vcf_file = optarg;
                break;
            case 1041:  // --species_shared_vcf
                species_shared_vcf_name = optarg;
                break;
            case 1042:  // --species_panel_mode
                species_panel_mode_str = optarg;
                break;
            case 1059:  // --allow_legacy_mixed_species_panel
                allow_legacy_mixed_species_panel = true;
                break;
            case 1043:  // --species_assignment_output
                species_assignment_output = true;
                break;
            case 1045:  // --species_counts_output
                species_counts_output = true;
                break;
            case 1044:  // --panel_metadata
                panel_metadata_file = optarg;
                break;
            case 1062:  // --species_ids
                species_idfile = optarg;
                break;
            case 1060:  // --species_doublet_rate
                species_doublet_rate = atof(optarg);
                break;
            case 1061:  // --species_reassign_from_prefix
                species_reassign_from_prefix = optarg;
                break;
            case 1050:  // --reuse_counts
                reuse_counts = true;
                break;
            case 1051:  // --force_recount
                force_recount = true;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Set global verbose flag
    g_verbose = verbose;
    
    // Set global debug flag
    g_debug = debug_mode;
    
    // Error check arguments
    if (reuse_counts && force_recount){
        fprintf(stderr, "ERROR: --reuse_counts and --force_recount are mutually exclusive\n");
        return 1;
    }
    if (!identity_prior_file.empty()){
        fprintf(stderr,
            "ERROR: --identity_prior is NOT IMPLEMENTED / NOT USED IN IDENTITY SCORING. "
            "Refusing to continue rather than silently ignoring the file.\n");
        return 1;
    }
    if (vq < 0){
        fprintf(stderr, "ERROR: variant quality must be >= 0\n");
        return 1;
    }
    if (output_prefix.length() == 0){
        fprintf(stderr, "ERROR: output_prefix/-o required\n");
        return 1;
    }
    if (is_dir(output_prefix)){
        fprintf(stderr, "ERROR: output_prefix %s is a directory\n", output_prefix.c_str());
        return 1;
    }
    if (doublet_rate < 0 || doublet_rate > 1){
        fprintf(stderr, "ERROR: doublet rate must be between 0 and 1\n");
        return 1;
    }
    if (idfile_given && idfile_doublet_given){
        fprintf(stderr, "ERROR: only one of -i/-I is allowed\n");
        return 1;
    }
    if (n_threads < 1){
        n_threads = 1;
    }
    if (species_doublet_rate < 0.0){
        species_doublet_rate = doublet_rate;
    }
    if (species_doublet_rate < 0.0 || species_doublet_rate > 1.0){
        fprintf(stderr, "ERROR: species doublet rate must be between 0 and 1\n");
        return 1;
    }

    // Counts-only native species reassignment is an isolated diagnostic mode.
    // Return before BAM/VCF/species-panel validation so it cannot trigger a recount.
    if (!species_reassign_from_prefix.empty()){
        return reassign_native_species_from_existing(
            species_reassign_from_prefix, output_prefix, species_doublet_rate,
            error_ref, error_alt, n_threads, n_target, n_runner_ups, close_threshold,
            barcode_group, cellranger, seurat, underscore, panel_metadata_file);
    }

    if (source_provenance_tag.size() != 2){
        fprintf(stderr, "ERROR: --source_provenance_tag must be exactly two characters (got '%s')\n", source_provenance_tag.c_str());
        return 1;
    }
    if (synthetic_id_tag.size() != 2){
        fprintf(stderr, "ERROR: --synthetic_id_tag must be exactly two characters (got '%s')\n", synthetic_id_tag.c_str());
        return 1;
    }
    if (dump_source_observations && dump_source_observations_prefix.empty()){
        fprintf(stderr, "ERROR: --dump_source_observations requires a non-empty prefix\n");
        return 1;
    }
    if (dump_source_observations && expected_synthetic_id.empty()){
        fprintf(stderr, "ERROR: --dump_source_observations requires --expected_synthetic_id so YS provenance can be validated\n");
        return 1;
    }
    if (!source_receiver_map.empty() && !dump_source_observations){
        fprintf(stderr, "ERROR: --source_receiver_map is meaningful only with --dump_source_observations\n");
        return 1;
    }
    if (source_reconciliation_mode && !dump_source_observations){
        fprintf(stderr, "ERROR: --source_reconciliation_mode requires --dump_source_observations\n");
        return 1;
    }
    if (source_reconciliation_mode && source_receiver_map.empty()){
        fprintf(stderr, "ERROR: --source_reconciliation_mode requires --source_receiver_map\n");
        return 1;
    }
    if (source_donor_site_audit && !dump_source_observations){
        fprintf(stderr, "ERROR: --source_donor_site_audit requires --dump_source_observations\n");
        return 1;
    }
    if (source_donor_site_audit && source_receiver_map.empty()){
        fprintf(stderr, "ERROR: --source_donor_site_audit requires --source_receiver_map\n");
        return 1;
    }
    if (source_donor_site_sample_mod < 1){
        fprintf(stderr, "ERROR: --source_donor_site_sample_mod must be >= 1\n");
        return 1;
    }
    if (!source_receiver_map.empty()){
        struct stat st;
        if (stat(source_receiver_map.c_str(), &st) != 0 || st.st_size <= 0){
            fprintf(stderr, "ERROR: --source_receiver_map missing or empty: %s\n", source_receiver_map.c_str());
            return 1;
        }
    }
    
    // Auto-calculate htslib threads per reader
    // With many workers (>=8), use 1 thread to avoid oversubscription
    // With few workers (<8), use 2 threads since decompression may bottleneck
    if (htslib_threads == 0){
        htslib_threads = (n_threads >= 8) ? 1 : 2;
    }
    
    fprintf(stderr, "demux_parallel: Parallel VCF-based demultiplexing\n");
    fprintf(stderr, "Using %d threads (%d htslib threads per reader)\n", n_threads, htslib_threads);

    // Validate new flags
    bool atac_mode = (atac_bamfile.length() > 0);
    if (atac_mode){
        if (atac_vcf_file.length() == 0 && atac_shared_vcf_name.length() == 0){
            fprintf(stderr, "ERROR: --atac_vcf or --atac_shared_vcf required when --atac_bam is set\n");
            return 1;
        }
        if (atac_vcf_file.length() > 0 && atac_shared_vcf_name.length() > 0){
            fprintf(stderr, "ERROR: --atac_vcf and --atac_shared_vcf are mutually exclusive\n");
            return 1;
        }
    }
    if (atac_het_vcf_file.length() > 0 && !atac_mode){
        fprintf(stderr, "WARNING: --atac_het_vcf ignored without --atac_bam\n");
    }

    // Species panel mode validation
    SpeciesPanelMode species_mode = SpeciesPanelMode::NONE;
    bool has_species_vcf = (species_vcf_file.length() > 0 || species_shared_vcf_name.length() > 0);
    if (has_species_vcf){
        if (species_panel_mode_str.length() == 0){
            fprintf(stderr, "ERROR: --species_panel_mode required when --species_vcf is set\n");
            return 1;
        }
        if (species_panel_mode_str == "count_only" || species_panel_mode_str == "native" || species_panel_mode_str == "none") species_mode = SpeciesPanelMode::COUNT_ONLY;
        else if (species_panel_mode_str == "filter" || species_panel_mode_str == "augment" || species_panel_mode_str == "both"){
            if (!allow_legacy_mixed_species_panel){
                fprintf(stderr, "ERROR: --species_panel_mode %s is a legacy mixed species-panel mode.\n",
                    species_panel_mode_str.c_str());
                fprintf(stderr, "       Under V1_R3 path separation, use --species_panel_mode count_only.\n");
                fprintf(stderr, "       Use --allow_legacy_mixed_species_panel only for intentional legacy debugging.\n");
                return 1;
            }
            if (species_panel_mode_str == "filter") species_mode = SpeciesPanelMode::FILTER;
            else if (species_panel_mode_str == "augment") species_mode = SpeciesPanelMode::AUGMENT;
            else species_mode = SpeciesPanelMode::BOTH;
        }
        else{
            fprintf(stderr, "ERROR: --species_panel_mode must be count_only (or legacy filter/augment/both with --allow_legacy_mixed_species_panel)\n");
            return 1;
        }
        if (species_vcf_file.length() > 0 && species_shared_vcf_name.length() > 0){
            fprintf(stderr, "ERROR: --species_vcf and --species_shared_vcf are mutually exclusive\n");
            return 1;
        }
    }
    if (species_assignment_output && !has_species_vcf){
        fprintf(stderr, "ERROR: --species_assignment_output requires --species_vcf\n");
        return 1;
    }
    if ((species_counts_output || species_assignment_output) && panel_metadata_file.length() == 0){
        fprintf(stderr, "ERROR: --species_counts_output/--species_assignment_output requires --panel_metadata so .species_* artifacts are species-native\n");
        return 1;
    }
    const bool species_panel_may_be_empty =
        has_species_vcf && species_mode == SpeciesPanelMode::COUNT_ONLY &&
        !species_counts_output && !species_assignment_output;
    
    // BAM header/TID mapping is read directly through HTSlib below when needed.

    const string final_output_prefix = output_prefix;
    OutputTransaction output_transaction;
    output_prefix = output_transaction.stage_prefix(final_output_prefix);
    if (dump_pileup){
        dump_pileup_prefix = output_transaction.stage_prefix(dump_pileup_prefix);
    }
    if (dump_source_observations){
        dump_source_observations_prefix =
            output_transaction.stage_prefix(dump_source_observations_prefix);
    }

    // Check if we can load counts from a previous successful run. All new
    // artifacts are staged under temporary names and published only at success.
    bool load_counts = false;
    string countsfilename = final_output_prefix + ".counts";
    if (file_exists(countsfilename)){
        if (force_recount){
            fprintf(stderr, "Found existing counts file: %s\n", countsfilename.c_str());
            fprintf(stderr,
                "  --force_recount set: recounting into staged outputs; existing final files remain intact until success\n");
        } else if (reuse_counts){
            // Validate the existing file before trusting it
            fprintf(stderr, "Found existing counts file: %s\n", countsfilename.c_str());
            fprintf(stderr, "  --reuse_counts set: validating before loading\n");

            // Check .samples companion
            string samplesfile = final_output_prefix + ".samples";
            if (!file_exists(samplesfile)){
                fprintf(stderr, "ERROR: --reuse_counts: .samples file missing alongside .counts\n");
                fprintf(stderr, "  Delete %s and rerun, or use --force_recount\n", countsfilename.c_str());
                return 1;
            }

            // Check .counts is gzip-readable (not truncated)
            if (!validate_gzip_readable(countsfilename)){
                fprintf(stderr, "ERROR: --reuse_counts: .counts file appears truncated or corrupt: %s\n",
                    countsfilename.c_str());
                fprintf(stderr, "  Delete it and rerun, or use --force_recount\n");
                return 1;
            }

            if (accepted_weighted_conditional){
                const string condf_file = final_output_prefix + ".condf";
                const string condf_basis = final_output_prefix + ".condf_basis.tsv";
                if (!file_exists(condf_file) || !file_exists(condf_basis) ||
                    !text_file_contains(condf_basis, "weighting\taccepted_observation_weighted") ||
                    !text_file_contains(condf_basis, "truth_inputs_used\tfalse")){
                    fprintf(stderr,
                        "ERROR: --reuse_counts with --accepted_weighted_conditional requires "
                        "a truth-free accepted-observation-weighted .condf and .condf_basis.tsv\n");
                    fprintf(stderr, "  Use --force_recount to regenerate a matched-basis bundle.\n");
                    return 1;
                }
            }

            // If species output is requested, the species companion outputs must
            // be present and readable too.  Otherwise a reused main .counts can
            // silently coexist with stale or missing .species_counts/.species_condf.
            if (species_counts_output){
                string sp_counts = final_output_prefix + ".species_counts";
                string sp_condf = final_output_prefix + ".species_condf";
                if (!file_exists(sp_counts)){
                    fprintf(stderr, "ERROR: --reuse_counts with --species_counts_output: missing %s\n",
                        sp_counts.c_str());
                    fprintf(stderr, "  Use --force_recount to regenerate the full dual-panel output set.\n");
                    return 1;
                }
                if (!file_exists(sp_condf)){
                    fprintf(stderr, "ERROR: --reuse_counts with --species_counts_output: missing %s\n",
                        sp_condf.c_str());
                    fprintf(stderr, "  Use --force_recount to regenerate the full dual-panel output set.\n");
                    return 1;
                }
                if (accepted_weighted_conditional){
                    const string sp_basis = final_output_prefix + ".species_condf_basis.tsv";
                    if (!file_exists(sp_basis) ||
                        !text_file_contains(sp_basis, "weighting\taccepted_observation_weighted") ||
                        !text_file_contains(sp_basis, "truth_inputs_used\tfalse")){
                        fprintf(stderr,
                            "ERROR: reused species bundle lacks a truth-free "
                            "accepted-observation-weighted .species_condf_basis.tsv\n");
                        fprintf(stderr, "  Use --force_recount to regenerate the full matched-basis bundle.\n");
                        return 1;
                    }
                }
                if (!validate_gzip_readable(sp_counts)){
                    fprintf(stderr, "ERROR: --reuse_counts: .species_counts appears truncated or corrupt: %s\n",
                        sp_counts.c_str());
                    fprintf(stderr, "  Use --force_recount to regenerate the full dual-panel output set.\n");
                    return 1;
                }
            }

            load_counts = true;
            fprintf(stderr, "  Validation passed; loading existing counts\n");
        } else {
            // Default: refuse to silently reuse stale counts
            fprintf(stderr, "ERROR: existing .counts file found: %s\n", countsfilename.c_str());
            fprintf(stderr, "  A previous run may have left this file truncated or stale.\n");
            fprintf(stderr, "  Use --reuse_counts to load it after validation,\n");
            fprintf(stderr, "  or --force_recount to delete it and recount from BAM.\n");
            return 1;
        }
    }
    else{
        if (bamfile.length() == 0 && !dump_conditional){
            fprintf(stderr, "ERROR: bam file (--bam) required\n");
            return 1;
        }
        // BAM readability/header validation is performed by direct HTSlib helpers
        // when chromosome/TID mapping is needed.
    }

    if (accepted_weighted_conditional && disable_conditional){
        fprintf(stderr, "ERROR: --accepted_weighted_conditional is incompatible with --disable_conditional/-f\n");
        return 1;
    }
    if (accepted_weighted_conditional && dump_conditional){
        fprintf(stderr, "ERROR: --accepted_weighted_conditional requires a BAM recount and cannot be used with VCF-only --dump_conditional\n");
        return 1;
    }
    if (accepted_weighted_conditional && n_threads <= 1 && !load_counts){
        fprintf(stderr, "ERROR: --accepted_weighted_conditional requires --threads > 1 for a BAM recount\n");
        return 1;
    }
    // --dump_pileup fail-closed guards. The pileup sidecars are produced inside
    // the BAM counting pass, so a run that loads existing counts (or takes the
    // single-threaded fallback, which has no pileup support) would silently
    // publish a demux bundle without the promised pileup files. Refuse instead.
    if (dump_pileup && load_counts){
        fprintf(stderr, "ERROR: --dump_pileup requires a BAM counting pass; existing counts were loaded\n");
        fprintf(stderr, "  Rerun with --force_recount to produce the pileup sidecars\n");
        return 1;
    }
    if (dump_pileup && n_threads <= 1){
        fprintf(stderr, "ERROR: --dump_pileup requires --threads > 1 (parallel counting paths only)\n");
        return 1;
    }

    // Load sample names
    vector<string> samples;
    bool samples_from_vcf = false;

    if (load_counts){
        string samplesfile = final_output_prefix + ".samples";
        if (file_exists(samplesfile)){
           load_samples(samplesfile, samples); 
        }
        else{
            if (vcf_file == ""){
                fprintf(stderr, "ERROR: vcf file is required\n");
                return 1;
            }
            read_vcf_samples(vcf_file, samples);         
        }
    }
    else{
        if (vcf_file == "" && shared_vcf_name == ""){
            fprintf(stderr, "ERROR: vcf file is required\n");
            return 1;
        }
        if (vcf_file != ""){
            read_vcf_samples(vcf_file, samples);
        }
        samples_from_vcf = true;
    }
    
    // Parse allowed IDs - deferred until after VCF loading if using shared memory
    set<int> allowed_ids;
    set<int> allowed_ids2;
    
    // Process samples now if not using shared VCF (otherwise defer until after attach).
    // EXCEPTION: when load_counts is true, sample names are already populated from
    // the .samples companion file above (or read_vcf_samples), so we can and must
    // parse the idfile here even when a shared VCF is in use. The deferred parse
    // further below only runs inside the (!load_counts) branch, so without this the
    // combination (--reuse_counts + --shared_vcf + -i/-I) would leave allowed_ids
    // empty and silently disable the identity restriction, letting the assignment
    // step consider every possible singlet/doublet combination.
    if (shared_vcf_name.length() == 0 || load_counts){
        fprintf(stderr, "Number of individuals in VCF: %lu\n", samples.size());
        
        if (idfile_given){
            parse_idfile(idfile, samples, allowed_ids, allowed_ids2, true);
            if (allowed_ids.size() == 0){
                fprintf(stderr, "ERROR: No valid individual names found in %s; refusing to ignore -i\n", idfile.c_str());
                return 1;
            }
        }
        if (idfile_doublet_given){
            parse_idfile(idfile_doublet, samples, allowed_ids, allowed_ids2, false);
            if (allowed_ids.size() == 0){
                fprintf(stderr, "ERROR: No valid individual names found in %s; refusing to ignore -I\n", idfile_doublet.c_str());
                return 1;
            }
        }

        if (idfile_given || idfile_doublet_given){
            fprintf(stderr, "Identity restriction active: allowed_ids=%lu allowed_ids2=%lu\n",
                allowed_ids.size(), allowed_ids2.size());
        }
        
        if (samples_from_vcf){
            string samplesfile = output_prefix + ".samples";
            if (!write_samples_checked(samplesfile, samples)) return 1;
        }
    }
    
    // Handle --dump_conditional / -F mode: VCF-only .condf generation
    if (dump_conditional){
        if (vcf_file.length() == 0 && shared_vcf_name.length() == 0){
            fprintf(stderr, "ERROR: --vcf/-v or --shared_vcf/-S required for --dump_conditional/-F\n");
            return 1;
        }
        
        robin_hood::unordered_map<int, ChromSNPs> snpdat_optimized;
        
        if (shared_vcf_name.length() > 0){
            // Fast path: attach to shared memory daemon
            print_elapsed(start_time, "dump_conditional mode: attaching to shared VCF...");
            fprintf(stderr, "Attaching to shared VCF: %s\n", shared_vcf_name.c_str());
            
            if (!attach_shared_vcf(shared_vcf_name, snpdat_optimized, samples)){
                fprintf(stderr, "ERROR: Could not attach to shared VCF: %s\n", shared_vcf_name.c_str());
                return 1;
            }
            fprintf(stderr, "Attached. %lu chromosomes, %lu samples\n",
                snpdat_optimized.size(), samples.size());
        }
        else{
            // Disk path: load VCF from file
            print_elapsed(start_time, "dump_conditional mode: loading VCF from disk...");
            
            set<string> chroms_vcf;
            get_vcf_chroms(vcf_file, chroms_vcf);
            
            // Build synthetic seq2tid (no BAM header available)
            map<string, int> seq2tid_synthetic;
            int tid_counter = 0;
            for (const auto& c : chroms_vcf){
                seq2tid_synthetic[c] = tid_counter++;
            }
            
            int nloaded = read_vcf_chroms_optimized(vcf_file, chroms_vcf,
                seq2tid_synthetic, snpdat_optimized, vq);
            if (nloaded < 0) return 1;
            fprintf(stderr, "Loaded %d SNPs from %lu chromosomes\n", nloaded, chroms_vcf.size());
        }
        if (!validate_panel_without_bam(
                snpdat_optimized, "main genotype panel", false)) return 1;
        
        // Compute conditional match fractions. The parallel CONDF function now
        // performs genotype-only precompute internally, avoiding unnecessary
        // per-SNP pair_targets allocation in VCF-only mode.
        print_elapsed(start_time, "Computing conditional match fractions...");
        map<pair<int, int>, map<int, float>> conditional_match_fracs;
        
        compute_conditional_match_fracs_parallel(snpdat_optimized,
            conditional_match_fracs, samples.size(), n_threads);
        
        // Write to a staged name and publish only after the complete operation succeeds.
        string outname = output_prefix + ".condf";
        FILE* outf = fopen(outname.c_str(), "w");
        if (!outf){
            fprintf(stderr, "ERROR: could not open staged .condf for writing: %s\n",
                strerror(errno));
            return 1;
        }
        dump_exp_fracs(outf, conditional_match_fracs);
        if (fclose(outf) != 0){
            fprintf(stderr, "ERROR: failed closing staged .condf: %s\n", strerror(errno));
            return 1;
        }
        if (!publish_outputs(output_transaction, "conditional-fraction")) return 1;

        print_elapsed(start_time, "Done. Wrote .condf file.");
        fprintf(stderr, "Wrote conditional match fractions to %s.condf\n",
            final_output_prefix.c_str());
        return 0;
    }
    
    // Load cell barcodes
    set<unsigned long> cell_barcodes;
    if (cell_barcode){
        parse_barcode_file(cell_barcode_file, cell_barcodes);
        if (cell_barcodes.size() == 0){
            fprintf(stderr, "ERROR reading cell barcode list\n");
            return 1;
        }
        fprintf(stderr, "Loaded %lu cell barcodes\n", cell_barcodes.size());
    }

    // Main data structures
    robin_hood::unordered_map<unsigned long, CellCounts> cell_counts;
    
    map<pair<int, int>, map<int, float> > conditional_match_fracs;
    map<pair<int, int>, map<int, float> > conditional_match_tots;
    AcceptedSiteWeightMap accepted_site_weights_individual;
    AcceptedSiteWeightMap accepted_site_weights_species;
    ConditionalWeightStats individual_condf_stats;
    ConditionalWeightStats species_condf_stats;

    // Position set for species panel dedup (§5.5)
    // Populated from demux VCF snpdat before it goes out of scope;
    // used later to filter overlapping sites from the species VCF.
    // Key: (tid, position)
    set<pair<int, int>> demux_positions;
    
    // Species panel data - declared early so dual counting path can populate them
    robin_hood::unordered_map<int, ChromSNPs> species_snpdat;
    robin_hood::unordered_map<unsigned long, CellCounts> species_cell_counts_rna;
    robin_hood::unordered_map<unsigned long, CellCounts> species_cell_counts_atac;
    robin_hood::unordered_map<unsigned long, CellCounts> species_cell_counts_native;
    vector<string> species_samples_native;
    PanelMetadata panel_meta;
    bool panel_meta_loaded = false;
    bool species_counted_dual = false;  // Set by WP3 dual counting path

    // Main BAM seqname->TID mapping, populated once for BAM-backed runs.
    // Keep it outside the counting block because later disk-based species/het
    // VCF loaders need the same TID mapping.
    map<string, int> bam_seq2tid;
    
    if (load_counts){
        // Load counts from previous run
        print_elapsed(start_time, "Loading counts from existing file...");
        fprintf(stderr, "Loading allele counts from %s\n", countsfilename.c_str());
        int n_lines = load_cellcounts_optimized(countsfilename, cell_counts, samples.size());
        fprintf(stderr, "Loaded %d count records for %lu cells\n", n_lines, cell_counts.size());
        
        // Generate a staged .condf only when the successful final bundle lacks one.
        const string existing_condf_file = final_output_prefix + ".condf";
        string condf_file = output_prefix + ".condf";
        if (!disable_conditional && !file_exists(existing_condf_file)){
            if (vcf_file.length() == 0 && shared_vcf_name.length() == 0){
                fprintf(stderr, "WARNING: No VCF or shared VCF provided and no .condf file found. "
                                "Skipping conditional match fraction computation.\n");
            }
            else{
                robin_hood::unordered_map<int, ChromSNPs> snpdat_optimized;
                
                if (shared_vcf_name.length() > 0){
                    // Attach to shared memory daemon
                    print_elapsed(start_time, "Generating .condf: attaching to shared VCF...");
                    fprintf(stderr, "Attaching to shared VCF: %s\n", shared_vcf_name.c_str());
                    
                    vector<string> shm_samples;
                    if (!attach_shared_vcf(shared_vcf_name, snpdat_optimized, shm_samples)){
                        fprintf(stderr, "ERROR: Could not attach to shared VCF needed for .condf\n");
                        return 1;
                    }
                }
                else{
                    // Load from disk
                    print_elapsed(start_time, "Generating .condf from VCF on disk...");
                    
                    set<string> chroms_vcf;
                    get_vcf_chroms(vcf_file, chroms_vcf);
                    
                    map<string, int> seq2tid_synthetic;
                    int tid_counter = 0;
                    for (const auto& c : chroms_vcf){
                        seq2tid_synthetic[c] = tid_counter++;
                    }
                    
                    int nloaded = read_vcf_chroms_optimized(vcf_file, chroms_vcf,
                        seq2tid_synthetic, snpdat_optimized, vq);
                    if (nloaded < 0) return 1;
                    fprintf(stderr, "Loaded %d SNPs for .condf computation\n", nloaded);
                }
                if (!validate_panel_without_bam(
                        snpdat_optimized, "main genotype panel", false)) return 1;
                // The parallel CONDF function performs genotype-only precompute
                // internally, avoiding unnecessary pair_targets allocation here.
                compute_conditional_match_fracs_parallel(snpdat_optimized,
                    conditional_match_fracs, samples.size(), n_threads);
                
                FILE* outf = fopen(condf_file.c_str(), "w");
                if (!outf){
                    fprintf(stderr, "ERROR: could not open staged .condf for writing: %s\n",
                        strerror(errno));
                    return 1;
                }
                dump_exp_fracs(outf, conditional_match_fracs);
                if (fclose(outf) != 0){
                    fprintf(stderr, "ERROR: failed closing staged .condf: %s\n",
                        strerror(errno));
                    return 1;
                }
                fprintf(stderr, "Prepared replacement .condf for %s\n",
                    final_output_prefix.c_str());
            }
        }
    }
    else{
        // Determine chromosomes to process
        set<string> chroms_bam;
        set<string> chroms_vcf;
        set<string> chroms_to_process;
        
        fprintf(stderr, "Identifying shared chromosomes between BAM and VCF...\n");
        string bam_header_error;
        if (!get_bam_header_chroms_and_seq2tid(
                bamfile, chroms_bam, bam_seq2tid, &bam_header_error)){
            fprintf(stderr, "ERROR: %s\n", bam_header_error.c_str());
            return 1;
        }
        fprintf(stderr, "BAM header loaded: %lu chromosomes\n", chroms_bam.size());
        
        if (shared_vcf_name.length() > 0){
            // Using shared memory VCF - chroms will be determined after loading
        }
        else{
            get_vcf_chroms(vcf_file, chroms_vcf);
            
            for (auto& c : chroms_bam){
                if (chroms_vcf.find(c) != chroms_vcf.end()){
                    chroms_to_process.insert(c);
                }
            }
            
            if (vcf_chroms_given){
                set<string> user_chroms;
                ifstream chromfile(vcf_chroms_file.c_str());
                string chrom;
                while (chromfile >> chrom){
                    user_chroms.insert(chrom);
                }
                set<string> filtered_chroms;
                for (auto& c : chroms_to_process){
                    if (user_chroms.find(c) != user_chroms.end()){
                        filtered_chroms.insert(c);
                    }
                }
                chroms_to_process = filtered_chroms;
            }
            
            fprintf(stderr, "Found %lu chromosomes in BAM, %lu in VCF, %lu shared\n",
                chroms_bam.size(), chroms_vcf.size(), chroms_to_process.size());
            if (chroms_to_process.empty()){
                fprintf(stderr,
                    "ERROR: required main genotype panel has zero shared BAM/VCF contigs\n");
                return 1;
            }
        }
        
        // Load VCF data
        print_elapsed(start_time, "Starting VCF loading...");
        robin_hood::unordered_map<int, ChromSNPs> snpdat_optimized;
        
        if (shared_vcf_name.length() > 0){
            // Attach to shared memory VCF
            fprintf(stderr, "Attaching to shared VCF: %s\n", shared_vcf_name.c_str());
            if (!attach_shared_vcf(shared_vcf_name, snpdat_optimized, samples)){
                fprintf(stderr, "ERROR: Could not attach to shared VCF\n");
                return 1;
            }
            if (!validate_and_restrict_panel_to_bam(
                    snpdat_optimized, bam_seq2tid, "main genotype panel", false)) return 1;
            
            // Deferred samples processing (samples now populated from shared memory)
            fprintf(stderr, "Number of individuals in VCF: %lu\n", samples.size());
            
            if (idfile_given){
                parse_idfile(idfile, samples, allowed_ids, allowed_ids2, true);
                if (allowed_ids.size() == 0){
                    fprintf(stderr, "ERROR: No valid individual names found in %s; refusing to ignore -i\n", idfile.c_str());
                    return 1;
                }
            }
            if (idfile_doublet_given){
                parse_idfile(idfile_doublet, samples, allowed_ids, allowed_ids2, false);
                if (allowed_ids.size() == 0){
                    fprintf(stderr, "ERROR: No valid individual names found in %s; refusing to ignore -I\n", idfile_doublet.c_str());
                    return 1;
                }
            }

            if (idfile_given || idfile_doublet_given){
                fprintf(stderr, "Identity restriction active: allowed_ids=%lu allowed_ids2=%lu\n",
                    allowed_ids.size(), allowed_ids2.size());
            }
            
            // Write samples file
            string samplesfile = output_prefix + ".samples";
            if (!write_samples_checked(samplesfile, samples)) return 1;
            precompute_all_genotypes(snpdat_optimized, samples.size());
        }
        else{
            // Load VCF into memory (default behavior)
            // Use --no_preload to disable if memory is limited
            if (!no_preload){
                print_elapsed(start_time, "Loading VCF data into memory...");
                fprintf(stderr, "Loading VCF data into memory...\n");
                int nloaded = read_vcf_chroms_optimized(vcf_file, chroms_to_process, 
                    bam_seq2tid, snpdat_optimized, vq);
                if (nloaded < 0) return 1;
                print_elapsed(start_time, "VCF loading complete");
                fprintf(stderr, "Loaded %d SNPs from %lu chromosomes\n", 
                    nloaded, chroms_to_process.size());
                if (!validate_and_restrict_panel_to_bam(
                        snpdat_optimized, bam_seq2tid, "main genotype panel", false)) return 1;
                precompute_all_genotypes(snpdat_optimized, samples.size());
            }
            else{
                fprintf(stderr, "VCF preloading disabled. Use --shared_vcf for large datasets.\n");
                return 1;
            }
        }
        
        // Count alleles
        print_elapsed(start_time, "Starting allele counting...");
        fprintf(stderr, "Counting alleles in BAM file...\n");
        
        // WP3 optimization: detect when both interindividual and species panels
        // need counting from the same RNA BAM. If so, load species VCF early,
        // merge SNP sets, and do a single BAM pass via count_alleles_parallel_dual.
        bool use_dual_counting = false;
        if (n_threads > 1 && has_species_vcf 
            && (species_mode == SpeciesPanelMode::AUGMENT 
                || species_mode == SpeciesPanelMode::BOTH
                || species_counts_output)){
            use_dual_counting = true;
        }
        
        if (use_dual_counting){
            // === DUAL-PANEL PATH (WP3) ===
            // Step 1: Build demux_positions for dedup BEFORE species VCF load
            for (auto& kv : snpdat_optimized){
                for (auto& snp : kv.second.snps){
                    demux_positions.insert(make_pair(kv.first, snp.pos));
                }
            }
            fprintf(stderr, "Built demux position set: %lu sites for species dedup\n",
                demux_positions.size());
            
            // Step 2: Load species VCF early
            print_elapsed(start_time, "Loading species panel VCF for dual counting...");
            vector<string> species_samples_early;
            if (species_shared_vcf_name.length() > 0){
                if (!attach_shared_vcf(species_shared_vcf_name, species_snpdat, species_samples_early)){
                    fprintf(stderr, "ERROR: Could not attach to shared species VCF: %s\n",
                        species_shared_vcf_name.c_str());
                    return 1;
                }
                if (!validate_and_restrict_panel_to_bam(
                        species_snpdat, bam_seq2tid, "species genotype panel", false)) return 1;
            }
            else{
                read_vcf_samples(species_vcf_file, species_samples_early);
                set<string> sp_chroms;
                get_vcf_chroms(species_vcf_file, sp_chroms);
                set<string> sp_shared_chroms = shared_contig_names(sp_chroms, bam_seq2tid);
                if (sp_shared_chroms.empty()){
                    fprintf(stderr,
                        "ERROR: required species genotype panel has zero shared BAM/VCF contigs\n");
                    return 1;
                }
                int nloaded = read_vcf_chroms_optimized(species_vcf_file, sp_shared_chroms,
                    bam_seq2tid, species_snpdat, vq);
                if (nloaded < 0) return 1;
                fprintf(stderr, "Loaded %d species panel SNPs\n", nloaded);
            }
            
            // Step 3: Sample set alignment check
            if (species_samples_early.size() != samples.size()){
                fprintf(stderr, "ERROR: Species VCF has %lu samples but demux VCF has %lu samples. "
                    "Sample sets must match.\n", species_samples_early.size(), samples.size());
                return 1;
            }
            for (size_t si = 0; si < samples.size(); ++si){
                if (species_samples_early[si] != samples[si]){
                    fprintf(stderr, "ERROR: Species VCF sample[%lu]='%s' differs from demux VCF sample[%lu]='%s'. "
                        "Sample sets must match in same order.\n",
                        si, species_samples_early[si].c_str(), si, samples[si].c_str());
                    return 1;
                }
            }
            fprintf(stderr, "Species VCF sample set matches demux VCF (%lu samples)\n", samples.size());
            
            // Step 4: Dedup species SNPs against demux positions
            if ((species_mode == SpeciesPanelMode::AUGMENT || species_mode == SpeciesPanelMode::BOTH)
                && !demux_positions.empty()){
                int overlap_count = 0;
                int total_species_sites = 0;
                for (auto& sp_chrom : species_snpdat){
                    int tid = sp_chrom.first;
                    vector<SNPData> filtered_snps;
                    for (auto& sp_snp : sp_chrom.second.snps){
                        total_species_sites++;
                        if (demux_positions.count(make_pair(tid, sp_snp.pos)) > 0){
                            overlap_count++;
                        }
                        else{
                            filtered_snps.push_back(sp_snp);
                        }
                    }
                    sp_chrom.second.snps = filtered_snps;
                }
                if (overlap_count > 0){
                    fprintf(stderr, "Species/demux overlap: %d of %d species sites removed "
                        "(demux VCF retains priority)\n", overlap_count, total_species_sites);
                }
                else{
                    fprintf(stderr, "Species/demux overlap: 0 of %d species sites overlap\n",
                        total_species_sites);
                }
            }
            if (!validate_and_restrict_panel_to_bam(
                    species_snpdat, bam_seq2tid, "species genotype panel", false)) return 1;
            
            if (species_counts_output || species_assignment_output){
                load_panel_metadata_if_needed(
                    panel_metadata_file, samples, panel_meta, panel_meta_loaded);
                species_samples_native = panel_meta.species_list;
            }

            // Step 5: Tag species SNPs with panel_id = 1
            for (auto& kv : species_snpdat){
                for (auto& snp : kv.second.snps){
                    snp.panel_id = 1;
                }
            }
            
            // Step 6: Merge into combined SNP set
            robin_hood::unordered_map<int, ChromSNPs> combined_snpdat = snpdat_optimized;
            for (auto& kv : species_snpdat){
                int tid = kv.first;
                auto& combined_chrom = combined_snpdat[tid];
                for (auto& snp : kv.second.snps){
                    combined_chrom.snps.push_back(snp);
                }
                combined_chrom.sort_snps();
            }
            
            long n_interindiv = 0, n_species = 0;
            for (auto& kv : combined_snpdat){
                for (auto& snp : kv.second.snps){
                    if (snp.panel_id == 0) n_interindiv++;
                    else n_species++;
                }
            }
            fprintf(stderr, "Combined SNP set: %ld interindividual + %ld species = %ld total\n",
                n_interindiv, n_species, n_interindiv + n_species);
            if (n_interindiv <= 0){
                fprintf(stderr,
                    "ERROR: required main genotype panel contains zero usable SNPs in dual-panel path\n");
                return 1;
            }
            if (n_species <= 0){
                fprintf(stderr,
                    "ERROR: required species genotype panel contains zero usable SNPs in dual-panel path\n");
                return 1;
            }
            
            // Precompute individual-shaped targets for both panels. Native
            // species targets remain in a separate compact table so SNPData
            // and the shared-VCF memory layout are unchanged.
            NativeSpeciesTargetTable species_native_targets;
            precompute_all_genotypes(combined_snpdat, samples.size());
            if (species_counts_output || species_assignment_output){
                precompute_native_species_targets(
                    combined_snpdat, panel_meta, samples.size(),
                    species_native_targets, 1);
            }
            
            // Step 7: Single BAM pass with dual output
            robin_hood::unordered_map<unsigned long, AlignedCellCounts> dual_panel0;
            robin_hood::unordered_map<unsigned long, AlignedCellCounts> dual_panel1;
            robin_hood::unordered_map<unsigned long, AlignedCellCounts> dual_panel1_native;
            if (!count_alleles_parallel_dual(bamfile, combined_snpdat,
                dual_panel0, dual_panel1,
                cell_barcodes, samples.size(), samples, n_threads, htslib_threads,
                dump_source_observations, dump_source_observations_prefix, source_provenance_tag,
                synthetic_id_tag, expected_synthetic_id, source_receiver_map,
                source_reconciliation_mode, source_donor_site_audit,
                source_donor_site_sample_mod,
                accepted_weighted_conditional ? &accepted_site_weights_individual : nullptr,
                accepted_weighted_conditional ? &accepted_site_weights_species : nullptr,
                (species_counts_output || species_assignment_output)
                    ? &species_native_targets : nullptr,
                (species_counts_output || species_assignment_output)
                    ? &dual_panel1_native : nullptr,
                (species_counts_output || species_assignment_output)
                    ? (int)species_samples_native.size() : 0,
                dump_pileup, dump_pileup_prefix)){
                fprintf(stderr, "ERROR: dual-panel BAM counting failed; no outputs will be published\n");
                return 1;
            }
            
            // Step 8: Finalize each panel separately
            finalize_parallel_counts(dual_panel0, cell_counts);
            finalize_parallel_counts(dual_panel1, species_cell_counts_rna);
            if (species_counts_output || species_assignment_output){
                finalize_parallel_counts(dual_panel1_native, species_cell_counts_native);
            }
            
            fprintf(stderr,
                "Dual counting complete: %lu interindiv cells, %lu individual-shaped species cells, %lu native species cells\n",
                cell_counts.size(), species_cell_counts_rna.size(),
                species_cell_counts_native.size());
            
            // Mark species RNA counting as done so the later block skips it
            species_counted_dual = true;
        }
        else if (n_threads > 1){
            if (dump_source_observations){
                fprintf(stderr, "ERROR: --dump_source_observations currently requires the dual-panel counting path (--species_counts_output).\n");
                return 1;
            }
            // Parallel processing
            robin_hood::unordered_map<unsigned long, AlignedCellCounts> parallel_counts;
        
            if (!count_alleles_parallel(bamfile, snpdat_optimized, parallel_counts,
                cell_barcodes, samples.size(), n_threads, htslib_threads,
                dump_pileup, dump_pileup_prefix,
                accepted_weighted_conditional ? &accepted_site_weights_individual : nullptr)){
                fprintf(stderr, "ERROR: BAM allele counting failed; no outputs will be published\n");
                return 1;
            }
        
            // Finalize counts
            finalize_parallel_counts(parallel_counts, cell_counts);
        }
        else{
            // Single-threaded fallback
            if (!count_alleles_single_threaded(bamfile, snpdat_optimized, cell_counts,
                cell_barcodes, samples.size(), conditional_match_fracs,
                conditional_match_tots, !disable_conditional)){
                fprintf(stderr, "ERROR: BAM allele counting failed; no outputs will be published\n");
                return 1;
            }
        }
    
        // Compute conditional match fractions if parallel (didn't do it during counting)
        if (n_threads > 1 && !disable_conditional){
            if (accepted_weighted_conditional){
                if (accepted_site_weights_individual.empty()){
                    fprintf(stderr, "ERROR: no accepted individual-panel site weights were collected; cannot build accepted-weighted .condf\n");
                    return 1;
                }
                fprintf(stderr, "Computing accepted-observation-weighted conditional match fractions...\n");
                individual_condf_stats = compute_conditional_match_fracs_weighted(
                    snpdat_optimized, accepted_site_weights_individual,
                    conditional_match_fracs, samples.size(), n_threads);
            }
            else{
                fprintf(stderr, "Computing VCF-site-weighted conditional match fractions...\n");
                compute_conditional_match_fracs_parallel(snpdat_optimized,
                    conditional_match_fracs, samples.size(), n_threads);
            }
        }
    
        if (!disable_conditional){
            // For single-threaded path, condf was accumulated during counting
            // and still needs normalization. For parallel path, the new function
            // already normalized internally, so this is a no-op (tots map is empty).
            conditional_match_fracs_normalize(conditional_match_fracs,
                conditional_match_tots, samples.size());
        
            string outname = output_prefix + ".condf";
            FILE* outf = fopen(outname.c_str(), "w");
            if (!outf){
                fprintf(stderr, "ERROR: could not open %s for writing: %s\n",
                    outname.c_str(), strerror(errno));
                return 1;
            }
            dump_exp_fracs(outf, conditional_match_fracs);
            if (fclose(outf) != 0){
                fprintf(stderr, "ERROR: failed closing %s: %s\n",
                    outname.c_str(), strerror(errno));
                return 1;
            }

            string basis_name = output_prefix + ".condf_basis.tsv";
            FILE* basis_out = fopen(basis_name.c_str(), "w");
            if (!basis_out){
                fprintf(stderr, "ERROR: could not open %s for writing\n", basis_name.c_str());
                return 1;
            }
            fprintf(basis_out, "contract_version\tcondf_basis_V1_R1\n");
            fprintf(basis_out, "panel\tindividual\n");
            fprintf(basis_out, "weighting\t%s\n",
                accepted_weighted_conditional ? "accepted_observation_weighted" : "vcf_site_weighted");
            fprintf(basis_out, "truth_inputs_used\tfalse\n");
            fprintf(basis_out, "observed_sites\t%llu\n",
                (unsigned long long)individual_condf_stats.observed_sites);
            fprintf(basis_out, "accepted_weight\t%.17g\n", individual_condf_stats.accepted_weight);
            if (fclose(basis_out) != 0){
                fprintf(stderr, "ERROR: failed closing %s: %s\n",
                    basis_name.c_str(), strerror(errno));
                return 1;
            }
        }

        // Build position set from demux VCF for species panel dedup (§5.5)
        if (has_species_vcf){
            for (auto& kv : snpdat_optimized){
                for (auto& snp : kv.second.snps){
                    demux_positions.insert(make_pair(kv.first, snp.pos));
                }
            }
            fprintf(stderr, "Built demux position set: %lu sites for species dedup\n",
                demux_positions.size());
        }
    } // end else (not load_counts)
    
    // Write counts to disk (skip if we loaded from existing file)
    if (!load_counts){
        print_elapsed(start_time, "Writing allele counts to disk...");
        {
            string fname = output_prefix + ".counts";
            gzFile outf = gzopen(fname.c_str(), "w");
            if (!outf){
                fprintf(stderr, "ERROR: could not open %s for writing: %s\n",
                    fname.c_str(), strerror(errno));
                return 1;
            }
            fprintf(stderr, "Writing allele counts to staged output...\n");
            dump_cellcounts_optimized(outf, cell_counts, samples.size());
            if (gzclose(outf) != Z_OK){
                fprintf(stderr, "ERROR: failed closing %s\n", fname.c_str());
                return 1;
            }
            fprintf(stderr, "Done writing staged counts\n");
        }
    }

    // ================================================================
    // 2A: ATAC dual-modality counting
    // ================================================================
    robin_hood::unordered_map<unsigned long, CellCounts> atac_cell_counts;

    if (atac_mode && !load_counts){
        print_elapsed(start_time, "Starting ATAC allele counting...");
        robin_hood::unordered_map<int, ChromSNPs> atac_snpdat;
        map<string, int> seq2tid_atac;
        set<string> chroms_atac_header;
        string atac_header_error;
        if (!get_bam_header_chroms_and_seq2tid(
                atac_bamfile, chroms_atac_header, seq2tid_atac, &atac_header_error)){
            fprintf(stderr, "ERROR: ATAC BAM/header validation failed: %s\n",
                atac_header_error.c_str());
            return 1;
        }

        if (atac_shared_vcf_name.length() > 0){
            vector<string> atac_samples;
            if (!attach_shared_vcf(atac_shared_vcf_name, atac_snpdat, atac_samples)){
                fprintf(stderr, "ERROR: Could not attach to shared ATAC VCF: %s\n",
                    atac_shared_vcf_name.c_str());
                return 1;
            }
            if (!validate_and_restrict_panel_to_bam(
                    atac_snpdat, seq2tid_atac, "ATAC genotype panel", false)) return 1;
            fprintf(stderr, "Attached to shared ATAC VCF with %lu chromosomes\n", atac_snpdat.size());
        }
        else{
            set<string> atac_chroms;
            get_vcf_chroms(atac_vcf_file, atac_chroms);
            set<string> atac_shared_chroms = shared_contig_names(atac_chroms, seq2tid_atac);
            if (atac_shared_chroms.empty()){
                fprintf(stderr,
                    "ERROR: required ATAC genotype panel has zero shared BAM/VCF contigs\n");
                return 1;
            }
            int nloaded = read_vcf_chroms_optimized(atac_vcf_file, atac_shared_chroms,
                seq2tid_atac, atac_snpdat, vq);
            if (nloaded < 0) return 1;
            fprintf(stderr, "Loaded %d ATAC SNPs\n", nloaded);
            if (!validate_and_restrict_panel_to_bam(
                    atac_snpdat, seq2tid_atac, "ATAC genotype panel", false)) return 1;
        }

        precompute_all_genotypes(atac_snpdat, samples.size());

        if (n_threads > 1){
            robin_hood::unordered_map<unsigned long, AlignedCellCounts> atac_parallel_counts;
            if (!count_alleles_parallel(atac_bamfile, atac_snpdat, atac_parallel_counts,
                    cell_barcodes, samples.size(), n_threads, htslib_threads)){
                fprintf(stderr, "ERROR: ATAC allele counting failed\n");
                return 1;
            }
            finalize_parallel_counts(atac_parallel_counts, atac_cell_counts);
        }
        else{
            map<pair<int, int>, map<int, float>> atac_cond_fracs, atac_cond_tots;
            if (!count_alleles_single_threaded(atac_bamfile, atac_snpdat, atac_cell_counts,
                    cell_barcodes, samples.size(), atac_cond_fracs, atac_cond_tots, false)){
                fprintf(stderr, "ERROR: ATAC allele counting failed\n");
                return 1;
            }
        }

        // CB tag intersection check
        int n_intersect = 0;
        int n_rna = (int)cell_counts.size();
        int n_atac = (int)atac_cell_counts.size();
        for (auto& kv : atac_cell_counts){
            if (cell_counts.count(kv.first) > 0) n_intersect++;
        }
        int min_set = (n_rna < n_atac) ? n_rna : n_atac;
        if (min_set > 0){
            double overlap_frac = (double)n_intersect / (double)min_set;
            fprintf(stderr, "ATAC/RNA barcode overlap: %d / %d = %.3f\n",
                n_intersect, min_set, overlap_frac);
            if (overlap_frac < 0.10){
                fprintf(stderr, "ERROR: ATAC/RNA barcode overlap (%.3f) below threshold (0.10).\n"
                    "Check that the ATAC BAM CB tag holds RNA-aligned barcodes.\n", overlap_frac);
                return 1;
            }
        }

        // Write ATAC counts
        {
            string fname = output_prefix + ".atac.counts";
            gzFile outf = gzopen(fname.c_str(), "w");
            if (!outf){
                fprintf(stderr, "ERROR: could not open %s for writing: %s\n",
                    fname.c_str(), strerror(errno));
                return 1;
            }
            dump_cellcounts_optimized(outf, atac_cell_counts, samples.size());
            if (gzclose(outf) != Z_OK){
                fprintf(stderr, "ERROR: failed closing %s\n", fname.c_str());
                return 1;
            }
            fprintf(stderr, "Wrote staged ATAC counts for %lu cells\n", atac_cell_counts.size());
        }
    }

    // ================================================================
    // 2C: Species panel loading (needed before assignment)
    // ================================================================
    // NOTE: species_snpdat, species_cell_counts_rna, species_cell_counts_atac,
    // panel_meta, and panel_meta_loaded are declared earlier (before counting)
    // to support the WP3 dual counting optimization path.

    if (has_species_vcf){
        if (bam_seq2tid.empty() && !species_panel_may_be_empty){
            set<string> species_bam_chroms;
            string species_bam_header_error;
            if (!get_bam_header_chroms_and_seq2tid(
                    bamfile, species_bam_chroms, bam_seq2tid, &species_bam_header_error)){
                fprintf(stderr, "ERROR: species-panel BAM/header validation failed: %s\n",
                    species_bam_header_error.c_str());
                return 1;
            }
        }
        // WP3: if dual counting already loaded species VCF, validated samples,
        // performed dedup, and counted RNA, skip those steps here.
        if (!species_counted_dual){
            // Species panel counts are never cached in .counts, so this block must run
            // even when load_counts is true (main demux counts were loaded from disk).
            print_elapsed(start_time, "Loading species panel VCF...");
            vector<string> species_samples;  // Sample names from species VCF
            if (species_shared_vcf_name.length() > 0){
                if (!attach_shared_vcf(species_shared_vcf_name, species_snpdat, species_samples)){
                    fprintf(stderr, "ERROR: Could not attach to shared species VCF: %s\n",
                        species_shared_vcf_name.c_str());
                    return 1;
                }
                if (!bam_seq2tid.empty() &&
                    !validate_and_restrict_panel_to_bam(
                        species_snpdat, bam_seq2tid, "species genotype panel",
                        species_panel_may_be_empty)) return 1;
            }
            else if (!load_counts){
                // Disk-based species VCF requires the BAM reader for tid mapping,
                // which is only initialized when load_counts is false.
                read_vcf_samples(species_vcf_file, species_samples);
                set<string> sp_chroms;
                get_vcf_chroms(species_vcf_file, sp_chroms);
                set<string> sp_shared_chroms = shared_contig_names(sp_chroms, bam_seq2tid);
                if (sp_shared_chroms.empty() && !species_panel_may_be_empty){
                    fprintf(stderr,
                        "ERROR: required species genotype panel has zero shared BAM/VCF contigs\n");
                    return 1;
                }
                int nloaded = 0;
                if (!sp_shared_chroms.empty()){
                    nloaded = read_vcf_chroms_optimized(species_vcf_file, sp_shared_chroms,
                        bam_seq2tid, species_snpdat, vq);
                    if (nloaded < 0) return 1;
                }
                fprintf(stderr, "Loaded %d species panel SNPs\n", nloaded);
            }
            else{
                // load_counts + disk-based species VCF: not currently supported
                fprintf(stderr, "ERROR: Disk-based --species_vcf with cached .counts is not supported. "
                    "Use --species_shared_vcf or delete the .counts file to re-run from BAM.\n");
                return 1;
            }

            // §5.3: Sample set alignment check
            if (species_samples.size() != samples.size()){
                fprintf(stderr, "ERROR: Species VCF has %lu samples but demux VCF has %lu samples. "
                    "Sample sets must match.\n", species_samples.size(), samples.size());
                return 1;
            }
            for (size_t si = 0; si < samples.size(); ++si){
                if (species_samples[si] != samples[si]){
                    fprintf(stderr, "ERROR: Species VCF sample[%lu]='%s' differs from demux VCF sample[%lu]='%s'. "
                        "Sample sets must match in same order.\n",
                        si, species_samples[si].c_str(), si, samples[si].c_str());
                    return 1;
                }
            }
            fprintf(stderr, "Species VCF sample set matches demux VCF (%lu samples)\n", samples.size());

            // §5.5: SNP overlap dedup - remove species sites that overlap with demux VCF
            // Demux VCF retains priority at overlapping positions.
            // When load_counts is true, demux_positions is empty (demux SNP data was not
            // loaded), so dedup is skipped. This is safe: species and demux panels are
            // designed to be non-overlapping, and any rare overlap would only cause minor
            // double-counting in the species-only scoring pass.
            if ((species_mode == SpeciesPanelMode::AUGMENT || species_mode == SpeciesPanelMode::BOTH)
                && !demux_positions.empty()){
                int overlap_count = 0;
                int total_species_sites = 0;
                for (auto& sp_chrom : species_snpdat){
                    int tid = sp_chrom.first;
                    vector<SNPData> filtered_snps;
                    for (auto& sp_snp : sp_chrom.second.snps){
                        total_species_sites++;
                        if (demux_positions.count(make_pair(tid, sp_snp.pos)) > 0){
                            overlap_count++;
                        }
                        else{
                            filtered_snps.push_back(sp_snp);
                        }
                    }
                    sp_chrom.second.snps = filtered_snps;
                }
                if (overlap_count > 0){
                    fprintf(stderr, "Species/demux overlap: %d of %d species sites removed "
                        "(demux VCF retains priority)\n", overlap_count, total_species_sites);
                }
                else{
                    fprintf(stderr, "Species/demux overlap: 0 of %d species sites overlap\n",
                        total_species_sites);
                }
            }
            if (!bam_seq2tid.empty() &&
                !validate_and_restrict_panel_to_bam(
                    species_snpdat, bam_seq2tid, "species genotype panel",
                    species_panel_may_be_empty)) return 1;

            // Load panel metadata if needed
            if (panel_metadata_file.length() > 0){
                panel_meta = load_panel_metadata(panel_metadata_file, samples);
                panel_meta_loaded = true;
            }

            // Count reads at species-panel sites.
            // Always needed for augment/both (LLR augmentation).
            // Also needed for filter mode when --species_counts_output requests persisting counts.
            if (species_mode == SpeciesPanelMode::AUGMENT || species_mode == SpeciesPanelMode::BOTH
                || species_counts_output){
                NativeSpeciesTargetTable species_native_targets;
                precompute_all_genotypes(species_snpdat, samples.size());
                if (species_counts_output || species_assignment_output){
                    species_samples_native = panel_meta.species_list;
                    precompute_native_species_targets(
                        species_snpdat, panel_meta, samples.size(),
                        species_native_targets);
                }
                print_elapsed(start_time, "Counting alleles at species panel sites (RNA)...");
                if (n_threads > 1){
                    robin_hood::unordered_map<unsigned long, AlignedCellCounts> sp_parallel;
                    robin_hood::unordered_map<unsigned long, AlignedCellCounts> sp_native_parallel;
                    if (!count_alleles_parallel(bamfile, species_snpdat, sp_parallel,
                            cell_barcodes, samples.size(), n_threads, htslib_threads,
                            false, "",
                            accepted_weighted_conditional ? &accepted_site_weights_species : nullptr,
                            (species_counts_output || species_assignment_output)
                                ? &species_native_targets : nullptr,
                            (species_counts_output || species_assignment_output)
                                ? &sp_native_parallel : nullptr,
                            (species_counts_output || species_assignment_output)
                                ? (int)species_samples_native.size() : 0)){
                        fprintf(stderr, "ERROR: species-panel RNA allele counting failed\n");
                        return 1;
                    }
                    finalize_parallel_counts(sp_parallel, species_cell_counts_rna);
                    if (species_counts_output || species_assignment_output){
                        finalize_parallel_counts(
                            sp_native_parallel, species_cell_counts_native);
                    }
                }
                else{
                    map<pair<int, int>, map<int, float>> sp_cond, sp_tots;
                    if (!count_alleles_single_threaded(bamfile, species_snpdat, species_cell_counts_rna,
                            cell_barcodes, samples.size(), sp_cond, sp_tots, false,
                            (species_counts_output || species_assignment_output)
                                ? &species_native_targets : nullptr,
                            (species_counts_output || species_assignment_output)
                                ? &species_cell_counts_native : nullptr,
                            (species_counts_output || species_assignment_output)
                                ? (int)species_samples_native.size() : 0)){
                        fprintf(stderr, "ERROR: species-panel RNA allele counting failed\n");
                        return 1;
                    }
                }
                fprintf(stderr, "Species RNA counts for %lu cells\n", species_cell_counts_rna.size());

                if (atac_mode){
                    print_elapsed(start_time, "Counting alleles at species panel sites (ATAC)...");
                    if (n_threads > 1){
                        robin_hood::unordered_map<unsigned long, AlignedCellCounts> sp_atac_parallel;
                        if (!count_alleles_parallel(atac_bamfile, species_snpdat, sp_atac_parallel,
                                cell_barcodes, samples.size(), n_threads, htslib_threads)){
                            fprintf(stderr, "ERROR: species-panel ATAC allele counting failed\n");
                            return 1;
                        }
                        finalize_parallel_counts(sp_atac_parallel, species_cell_counts_atac);
                    }
                    else{
                        map<pair<int, int>, map<int, float>> sp_cond2, sp_tots2;
                        if (!count_alleles_single_threaded(atac_bamfile, species_snpdat, species_cell_counts_atac,
                                cell_barcodes, samples.size(), sp_cond2, sp_tots2, false)){
                            fprintf(stderr, "ERROR: species-panel ATAC allele counting failed\n");
                            return 1;
                        }
                    }
                }
            }
        }
        else{
            // Dual counting path already handled VCF loading, validation, dedup, 
            // and RNA counting. Only ATAC counting remains (if needed).
            fprintf(stderr, "Species RNA counting already done via dual-panel path\n");
            
            if (atac_mode && (species_mode == SpeciesPanelMode::AUGMENT
                || species_mode == SpeciesPanelMode::BOTH || species_counts_output)){
                print_elapsed(start_time, "Counting alleles at species panel sites (ATAC)...");
                if (n_threads > 1){
                    robin_hood::unordered_map<unsigned long, AlignedCellCounts> sp_atac_parallel;
                    if (!count_alleles_parallel(atac_bamfile, species_snpdat, sp_atac_parallel,
                            cell_barcodes, samples.size(), n_threads, htslib_threads)){
                        fprintf(stderr, "ERROR: species-panel ATAC allele counting failed\n");
                        return 1;
                    }
                    finalize_parallel_counts(sp_atac_parallel, species_cell_counts_atac);
                }
                else{
                    map<pair<int, int>, map<int, float>> sp_cond2, sp_tots2;
                    if (!count_alleles_single_threaded(atac_bamfile, species_snpdat, species_cell_counts_atac,
                            cell_barcodes, samples.size(), sp_cond2, sp_tots2, false)){
                        fprintf(stderr, "ERROR: species-panel ATAC allele counting failed\n");
                        return 1;
                    }
                }
            }
        }
    }

    // ================================================================
    // Write native species counts / condf / samples to disk (if requested)
    // ================================================================
    if (species_counts_output || species_assignment_output){
        load_panel_metadata_if_needed(
            panel_metadata_file, samples, panel_meta, panel_meta_loaded);
        species_samples_native = panel_meta.species_list;
        if (species_cell_counts_native.empty()){
            fprintf(stderr,
                "ERROR: native species output requested but no per-site normalized native species counts were produced\n");
            return 1;
        }
    }

    if (species_counts_output && !species_cell_counts_native.empty()){
        {
            string fname = output_prefix + ".species_counts";
            gzFile outf = gzopen(fname.c_str(), "w");
            if (!outf){
                fprintf(stderr, "ERROR: could not open %s for writing: %s\n",
                    fname.c_str(), strerror(errno));
                return 1;
            }
            dump_cellcounts_optimized(outf, species_cell_counts_native,
                species_samples_native.size());
            if (gzclose(outf) != Z_OK){
                fprintf(stderr, "ERROR: failed closing %s\n", fname.c_str());
                return 1;
            }
            fprintf(stderr, "Wrote staged native species counts for %lu cells (%lu species columns)\n",
                species_cell_counts_native.size(), species_samples_native.size());
        }

        {
            map<pair<int, int>, map<int, float>> sp_condf;
            species_condf_stats = compute_species_condf_native(
                species_snpdat, sp_condf, panel_meta, samples.size(),
                accepted_weighted_conditional ? &accepted_site_weights_species : nullptr);
            string fname = output_prefix + ".species_condf";
            FILE* outf = fopen(fname.c_str(), "w");
            if (!outf){
                fprintf(stderr, "ERROR: could not open %s for writing: %s\n",
                    fname.c_str(), strerror(errno));
                return 1;
            }
            dump_exp_fracs(outf, sp_condf);
            if (fclose(outf) != 0){
                fprintf(stderr, "ERROR: failed closing %s: %s\n",
                    fname.c_str(), strerror(errno));
                return 1;
            }
            fprintf(stderr, "Wrote staged native species conditional match fracs\n");

            string basis_name = output_prefix + ".species_condf_basis.tsv";
            FILE* basis_out = fopen(basis_name.c_str(), "w");
            if (!basis_out){
                fprintf(stderr, "ERROR: could not open %s for writing\n", basis_name.c_str());
                return 1;
            }
            fprintf(basis_out, "contract_version\tcondf_basis_V1_R1\n");
            fprintf(basis_out, "panel\tspecies\n");
            fprintf(basis_out, "weighting\t%s\n",
                accepted_weighted_conditional ? "accepted_observation_weighted" : "vcf_site_weighted");
            fprintf(basis_out, "truth_inputs_used\tfalse\n");
            fprintf(basis_out, "observed_sites\t%llu\n",
                (unsigned long long)species_condf_stats.observed_sites);
            fprintf(basis_out, "accepted_weight\t%.17g\n", species_condf_stats.accepted_weight);
            if (fclose(basis_out) != 0){
                fprintf(stderr, "ERROR: failed closing %s: %s\n",
                    basis_name.c_str(), strerror(errno));
                return 1;
            }
        }

        {
            string fname = output_prefix + ".species_samples";
            if (!write_samples_checked(fname, species_samples_native)) return 1;
            fprintf(stderr, "Wrote native species sample list to %s\n", fname.c_str());
        }
    }

    // ================================================================
    // Step 0a: --skip_assignment bypass
    // ================================================================
    if (skip_assignment){
        print_elapsed(start_time, "Skipping identity assignment (--skip_assignment set)");
        if (!publish_outputs(output_transaction, "counts-only")) return 1;
        fprintf(stderr, "Published %s.counts; exiting (no assignment performed).\n",
                final_output_prefix.c_str());
        return 0;
    }

    // Identity priors are deliberately unavailable in the active maximin
    // scoring path. The command-line parser rejects --identity_prior above.
    const double z_doublet_prior = 0.0;

    // Assign identities
    print_elapsed(start_time, "Starting identity assignment (round 1)...");
    robin_hood::unordered_map<unsigned long, int> assn;
    robin_hood::unordered_map<unsigned long, double> assn_llr;
    
    fprintf(stderr, "Finding likeliest identities of cells...\n");
    map<int, double> prior_weights;
    
    // First round of assignments
    // NOTE: allowed_ids contains all singlets + doublets (including singlet components)
    //       allowed_ids2 contains ONLY what was in the original ID file
    // When -I is used with doublet combinations, allowed_ids2 has only doublets,
    // which causes singlets to be disallowed in the final assignment step.
    if (!assign_ids_parallel(cell_counts, samples, assn, assn_llr,
            allowed_ids, allowed_ids2, doublet_rate, error_ref, error_alt,
            false, prior_weights, n_threads, n_target)) {
        fprintf(stderr, "ERROR: identity scoring failed in round 1\n");
        return 1;
    }
    
    // Re-estimate error rates
    print_elapsed(start_time, "Estimating error rates...");
    fprintf(stderr, "Finding likeliest error rates...\n");
    pair<double, double> err_new = infer_error_rates_optimized(cell_counts, samples.size(),
        assn, assn_llr, error_ref, error_alt, error_sigma);
    
    double error_ref_posterior = err_new.first;
    double error_alt_posterior = err_new.second; 
    
    fprintf(stderr, "Posterior error rates:\n");
    fprintf(stderr, "\tref mismatch: %f\n", error_ref_posterior);
    fprintf(stderr, "\talt mismatch: %f\n", error_alt_posterior);
    
    // Re-assign with posterior error rates
    print_elapsed(start_time, "Re-inferring identities (round 2)...");
    fprintf(stderr, "Re-inferring identities of cells...\n");
    if (!assign_ids_parallel(cell_counts, samples, assn, assn_llr,
            allowed_ids, allowed_ids2, doublet_rate,
            error_ref_posterior, error_alt_posterior,
            false, prior_weights, n_threads, n_target)) {
        fprintf(stderr, "ERROR: identity scoring failed in round 2\n");
        return 1;
    }

    if (idfile_doublet_given && allowed_ids.size() > allowed_ids2.size()) {
        fprintf(stderr,
            "Identity restriction: undeclared component singlets retained internally "
            "for declared fusion likelihoods and remain ineligible as final outputs.\n");
    }

    // NEW: Final assignment with diagnostic collection
    robin_hood::unordered_map<unsigned long, CellDiagnostics> cell_diagnostics;
    robin_hood::unordered_map<unsigned long, vector<RunnerUp> > cell_runner_ups;
    
    // Het data structures
    robin_hood::unordered_map<unsigned long, CellHetData> het_data;
    robin_hood::unordered_map<int, ChromSNPs> het_snpdat;
    vector<pair<int, int>> idx_to_site;  // Site index to (tid, pos) for persite method
    
    bool het_vcf_available = (het_vcf_file.length() > 0 || shared_het_vcf_name.length() > 0);
    int n_het_loaded = 0;
    
    if (het_vcf_available && write_diagnostics){
        if (load_counts && bamfile.length() == 0){
            fprintf(stderr, "WARNING: Het VCF diagnostics require a BAM file. Skipping het balance "
                            "computation in counts-only mode.\n");
            fprintf(stderr, "Re-run with --bam to include ploidy diagnostics.\n");
        }
        else{
            vector<string> het_samples;  // Not used but needed for attach_shared_vcf
            
            if (shared_het_vcf_name.length() > 0){
                // Load from shared memory
                print_elapsed(start_time, "Attaching to shared het VCF for ploidy diagnostics...");
                fprintf(stderr, "Attaching to shared het VCF: %s\n", shared_het_vcf_name.c_str());
                
                if (!attach_shared_vcf(shared_het_vcf_name, het_snpdat, het_samples)){
                    fprintf(stderr, "ERROR: could not attach to required shared het VCF: %s\n",
                        shared_het_vcf_name.c_str());
                    return 1;
                }
                else{
                    // Count SNPs loaded
                    for (auto& kv : het_snpdat){
                        n_het_loaded += kv.second.snps.size();
                    }
                    fprintf(stderr, "Attached to shared het VCF with %d sites\n", n_het_loaded);
                }
            }
            else{
                // Load from file
                print_elapsed(start_time, "Loading het VCF for ploidy diagnostics...");
                fprintf(stderr, "Loading het VCF: %s\n", het_vcf_file.c_str());
                
                map<string, int> seq2tid_het;
                
                // Get chromosome names from BAM header
                set<string> chroms_for_het;
                {
                    string header_error;
                    if (!get_bam_header_chroms_and_seq2tid(
                            bamfile, chroms_for_het, seq2tid_het, &header_error)){
                        fprintf(stderr, "ERROR: BAM/header validation for het diagnostics failed: %s\n",
                            header_error.c_str());
                        return 1;
                    }
                }
                
                n_het_loaded = load_het_vcf(het_vcf_file, chroms_for_het, seq2tid_het, het_snpdat, vq);
                if (n_het_loaded < 0){
                    return 1;
                }
                fprintf(stderr, "Loaded %d het sites\n", n_het_loaded);
            }
            
            if (n_het_loaded > 0){
                const char* method_name = (het_method == HetBalanceMethod::PERSITE) ? "per-site" : "Welford";
                fprintf(stderr, "Using %s method (min_het_sites=%d", method_name, min_het_sites);
                if (het_method == HetBalanceMethod::PERSITE) {
                    fprintf(stderr, ", min_het_depth=%.1f", min_het_depth);
                }
                fprintf(stderr, ")\n");
                print_elapsed(start_time, "Counting alleles at het sites...");
                if (!count_het_alleles_extended(bamfile, het_snpdat, het_data, idx_to_site,
                        cell_barcodes, samples.size(), n_threads, htslib_threads, het_method)){
                    fprintf(stderr, "ERROR: het/ploidy allele counting failed\n");
                    return 1;
                }
            }
            else{
                fprintf(stderr, "WARNING: No het sites available\n");
            }
        }
    }
    else if (write_diagnostics && !het_vcf_available){
        fprintf(stderr, "WARNING: --het_vcf/--shared_het_vcf not provided. Ploidy-related diagnostics "
                        "(het_balance_var, n_het_sites, het_total_depth) will not be computed.\n"
                        "For full tetraploid analysis, provide het VCF from downsample_vcf_parallel.\n");
    }

    // ================================================================
    // 3.6: ATAC het balance diagnostics
    // ================================================================
    robin_hood::unordered_map<unsigned long, CellHetData> atac_het_data;
    robin_hood::unordered_map<int, ChromSNPs> atac_het_snpdat;
    vector<pair<int, int>> atac_idx_to_site;
    bool atac_het_available = atac_mode &&
        (atac_het_vcf_file.length() > 0 || atac_shared_het_vcf_name.length() > 0);
    int n_atac_het_loaded = 0;

    if (atac_het_available && write_diagnostics && !load_counts){
        if (atac_shared_het_vcf_name.length() > 0){
            print_elapsed(start_time, "Attaching to shared ATAC het VCF...");
            vector<string> atac_het_samples;
            if (!attach_shared_vcf(atac_shared_het_vcf_name, atac_het_snpdat, atac_het_samples)){
                fprintf(stderr, "ERROR: could not attach to required shared ATAC het VCF: %s\n",
                    atac_shared_het_vcf_name.c_str());
                return 1;
            }
            else{
                for (auto& kv : atac_het_snpdat){
                    n_atac_het_loaded += kv.second.snps.size();
                }
                fprintf(stderr, "Attached to shared ATAC het VCF with %d sites\n", n_atac_het_loaded);
            }
        }
        else{
            print_elapsed(start_time, "Loading ATAC het VCF...");
            map<string, int> seq2tid_atac_het;
            set<string> chroms_atac_het;
            {
                string header_error;
                if (!get_bam_header_chroms_and_seq2tid(
                        atac_bamfile, chroms_atac_het, seq2tid_atac_het, &header_error)){
                    fprintf(stderr, "ERROR: ATAC BAM/header validation for het diagnostics failed: %s\n",
                        header_error.c_str());
                    return 1;
                }
            }
            n_atac_het_loaded = load_het_vcf(atac_het_vcf_file, chroms_atac_het,
                seq2tid_atac_het, atac_het_snpdat, vq);
            if (n_atac_het_loaded < 0){
                return 1;
            }
            fprintf(stderr, "Loaded %d ATAC het sites\n", n_atac_het_loaded);
        }

        if (n_atac_het_loaded > 0){
            print_elapsed(start_time, "Counting alleles at ATAC het sites...");
            if (!count_het_alleles_extended(atac_bamfile, atac_het_snpdat, atac_het_data,
                    atac_idx_to_site, cell_barcodes, samples.size(), n_threads,
                    htslib_threads, het_method)){
                fprintf(stderr, "ERROR: ATAC het/ploidy allele counting failed\n");
                return 1;
            }
            fprintf(stderr, "ATAC het data for %lu cells\n", atac_het_data.size());
        }
    }
    
    // Prepare extended evidence pointers for diagnostic calls
    robin_hood::unordered_map<unsigned long, CellCounts>* atac_ptr =
        atac_mode ? &atac_cell_counts : nullptr;
    const map<int, double>* prior_ptr = nullptr;
    robin_hood::unordered_map<unsigned long, CellCounts>* sp_rna_ptr =
        (species_mode == SpeciesPanelMode::AUGMENT || species_mode == SpeciesPanelMode::BOTH)
        ? &species_cell_counts_rna : nullptr;
    robin_hood::unordered_map<unsigned long, CellCounts>* sp_atac_ptr =
        (atac_mode && sp_rna_ptr != nullptr) ? &species_cell_counts_atac : nullptr;

    // Do final assignment with diagnostics if requested
    if (write_diagnostics){
        print_elapsed(start_time, "Final assignment with diagnostic collection...");
        
        if (het_data.empty()){
            // No het data - use original function
            if (!assign_ids_parallel_with_diagnostics(
                    cell_counts, samples, assn, assn_llr,
                    allowed_ids, allowed_ids2, doublet_rate,
                    error_ref_posterior, error_alt_posterior,
                    false, prior_weights, n_threads, n_target,
                    true,  // compute_diagnostics
                    n_runner_ups,
                    close_threshold,
                    NULL,  // no het counts
                    cell_diagnostics,
                    cell_runner_ups,
                    atac_ptr, prior_ptr, z_doublet_prior,
                    sp_rna_ptr, sp_atac_ptr)) {
                fprintf(stderr, "ERROR: final identity scoring with diagnostics failed\n");
                return 1;
            }
        }
        else{
            // Use extended method (Welford or per-site)
            if (!assign_ids_parallel_with_diagnostics_extended(
                    cell_counts, samples, assn, assn_llr,
                    allowed_ids, allowed_ids2, doublet_rate,
                    error_ref_posterior, error_alt_posterior,
                    false, prior_weights, n_threads, n_target,
                    true,  // compute_diagnostics
                    n_runner_ups,
                    close_threshold,
                    &het_data, &het_snpdat, &idx_to_site,
                    het_method, min_het_sites, min_het_depth,
                    cell_diagnostics,
                    cell_runner_ups,
                    atac_ptr, prior_ptr, z_doublet_prior,
                    sp_rna_ptr, sp_atac_ptr)) {
                fprintf(stderr, "ERROR: final identity scoring with het diagnostics failed\n");
                return 1;
            }
        }
    }

    // QC
    print_elapsed(start_time, "Running QC...");
    map<int, double> p_ncell;
    map<int, double> p_llr;
    id_qc(assn, assn_llr, p_ncell, p_llr);

    // Write assignments to staged outputs.
    print_elapsed(start_time, "Writing staged outputs...");
    {
        string fname = output_prefix + ".assignments";
        FILE* outf = fopen(fname.c_str(), "w");
        if (!outf){
            fprintf(stderr, "ERROR: could not open %s for writing: %s\n",
                fname.c_str(), strerror(errno));
            return 1;
        }
        fprintf(stderr, "Writing cell-individual assignments to staged output...\n");
        dump_assignments(outf, assn, assn_llr, samples, barcode_group, cellranger, seurat, underscore);
        if (fclose(outf) != 0){
            fprintf(stderr, "ERROR: failed closing %s: %s\n",
                fname.c_str(), strerror(errno));
            return 1;
        }
    }

    // The summary reports the user-visible final prefix, never the temporary name.
    {
        string fname = output_prefix + ".summary";
        FILE* outf = fopen(fname.c_str(), "w");
        if (!outf){
            fprintf(stderr, "ERROR: could not open %s for writing: %s\n",
                fname.c_str(), strerror(errno));
            return 1;
        }
        string summary_prefix = final_output_prefix;
        write_summary(outf, summary_prefix, assn, samples, error_ref,
            error_alt, error_sigma, error_ref_posterior,
            error_alt_posterior, vcf_file, vq, doublet_rate,
            p_ncell, p_llr);
        if (fclose(outf) != 0){
            fprintf(stderr, "ERROR: failed closing %s: %s\n",
                fname.c_str(), strerror(errno));
            return 1;
        }
    }
    
    // NEW: Write diagnostic files
    if (write_diagnostics){
        // Write .diagnostics.gz
        {
            string fname = output_prefix + ".diagnostics.gz";
            if (!write_diagnostics_gz(fname, samples, assn, assn_llr, cell_diagnostics,
                    samples.size(), barcode_group, cellranger, seurat, underscore,
                    atac_het_available ? &atac_het_data : NULL,
                    atac_het_available ? &atac_het_snpdat : NULL,
                    atac_het_available ? &atac_idx_to_site : NULL,
                    het_method, min_het_sites, min_het_depth,
                    dump_selection_audit)){
                return 1;
            }
        }
        
        // Write .runner_ups.gz
        {
            string fname = output_prefix + ".runner_ups.gz";
            if (!write_runner_ups_gz(fname, samples, assn, cell_diagnostics,
                    cell_runner_ups, samples.size(), barcode_group,
                    cellranger, seurat, underscore)){
                return 1;
            }
        }
        
        fprintf(stderr, "Diagnostic output complete.\n");
        if (!het_vcf_available){
            fprintf(stderr, "Note: unavailable het/ploidy values are written as NA with an explicit state.\n");
        }
    }
    
    // ================================================================
    // 2C: Native species assignment output
    // ================================================================
    // This is deliberately independent from .assignments.  .assignments remains
    // the individual-native demux call from the individual SNP panel.  The species
    // panel writes its own species-native assignment file with species labels only.
    if (species_assignment_output){
        if (species_cell_counts_native.empty()){
            fprintf(stderr,
                "ERROR: --species_assignment_output requested but native species counts are unavailable\n");
            return 1;
        }
        if (!write_native_species_assignments(output_prefix, species_cell_counts_native,
                species_samples_native, species_idfile,
                species_doublet_rate, error_ref_posterior, error_alt_posterior,
                n_threads, n_target, barcode_group, cellranger, seurat, underscore)){
            return 1;
        }
    }

    if (!publish_outputs(output_transaction, "demultiplexing")) return 1;
    print_elapsed(start_time, "Complete!");
    fprintf(stderr, "Published complete output bundle at prefix %s\n",
        final_output_prefix.c_str());

    return 0;
}
