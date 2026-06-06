// ============================================================================
// tet_contam_estimate.cpp
// Tetraploid contamination estimator with panel selection and warm start
//
// Joint estimation of allele ratio (r) and contamination (c) for heterotypic
// tetraploid cells. Diploid singlets and homotypic tetraploids use the standard
// two-component model. Outputs standard .contam_rate/.contam_prof files plus
// a .allele_ratio file with per-cell genome A fraction estimates.
//
// Replaces quant3_contam_ap with a cleaner interface: explicit panel selection
// (--interindividual / --interspecies), unified warm start (--warm_start),
// and removal of degenerate/harmful flags.
// ============================================================================

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
#include <map>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <regex>
#include <math.h>
#include <cctype>
#include <zlib.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include <htswrapper/mex.h>
#include "common.h"
#include "ambient_rna_three_ap.h"
#include "ambient_rna_gex.h"
#include "demux_vcf_io.h"

using std::cout;
using std::endl;
using namespace std;

// Version tracking
const string QC_VERSION = "3.0-tet";
const string QC_VERSION_MSG = "Tetraploid contamination estimator with panel selection and warm start";
const string TOOL_NAME = "tet_contam_estimate";


static string trim_copy(const string& x){
    size_t a = x.find_first_not_of(" \t\r\n");
    if (a == string::npos) return "";
    size_t b = x.find_last_not_of(" \t\r\n");
    return x.substr(a, b - a + 1);
}

static bool is_chinobo_or_hybrid_label(const string& x){
    string y = x;
    transform(y.begin(), y.end(), y.begin(), ::tolower);
    return y == "hy" || y == "hybrid" || y.find("chinobo") != string::npos;
}

static vector<string> split_plus_tokens(const string& label){
    vector<string> out;
    string cur;
    stringstream ss(label);
    while (getline(ss, cur, '+')){
        cur = trim_copy(cur);
        if (!cur.empty()) out.push_back(cur);
    }
    return out;
}

static void normalize_species_composition(map<string, double>& comp){
    double s = 0.0;
    for (auto it = comp.begin(); it != comp.end(); ++it){
        if (it->second > 0.0) s += it->second;
    }
    if (s <= 0.0) return;
    for (auto it = comp.begin(); it != comp.end(); ++it){
        it->second /= s;
    }
}

static map<string, double> token_species_composition(
    const string& token,
    const vector<string>& indiv_samples,
    const map<string, int>& indiv_name_to_idx,
    const PanelMetadata& panel,
    const set<string>& native_species){

    map<string, double> comp;
    string tok = trim_copy(token);
    if (tok.empty()) return comp;

    if (native_species.count(tok) > 0){
        comp[tok] = 1.0;
        return comp;
    }

    if (is_chinobo_or_hybrid_label(tok)){
        comp["B"] += 0.5;
        comp["C"] += 0.5;
        return comp;
    }

    auto ni = indiv_name_to_idx.find(tok);
    if (ni != indiv_name_to_idx.end()){
        int idx = ni->second;
        for (const string& sp : panel.species_list){
            auto sit = panel.species_to_sample_indices.find(sp);
            if (sit == panel.species_to_sample_indices.end()) continue;
            for (int member_idx : sit->second){
                if (member_idx == idx){
                    comp[sp] += panel.get_weight(sp, idx);
                }
            }
        }
        if (!comp.empty()){
            normalize_species_composition(comp);
            return comp;
        }
    }

    auto sp = panel.indiv_to_species.find(tok);
    if (sp != panel.indiv_to_species.end()){
        if (native_species.count(sp->second) > 0){
            comp[sp->second] = 1.0;
        } else if (is_chinobo_or_hybrid_label(sp->second)){
            comp["B"] += 0.5;
            comp["C"] += 0.5;
        }
    }
    normalize_species_composition(comp);
    return comp;
}

static map<string, double> identity_species_composition(
    const string& identity,
    const vector<string>& indiv_samples,
    const map<string, int>& indiv_name_to_idx,
    const PanelMetadata& panel,
    const set<string>& native_species){

    vector<string> toks = split_plus_tokens(identity);
    map<string, double> comp;
    if (toks.empty()) return comp;

    double component_weight = 1.0 / (double)toks.size();
    for (const string& tok : toks){
        map<string, double> tc = token_species_composition(tok, indiv_samples,
            indiv_name_to_idx, panel, native_species);
        for (auto it = tc.begin(); it != tc.end(); ++it){
            comp[it->first] += component_weight * it->second;
        }
    }
    normalize_species_composition(comp);
    return comp;
}

static bool identity_contains_hybrid_component(const string& identity){
    vector<string> toks = split_plus_tokens(identity);
    for (const string& tok : toks){
        if (is_chinobo_or_hybrid_label(tok)) return true;
    }
    return false;
}

static string assignment_label_from_idx(int idx, const vector<string>& samples){
    if (idx >= 0 && idx < (int)samples.size()) return samples[idx];
    if (idx >= (int)samples.size()){
        pair<int, int> c = idx_to_hap_comb(idx, samples.size());
        if (c.first >= 0 && c.first < (int)samples.size() &&
            c.second >= 0 && c.second < (int)samples.size()){
            return samples[c.first] + "+" + samples[c.second];
        }
    }
    return "";
}

static map<unsigned long, map<int, double> > build_weighted_species_composition_overrides(
    const string& individual_assignments_file,
    const string& individual_samples_file,
    const string& panel_metadata_file,
    const vector<string>& species_samples){

    map<unsigned long, map<int, double> > overrides;

    fprintf(stderr, "Weighted species-composition override loader: START\n");
    fprintf(stderr, "  individual assignments: %s\n",
        individual_assignments_file.empty() ? "<empty>" : individual_assignments_file.c_str());
    fprintf(stderr, "  individual samples:     %s\n",
        individual_samples_file.empty() ? "<empty>" : individual_samples_file.c_str());
    fprintf(stderr, "  panel metadata:         %s\n",
        panel_metadata_file.empty() ? "<empty>" : panel_metadata_file.c_str());
    fprintf(stderr, "  species sample universe (%lu):", species_samples.size());
    for (int i = 0; i < (int)species_samples.size(); ++i){
        fprintf(stderr, " %s", species_samples[i].c_str());
    }
    fprintf(stderr, "\n");

    if (individual_assignments_file.empty() || individual_samples_file.empty() ||
        panel_metadata_file.empty()){
        fprintf(stderr, "WARNING: weighted species-composition overrides skipped because one or more required paths are empty.\n");
        return overrides;
    }

    bool have_indiv_assignments = file_exists(individual_assignments_file);
    bool have_indiv_samples = file_exists(individual_samples_file);
    bool have_panel_metadata = file_exists(panel_metadata_file);
    fprintf(stderr, "  file_exists(assignments)=%s file_exists(samples)=%s file_exists(panel_metadata)=%s\n",
        have_indiv_assignments ? "YES" : "NO",
        have_indiv_samples ? "YES" : "NO",
        have_panel_metadata ? "YES" : "NO");

    if (!have_indiv_assignments || !have_indiv_samples || !have_panel_metadata){
        fprintf(stderr, "WARNING: weighted species-composition overrides skipped because one or more required files are missing.\n");
        return overrides;
    }

    string individual_samples_path = individual_samples_file;
    vector<string> indiv_samples;
    load_samples(individual_samples_path, indiv_samples);
    fprintf(stderr, "  loaded %lu individual samples for weighted composition lookup\n",
        indiv_samples.size());
    if (indiv_samples.empty()){
        fprintf(stderr, "WARNING: weighted species-composition overrides skipped because individual samples file loaded zero samples.\n");
        return overrides;
    }

    map<string, int> indiv_name_to_idx;
    for (int i = 0; i < (int)indiv_samples.size(); ++i){
        indiv_name_to_idx[indiv_samples[i]] = i;
    }

    PanelMetadata panel = load_panel_metadata(panel_metadata_file, indiv_samples, true);
    fprintf(stderr, "  loaded panel metadata with %lu species labels\n",
        panel.species_list.size());

    set<string> native_species(species_samples.begin(), species_samples.end());
    map<string, int> species_name_to_idx;
    for (int i = 0; i < (int)species_samples.size(); ++i){
        species_name_to_idx[species_samples[i]] = i;
    }

    robin_hood::unordered_map<unsigned long, int> indiv_assn;
    robin_hood::unordered_map<unsigned long, double> indiv_llr;
    string individual_assignments_path = individual_assignments_file;
    load_assignments_from_file(individual_assignments_path, indiv_assn, indiv_llr, indiv_samples);
    fprintf(stderr, "  loaded %lu original individual assignments for weighted composition lookup\n",
        indiv_assn.size());
    if (indiv_assn.empty()){
        fprintf(stderr, "WARNING: weighted species-composition overrides skipped because original individual assignments loaded zero rows.\n");
        return overrides;
    }

    int n_loaded = 0;
    int n_weighted_multi = 0;
    int n_seen_hybrid = 0;
    int n_unmapped_hybrid = 0;
    for (auto it = indiv_assn.begin(); it != indiv_assn.end(); ++it){
        string label = assignment_label_from_idx(it->second, indiv_samples);
        if (label.empty()) continue;
        // Preserve the existing two-component/r-feedback model for ordinary
        // species fusions.  Override only identities involving the real
        // bonobo-chimp hybrid, where no single native pair identity correctly
        // represents all species components and weights.
        if (!identity_contains_hybrid_component(label)) continue;
        n_seen_hybrid++;

        map<string, double> sp_comp = identity_species_composition(label, indiv_samples,
            indiv_name_to_idx, panel, native_species);
        if (sp_comp.empty()){
            n_unmapped_hybrid++;
            if (n_unmapped_hybrid <= 10){
                fprintf(stderr, "WARNING: could not derive species composition for hybrid-containing identity: %s\n",
                    label.c_str());
            }
            continue;
        }

        map<int, double> idx_comp;
        for (auto ci = sp_comp.begin(); ci != sp_comp.end(); ++ci){
            auto si = species_name_to_idx.find(ci->first);
            if (si != species_name_to_idx.end() && ci->second > 0.0){
                idx_comp[si->second] += ci->second;
            }
        }
        double s = 0.0;
        for (auto ci = idx_comp.begin(); ci != idx_comp.end(); ++ci) s += ci->second;
        if (s <= 0.0) continue;
        for (auto ci = idx_comp.begin(); ci != idx_comp.end(); ++ci) ci->second /= s;

        overrides[it->first] = idx_comp;
        n_loaded++;
        if (idx_comp.size() > 1) n_weighted_multi++;
    }

    fprintf(stderr, "Weighted species-composition overrides loaded: %d cells from %s (%d multi-species compositions)\n",
        n_loaded, individual_assignments_file.c_str(), n_weighted_multi);
    fprintf(stderr, "  Hybrid-containing original assignments seen: %d; unmapped hybrid identities: %d\n",
        n_seen_hybrid, n_unmapped_hybrid);
    fprintf(stderr, "  Composition model: each diploid identity maps to native species weights; fusions average component vectors.\n");
    fprintf(stderr, "  Example: Chinobo-mCherry -> B:0.5,C:0.5; Chinobo-mCherry+JOS3C1 -> B:0.25,C:0.25,O:0.5\n");
    if (n_loaded == 0){
        fprintf(stderr, "WARNING: weighted species-composition loader completed but produced zero usable overrides.\n");
    }

    return overrides;
}

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "%s v%s\n", TOOL_NAME.c_str(), QC_VERSION.c_str());
    fprintf(stderr, "%s\n\n", QC_VERSION_MSG.c_str());
    fprintf(stderr, "Usage: %s [OPTIONS]\n\n", TOOL_NAME.c_str());
    fprintf(stderr, "Estimates per-cell ambient RNA contamination rate (c) and allele ratio (r)\n");
    fprintf(stderr, "for tetraploid cell fusions, given output of a demux_vcf or demux_parallel run.\n\n");

    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --output_prefix -o  Output prefix. Interindividual mode finds .counts/.condf/.assignments/.samples;\n");
    fprintf(stderr, "                        interspecies mode finds .species_counts/.species_condf/.species_assignments/.species_samples.\n\n");

    fprintf(stderr, "===== SNP PANEL (exactly one required, mutually exclusive) =====\n");
    fprintf(stderr, "    --interindividual   Native individual path: use .counts/.condf/.assignments/.samples\n");
    fprintf(stderr, "    --interspecies      Native species path: use .species_counts/.species_condf/.species_assignments/.species_samples.\n");
    fprintf(stderr, "                        With --panel_metadata plus .assignments/.samples, derives weighted native species composition\n");
    fprintf(stderr, "                        from original individual calls (e.g. Chinobo+O -> B:0.25,C:0.25,O:0.5).\n");
    fprintf(stderr, "    --condf FILE        Override .condf file path (interindividual mode)\n");
    fprintf(stderr, "    --species_condf FILE  Override .species_condf file path (interspecies mode)\n\n");

    fprintf(stderr, "===== AMBIENT PROFILE INITIALIZATION (at most one) =====\n");
    fprintf(stderr, "    --warm_start -W FILE   Load ambient profile as solver starting point.\n");
    fprintf(stderr, "                           FILE is .contam_prof_empty or .species_prof_empty.\n");
    fprintf(stderr, "                           Profile is freely refined during estimation.\n");
    fprintf(stderr, "    --fixed_ambient -A FILE  Lock ambient profile (no estimation). Use with caution.\n\n");

    fprintf(stderr, "===== TETRAPLOID-AWARE =====\n");
    fprintf(stderr, "    --expected_lines -X FILE  Expected lines file (locks combo identities)\n\n");

    fprintf(stderr, "===== SOLVER OPTIONS =====\n");
    fprintf(stderr, "    --loo               Leave-one-out ambient profiles\n");
    fprintf(stderr, "    --r_feedback        Feed per-cell allele ratios into ambient profile estimation\n");
    fprintf(stderr, "    --adaptive_prior    Auto-detect pathological distributions, apply fixed prior fallback\n");
    fprintf(stderr, "    --run_once -r       Single pass, no iterative convergence\n");
    fprintf(stderr, "    --num_threads -T    Parallel threads (default: 1)\n");
    fprintf(stderr, "    --bootstrap -b      Bootstrap replicates (default: 100)\n");
    fprintf(stderr, "    --thorough_multistart  Always run r=0.2/0.8 alternate BFGS starts for\n");
    fprintf(stderr, "                           heterotypic cells (slower, maximum robustness).\n");
    fprintf(stderr, "    --no_adaptive_multistart  Disable boundary-triggered alternate BFGS starts.\n\n");

    fprintf(stderr, "===== SPECIES REGULARIZATION (Mode 1/3 enhancement) =====\n");
    fprintf(stderr, "    --species_regularize  LEGACY/PARKED mixed-resolution regularizer. Not native V1_R3 joint evidence.\n");
    fprintf(stderr, "                          Requires --allow_legacy_species_regularize and --panel_metadata. Incompatible with --interspecies.\n");
    fprintf(stderr, "    --allow_legacy_species_regularize  Permit the parked legacy regularizer intentionally.\n");
    fprintf(stderr, "    --panel_metadata -P FILE  TSV mapping individual to species\n");
    fprintf(stderr, "    --species_counts FILE  Species-diagnostic counts for species-level solver signal.\n");
    fprintf(stderr, "                           Requires companion .species_condf file alongside.\n\n");

    fprintf(stderr, "===== STANDARD OPTIONS =====\n");
    fprintf(stderr, "    --error_ref -e      Ref->alt error rate (default: 0.001)\n");
    fprintf(stderr, "    --error_alt -E      Alt->ref error rate (default: 0.001)\n");
    fprintf(stderr, "    --ids -i            Filtered individual list (one name per line)\n");
    fprintf(stderr, "    --ids_doublet -I    Filtered doublet identity list\n");
    fprintf(stderr, "    --doublet_rate -D   Expected doublet rate (default: no expectation)\n");
    fprintf(stderr, "    --no_weights -w     Disable LLR-based weighting\n");
    fprintf(stderr, "    --libname -n        Library name for barcode formatting\n");
    fprintf(stderr, "    --cellranger -C     CellRanger barcode format (-1 suffix)\n");
    fprintf(stderr, "    --seurat -S         Seurat barcode format (libname_ prefix)\n");
    fprintf(stderr, "    --underscore -U     Underscore separator for libname\n");
    fprintf(stderr, "    --llr -l            LLR filter cutoff for assignments (default: 0)\n\n");

    fprintf(stderr, "===== GEX DECONTAMINATION =====\n");
    fprintf(stderr, "    --barcodes -B       Barcodes file (MEX format, optionally gzipped)\n");
    fprintf(stderr, "    --features -F       Features file (MEX format, optionally gzipped)\n");
    fprintf(stderr, "    --matrix -M         Matrix file (MEX format, optionally gzipped)\n");
    fprintf(stderr, "    --feature_type -t   Feature type filter (e.g. \"Gene Expression\")\n");
    fprintf(stderr, "    --clusts -c         Cell-cluster assignments file\n");
    fprintf(stderr, "    --noround -R        Output unrounded decimal counts\n");
    fprintf(stderr, "    --skip_genes -g     File of gene names to exclude\n");
    fprintf(stderr, "    --skip_genes_regex -G  Pipe-separated regex for genes to exclude (default: \"^MT-\")\n\n");

    fprintf(stderr, "    --help -h           Display this message and exit.\n");
    exit(code);
}

int parse_clustfile(string& filename,
    robin_hood::unordered_map<unsigned long, int>& clusts,
    vector<string>& clustnames){

    ifstream inf(filename.c_str());
    string bcstr;
    string clustn;
    set<string> clustnames_set;
    robin_hood::unordered_map<unsigned long, string> ctmp;
    while (inf >> bcstr >> clustn){
        clustnames_set.insert(clustn);
        unsigned long ul = bc_ul(bcstr);
        ctmp.emplace(ul, clustn);
    }

    map<string, int> clust2idx;
    int idx = 0;
    for (set<string>::iterator s = clustnames_set.begin(); s != clustnames_set.end(); ++s){
        clust2idx.insert(make_pair(*s, idx));
        clustnames.push_back(*s);
        idx++;
    }
    for (robin_hood::unordered_map<unsigned long, string>::iterator t = ctmp.begin();
        t != ctmp.end(); ++t){
        clusts.emplace(t->first, clust2idx[t->second]);
    }
    return clustnames.size();
}

static void intersect_int_sets(set<int>& target, const set<int>& restriction){
    set<int> keep;
    for (set<int>::iterator it = target.begin(); it != target.end(); ++it){
        if (restriction.count(*it) > 0){
            keep.insert(*it);
        }
    }
    target = keep;
}

static int count_singlet_ids(const set<int>& ids, int n_samples){
    int n = 0;
    for (set<int>::iterator it = ids.begin(); it != ids.end(); ++it){
        if (*it >= 0 && *it < n_samples){
            n++;
        }
    }
    return n;
}

enum class NativePanelMode { INTERINDIVIDUAL_NATIVE, INTERSPECIES_NATIVE };

static const char* panel_mode_name(NativePanelMode mode){
    return mode == NativePanelMode::INTERSPECIES_NATIVE
        ? "INTERSPECIES_NATIVE" : "INTERINDIVIDUAL_NATIVE";
}

static void build_species_prof_from_native_contam(
    const map<int, double>& contam_prof,
    const map<int, double>& contam_prof_conc,
    const vector<string>& samples,
    map<string, double>& species_prof,
    map<string, double>& species_prof_conc){

    species_prof.clear();
    species_prof_conc.clear();
    for (const auto& kv : contam_prof){
        if (kv.first >= 0 && kv.first < (int)samples.size()){
            species_prof[samples[kv.first]] = kv.second;
        }
    }
    for (const auto& kv : contam_prof_conc){
        if (kv.first >= 0 && kv.first < (int)samples.size()){
            species_prof_conc[samples[kv.first]] = kv.second;
        }
    }
}

static int max_index_in_counts_file(const string& filename){
    gzFile gz = gzopen(filename.c_str(), "r");
    if (!gz){
        fprintf(stderr, "ERROR: could not open counts file for dimension guard: %s\n", filename.c_str());
        exit(1);
    }
    char buf[1 << 16];
    int max_idx = -1;
    while (gzgets(gz, buf, sizeof(buf)) != Z_NULL){
        unsigned long cell = 0;
        int i = -1, g1 = -1, j = -1, g2 = -1;
        float ref = 0.0f, alt = 0.0f;
        int n = sscanf(buf, "%lu\t%d\t%d\t%d\t%d\t%f\t%f", &cell, &i, &g1, &j, &g2, &ref, &alt);
        if (n != 7){
            n = sscanf(buf, "%lu %d %d %d %d %f %f", &cell, &i, &g1, &j, &g2, &ref, &alt);
        }
        if (n != 7) continue;
        if (i > max_idx) max_idx = i;
        if (j > max_idx) max_idx = j;
    }
    gzclose(gz);
    return max_idx;
}

static int max_index_in_condf_file(const string& filename){
    ifstream infile(filename.c_str());
    if (!infile.good()){
        fprintf(stderr, "ERROR: could not open condf file for dimension guard: %s\n", filename.c_str());
        exit(1);
    }
    int idx1 = -1, type = -1, idx2 = -1;
    float frac = 0.0f;
    int max_idx = -1;
    while (infile >> idx1 >> type >> idx2 >> frac){
        if (idx1 > max_idx) max_idx = idx1;
        if (idx2 > max_idx) max_idx = idx2;
    }
    return max_idx;
}

static void enforce_native_dimensions(
    NativePanelMode mode,
    const string& counts_name,
    const string& condf_name,
    const vector<string>& samples){

    int n = (int)samples.size();
    int max_counts = max_index_in_counts_file(counts_name);
    int max_condf = max_index_in_condf_file(condf_name);
    fprintf(stderr, "Dimension guard [%s]: n_samples=%d max_counts_index=%d max_condf_index=%d\n",
        panel_mode_name(mode), n, max_counts, max_condf);
    if (max_counts >= n){
        fprintf(stderr, "ERROR: %s has index %d but %s loaded only %d samples.\n",
            counts_name.c_str(), max_counts, panel_mode_name(mode), n);
        fprintf(stderr, "       This indicates cross-resolution input, e.g. an individual-shaped species file.\n");
        exit(1);
    }
    if (max_condf >= n){
        fprintf(stderr, "ERROR: %s has index %d but %s loaded only %d samples.\n",
            condf_name.c_str(), max_condf, panel_mode_name(mode), n);
        fprintf(stderr, "       This indicates cross-resolution input, e.g. an individual-shaped species condf.\n");
        exit(1);
    }
}

void infer_from_genotypes(string& output_prefix,
    robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    string& idfile,
    bool idfile_given,
    string& idfile_doublet,
    bool idfile_doublet_given,
    vector<string>& samples,
    map<int, double>& contam_prof,
    robin_hood::unordered_map<unsigned long, double>& contam_rate,
    double doublet_rate,
    int num_threads,
    double error_ref,
    double error_alt,
    bool weight,
    bool run_once,
    int bootstrap,
    string& libname,
    bool seurat,
    bool cellranger,
    bool underscore,
    // Tetraploid-aware parameters
    const set<int>& locked_identities,
    const set<int>& safe_singlets,
    bool tetraploid_aware,
    bool ids_restricted,
    const set<int>& expected_allowed_ids,
    const set<int>& expected_allowed_ids2,
    // Solver options
    bool leave_one_out,
    bool r_feedback,
    bool adaptive_prior,
    bool thorough_multistart,
    bool adaptive_multistart,
    // Warm start / fixed ambient
    const string& warm_start_file,
    const string& fixed_ambient_file,
    // Species regularization (Mode 1/3)
    bool species_regularize,
    const string& panel_metadata_file,
    const string& species_counts_file,
    // Panel selection
    bool use_interspecies,
    // Explicit condf path overrides
    const string& user_condf_file,
    const string& user_species_condf_file,
    // Optional native species composition overrides derived from original individual calls
    const map<unsigned long, map<int, double> >& species_composition_overrides){

    // Hardcoded solver parameters (no longer CLI-configurable)
    int n_mixprop_trials = 10;
    double min_signal_gap = 0.10;
    double adaptive_mean = 0.05;
    double adaptive_init_var = 0.04;
    double adaptive_thresh = 0.20;

    // Determine which counts and condf files to load
    string counts_name, condf_name;
    if (use_interspecies){
        counts_name = output_prefix + ".species_counts";
        condf_name = !user_species_condf_file.empty()
                   ? user_species_condf_file
                   : output_prefix + ".species_condf";
        fprintf(stderr, "Interspecies mode: using species-diagnostic SNP panel\n");
        fprintf(stderr, "  Counts: %s\n", counts_name.c_str());
        fprintf(stderr, "  Condf:  %s\n", condf_name.c_str());
    } else {
        counts_name = output_prefix + ".counts";
        condf_name = !user_condf_file.empty()
                   ? user_condf_file
                   : output_prefix + ".condf";
    }

    // Load conditional matching probabilities
    map<pair<int, int>, map<int, float> > exp_match_fracs;
    if (file_exists(condf_name)){
        load_exp_fracs(condf_name, exp_match_fracs);
    }
    else{
        fprintf(stderr, "ERROR: no conditional matching probability file found: %s\n",
            condf_name.c_str());
        fprintf(stderr, "Please ensure demux_parallel was run with the appropriate VCF\n"
            "and that .condf (or .species_condf) files are present.\n");
        exit(1);
    }

    // Load filtered ID list(s), if given
    set<int> allowed_ids;
    set<int> allowed_ids2;

    if (idfile_given){
        parse_idfile(idfile, samples, allowed_ids, allowed_ids2, true);
        if (allowed_ids.size() == 0){
            fprintf(stderr, "No valid individual names found in file %s; allowing "
                "all possible individuals\n", idfile.c_str());
        }
    }
    if (idfile_doublet_given){
        parse_idfile(idfile_doublet, samples, allowed_ids, allowed_ids2, false);
        if (allowed_ids.size() == 0){
            fprintf(stderr, "No valid individual names found in file %s; allowing "
                "all possible individuals\n", idfile_doublet.c_str());
        }
    }

    // Expected-lines restriction: --expected_lines defines the legal library pool.
    // Use it to keep absent individuals/species out of the ambient candidate set
    // and out of reassignment.  This is separate from the tetraploid lock/safe
    // singlet logic in main().
    if (!expected_allowed_ids.empty()){
        if (allowed_ids.empty()){
            allowed_ids = expected_allowed_ids;
        }
        else{
            intersect_int_sets(allowed_ids, expected_allowed_ids);
        }

        if (allowed_ids2.empty()){
            allowed_ids2 = expected_allowed_ids2;
        }
        else if (!expected_allowed_ids2.empty()){
            intersect_int_sets(allowed_ids2, expected_allowed_ids2);
        }

        fprintf(stderr, "Expected-lines restriction active: %lu allowed identities, "
            "%lu exact identities, %d ambient source individuals\n",
            allowed_ids.size(), allowed_ids2.size(),
            count_singlet_ids(allowed_ids, samples.size()));

        if (allowed_ids.empty()){
            fprintf(stderr, "ERROR: expected-lines restriction left no allowed identities.\n");
            exit(1);
        }
    }

    enforce_native_dimensions(use_interspecies ? NativePanelMode::INTERSPECIES_NATIVE : NativePanelMode::INTERINDIVIDUAL_NATIVE,
        counts_name, condf_name, samples);

    // Load stored allele counts
    robin_hood::unordered_map<unsigned long, map<pair<int, int>,
        map<pair<int, int>, pair<float, float> > > > indv_allelecounts;
    if (file_exists(counts_name)){
        fprintf(stderr, "Loading counts from %s...\n", counts_name.c_str());
        load_counts_from_file(indv_allelecounts, samples, counts_name, allowed_ids);
    }
    else{
        fprintf(stderr, "ERROR: no counts found: %s\n", counts_name.c_str());
        exit(1);
    }

    double llprev = 0.0;
    double delta = 999;
    double delta_thresh = 0.1;

    fprintf(stderr, "%s v%s: %s\n", TOOL_NAME.c_str(), QC_VERSION.c_str(), QC_VERSION_MSG.c_str());

    map<int, double> contam_prof_conc;
    robin_hood::unordered_map<unsigned long, double> contam_rate_se;
    robin_hood::unordered_map<unsigned long, double> allele_ratio;
    robin_hood::unordered_map<unsigned long, double> allele_ratio_se;
    map<string, double> species_contam_prof_out;
    map<string, double> species_contam_prof_conc_out;

    // Fixed ambient profile is immutable state.  Native species mode uses the
    // same integer-keyed contam profile path as interindividual mode because
    // samples == species labels in INTERSPECIES_NATIVE.
    map<int, double> fixed_indiv_prof;
    map<int, double> fixed_indiv_conc;
    bool fixed_ambient_loaded = false;
    if (!fixed_ambient_file.empty()){
        load_contam_prof(fixed_ambient_file, fixed_indiv_prof, fixed_indiv_conc,
            samples, true);
        if (fixed_indiv_prof.empty()){
            fprintf(stderr, "ERROR: fixed ambient file %s had no matching entries for current %s sample set.\n",
                fixed_ambient_file.c_str(), use_interspecies ? "species" : "individual");
            exit(1);
        }
        fixed_ambient_loaded = true;
    }

    int nits = 0;
    while (delta > delta_thresh){
        fprintf(stderr, "===== ITERATION %d =====\n", nits+1);
        contamFinder3 cf(indv_allelecounts, assn, assn_llr, exp_match_fracs, samples.size(),
            allowed_ids, allowed_ids2);
        if (!species_composition_overrides.empty()){
            cf.set_cell_composition_overrides(species_composition_overrides);
            if (nits == 0){
                fprintf(stderr, "Weighted native species-composition model active for %lu cells\n",
                    species_composition_overrides.size());
            }
        }
        cf.set_doublet_rate(doublet_rate);
        cf.set_num_threads(num_threads);
        cf.set_thorough_multistart(thorough_multistart);
        cf.set_adaptive_multistart(adaptive_multistart);
        if (run_once){
            cf.no_reassign();
        }

        // Tetraploid-aware mode
        if (tetraploid_aware){
            cf.set_tetraploid_aware(true);
            cf.set_locked_identities(locked_identities);
            cf.set_safe_singlets(safe_singlets);
            cf.set_min_signal_gap(min_signal_gap);
            cf.set_ids_restricted(ids_restricted);
            if (nits == 0){
                fprintf(stderr, "Tetraploid-aware mode: %lu locked identities, %lu safe singlets, "
                    "min_signal_gap=%.2f\n", locked_identities.size(), safe_singlets.size(),
                    min_signal_gap);
            }
        }

        // Solver options
        if (leave_one_out){
            cf.set_use_loo(true);
            if (nits == 0){
                fprintf(stderr, "Leave-one-out ambient profile enabled\n");
            }
        }

        // Ambient profile refinements: r-feedback and adaptive prior
        if (r_feedback){
            cf.set_r_feedback(true);
        }
        if (adaptive_prior){
            cf.set_adaptive_prior(true, adaptive_mean, adaptive_init_var, 0.001, adaptive_thresh);
        }

        // Initialize to whatever was the final estimate last time,
        // but only if the individual set hasn't changed between iterations.
        // reclassify_cells() can add/remove individuals from assignments,
        // causing idx2samp to differ from the previous contam_prof keys.
        if (nits > 0){
            // Build the set of singlet IDs that the new contamFinder will use
            set<int> new_singlets;
            for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
                a != assn.end(); ++a){
                if (a->second < (int)samples.size()){
                    new_singlets.insert(a->second);
                }
                else{
                    pair<int, int> combo = idx_to_hap_comb(a->second, samples.size());
                    new_singlets.insert(combo.first);
                    new_singlets.insert(combo.second);
                }
            }
            // Also include any from allowed_ids/allowed_ids2
            for (set<int>::iterator ai = allowed_ids.begin(); ai != allowed_ids.end(); ++ai){
                if (*ai < (int)samples.size()){
                    new_singlets.insert(*ai);
                }
            }
            for (set<int>::iterator ai = allowed_ids2.begin(); ai != allowed_ids2.end(); ++ai){
                if (*ai < (int)samples.size()){
                    new_singlets.insert(*ai);
                }
            }

            // Check if old contam_prof keys (excluding -1) match new singlets
            set<int> old_keys;
            for (map<int, double>::iterator cp = contam_prof.begin();
                cp != contam_prof.end(); ++cp){
                if (cp->first >= 0){
                    old_keys.insert(cp->first);
                }
            }
            if (old_keys == new_singlets){
                cf.set_init_contam_prof(contam_prof);
            }
            else{
                fprintf(stderr, "  Individual set changed (%lu -> %lu singlets); "
                    "re-initializing contamination profile\n",
                    old_keys.size(), new_singlets.size());
            }
        }
        if (nits > 0){
            double meanc = 0.0;
            double cfrac = 1.0/(double)contam_rate.size();
            for (robin_hood::unordered_map<unsigned long, double>::iterator c = contam_rate.begin();
                c != contam_rate.end(); ++c){
                meanc += cfrac * c->second;
            }
            cf.set_init_c(meanc);
        }
        // Do standard initialization
        cf.set_error_rates(error_ref, error_alt);
        cf.set_mixprop_trials(n_mixprop_trials);
        if (weight){
            cf.use_weights();
        }

        // Warm start. Native species mode uses .species_samples, so species
        // profile files load through the same index-keyed path as individual files.
        bool contam_prof_initialized_from_warm_start = false;
        if (nits == 0 && !warm_start_file.empty()){
            map<int, double> loaded_prof;
            map<int, double> loaded_conc;
            load_contam_prof(warm_start_file, loaded_prof, loaded_conc, samples, false);

            if (!loaded_prof.empty()){
                cf.set_init_contam_prof(loaded_prof);
                contam_prof_initialized_from_warm_start = true;
                fprintf(stderr, "Warm start loaded from: %s (%lu entries, %s)\n",
                    warm_start_file.c_str(), loaded_prof.size(),
                    use_interspecies ? "species-native" : "individual-native");
            } else {
                fprintf(stderr, "WARNING: warm start file %s had no matching entries for current sample set. Ignoring warm start.\n",
                    warm_start_file.c_str());
            }
        }

        // Species-level warm start fallback for --species_regularize case.
        // This fires when samples = individual names but warm_start file has
        // species labels, so load_contam_prof() found nothing.
        if (nits == 0 && !warm_start_file.empty() && species_regularize
            && !contam_prof_initialized_from_warm_start){
            map<string, double> sp_init, sp_init_conc;
            load_species_prior(warm_start_file, sp_init, sp_init_conc);
            if (!sp_init.empty()){
                cf.set_species_init(sp_init);
                fprintf(stderr, "Warm start loaded (species-level for regularizer) from: %s (%lu species)\n",
                    warm_start_file.c_str(), sp_init.size());
            }
        }

        // Fixed ambient profile.  Apply immutable fixed profile on every
        // iteration because each iteration constructs a new contamFinder3 object.
        if (fixed_ambient_loaded){
            cf.set_init_contam_prof(fixed_indiv_prof);
            cf.set_fixed_amb_prof(true);
            if (nits == 0){
                fprintf(stderr, "Fixed ambient profile loaded from: %s (%lu entries, %s)\n",
                    fixed_ambient_file.c_str(), fixed_indiv_prof.size(),
                    use_interspecies ? "species-native" : "individual-native");
            } else {
                fprintf(stderr, "Fixed ambient profile reapplied for iteration %d\n", nits + 1);
            }
        }

        // Species regularization (Mode 1/3 enhancement)
        // Do NOT call set_species_mode() for --interspecies (Mode 2/4).
        if (species_regularize){
            PanelMetadata pm = load_panel_metadata(panel_metadata_file, samples);
            cf.set_species_mode(pm);
            fprintf(stderr, "Species regularization active (%lu species: ",
                pm.species_list.size());
            for (int si = 0; si < (int)pm.species_list.size(); si++){
                if (si > 0) fprintf(stderr, ", ");
                fprintf(stderr, "%s", pm.species_list[si].c_str());
            }
            fprintf(stderr, ")\n");
        }

        // Species-diagnostic counts for species-level solver
        if (!species_counts_file.empty() && species_regularize){
            robin_hood::unordered_map<unsigned long, map<pair<int, int>,
                map<pair<int, int>, pair<float, float> > > > sp_counts;
            string sp_counts_name = species_counts_file;
            if (file_exists(sp_counts_name)){
                fprintf(stderr, "Loading species-diagnostic counts from %s...\n",
                    sp_counts_name.c_str());
                load_counts_from_file(sp_counts, samples, sp_counts_name, allowed_ids);
                fprintf(stderr, "  Loaded species counts for %lu cells\n", sp_counts.size());
            } else {
                fprintf(stderr, "ERROR: species counts file not found: %s\n",
                    sp_counts_name.c_str());
                exit(1);
            }

            // Load companion .species_condf
            string sp_condf_name = species_counts_file;
            size_t pos = sp_condf_name.rfind(".species_counts");
            if (pos != string::npos){
                sp_condf_name = sp_condf_name.substr(0, pos) + ".species_condf";
            } else {
                sp_condf_name = sp_condf_name + ".species_condf";
            }
            map<pair<int, int>, map<int, float> > sp_expfracs;
            if (file_exists(sp_condf_name)){
                load_exp_fracs(sp_condf_name, sp_expfracs);
                fprintf(stderr, "  Loaded species condf from %s\n", sp_condf_name.c_str());
            } else {
                fprintf(stderr, "ERROR: species condf file not found: %s\n",
                    sp_condf_name.c_str());
                fprintf(stderr, "  (expected alongside %s)\n", species_counts_file.c_str());
                exit(1);
            }

            cf.set_species_counts(sp_counts, sp_expfracs);
        }

        cf.fit();

        for (int i = 0; i < (int)samples.size(); ++i){
            if (cf.contam_prof.count(i) > 0){
                fprintf(stderr, "%s) %f\n", samples[i].c_str(), cf.contam_prof[i]);
            }
        }

        double ll = cf.compute_ll();

        bool last_iter_accepted = true;
        if (run_once){
            // Break out of cycle
            assn = cf.assn;
            assn_llr = cf.assn_llr;
            contam_prof = cf.contam_prof;
            contam_rate = cf.contam_rate;
            contam_rate_se = cf.contam_rate_se;
            allele_ratio = cf.allele_ratio;
            allele_ratio_se = cf.allele_ratio_se;
            species_contam_prof_out = cf.species_contam_prof;
            species_contam_prof_conc_out = cf.species_contam_prof_conc;

            delta = 0;
        }
        else{
            last_iter_accepted = false;
            if (llprev == 0 || ll > llprev){
                assn = cf.assn;
                assn_llr = cf.assn_llr;
                contam_prof = cf.contam_prof;
                contam_rate = cf.contam_rate;
                contam_rate_se = cf.contam_rate_se;
                allele_ratio = cf.allele_ratio;
                allele_ratio_se = cf.allele_ratio_se;
                species_contam_prof_out = cf.species_contam_prof;
                species_contam_prof_conc_out = cf.species_contam_prof_conc;
                last_iter_accepted = true;
            }
            else{
                fprintf(stderr, "  Rejecting iteration %d update because likelihood did not improve "
                    "(current=%f, previous_best=%f). Keeping previous accepted state.\n",
                    nits + 1, ll, llprev);
            }
            fprintf(stderr, " -- Log likelihood: %f", ll);
            if (llprev != 0){
                delta = ll - llprev;
                fprintf(stderr, " delta = %f\n", delta);
            }
            else{
                fprintf(stderr, "\n");
            }
            if (last_iter_accepted){
                llprev = ll;
            }
            nits++;
        }

        // After iteration 1, write single-pass results to .pass1.* files
        if (nits == 1){
            {
                string fname = output_prefix + ".pass1.contam_prof";
                FILE* outf = fopen(fname.c_str(), "w");
                fprintf(stderr, "Writing pass1 contamination profile...\n");
                map<int, double> empty_conc;
                dump_contam_prof(outf, contam_prof, empty_conc, samples);
                fclose(outf);
            }
            if ((species_regularize || use_interspecies) && !species_contam_prof_out.empty()){
                string fname = output_prefix + ".pass1.species_prof";
                FILE* outf = fopen(fname.c_str(), "w");
                fprintf(stderr, "Writing pass1 species-level profile...\n");
                dump_species_prof(outf, species_contam_prof_out, species_contam_prof_conc_out);
                fclose(outf);
            }
            {
                string fname = output_prefix + ".pass1.contam_rate";
                FILE* outf = fopen(fname.c_str(), "w");
                dump_contam_rates(outf, contam_rate, contam_rate_se, samples,
                    libname, cellranger, seurat, underscore);
                fclose(outf);
            }
            {
                string fname = output_prefix + ".pass1.decontam.assignments";
                FILE* outf = fopen(fname.c_str(), "w");
                dump_assignments(outf, assn, assn_llr, samples, libname,
                    cellranger, seurat, underscore);
                fclose(outf);
            }
            fprintf(stderr, "Pass1 results written to %s.pass1.*\n", output_prefix.c_str());
        }
        if (delta <= delta_thresh && bootstrap > 0){
            if (fixed_ambient_loaded){
                // Fixed ambient proportions are externally supplied.  Bootstrapping
                // them as if they were MLE mixture proportions can produce invalid
                // zero/NaN starts and is not meaningful.  Preserve concentrations
                // from the fixed input file when available and skip bootstrap.
                fprintf(stderr, "Skipping bootstrap for --fixed_ambient; using fixed-profile "
                    "proportions/concentrations as supplied.\n");
                if (!fixed_indiv_conc.empty()){
                    contam_prof_conc = fixed_indiv_conc;
                }
            }
            else if (!last_iter_accepted){
                fprintf(stderr, "Skipping bootstrap because final iteration was rejected. "
                    "Bootstrap must run only on an accepted state.\n");
            }
            else{
                // Do bootstrapping using the current accepted iteration's contamFinder.
                cf.assn = assn;
                cf.assn_llr = assn_llr;
                cf.contam_rate = contam_rate;
                cf.contam_rate_se = contam_rate_se;
                fprintf(stderr, "Computing Dirichlet concentration parameters "
                    "on mixture proportions...\n");
                cf.bootstrap_amb_prof(bootstrap, contam_prof_conc);
                if (species_regularize && !cf.species_contam_prof_conc.empty()){
                    species_contam_prof_conc_out = cf.species_contam_prof_conc;
                }
            }
        }
    }

    // Write contamination profile to disk
    {
        string fname = output_prefix + ".contam_prof";
        FILE* outf = fopen(fname.c_str(), "w");
        fprintf(stderr, "Writing contamination profile to disk...\n");
        dump_contam_prof(outf, contam_prof, contam_prof_conc, samples);
        fclose(outf);
    }
    // Write species-level profile when any species-level solver was active
    if ((species_regularize || use_interspecies) && !species_contam_prof_out.empty()){
        string fname = output_prefix + ".species_prof";
        FILE* outf = fopen(fname.c_str(), "w");
        fprintf(stderr, "Writing species-level profile to disk...\n");
        dump_species_prof(outf, species_contam_prof_out, species_contam_prof_conc_out);
        fclose(outf);
    }
    // Write contamination rate (and standard error) per cell to disk
    {
        string fname = output_prefix + ".contam_rate";
        FILE* outf = fopen(fname.c_str(), "w");
        dump_contam_rates(outf, contam_rate, contam_rate_se, samples,
            libname, cellranger, seurat, underscore);
        fclose(outf);
    }
    // Write allele ratio (genome A fraction) per cell to disk
    {
        string fname = output_prefix + ".allele_ratio";
        FILE* outf = fopen(fname.c_str(), "w");
        fprintf(stderr, "Writing allele ratios to disk (%lu heterotypic cells)...\n",
            allele_ratio.size());
        // Header
        fprintf(outf, "barcode\tallele_ratio_r\tallele_ratio_se\tcontam_rate_c\tcontam_rate_se\tidentity\n");
        for (robin_hood::unordered_map<unsigned long, double>::iterator ar = allele_ratio.begin();
            ar != allele_ratio.end(); ++ar){
            string bc_str = bc2str(ar->first);
            // Apply barcode formatting
            string bc_out;
            if (seurat && libname != ""){
                bc_out = libname + "_" + bc_str;
            } else if (cellranger){
                bc_out = bc_str + "-1";
            } else if (libname != ""){
                if (underscore){
                    bc_out = bc_str + "_" + libname;
                } else {
                    bc_out = bc_str + "-" + libname;
                }
            } else {
                bc_out = bc_str;
            }
            double r_val = ar->second;
            double r_se_val = allele_ratio_se.count(ar->first) > 0 ? allele_ratio_se[ar->first] : 0.0;
            double c_val = contam_rate.count(ar->first) > 0 ? contam_rate[ar->first] : -1.0;
            double c_se_val = contam_rate_se.count(ar->first) > 0 ? contam_rate_se[ar->first] : 0.0;
            string ident = "NA";
            if (assn.count(ar->first) > 0){
                ident = idx2name(assn[ar->first], samples);
            }
            fprintf(outf, "%s\t%f\t%f\t%f\t%f\t%s\n",
                bc_out.c_str(), r_val, r_se_val, c_val, c_se_val, ident.c_str());
        }
        fclose(outf);
    }

    // Write refined assignments to disk
    {
        string fname = output_prefix + ".decontam.assignments";
        FILE* outf = fopen(fname.c_str(), "w");
        dump_assignments(outf, assn, assn_llr, samples, libname,
            cellranger, seurat, underscore);
        fclose(outf);
    }
}

void parse_rates(string& filename, robin_hood::unordered_map<unsigned long, double>& contam_rates){
    ifstream inf(filename);
    string bc_str;
    double rate;
    double rate_se;
    while (inf >> bc_str >> rate >> rate_se){
        unsigned long bc_key = bc_ul(bc_str);
        contam_rates.emplace(bc_key, rate);
    }
}

void parse_prof(string& filename, map<int, double>& contam_prof, vector<string>& samples){

    // Note: sometimes contam_prof can contain an "other_species" entry, keyed to -1.
    // For this purpose, we are loading contam prof to help guide the GEX profiling.
    // That can't make use of individuals that aren't in the data set, so we will
    // exclude this individual if it exists and normalize the profile so it sums to 1.

    map<string, int> samp2idx;
    for (int i = 0; i < (int)samples.size(); ++i){
        samp2idx.insert(make_pair(samples[i], i));
    }

    ifstream inf(filename);
    string line;
    double proptot = 0.0;
    while(getline(inf, line)){
        istringstream splitter(line);
        string field;
        int fld_idx = 0;
        int samp_idx = -1;
        double samp_prop = 0.0;
        while(getline(splitter, field, '\t')){
            if (fld_idx == 0){
                if (samp2idx.count(field) > 0){
                    samp_idx = samp2idx[field];
                }
            }
            else if (fld_idx == 1){
                samp_prop = atof(field.c_str());
            }
            ++fld_idx;
        }

        if (samp_idx != -1){
            contam_prof.insert(make_pair(samp_idx, samp_prop));
            proptot += samp_prop;
        }
    }
    for (map<int, double>::iterator p = contam_prof.begin(); p != contam_prof.end();
        ++p){
        p->second /= proptot;
    }
}

void parse_skip_genes(string& skipgenesfile, set<string>& skipgenes){
    ifstream inf(skipgenesfile);
    string line;
    while (inf >> line){
        fprintf(stderr, "Skipping gene (raw expr output): %s\n", line.c_str());
        skipgenes.insert(line);
    }
}

void process_gex_data(string& output_prefix,
    string& barcodesfile,
    string& featuresfile,
    string& matrixfile,
    string& feature_type,
    string& clustfile,
    bool idfile_doublet_given,
    robin_hood::unordered_map<unsigned long, int>& assn,
    vector<string>& samples,
    robin_hood::unordered_map<unsigned long, double>& contam_rate,
    map<int, double>& contam_prof,
    bool round,
    int num_threads,
    string& libname,
    bool seurat,
    bool cellranger,
    bool underscore,
    string& skipgenesfile,
    string& skip_genes_regex){

    robin_hood::unordered_map<unsigned long, map<int, long int> > mtx;
    vector<string> features;
    fprintf(stderr, "Loading gene expression data...\n");
    bool success = parse_mex(barcodesfile, featuresfile, matrixfile, mtx, features, feature_type);
    if (!success){
        exit(1);
    }
    robin_hood::unordered_map<unsigned long, int> clusts;
    int nclusts = 0;
    vector<string> clustnames;
    if (clustfile != ""){
        nclusts = parse_clustfile(clustfile, clusts, clustnames);
    }
    else{
        fprintf(stderr, "Using cell identities as clusters\n");

        if (!idfile_doublet_given){
            clustnames = samples;
            for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
                a != assn.end(); ++a){
                if (a->second < (int)samples.size()){
                    clusts.emplace(a->first, a->second);
                }
            }
            nclusts = samples.size();
        }
        else{
            set<pair<string, int> > name_set;
            for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
                a != assn.end(); ++a){
                name_set.insert(make_pair(idx2name(a->second, samples), a->second));
            }
            nclusts = name_set.size();
            int idx = 0;
            map<int, int> idx2idx;
            for (set<pair<string, int> >::iterator n = name_set.begin(); n != name_set.end(); ++n){
                idx2idx.insert(make_pair(n->second, idx));
                clustnames.push_back(n->first);
                ++idx;
            }
            for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
                a != assn.end(); ++a){
                clusts.emplace(a->first, idx2idx[a->second]);
            }
        }
    }
    contam_profiler_gex contam_profiler(contam_rate, contam_prof,
        assn, samples.size(), idfile_doublet_given);
    contam_profiler.set_threads(num_threads);
    contam_profiler.set_mtx(mtx, features.size());
    contam_profiler.set_clusts(clusts, nclusts);
    if (round){
        contam_profiler.round_counts();
    }
    if (skipgenesfile != "" || skip_genes_regex != ""){
        set<string> skip_genes_txt;
        if (skipgenesfile != ""){
            parse_skip_genes(skipgenesfile, skip_genes_txt);
        }
        if (skip_genes_regex != ""){
            const std::regex skip_regex(skip_genes_regex);
            for (int i = 0; i < (int)features.size(); ++i){
                smatch matches;
                if (regex_search(features[i], matches, skip_regex)){
                    fprintf(stderr, "Skipping gene (raw expr output): %s\n", features[i].c_str());
                    skip_genes_txt.insert(features[i]);
                }
            }
        }

        if (skip_genes_txt.size() > 0){
            fprintf(stderr, "To avoid skipping genes, omit -g and set option -G \"\"\n");
        }
        map<string, int> gene2idx;
        for (int i = 0; i < (int)features.size(); ++i){
            gene2idx.insert(make_pair(features[i], i));
        }
        set<int> skipgenes;
        for (set<string>::iterator sgt = skip_genes_txt.begin(); sgt != skip_genes_txt.end();
            ++sgt){
            if (gene2idx.count(*sgt) > 0){
                skipgenes.insert(gene2idx[*sgt]);
            }
        }
        contam_profiler.skip_genes(skipgenes);
    }

    // Infer ambient RNA profile
    contam_profiler.get_profile();

    // Write to disk
    {
        string fname = output_prefix + ".gex_profile";
        FILE* outf = fopen(fname.c_str(), "w");
        fprintf(outf, "gene\tambient");
        for (int i = 0; i < nclusts; ++i){
            fprintf(outf, "\t%s", clustnames[i].c_str());
        }
        fprintf(outf, "\n");
        for (int i = 0; i < (int)features.size(); ++i){
            fprintf(outf, "%s\t%e", features[i].c_str(), contam_profiler.prof_ambient[i]);
            for (int j = 0; j < nclusts; ++j){
                fprintf(outf, "\t%e", contam_profiler.prof_clusts[j][i]);
            }
            fprintf(outf, "\n");
        }
        fclose(outf);
    }

    // Remove ambient RNA profile
    contam_profiler.decontam();

    // Write cleaned up matrix to disk
    string out_mtx = output_prefix + "_mtx";
    write_mex(out_mtx, contam_profiler.mtx_decontam,
        features, round, libname, cellranger, seurat, underscore);
}

int main(int argc, char *argv[]) {

    static struct option long_options[] = {
       {"output_prefix", required_argument, 0, 'o'},
       {"error_ref", required_argument, 0, 'e'},
       {"error_alt", required_argument, 0, 'E'},
       {"doublet_rate", required_argument, 0, 'D'},
       {"llr", required_argument, 0, 'l'},
       {"no_weights", no_argument, 0, 'w'},
       {"ids", required_argument, 0, 'i'},
       {"ids_doublet", required_argument, 0, 'I'},
       {"libname", required_argument, 0, 'n'},
       {"cellranger", no_argument, 0, 'C'},
       {"seurat", no_argument, 0, 'S'},
       {"underscore", no_argument, 0, 'U'},
       {"run_once", no_argument, 0, 'r'},
       {"bootstrap", required_argument, 0, 'b'},
       {"barcodes", required_argument, 0, 'B'},
       {"features", required_argument, 0, 'F'},
       {"matrix", required_argument, 0, 'M'},
       {"feature_type", required_argument, 0, 't'},
       {"clusts", required_argument, 0, 'c'},
       {"skip_genes", required_argument, 0, 'g'},
       {"skip_genes_regex", required_argument, 0, 'G'},
       {"num_threads", required_argument, 0, 'T'},
       {"noround", no_argument, 0, 'R'},
       {"expected_lines", required_argument, 0, 'X'},
       {"warm_start", required_argument, 0, 'W'},
       {"fixed_ambient", required_argument, 0, 'A'},
       {"panel_metadata", required_argument, 0, 'P'},
       {"help", no_argument, 0, 'h'},
       // Long-only options
       {"interindividual", no_argument, 0, 2001},
       {"interspecies", no_argument, 0, 2002},
       {"species_regularize", no_argument, 0, 2003},
       {"loo", no_argument, 0, 2004},
       {"species_counts", required_argument, 0, 2005},
       {"allow_legacy_species_regularize", no_argument, 0, 2012},
       {"r_feedback", no_argument, 0, 2006},
       {"adaptive_prior", no_argument, 0, 2007},
       {"condf", required_argument, 0, 2008},
       {"species_condf", required_argument, 0, 2009},
       {"thorough_multistart", no_argument, 0, 2010},
       {"no_adaptive_multistart", no_argument, 0, 2011},
       {0, 0, 0, 0}
    };

    // Set default values
    string output_prefix = "";
    double error_ref = 0.001;
    double error_alt = 0.001;
    double llr = 0.0;
    bool weight = true;
    string idfile;
    bool idfile_given = false;
    string idfile_doublet;
    bool idfile_doublet_given = false;
    string libname = "";
    bool cellranger = false;
    bool seurat = false;
    bool underscore = false;
    bool run_once = false;
    int bootstrap = 100;
    double doublet_rate = -1.0;
    int num_threads = 0;
    string skipgenesfile = "";

    string skip_genes_regex = R"(^MT-)";

    string barcodesfile = "";
    string featuresfile = "";
    string matrixfile = "";
    string feature_type = "";
    string clustfile = "";
    bool round = true;

    // Tetraploid-aware mode
    string expected_lines_file = "";

    // Panel selection
    bool use_interindividual = false;
    bool use_interspecies = false;

    // Warm start / fixed ambient
    string warm_start_file = "";
    string fixed_ambient_file = "";

    // Species regularization
    bool species_regularize = false;
    bool allow_legacy_species_regularize = false;
    string panel_metadata_file = "";
    string species_counts_file = "";

    // Solver options
    bool leave_one_out = false;
    bool r_feedback = false;
    bool adaptive_prior = false;
    bool thorough_multistart = false;
    bool adaptive_multistart = true;

    // Explicit condf/species_condf paths (override prefix-derived defaults)
    string user_condf_file = "";
    string user_species_condf_file = "";

    int option_index = 0;
    int ch;

    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "o:e:g:G:E:l:i:I:n:b:D:B:F:M:t:c:T:X:W:A:P:RrCSUwh", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                break;
            case 'h':
                help(0);
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 'X':
                expected_lines_file = optarg;
                break;
            case 'W':
                warm_start_file = optarg;
                break;
            case 'A':
                fixed_ambient_file = optarg;
                break;
            case 'P':
                panel_metadata_file = optarg;
                break;
            case 2001:
                use_interindividual = true;
                break;
            case 2002:
                use_interspecies = true;
                break;
            case 2003:
                species_regularize = true;
                break;
            case 2004:
                leave_one_out = true;
                break;
            case 2005:
                species_counts_file = optarg;
                break;
            case 2012:
                allow_legacy_species_regularize = true;
                break;
            case 2006:
                r_feedback = true;
                break;
            case 2007:
                adaptive_prior = true;
                break;
            case 2008:
                user_condf_file = optarg;
                break;
            case 2009:
                user_species_condf_file = optarg;
                break;
            case 2010:
                thorough_multistart = true;
                break;
            case 2011:
                adaptive_multistart = false;
                break;
            case 'g':
                skipgenesfile = optarg;
                break;
            case 'G':
                skip_genes_regex = optarg;
                break;
            case 'i':
                idfile = optarg;
                idfile_given = true;
                break;
            case 'I':
                idfile_doublet = optarg;
                idfile_doublet_given = true;
                break;
            case 'D':
                doublet_rate = atof(optarg);
                break;
            case 'e':
                error_ref = atof(optarg);
                break;
            case 'E':
                error_alt = atof(optarg);
                break;
            case 'l':
                llr = atof(optarg);
                break;
            case 'w':
                weight = false;
                break;
            case 'n':
                libname = optarg;
                break;
            case 'C':
                cellranger = true;
                break;
            case 'S':
                seurat = true;
                break;
            case 'U':
                underscore = true;
                break;
            case 'r':
                run_once = true;
                break;
            case 'b':
                bootstrap = atoi(optarg);
                break;
            case 'B':
                barcodesfile = optarg;
                break;
            case 'F':
                featuresfile = optarg;
                break;
            case 'M':
                matrixfile = optarg;
                break;
            case 't':
                feature_type = optarg;
                break;
            case 'c':
                clustfile = optarg;
                break;
            case 'R':
                round = false;
                break;
            case 'T':
                num_threads = atoi(optarg);
                break;
            default:
                help(0);
                break;
        }
    }

    // ========================================================================
    // Validation
    // ========================================================================
    if (output_prefix == ""){
        fprintf(stderr, "ERROR: output_prefix required\n");
        exit(1);
    }
    if (!use_interindividual && !use_interspecies){
        fprintf(stderr, "ERROR: exactly one of --interindividual or --interspecies is required\n");
        exit(1);
    }
    if (use_interindividual && use_interspecies){
        fprintf(stderr, "ERROR: --interindividual and --interspecies are mutually exclusive\n");
        exit(1);
    }
    if (species_regularize && use_interspecies){
        fprintf(stderr, "ERROR: --species_regularize is incompatible with --interspecies\n");
        exit(1);
    }
    if (species_regularize && panel_metadata_file.empty()){
        fprintf(stderr, "ERROR: --species_regularize requires --panel_metadata\n");
        exit(1);
    }
    if (!warm_start_file.empty() && !fixed_ambient_file.empty()){
        fprintf(stderr, "ERROR: --warm_start and --fixed_ambient are mutually exclusive\n");
        exit(1);
    }
    if (error_ref <= 0 || error_ref >= 1.0 || error_alt <= 0 || error_alt >= 1.0){
        fprintf(stderr, "ERROR: error rates must be between 0 and 1, exclusive.\n");
        exit(1);
    }
    if (idfile_given && idfile_doublet_given){
        fprintf(stderr, "ERROR: only one of -i and -I is allowed.\n");
        exit(1);
    }
    if (bootstrap <= 0){
        fprintf(stderr, "WARNING: bootstrapping disabled. Ambient RNA pool proportions will "
            "be reported without concentration parameters (variance will be unknown).\n");
    }
    if (doublet_rate != -1 && (doublet_rate < 0 || doublet_rate > 1)){
        fprintf(stderr, "ERROR: --doublet_rate/-D must be between 0 and 1, inclusive.\n");
        exit(1);
    }
    if ((barcodesfile != "" || featuresfile != "" || matrixfile != "") &&
        !(barcodesfile != "" && featuresfile != "" && matrixfile != "")){
        fprintf(stderr, "ERROR: if inferring gene expression profile, you must provide all\n");
        fprintf(stderr, "three of --barcodes/-B, --features/-F, and --matrix/-M\n");
        exit(1);
    }
    if (barcodesfile != "" && clustfile == ""){
        fprintf(stderr, "WARNING: inferring expression profile of contamination without "
            "cluster information. Assuming one default expression profile for each individual "
            "(results may be less accurate if there is much cell type heterogeneity).\n");
    }
    if (clustfile != "" && barcodesfile == ""){
        fprintf(stderr, "ERROR: --clusters/-c only applicable when loading gene expression data\n");
        exit(1);
    }
    if (!species_counts_file.empty() && !species_regularize){
        fprintf(stderr, "ERROR: --species_counts requires --species_regularize\n");
        exit(1);
    }

    // ========================================================================
    // Load native sample set
    // ========================================================================
    NativePanelMode panel_mode = use_interspecies
        ? NativePanelMode::INTERSPECIES_NATIVE
        : NativePanelMode::INTERINDIVIDUAL_NATIVE;

    string sample_name = output_prefix + (use_interspecies ? ".species_samples" : ".samples");
    vector<string> samples;
    if (file_exists(sample_name)){
        load_samples(sample_name, samples);
    }
    else{
        fprintf(stderr, "ERROR: no %s samples file found for %s. Expected %s\n",
            use_interspecies ? "species-native" : "individual-native",
            output_prefix.c_str(), sample_name.c_str());
        exit(1);
    }
    fprintf(stderr, "Panel mode: %s (%lu samples loaded from %s)\n",
        panel_mode_name(panel_mode), samples.size(), sample_name.c_str());

    // 1 thread means 0 threads (don't launch extra processes)
    if (num_threads <= 1){
        num_threads = 0;
    }

    // Map cell barcodes to numeric IDs of best individual assignments
    robin_hood::unordered_map<unsigned long, int> assn;

    // Map cell barcodes to log likelihood ratio of best individual assignments
    robin_hood::unordered_map<unsigned long, double> assn_llr;

    string assn_name = output_prefix + (use_interspecies ? ".species_assignments" : ".assignments");
    if (file_exists(assn_name)){
        fprintf(stderr, "Loading %s assignments from %s...\n",
            use_interspecies ? "species-native" : "individual-native", assn_name.c_str());
        load_assignments_from_file(assn_name, assn, assn_llr, samples);
        if (llr > 0.0){
            // Filter assignments.
            for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
                a != assn.end();){
                if (assn_llr[a->first] <= llr){
                    assn_llr.erase(a->first);
                    a = assn.erase(a);
                }
                else{
                    ++a;
                }
            }
            if (assn.size() == 0){
                fprintf(stderr, "ERROR: LLR filter too high; no assignments left to use.\n");
                exit(1);
            }
        }
    }
    else{
        fprintf(stderr, "ERROR: no assignments found for %s. Please run demux_vcf with same\n",
            output_prefix.c_str());
        fprintf(stderr, "output prefix.\n");
        exit(1);
    }

    // Optional weighted native species-composition model for --interspecies.
    // This uses the original individual assignment file as the biological identity
    // source, while still fitting contamination against native B/C/H/O species counts.
    map<unsigned long, map<int, double> > species_composition_overrides;
    if (use_interspecies){
        string indiv_assn_name = output_prefix + ".assignments";
        string indiv_samples_name = output_prefix + ".samples";
        fprintf(stderr, "Weighted species-composition activation check:\n");
        fprintf(stderr, "  use_interspecies=YES\n");
        fprintf(stderr, "  panel_metadata=%s\n",
            panel_metadata_file.empty() ? "<empty>" : panel_metadata_file.c_str());
        fprintf(stderr, "  individual assignments path=%s exists=%s\n",
            indiv_assn_name.c_str(), file_exists(indiv_assn_name) ? "YES" : "NO");
        fprintf(stderr, "  individual samples path=%s exists=%s\n",
            indiv_samples_name.c_str(), file_exists(indiv_samples_name) ? "YES" : "NO");

        if (!panel_metadata_file.empty()){
            species_composition_overrides = build_weighted_species_composition_overrides(
                indiv_assn_name, indiv_samples_name, panel_metadata_file, samples);
            if (species_composition_overrides.empty()){
                fprintf(stderr, "WARNING: no weighted species-composition overrides loaded; native species estimator will use .species_assignments only.\n");
                fprintf(stderr, "         Expected files: %s and %s; panel_metadata: %s\n",
                    indiv_assn_name.c_str(), indiv_samples_name.c_str(), panel_metadata_file.c_str());
            }
        }
        else{
            fprintf(stderr, "WARNING: weighted species-composition overrides disabled because --panel_metadata/-P was not supplied.\n");
        }
    }

    // Parse expected_lines for tetraploid-aware mode
    set<int> locked_identities;
    set<int> safe_singlets;
    set<int> expected_allowed_ids;
    set<int> expected_allowed_ids2;
    bool tetraploid_aware = false;
    bool ids_restricted = idfile_given || idfile_doublet_given;

    if (expected_lines_file != ""){
        tetraploid_aware = true;
        fprintf(stderr, "Tetraploid-aware mode: parsing %s\n", expected_lines_file.c_str());

        // Build name->index mapping from samples (which may now be species names
        // if --interspecies was set)
        map<string, int> name_to_idx;
        for (int i = 0; i < (int)samples.size(); i++){
            name_to_idx[samples[i]] = i;
        }

        set<string> combo_individuals;
        set<string> singlet_only;
        set<string> all_singlets;

        ifstream elf(expected_lines_file);
        if (!elf.is_open()){
            fprintf(stderr, "ERROR: cannot open expected_lines file: %s\n", expected_lines_file.c_str());
            exit(1);
        }
        string el_line;
        while (getline(elf, el_line)){
            if (el_line.empty() || el_line[0] == '#') continue;
            while (!el_line.empty() && (el_line.back() == '\r' || el_line.back() == '\n' || el_line.back() == ' '))
                el_line.pop_back();
            if (el_line.empty()) continue;

            size_t plus = el_line.find('+');
            if (plus != string::npos){
                string id1 = el_line.substr(0, plus);
                string id2 = el_line.substr(plus + 1);
                if (id1 == id2){
                    // Homotypic tetraploid (A+A): genotype-identical to singlet A.
                    // Do NOT create a combo identity; treat as a pure singlet.
                    all_singlets.insert(id1);
                    if (name_to_idx.count(id1) > 0){
                        int idx = name_to_idx[id1];
                        expected_allowed_ids.insert(idx);
                        expected_allowed_ids2.insert(idx);
                    }
                    fprintf(stderr, "  Homotypic entry %s+%s treated as singlet %s\n",
                        id1.c_str(), id2.c_str(), id1.c_str());
                    continue;
                }
                combo_individuals.insert(id1);
                combo_individuals.insert(id2);
                // Lock the combo identity index
                if (name_to_idx.count(id1) > 0 && name_to_idx.count(id2) > 0){
                    int idx1 = name_to_idx[id1];
                    int idx2 = name_to_idx[id2];
                    int combo_idx;
                    if (idx1 < idx2) combo_idx = hap_comb_to_idx(idx1, idx2, samples.size());
                    else combo_idx = hap_comb_to_idx(idx2, idx1, samples.size());
                    locked_identities.insert(combo_idx);
                    expected_allowed_ids.insert(combo_idx);
                    expected_allowed_ids2.insert(combo_idx);
                    expected_allowed_ids.insert(idx1);
                    expected_allowed_ids.insert(idx2);
                }
            } else {
                all_singlets.insert(el_line);
                if (name_to_idx.count(el_line) > 0){
                    int idx = name_to_idx[el_line];
                    expected_allowed_ids.insert(idx);
                    expected_allowed_ids2.insert(idx);
                }
            }
        }

        // Build locked and safe singlet sets
        for (const auto& s : all_singlets){
            if (name_to_idx.count(s) > 0){
                int idx = name_to_idx[s];
                if (combo_individuals.count(s) > 0){
                    locked_identities.insert(idx);
                } else {
                    safe_singlets.insert(idx);
                }
            }
        }

        // Also lock singlet indices for individuals in combos but NOT in singlet list
        for (const auto& ci : combo_individuals){
            if (name_to_idx.count(ci) > 0 && all_singlets.count(ci) == 0){
                locked_identities.insert(name_to_idx[ci]);
            }
        }

        if (!expected_allowed_ids.empty()){
            ids_restricted = true;
        }

        fprintf(stderr, "  Expected-lines legal identities: %lu exact, %lu with component singlets; "
            "%d ambient source individuals\n",
            expected_allowed_ids2.size(), expected_allowed_ids.size(),
            count_singlet_ids(expected_allowed_ids, samples.size()));
        fprintf(stderr, "  Locked identities: %lu (combos + ambiguous singlets)\n", locked_identities.size());
        fprintf(stderr, "  Safe singlets: %lu (pure diploids, reassignment allowed)\n", safe_singlets.size());
    }

    // ========================================================================
    // Check for existing outputs / run estimation
    // ========================================================================
    string prof_name = output_prefix + ".contam_prof";
    string rate_name = output_prefix + ".contam_rate";

    robin_hood::unordered_map<unsigned long, double> contam_rate;
    map<int, double> contam_prof;

    bool load_gex = (barcodesfile != "" && featuresfile != "" && matrixfile != "");

    if (file_exists(prof_name) && file_exists(rate_name)){
        if (!load_gex){
            fprintf(stderr, "ERROR: previous run detected, and no gene expression data "
                "provided. Nothing to do.\n");
            fprintf(stderr, "To repeat the run, remove %s and %s.\n", prof_name.c_str(),
                rate_name.c_str());
            exit(1);
        }
        fprintf(stderr, "Previous run detected. Loading profile\n");
        fprintf(stderr, "To prevent this behavior, change --output_prefix/-o or remove "
            "%s and %s.\n", prof_name.c_str(), rate_name.c_str());
        parse_rates(rate_name, contam_rate);
        parse_prof(prof_name, contam_prof, samples);
    }
    else{
        infer_from_genotypes(output_prefix,
            assn,
            assn_llr,
            idfile,
            idfile_given,
            idfile_doublet,
            idfile_doublet_given,
            samples,
            contam_prof,
            contam_rate,
            doublet_rate,
            num_threads,
            error_ref,
            error_alt,
            weight,
            run_once,
            bootstrap,
            libname,
            seurat,
            cellranger,
            underscore,
            locked_identities,
            safe_singlets,
            tetraploid_aware,
            ids_restricted,
            expected_allowed_ids,
            expected_allowed_ids2,
            leave_one_out,
            r_feedback,
            adaptive_prior,
            thorough_multistart,
            adaptive_multistart,
            warm_start_file,
            fixed_ambient_file,
            species_regularize,
            panel_metadata_file,
            species_counts_file,
            use_interspecies,
            user_condf_file,
            user_species_condf_file,
            species_composition_overrides);
    }

    if (load_gex){
        process_gex_data(output_prefix,
            barcodesfile,
            featuresfile,
            matrixfile,
            feature_type,
            clustfile,
            idfile_doublet_given,
            assn,
            samples,
            contam_rate,
            contam_prof,
            round,
            num_threads,
            libname,
            seurat,
            cellranger,
            underscore,
            skipgenesfile,
            skip_genes_regex);
    }
    return 0;
}
