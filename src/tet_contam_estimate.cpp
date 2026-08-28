#include <limits>
#include <cfloat>
#include <cerrno>
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
#include <iomanip>
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
#include <stdexcept>
#include <zlib.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include <htswrapper/mex.h>
#include "common.h"
#include "ambient_rna_three_ap.h"
#include "ambient_rna_gex.h"
#include "io.h"

using std::cout;
using std::endl;
using namespace std;

// Version tracking
const string QC_VERSION = "4.7-fixed-assignment-comparison";
const string QC_VERSION_MSG = "Adds iterative fixed-assignment fitting for controlled identity/ambient-roster comparisons while preserving historical behavior unless --freeze_assignments is requested";
const string TOOL_NAME = "tet_contam_estimate";

static string assignment_update_mode_name(bool run_once, bool freeze_assignments){
    if (run_once) return "single_pass_frozen";
    if (freeze_assignments) return "iterative_frozen";
    return "iterative_reclassification";
}

static void validate_assignment_update_mode(bool run_once, bool freeze_assignments){
    if (run_once && freeze_assignments){
        fprintf(stderr,
            "ERROR: --run_once and --freeze_assignments are mutually exclusive; "
            "use --freeze_assignments for iterative fixed-assignment fitting.\n");
        exit(1);
    }
}

static string json_escape(const string& value){
    string out;
    out.reserve(value.size() + 8);
    for (char ch : value){
        switch(ch){
            case '\\': out += "\\\\"; break;
            case '"': out += "\\\""; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default: out.push_back(ch); break;
        }
    }
    return out;
}

static long long file_size_or_minus_one(const string& path){
    if (path.empty()) return -1;
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return -1;
    return (long long)st.st_size;
}

static void write_run_contract_json(
    const string& path,
    const string& output_prefix,
    const string& run_class,
    bool production_contract_pass,
    const string& production_contract_reason,
    bool use_interspecies,
    const string& counts_path,
    const string& condf_path,
    const string& assignments_path,
    const string& samples_path,
    const string& expected_lines_path,
    const string& ambient_candidates_path,
    const string& warm_start_path,
    const string& fixed_ambient_path,
    const string& fix_r_path,
    const string& assignments_basis,
    const string& expected_lines_basis,
    const string& ambient_candidates_basis,
    const string& warm_start_basis,
    const string& fixed_ambient_basis,
    const string& fix_r_basis,
    bool strict_condf,
    bool run_once,
    bool freeze_assignments,
    const string& condition_key,
    const string& synthetic_id,
    double source_exclusion_strength,
    bool source_exclusion_explicit,
    const string& profile_holdout_path,
    const string& profile_holdout_basis,
    unsigned long profile_holdout_count,
    const string& surface_selector_path,
    const string& surface_output_path){

    ofstream out(path.c_str());
    if (!out.is_open()){
        fprintf(stderr, "ERROR: could not write run contract %s\n", path.c_str());
        exit(1);
    }
    auto field = [&](const string& key, const string& value, bool comma=true){
        out << "  \"" << json_escape(key) << "\": \"" << json_escape(value) << "\"";
        if (comma) out << ",";
        out << "\n";
    };
    auto path_field = [&](const string& key, const string& value, bool comma=true){
        out << "  \"" << json_escape(key) << "\": {\"path\": \""
            << json_escape(value) << "\", \"size_bytes\": " << file_size_or_minus_one(value) << "}";
        if (comma) out << ",";
        out << "\n";
    };
    out << "{\n";
    field("contract_version", "tet_contam_estimate_run_contract_V1_R4");
    field("tool", TOOL_NAME);
    field("tool_version", QC_VERSION);
    field("run_class", run_class);
    out << "  \"production_contract_pass\": " << (production_contract_pass ? "true" : "false") << ",\n";
    field("production_contract_reason", production_contract_reason);
    field("panel_mode", use_interspecies ? "interspecies" : "interindividual");
    field("output_prefix", output_prefix);
    out << "  \"strict_condf\": " << (strict_condf ? "true" : "false") << ",\n";
    out << "  \"run_once\": " << (run_once ? "true" : "false") << ",\n";
    out << "  \"freeze_assignments\": "
        << (freeze_assignments ? "true" : "false") << ",\n";
    field("assignment_update_mode",
        assignment_update_mode_name(run_once, freeze_assignments));
    field("assignments_basis", assignments_basis);
    field("expected_lines_basis", expected_lines_basis);
    field("ambient_candidates_basis", ambient_candidates_basis);
    field("warm_start_basis", warm_start_basis);
    field("fixed_ambient_basis", fixed_ambient_basis);
    field("fixed_r_basis", fix_r_basis);
    field("fix_r_basis", fix_r_basis);
    field("condition_key", condition_key);
    field("synthetic_id", synthetic_id);
    out << "  \"source_exclusion_strength\": " << std::setprecision(17)
        << source_exclusion_strength << ",\n";
    out << "  \"source_exclusion_explicit\": "
        << (source_exclusion_explicit ? "true" : "false") << ",\n";
    field("profile_holdout_basis", profile_holdout_basis);
    out << "  \"profile_holdout_count\": " << profile_holdout_count << ",\n";
    out << "  \"all_source_columns_retained\": "
        << (source_exclusion_strength <= 1e-12 ? "true" : "false") << ",\n";
    out << "  \"fixed_r_enabled\": " << (!fix_r_path.empty() ? "true" : "false") << ",\n";
    out << "  \"fixed_ambient_enabled\": "
        << (!fixed_ambient_path.empty() ? "true" : "false") << ",\n";
    out << "  \"truth_assisted\": "
        << ((run_class == "oracle") ? "true" : "false") << ",\n";
    path_field("counts", counts_path);
    path_field("condf", condf_path);
    path_field("assignments", assignments_path);
    path_field("samples", samples_path);
    path_field("expected_lines", expected_lines_path);
    path_field("ambient_candidates", ambient_candidates_path);
    path_field("warm_start", warm_start_path);
    path_field("fixed_ambient", fixed_ambient_path);
    path_field("fix_r", fix_r_path);
    path_field("profile_holdout_barcodes", profile_holdout_path);
    path_field("r_c_surface_selector", surface_selector_path);
    path_field("r_c_surface_out", surface_output_path, false);
    out << "}\n";
    out.close();
}

static string trim_copy(const string& x){
    size_t a = x.find_first_not_of(" \t\r\n");
    if (a == string::npos) return "";
    size_t b = x.find_last_not_of(" \t\r\n");
    return x.substr(a, b - a + 1);
}

static set<unsigned long> load_profile_holdout_barcodes(
    const string& path){
    set<unsigned long> out;
    if (path.empty()) return out;
    ifstream in(path.c_str());
    if (!in.is_open()){
        fprintf(stderr, "ERROR: could not open --profile_holdout_barcodes file %s\n",
            path.c_str());
        exit(1);
    }
    string line;
    while (getline(in, line)){
        line = trim_copy(line);
        if (line.empty() || line[0] == '#') continue;
        istringstream iss(line);
        string barcode;
        iss >> barcode;
        if (barcode.empty() || barcode == "barcode") continue;
        out.insert(bc_ul(barcode));
    }
    if (out.empty()){
        fprintf(stderr, "ERROR: --profile_holdout_barcodes %s contained no barcodes.\n",
            path.c_str());
        exit(1);
    }
    return out;
}

static set<string> load_r_c_surface_selector(
    const string& path, const string& synthetic_id){
    set<string> selected;
    if (path.empty()) return selected;
    ifstream in(path.c_str());
    if (!in.is_open()){
        fprintf(stderr, "ERROR: could not open --r_c_surface_selector %s\n", path.c_str());
        exit(1);
    }
    string line;
    if (!getline(in, line)){
        fprintf(stderr, "ERROR: empty --r_c_surface_selector %s\n", path.c_str());
        exit(1);
    }
    vector<string> header;
    string field;
    stringstream hs(line);
    while (getline(hs, field, '\t')) header.push_back(trim_copy(field));
    int sid_col = -1, barcode_col = -1;
    for (size_t i = 0; i < header.size(); ++i){
        if (header[i] == "synthetic_id") sid_col = (int)i;
        if (header[i] == "selected_barcode" || header[i] == "barcode") barcode_col = (int)i;
    }
    if (sid_col < 0 || barcode_col < 0){
        fprintf(stderr, "ERROR: %s must contain synthetic_id and selected_barcode columns\n",
            path.c_str());
        exit(1);
    }
    int line_no = 1;
    while (getline(in, line)){
        ++line_no;
        if (trim_copy(line).empty()) continue;
        vector<string> fields;
        stringstream ss(line);
        while (getline(ss, field, '\t')) fields.push_back(trim_copy(field));
        if ((size_t)std::max(sid_col, barcode_col) >= fields.size()){
            fprintf(stderr, "ERROR: malformed selector row %d in %s\n", line_no, path.c_str());
            exit(1);
        }
        if (fields[sid_col] == synthetic_id && !fields[barcode_col].empty()){
            selected.insert(fields[barcode_col]);
        }
    }
    fprintf(stderr, "Loaded %lu r-C surface barcode(s) for synthetic_id=%s from %s\n",
        selected.size(), synthetic_id.c_str(), path.c_str());
    return selected;
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
    fprintf(stderr, "    --species_condf FILE  Override .species_condf file path (interspecies mode)\n");
    fprintf(stderr, "    --strict_condf       Fail before fitting if any observed category x active donor conditional is missing.\n\n");

    fprintf(stderr, "===== AMBIENT PROFILE ROW DESIGN (individual path only) =====\n");
    fprintf(stderr, "    --candidate_keyed_rows      Also fit the ambient profile on the marginal allele-count bins of every\n");
    fprintf(stderr, "                                ambient candidate, not only the bins keyed on each cell's own assigned\n");
    fprintf(stderr, "                                identity. The assignment-keyed design alone is rank deficient by 13 to 16\n");
    fprintf(stderr, "                                of 23 free profile directions, so those donors cannot be separated at any\n");
    fprintf(stderr, "                                depth. Adds no new data: the bins are already loaded.\n");
    fprintf(stderr, "    --candidate_keyed_split F   Fraction of the per-read weight budget kept by the assignment-keyed block,\n");
    fprintf(stderr, "                                default 0.05, measured. The candidate-keyed block receives 1-F over the candidate\n");
    fprintf(stderr, "                                count, so total evidence per read is unchanged. 1.0 reproduces the historical\n");
    fprintf(stderr, "                                design; 0.0 uses candidate-keyed marginals alone.\n\n");

    fprintf(stderr, "===== AMBIENT PROFILE INITIALIZATION (at most one) =====\n");
    fprintf(stderr, "    --warm_start -W FILE   Load ambient profile as solver starting point.\n");
    fprintf(stderr, "                           FILE is .contam_prof_empty or .species_prof_empty.\n");
    fprintf(stderr, "                           Profile is freely refined during estimation.\n");
    fprintf(stderr, "    --fixed_ambient -A FILE  Lock ambient profile (oracle/diagnostic only).\n");
    fprintf(stderr, "    --run_class CLASS  unspecified, production, or oracle. Production rejects answer-key inputs.\n");
    fprintf(stderr, "    --run_contract FILE  Write machine-readable input/run contract (default: <prefix>.run_contract.json).\n\n");

    fprintf(stderr, "===== TETRAPLOID-AWARE =====\n");
    fprintf(stderr, "    --expected_lines -X FILE  Expected receiver lines file (locks/restricts receiver identities)\n");
    fprintf(stderr, "    --ambient_candidates FILE  Legal ambient source singlets, separate from receiver identities\n");
    fprintf(stderr, "    --assignments_basis BASIS  library, truth, or unspecified (required=library for production).\n");
    fprintf(stderr, "    --expected_lines_basis BASIS  library, truth, or unspecified.\n");
    fprintf(stderr, "    --ambient_candidates_basis BASIS  library, truth, or unspecified.\n");
    fprintf(stderr, "    --warm_start_basis BASIS  library, truth, or unspecified.\n");
    fprintf(stderr, "    --fixed_ambient_basis BASIS  library, truth, or unspecified.\n");
    fprintf(stderr, "    --fix_r_basis BASIS  library, truth, or unspecified.\n\n");

    fprintf(stderr, "===== SOLVER OPTIONS =====\n");
    fprintf(stderr, "    --leave_source_out  Diagnostic: exact legacy hard source exclusion (lambda=1)\n");
    fprintf(stderr, "    --loo               Deprecated alias for --leave_source_out\n");
    fprintf(stderr, "    --source_exclusion_strength FLOAT  Continuous source exclusion in [0,1]; 0 is exact global profile, 1 is exact hard exclusion\n");
    fprintf(stderr, "    --profile_holdout_barcodes FILE  Exclude listed receiver cells from ambient-profile training only; all cells remain scored\n");
    fprintf(stderr, "    --profile_holdout_basis BASIS  library, truth, or unspecified; production cross-fitting requires library\n");
    fprintf(stderr, "    --r_c_surface_selector FILE  Signed selector TSV with synthetic_id and selected_barcode\n");
    fprintf(stderr, "    --r_c_surface_out FILE  Write selector-limited r-C likelihood surface TSV\n");
    fprintf(stderr, "    --condition_key KEY  Condition provenance written to surfaces/contracts\n");
    fprintf(stderr, "    --synthetic_id ID  Unit provenance and selector filter\n");
    fprintf(stderr, "    --r_feedback        Feed per-cell allele ratios into ambient profile estimation\n");
    fprintf(stderr, "    --fix_r FILE        Do NOT fit r. Read a fixed allele ratio per identity from\n");
    fprintf(stderr, "                        FILE (2 cols: identity<TAB>r, e.g. C40670+CongoA4B<TAB>0.4941)\n");
    fprintf(stderr, "                        and solve c alone. Identities absent from FILE keep the\n");
    fprintf(stderr, "                        normal free-r fit. The pooling policy is FILE's business.\n");
    fprintf(stderr, "    --adaptive_prior    Auto-detect pathological distributions, apply fixed prior fallback\n");
    fprintf(stderr, "    --run_once -r       Single pass, no iterative convergence; supplied identities remain fixed\n");
    fprintf(stderr, "                        Mutually exclusive with --freeze_assignments.\n");
    fprintf(stderr, "    --freeze_assignments  Fit profile and per-cell c/r to iterative convergence without changing supplied identities\n");
    fprintf(stderr, "                        Mutually exclusive with --run_once.\n");
    fprintf(stderr, "    --num_threads -T    Parallel threads (default: 1)\n");
    fprintf(stderr, "    --bootstrap -b      Bootstrap replicates (default: 100)\n");
    fprintf(stderr, "    --heterotypic_start_mode MODE  single, topk (default), or exhaustive\n");
    fprintf(stderr, "    --heterotypic_start_top_k N   BFGS-refine top N deterministic basins (default: 5)\n");
    fprintf(stderr, "    --thorough_multistart  Deprecated alias for exhaustive start mode\n");
    fprintf(stderr, "    --no_profile_intervals  Disable adaptive ridge-aware profile intervals\n");
    fprintf(stderr, "    --contam_prior_mode MODE  none, global, heterotypic, or fusion (default)\n");
    fprintf(stderr, "    --contam_prior_min_cells N  Minimum supported cells per prior group (default: 20)\n");
    fprintf(stderr, "    --contam_prior_max_ci_width X  Maximum MLE c interval width for prior training (default: 0.50)\n");
    fprintf(stderr, "    --contam_prior_min_weight X  Minimum informative allele weight (default: 10)\n");
    fprintf(stderr, "    --contam_prior_max_gradient X  DEPRECATED compatibility option; ignored\n");
    fprintf(stderr, "    --profile_restarts N  Exact total ambient-profile starts (>=1; supplied/warm,\n");
    fprintf(stderr, "                           uniform, then deterministic Dirichlet starts).\n");
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

static void restrict_profile_to_explicit_candidates(
    map<int, double>& profile,
    map<int, double>& concentration,
    const set<int>& candidates,
    const string& label){
    if (candidates.empty()) return;
    int removed = 0;
    for (auto it = profile.begin(); it != profile.end(); ){
        if (candidates.count(it->first) == 0){
            concentration.erase(it->first);
            it = profile.erase(it);
            removed++;
        }
        else ++it;
    }
    double total = 0.0;
    for (const auto& kv : profile){
        if (std::isfinite(kv.second) && kv.second > 0.0) total += kv.second;
    }
    if (!(total > 0.0)){
        fprintf(stderr, "ERROR: %s has no positive mass on the explicit ambient candidate set.\n", label.c_str());
        exit(1);
    }
    for (auto& kv : profile){
        kv.second = std::max(0.0, kv.second) / total;
    }
    if (removed > 0){
        fprintf(stderr, "WARNING: removed %d %s entries outside the explicit ambient candidate set.\n",
            removed, label.c_str());
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
    bool freeze_assignments,
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
    const set<int>& ambient_candidate_ids,
    // Solver options
    bool leave_one_out,
    double source_exclusion_strength,
    bool r_feedback,
    const std::string& fix_r_file,
    const std::string& fixed_r_basis,
    const std::string& fixed_ambient_basis,
    bool truth_assisted_condition,
    const std::set<std::string>& surface_selector_barcodes,
    const std::string& surface_output_file,
    const std::string& condition_key,
    const std::string& synthetic_id,
    bool adaptive_prior,
    bool thorough_multistart,
    bool adaptive_multistart,
    int profile_restarts,
    const string& heterotypic_start_mode,
    int heterotypic_start_top_k,
    bool adaptive_profile_intervals,
    const string& contam_prior_mode,
    int contam_prior_min_cells,
    double contam_prior_max_ci_width,
    double contam_prior_min_weight,
    double contam_prior_max_gradient,
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
    bool strict_condf,
    // Optional native species composition overrides derived from original individual calls
    const map<unsigned long, map<int, double> >& species_composition_overrides,
    // Receiver-data cross-fitting: exclude these cells only from profile training.
    const set<unsigned long>& profile_holdout_barcodes,
    const string& profile_holdout_basis,
    // O53 step 2: candidate-keyed ambient-profile rows
    bool candidate_keyed_rows,
    double candidate_keyed_split){

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

    // Ambient source candidates are a distinct model dimension.  When the new
    // file is supplied, it controls only the profile simplex; expected-lines
    // continues to control receiver locking/reassignment.  Legacy runs without
    // --ambient_candidates retain the old expected-lines behavior.
    set<int> profile_allowed_ids = ambient_candidate_ids.empty() ? allowed_ids : ambient_candidate_ids;
    if (!ambient_candidate_ids.empty()){
        fprintf(stderr, "Separate ambient candidate restriction active: %lu source singlets; "
            "receiver reassignment remains restricted to %lu identities.\n",
            profile_allowed_ids.size(), allowed_ids.size());
    }
    if (profile_allowed_ids.empty()){
        fprintf(stderr, "ERROR: ambient source candidate set is empty.\n");
        exit(1);
    }

    // Counts must retain both receiver identity rows and source singlet rows.
    set<int> count_allowed_ids = allowed_ids;
    count_allowed_ids.insert(profile_allowed_ids.begin(), profile_allowed_ids.end());

    enforce_native_dimensions(use_interspecies ? NativePanelMode::INTERSPECIES_NATIVE : NativePanelMode::INTERINDIVIDUAL_NATIVE,
        counts_name, condf_name, samples);

    // Load stored allele counts
    robin_hood::unordered_map<unsigned long, map<pair<int, int>,
        map<pair<int, int>, pair<float, float> > > > indv_allelecounts;
    if (file_exists(counts_name)){
        fprintf(stderr, "Loading counts from %s...\n", counts_name.c_str());
        load_counts_from_file(indv_allelecounts, samples, counts_name, count_allowed_ids);
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

    // Carried out of the iteration loop with contam_prof, because the
    // contamFinder is scoped to one iteration and goes out of scope before the
    // outputs are written.
    double final_loglik_out = 0.0;
    bool final_loglik_valid_out = false;
    robin_hood::unordered_map<unsigned long, double> contam_rate_se;
    robin_hood::unordered_map<unsigned long, double> allele_ratio;
    robin_hood::unordered_map<unsigned long, double> allele_ratio_se;
    map<string, double> species_contam_prof_out;
    map<string, double> species_contam_prof_conc_out;
    int profile_successful_starts_out = 0;
    int profile_near_optimal_count_out = 0;
    double profile_best_ll_out = -DBL_MAX;
    double profile_second_best_ll_out = -DBL_MAX;
    double profile_near_optimal_l1_spread_out = std::numeric_limits<double>::quiet_NaN();
    // Multistart diagnostics are carried separately from the refinement-solve
    // diagnostics above. See ambient_rna_three_ap.h for why the two differ.
    bool multistart_attempted_out = false;
    int multistart_configured_starts_out = 0;
    int multistart_successful_starts_out = 0;
    int multistart_near_optimal_count_out = 0;
    double multistart_best_ll_out = -DBL_MAX;
    double multistart_second_best_ll_out = -DBL_MAX;
    double multistart_near_optimal_l1_spread_out = std::numeric_limits<double>::quiet_NaN();

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
        restrict_profile_to_explicit_candidates(
            fixed_indiv_prof, fixed_indiv_conc, ambient_candidate_ids, "fixed ambient profile");
        fixed_ambient_loaded = true;
    }

    if (freeze_assignments && !run_once){
        fprintf(stderr, "Assignment update mode: iterative frozen assignments; "
            "profile and per-cell c/r will converge without reclassification.\n");
    }

    int nits = 0;
    while (delta > delta_thresh){
        fprintf(stderr, "===== ITERATION %d =====\n", nits+1);
        contamFinder3 cf(indv_allelecounts, assn, assn_llr, exp_match_fracs, samples.size(),
            profile_allowed_ids, allowed_ids, allowed_ids2);
        if (nits == 0){
            contamFinder3::CondfCoverageReport condf_report = cf.write_condf_coverage_report(
                output_prefix + ".condf_coverage.tsv", samples);
            if (strict_condf && condf_report.missing_lookups > 0){
                fprintf(stderr,
                    "ERROR: --strict_condf rejected %llu missing conditional-fraction lookups "
                    "covering %.6g weighted observations. See %s.condf_coverage.tsv\n",
                    condf_report.missing_lookups, condf_report.missing_weight, output_prefix.c_str());
                exit(1);
            }
        }
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
        cf.set_heterotypic_start_mode(thorough_multistart ? "exhaustive" : heterotypic_start_mode,
            heterotypic_start_top_k);
        cf.set_adaptive_multistart(adaptive_multistart);
        cf.set_adaptive_profile_intervals(adaptive_profile_intervals);
        cf.set_contam_prior_mode(contam_prior_mode);
        cf.set_contam_prior_support(contam_prior_min_cells, contam_prior_max_ci_width,
            contam_prior_min_weight, contam_prior_max_gradient);
        if (profile_restarts > 0){
            cf.set_profile_total_starts(profile_restarts);
        }
        if (run_once || freeze_assignments){
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

        // O53 step 2: candidate-keyed ambient-profile rows.
        //
        // Individual path only. Candidate-keyed expansion of the species path
        // was measured and rejected, worse on 20 of 20 native species units, so
        // requesting it there is refused rather than silently honoured.
        if (candidate_keyed_rows && use_interspecies){
            fprintf(stderr, "ERROR: --candidate_keyed_rows applies to the individual path "
                "only.\n");
            fprintf(stderr, "       Candidate-keyed expansion of the species path was measured "
                "and rejected.\n");
            exit(1);
        }
        if (candidate_keyed_rows){
            cf.set_candidate_keyed_rows(true, candidate_keyed_split);
            if (nits == 0){
                fprintf(stderr, "Candidate-keyed profile rows enabled: assignment-keyed weight "
                    "share %.3f, candidate-keyed share %.3f split over the ambient candidate "
                    "roster\n", candidate_keyed_split, 1.0 - candidate_keyed_split);
            }
        }

        // Source exclusion. The explicit strength setter owns both endpoints;
        // --loo/--leave_source_out are resolved to lambda=1 in main().
        cf.set_source_exclusion_strength(source_exclusion_strength);
        if (!profile_holdout_barcodes.empty()){
            cf.set_profile_holdout_barcodes(profile_holdout_barcodes, profile_holdout_basis);
            if (nits == 0){
                fprintf(stderr, "Cross-fitted profile training excludes %lu held-out cells; "
                    "all cells remain in C/r scoring.\n",
                    (unsigned long)profile_holdout_barcodes.size());
            }
        }
        if (source_exclusion_strength > 0.0 && nits == 0){
            fprintf(stderr, "Source exclusion enabled at lambda=%.6g%s\n",
                source_exclusion_strength, leave_one_out ? " (legacy hard-exclusion request)" : "");
        }

        // Ambient profile refinements: r-feedback and adaptive prior
        // ---- Fixed-r mode (pooled-r experiment) --------------------------
        // r is supplied, not fitted. See --fix_r in the usage text. The file
        // maps an identity string, exactly as written into the `identity`
        // column of .allele_ratio, to a fixed allele ratio.
        if (fix_r_file != ""){
            map<int, double> fixed_r;
            set<int> idents_seen;
            for (robin_hood::unordered_map<unsigned long, int>::iterator ai = assn.begin();
                ai != assn.end(); ++ai){
                idents_seen.insert(ai->second);
            }
            map<string, int> name2idx;
            int n_heterotypic_identities = 0;
            for (set<int>::iterator ii = idents_seen.begin(); ii != idents_seen.end(); ++ii){
                name2idx.insert(make_pair(idx2name(*ii, samples), *ii));
                if (*ii >= (int)samples.size()){
                    pair<int,int> combo = idx_to_hap_comb(*ii, (int)samples.size());
                    if (combo.first != combo.second) n_heterotypic_identities++;
                }
            }
            ifstream frf(fix_r_file.c_str());
            if (!frf.is_open()){
                fprintf(stderr, "ERROR: could not open --fix_r file %s\n", fix_r_file.c_str());
                exit(1);
            }
            string line;
            int nline = 0;
            int nmatched = 0;
            while (getline(frf, line)){
                nline++;
                if (line.size() == 0 || line[0] == '#') continue;
                size_t tab = line.find('\t');
                if (tab == string::npos) continue;
                string ident = line.substr(0, tab);
                string rstr = line.substr(tab + 1);
                if (ident == "identity") continue;   // tolerate a header row
                double rv = atof(rstr.c_str());
                if (rv < 0.0 || rv > 1.0){
                    fprintf(stderr, "ERROR: --fix_r line %d: r=%f is outside [0,1]\n", nline, rv);
                    exit(1);
                }
                map<string, int>::iterator ni = name2idx.find(ident);
                if (ni != name2idx.end()){
                    fixed_r.insert(make_pair(ni->second, rv));
                    nmatched++;
                }
            }
            frf.close();
            fprintf(stderr, "--fix_r: %d identity/identities matched in this library.\n", nmatched);
            if (nmatched == 0 && n_heterotypic_identities > 0){
                fprintf(stderr, "ERROR: --fix_r matched no heterotypic identity present in this library. "
                    "Nothing would be fixed, which is silently not the requested experiment.\n");
                exit(1);
            }
            if (nmatched == 0){
                fprintf(stderr, "--fix_r: this unit contains no heterotypic identities; "
                    "the signed header-only artifact leaves diploid fits unchanged.\n");
            }
            cf.set_fixed_r(fixed_r);
        }

        if (r_feedback){
            cf.set_r_feedback(true);
        }
        if (adaptive_prior){
            cf.set_adaptive_prior(true, adaptive_mean, adaptive_init_var, 0.001, adaptive_thresh);
        }

        // Initialize to the final profile from the previous outer iteration.
        //
        // With --ambient_candidates, the ambient-profile simplex is explicitly
        // fixed for the entire run and is intentionally independent of receiver
        // assignments.  Reclassification can change assn/allowed receiver IDs,
        // but it must not be interpreted as a change in the source-candidate
        // universe.  set_init_contam_prof() already filters/fills against the
        // constructor's profile_allowed_ids, so carry the profile forward
        // unconditionally in this mode.
        //
        // Legacy runs without --ambient_candidates retain the historical
        // assignment-derived compatibility check because their source universe
        // is still coupled to active/allowed receiver identities.
        if (nits > 0 && !contam_prof.empty()){
            if (!ambient_candidate_ids.empty()){
                cf.set_init_contam_prof(contam_prof);
                fprintf(stderr, "  Carrying forward ambient profile across fixed explicit "
                    "source candidate set (%lu singlets)\n",
                    profile_allowed_ids.size());
            }
            else{
                // Build the set of singlet IDs that the new contamFinder will use.
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
                // Also include any from allowed_ids/allowed_ids2.
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

                // Check if old contam_prof keys (excluding -1) match new singlets.
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
                restrict_profile_to_explicit_candidates(
                    loaded_prof, loaded_conc, ambient_candidate_ids, "warm-start ambient profile");
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

        cf.set_diagnostic_context(fixed_r_basis, fixed_ambient_basis,
            truth_assisted_condition);

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
            final_loglik_out = cf.final_loglik;
            final_loglik_valid_out = cf.final_loglik_valid;
            contam_rate = cf.contam_rate;
            contam_rate_se = cf.contam_rate_se;
            allele_ratio = cf.allele_ratio;
            allele_ratio_se = cf.allele_ratio_se;
            species_contam_prof_out = cf.species_contam_prof;
            species_contam_prof_conc_out = cf.species_contam_prof_conc;
            profile_successful_starts_out = cf.profile_successful_starts;
            profile_near_optimal_count_out = cf.profile_near_optimal_count;
            profile_best_ll_out = cf.profile_best_ll;
            profile_second_best_ll_out = cf.profile_second_best_ll;
            profile_near_optimal_l1_spread_out = cf.profile_near_optimal_l1_spread;
            multistart_attempted_out = cf.multistart_attempted;
            multistart_configured_starts_out = cf.multistart_configured_starts;
            multistart_successful_starts_out = cf.multistart_successful_starts;
            multistart_near_optimal_count_out = cf.multistart_near_optimal_count;
            multistart_best_ll_out = cf.multistart_best_ll;
            multistart_second_best_ll_out = cf.multistart_second_best_ll;
            multistart_near_optimal_l1_spread_out = cf.multistart_near_optimal_l1_spread;

            delta = 0;
        }
        else{
            last_iter_accepted = false;
            if (llprev == 0 || ll > llprev){
                assn = cf.assn;
                assn_llr = cf.assn_llr;
                contam_prof = cf.contam_prof;
                final_loglik_out = cf.final_loglik;
                final_loglik_valid_out = cf.final_loglik_valid;
                contam_rate = cf.contam_rate;
                contam_rate_se = cf.contam_rate_se;
                allele_ratio = cf.allele_ratio;
                allele_ratio_se = cf.allele_ratio_se;
                species_contam_prof_out = cf.species_contam_prof;
                species_contam_prof_conc_out = cf.species_contam_prof_conc;
                profile_successful_starts_out = cf.profile_successful_starts;
                profile_near_optimal_count_out = cf.profile_near_optimal_count;
                profile_best_ll_out = cf.profile_best_ll;
                profile_second_best_ll_out = cf.profile_second_best_ll;
                profile_near_optimal_l1_spread_out = cf.profile_near_optimal_l1_spread;
                multistart_attempted_out = cf.multistart_attempted;
                multistart_configured_starts_out = cf.multistart_configured_starts;
                multistart_successful_starts_out = cf.multistart_successful_starts;
                multistart_near_optimal_count_out = cf.multistart_near_optimal_count;
                multistart_best_ll_out = cf.multistart_best_ll;
                multistart_second_best_ll_out = cf.multistart_second_best_ll;
                multistart_near_optimal_l1_spread_out = cf.multistart_near_optimal_l1_spread;
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

        if (last_iter_accepted){
            cf.write_cell_fit_diagnostics(output_prefix + ".contam_diagnostics.tsv",
                samples, libname, cellranger, seurat, underscore);
            cf.write_cell_source_profile(output_prefix + ".cell_source_profile.tsv",
                samples, libname, cellranger, seurat, underscore);
            cf.write_class_residual_report(output_prefix + ".class_residuals.tsv", samples);
            if (!surface_output_file.empty()){
                cf.write_r_c_likelihood_surface(surface_output_file,
                    surface_selector_barcodes, samples, libname, cellranger,
                    seurat, underscore, condition_key, synthetic_id);
            }
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
    // Write the full-model objective to disk.
    //
    // Running one condition with a fitted profile and one with a supplied
    // profile, on the same unit and the same count bundle, and differencing
    // this value tells whether the likelihood maximum sits at the supplied
    // profile.
    //
    // ⚠️ Both conditions must use the same prior mode. compute_ll() adds a
    // per-cell beta prior term, so a fitted run under --contam_prior_mode none
    // and a fixed run under any other mode are sums over different objectives
    // and are not comparable.
    //
    // ⚠️ The profile is fit over compile_amb_prof_dat rows and this is
    // evaluated over compile_data rows, so a fitted profile has no obligation
    // to beat a fixed one here. A negative difference is a finding.
    {
        string fname = output_prefix + ".model_fit.tsv";
        FILE* outf = fopen(fname.c_str(), "w");
        if (outf == NULL){
            fprintf(stderr, "ERROR: could not write model fit summary %s\n", fname.c_str());
            exit(1);
        }
        fprintf(outf, "metric\tvalue\n");
        fprintf(outf, "tool_version\t%s\n", QC_VERSION.c_str());
        fprintf(outf, "panel_mode\t%s\n", use_interspecies ? "interspecies" : "interindividual");
        fprintf(outf, "assignment_update_mode\t%s\n",
            assignment_update_mode_name(run_once, freeze_assignments).c_str());
        fprintf(outf, "run_once\t%s\n", run_once ? "true" : "false");
        fprintf(outf, "freeze_assignments\t%s\n",
            freeze_assignments ? "true" : "false");
        fprintf(outf, "final_loglik\t%.10f\n", final_loglik_out);
        fprintf(outf, "final_loglik_valid\t%s\n", final_loglik_valid_out ? "true" : "false");
        fprintf(outf, "ambient_profile_fixed\t%s\n", fixed_ambient_file.empty() ? "false" : "true");
        fprintf(outf, "candidate_keyed_rows\t%s\n", candidate_keyed_rows ? "true" : "false");
        fprintf(outf, "candidate_keyed_split\t%.6f\n", candidate_keyed_split);
        fprintf(outf, "n_ambient_candidates\t%lu\n", (unsigned long)contam_prof.size());
        fclose(outf);
        fprintf(stderr, "Wrote model fit summary to %s\n", fname.c_str());
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

    // Write ambient-profile optimizer stability diagnostics.  These values make
    // flat/non-unique mixture fits explicit instead of presenting one arbitrary
    // simplex point as a uniquely identified biological source profile.
    {
        string fname = output_prefix + ".profile_fit_diagnostics.tsv";
        FILE* outf = fopen(fname.c_str(), "w");
        if (outf == NULL){
            fprintf(stderr, "ERROR: could not open profile diagnostics output: %s\n", fname.c_str());
            exit(1);
        }
        // Columns 1-6 keep their existing names, order, and meaning: they
        // describe the final (refinement) profile solve exactly as before.
        //
        // COLUMN 7 CHANGES MEANING. profile_nonunique_flag previously derived
        // from profile_near_optimal_count_out, which the refinement path
        // hardcodes to 1, so the "> 1" test could never fire and the column was
        // a constant 0. It now derives from the multistart, which is the only
        // solve that compares starting points, and can therefore read 1.
        // Consumers diffing this file across builds must expect column 7 to
        // differ. Estimate outputs (.contam_prof, .contam_rate, .allele_ratio)
        // are unaffected by this change and remain byte-identical.
        //
        // Columns 8-14 are new and describe the multistart solve.
        //
        // Column 1 (profile_restarts, the CLI request) and column 9
        // (multistart_configured_starts, what the loop actually stood up)
        // should agree. A disagreement means --profile_restarts is not
        // reaching the trial count, so warn rather than let it pass silently.
        fprintf(outf, "configured_starts\tsuccessful_starts\tnear_optimal_starts"
            "\tbest_log_likelihood\tsecond_best_log_likelihood\tnear_optimal_l1_spread"
            "\tprofile_nonunique_flag\tmultistart_attempted\tmultistart_configured_starts"
            "\tmultistart_successful_starts\tmultistart_near_optimal_starts"
            "\tmultistart_best_log_likelihood\tmultistart_second_best_log_likelihood"
            "\tmultistart_near_optimal_l1_spread\n");
        if (multistart_attempted_out && profile_restarts > 0 &&
            multistart_configured_starts_out != profile_restarts){
            fprintf(stderr, "WARNING: requested --profile_restarts %d but the ambient "
                "profile multistart configured %d starts; the restart flag is not "
                "reaching the trial count.\n",
                profile_restarts, multistart_configured_starts_out);
        }
        bool nonunique = multistart_attempted_out &&
            multistart_near_optimal_count_out > 1 &&
            std::isfinite(multistart_near_optimal_l1_spread_out) &&
            multistart_near_optimal_l1_spread_out > 0.10;
        fprintf(outf, "%d\t%d\t%d\t%.17g\t%.17g\t%.17g\t%d"
            "\t%d\t%d\t%d\t%d\t%.17g\t%.17g\t%.17g\n",
            profile_restarts, profile_successful_starts_out, profile_near_optimal_count_out,
            profile_best_ll_out, profile_second_best_ll_out,
            profile_near_optimal_l1_spread_out, nonunique ? 1 : 0,
            multistart_attempted_out ? 1 : 0,
            multistart_configured_starts_out, multistart_successful_starts_out,
            multistart_near_optimal_count_out, multistart_best_ll_out,
            multistart_second_best_ll_out, multistart_near_optimal_l1_spread_out);
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
       {"fix_r", required_argument, 0, 2101},
       {"adaptive_prior", no_argument, 0, 2007},
       {"condf", required_argument, 0, 2008},
       {"species_condf", required_argument, 0, 2009},
       {"thorough_multistart", no_argument, 0, 2010},
       {"no_adaptive_multistart", no_argument, 0, 2011},
       {"ambient_candidates", required_argument, 0, 2013},
       {"profile_restarts", required_argument, 0, 2014},
       {"leave_source_out", no_argument, 0, 2015},
       {"heterotypic_start_mode", required_argument, 0, 2016},
       {"heterotypic_start_top_k", required_argument, 0, 2017},
       {"no_profile_intervals", no_argument, 0, 2018},
       {"contam_prior_mode", required_argument, 0, 2019},
       {"contam_prior_min_cells", required_argument, 0, 2020},
       {"contam_prior_max_ci_width", required_argument, 0, 2021},
       {"contam_prior_min_weight", required_argument, 0, 2022},
       {"contam_prior_max_gradient", required_argument, 0, 2023},
       {"strict_condf", no_argument, 0, 2024},
       {"run_class", required_argument, 0, 2025},
       {"run_contract", required_argument, 0, 2026},
       {"assignments_basis", required_argument, 0, 2027},
       {"expected_lines_basis", required_argument, 0, 2028},
       {"ambient_candidates_basis", required_argument, 0, 2029},
       {"warm_start_basis", required_argument, 0, 2030},
       {"fixed_ambient_basis", required_argument, 0, 2031},
       {"fix_r_basis", required_argument, 0, 2032},
       {"candidate_keyed_rows", no_argument, 0, 2033},
       {"candidate_keyed_split", required_argument, 0, 2034},
       {"source_exclusion_strength", required_argument, 0, 2035},
       {"r_c_surface_selector", required_argument, 0, 2036},
       {"r_c_surface_out", required_argument, 0, 2037},
       {"condition_key", required_argument, 0, 2038},
       {"synthetic_id", required_argument, 0, 2039},
       {"profile_holdout_barcodes", required_argument, 0, 2040},
       {"profile_holdout_basis", required_argument, 0, 2041},
       {"freeze_assignments", no_argument, 0, 2042},
       {0, 0, 0, 0}
    };

    // Set default values
    // O53 step 2: candidate-keyed ambient-profile rows, off by default so the
    // historical design is reproduced byte for byte unless requested.
    bool candidate_keyed_rows = false;
    double candidate_keyed_split = 0.05;
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
    bool freeze_assignments = false;
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
    string ambient_candidates_file = "";

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
    double source_exclusion_strength = 0.0;
    bool source_exclusion_explicit = false;
    string r_c_surface_selector_file = "";
    string r_c_surface_out_file = "";
    string condition_key = "";
    string synthetic_id = "";
    string profile_holdout_barcodes_file = "";
    string profile_holdout_basis = "unspecified";
    bool r_feedback = false;
    string fix_r_file = "";
    bool adaptive_prior = false;
    bool thorough_multistart = false;
    bool adaptive_multistart = true;
    int profile_restarts = -1;
    string heterotypic_start_mode = "topk";
    int heterotypic_start_top_k = 5;
    bool adaptive_profile_intervals = true;
    string contam_prior_mode = "fusion";
    int contam_prior_min_cells = 20;
    double contam_prior_max_ci_width = 0.50;
    double contam_prior_min_weight = 10.0;
    double contam_prior_max_gradient = 1e-3;
    bool contam_prior_max_gradient_supplied = false;

    // Explicit condf/species_condf paths (override prefix-derived defaults)
    string user_condf_file = "";
    string user_species_condf_file = "";
    bool strict_condf = false;

    // Run/input contract.  Existing users remain in unspecified mode; the
    // benchmark explicitly requests production or oracle and supplies provenance.
    string run_class = "unspecified";
    string run_contract_file = "";
    string assignments_basis = "unspecified";
    string expected_lines_basis = "unspecified";
    string ambient_candidates_basis = "unspecified";
    string warm_start_basis = "unspecified";
    string fixed_ambient_basis = "unspecified";
    string fix_r_basis = "unspecified";

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
                fprintf(stderr, "WARNING: --loo is deprecated; use --leave_source_out.\n");
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
            case 2101:
                fix_r_file = optarg;
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
            case 2013:
                ambient_candidates_file = optarg;
                break;
            case 2014:
                profile_restarts = atoi(optarg);
                if (profile_restarts < 1){
                    fprintf(stderr, "ERROR: --profile_restarts must be >= 1.\n");
                    exit(1);
                }
                break;
            case 2015:
                leave_one_out = true;
                break;
            case 2016:
                heterotypic_start_mode = optarg;
                if (heterotypic_start_mode != "single" && heterotypic_start_mode != "topk" &&
                    heterotypic_start_mode != "exhaustive") {
                    fprintf(stderr, "ERROR: --heterotypic_start_mode must be single, topk, or exhaustive.\n");
                    exit(1);
                }
                break;
            case 2017:
                heterotypic_start_top_k = atoi(optarg);
                if (heterotypic_start_top_k < 1 || heterotypic_start_top_k > 25){
                    fprintf(stderr, "ERROR: --heterotypic_start_top_k must be 1..25.\n");
                    exit(1);
                }
                break;
            case 2018:
                adaptive_profile_intervals = false;
                break;
            case 2019:
                contam_prior_mode = optarg;
                if (contam_prior_mode != "none" && contam_prior_mode != "global" &&
                    contam_prior_mode != "heterotypic" && contam_prior_mode != "fusion") {
                    fprintf(stderr, "ERROR: --contam_prior_mode must be none, global, heterotypic, or fusion.\n");
                    exit(1);
                }
                break;
            case 2020:
                contam_prior_min_cells = atoi(optarg);
                if (contam_prior_min_cells < 2){
                    fprintf(stderr, "ERROR: --contam_prior_min_cells must be >=2.\n");
                    exit(1);
                }
                break;
            case 2021:
                contam_prior_max_ci_width = atof(optarg);
                break;
            case 2022:
                contam_prior_min_weight = atof(optarg);
                break;
            case 2023:
                contam_prior_max_gradient = atof(optarg);
                contam_prior_max_gradient_supplied = true;
                break;
            case 2024:
                strict_condf = true;
                break;
            case 2025:
                run_class = optarg;
                break;
            case 2026:
                run_contract_file = optarg;
                break;
            case 2027:
                assignments_basis = optarg;
                break;
            case 2028:
                expected_lines_basis = optarg;
                break;
            case 2029:
                ambient_candidates_basis = optarg;
                break;
            case 2030:
                warm_start_basis = optarg;
                break;
            case 2031:
                fixed_ambient_basis = optarg;
                break;
            case 2032:
                fix_r_basis = optarg;
                break;
            case 2033:
                candidate_keyed_rows = true;
                break;
            case 2034:
                candidate_keyed_split = atof(optarg);
                break;
            case 2035: {
                errno = 0;
                char* parse_end = NULL;
                const double parsed_strength = strtod(optarg, &parse_end);
                if (parse_end == optarg || parse_end == NULL || *parse_end != '\0' ||
                    errno == ERANGE || !std::isfinite(parsed_strength)){
                    fprintf(stderr,
                        "ERROR: --source_exclusion_strength requires a finite numeric value in [0,1]; received %s.\n",
                        optarg);
                    exit(1);
                }
                source_exclusion_strength = parsed_strength;
                source_exclusion_explicit = true;
                break;
            }
            case 2036:
                r_c_surface_selector_file = optarg;
                break;
            case 2037:
                r_c_surface_out_file = optarg;
                break;
            case 2038:
                condition_key = optarg;
                break;
            case 2039:
                synthetic_id = optarg;
                break;
            case 2040:
                profile_holdout_barcodes_file = optarg;
                break;
            case 2041:
                profile_holdout_basis = optarg;
                break;
            case 2042:
                freeze_assignments = true;
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
    validate_assignment_update_mode(run_once, freeze_assignments);
    if (source_exclusion_explicit){
        if (!std::isfinite(source_exclusion_strength) ||
            source_exclusion_strength < 0.0 || source_exclusion_strength > 1.0){
            fprintf(stderr, "ERROR: --source_exclusion_strength must be finite and in [0,1].\n");
            exit(1);
        }
        if (leave_one_out && fabs(source_exclusion_strength - 1.0) > 1e-12){
            fprintf(stderr, "ERROR: --loo/--leave_source_out conflicts with explicit --source_exclusion_strength %.17g; legacy flags require exactly 1.\n",
                source_exclusion_strength);
            exit(1);
        }
    } else {
        source_exclusion_strength = leave_one_out ? 1.0 : 0.0;
    }
    if (!r_c_surface_selector_file.empty() != !r_c_surface_out_file.empty()){
        fprintf(stderr, "ERROR: --r_c_surface_selector and --r_c_surface_out must be supplied together.\n");
        exit(1);
    }
    if (!r_c_surface_selector_file.empty() && synthetic_id.empty()){
        fprintf(stderr, "ERROR: --synthetic_id is required with --r_c_surface_selector.\n");
        exit(1);
    }
    if (!r_c_surface_selector_file.empty() && condition_key.empty()){
        fprintf(stderr, "ERROR: --condition_key is required with --r_c_surface_selector.\n");
        exit(1);
    }

    if (!(contam_prior_max_ci_width > 0.0 && contam_prior_max_ci_width <= 1.0)){
        fprintf(stderr, "ERROR: --contam_prior_max_ci_width must be in (0,1].\n");
        exit(1);
    }
    if (contam_prior_min_weight < 0.0){
        fprintf(stderr, "ERROR: --contam_prior_min_weight must be >=0.\n");
        exit(1);
    }
    if (contam_prior_max_gradient_supplied){
        fprintf(stderr, "WARNING: --contam_prior_max_gradient is deprecated and ignored; "
            "prior training now uses scale-invariant profile-likelihood evidence.\n");
    }
    const set<string> valid_run_classes = {"unspecified", "production", "oracle"};
    const set<string> valid_bases = {"unspecified", "library", "truth"};
    if (valid_run_classes.count(run_class) == 0){
        fprintf(stderr, "ERROR: --run_class must be unspecified, production, or oracle.\n");
        exit(1);
    }
    for (const auto& named_basis : vector<pair<string,string> >{
            {"assignments_basis", assignments_basis},
            {"expected_lines_basis", expected_lines_basis},
            {"ambient_candidates_basis", ambient_candidates_basis},
            {"warm_start_basis", warm_start_basis},
            {"fixed_ambient_basis", fixed_ambient_basis},
            {"fix_r_basis", fix_r_basis},
            {"profile_holdout_basis", profile_holdout_basis}}){
        if (valid_bases.count(named_basis.second) == 0){
            fprintf(stderr, "ERROR: --%s must be unspecified, library, or truth.\n",
                named_basis.first.c_str());
            exit(1);
        }
    }
    bool production_contract_pass = true;
    string production_contract_reason = "not_requested";
    if (run_class == "production"){
        vector<string> violations;
        if (!fixed_ambient_file.empty() && fixed_ambient_basis != "library")
            violations.push_back("fixed_ambient_basis_not_library");
        if (!fix_r_file.empty() && fix_r_basis != "library")
            violations.push_back("fix_r_basis_not_library");
        if (assignments_basis != "library") violations.push_back("assignments_basis_not_library");
        if (!expected_lines_file.empty() && expected_lines_basis != "library")
            violations.push_back("expected_lines_basis_not_library");
        if (!ambient_candidates_file.empty() && ambient_candidates_basis != "library")
            violations.push_back("ambient_candidates_basis_not_library");
        if (!warm_start_file.empty() && warm_start_basis != "library")
            violations.push_back("warm_start_basis_not_library");
        if (!profile_holdout_barcodes_file.empty() && profile_holdout_basis != "library")
            violations.push_back("profile_holdout_basis_not_library");
        if (!strict_condf) violations.push_back("strict_condf_required");
        if (!violations.empty()){
            production_contract_pass = false;
            production_contract_reason.clear();
            for (size_t i = 0; i < violations.size(); ++i){
                if (i) production_contract_reason += ";";
                production_contract_reason += violations[i];
            }
            fprintf(stderr, "ERROR: production input contract failed: %s\n",
                production_contract_reason.c_str());
            exit(1);
        }
        production_contract_reason = "ok";
    } else if (run_class == "oracle"){
        production_contract_pass = false;
        production_contract_reason = "oracle_run_not_eligible_for_production_scoring";
    }
    if (output_prefix == ""){
        fprintf(stderr, "ERROR: output_prefix required\n");
        exit(1);
    }
    const string effective_fixed_r_basis = fix_r_file.empty() ? "none" : fix_r_basis;
    const string effective_fixed_ambient_basis = fixed_ambient_file.empty() ? "none" : fixed_ambient_basis;
    const bool truth_assisted_condition = run_class == "oracle";
    const set<string> r_c_surface_barcodes = load_r_c_surface_selector(
        r_c_surface_selector_file, synthetic_id);
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

    const set<unsigned long> profile_holdout_barcodes =
        load_profile_holdout_barcodes(profile_holdout_barcodes_file);
    if (!profile_holdout_barcodes.empty()){
        if (use_interspecies){
            fprintf(stderr, "ERROR: --profile_holdout_barcodes is currently supported only "
                "for the interindividual CK path.\n");
            exit(1);
        }
        unsigned long missing = 0;
        for (set<unsigned long>::const_iterator it = profile_holdout_barcodes.begin();
            it != profile_holdout_barcodes.end(); ++it){
            if (assn.count(*it) == 0) missing++;
        }
        if (missing > 0){
            fprintf(stderr, "ERROR: --profile_holdout_barcodes contains %lu barcodes "
                "not present after assignment/LLR filtering.\n", missing);
            exit(1);
        }
        if (profile_holdout_barcodes.size() >= assn.size()){
            fprintf(stderr, "ERROR: profile holdout contains all %lu scored cells; "
                "profile training would be empty.\n", (unsigned long)assn.size());
            exit(1);
        }
        if (source_exclusion_strength > 1e-12){
            fprintf(stderr, "ERROR: receiver-data profile cross-fitting must retain all source "
                "columns; combine --profile_holdout_barcodes only with "
                "--source_exclusion_strength 0.\n");
            exit(1);
        }
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
    if (candidate_keyed_split < 0.0 || candidate_keyed_split > 1.0){
        fprintf(stderr, "ERROR: --candidate_keyed_split must lie in [0, 1]; got %f\n",
            candidate_keyed_split);
        exit(1);
    }

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

    set<int> ambient_candidate_ids;
    if (!ambient_candidates_file.empty()){
        ifstream acf(ambient_candidates_file);
        if (!acf.is_open()){
            fprintf(stderr, "ERROR: cannot open ambient candidate file: %s\n", ambient_candidates_file.c_str());
            exit(1);
        }
        map<string, int> name_to_idx;
        for (int i = 0; i < (int)samples.size(); i++) name_to_idx[samples[i]] = i;
        string line;
        vector<string> unknown_candidates;
        while (getline(acf, line)){
            while (!line.empty() && (line.back() == '\r' || line.back() == '\n' || line.back() == ' ' || line.back() == '\t')) line.pop_back();
            size_t start = line.find_first_not_of(" \t");
            if (start == string::npos) continue;
            line = line.substr(start);
            if (line.empty() || line[0] == '#') continue;
            if (line.find('+') != string::npos){
                fprintf(stderr, "ERROR: ambient candidate entries must be singlet source IDs, not combinations: %s\n", line.c_str());
                exit(1);
            }
            auto it = name_to_idx.find(line);
            if (it == name_to_idx.end()) unknown_candidates.push_back(line);
            else ambient_candidate_ids.insert(it->second);
        }
        if (!unknown_candidates.empty()){
            fprintf(stderr, "ERROR: %lu ambient candidates were not present in the active sample universe:\n", unknown_candidates.size());
            for (const string& x : unknown_candidates) fprintf(stderr, "  %s\n", x.c_str());
            exit(1);
        }
        if (ambient_candidate_ids.empty()){
            fprintf(stderr, "ERROR: ambient candidate file %s contained no valid source IDs.\n", ambient_candidates_file.c_str());
            exit(1);
        }
        fprintf(stderr, "Loaded %lu explicit ambient source candidates from %s\n",
            ambient_candidate_ids.size(), ambient_candidates_file.c_str());
    }

    // ========================================================================
    // Check for existing outputs / run estimation
    // ========================================================================
    string prof_name = output_prefix + ".contam_prof";
    string rate_name = output_prefix + ".contam_rate";

    robin_hood::unordered_map<unsigned long, double> contam_rate;
    map<int, double> contam_prof;

    bool load_gex = (barcodesfile != "" && featuresfile != "" && matrixfile != "");
    bool genotype_inference_ran = false;

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
        genotype_inference_ran = true;
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
            freeze_assignments,
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
            ambient_candidate_ids,
            leave_one_out,
            source_exclusion_strength,
            r_feedback,
            fix_r_file,
            effective_fixed_r_basis,
            effective_fixed_ambient_basis,
            truth_assisted_condition,
            r_c_surface_barcodes,
            r_c_surface_out_file,
            condition_key,
            synthetic_id,
            adaptive_prior,
            thorough_multistart,
            adaptive_multistart,
            profile_restarts,
            heterotypic_start_mode,
            heterotypic_start_top_k,
            adaptive_profile_intervals,
            contam_prior_mode,
            contam_prior_min_cells,
            contam_prior_max_ci_width,
            contam_prior_min_weight,
            contam_prior_max_gradient,
            warm_start_file,
            fixed_ambient_file,
            species_regularize,
            panel_metadata_file,
            species_counts_file,
            use_interspecies,
            user_condf_file,
            user_species_condf_file,
            strict_condf,
            species_composition_overrides,
            profile_holdout_barcodes,
            profile_holdout_basis,
            candidate_keyed_rows,
            candidate_keyed_split);
    }

    if (genotype_inference_ran){
        if (run_contract_file.empty()) run_contract_file = output_prefix + ".run_contract.json";
        const string counts_contract_path = output_prefix + (use_interspecies ? ".species_counts" : ".counts");
        const string condf_contract_path = !use_interspecies
            ? (user_condf_file.empty() ? output_prefix + ".condf" : user_condf_file)
            : (user_species_condf_file.empty() ? output_prefix + ".species_condf" : user_species_condf_file);
        write_run_contract_json(
            run_contract_file, output_prefix, run_class, production_contract_pass,
            production_contract_reason, use_interspecies, counts_contract_path,
            condf_contract_path, assn_name, sample_name, expected_lines_file,
            ambient_candidates_file, warm_start_file, fixed_ambient_file, fix_r_file,
            assignments_basis, expected_lines_basis, ambient_candidates_basis,
            warm_start_basis, effective_fixed_ambient_basis, effective_fixed_r_basis,
            strict_condf, run_once, freeze_assignments,
            condition_key, synthetic_id, source_exclusion_strength,
            source_exclusion_explicit, profile_holdout_barcodes_file,
            profile_holdout_basis, (unsigned long)profile_holdout_barcodes.size(),
            r_c_surface_selector_file, r_c_surface_out_file);
    }
    else{
        fprintf(stderr, "Existing contamination estimates reused for GEX processing; "
            "estimator run contract was not rewritten.\n");
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
