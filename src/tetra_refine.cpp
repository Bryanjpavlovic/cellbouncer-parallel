// tetra_refine v3.2
//
// Revision history:
//   V2_R1 - Complete rewrite separating ploidy from droplet contents (chat 543ca616)
//   V3_R1 - Add DNLLR, margin_ratio, posterior/entropy carry-through,
//           demote het_var to annotation, add --contam_rate and --external_ploidy inputs,
//           backward-compatible 11-col and 13-col diagnostics parsing

#include <getopt.h>
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
#include <map>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <float.h>
#include <zlib.h>
#include <htswrapper/gzreader.h>
#include <iomanip>

using std::cout;
using std::endl;
using namespace std;

// ============================================================================
// TERMINOLOGY
// ============================================================================
// 
// PLOIDY (property of a CELL):
//   - diploid: cell has 2 copies of genome (normal cell)
//   - tetraploid: cell has 4 copies of genome (fused cell)
//     - heterotypic tetraploid: fused from two different individuals (A x B)
//     - homotypic tetraploid: fused from same individual (A x A)
//
// DROPLET CONTENTS:
//   - 1 cell in droplet: singlet (can be diploid OR tetraploid)
//   - 2 cells in droplet: doublet (both could be diploid, both tetraploid, or mixed)
//
// demux_parallel OUTPUT:
//   - "S" (singlet): assigned to one individual identity
//   - "D" (doublet): assigned to two individual identities (A+B)
//
// This tool's job:
//   1. Recognize that "D" assignments (A+B) are tetraploid SINGLETS, not doublets
//   2. Flag candidate homotypic tetraploids among "S" assignments
//   3. Detect actual doublets of tetraploids (quads) via runner-up patterns
//   4. Compute cross-library confidence metrics (DNLLR, margin_ratio)
//   5. Carry through posterior/entropy from demux_parallel and contam_rate
//      from quant_contam
// ============================================================================

// ============================================================================
// Data Structures
// ============================================================================

struct Assignment {
    string barcode;
    string identity;
    char type;  // 'S' singlet, 'D' doublet (from demux)
    double llr;
};

struct Diagnostics {
    string barcode;
    string assignment;
    char singlet_doublet;
    double llr;
    double min_margin;
    string worst_competitor;
    int n_close;
    double total_depth;
    double het_balance_var;
    int n_het_sites;
    double het_total_depth;
    // v3: optional columns from updated demux_parallel (13-col format)
    double posterior;   // -1 if not available (11-col format)
    double entropy;     // -1 if not available (11-col format)
};

struct RunnerUp {
    string barcode;
    int rank;
    string identity;
    double llr_vs_winner;
    double min_margin;
};

struct RefinedCell {
    string barcode;
    
    // Original demux output
    string original_assignment;
    char original_type;
    double llr;
    
    // Refined assignment
    string refined_assignment;
    
    // PLOIDY: diploid vs tetraploid (property of the cell)
    string ploidy;              // "diploid", "tetraploid", "unknown"
    string ploidy_method;       // "heterotypic", "expected_lines", "external", "ambiguous", "no_data"
    string ploidy_confidence;   // "HIGH", "MEDIUM", "LOW"
    
    // DROPLET: how many cells in droplet
    int cells_in_droplet;       // 1 = singlet, 2 = doublet
    string droplet_flag;        // "NONE", "POSSIBLE_DOUBLET", "LIKELY_DOUBLET"
    string droplet_candidates;
    
    // Confidence metrics (v3)
    double dnllr;               // LLR / total_depth
    double margin_ratio;        // winner_min_margin / runner_up_rank1_min_margin
    double posterior;            // from demux_parallel (-1 if unavailable)
    double entropy;             // from demux_parallel (-1 if unavailable)
    
    // Contamination (v3)
    double contam_rate;         // from quant_contam (-1 if unavailable)
    
    // Quad detection score (v3.2)
    double quad_pattern_score;  // 0-1, how many expected runner-up combos found (-1 if N/A)
    
    // Het balance annotation (v3: demoted from decision-maker)
    string het_var_signal;      // "suggests_diploid", "suggests_tetraploid", "ambiguous", "no_data"
    double het_balance_var;     // raw value carried through
    
    // Overall
    string overall_confidence;
    bool changed;
};

struct Thresholds {
    double min_margin_threshold;       // strict threshold (LIKELY_DOUBLET)
    double min_margin_possible;        // lenient threshold (POSSIBLE_DOUBLET)
    double close_count_threshold;
    double het_var_diploid;
    double het_var_tetraploid;
    double depth_ratio_threshold;
    string source;
};

struct SummaryStats {
    int total_cells;
    int original_singlets;
    int original_doublets;
    int ploidy_diploid;
    int ploidy_tetraploid;
    int ploidy_unknown;
    int tetraploid_heterotypic;
    int tetraploid_homotypic;
    int tetraploid_from_expected;
    int tetraploid_from_external;
    int droplet_singlet;
    int droplet_possible;
    int droplet_likely;
    int assignment_changed;
    int ploidy_inferred;
    int cells_changed;
};

// ============================================================================
// Expected Lines Parsing
// ============================================================================

struct ExpectedLines {
    set<string> diploid_singlets;
    set<string> tetraploid_heterotypic;
    set<string> tetraploid_homotypic;
    set<string> all_individuals;
    map<string, bool> singlet_has_diploid;
    map<string, bool> singlet_has_homotypic;
};

void parse_expected_lines(const string& filename, ExpectedLines& expected) {
    ifstream inf(filename);
    if (!inf.is_open()) {
        fprintf(stderr, "ERROR: Cannot open expected lines file: %s\n", filename.c_str());
        exit(1);
    }
    
    string line;
    while (getline(inf, line)) {
        if (line.empty() || line[0] == '#') continue;
        // Strip trailing whitespace/CR
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n' || line.back() == ' '))
            line.pop_back();
        if (line.empty()) continue;
        
        size_t plus_pos = line.find('+');
        if (plus_pos != string::npos) {
            string id1 = line.substr(0, plus_pos);
            string id2 = line.substr(plus_pos + 1);
            expected.all_individuals.insert(id1);
            expected.all_individuals.insert(id2);
            if (id1 == id2) {
                expected.tetraploid_homotypic.insert(line);
                expected.singlet_has_homotypic[id1] = true;
            } else {
                expected.tetraploid_heterotypic.insert(line);
            }
        } else {
            expected.diploid_singlets.insert(line);
            expected.all_individuals.insert(line);
            expected.singlet_has_diploid[line] = true;
        }
    }
    
    for (const auto& indv : expected.all_individuals) {
        if (expected.singlet_has_diploid.find(indv) == expected.singlet_has_diploid.end()) {
            expected.singlet_has_diploid[indv] = false;
        }
        if (expected.singlet_has_homotypic.find(indv) == expected.singlet_has_homotypic.end()) {
            expected.singlet_has_homotypic[indv] = false;
        }
    }
}

// ============================================================================
// File Parsing Functions
// ============================================================================

void parse_assignments_file(const string& filename, map<string, Assignment>& assignments) {
    ifstream inf(filename);
    if (!inf.is_open()) {
        fprintf(stderr, "ERROR: Cannot open assignments file: %s\n", filename.c_str());
        exit(1);
    }
    
    string line;
    while (getline(inf, line)) {
        if (line.empty()) continue;
        
        istringstream splitter(line);
        string field;
        Assignment a;
        int idx = 0;
        
        while (getline(splitter, field, '\t')) {
            if (idx == 0) a.barcode = field;
            else if (idx == 1) a.identity = field;
            else if (idx == 2 && field.length() > 0) a.type = field[0];
            else if (idx == 3) a.llr = atof(field.c_str());
            idx++;
        }
        
        if (idx >= 4) {
            assignments[a.barcode] = a;
        }
    }
}

// Backward-compatible: handles both 11-col (v2) and 13-col (v3+) diagnostics
void parse_diagnostics_file(const string& filename, map<string, Diagnostics>& diagnostics,
                            bool& has_posterior_data) {
    gzreader reader(filename);
    bool header_read = false;
    has_posterior_data = false;
    int n_cols_detected = 0;
    
    while (reader.next()) {
        if (!header_read) {
            // Count header columns to detect format
            istringstream hdr(reader.line);
            string tok;
            while (getline(hdr, tok, '\t')) n_cols_detected++;
            if (n_cols_detected >= 13) has_posterior_data = true;
            header_read = true;
            continue;
        }
        
        istringstream splitter(reader.line);
        string field;
        Diagnostics d;
        d.posterior = -1.0;
        d.entropy = -1.0;
        int idx = 0;
        
        while (getline(splitter, field, '\t')) {
            if (idx == 0) d.barcode = field;
            else if (idx == 1) d.assignment = field;
            else if (idx == 2 && field.length() > 0) d.singlet_doublet = field[0];
            else if (idx == 3) d.llr = atof(field.c_str());
            else if (idx == 4) d.min_margin = atof(field.c_str());
            else if (idx == 5) d.worst_competitor = field;
            else if (idx == 6) d.n_close = atoi(field.c_str());
            else if (idx == 7) d.total_depth = atof(field.c_str());
            else if (idx == 8) d.het_balance_var = atof(field.c_str());
            else if (idx == 9) d.n_het_sites = atoi(field.c_str());
            else if (idx == 10) d.het_total_depth = atof(field.c_str());
            else if (idx == 11) d.posterior = atof(field.c_str());
            else if (idx == 12) d.entropy = atof(field.c_str());
            idx++;
        }
        
        if (idx >= 8) {
            diagnostics[d.barcode] = d;
        }
    }
}

void parse_runner_ups_file(const string& filename, map<string, vector<RunnerUp>>& runner_ups) {
    gzreader reader(filename);
    bool header_read = false;
    
    while (reader.next()) {
        if (!header_read) {
            header_read = true;
            continue;
        }
        
        istringstream splitter(reader.line);
        string field;
        RunnerUp ru;
        int idx = 0;
        
        while (getline(splitter, field, '\t')) {
            if (idx == 0) ru.barcode = field;
            else if (idx == 1) ru.rank = atoi(field.c_str());
            else if (idx == 2) ru.identity = field;
            else if (idx == 3) ru.llr_vs_winner = atof(field.c_str());
            else if (idx == 4) ru.min_margin = atof(field.c_str());
            idx++;
        }
        
        if (idx >= 5) {
            runner_ups[ru.barcode].push_back(ru);
        }
    }
}

// Parse quant_contam .contam_rate file: barcode \t rate \t se
void parse_contam_rate_file(const string& filename, map<string, double>& contam_rates) {
    ifstream inf(filename);
    if (!inf.is_open()) {
        fprintf(stderr, "ERROR: Cannot open contam_rate file: %s\n", filename.c_str());
        exit(1);
    }
    
    string line;
    while (getline(inf, line)) {
        if (line.empty()) continue;
        istringstream splitter(line);
        string barcode, rate_str;
        if (getline(splitter, barcode, '\t') && getline(splitter, rate_str, '\t')) {
            contam_rates[barcode] = atof(rate_str.c_str());
        }
    }
    
    fprintf(stderr, "  Loaded contamination rates for %lu cells\n", contam_rates.size());
}

// Parse external ploidy file: barcode \t ploidy_call [\t ploidy_probability]
struct ExternalPloidy {
    string ploidy_call;   // "diploid" or "tetraploid"
    double probability;   // confidence, -1 if not provided
};

void parse_external_ploidy_file(const string& filename, map<string, ExternalPloidy>& ext_ploidy,
                                const string& filter_library) {
    ifstream inf(filename);
    if (!inf.is_open()) {
        fprintf(stderr, "ERROR: Cannot open external ploidy file: %s\n", filename.c_str());
        exit(1);
    }
    
    string line;
    bool header_skipped = false;
    int lib_col = -1;  // column index for library (-1 = not present)
    int call_col = 1;  // default: barcode \t call \t prob (old format)
    int prob_col = 2;
    size_t total_rows = 0;
    size_t filtered_rows = 0;
    
    while (getline(inf, line)) {
        if (line.empty()) continue;
        // Parse header to detect column layout
        if (!header_skipped && (line.find("barcode") != string::npos || line.find("Barcode") != string::npos)) {
            header_skipped = true;
            // Detect column positions from header
            istringstream hsplit(line);
            string col;
            int idx = 0;
            while (getline(hsplit, col, '\t')) {
                if (col == "library") lib_col = idx;
                else if (col == "ploidy_call") call_col = idx;
                else if (col == "ploidy_probability") prob_col = idx;
                idx++;
            }
            continue;
        }
        header_skipped = true;
        
        // Split line into fields
        vector<string> fields;
        istringstream splitter(line);
        string field;
        while (getline(splitter, field, '\t')) {
            fields.push_back(field);
        }
        
        if ((int)fields.size() <= call_col) continue;
        total_rows++;
        
        // Filter by library if requested
        if (!filter_library.empty() && lib_col >= 0 && lib_col < (int)fields.size()) {
            if (fields[lib_col] != filter_library) {
                filtered_rows++;
                continue;
            }
        }
        
        string barcode = fields[0];
        ExternalPloidy ep;
        ep.ploidy_call = fields[call_col];
        ep.probability = (prob_col < (int)fields.size()) ? atof(fields[prob_col].c_str()) : -1.0;
        
        ext_ploidy[barcode] = ep;
    }
    
    if (!filter_library.empty() && lib_col >= 0) {
        fprintf(stderr, "  Loaded external ploidy: %lu cells for library %s (filtered %lu of %lu total)\n",
                ext_ploidy.size(), filter_library.c_str(), filtered_rows, total_rows);
    } else if (!filter_library.empty() && lib_col < 0) {
        fprintf(stderr, "  WARNING: --external_ploidy_library %s specified but no 'library' column found in ploidy file\n",
                filter_library.c_str());
        fprintf(stderr, "  Loaded all %lu external ploidy calls (no filtering)\n", ext_ploidy.size());
    } else {
        fprintf(stderr, "  Loaded external ploidy calls for %lu cells\n", ext_ploidy.size());
    }
}

// ============================================================================
// Statistical Functions
// ============================================================================

double compute_percentile(vector<double>& values, double percentile) {
    if (values.empty()) return 0.0;
    sort(values.begin(), values.end());
    double idx = (percentile / 100.0) * (values.size() - 1);
    size_t lower = (size_t)floor(idx);
    size_t upper = (size_t)ceil(idx);
    if (lower == upper || upper >= values.size()) return values[lower];
    double frac = idx - lower;
    return values[lower] * (1.0 - frac) + values[upper] * frac;
}

double compute_median(vector<double>& values) {
    return compute_percentile(values, 50.0);
}

// ============================================================================
// Identity Parsing Helpers
// ============================================================================

bool parse_doublet_identity(const string& identity, string& id1, string& id2) {
    size_t plus_pos = identity.find('+');
    if (plus_pos == string::npos) return false;
    id1 = identity.substr(0, plus_pos);
    id2 = identity.substr(plus_pos + 1);
    return true;
}

bool is_homotypic(const string& identity) {
    string id1, id2;
    if (!parse_doublet_identity(identity, id1, id2)) return false;
    return id1 == id2;
}

string make_homotypic(const string& singlet) {
    return singlet + "+" + singlet;
}

// ============================================================================
// Confidence Metric Computation
// ============================================================================

double compute_dnllr(double llr, double total_depth) {
    if (total_depth <= 0) return -1.0;
    return llr / total_depth;
}

double compute_margin_ratio(double winner_min_margin, const vector<RunnerUp>* runner_ups) {
    if (runner_ups == NULL || runner_ups->empty()) return 1e6;
    // rank 1 runner-up (first entry, should be sorted by rank)
    double ru_margin = (*runner_ups)[0].min_margin;
    if (ru_margin <= 0) return 1e6;  // runner-up has no support
    return winner_min_margin / ru_margin;
}

// ============================================================================
// Het Var Signal Annotation (v3: no longer drives assignments)
// ============================================================================

string compute_het_var_signal(double het_balance_var, const Thresholds& thresholds, bool has_het_data) {
    if (het_balance_var < 0 || !has_het_data) return "no_data";
    if (thresholds.het_var_diploid <= 0 && thresholds.het_var_tetraploid <= 0) return "no_data";
    
    double threshold = (thresholds.het_var_diploid + thresholds.het_var_tetraploid) / 2.0;
    if (threshold <= 0) return "no_data";
    
    double dist = fabs(het_balance_var - threshold);
    double relative_dist = dist / threshold;
    
    if (het_balance_var < threshold) {
        return (relative_dist > 0.3) ? "suggests_tetraploid" : "ambiguous";
    } else {
        return (relative_dist > 0.3) ? "suggests_diploid" : "ambiguous";
    }
}

// ============================================================================
// Droplet Detection (Quad Problem)
// ============================================================================

bool sets_are_disjoint(const set<string>& s1, const set<string>& s2) {
    for (const auto& elem : s1) {
        if (s2.count(elem) > 0) return false;
    }
    return true;
}

string make_doublet_key(const string& id1, const string& id2) {
    if (id1 < id2) return id1 + "+" + id2;
    return id2 + "+" + id1;
}

string detect_droplet_doublet(const Diagnostics& diag,
                              const vector<RunnerUp>& runner_ups,
                              const Thresholds& thresholds,
                              string& candidates,
                              double& out_pattern_score) {
    candidates = ".";
    out_pattern_score = -1.0;
    
    if (diag.singlet_doublet != 'D') return "NONE";
    
    string winner_id1, winner_id2;
    if (!parse_doublet_identity(diag.assignment, winner_id1, winner_id2)) return "NONE";
    if (winner_id1 == winner_id2) return "NONE";
    
    set<string> winner_ids;
    winner_ids.insert(winner_id1);
    winner_ids.insert(winner_id2);
    
    set<string> runner_up_keys;
    map<string, set<string>> runner_up_components;
    
    for (const auto& ru : runner_ups) {
        string ru_id1, ru_id2;
        if (parse_doublet_identity(ru.identity, ru_id1, ru_id2)) {
            string key = make_doublet_key(ru_id1, ru_id2);
            runner_up_keys.insert(key);
            set<string> components;
            components.insert(ru_id1);
            components.insert(ru_id2);
            runner_up_components[key] = components;
        }
    }
    
    string best_complement_key = "";
    set<string> complement_ids;
    
    for (const auto& ru_pair : runner_up_components) {
        if (sets_are_disjoint(winner_ids, ru_pair.second)) {
            best_complement_key = ru_pair.first;
            complement_ids = ru_pair.second;
            break;
        }
    }
    
    if (best_complement_key.empty()) return "NONE";
    
    int cross_combos_found = 0;
    for (const auto& w_id : winner_ids) {
        for (const auto& c_id : complement_ids) {
            string cross_key = make_doublet_key(w_id, c_id);
            if (runner_up_keys.count(cross_key) > 0) cross_combos_found++;
        }
    }
    
    double pattern_score = (1.0 + cross_combos_found) / 5.0;
    out_pattern_score = pattern_score;
    
    vector<string> all_four;
    for (const auto& id : winner_ids) all_four.push_back(id);
    for (const auto& id : complement_ids) all_four.push_back(id);
    sort(all_four.begin(), all_four.end());
    
    candidates = "";
    for (size_t i = 0; i < all_four.size(); i++) {
        if (i > 0) candidates += "+";
        candidates += all_four[i];
    }
    
    bool strict_margin = (diag.min_margin < thresholds.min_margin_threshold);
    bool lenient_margin = (diag.min_margin < thresholds.min_margin_possible);
    bool has_close = (diag.n_close >= 1);
    bool strong_pattern = (pattern_score >= 0.6);  // complement + 2 cross-combos
    bool weak_pattern = (pattern_score >= 0.4);    // complement + 1 cross-combo
    
    // LIKELY: strict margin + close competitors + strong pattern
    if (strict_margin && has_close && strong_pattern) return "LIKELY_DOUBLET";
    
    // POSSIBLE: lenient margin + weak pattern (use quad_pattern_score to filter downstream)
    if (lenient_margin && weak_pattern) return "POSSIBLE_DOUBLET";
    
    // Also flag if has_close + strong pattern even without low margin
    if (has_close && strong_pattern) return "POSSIBLE_DOUBLET";
    
    // Not flagged, but preserve candidates and score for downstream filtering
    // (candidates and out_pattern_score are already set above)
    return "NONE";
}

// ============================================================================
// Ploidy Classification
// ============================================================================

void classify_cell(const Assignment& assn,
                   const Diagnostics* diag,
                   const vector<RunnerUp>* runner_ups,
                   const ExternalPloidy* ext_ploidy,
                   const ExpectedLines& expected,
                   const Thresholds& thresholds,
                   double median_depth,
                   bool has_het_data,
                   RefinedCell& cell) {
    
    cell.barcode = assn.barcode;
    cell.original_assignment = assn.identity;
    cell.original_type = assn.type;
    cell.llr = assn.llr;
    cell.refined_assignment = assn.identity;
    cell.changed = false;
    cell.cells_in_droplet = 1;
    cell.droplet_flag = "NONE";
    cell.droplet_candidates = ".";
    cell.het_var_signal = "no_data";
    cell.het_balance_var = -1.0;
    cell.quad_pattern_score = -1.0;
    
    // Compute het_var annotation (never drives decisions now)
    if (diag != NULL && diag->het_balance_var >= 0 && diag->singlet_doublet == 'S') {
        cell.het_balance_var = diag->het_balance_var;
        cell.het_var_signal = compute_het_var_signal(diag->het_balance_var, thresholds, has_het_data);
    }
    
    // =========================================================================
    // CASE 1: Heterotypic tetraploid (demux assigned A+B where A != B)
    // =========================================================================
    if (assn.type == 'D') {
        string id1, id2;
        parse_doublet_identity(assn.identity, id1, id2);
        
        if (id1 != id2) {
            cell.ploidy = "tetraploid";
            cell.ploidy_method = "heterotypic";
            cell.ploidy_confidence = "HIGH";
            cell.cells_in_droplet = 1;
            
            if (diag != NULL && runner_ups != NULL) {
                cell.droplet_flag = detect_droplet_doublet(*diag, *runner_ups,
                                                           thresholds, cell.droplet_candidates,
                                                           cell.quad_pattern_score);
                if (cell.droplet_flag != "NONE") cell.cells_in_droplet = 2;
            }
            
            cell.overall_confidence = (cell.droplet_flag == "NONE") ? "HIGH" : "LOW";
            return;
        } else {
            // demux output A+A as "D" (shouldn't happen but handle it)
            cell.ploidy = "tetraploid";
            cell.ploidy_method = "homotypic_direct";
            cell.ploidy_confidence = "HIGH";
            cell.cells_in_droplet = 1;
            cell.overall_confidence = "HIGH";
            return;
        }
    }
    
    // =========================================================================
    // CASE 2: Singlet assignment - determine if diploid or homotypic tetraploid
    // =========================================================================
    
    string indv = assn.identity;
    bool has_diploid_in_pool = expected.singlet_has_diploid.count(indv) &&
                                expected.singlet_has_diploid.at(indv);
    bool has_homotypic_in_pool = expected.singlet_has_homotypic.count(indv) &&
                                  expected.singlet_has_homotypic.at(indv);
    
    // -------------------------------------------------------------------------
    // CASE 2a: Only homotypic tetraploid in pool - CHANGE ASSIGNMENT
    // -------------------------------------------------------------------------
    if (has_homotypic_in_pool && !has_diploid_in_pool) {
        cell.refined_assignment = make_homotypic(indv);
        cell.ploidy = "tetraploid";
        cell.ploidy_method = "expected_lines";
        cell.ploidy_confidence = "HIGH";
        cell.cells_in_droplet = 1;
        cell.overall_confidence = "HIGH";
        cell.changed = true;
        return;
    }
    
    // -------------------------------------------------------------------------
    // CASE 2b: Only diploid in pool - confirm as diploid
    // -------------------------------------------------------------------------
    if (has_diploid_in_pool && !has_homotypic_in_pool) {
        cell.ploidy = "diploid";
        cell.ploidy_method = "expected_lines";
        cell.ploidy_confidence = "HIGH";
        cell.cells_in_droplet = 1;
        cell.overall_confidence = "HIGH";
        return;
    }
    
    // -------------------------------------------------------------------------
    // CASE 2c: Both or neither in pool - use external ploidy if available
    // -------------------------------------------------------------------------
    if (ext_ploidy != NULL) {
        cell.ploidy = ext_ploidy->ploidy_call;
        cell.ploidy_method = "external";
        
        if (ext_ploidy->probability >= 0.9) {
            cell.ploidy_confidence = "HIGH";
        } else if (ext_ploidy->probability >= 0.7) {
            cell.ploidy_confidence = "MEDIUM";
        } else if (ext_ploidy->probability > 0) {
            cell.ploidy_confidence = "LOW";
        } else {
            // probability not provided, default MEDIUM
            cell.ploidy_confidence = "MEDIUM";
        }
        
        if (cell.ploidy == "tetraploid" && has_homotypic_in_pool) {
            cell.refined_assignment = make_homotypic(indv);
            cell.changed = true;
        }
        
        cell.cells_in_droplet = 1;
        cell.overall_confidence = cell.ploidy_confidence;
        return;
    }
    
    // -------------------------------------------------------------------------
    // CASE 2d: No external ploidy - mark as unknown (het_var is annotation only)
    // -------------------------------------------------------------------------
    cell.ploidy = "unknown";
    cell.ploidy_method = "ambiguous";
    cell.ploidy_confidence = "LOW";
    cell.cells_in_droplet = 1;
    cell.overall_confidence = "LOW";
}

// ============================================================================
// Output Functions
// ============================================================================

void write_refined_assignments(const string& filename, const vector<RefinedCell>& cells) {
    ofstream outf(filename);
    if (!outf.is_open()) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", filename.c_str());
        exit(1);
    }
    
    outf << "barcode\t"
         << "original_assignment\t"
         << "original_type\t"
         << "refined_assignment\t"
         << "llr\t"
         << "ploidy\t"
         << "ploidy_method\t"
         << "ploidy_confidence\t"
         << "cells_in_droplet\t"
         << "droplet_flag\t"
         << "droplet_candidates\t"
         << "overall_confidence\t"
         << "changed\t"
         << "dnllr\t"
         << "margin_ratio\t"
         << "posterior\t"
         << "entropy\t"
         << "contam_rate\t"
         << "het_var_signal\t"
         << "het_balance_var\t"
         << "quad_pattern_score" << endl;
    
    for (const auto& cell : cells) {
        outf << cell.barcode << "\t"
             << cell.original_assignment << "\t"
             << cell.original_type << "\t"
             << cell.refined_assignment << "\t"
             << cell.llr << "\t"
             << cell.ploidy << "\t"
             << cell.ploidy_method << "\t"
             << cell.ploidy_confidence << "\t"
             << cell.cells_in_droplet << "\t"
             << cell.droplet_flag << "\t"
             << cell.droplet_candidates << "\t"
             << cell.overall_confidence << "\t"
             << (cell.changed ? "TRUE" : "FALSE") << "\t";
        
        // DNLLR
        if (cell.dnllr >= 0) outf << fixed << setprecision(6) << cell.dnllr;
        else outf << ".";
        outf << "\t";
        
        // margin_ratio
        if (cell.margin_ratio < 1e5) outf << fixed << setprecision(2) << cell.margin_ratio;
        else outf << "INF";
        outf << "\t";
        
        // posterior
        if (cell.posterior >= 0) outf << fixed << setprecision(6) << cell.posterior;
        else outf << ".";
        outf << "\t";
        
        // entropy
        if (cell.entropy >= 0) outf << fixed << setprecision(4) << cell.entropy;
        else outf << ".";
        outf << "\t";
        
        // contam_rate
        if (cell.contam_rate >= 0) outf << fixed << setprecision(6) << cell.contam_rate;
        else outf << ".";
        outf << "\t";
        
        // het_var_signal
        outf << cell.het_var_signal << "\t";
        
        // het_balance_var
        if (cell.het_balance_var >= 0) outf << fixed << setprecision(6) << cell.het_balance_var;
        else outf << ".";
        outf << "\t";
        
        // quad_pattern_score
        if (cell.quad_pattern_score >= 0) outf << fixed << setprecision(3) << cell.quad_pattern_score;
        else outf << ".";
        
        outf << endl;
    }
}

void write_simple_assignments(const string& filename, const vector<RefinedCell>& cells) {
    ofstream outf(filename);
    if (!outf.is_open()) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", filename.c_str());
        exit(1);
    }
    
    outf << "barcode\tassignment\tploidy\tcells_in_droplet\tllr" << endl;
    
    for (const auto& cell : cells) {
        outf << cell.barcode << "\t"
             << cell.refined_assignment << "\t"
             << cell.ploidy << "\t"
             << cell.cells_in_droplet << "\t"
             << cell.llr << endl;
    }
}

void write_summary(const string& filename,
                   const string& assignments_file,
                   const string& diagnostics_file,
                   const string& runner_ups_file,
                   const string& expected_lines_file,
                   const string& contam_rate_file,
                   const string& external_ploidy_file,
                   const string& external_ploidy_library,
                   const ExpectedLines& expected,
                   const Thresholds& thresholds,
                   bool has_het_data,
                   bool has_posterior_data,
                   double median_depth,
                   const SummaryStats& stats) {
    
    ofstream outf(filename);
    if (!outf.is_open()) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", filename.c_str());
        exit(1);
    }
    
    outf << "tetra_refine v3.2 summary" << endl;
    outf << "=======================" << endl;
    outf << endl;
    
    outf << "Input files:" << endl;
    outf << "  Assignments: " << assignments_file << endl;
    outf << "  Diagnostics: " << diagnostics_file << endl;
    outf << "  Runner-ups: " << runner_ups_file << endl;
    outf << "  Expected lines: " << expected_lines_file << endl;
    if (!contam_rate_file.empty())
        outf << "  Contam rate: " << contam_rate_file << endl;
    if (!external_ploidy_file.empty()) {
        outf << "  External ploidy: " << external_ploidy_file << endl;
        if (!external_ploidy_library.empty())
            outf << "  External ploidy library filter: " << external_ploidy_library << endl;
    }
    outf << endl;
    
    outf << "Expected lines summary:" << endl;
    outf << "  Diploid singlets: " << expected.diploid_singlets.size() << endl;
    outf << "  Heterotypic tetraploids: " << expected.tetraploid_heterotypic.size() << endl;
    outf << "  Homotypic tetraploids: " << expected.tetraploid_homotypic.size() << endl;
    outf << "  Total individuals: " << expected.all_individuals.size() << endl;
    outf << endl;
    
    outf << "Data availability:" << endl;
    outf << "  Het variance data: " << (has_het_data ? "YES" : "NO") << endl;
    outf << "  Posterior/entropy data: " << (has_posterior_data ? "YES" : "NO") << endl;
    outf << "  Contamination rate data: " << (!contam_rate_file.empty() ? "YES" : "NO") << endl;
    outf << "  External ploidy data: " << (!external_ploidy_file.empty() ? "YES" : "NO") << endl;
    outf << "  Median depth: " << median_depth << endl;
    outf << endl;
    
    outf << "Thresholds:" << endl;
    outf << "  min_margin_threshold (strict): " << thresholds.min_margin_threshold << endl;
    outf << "  min_margin_possible (lenient): " << thresholds.min_margin_possible << endl;
    outf << "  het_var_diploid: " << thresholds.het_var_diploid << " (annotation only)" << endl;
    outf << "  het_var_tetraploid: " << thresholds.het_var_tetraploid << " (annotation only)" << endl;
    outf << "  Source: " << thresholds.source << endl;
    outf << endl;
    
    outf << "Cell counts:" << endl;
    outf << "  Total cells: " << stats.total_cells << endl;
    outf << "  Original singlets (S): " << stats.original_singlets << endl;
    outf << "  Original doublets (D): " << stats.original_doublets << endl;
    outf << endl;
    
    outf << "Ploidy classification:" << endl;
    outf << "  Diploid: " << stats.ploidy_diploid << endl;
    outf << "  Tetraploid: " << stats.ploidy_tetraploid << endl;
    outf << "    - Heterotypic: " << stats.tetraploid_heterotypic << endl;
    outf << "    - Homotypic (from expected_lines): " << stats.tetraploid_from_expected << endl;
    outf << "    - Homotypic (from external ploidy): " << stats.tetraploid_from_external << endl;
    outf << "  Unknown: " << stats.ploidy_unknown << endl;
    outf << endl;
    
    outf << "Droplet classification:" << endl;
    outf << "  Singlet (1 cell): " << stats.droplet_singlet << endl;
    outf << "  Possible doublet: " << stats.droplet_possible << endl;
    outf << "  Likely doublet: " << stats.droplet_likely << endl;
    outf << endl;
    
    outf << "Refinement summary:" << endl;
    outf << "  Cells with changed assignment: " << stats.assignment_changed << endl;
    outf << "  Total cells changed (any field): " << stats.cells_changed << endl;
    double pct = (stats.total_cells > 0) ? 100.0 * stats.cells_changed / stats.total_cells : 0.0;
    outf << "  Percent changed: " << fixed << setprecision(1) << pct << "%" << endl;
}

// ============================================================================
// Help Function
// ============================================================================

void help(int code) {
    fprintf(stderr, "tetra_refine v3.2 [OPTIONS]\n");
    fprintf(stderr, "Refines demux_parallel assignments for tetraploid pools.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "This tool:\n");
    fprintf(stderr, "  1. Recognizes heterotypic tetraploids (A+B) as single cells, not doublets\n");
    fprintf(stderr, "  2. Reclassifies homotypic tetraploids (A+A) from expected_lines\n");
    fprintf(stderr, "  3. Flags potential doublets of tetraploids (quads)\n");
    fprintf(stderr, "  4. Computes cross-library confidence metrics (DNLLR, margin_ratio)\n");
    fprintf(stderr, "  5. Carries through posterior, entropy, and contamination rate\n");
    fprintf(stderr, "  6. Accepts external ploidy classifications for ambiguous cells\n");
    fprintf(stderr, "\n===== REQUIRED =====\n");
    fprintf(stderr, "    --assignments -a FILE    .assignments file from demux_parallel\n");
    fprintf(stderr, "    --diagnostics -d FILE    .diagnostics.gz from demux_parallel (11 or 13 col)\n");
    fprintf(stderr, "    --runner_ups -r FILE     .runner_ups.gz from demux_parallel\n");
    fprintf(stderr, "    --expected -e FILE       Expected lines file (same format as -I flag)\n");
    fprintf(stderr, "    --output -o PREFIX       Output file prefix\n");
    fprintf(stderr, "\n===== OPTIONAL - Additional data =====\n");
    fprintf(stderr, "    --contam_rate FILE       .contam_rate file from quant_contam\n");
    fprintf(stderr, "    --external_ploidy FILE   External ploidy calls (barcode\\tploidy[\\tprob])\n");
    fprintf(stderr, "    --external_ploidy_library STR  Filter external ploidy to this library number\n");
    fprintf(stderr, "\n===== OPTIONAL - Threshold overrides =====\n");
    fprintf(stderr, "    --het_var_diploid F      Het variance annotation threshold for diploid\n");
    fprintf(stderr, "    --het_var_tetraploid F   Het variance annotation threshold for tetraploid\n");
    fprintf(stderr, "    --min_margin F           Min margin threshold for doublet detection\n");
    fprintf(stderr, "\n===== OPTIONAL - Output =====\n");
    fprintf(stderr, "    --write_changed_only     Also write separate file with only changed cells\n");
    fprintf(stderr, "    --write_simple           Also write simple .assignments format file\n");
    fprintf(stderr, "    --verbose -V             Detailed progress messages\n");
    fprintf(stderr, "    --help -h                Display this message and exit\n");
    exit(code);
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char* argv[]) {
    
    static struct option long_options[] = {
        {"assignments", required_argument, 0, 'a'},
        {"diagnostics", required_argument, 0, 'd'},
        {"runner_ups", required_argument, 0, 'r'},
        {"expected", required_argument, 0, 'e'},
        {"output", required_argument, 0, 'o'},
        {"contam_rate", required_argument, 0, 1006},
        {"external_ploidy", required_argument, 0, 1007},
        {"external_ploidy_library", required_argument, 0, 1008},
        {"het_var_diploid", required_argument, 0, 1001},
        {"het_var_tetraploid", required_argument, 0, 1002},
        {"min_margin", required_argument, 0, 1003},
        {"write_changed_only", no_argument, 0, 1004},
        {"write_simple", no_argument, 0, 1005},
        {"verbose", no_argument, 0, 'V'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    string assignments_file = "";
    string diagnostics_file = "";
    string runner_ups_file = "";
    string expected_lines_file = "";
    string output_prefix = "";
    string contam_rate_file = "";
    string external_ploidy_file = "";
    string external_ploidy_library = "";
    
    double het_var_diploid_override = -1.0;
    double het_var_tetraploid_override = -1.0;
    double min_margin_override = -1.0;
    
    bool write_changed_only = false;
    bool write_simple = false;
    bool verbose = false;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1) {
        help(0);
    }
    
    while ((ch = getopt_long(argc, argv, "a:d:r:e:o:Vh", long_options, &option_index)) != -1) {
        switch (ch) {
            case 'a': assignments_file = optarg; break;
            case 'd': diagnostics_file = optarg; break;
            case 'r': runner_ups_file = optarg; break;
            case 'e': expected_lines_file = optarg; break;
            case 'o': output_prefix = optarg; break;
            case 1001: het_var_diploid_override = atof(optarg); break;
            case 1002: het_var_tetraploid_override = atof(optarg); break;
            case 1003: min_margin_override = atof(optarg); break;
            case 1004: write_changed_only = true; break;
            case 1005: write_simple = true; break;
            case 1006: contam_rate_file = optarg; break;
            case 1007: external_ploidy_file = optarg; break;
            case 1008: external_ploidy_library = optarg; break;
            case 'V': verbose = true; break;
            case 'h': help(0); break;
            default: help(1);
        }
    }
    
    // Validate required arguments
    if (assignments_file.empty()) { fprintf(stderr, "ERROR: --assignments/-a is required\n"); exit(1); }
    if (diagnostics_file.empty()) { fprintf(stderr, "ERROR: --diagnostics/-d is required\n"); exit(1); }
    if (runner_ups_file.empty()) { fprintf(stderr, "ERROR: --runner_ups/-r is required\n"); exit(1); }
    if (expected_lines_file.empty()) { fprintf(stderr, "ERROR: --expected/-e is required\n"); exit(1); }
    if (output_prefix.empty()) { fprintf(stderr, "ERROR: --output/-o is required\n"); exit(1); }
    
    fprintf(stderr, "tetra_refine v3.2\n");
    fprintf(stderr, "===============\n");
    
    // ========================================================================
    // Load data
    // ========================================================================
    
    if (verbose) fprintf(stderr, "Loading expected lines from %s...\n", expected_lines_file.c_str());
    ExpectedLines expected;
    parse_expected_lines(expected_lines_file, expected);
    
    if (verbose) {
        fprintf(stderr, "  Diploid singlets: %lu\n", expected.diploid_singlets.size());
        fprintf(stderr, "  Heterotypic tetraploids: %lu\n", expected.tetraploid_heterotypic.size());
        fprintf(stderr, "  Homotypic tetraploids: %lu\n", expected.tetraploid_homotypic.size());
    }
    
    if (verbose) fprintf(stderr, "Loading assignments from %s...\n", assignments_file.c_str());
    map<string, Assignment> assignments;
    parse_assignments_file(assignments_file, assignments);
    if (verbose) fprintf(stderr, "  Loaded %lu assignments\n", assignments.size());
    
    if (verbose) fprintf(stderr, "Loading diagnostics from %s...\n", diagnostics_file.c_str());
    map<string, Diagnostics> diagnostics;
    bool has_posterior_data = false;
    parse_diagnostics_file(diagnostics_file, diagnostics, has_posterior_data);
    if (verbose) {
        fprintf(stderr, "  Loaded %lu diagnostics\n", diagnostics.size());
        fprintf(stderr, "  Posterior/entropy columns: %s\n", has_posterior_data ? "YES" : "NO");
    }
    
    if (verbose) fprintf(stderr, "Loading runner-ups from %s...\n", runner_ups_file.c_str());
    map<string, vector<RunnerUp>> runner_ups;
    parse_runner_ups_file(runner_ups_file, runner_ups);
    if (verbose) fprintf(stderr, "  Loaded runner-ups for %lu cells\n", runner_ups.size());
    
    // Optional: contamination rates
    map<string, double> contam_rates;
    if (!contam_rate_file.empty()) {
        if (verbose) fprintf(stderr, "Loading contamination rates from %s...\n", contam_rate_file.c_str());
        parse_contam_rate_file(contam_rate_file, contam_rates);
    }
    
    // Optional: external ploidy
    map<string, ExternalPloidy> ext_ploidy_map;
    if (!external_ploidy_file.empty()) {
        if (verbose) fprintf(stderr, "Loading external ploidy from %s...\n", external_ploidy_file.c_str());
        parse_external_ploidy_file(external_ploidy_file, ext_ploidy_map, external_ploidy_library);
    }
    
    // ========================================================================
    // Compute thresholds from data
    // ========================================================================
    
    vector<double> het_vars_singlet;
    vector<double> min_margins_d;   // D-type only for quad detection thresholds
    vector<double> depths;
    bool has_het_data = false;
    
    for (const auto& dp : diagnostics) {
        depths.push_back(dp.second.total_depth);
        
        if (dp.second.singlet_doublet == 'D') {
            min_margins_d.push_back(dp.second.min_margin);
        }
        
        if (dp.second.singlet_doublet == 'S' && dp.second.het_balance_var >= 0) {
            het_vars_singlet.push_back(dp.second.het_balance_var);
            has_het_data = true;
        }
    }
    
    double median_depth = !depths.empty() ? compute_median(depths) : 0.0;
    
    Thresholds thresholds;
    thresholds.source = "";
    
    // Min margin thresholds for quad detection (D-type only)
    if (min_margin_override > 0) {
        thresholds.min_margin_threshold = min_margin_override;
        thresholds.min_margin_possible = min_margin_override * 2.5;
        thresholds.source += "min_margin: user (strict), user*2.5 (lenient); ";
    } else if (!min_margins_d.empty()) {
        thresholds.min_margin_threshold = compute_percentile(min_margins_d, 10.0);
        thresholds.min_margin_possible = compute_percentile(min_margins_d, 25.0);
        thresholds.source += "min_margin: D-type p10 (strict), p25 (lenient); ";
    } else {
        thresholds.min_margin_threshold = 50.0;
        thresholds.min_margin_possible = 125.0;
        thresholds.source += "min_margin: default; ";
    }
    
    // Het variance thresholds (annotation only in v3)
    if (het_var_diploid_override > 0 && het_var_tetraploid_override > 0) {
        thresholds.het_var_diploid = het_var_diploid_override;
        thresholds.het_var_tetraploid = het_var_tetraploid_override;
        thresholds.source += "het_var: user (annotation)";
    } else if (has_het_data && !het_vars_singlet.empty()) {
        double p25 = compute_percentile(het_vars_singlet, 25.0);
        double p75 = compute_percentile(het_vars_singlet, 75.0);
        thresholds.het_var_tetraploid = p25;
        thresholds.het_var_diploid = p75;
        thresholds.source += "het_var: 25th/75th percentile (annotation)";
    } else {
        thresholds.het_var_diploid = 0.0;
        thresholds.het_var_tetraploid = 0.0;
        thresholds.source += "het_var: N/A";
    }
    
    thresholds.close_count_threshold = 1;
    thresholds.depth_ratio_threshold = 1.5;
    
    if (verbose) {
        fprintf(stderr, "\nThresholds:\n");
        fprintf(stderr, "  min_margin_threshold (strict): %.2f\n", thresholds.min_margin_threshold);
        fprintf(stderr, "  min_margin_possible (lenient): %.2f\n", thresholds.min_margin_possible);
        fprintf(stderr, "  het_var_tetraploid: %.6f (annotation only)\n", thresholds.het_var_tetraploid);
        fprintf(stderr, "  het_var_diploid: %.6f (annotation only)\n", thresholds.het_var_diploid);
        fprintf(stderr, "  median_depth: %.1f\n", median_depth);
        fprintf(stderr, "  Source: %s\n", thresholds.source.c_str());
    }
    
    // ========================================================================
    // Process each cell
    // ========================================================================
    
    if (verbose) fprintf(stderr, "\nProcessing cells...\n");
    
    vector<RefinedCell> refined_cells;
    SummaryStats stats = {};
    
    for (const auto& ap : assignments) {
        const string& barcode = ap.first;
        const Assignment& assn = ap.second;
        
        // Get diagnostics and runner-ups if available
        Diagnostics* diag_ptr = NULL;
        vector<RunnerUp>* runner_ptr = NULL;
        
        if (diagnostics.find(barcode) != diagnostics.end()) {
            diag_ptr = &diagnostics[barcode];
        }
        if (runner_ups.find(barcode) != runner_ups.end()) {
            runner_ptr = &runner_ups[barcode];
        }
        
        // Get external ploidy if available (only for singlets in Case 2c)
        ExternalPloidy* ext_ptr = NULL;
        if (!ext_ploidy_map.empty() && ext_ploidy_map.find(barcode) != ext_ploidy_map.end()) {
            ext_ptr = &ext_ploidy_map[barcode];
        }
        
        RefinedCell cell;
        classify_cell(assn, diag_ptr, runner_ptr, ext_ptr, expected, thresholds,
                      median_depth, has_het_data, cell);
        
        // Compute confidence metrics
        if (diag_ptr != NULL) {
            cell.dnllr = compute_dnllr(assn.llr, diag_ptr->total_depth);
            cell.margin_ratio = compute_margin_ratio(diag_ptr->min_margin, runner_ptr);
            cell.posterior = diag_ptr->posterior;
            cell.entropy = diag_ptr->entropy;
        } else {
            cell.dnllr = -1.0;
            cell.margin_ratio = 1e6;
            cell.posterior = -1.0;
            cell.entropy = -1.0;
        }
        
        // Carry through contamination rate
        if (!contam_rates.empty() && contam_rates.find(barcode) != contam_rates.end()) {
            cell.contam_rate = contam_rates[barcode];
        } else {
            cell.contam_rate = -1.0;
        }
        
        // Update stats
        stats.total_cells++;
        
        if (assn.type == 'S') stats.original_singlets++;
        else stats.original_doublets++;
        
        if (cell.ploidy == "diploid") stats.ploidy_diploid++;
        else if (cell.ploidy == "tetraploid") stats.ploidy_tetraploid++;
        else stats.ploidy_unknown++;
        
        if (cell.ploidy_method == "heterotypic") stats.tetraploid_heterotypic++;
        else if (cell.ploidy_method == "expected_lines") stats.tetraploid_from_expected++;
        else if (cell.ploidy_method == "external") stats.tetraploid_from_external++;
        
        if (cell.droplet_flag == "NONE") stats.droplet_singlet++;
        else if (cell.droplet_flag == "POSSIBLE_DOUBLET") stats.droplet_possible++;
        else stats.droplet_likely++;
        
        if (cell.original_assignment != cell.refined_assignment) stats.assignment_changed++;
        if (cell.changed) stats.cells_changed++;
        
        refined_cells.push_back(cell);
    }
    
    // ========================================================================
    // Write outputs
    // ========================================================================
    
    string refined_file = output_prefix + ".refined_assignments";
    if (verbose) fprintf(stderr, "Writing refined assignments to %s...\n", refined_file.c_str());
    write_refined_assignments(refined_file, refined_cells);
    
    if (write_changed_only) {
        string changed_file = output_prefix + ".refined_changed";
        if (verbose) fprintf(stderr, "Writing changed cells to %s...\n", changed_file.c_str());
        vector<RefinedCell> changed_cells;
        for (const auto& cell : refined_cells) {
            if (cell.changed) changed_cells.push_back(cell);
        }
        write_refined_assignments(changed_file, changed_cells);
    }
    
    if (write_simple) {
        string simple_file = output_prefix + ".assignments_refined";
        if (verbose) fprintf(stderr, "Writing simple assignments to %s...\n", simple_file.c_str());
        write_simple_assignments(simple_file, refined_cells);
    }
    
    string summary_file = output_prefix + ".refine_summary";
    if (verbose) fprintf(stderr, "Writing summary to %s...\n", summary_file.c_str());
    write_summary(summary_file, assignments_file, diagnostics_file, runner_ups_file,
                  expected_lines_file, contam_rate_file, external_ploidy_file,
                  external_ploidy_library,
                  expected, thresholds, has_het_data, has_posterior_data, median_depth, stats);
    
    // Print summary to stderr
    fprintf(stderr, "\ntetra_refine v3.2 complete\n");
    fprintf(stderr, "========================\n");
    fprintf(stderr, "Total cells: %d\n", stats.total_cells);
    fprintf(stderr, "\nPloidy classification:\n");
    fprintf(stderr, "  Diploid: %d\n", stats.ploidy_diploid);
    fprintf(stderr, "  Tetraploid: %d\n", stats.ploidy_tetraploid);
    fprintf(stderr, "    - Heterotypic (A+B): %d\n", stats.tetraploid_heterotypic);
    fprintf(stderr, "    - Homotypic from expected_lines: %d\n", stats.tetraploid_from_expected);
    fprintf(stderr, "    - From external ploidy: %d\n", stats.tetraploid_from_external);
    fprintf(stderr, "  Unknown: %d\n", stats.ploidy_unknown);
    fprintf(stderr, "\nDroplet status:\n");
    fprintf(stderr, "  Singlets (1 cell): %d\n", stats.droplet_singlet);
    fprintf(stderr, "  Possible doublets: %d\n", stats.droplet_possible);
    fprintf(stderr, "  Likely doublets: %d\n", stats.droplet_likely);
    fprintf(stderr, "\nAssignments changed: %d (%.1f%%)\n",
            stats.assignment_changed,
            stats.total_cells > 0 ? 100.0 * stats.assignment_changed / stats.total_cells : 0.0);
    fprintf(stderr, "\nData sources:\n");
    fprintf(stderr, "  Posterior/entropy: %s\n", has_posterior_data ? "available" : "not in diagnostics");
    fprintf(stderr, "  Contamination rate: %s\n", !contam_rate_file.empty() ? "loaded" : "not provided");
    fprintf(stderr, "  External ploidy: %s\n", !external_ploidy_file.empty() ? "loaded" : "not provided");
    fprintf(stderr, "\nOutput files:\n");
    fprintf(stderr, "  %s\n", refined_file.c_str());
    if (write_changed_only) fprintf(stderr, "  %s.refined_changed\n", output_prefix.c_str());
    if (write_simple) fprintf(stderr, "  %s.assignments_refined\n", output_prefix.c_str());
    fprintf(stderr, "  %s\n", summary_file.c_str());
    
    return 0;
}
