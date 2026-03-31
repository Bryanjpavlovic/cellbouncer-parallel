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
// TERMINOLOGY - IMPORTANT!
// ============================================================================
// 
// PLOIDY (property of a CELL):
//   - diploid: cell has 2 copies of genome (normal cell)
//   - tetraploid: cell has 4 copies of genome (fused cell)
//     - heterotypic tetraploid: fused from two different individuals (A×B)
//     - homotypic tetraploid: fused from same individual (A×A)
//
// DROPLET CONTENTS:
//   - 1 cell in droplet: singlet (can be diploid OR tetraploid)
//   - 2 cells in droplet: doublet (both could be diploid, both tetraploid, or mixed)
//
// demux_parallel OUTPUT:
//   - "S" (singlet): assigned to one individual identity
//   - "D" (doublet): assigned to two individual identities (A+B)
//
// THE TRICK: We represent heterotypic tetraploids as "D" in demux because they
// have alleles from two individuals. But they're actually SINGLETS (one cell).
//
// This tool's job:
//   1. Recognize that "D" assignments are tetraploid SINGLETS, not diploid doublets
//   2. Detect homotypic tetraploids among "S" assignments using het_balance_var
//   3. Detect actual doublets (2 cells in droplet) - the "quad" problem
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
};

struct RunnerUp {
    string barcode;
    int rank;
    string identity;
    double llr_vs_winner;
    double min_margin;
};

// NEW: Clear output structure with proper terminology
struct RefinedCell {
    string barcode;
    
    // Original demux output
    string original_assignment;  // e.g., "H20961" or "H20961+H27322"
    char original_type;          // 'S' or 'D' from demux
    double llr;
    
    // Refined assignment
    string refined_assignment;   // e.g., "H20961+H20961" if homotypic tetraploid detected
    
    // PLOIDY: diploid vs tetraploid (property of the cell)
    string ploidy;              // "diploid", "tetraploid", "unknown"
    string ploidy_method;       // How we determined ploidy: "heterotypic", "expected_lines", "het_var", "depth", "ambiguous"
    string ploidy_confidence;   // "HIGH", "MEDIUM", "LOW"
    
    // DROPLET: how many cells in droplet
    int cells_in_droplet;       // 1 = singlet, 2 = doublet (quad = doublet of tetraploids)
    string droplet_flag;        // "NONE", "POSSIBLE_DOUBLET", "LIKELY_DOUBLET"
    string droplet_candidates;  // If flagged: the suspected identities
    
    // Overall
    string overall_confidence;
    bool changed;               // TRUE if anything differs from original interpretation
};

struct Thresholds {
    double min_margin_threshold;
    double close_count_threshold;
    double het_var_diploid;
    double het_var_tetraploid;
    double depth_ratio_threshold;  // Depth ratio for tetraploid detection
    string source;
};

struct SummaryStats {
    int total_cells;
    
    // Original assignment types
    int original_singlets;
    int original_doublets;
    
    // Ploidy classification
    int ploidy_diploid;
    int ploidy_tetraploid;
    int ploidy_unknown;
    
    // Tetraploid subtypes
    int tetraploid_heterotypic;
    int tetraploid_homotypic;
    int tetraploid_from_expected;
    int tetraploid_from_het_var;
    
    // Droplet classification
    int droplet_singlet;
    int droplet_possible;
    int droplet_likely;
    
    // Changes
    int assignment_changed;
    int ploidy_inferred;
    int cells_changed;
};

// ============================================================================
// Expected Lines Parsing
// ============================================================================

struct ExpectedLines {
    set<string> diploid_singlets;      // Individual names that appear alone (diploid)
    set<string> tetraploid_heterotypic; // A+B combinations (heterotypic tetraploid)
    set<string> tetraploid_homotypic;   // A+A combinations (homotypic tetraploid)
    set<string> all_individuals;        // All individual names seen
    
    // For quick lookup
    map<string, bool> singlet_has_diploid;    // Does individual X have diploid form in pool?
    map<string, bool> singlet_has_homotypic;  // Does individual X have homotypic tetraploid form?
};

void parse_expected_lines(const string& filename, ExpectedLines& expected) {
    ifstream inf(filename);
    if (!inf.is_open()) {
        fprintf(stderr, "ERROR: Cannot open expected lines file: %s\n", filename.c_str());
        exit(1);
    }
    
    string line;
    while (getline(inf, line)) {
        if (line.empty()) continue;
        
        // Check if it's a doublet (A+B format)
        size_t plus_pos = line.find('+');
        if (plus_pos != string::npos) {
            string id1 = line.substr(0, plus_pos);
            string id2 = line.substr(plus_pos + 1);
            
            expected.all_individuals.insert(id1);
            expected.all_individuals.insert(id2);
            
            if (id1 == id2) {
                // Homotypic tetraploid (A+A)
                expected.tetraploid_homotypic.insert(line);
                expected.singlet_has_homotypic[id1] = true;
            } else {
                // Heterotypic tetraploid (A+B)
                expected.tetraploid_heterotypic.insert(line);
            }
        } else {
            // Diploid singlet
            expected.diploid_singlets.insert(line);
            expected.all_individuals.insert(line);
            expected.singlet_has_diploid[line] = true;
        }
    }
    
    // Fill in missing entries as false
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

void parse_diagnostics_file(const string& filename, map<string, Diagnostics>& diagnostics) {
    gzreader reader(filename);
    bool header_read = false;
    
    while (reader.next()) {
        if (!header_read) {
            header_read = true;
            continue;
        }
        
        istringstream splitter(reader.line);
        string field;
        Diagnostics d;
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
// Droplet Detection (Doublet of Tetraploids - the "Quad" Problem)
// ============================================================================
//
// A true quad (two tetraploid cells A×B and C×D in one droplet) should show:
// 1. Winner is some doublet (A+B)
// 2. There exists a complementary pair C+D in runner-ups where {C,D} ∩ {A,B} = ∅
// 3. Cross-combinations (A+C, A+D, B+C, B+D) appear in runner-ups
// 4. Low min_margin (winner barely beats alternatives)
// 5. Multiple close competitors (n_close >= 1)
//
// Scoring: count how many of the 6 expected combinations are present
// (winner A+B, complement C+D, plus 4 cross-combos)

// Helper: check if two identity sets are disjoint
bool sets_are_disjoint(const set<string>& s1, const set<string>& s2) {
    for (const auto& elem : s1) {
        if (s2.count(elem) > 0) return false;
    }
    return true;
}

// Helper: create normalized doublet string (alphabetically sorted)
string make_doublet_key(const string& id1, const string& id2) {
    if (id1 < id2) return id1 + "+" + id2;
    return id2 + "+" + id1;
}

string detect_droplet_doublet(const Diagnostics& diag,
                              const vector<RunnerUp>& runner_ups,
                              const Thresholds& thresholds,
                              string& candidates) {
    
    candidates = ".";
    
    // Only check cells assigned as heterotypic tetraploids (A+B where A != B)
    if (diag.singlet_doublet != 'D') {
        return "NONE";
    }
    
    // Parse winner identity
    string winner_id1, winner_id2;
    if (!parse_doublet_identity(diag.assignment, winner_id1, winner_id2)) {
        return "NONE";
    }
    
    // Skip homotypic (shouldn't happen but check anyway)
    if (winner_id1 == winner_id2) {
        return "NONE";
    }
    
    set<string> winner_ids;
    winner_ids.insert(winner_id1);
    winner_ids.insert(winner_id2);
    
    // Build set of all runner-up identities (normalized keys)
    set<string> runner_up_keys;
    map<string, set<string>> runner_up_components;  // key -> {id1, id2}
    
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
    
    // Look for a complementary pair: C+D where {C,D} is disjoint from {A,B}
    string best_complement_key = "";
    set<string> complement_ids;
    
    for (const auto& ru_pair : runner_up_components) {
        if (sets_are_disjoint(winner_ids, ru_pair.second)) {
            // Found a disjoint pair - this is a potential complement
            best_complement_key = ru_pair.first;
            complement_ids = ru_pair.second;
            break;  // Take the first (highest ranked) one
        }
    }
    
    if (best_complement_key.empty()) {
        // No complementary pair found - not a quad
        return "NONE";
    }
    
    // Now check for cross-combinations
    // For quad A+B+C+D, expect: A+C, A+D, B+C, B+D (4 cross-combos)
    int cross_combos_found = 0;
    int cross_combos_expected = 4;
    
    for (const auto& w_id : winner_ids) {
        for (const auto& c_id : complement_ids) {
            string cross_key = make_doublet_key(w_id, c_id);
            if (runner_up_keys.count(cross_key) > 0) {
                cross_combos_found++;
            }
        }
    }
    
    // Calculate quad pattern score
    // Total expected: complement (1) + cross-combos (4) = 5
    // (winner is implicit, not counted)
    double pattern_score = (1.0 + cross_combos_found) / 5.0;
    
    // Build the quad candidates string (all 4 individuals)
    vector<string> all_four;
    for (const auto& id : winner_ids) all_four.push_back(id);
    for (const auto& id : complement_ids) all_four.push_back(id);
    sort(all_four.begin(), all_four.end());
    
    candidates = "";
    for (size_t i = 0; i < all_four.size(); i++) {
        if (i > 0) candidates += "+";
        candidates += all_four[i];
    }
    
    // Decision logic based on multiple signals
    bool low_margin = (diag.min_margin < thresholds.min_margin_threshold);
    bool has_close = (diag.n_close >= 1);
    bool strong_pattern = (pattern_score >= 0.6);  // At least complement + 2 cross-combos
    bool weak_pattern = (pattern_score >= 0.4);    // At least complement + 1 cross-combo
    
    // HIGH confidence: low margin AND close competitors AND strong pattern
    if (low_margin && has_close && strong_pattern) {
        return "LIKELY_DOUBLET";
    }
    
    // MEDIUM confidence: need at least 2 of 3 signals, plus weak pattern
    int signals = (low_margin ? 1 : 0) + (has_close ? 1 : 0) + (strong_pattern ? 1 : 0);
    if (signals >= 2 && weak_pattern) {
        return "POSSIBLE_DOUBLET";
    }
    
    // Not enough evidence
    candidates = ".";
    return "NONE";
}

// ============================================================================
// Ploidy Classification
// ============================================================================

void classify_cell(const Assignment& assn,
                   const Diagnostics* diag,
                   const vector<RunnerUp>* runner_ups,
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
    cell.cells_in_droplet = 1;  // Default: assume 1 cell
    cell.droplet_flag = "NONE";
    cell.droplet_candidates = ".";
    
    // =========================================================================
    // CASE 1: Heterotypic tetraploid (demux assigned A+B where A != B)
    // =========================================================================
    if (assn.type == 'D') {
        string id1, id2;
        parse_doublet_identity(assn.identity, id1, id2);
        
        if (id1 != id2) {
            // Heterotypic tetraploid - this is DEFINITELY a tetraploid
            cell.ploidy = "tetraploid";
            cell.ploidy_method = "heterotypic";
            cell.ploidy_confidence = "HIGH";
            cell.cells_in_droplet = 1;  // One tetraploid cell, not two diploid cells!
            
            // Check for droplet doublet (quad)
            if (diag != NULL && runner_ups != NULL) {
                cell.droplet_flag = detect_droplet_doublet(*diag, *runner_ups, 
                                                           thresholds, cell.droplet_candidates);
                if (cell.droplet_flag != "NONE") {
                    cell.cells_in_droplet = 2;  // Two tetraploid cells
                }
            }
            
            cell.overall_confidence = (cell.droplet_flag == "NONE") ? "HIGH" : "LOW";
            return;
        } else {
            // This shouldn't happen (demux shouldn't output A+A as "D")
            // But if it does, treat as homotypic tetraploid
            cell.ploidy = "tetraploid";
            cell.ploidy_method = "homotypic_direct";
            cell.ploidy_confidence = "HIGH";
            cell.cells_in_droplet = 1;
            cell.overall_confidence = "HIGH";
            return;
        }
    }
    
    // =========================================================================
    // CASE 2: Singlet assignment - need to determine if diploid or homotypic tetraploid
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
    // CASE 2c: Both or neither in pool - use het_balance_var if available
    // -------------------------------------------------------------------------
    if (diag != NULL && diag->het_balance_var >= 0 && has_het_data) {
        // Use het_balance_var to distinguish
        // Lower variance = tetraploid (sampling from 4 alleles)
        // Higher variance = diploid (sampling from 2 alleles)
        
        double threshold = (thresholds.het_var_diploid + thresholds.het_var_tetraploid) / 2.0;
        
        if (diag->het_balance_var < threshold) {
            // Low variance suggests tetraploid
            if (has_homotypic_in_pool) {
                cell.refined_assignment = make_homotypic(indv);
                cell.changed = true;
            }
            cell.ploidy = "tetraploid";
            cell.ploidy_method = "het_var";
            
            // Confidence based on how far from threshold
            double dist_from_thresh = threshold - diag->het_balance_var;
            if (dist_from_thresh > threshold * 0.3) {
                cell.ploidy_confidence = "HIGH";
            } else {
                cell.ploidy_confidence = "MEDIUM";
            }
        } else {
            // High variance suggests diploid
            cell.ploidy = "diploid";
            cell.ploidy_method = "het_var";
            
            double dist_from_thresh = diag->het_balance_var - threshold;
            if (dist_from_thresh > threshold * 0.3) {
                cell.ploidy_confidence = "HIGH";
            } else {
                cell.ploidy_confidence = "MEDIUM";
            }
        }
        
        cell.cells_in_droplet = 1;
        cell.overall_confidence = cell.ploidy_confidence;
        return;
    }
    
    // -------------------------------------------------------------------------
    // CASE 2d: No het_var data - use depth as weak signal
    // -------------------------------------------------------------------------
    if (diag != NULL && median_depth > 0) {
        bool depth_suggests_tetraploid = (diag->total_depth > median_depth * 1.5);
        
        if (depth_suggests_tetraploid && has_homotypic_in_pool) {
            cell.refined_assignment = make_homotypic(indv);
            cell.ploidy = "tetraploid";
            cell.ploidy_method = "depth";
            cell.ploidy_confidence = "LOW";
            cell.changed = true;
        } else if (!depth_suggests_tetraploid && has_diploid_in_pool) {
            cell.ploidy = "diploid";
            cell.ploidy_method = "depth";
            cell.ploidy_confidence = "LOW";
        } else {
            cell.ploidy = "unknown";
            cell.ploidy_method = "ambiguous";
            cell.ploidy_confidence = "LOW";
        }
        
        cell.cells_in_droplet = 1;
        cell.overall_confidence = "LOW";
        return;
    }
    
    // -------------------------------------------------------------------------
    // CASE 2e: No diagnostics at all - unknown
    // -------------------------------------------------------------------------
    cell.ploidy = "unknown";
    cell.ploidy_method = "no_data";
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
    
    // Header
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
         << "changed" << endl;
    
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
             << (cell.changed ? "TRUE" : "FALSE") << endl;
    }
}

// Write a simple assignments file format (barcode, identity, ploidy, cells_in_droplet, llr)
// Clear columns, no ambiguous single-letter codes
void write_simple_assignments(const string& filename, const vector<RefinedCell>& cells) {
    ofstream outf(filename);
    if (!outf.is_open()) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", filename.c_str());
        exit(1);
    }
    
    // Header
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
                   const ExpectedLines& expected,
                   const Thresholds& thresholds,
                   bool has_het_data,
                   double median_depth,
                   const SummaryStats& stats) {
    
    ofstream outf(filename);
    if (!outf.is_open()) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", filename.c_str());
        exit(1);
    }
    
    outf << "tetra_refine v2 summary" << endl;
    outf << "=======================" << endl;
    outf << endl;
    
    outf << "Input files:" << endl;
    outf << "  Assignments: " << assignments_file << endl;
    outf << "  Diagnostics: " << diagnostics_file << endl;
    outf << "  Runner-ups: " << runner_ups_file << endl;
    outf << "  Expected lines: " << expected_lines_file << endl;
    outf << endl;
    
    outf << "Expected lines summary:" << endl;
    outf << "  Diploid singlets: " << expected.diploid_singlets.size() << endl;
    outf << "  Heterotypic tetraploids: " << expected.tetraploid_heterotypic.size() << endl;
    outf << "  Homotypic tetraploids: " << expected.tetraploid_homotypic.size() << endl;
    outf << "  Total individuals: " << expected.all_individuals.size() << endl;
    outf << endl;
    
    outf << "Data availability:" << endl;
    outf << "  Het variance data: " << (has_het_data ? "YES" : "NO") << endl;
    outf << "  Median depth: " << median_depth << endl;
    outf << endl;
    
    outf << "Thresholds:" << endl;
    outf << "  min_margin_threshold: " << thresholds.min_margin_threshold << endl;
    outf << "  het_var_diploid: " << thresholds.het_var_diploid << endl;
    outf << "  het_var_tetraploid: " << thresholds.het_var_tetraploid << endl;
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
    outf << "    - Homotypic (from het_var): " << stats.tetraploid_from_het_var << endl;
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
    fprintf(stderr, "tetra_refine v2 [OPTIONS]\n");
    fprintf(stderr, "Refines demux_parallel assignments for tetraploid pools.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "This tool:\n");
    fprintf(stderr, "  1. Recognizes heterotypic tetraploids (A+B) as single cells, not doublets\n");
    fprintf(stderr, "  2. Detects homotypic tetraploids (A+A) among singlet assignments\n");
    fprintf(stderr, "  3. Flags potential doublets of tetraploids (quads)\n");
    fprintf(stderr, "\n===== REQUIRED =====\n");
    fprintf(stderr, "    --assignments -a FILE    .assignments file from demux_parallel\n");
    fprintf(stderr, "    --diagnostics -d FILE    .diagnostics.gz file from demux_parallel\n");
    fprintf(stderr, "    --runner_ups -r FILE     .runner_ups.gz file from demux_parallel\n");
    fprintf(stderr, "    --expected -e FILE       Expected lines file (same format as demux_parallel -I)\n");
    fprintf(stderr, "    --output -o PREFIX       Output file prefix\n");
    fprintf(stderr, "\n===== OPTIONAL - Thresholds =====\n");
    fprintf(stderr, "    --het_var_diploid F      Het variance threshold for diploid (default: auto)\n");
    fprintf(stderr, "    --het_var_tetraploid F   Het variance threshold for tetraploid (default: auto)\n");
    fprintf(stderr, "    --min_margin F           Min margin threshold for doublet detection (default: 10th percentile)\n");
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
            case 'a':
                assignments_file = optarg;
                break;
            case 'd':
                diagnostics_file = optarg;
                break;
            case 'r':
                runner_ups_file = optarg;
                break;
            case 'e':
                expected_lines_file = optarg;
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 1001:
                het_var_diploid_override = atof(optarg);
                break;
            case 1002:
                het_var_tetraploid_override = atof(optarg);
                break;
            case 1003:
                min_margin_override = atof(optarg);
                break;
            case 1004:
                write_changed_only = true;
                break;
            case 1005:
                write_simple = true;
                break;
            case 'V':
                verbose = true;
                break;
            case 'h':
                help(0);
                break;
            default:
                help(1);
        }
    }
    
    // Validate required arguments
    if (assignments_file.empty()) {
        fprintf(stderr, "ERROR: --assignments/-a is required\n");
        exit(1);
    }
    if (diagnostics_file.empty()) {
        fprintf(stderr, "ERROR: --diagnostics/-d is required\n");
        exit(1);
    }
    if (runner_ups_file.empty()) {
        fprintf(stderr, "ERROR: --runner_ups/-r is required\n");
        exit(1);
    }
    if (expected_lines_file.empty()) {
        fprintf(stderr, "ERROR: --expected/-e is required\n");
        exit(1);
    }
    if (output_prefix.empty()) {
        fprintf(stderr, "ERROR: --output/-o is required\n");
        exit(1);
    }
    
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
    parse_diagnostics_file(diagnostics_file, diagnostics);
    if (verbose) fprintf(stderr, "  Loaded %lu diagnostics\n", diagnostics.size());
    
    if (verbose) fprintf(stderr, "Loading runner-ups from %s...\n", runner_ups_file.c_str());
    map<string, vector<RunnerUp>> runner_ups;
    parse_runner_ups_file(runner_ups_file, runner_ups);
    if (verbose) fprintf(stderr, "  Loaded runner-ups for %lu cells\n", runner_ups.size());
    
    // ========================================================================
    // Compute thresholds from data
    // ========================================================================
    
    vector<double> het_vars_singlet;
    vector<double> min_margins;
    vector<double> depths;
    bool has_het_data = false;
    
    for (const auto& dp : diagnostics) {
        min_margins.push_back(dp.second.min_margin);
        depths.push_back(dp.second.total_depth);
        
        // Only collect het_var from singlets (doublets have -1)
        if (dp.second.singlet_doublet == 'S' && dp.second.het_balance_var >= 0) {
            het_vars_singlet.push_back(dp.second.het_balance_var);
            has_het_data = true;
        }
    }
    
    double median_depth = !depths.empty() ? compute_median(depths) : 0.0;
    
    Thresholds thresholds;
    thresholds.source = "";
    
    // Min margin threshold
    if (min_margin_override > 0) {
        thresholds.min_margin_threshold = min_margin_override;
        thresholds.source += "min_margin: user; ";
    } else if (!min_margins.empty()) {
        thresholds.min_margin_threshold = compute_percentile(min_margins, 10.0);
        thresholds.source += "min_margin: 10th percentile; ";
    } else {
        thresholds.min_margin_threshold = 50.0;  // Default
        thresholds.source += "min_margin: default; ";
    }
    
    // Het variance thresholds
    if (het_var_diploid_override > 0 && het_var_tetraploid_override > 0) {
        thresholds.het_var_diploid = het_var_diploid_override;
        thresholds.het_var_tetraploid = het_var_tetraploid_override;
        thresholds.source += "het_var: user";
    } else if (has_het_data && !het_vars_singlet.empty()) {
        // Use percentiles from singlet het_var distribution
        // Tetraploids should have LOWER variance than diploids
        double p25 = compute_percentile(het_vars_singlet, 25.0);
        double p50 = compute_percentile(het_vars_singlet, 50.0);
        double p75 = compute_percentile(het_vars_singlet, 75.0);
        
        // Heuristic: below 25th percentile → tetraploid, above 75th → diploid
        thresholds.het_var_tetraploid = p25;
        thresholds.het_var_diploid = p75;
        thresholds.source += "het_var: 25th/75th percentile";
    } else {
        thresholds.het_var_diploid = 0.0;
        thresholds.het_var_tetraploid = 0.0;
        thresholds.source += "het_var: N/A";
    }
    
    thresholds.close_count_threshold = 1;
    thresholds.depth_ratio_threshold = 1.5;
    
    if (verbose) {
        fprintf(stderr, "\nThresholds:\n");
        fprintf(stderr, "  min_margin_threshold: %.2f\n", thresholds.min_margin_threshold);
        fprintf(stderr, "  het_var_tetraploid: %.6f (below this → tetraploid)\n", thresholds.het_var_tetraploid);
        fprintf(stderr, "  het_var_diploid: %.6f (above this → diploid)\n", thresholds.het_var_diploid);
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
        
        RefinedCell cell;
        classify_cell(assn, diag_ptr, runner_ptr, expected, thresholds, 
                      median_depth, has_het_data, cell);
        
        // Update stats
        stats.total_cells++;
        
        if (assn.type == 'S') stats.original_singlets++;
        else stats.original_doublets++;
        
        if (cell.ploidy == "diploid") stats.ploidy_diploid++;
        else if (cell.ploidy == "tetraploid") stats.ploidy_tetraploid++;
        else stats.ploidy_unknown++;
        
        if (cell.ploidy_method == "heterotypic") stats.tetraploid_heterotypic++;
        else if (cell.ploidy_method == "expected_lines") stats.tetraploid_from_expected++;
        else if (cell.ploidy_method == "het_var") stats.tetraploid_from_het_var++;
        
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
            if (cell.changed) {
                changed_cells.push_back(cell);
            }
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
                  expected_lines_file, expected, thresholds, has_het_data, median_depth, stats);
    
    // Print summary to stderr
    fprintf(stderr, "\n");
    fprintf(stderr, "tetra_refine v2 complete\n");
    fprintf(stderr, "========================\n");
    fprintf(stderr, "Total cells: %d\n", stats.total_cells);
    fprintf(stderr, "\nPloidy classification:\n");
    fprintf(stderr, "  Diploid: %d\n", stats.ploidy_diploid);
    fprintf(stderr, "  Tetraploid: %d\n", stats.ploidy_tetraploid);
    fprintf(stderr, "    - Heterotypic (A+B): %d\n", stats.tetraploid_heterotypic);
    fprintf(stderr, "    - Homotypic from expected_lines: %d\n", stats.tetraploid_from_expected);
    fprintf(stderr, "    - Homotypic from het_var: %d\n", stats.tetraploid_from_het_var);
    fprintf(stderr, "  Unknown: %d\n", stats.ploidy_unknown);
    fprintf(stderr, "\nDroplet status:\n");
    fprintf(stderr, "  Singlets (1 cell): %d\n", stats.droplet_singlet);
    fprintf(stderr, "  Possible doublets: %d\n", stats.droplet_possible);
    fprintf(stderr, "  Likely doublets: %d\n", stats.droplet_likely);
    fprintf(stderr, "\nAssignments changed: %d (%.1f%%)\n", 
            stats.assignment_changed,
            stats.total_cells > 0 ? 100.0 * stats.assignment_changed / stats.total_cells : 0.0);
    fprintf(stderr, "\nOutput files:\n");
    fprintf(stderr, "  %s\n", refined_file.c_str());
    if (write_changed_only) {
        fprintf(stderr, "  %s.refined_changed\n", output_prefix.c_str());
    }
    if (write_simple) {
        fprintf(stderr, "  %s.assignments_refined\n", output_prefix.c_str());
    }
    fprintf(stderr, "  %s\n", summary_file.c_str());
    
    return 0;
}
