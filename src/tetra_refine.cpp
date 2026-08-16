// tetra_refine v3.5
//
// Revision history:
//   V2_R1 - Complete rewrite separating ploidy from droplet contents
//   V3_R1 - Add depth-normalized runner-up contrast, margin_ratio,
//           margin-softmax/entropy carry-through,
//           demote het_var to annotation, add --contam_rate and --external_ploidy inputs,
//           backward-compatible legacy and versioned name-addressed diagnostics parsing
//   V3_R3 - Add --scoring_only pass-through mode: emit scoring metrics only,
//           no ploidy reclassification, no droplet flagging, no threshold computation.
//           --external_ploidy is ignored in this mode (with warning).
//   V3_R4 - Version bump to distinguish this name-addressed build from the
//           position-parsed v3.3 build that shared the same version string.
//           This build reads .diagnostics.gz and .runner_ups.gz by COLUMN NAME
//           with alias tolerance: it accepts both the currently deployed demux
//           names (llr, posterior, entropy) and the next-version demux names
//           (llr_vs_runner_up, margin_softmax_score, margin_entropy), so it parses
//           either schema without column-position assumptions. Unavailable fields
//           are tokenized as NA rather than read as 0.0. Consumes the new
//           runnerup/worst comparison_state fields and gates DNLLR and margin_ratio
//           on comparisons that are actually present, so runner-ups with no
//           shared-site support no longer distort the derived scores.
//           Required to consume the next demux diagnostics schema
//           (schema_version = demux_parallel_diagnostics_v3); do NOT pair the
//           position-parsed v3.3 build with that schema, because its atof() parser
//           reads the NA fields as 0.0. Ploidy classification, the --external_ploidy
//           relabel gate (default 0.90), and --scoring_only are unchanged from v3.3.
//   V3_R5 - Homotypic entries in expected_lines are eligibility/candidate information
//           only. A demux singlet A is never converted to A+A from line expectation
//           alone. Relabeling A -> A+A now requires cell-specific external ploidy
//           evidence calling tetraploid at --external_ploidy_min_prob (default 0.90),
//           while A+A must still be expected for that identity. het_balance_var remains
//           carried through as orthogonal annotation for evaluation rather than an
//           unvalidated hard gate. Also retains the finite-subnormal diagnostics parser
//           fix: finite strtod() results are accepted even when underflow sets ERANGE.

#include <getopt.h>
#include <string>
#include <algorithm>
#include <cctype>
#include <cerrno>
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
#include <initializer_list>
#include <limits>

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
//   4. Compute cross-library diagnostics from the direct winner/runner-up contrast
//   5. Carry through margin_softmax_score/margin_entropy from demux_parallel and contam_rate
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
    char singlet_doublet = 'U';
    double llr_vs_runner_up = NAN;
    string runnerup_comparison_state = "unavailable";
    double min_margin = NAN;
    string worst_competitor;
    int n_close = 0;
    double total_depth = NAN;
    double het_balance_var = NAN;
    int n_het_sites = 0;
    double het_total_depth = NAN;
    double margin_softmax_score = NAN;
    double margin_entropy = NAN;
};

struct RunnerUp {
    string barcode;
    int rank = 0;
    string identity;
    double llr_vs_winner = NAN;
    string comparison_state = "unavailable";
    double min_margin = NAN;
};

struct RefinedCell {
    string barcode;
    
    // Original demux output
    string original_assignment;
    char original_type;
    double assignment_score;
    
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
    
    // Confidence/diagnostic metrics (v3)
    double llr_vs_runner_up;    // direct winner-versus-runner comparison
    string runnerup_comparison_state;
    double depth_normalized_llr_vs_runner_up;
    double margin_ratio;        // winner_min_margin / runner_up_rank1_min_margin
    double margin_softmax_score; // descriptive score; not a calibrated posterior
    double margin_entropy;
    
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
    int passthrough_count;   // v3.3: scoring-only mode count
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

static string trim_field(const string& raw) {
    const size_t first = raw.find_first_not_of(" \t\r\n");
    if (first == string::npos) return "";
    const size_t last = raw.find_last_not_of(" \t\r\n");
    return raw.substr(first, last - first + 1);
}

static string lowercase_field(string value) {
    transform(value.begin(), value.end(), value.begin(),
        [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    return value;
}

static vector<string> split_tsv_strict(const string& line) {
    vector<string> fields;
    size_t start = 0;
    while (true) {
        const size_t pos = line.find('\t', start);
        if (pos == string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return fields;
}

static void schema_error(const string& message) {
    fprintf(stderr, "ERROR: %s\n", message.c_str());
    exit(1);
}

static bool is_missing_diagnostic_token(const string& raw) {
    const string token = lowercase_field(trim_field(raw));
    return token.empty() || token == "." || token == "na" ||
           token == "unavailable" || token == "partial_support" ||
           token == "not_applicable";
}

static double parse_optional_diagnostic_number(
    const string& raw, const string& context) {
    if (is_missing_diagnostic_token(raw)) return NAN;
    errno = 0;
    char* end = NULL;
    const double value = strtod(raw.c_str(), &end);
    if (end == raw.c_str() || *end != '\0' || !isfinite(value)) {
        schema_error(context + ": expected a finite number or an explicit missing state, saw '" + raw + "'");
    }
    return value;
}

static long parse_required_integer(const string& raw, const string& context) {
    errno = 0;
    char* end = NULL;
    const long value = strtol(raw.c_str(), &end, 10);
    if (errno != 0 || end == raw.c_str() || *end != '\0') {
        schema_error(context + ": expected an integer, saw '" + raw + "'");
    }
    return value;
}

static int parse_optional_integer(const string& raw, const string& context) {
    if (is_missing_diagnostic_token(raw)) return 0;
    const long value = parse_required_integer(raw, context);
    if (value < numeric_limits<int>::min() || value > numeric_limits<int>::max()) {
        schema_error(context + ": integer is out of range");
    }
    return static_cast<int>(value);
}

static map<string, int> make_header_index(
    const vector<string>& header, const string& filename) {
    map<string, int> index;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) {
        if (header[i].empty()) schema_error(filename + ": empty column name");
        if (index.count(header[i]) != 0) {
            schema_error(filename + ": duplicate column name '" + header[i] + "'");
        }
        index[header[i]] = i;
    }
    return index;
}

static int optional_column(
    const map<string, int>& index,
    initializer_list<const char*> names) {
    for (const char* name : names) {
        const auto it = index.find(name);
        if (it != index.end()) return it->second;
    }
    return -1;
}

static int required_column(
    const map<string, int>& index,
    initializer_list<const char*> names,
    const string& filename) {
    const int column = optional_column(index, names);
    if (column >= 0) return column;
    string wanted;
    for (const char* name : names) {
        if (!wanted.empty()) wanted += " or ";
        wanted += name;
    }
    schema_error(filename + ": required column missing: " + wanted);
    return -1;
}

static bool comparison_state_is_present(const string& state) {
    return state == "present_nonzero" || state == "present_zero";
}

static string parse_comparison_state(
    const string& raw_state,
    double numeric_value,
    bool explicit_state,
    const string& context) {
    if (!explicit_state) {
        if (!isfinite(numeric_value)) return "unavailable";
        return numeric_value == 0.0 ? "present_zero" : "present_nonzero";
    }
    const string state = lowercase_field(trim_field(raw_state));
    static const set<string> supported = {
        "present_nonzero", "present_zero", "unavailable",
        "partial_support", "not_applicable"
    };
    if (supported.count(state) == 0) {
        schema_error(context + ": unsupported comparison state '" + raw_state + "'");
    }
    if (comparison_state_is_present(state)) {
        if (!isfinite(numeric_value)) {
            schema_error(context + ": state " + state + " requires a numeric direct comparison");
        }
        if (state == "present_zero" && numeric_value != 0.0) {
            schema_error(context + ": present_zero requires numeric zero");
        }
        if (state == "present_nonzero" && numeric_value == 0.0) {
            schema_error(context + ": present_nonzero cannot carry numeric zero");
        }
    } else if (isfinite(numeric_value)) {
        schema_error(context + ": missing comparison state " + state +
                     " must not carry a numeric value");
    }
    return state;
}

static bool supported_schema(const string& schema, const string& prefix) {
    return schema == prefix + "_v2" || schema == prefix + "_v3";
}

// Parse both legacy name-addressed diagnostics and the versioned corrected
// schemas. Column order is never used as a contract.
void parse_diagnostics_file(const string& filename,
                            map<string, Diagnostics>& diagnostics,
                            bool& has_margin_softmax_data) {
    gzreader reader(filename);
    if (!reader.next()) schema_error(filename + ": empty diagnostics file");

    const vector<string> header = split_tsv_strict(reader.line);
    const map<string, int> index = make_header_index(header, filename);
    const int barcode_col = required_column(index, {"barcode"}, filename);
    const int assignment_col = required_column(index, {"assignment"}, filename);
    const int sd_col = required_column(index, {"singlet_doublet"}, filename);
    const int direct_col = required_column(index, {"llr_vs_runner_up", "llr"}, filename);
    const int min_margin_col = required_column(index, {"min_margin"}, filename);
    const int worst_col = required_column(index, {"worst_competitor"}, filename);
    const int n_close_col = required_column(index, {"n_close"}, filename);
    const int total_depth_col = required_column(index, {"total_depth"}, filename);
    const int het_var_col = optional_column(index, {"het_balance_var"});
    const int n_het_col = optional_column(index, {"n_het_sites"});
    const int het_depth_col = optional_column(index, {"het_total_depth"});
    const int state_col = optional_column(index,
        {"runnerup_comparison_state", "runner_up_comparison_state"});
    const int softmax_col = optional_column(index,
        {"margin_softmax_score", "posterior"});
    const int entropy_col = optional_column(index,
        {"margin_entropy", "entropy"});
    const int schema_col = optional_column(index, {"schema_version"});
    const bool versioned = schema_col >= 0;

    if (versioned) {
        if (index.count("llr_vs_runner_up") == 0) {
            schema_error(filename + ": versioned diagnostics require llr_vs_runner_up");
        }
        if (state_col < 0) {
            schema_error(filename + ": versioned diagnostics require runnerup_comparison_state");
        }
        if (index.count("margin_softmax_score") == 0) {
            schema_error(filename + ": versioned diagnostics require margin_softmax_score");
        }
        if (index.count("margin_entropy") == 0) {
            schema_error(filename + ": versioned diagnostics require margin_entropy");
        }
    }
    has_margin_softmax_data = softmax_col >= 0;

    long line_no = 1;
    while (reader.next()) {
        ++line_no;
        if (reader.line == NULL || reader.line[0] == '\0') continue;
        const vector<string> fields = split_tsv_strict(reader.line);
        const string context = filename + ": line " + to_string(line_no);
        if (fields.size() != header.size()) {
            schema_error(context + ": malformed row has " +
                to_string(fields.size()) + " fields; expected " +
                to_string(header.size()));
        }
        if (versioned && !supported_schema(
                fields[schema_col], "demux_parallel_diagnostics")) {
            schema_error(context + ": unsupported diagnostics schema '" +
                fields[schema_col] + "'");
        }

        Diagnostics d;
        d.barcode = fields[barcode_col];
        d.assignment = fields[assignment_col];
        if (fields[sd_col].empty()) {
            schema_error(context + ": empty singlet_doublet value");
        }
        d.singlet_doublet = fields[sd_col][0];
        d.llr_vs_runner_up = parse_optional_diagnostic_number(
            fields[direct_col], context + " llr_vs_runner_up");
        d.runnerup_comparison_state = parse_comparison_state(
            state_col >= 0 ? fields[state_col] : "",
            d.llr_vs_runner_up, state_col >= 0,
            context + " runner-up comparison");
        d.min_margin = parse_optional_diagnostic_number(
            fields[min_margin_col], context + " min_margin");
        d.worst_competitor = fields[worst_col];
        d.n_close = static_cast<int>(parse_required_integer(
            fields[n_close_col], context + " n_close"));
        d.total_depth = parse_optional_diagnostic_number(
            fields[total_depth_col], context + " total_depth");
        if (het_var_col >= 0) {
            d.het_balance_var = parse_optional_diagnostic_number(
                fields[het_var_col], context + " het_balance_var");
        }
        if (n_het_col >= 0) {
            d.n_het_sites = parse_optional_integer(
                fields[n_het_col], context + " n_het_sites");
        }
        if (het_depth_col >= 0) {
            d.het_total_depth = parse_optional_diagnostic_number(
                fields[het_depth_col], context + " het_total_depth");
        }
        if (softmax_col >= 0) {
            d.margin_softmax_score = parse_optional_diagnostic_number(
                fields[softmax_col], context + " margin_softmax_score");
        }
        if (entropy_col >= 0) {
            d.margin_entropy = parse_optional_diagnostic_number(
                fields[entropy_col], context + " margin_entropy");
        }
        diagnostics[d.barcode] = d;
    }
}

void parse_runner_ups_file(
    const string& filename,
    map<string, vector<RunnerUp>>& runner_ups) {
    gzreader reader(filename);
    if (!reader.next()) schema_error(filename + ": empty runner-up file");

    const vector<string> header = split_tsv_strict(reader.line);
    const map<string, int> index = make_header_index(header, filename);
    const int barcode_col = required_column(index, {"barcode"}, filename);
    const int rank_col = required_column(index, {"rank"}, filename);
    const int identity_col = required_column(index, {"identity"}, filename);
    const int direct_col = required_column(index, {"llr_vs_winner"}, filename);
    const int min_margin_col = required_column(index, {"min_margin"}, filename);
    const int state_col = optional_column(index, {"comparison_state"});
    const int schema_col = optional_column(index, {"schema_version"});
    const bool versioned = schema_col >= 0;
    if (versioned && state_col < 0) {
        schema_error(filename + ": versioned runner-ups require comparison_state");
    }

    long line_no = 1;
    while (reader.next()) {
        ++line_no;
        if (reader.line == NULL || reader.line[0] == '\0') continue;
        const vector<string> fields = split_tsv_strict(reader.line);
        const string context = filename + ": line " + to_string(line_no);
        if (fields.size() != header.size()) {
            schema_error(context + ": malformed row has " +
                to_string(fields.size()) + " fields; expected " +
                to_string(header.size()));
        }
        if (versioned && !supported_schema(
                fields[schema_col], "demux_parallel_runner_ups")) {
            schema_error(context + ": unsupported runner-up schema '" +
                fields[schema_col] + "'");
        }

        RunnerUp ru;
        ru.barcode = fields[barcode_col];
        ru.rank = static_cast<int>(parse_required_integer(
            fields[rank_col], context + " rank"));
        if (ru.rank <= 0) schema_error(context + ": runner-up rank must be positive");
        ru.identity = fields[identity_col];
        if (ru.identity.empty()) schema_error(context + ": empty runner-up identity");
        ru.llr_vs_winner = parse_optional_diagnostic_number(
            fields[direct_col], context + " llr_vs_winner");
        ru.comparison_state = parse_comparison_state(
            state_col >= 0 ? fields[state_col] : "",
            ru.llr_vs_winner, state_col >= 0,
            context + " runner-up comparison");
        ru.min_margin = parse_optional_diagnostic_number(
            fields[min_margin_col], context + " min_margin");
        runner_ups[ru.barcode].push_back(ru);
    }

    for (auto& entry : runner_ups) {
        sort(entry.second.begin(), entry.second.end(),
             [](const RunnerUp& lhs, const RunnerUp& rhs) {
                 if (lhs.rank != rhs.rank) return lhs.rank < rhs.rank;
                 return lhs.identity < rhs.identity;
             });
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

double compute_depth_normalized_llr_vs_runner_up(
    double llr_vs_runner_up,
    const string& comparison_state,
    double total_depth) {
    if (!comparison_state_is_present(comparison_state) ||
        !isfinite(llr_vs_runner_up) || !isfinite(total_depth) || total_depth <= 0) {
        return NAN;
    }
    return llr_vs_runner_up / total_depth;
}

double compute_margin_ratio(double winner_min_margin, const vector<RunnerUp>* runner_ups) {
    if (!isfinite(winner_min_margin) || runner_ups == NULL || runner_ups->empty()) {
        return NAN;
    }
    const RunnerUp* complete_runner = NULL;
    for (const RunnerUp& runner : *runner_ups) {
        if (comparison_state_is_present(runner.comparison_state) &&
            isfinite(runner.min_margin)) {
            complete_runner = &runner;
            break;
        }
    }
    if (complete_runner == NULL) return NAN;
    // Both maximin margins can be negative. Use absolute values while keeping
    // missing/partial direct comparisons out of the derived metric.
    const double abs_winner = fabs(winner_min_margin);
    const double abs_runner = fabs(complete_runner->min_margin);
    if (abs_runner < 1e-10) return 1e6;  // genuine numeric zero support
    return abs_winner / abs_runner;
}

// ============================================================================
// Het Var Signal Annotation (v3: no longer drives assignments)
// ============================================================================

string compute_het_var_signal(double het_balance_var, const Thresholds& thresholds, bool has_het_data) {
    if (!isfinite(het_balance_var) || !has_het_data) return "no_data";
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
        // Unavailable/partial/not-applicable runner edges are missing evidence,
        // not zero-strength biological support for a quad pattern.
        if (!comparison_state_is_present(ru.comparison_state)) continue;
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
                   double external_ploidy_min_prob,
                   bool has_het_data,
                   RefinedCell& cell) {
    
    cell.barcode = assn.barcode;
    cell.original_assignment = assn.identity;
    cell.original_type = assn.type;
    cell.assignment_score = assn.llr;
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
    // CASE 2a: Homotypic tetraploid is expected - candidate only.
    //
    // expected_lines establishes biological eligibility (A+A is possible), but
    // it is not cell-level evidence that this particular demux singlet is
    // tetraploid. Do not relabel here. If external ploidy evidence is available,
    // the shared external-evidence branch below decides whether A -> A+A.
    // -------------------------------------------------------------------------
    
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
    // CASE 2c: Use cell-specific external ploidy when available
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
            // probability not provided; annotate but do not relabel by default.
            cell.ploidy_confidence = "LOW";
        }
        
        bool external_confident_enough = (ext_ploidy->probability >= external_ploidy_min_prob);
        if (cell.ploidy == "tetraploid" && has_homotypic_in_pool && external_confident_enough) {
            cell.refined_assignment = make_homotypic(indv);
            cell.changed = true;
        } else if (cell.ploidy == "tetraploid" && has_homotypic_in_pool && !external_confident_enough) {
            // Preserve the external ploidy annotation, but avoid biological relabeling
            // from a borderline NN call.  The audit can report this as optional
            // evidence without changing assignment.
            cell.ploidy_method = "external_low_confidence";
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
// Scoring-only pass-through (v3.3)
// ============================================================================

void score_only_cell(const Assignment& assn,
                     const Diagnostics* diag,
                     const vector<RunnerUp>* runner_ups,
                     const ExpectedLines& expected,
                     const Thresholds& thresholds,
                     bool has_het_data,
                     RefinedCell& cell) {

    (void)expected;  // unused; signature parallels classify_cell()

    cell.barcode = assn.barcode;
    cell.original_assignment = assn.identity;
    cell.original_type = assn.type;
    cell.assignment_score = assn.llr;

    // Pass-through: no reclassification
    cell.refined_assignment = assn.identity;
    cell.changed = false;

    // Fixed tokens for the ploidy/droplet columns
    cell.ploidy = "passthrough";
    cell.ploidy_method = "scoring_only";
    cell.ploidy_confidence = "N/A";
    cell.cells_in_droplet = 1;
    cell.droplet_flag = "NONE";
    cell.droplet_candidates = ".";

    // Het diagnostics annotation only
    cell.het_var_signal = "no_data";
    cell.het_balance_var = -1.0;
    if (diag != NULL && diag->het_balance_var >= 0 && diag->singlet_doublet == 'S') {
        cell.het_balance_var = diag->het_balance_var;
        cell.het_var_signal = compute_het_var_signal(diag->het_balance_var,
                                                     thresholds, has_het_data);
    }

    // Quad pattern score is informative; still compute it but don't act on it
    cell.quad_pattern_score = -1.0;
    if (assn.type == 'D' && diag != NULL && runner_ups != NULL) {
        string id1, id2;
        parse_doublet_identity(assn.identity, id1, id2);
        if (id1 != id2) {
            string dummy_candidates;
            // detect_droplet_doublet writes the pattern score as a side effect;
            // we throw away the returned label
            (void)detect_droplet_doublet(*diag, *runner_ups, thresholds,
                                         dummy_candidates, cell.quad_pattern_score);
        }
    }

    // The margin softmax is descriptive, not a calibrated posterior. In
    // scoring-only mode report the state/sign of the direct winner-runner
    // comparison instead of converting that score into confidence classes.
    if (diag == NULL) {
        cell.overall_confidence = "DIRECT_CONTRAST_UNAVAILABLE";
    } else if (!comparison_state_is_present(diag->runnerup_comparison_state)) {
        cell.overall_confidence = "DIRECT_CONTRAST_";
        string state = diag->runnerup_comparison_state;
        transform(state.begin(), state.end(), state.begin(),
            [](unsigned char c){ return static_cast<char>(std::toupper(c)); });
        cell.overall_confidence += state;
    } else if (diag->llr_vs_runner_up > 0.0) {
        cell.overall_confidence = "DIRECT_CONTRAST_POSITIVE";
    } else if (diag->llr_vs_runner_up == 0.0) {
        cell.overall_confidence = "DIRECT_CONTRAST_ZERO";
    } else {
        cell.overall_confidence = "DIRECT_CONTRAST_NONPOSITIVE";
    }
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
         << "assignment_score\t"
         << "ploidy\t"
         << "ploidy_method\t"
         << "ploidy_confidence\t"
         << "cells_in_droplet\t"
         << "droplet_flag\t"
         << "droplet_candidates\t"
         << "overall_confidence\t"
         << "changed\t"
         << "llr_vs_runner_up\t"
         << "runnerup_comparison_state\t"
         << "depth_normalized_llr_vs_runner_up\t"
         << "margin_ratio\t"
         << "margin_softmax_score\t"
         << "margin_entropy\t"
         << "contam_rate\t"
         << "het_var_signal\t"
         << "het_balance_var\t"
         << "quad_pattern_score" << endl;

    for (const auto& cell : cells) {
        outf << cell.barcode << "\t"
             << cell.original_assignment << "\t"
             << cell.original_type << "\t"
             << cell.refined_assignment << "\t"
             << cell.assignment_score << "\t"
             << cell.ploidy << "\t"
             << cell.ploidy_method << "\t"
             << cell.ploidy_confidence << "\t"
             << cell.cells_in_droplet << "\t"
             << cell.droplet_flag << "\t"
             << cell.droplet_candidates << "\t"
             << cell.overall_confidence << "\t"
             << (cell.changed ? "TRUE" : "FALSE") << "\t";

        if (isfinite(cell.llr_vs_runner_up)) {
            outf << fixed << setprecision(6) << cell.llr_vs_runner_up;
        } else {
            outf << "NA";
        }
        outf << "\t" << cell.runnerup_comparison_state << "\t";

        if (isfinite(cell.depth_normalized_llr_vs_runner_up)) {
            outf << fixed << setprecision(6)
                 << cell.depth_normalized_llr_vs_runner_up;
        } else {
            outf << "NA";
        }
        outf << "\t";

        if (isfinite(cell.margin_ratio)) {
            if (cell.margin_ratio < 1e5) {
                outf << fixed << setprecision(2) << cell.margin_ratio;
            } else {
                outf << "INF";
            }
        } else {
            outf << "NA";
        }
        outf << "\t";

        if (isfinite(cell.margin_softmax_score)) {
            outf << fixed << setprecision(6) << cell.margin_softmax_score;
        } else {
            outf << "NA";
        }
        outf << "\t";

        if (isfinite(cell.margin_entropy)) {
            outf << fixed << setprecision(4) << cell.margin_entropy;
        } else {
            outf << "NA";
        }
        outf << "\t";

        if (cell.contam_rate >= 0) {
            outf << fixed << setprecision(6) << cell.contam_rate;
        } else {
            outf << "NA";
        }
        outf << "\t" << cell.het_var_signal << "\t";

        if (isfinite(cell.het_balance_var) && cell.het_balance_var >= 0) {
            outf << fixed << setprecision(6) << cell.het_balance_var;
        } else {
            outf << "NA";
        }
        outf << "\t";

        if (cell.quad_pattern_score >= 0) {
            outf << fixed << setprecision(3) << cell.quad_pattern_score;
        } else {
            outf << "NA";
        }
        outf << endl;
    }
}

void write_simple_assignments(const string& filename, const vector<RefinedCell>& cells) {
    ofstream outf(filename);
    if (!outf.is_open()) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", filename.c_str());
        exit(1);
    }
    
    // Demux-compatible assignment-like output: no header, four columns.
    // Rich biological annotations remain in *.refined_assignments.
    for (const auto& cell : cells) {
        char sd = (cell.refined_assignment.find('+') != string::npos) ? 'D' : 'S';
        outf << cell.barcode << "\t"
             << cell.refined_assignment << "\t"
             << sd << "\t"
             << cell.assignment_score << endl;
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
                   bool has_margin_softmax_data,
                   double median_depth,
                   const SummaryStats& stats,
                   bool scoring_only) {
    
    ofstream outf(filename);
    if (!outf.is_open()) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", filename.c_str());
        exit(1);
    }
    
    outf << "tetra_refine v3.5 summary" << endl;
    outf << "=======================" << endl;
    if (scoring_only) {
        outf << "Mode: scoring-only (pass-through, no reclassification)" << endl;
    }
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
    outf << "  margin_softmax_score/margin_entropy data: " << (has_margin_softmax_data ? "YES" : "NO") << endl;
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
    
    if (!scoring_only) {
        outf << "Ploidy classification:" << endl;
        outf << "  Diploid: " << stats.ploidy_diploid << endl;
        outf << "  Tetraploid: " << stats.ploidy_tetraploid << endl;
        outf << "    - Heterotypic: " << stats.tetraploid_heterotypic << endl;
        outf << "    - Homotypic (from expected_lines): " << stats.tetraploid_from_expected << endl;
        outf << "    - Homotypic (from external ploidy): " << stats.tetraploid_from_external << endl;
        outf << "  Unknown: " << stats.ploidy_unknown << endl;
        outf << endl;
    } else {
        outf << "Pass-through cells: " << stats.passthrough_count << endl;
        outf << endl;
    }
    
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
    fprintf(stderr, "tetra_refine v3.5 [OPTIONS]\n");
    fprintf(stderr, "Refines demux_parallel assignments for tetraploid pools.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "This tool:\n");
    fprintf(stderr, "  1. Recognizes heterotypic tetraploids (A+B) as single cells, not doublets\n");
    fprintf(stderr, "  2. Reclassifies homotypic tetraploids (A+A) from expected_lines\n");
    fprintf(stderr, "  3. Flags potential doublets of tetraploids (quads)\n");
    fprintf(stderr, "  4. Computes direct winner/runner-up contrast diagnostics and margin_ratio\n");
    fprintf(stderr, "  5. Carries through margin_softmax_score, margin_entropy, and contamination rate\n");
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
    fprintf(stderr, "    --external_ploidy_min_prob F  Min confidence to relabel A -> A+A from external ploidy [0.90]\n");
    fprintf(stderr, "\n===== OPTIONAL - Threshold overrides =====\n");
    fprintf(stderr, "    --het_var_diploid F      Het variance annotation threshold for diploid\n");
    fprintf(stderr, "    --het_var_tetraploid F   Het variance annotation threshold for tetraploid\n");
    fprintf(stderr, "    --min_margin F           Min margin threshold for doublet detection\n");
    fprintf(stderr, "\n===== OPTIONAL - Output =====\n");
    fprintf(stderr, "    --write_changed_only     Also write separate file with only changed cells\n");
    fprintf(stderr, "    --write_simple           Also write simple .assignments format file\n");
    fprintf(stderr, "\n===== OPTIONAL - Scoring-only mode =====\n");
    fprintf(stderr, "    --scoring_only           Pass through assignments unchanged; emit only scoring\n");
    fprintf(stderr, "                             columns (llr_vs_runner_up, depth-normalized contrast,\n");
    fprintf(stderr, "                             margin_ratio, margin_softmax_score, margin_entropy,\n");
    fprintf(stderr, "                             contam_rate, het_balance_var, quad_pattern_score).\n");
    fprintf(stderr, "                             No ploidy reclassification, no doublet flagging.\n");
    fprintf(stderr, "                             --external_ploidy is ignored in this mode.\n");
    fprintf(stderr, "                             --expected is still required (used for pool-membership\n");
    fprintf(stderr, "                             annotations only).\n");
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
        {"external_ploidy_min_prob", required_argument, 0, 1011},
        {"het_var_diploid", required_argument, 0, 1001},
        {"het_var_tetraploid", required_argument, 0, 1002},
        {"min_margin", required_argument, 0, 1003},
        {"write_changed_only", no_argument, 0, 1004},
        {"write_simple", no_argument, 0, 1005},
        {"scoring_only", no_argument, 0, 1010},
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
    double external_ploidy_min_prob = 0.90;
    
    double het_var_diploid_override = -1.0;
    double het_var_tetraploid_override = -1.0;
    double min_margin_override = -1.0;
    
    bool write_changed_only = false;
    bool write_simple = false;
    bool verbose = false;
    bool scoring_only = false;
    
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
            case 1011: external_ploidy_min_prob = atof(optarg); break;
            case 1010: scoring_only = true; break;
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
    if (external_ploidy_min_prob < 0.0 || external_ploidy_min_prob > 1.0) {
        fprintf(stderr, "ERROR: --external_ploidy_min_prob must be between 0 and 1\n");
        exit(1);
    }

    // Scoring-only mode: --external_ploidy is ignored (warn and clear)
    if (scoring_only && !external_ploidy_file.empty()) {
        fprintf(stderr, "WARNING: --external_ploidy ignored in --scoring_only mode\n");
        external_ploidy_file = "";
    }

    fprintf(stderr, "tetra_refine v3.5\n");
    fprintf(stderr, "===============\n");
    if (verbose && scoring_only) {
        fprintf(stderr, "Mode: scoring-only (no reclassification)\n");
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
    bool has_margin_softmax_data = false;
    parse_diagnostics_file(diagnostics_file, diagnostics, has_margin_softmax_data);
    if (verbose) {
        fprintf(stderr, "  Loaded %lu diagnostics\n", diagnostics.size());
        fprintf(stderr, "  margin_softmax_score/margin_entropy columns: %s\n", has_margin_softmax_data ? "YES" : "NO");
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
        if (isfinite(dp.second.total_depth)) {
            depths.push_back(dp.second.total_depth);
        }

        if (dp.second.singlet_doublet == 'D' && isfinite(dp.second.min_margin)) {
            min_margins_d.push_back(dp.second.min_margin);
        }

        if (dp.second.singlet_doublet == 'S' &&
            isfinite(dp.second.het_balance_var) && dp.second.het_balance_var >= 0) {
            het_vars_singlet.push_back(dp.second.het_balance_var);
            has_het_data = true;
        }
    }
    
    double median_depth = !depths.empty() ? compute_median(depths) : 0.0;
    
    Thresholds thresholds;
    thresholds.source = "";

    if (scoring_only) {
        thresholds.min_margin_threshold = 0.0;
        thresholds.min_margin_possible = 0.0;
        thresholds.het_var_diploid = 0.0;
        thresholds.het_var_tetraploid = 0.0;
        thresholds.close_count_threshold = 1;
        thresholds.depth_ratio_threshold = 1.5;
        thresholds.source = "scoring_only";
    } else {
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
    }
    
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
        if (scoring_only) {
            score_only_cell(assn, diag_ptr, runner_ptr, expected, thresholds,
                            has_het_data, cell);
        } else {
            classify_cell(assn, diag_ptr, runner_ptr, ext_ptr, expected, thresholds,
                          median_depth, external_ploidy_min_prob, has_het_data, cell);
        }
        
        // Compute confidence/diagnostic metrics from the direct runner-up
        // comparison. margin_softmax_score is carried through descriptively and
        // is never treated as a calibrated posterior.
        if (diag_ptr != NULL) {
            cell.llr_vs_runner_up = diag_ptr->llr_vs_runner_up;
            cell.runnerup_comparison_state = diag_ptr->runnerup_comparison_state;
            cell.depth_normalized_llr_vs_runner_up =
                compute_depth_normalized_llr_vs_runner_up(
                    diag_ptr->llr_vs_runner_up,
                    diag_ptr->runnerup_comparison_state,
                    diag_ptr->total_depth);
            cell.margin_ratio = compute_margin_ratio(diag_ptr->min_margin, runner_ptr);
            cell.margin_softmax_score = diag_ptr->margin_softmax_score;
            cell.margin_entropy = diag_ptr->margin_entropy;
        } else {
            cell.llr_vs_runner_up = NAN;
            cell.runnerup_comparison_state = "unavailable";
            cell.depth_normalized_llr_vs_runner_up = NAN;
            cell.margin_ratio = NAN;
            cell.margin_softmax_score = NAN;
            cell.margin_entropy = NAN;
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
        
        if (cell.ploidy == "passthrough") stats.passthrough_count++;
        else if (cell.ploidy == "diploid") stats.ploidy_diploid++;
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
                  expected, thresholds, has_het_data, has_margin_softmax_data, median_depth, stats,
                  scoring_only);
    
    // Print summary to stderr
    fprintf(stderr, "\ntetra_refine v3.5 complete\n");
    fprintf(stderr, "========================\n");
    fprintf(stderr, "Total cells: %d\n", stats.total_cells);
    if (!scoring_only) {
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
    } else {
        fprintf(stderr, "\nPass-through cells: %d\n", stats.passthrough_count);
    }
    fprintf(stderr, "\nData sources:\n");
    fprintf(stderr, "  margin_softmax_score/margin_entropy: %s\n", has_margin_softmax_data ? "available" : "not in diagnostics");
    fprintf(stderr, "  Contamination rate: %s\n", !contam_rate_file.empty() ? "loaded" : "not provided");
    fprintf(stderr, "  External ploidy: %s\n", !external_ploidy_file.empty() ? "loaded" : "not provided");
    fprintf(stderr, "\nOutput files:\n");
    fprintf(stderr, "  %s\n", refined_file.c_str());
    if (write_changed_only) fprintf(stderr, "  %s.refined_changed\n", output_prefix.c_str());
    if (write_simple) fprintf(stderr, "  %s.assignments_refined\n", output_prefix.c_str());
    fprintf(stderr, "  %s\n", summary_file.c_str());
    
    return 0;
}
