#include <algorithm>
#include <cctype>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <limits>
#include <initializer_list>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>

// CellBouncer / htswrapper barcode encoder used by .counts files.
#include <htswrapper/bc.h>

using namespace std;

struct Identity {
    string name;
    int a = -1;
    int b = -1;   // -1 for singlet; a==b for homotypic
    string type = "S";
};

struct ScoreAccum {
    double abs_error_weighted = 0.0;
    double depth = 0.0;
    long bins = 0;
};

struct Diagnostics {
    double llr_vs_runner_up = NAN;
    string runnerup_comparison_state = "unavailable";
    double total_depth = NAN;
    long n_close = 0;
    double margin_softmax_score = NAN;
};

struct CellInfo {
    string barcode;
    unsigned long bc_ul_val = 0;
    Identity assignment;
    double assignment_llr = NAN;
    vector<Identity> runnerups;
    vector<string> runnerup_names;
    vector<string> runnerup_comparison_states;
    ScoreAccum assigned_acc;
    vector<ScoreAccum> runner_acc;
    Diagnostics diag;

    bool has_expected_species = false;
    Identity expected_species_identity;
    ScoreAccum species_expected_acc;
    vector<ScoreAccum> species_candidate_acc;
};

static string trim(const string& s){
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b-a+1);
}

static void die(const string& msg);

static vector<string> split(const string& s, char delim){
    vector<string> out;
    string item;
    stringstream ss(s);
    while (getline(ss, item, delim)) out.push_back(item);
    return out;
}

static vector<string> split_tsv_strict(const string& s){
    vector<string> out;
    size_t start = 0;
    while (true){
        size_t pos = s.find('\t', start);
        if (pos == string::npos){
            out.push_back(s.substr(start));
            break;
        }
        out.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return out;
}

static string lowercase(string s){
    transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    return s;
}

static bool is_missing_diagnostic_token(const string& raw){
    const string token = lowercase(trim(raw));
    return token.empty() || token == "." || token == "na" ||
           token == "unavailable" || token == "partial_support" ||
           token == "not_applicable";
}

static double parse_optional_diagnostic_number(
    const string& raw, const string& context){
    if (is_missing_diagnostic_token(raw)) return NAN;
    errno = 0;
    char* end = NULL;
    const double value = strtod(raw.c_str(), &end);
    if (errno != 0 || end == raw.c_str() || *end != '\0' || !isfinite(value)){
        die(context + ": expected a finite number or an explicit missing state, saw '" + raw + "'");
    }
    return value;
}

static long parse_required_long(const string& raw, const string& context){
    errno = 0;
    char* end = NULL;
    const long value = strtol(raw.c_str(), &end, 10);
    if (errno != 0 || end == raw.c_str() || *end != '\0'){
        die(context + ": expected an integer, saw '" + raw + "'");
    }
    return value;
}

static map<string,int> diagnostic_header_index(
    const vector<string>& header, const string& path){
    map<string,int> index;
    for (int i = 0; i < (int)header.size(); ++i){
        if (header[i].empty()) die(path + ": empty diagnostic column name");
        if (index.count(header[i]) > 0) die(path + ": duplicate diagnostic column: " + header[i]);
        index[header[i]] = i;
    }
    return index;
}

static int optional_column(const map<string,int>& index,
                           initializer_list<const char*> names){
    for (const char* name : names){
        auto it = index.find(name);
        if (it != index.end()) return it->second;
    }
    return -1;
}

static int required_column(const map<string,int>& index,
                           initializer_list<const char*> names,
                           const string& path){
    const int col = optional_column(index, names);
    if (col >= 0) return col;
    string wanted;
    for (const char* name : names){
        if (!wanted.empty()) wanted += " or ";
        wanted += name;
    }
    die(path + ": required diagnostic column missing: " + wanted);
    return -1;
}

static bool comparison_state_is_present(const string& state){
    return state == "present_nonzero" || state == "present_zero";
}

static string parse_comparison_state(
    const string& raw_state,
    double numeric_value,
    bool explicit_state,
    const string& context){
    if (!explicit_state){
        if (!isfinite(numeric_value)) return "unavailable";
        return numeric_value == 0.0 ? "present_zero" : "present_nonzero";
    }
    const string state = lowercase(trim(raw_state));
    static const set<string> supported = {
        "present_nonzero", "present_zero", "unavailable",
        "partial_support", "not_applicable"
    };
    if (supported.count(state) == 0){
        die(context + ": unsupported comparison state '" + raw_state + "'");
    }
    if (comparison_state_is_present(state)){
        if (!isfinite(numeric_value)){
            die(context + ": state " + state + " requires a numeric direct comparison");
        }
        if (state == "present_zero" && numeric_value != 0.0){
            die(context + ": present_zero requires numeric zero");
        }
        if (state == "present_nonzero" && numeric_value == 0.0){
            die(context + ": present_nonzero cannot carry numeric zero");
        }
    } else if (isfinite(numeric_value)){
        die(context + ": missing comparison state " + state +
            " must not carry a numeric value");
    }
    return state;
}

static bool supported_schema(const string& schema,
                             const string& prefix){
    return schema == prefix + "_v2" || schema == prefix + "_v3";
}

static bool file_exists(const string& path){
    if (path.empty()) return false;
    ifstream in(path.c_str());
    return (bool)in;
}

static void die(const string& msg){
    fprintf(stderr, "ERROR: %s\n", msg.c_str());
    exit(1);
}

static string canonical_pair_name(const string& a, const string& b){
    if (a == b) return a + "+" + b;
    return (a < b) ? (a + "+" + b) : (b + "+" + a);
}

static Identity make_identity_from_indices(const string& a_name, int a_idx, const string& b_name="", int b_idx=-1){
    Identity id;
    if (b_idx < 0 || b_name.empty() || a_name == b_name){
        id.name = a_name;
        id.a = a_idx;
        id.b = -1;
        id.type = "S";
        return id;
    }
    if (a_name < b_name){
        id.name = a_name + "+" + b_name;
        id.a = a_idx;
        id.b = b_idx;
    } else {
        id.name = b_name + "+" + a_name;
        id.a = b_idx;
        id.b = a_idx;
    }
    id.type = "D";
    return id;
}

static Identity parse_identity(const string& raw, const unordered_map<string,int>& sample2idx){
    string s = trim(raw);
    s.erase(remove_if(s.begin(), s.end(), [](unsigned char ch){ return std::isspace(ch); }), s.end());
    if (s.empty()) die("empty identity");
    size_t p = s.find('+');
    Identity id;
    if (p == string::npos){
        auto it = sample2idx.find(s);
        if (it == sample2idx.end()) die("identity component not found in sample vector: " + s);
        id.name = s;
        id.a = it->second;
        id.b = -1;
        id.type = "S";
        return id;
    }
    if (s.find('+', p+1) != string::npos) die("malformed identity with multiple '+': " + s);
    string x = s.substr(0, p);
    string y = s.substr(p+1);
    if (x.empty() || y.empty()) die("malformed identity with empty component: " + s);
    auto ix = sample2idx.find(x);
    auto iy = sample2idx.find(y);
    if (ix == sample2idx.end()) die("identity component not found in sample vector: " + x);
    if (iy == sample2idx.end()) die("identity component not found in sample vector: " + y);
    if (x == y){
        id.name = x + "+" + y;
        id.a = ix->second;
        id.b = ix->second;
        id.type = "H";
    }
    else if (x < y){
        id.name = x + "+" + y;
        id.a = ix->second;
        id.b = iy->second;
        id.type = "D";
    }
    else{
        id.name = y + "+" + x;
        id.a = iy->second;
        id.b = ix->second;
        id.type = "D";
    }
    return id;
}

static vector<string> load_samples(const string& path){
    ifstream in(path.c_str());
    if (!in) die("could not open samples: " + path);
    vector<string> samples;
    string s;
    while (in >> s) samples.push_back(s);
    if (samples.empty()) die("empty samples file: " + path);
    return samples;
}

static unordered_map<string,string> load_panel(const string& path){
    unordered_map<string,string> panel;
    if (path.empty()) return panel;
    ifstream in(path.c_str());
    if (!in) return panel;
    string header;
    if (!getline(in, header)) return panel;
    vector<string> h = split(header, '\t');
    int id_col = -1, sp_col = -1;
    for (int i=0; i<(int)h.size(); ++i){
        if (h[i] == "indiv_id" || h[i] == "VCF_ID" || h[i] == "sample" || h[i] == "Sample") id_col = i;
        if (h[i] == "species" || h[i] == "Species") sp_col = i;
    }
    if (id_col < 0) id_col = 0;
    if (sp_col < 0) return panel;
    string line;
    while (getline(in, line)){
        if (line.empty()) continue;
        vector<string> f = split(line, '\t');
        if ((int)f.size() > max(id_col, sp_col)) panel[trim(f[id_col])] = trim(f[sp_col]);
    }
    return panel;
}

static void add_expanded_species(set<string>& sp, const string& label){
    if (label == "Hy" || label == "Chinobo-mCherry"){
        sp.insert("B");
        sp.insert("C");
    } else if (!label.empty() && label != "NA") {
        sp.insert(label);
    }
}

static string join_species_set(const set<string>& sp){
    string out;
    for (auto it=sp.begin(); it!=sp.end(); ++it){
        if (!out.empty()) out += ",";
        out += *it;
    }
    return out.empty() ? "NA" : out;
}

static set<string> parse_species_set_string(const string& s){
    set<string> out;
    string token;
    auto flush = [&](){
        string t = trim(token);
        if (!t.empty() && t != "NA") add_expanded_species(out, t);
        token.clear();
    };
    for (char ch : s){
        if (ch == ',' || ch == '+' || ch == ';') flush();
        else token.push_back(ch);
    }
    flush();
    return out;
}

static bool is_subset_of(const set<string>& a, const set<string>& b){
    for (const auto& x : a){
        if (b.find(x) == b.end()) return false;
    }
    return true;
}

static bool intersects_set(const set<string>& a, const set<string>& b){
    for (const auto& x : a){
        if (b.find(x) != b.end()) return true;
    }
    return false;
}

static string species_relation(const string& expected, const string& observed){
    set<string> exp = parse_species_set_string(expected);
    set<string> obs = parse_species_set_string(observed);
    if (exp.empty() || obs.empty()) return "missing_species_evidence";
    if (exp == obs) return "exact_match";
    if (is_subset_of(obs, exp)) return "expected_subset_only_component_missing";
    if (is_subset_of(exp, obs)) return "expected_superset_with_extra_species";
    if (intersects_set(exp, obs)) return "partial_overlap_with_extra_and_missing";
    return "disjoint_wrong_species";
}

static string expected_species_set(const Identity& id, const vector<string>& samples, const unordered_map<string,string>& panel){
    if (panel.empty()) return "NA";
    set<string> sp;
    auto add = [&](int idx){
        if (idx >= 0 && idx < (int)samples.size()){
            auto it = panel.find(samples[idx]);
            string label = it == panel.end() ? (samples[idx] == "Chinobo-mCherry" ? "Chinobo-mCherry" : "UNKNOWN") : it->second;
            add_expanded_species(sp, label);
        }
    };
    add(id.a);
    if (id.b >= 0) add(id.b);
    return join_species_set(sp);
}

static bool make_expected_species_identity(const Identity& indiv_id,
                                           const vector<string>& indiv_samples,
                                           const unordered_map<string,string>& panel,
                                           const unordered_map<string,int>& species2idx,
                                           Identity& out){
    if (panel.empty() || species2idx.empty()) return false;
    set<string> sp_set;
    auto add = [&](int idx){
        if (idx < 0 || idx >= (int)indiv_samples.size()) return;
        auto it = panel.find(indiv_samples[idx]);
        string sp = (it == panel.end()) ? (indiv_samples[idx] == "Chinobo-mCherry" ? "Chinobo-mCherry" : "UNKNOWN") : it->second;
        add_expanded_species(sp_set, sp);
    };
    add(indiv_id.a);
    if (indiv_id.b >= 0) add(indiv_id.b);
    if (sp_set.empty()) return false;
    vector<string> species(sp_set.begin(), sp_set.end());
    if (species.size() == 1){
        auto it = species2idx.find(species[0]);
        if (it == species2idx.end()) return false;
        out = make_identity_from_indices(species[0], it->second);
        return true;
    }
    if (species.size() == 2){
        auto ia = species2idx.find(species[0]);
        auto ib = species2idx.find(species[1]);
        if (ia == species2idx.end() || ib == species2idx.end()) return false;
        out = make_identity_from_indices(species[0], ia->second, species[1], ib->second);
        return true;
    }
    return false;
}

static vector<Identity> build_species_candidates(const vector<string>& species_samples){
    vector<Identity> out;
    for (int i=0; i<(int)species_samples.size(); ++i){
        out.push_back(make_identity_from_indices(species_samples[i], i));
    }
    for (int i=0; i<(int)species_samples.size(); ++i){
        for (int j=i+1; j<(int)species_samples.size(); ++j){
            out.push_back(make_identity_from_indices(species_samples[i], i, species_samples[j], j));
        }
    }
    return out;
}

static string infer_species_samples_path(const string& species_counts){
    const string suffix = ".species_counts";
    if (species_counts.size() >= suffix.size() && species_counts.substr(species_counts.size()-suffix.size()) == suffix){
        return species_counts.substr(0, species_counts.size()-suffix.size()) + ".species_samples";
    }
    return "";
}

static void load_assignments(const string& path, const unordered_map<string,int>& sample2idx,
                             unordered_map<unsigned long, CellInfo>& cells,
                             vector<unsigned long>& order){
    ifstream in(path.c_str());
    if (!in) die("could not open assignments: " + path);
    string bc, ident, sd;
    double llr;
    while (in >> bc >> ident >> sd >> llr){
        unsigned long ul = bc_ul(bc);
        CellInfo ci;
        ci.barcode = bc;
        ci.bc_ul_val = ul;
        ci.assignment = parse_identity(ident, sample2idx);
        ci.assignment_llr = llr;
        cells[ul] = ci;
        order.push_back(ul);
    }
}

static void load_diagnostics(const string& path, unordered_map<unsigned long, CellInfo>& cells){
    gzFile gz = gzopen(path.c_str(), "rb");
    if (!gz) die("could not open diagnostics: " + path);
    char buf[1<<20];
    if (!gzgets(gz, buf, sizeof(buf))){
        gzclose(gz);
        die(path + ": empty diagnostics file");
    }
    string header_line(buf);
    header_line.erase(remove(header_line.begin(), header_line.end(), '\n'), header_line.end());
    header_line.erase(remove(header_line.begin(), header_line.end(), '\r'), header_line.end());
    const vector<string> header = split_tsv_strict(header_line);
    const map<string,int> index = diagnostic_header_index(header, path);
    const int bc_i = required_column(index, {"barcode"}, path);
    const int llr_i = required_column(index, {"llr_vs_runner_up", "llr"}, path);
    const int nc_i = required_column(index, {"n_close"}, path);
    const int td_i = required_column(index, {"total_depth"}, path);
    const int state_i = optional_column(index,
        {"runnerup_comparison_state", "runner_up_comparison_state"});
    const int softmax_i = optional_column(index,
        {"margin_softmax_score", "posterior"});
    const int schema_i = optional_column(index, {"schema_version"});
    const bool versioned = schema_i >= 0;
    if (versioned){
        if (index.count("llr_vs_runner_up") == 0){
            gzclose(gz);
            die(path + ": versioned diagnostics require llr_vs_runner_up");
        }
        if (state_i < 0){
            gzclose(gz);
            die(path + ": versioned diagnostics require runnerup_comparison_state");
        }
        if (index.count("margin_softmax_score") == 0){
            gzclose(gz);
            die(path + ": versioned diagnostics require margin_softmax_score");
        }
    }

    long line_no = 1;
    while (gzgets(gz, buf, sizeof(buf))){
        ++line_no;
        string line(buf);
        line.erase(remove(line.begin(), line.end(), '\n'), line.end());
        line.erase(remove(line.begin(), line.end(), '\r'), line.end());
        if (line.empty()) continue;
        const vector<string> f = split_tsv_strict(line);
        const string context = path + ": line " + to_string(line_no);
        if (f.size() != header.size()){
            gzclose(gz);
            die(context + ": malformed row has " + to_string(f.size()) +
                " fields; expected " + to_string(header.size()));
        }
        if (versioned && !supported_schema(
                f[schema_i], "demux_parallel_diagnostics")){
            gzclose(gz);
            die(context + ": unsupported diagnostics schema '" + f[schema_i] + "'");
        }
        const double direct = parse_optional_diagnostic_number(f[llr_i],
            context + " llr_vs_runner_up");
        const string state = parse_comparison_state(
            state_i >= 0 ? f[state_i] : "", direct, state_i >= 0,
            context + " runner-up comparison");
        const long n_close = parse_required_long(f[nc_i], context + " n_close");
        const double total_depth = parse_optional_diagnostic_number(
            f[td_i], context + " total_depth");
        const double softmax = softmax_i >= 0
            ? parse_optional_diagnostic_number(
                f[softmax_i], context + " margin_softmax_score")
            : NAN;

        string barcode = f[bc_i];
        const unsigned long ul = bc_ul(barcode);
        auto it = cells.find(ul);
        if (it == cells.end()) continue;
        it->second.diag.llr_vs_runner_up = direct;
        it->second.diag.runnerup_comparison_state = state;
        it->second.diag.n_close = n_close;
        it->second.diag.total_depth = total_depth;
        it->second.diag.margin_softmax_score = softmax;
    }
    if (gzclose(gz) != Z_OK) die("failed closing diagnostics: " + path);
}

static void load_runnerups(const string& path, const unordered_map<string,int>& sample2idx,
                           unordered_map<unsigned long, CellInfo>& cells){
    if (path.empty()) return;
    gzFile gz = gzopen(path.c_str(), "rb");
    if (!gz) die("could not open runner-ups: " + path);
    char buf[1<<20];
    if (!gzgets(gz, buf, sizeof(buf))){
        gzclose(gz);
        die(path + ": empty runner-up file");
    }
    string header_line(buf);
    header_line.erase(remove(header_line.begin(), header_line.end(), '\n'), header_line.end());
    header_line.erase(remove(header_line.begin(), header_line.end(), '\r'), header_line.end());
    const vector<string> header = split_tsv_strict(header_line);
    const map<string,int> index = diagnostic_header_index(header, path);
    const int bc_i = required_column(index, {"barcode"}, path);
    const int id_i = required_column(index, {"identity"}, path);
    const int rank_i = required_column(index, {"rank"}, path);
    const int value_i = required_column(index, {"llr_vs_winner"}, path);
    const int state_i = optional_column(index, {"comparison_state"});
    const int schema_i = optional_column(index, {"schema_version"});
    const bool versioned = schema_i >= 0;
    if (versioned && state_i < 0){
        gzclose(gz);
        die(path + ": versioned runner-ups require comparison_state");
    }

    long line_no = 1;
    while (gzgets(gz, buf, sizeof(buf))){
        ++line_no;
        string line(buf);
        line.erase(remove(line.begin(), line.end(), '\n'), line.end());
        line.erase(remove(line.begin(), line.end(), '\r'), line.end());
        if (line.empty()) continue;
        const vector<string> f = split_tsv_strict(line);
        const string context = path + ": line " + to_string(line_no);
        if (f.size() != header.size()){
            gzclose(gz);
            die(context + ": malformed row has " + to_string(f.size()) +
                " fields; expected " + to_string(header.size()));
        }
        if (versioned && !supported_schema(
                f[schema_i], "demux_parallel_runner_ups")){
            gzclose(gz);
            die(context + ": unsupported runner-up schema '" + f[schema_i] + "'");
        }
        const long rank = parse_required_long(f[rank_i], context + " rank");
        if (rank <= 0){
            gzclose(gz);
            die(context + ": runner-up rank must be positive");
        }
        const double direct = parse_optional_diagnostic_number(
            f[value_i], context + " llr_vs_winner");
        const string state = parse_comparison_state(
            state_i >= 0 ? f[state_i] : "", direct, state_i >= 0,
            context + " runner-up comparison");
        Identity r;
        try{
            r = parse_identity(f[id_i], sample2idx);
        } catch (const exception& e) {
            gzclose(gz);
            die(context + ": invalid runner-up identity '" + f[id_i] + "': " + e.what());
        }

        string barcode = f[bc_i];
        const unsigned long ul = bc_ul(barcode);
        auto it = cells.find(ul);
        if (it == cells.end()) continue;
        if (r.name == it->second.assignment.name) continue;
        bool dup = false;
        for (const auto& existing : it->second.runnerups){
            if (existing.name == r.name){ dup = true; break; }
        }
        if (!dup){
            it->second.runnerups.push_back(r);
            it->second.runnerup_names.push_back(r.name);
            it->second.runnerup_comparison_states.push_back(state);
        }
    }
    if (gzclose(gz) != Z_OK) die("failed closing runner-ups: " + path);
    for (auto& kv : cells) kv.second.runner_acc.resize(kv.second.runnerups.size());
}

static bool expected_for_record(const Identity& id, int indv1, int type1, int indv2, int type2, double& p){
    if (id.b < 0 || id.a == id.b){
        if (indv2 == -1 && indv1 == id.a){
            p = type1 / 2.0;
            return true;
        }
        return false;
    }
    if (indv2 == -1) return false;
    if (indv1 == id.a && indv2 == id.b){
        p = (type1 + type2) / 4.0;
        return true;
    }
    if (indv1 == id.b && indv2 == id.a){
        p = (type1 + type2) / 4.0;
        return true;
    }
    return false;
}

static double adjust_expected(double p, double e_ref, double e_alt){
    return p - p * e_alt + (1.0 - p) * e_ref;
}

static void update_score(ScoreAccum& acc, double expected_p, double ref, double alt, double e_ref, double e_alt){
    double d = ref + alt;
    if (d <= 0) return;
    double obs = alt / d;
    double pe = adjust_expected(expected_p, e_ref, e_alt);
    acc.abs_error_weighted += d * fabs(obs - pe);
    acc.depth += d;
    acc.bins += 1;
}

static void stream_counts(const string& path, unordered_map<unsigned long, CellInfo>& cells,
                          double e_ref, double e_alt){
    gzFile gz = gzopen(path.c_str(), "rb");
    if (!gz) die("could not open counts: " + path);
    char buf[1<<20];
    while (gzgets(gz, buf, sizeof(buf))){
        string line(buf);
        if (line.empty()) continue;
        vector<string> f = split(line, '\t');
        if (f.size() < 7) continue;
        unsigned long cell = strtoul(f[0].c_str(), nullptr, 10);
        auto it = cells.find(cell);
        if (it == cells.end()) continue;
        int indv1 = atoi(f[1].c_str());
        int type1 = atoi(f[2].c_str());
        int indv2 = atoi(f[3].c_str());
        int type2 = atoi(f[4].c_str());
        double ref = atof(f[5].c_str());
        double alt = atof(f[6].c_str());
        double p = 0.0;
        if (expected_for_record(it->second.assignment, indv1, type1, indv2, type2, p)){
            update_score(it->second.assigned_acc, p, ref, alt, e_ref, e_alt);
        }
        for (size_t i=0; i<it->second.runnerups.size(); ++i){
            if (expected_for_record(it->second.runnerups[i], indv1, type1, indv2, type2, p)){
                update_score(it->second.runner_acc[i], p, ref, alt, e_ref, e_alt);
            }
        }
    }
    gzclose(gz);
}

static void stream_species_counts(const string& path,
                                  unordered_map<unsigned long, CellInfo>& cells,
                                  const vector<Identity>& candidates,
                                  int n_species,
                                  double e_ref,
                                  double e_alt){
    gzFile gz = gzopen(path.c_str(), "rb");
    if (!gz) die("could not open species_counts: " + path);
    char buf[1<<20];
    long line_no = 0;
    while (gzgets(gz, buf, sizeof(buf))){
        ++line_no;
        string line(buf);
        if (line.empty()) continue;
        vector<string> f = split(line, '\t');
        if (f.size() < 7) continue;
        unsigned long cell = strtoul(f[0].c_str(), nullptr, 10);
        int indv1 = atoi(f[1].c_str());
        int type1 = atoi(f[2].c_str());
        int indv2 = atoi(f[3].c_str());
        int type2 = atoi(f[4].c_str());
        if (indv1 < 0 || indv1 >= n_species || (indv2 >= 0 && indv2 >= n_species)){
            stringstream ss;
            ss << "species_counts dimensional guard failed at line " << line_no
               << ": saw index " << indv1 << "," << indv2
               << " but n_species=" << n_species << ". This is not a native species-shaped file.";
            die(ss.str());
        }
        auto it = cells.find(cell);
        if (it == cells.end()) continue;
        double ref = atof(f[5].c_str());
        double alt = atof(f[6].c_str());
        double p = 0.0;
        if (it->second.has_expected_species && expected_for_record(it->second.expected_species_identity, indv1, type1, indv2, type2, p)){
            update_score(it->second.species_expected_acc, p, ref, alt, e_ref, e_alt);
        }
        if (!candidates.empty()){
            if (it->second.species_candidate_acc.size() != candidates.size()){
                it->second.species_candidate_acc.resize(candidates.size());
            }
            for (size_t i=0; i<candidates.size(); ++i){
                if (expected_for_record(candidates[i], indv1, type1, indv2, type2, p)){
                    update_score(it->second.species_candidate_acc[i], p, ref, alt, e_ref, e_alt);
                }
            }
        }
    }
    gzclose(gz);
}

static double concordance(const ScoreAccum& acc, long min_evidence){
    if (acc.depth < min_evidence || acc.depth <= 0) return NAN;
    double err = acc.abs_error_weighted / acc.depth;
    double c = 1.0 - err;
    if (c < 0) c = 0;
    if (c > 1) c = 1;
    return c;
}

static string fmt(double x){
    if (!isfinite(x)) return "NA";
    char b[64];
    snprintf(b, sizeof(b), "%.6g", x);
    return string(b);
}

static string join_flags(const vector<string>& flags){
    string out;
    for (size_t i=0; i<flags.size(); ++i){
        if (i) out += ",";
        out += flags[i];
    }
    return out;
}

static double median(vector<double> vals){
    vector<double> v;
    for (double x : vals) if (isfinite(x)) v.push_back(x);
    if (v.empty()) return NAN;
    sort(v.begin(), v.end());
    size_t n = v.size();
    if (n % 2) return v[n/2];
    return 0.5 * (v[n/2 - 1] + v[n/2]);
}

static void usage(){
    fprintf(stderr,
        "tetra_score_calls --counts FILE --samples FILE --assignments FILE --diagnostics FILE --output FILE [options]\n"
        "Options:\n"
        "  --runner_ups FILE              enables constrained dosage gap\n"
        "  --panel_metadata FILE          maps individual assignments to expected species labels\n"
        "  --species_counts FILE          native species-shaped counts; enables real species support scoring\n"
        "  --species_condf FILE           accepted for manifest compatibility; dimensional checks use species_counts\n"
        "  --species_samples FILE         species sample list; default inferred from --species_counts prefix\n"
        "  --species_support_threshold X  default 0.70; used for species support QC\n"
        "  --condf FILE                   accepted for manifest compatibility\n"
        "  --libname STR                  output library label\n"
        "  --threads N                    accepted; current streaming scorer is deterministic single-writer\n"
        "  --min_evidence INT             default 10\n"
        "  --error_ref FLOAT              default 0.001\n"
        "  --error_alt FLOAT              default 0.001\n"
        "  --strict | --best-effort       default best-effort\n");
}

int main(int argc, char** argv){
    string counts, samples_path, assignments, diagnostics, output, runnerups, panel_path, libname="NA";
    string species_counts, species_condf, species_samples_path, condf;
    long min_evidence = 10;
    double e_ref = 0.001, e_alt = 0.001;
    double species_support_threshold = 0.70;
    bool strict = false;
    for (int i=1; i<argc; ++i){
        string a = argv[i];
        auto need = [&](string& dest){ if (i+1 >= argc) die("missing value after " + a); dest = argv[++i]; };
        if (a == "--counts") need(counts);
        else if (a == "--samples") need(samples_path);
        else if (a == "--assignments") need(assignments);
        else if (a == "--diagnostics") need(diagnostics);
        else if (a == "--output") need(output);
        else if (a == "--runner_ups") need(runnerups);
        else if (a == "--panel_metadata") need(panel_path);
        else if (a == "--species_counts") need(species_counts);
        else if (a == "--species_condf") need(species_condf);
        else if (a == "--species_samples") need(species_samples_path);
        else if (a == "--condf") need(condf);
        else if (a == "--libname") need(libname);
        else if (a == "--threads") { string tmp; need(tmp); }
        else if (a == "--min_evidence") { string tmp; need(tmp); min_evidence = atol(tmp.c_str()); }
        else if (a == "--error_ref") { string tmp; need(tmp); e_ref = atof(tmp.c_str()); }
        else if (a == "--error_alt") { string tmp; need(tmp); e_alt = atof(tmp.c_str()); }
        else if (a == "--species_support_threshold") { string tmp; need(tmp); species_support_threshold = atof(tmp.c_str()); }
        else if (a == "--strict") strict = true;
        else if (a == "--best-effort") strict = false;
        else if (a == "--help" || a == "-h") { usage(); return 0; }
        else die("unknown argument: " + a);
    }
    if (counts.empty() || samples_path.empty() || assignments.empty() || diagnostics.empty() || output.empty()){
        usage();
        return 1;
    }

    vector<string> samples = load_samples(samples_path);
    unordered_map<string,int> sample2idx;
    for (int i=0; i<(int)samples.size(); ++i) sample2idx[samples[i]] = i;
    unordered_map<string,string> panel = load_panel(panel_path);
    unordered_map<unsigned long, CellInfo> cells;
    vector<unsigned long> order;
    load_assignments(assignments, sample2idx, cells, order);
    load_diagnostics(diagnostics, cells);
    if (!runnerups.empty()) load_runnerups(runnerups, sample2idx, cells);
    stream_counts(counts, cells, e_ref, e_alt);

    bool species_scoring_requested = !species_counts.empty() || !species_condf.empty() || !species_samples_path.empty();
    bool species_scoring_enabled = false;
    vector<string> species_samples;
    unordered_map<string,int> species2idx;
    vector<Identity> species_candidates;
    string species_disable_reason;
    if (species_scoring_requested){
        if (species_counts.empty()){
            species_disable_reason = "SPECIES_COUNTS_MISSING";
            if (strict) die("--species_counts is required when species scoring inputs are requested");
        } else {
            if (species_samples_path.empty()) species_samples_path = infer_species_samples_path(species_counts);
            if (species_samples_path.empty() || !file_exists(species_samples_path)){
                species_disable_reason = "SPECIES_SAMPLES_MISSING";
                if (strict) die("native species scoring requires --species_samples or an inferable .species_samples file");
            } else if (panel_path.empty() || panel.empty()){
                species_disable_reason = "PANEL_METADATA_MISSING";
                if (strict) die("species support scoring requires --panel_metadata to map individual assignments to species labels");
            } else {
                species_samples = load_samples(species_samples_path);
                for (int i=0; i<(int)species_samples.size(); ++i){
                    species2idx[species_samples[i]] = i;
                }
                species_candidates = build_species_candidates(species_samples);
                for (auto& kv : cells){
                    Identity spid;
                    if (make_expected_species_identity(kv.second.assignment, samples, panel, species2idx, spid)){
                        kv.second.has_expected_species = true;
                        kv.second.expected_species_identity = spid;
                    }
                    kv.second.species_candidate_acc.resize(species_candidates.size());
                }
                stream_species_counts(species_counts, cells, species_candidates, (int)species_samples.size(), e_ref, e_alt);
                species_scoring_enabled = true;
            }
        }
    }

    gzFile out = gzopen(output.c_str(), "wb");
    if (!out) die("could not open output: " + output);
    string header = "barcode\tlibname\tassignment\tassignment_type\tploidy_status\tllr_vs_runner_up\trunnerup_comparison_state\tmargin_softmax_score\ttotal_depth\tn_informative_bins\tn_informative_depth\tn_close\tdepth_normalized_llr_vs_runner_up\tdosage_concordance\tdosage_runnerup_identity\tdosage_runnerup_comparison_state\trunnerup_dosage_concordance\tdosage_gap_constrained\tresidual_mismatch\texpected_species_set\tspecies_support_expected\tspecies_conflict_flag\tspecies_relation\tspecies_missing_expected_component\tspecies_has_unexpected_component\tspecies_disjoint_wrong_species\tspecies_best_identity\tspecies_best_support\tspecies_gap\tcall_qc_flags\twarnings\n";
    gzwrite(out, header.c_str(), header.size());

    vector<double> species_supports;
    long n_species_conflict = 0;
    long n_species_evidence = 0;
    map<string,long> species_relation_counts;
    set<string> expected_species_sets;
    map<string,long> observed_species_best_counts;

    for (unsigned long ul : order){
        auto it = cells.find(ul);
        if (it == cells.end()) continue;
        CellInfo& c = it->second;
        double conc = concordance(c.assigned_acc, min_evidence);
        double best_runner = NAN;
        string best_runner_name = "NA";
        string best_runner_state = "not_applicable";
        if (!c.runnerups.empty()) {
            best_runner_name = c.runnerups[0].name;
            best_runner_state = !c.runnerup_comparison_states.empty()
                ? c.runnerup_comparison_states[0] : "unavailable";
        }
        for (size_t i=0; i<c.runner_acc.size(); ++i){
            const string state = i < c.runnerup_comparison_states.size()
                ? c.runnerup_comparison_states[i] : "unavailable";
            // Missing/partial comparison states are not usable confidence
            // comparators. Preserve their explicit state for reporting when no
            // complete runner exists, but never turn them into a numeric gap.
            if (!comparison_state_is_present(state)) continue;
            double rc = concordance(c.runner_acc[i], min_evidence);
            if (isfinite(rc) && (!isfinite(best_runner) || rc > best_runner)){
                best_runner = rc;
                best_runner_name = c.runnerups[i].name;
                best_runner_state = state;
            }
        }
        double gap = (isfinite(conc) && isfinite(best_runner)) ? (conc - best_runner) : NAN;
        vector<string> flags;
        vector<string> warnings;
        if (!isfinite(conc)){
            flags.push_back("LOW_EVIDENCE");
            warnings.push_back("LOW_EVIDENCE");
        }
        if (runnerups.empty()) warnings.push_back("NO_RUNNER_UPS_AVAILABLE");
        if (isfinite(gap) && gap < 0) flags.push_back("NEG_GAP");
        if (isfinite(conc) && conc < 0.70) flags.push_back("LOW_CONCORDANCE");

        string species = expected_species_set(c.assignment, samples, panel);
        if (species != "NA") expected_species_sets.insert(species);

        double species_support = NAN;
        string species_conflict = "NA";
        string sp_relation = "NA";
        string sp_missing_expected = "NA";
        string sp_has_unexpected = "NA";
        string sp_disjoint_wrong = "NA";
        string species_best = "NA";
        double species_best_support = NAN;
        double species_gap = NAN;
        if (!species_scoring_enabled){
            if (!species_scoring_requested) warnings.push_back("NO_SPECIES_INPUTS");
            else warnings.push_back(species_disable_reason.empty() ? "SPECIES_SCORING_DISABLED" : species_disable_reason);
        } else if (!c.has_expected_species){
            warnings.push_back("EXPECTED_SPECIES_UNRESOLVED");
        } else {
            species_support = concordance(c.species_expected_acc, min_evidence);
            for (size_t i=0; i<c.species_candidate_acc.size(); ++i){
                double sc = concordance(c.species_candidate_acc[i], min_evidence);
                if (isfinite(sc) && (!isfinite(species_best_support) || sc > species_best_support)){
                    species_best_support = sc;
                    species_best = species_candidates[i].name;
                }
            }
            if (isfinite(species_support)){
                ++n_species_evidence;
                species_supports.push_back(species_support);
                if (isfinite(species_best_support)) species_gap = species_support - species_best_support;
                sp_relation = species_relation(species, species_best);
                bool missing_expected = (sp_relation == "expected_subset_only_component_missing" || sp_relation == "partial_overlap_with_extra_and_missing" || sp_relation == "missing_species_evidence");
                bool has_unexpected = (sp_relation == "expected_superset_with_extra_species" || sp_relation == "partial_overlap_with_extra_and_missing" || sp_relation == "disjoint_wrong_species");
                bool disjoint_wrong = (sp_relation == "disjoint_wrong_species");
                sp_missing_expected = missing_expected ? "1" : "0";
                sp_has_unexpected = has_unexpected ? "1" : "0";
                sp_disjoint_wrong = disjoint_wrong ? "1" : "0";
                bool conflict = has_unexpected || disjoint_wrong || species_support < species_support_threshold;
                species_conflict = conflict ? "1" : "0";
                species_relation_counts[sp_relation]++;
                if (conflict) ++n_species_conflict;
                if (species_best != "NA") observed_species_best_counts[species_best]++;
            } else {
                warnings.push_back("LOW_SPECIES_EVIDENCE");
            }
        }

        string assignment_type = (c.assignment.b < 0 ? "S" : "D");
        string ploidy = "singlet";
        if (c.assignment.b >= 0 && c.assignment.a == c.assignment.b) ploidy = "unresolved_by_SNPs";
        else if (c.assignment.b >= 0) ploidy = "heterotypic";
        double dnllr = (comparison_state_is_present(c.diag.runnerup_comparison_state) &&
                        isfinite(c.diag.llr_vs_runner_up) && c.diag.total_depth > 0)
            ? c.diag.llr_vs_runner_up / c.diag.total_depth : NAN;
        if (!comparison_state_is_present(c.diag.runnerup_comparison_state)){
            warnings.push_back("RUNNER_UP_COMPARISON_" +
                lowercase(c.diag.runnerup_comparison_state));
        }
        string line;
        line += c.barcode + "\t" + libname + "\t" + c.assignment.name + "\t" + assignment_type + "\t" + ploidy + "\t";
        line += fmt(c.diag.llr_vs_runner_up) + "\t" + c.diag.runnerup_comparison_state + "\t" +
                fmt(c.diag.margin_softmax_score) + "\t" + fmt(c.diag.total_depth) + "\t" +
                to_string(c.assigned_acc.bins) + "\t" + fmt(c.assigned_acc.depth) + "\t" +
                to_string(c.diag.n_close) + "\t" + fmt(dnllr) + "\t";
        line += fmt(conc) + "\t" + best_runner_name + "\t" + best_runner_state + "\t" +
                fmt(best_runner) + "\t" + fmt(gap) + "\t" +
                (isfinite(conc) ? fmt(1.0-conc) : "NA") + "\t";
        line += species + "\t" + fmt(species_support) + "\t" + species_conflict + "\t" + sp_relation + "\t" + sp_missing_expected + "\t" + sp_has_unexpected + "\t" + sp_disjoint_wrong + "\t" + species_best + "\t" + fmt(species_best_support) + "\t" + fmt(species_gap) + "\t" + join_flags(flags) + "\t" + join_flags(warnings) + "\n";
        gzwrite(out, line.c_str(), line.size());
    }
    gzclose(out);

    string side = output;
    string suffix = ".call_qc.tsv.gz";
    if (side.size() >= suffix.size() && side.substr(side.size()-suffix.size()) == suffix){
        side = side.substr(0, side.size()-suffix.size()) + ".species_qc.tsv";
        ofstream sp(side.c_str());
        if (sp){
            sp << "library\tn_cells_with_species_evidence\tmedian_species_support_expected\tfrac_cells_species_conflict\tfrac_cells_species_exact_match\tfrac_cells_species_component_missing\tfrac_cells_species_unexpected_extra\tfrac_cells_species_partial_overlap_extra\tfrac_cells_species_disjoint_wrong\tfrac_cells_species_unexpected_or_disjoint\texpected_species_set\tobserved_species_evidence\twarnings\n";
            string expected_join;
            for (auto it=expected_species_sets.begin(); it!=expected_species_sets.end(); ++it){
                if (!expected_join.empty()) expected_join += ";";
                expected_join += *it;
            }
            if (expected_join.empty()) expected_join = "NA";
            string observed;
            long total_best = 0;
            for (auto& kv : observed_species_best_counts) total_best += kv.second;
            for (auto& kv : observed_species_best_counts){
                if (!observed.empty()) observed += ",";
                observed += kv.first + ":" + fmt(total_best ? ((double)kv.second / (double)total_best) : NAN);
            }
            if (observed.empty()) observed = "NA";
            vector<string> side_warn;
            if (!species_scoring_enabled){
                if (!species_scoring_requested) side_warn.push_back("NO_SPECIES_INPUTS");
                else side_warn.push_back(species_disable_reason.empty() ? "SPECIES_SCORING_DISABLED" : species_disable_reason);
            }
            auto frac_rel = [&](const string& key)->string {
                return n_species_evidence ? fmt((double)species_relation_counts[key] / (double)n_species_evidence) : string("NA");
            };
            long unexpected_or_disjoint = species_relation_counts["expected_superset_with_extra_species"] + species_relation_counts["partial_overlap_with_extra_and_missing"] + species_relation_counts["disjoint_wrong_species"];
            sp << libname << "\t" << n_species_evidence << "\t" << fmt(median(species_supports)) << "\t";
            sp << (n_species_evidence ? fmt((double)n_species_conflict / (double)n_species_evidence) : string("NA")) << "\t";
            sp << frac_rel("exact_match") << "\t";
            sp << frac_rel("expected_subset_only_component_missing") << "\t";
            sp << frac_rel("expected_superset_with_extra_species") << "\t";
            sp << frac_rel("partial_overlap_with_extra_and_missing") << "\t";
            sp << frac_rel("disjoint_wrong_species") << "\t";
            sp << (n_species_evidence ? fmt((double)unexpected_or_disjoint / (double)n_species_evidence) : string("NA")) << "\t";
            sp << expected_join << "\t" << observed << "\t" << join_flags(side_warn) << "\n";
        }
    }
    fprintf(stderr, "Wrote %s for %lu cells\n", output.c_str(), (unsigned long)order.size());
    if (species_scoring_enabled){
        fprintf(stderr, "Native species support scoring enabled: n_species=%lu, species_counts=%s, species_samples=%s\n",
                (unsigned long)species_samples.size(), species_counts.c_str(), species_samples_path.c_str());
    } else if (species_scoring_requested){
        fprintf(stderr, "Native species support scoring disabled: %s\n", species_disable_reason.c_str());
    }
    return 0;
}
