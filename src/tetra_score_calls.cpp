#include <algorithm>
#include <cctype>
#include <cerrno>
#include <climits>
#include <clocale>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <limits>
#include <initializer_list>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <queue>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <dirent.h>
#include <float.h>
#include <sys/stat.h>
#include <unistd.h>
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



struct CandidateHypothesis {
    string library;
    string barcode;
    string hypothesis_id;
    string state_notation;
    string donor_genotype;
    string current_donor_genotype;
    string score_pair_id;
    string score_pair_role;
    string score_scope_contract;
    string schema_version;
    string candidate_origin;
    string expected_genotype_status;
    string project_genotype_status;
    string biological_admissibility;
    string score_pair_source;
    string score_population_scope;
    string population_votes_in_authoritative_event;
    string supported_event_key;
    string selected_supported_event_id;
    string selected_supported_event_proposal;
    string reconciliation_event_id;
    string reconciliation_event_class;
    string reconciliation_event_confidence;
    string reconciliation_final_action;
    string reconciliation_decision_confidence;
    string reconciliation_reassignment_applied;
    string reconciliation_current_refined_assignment;
    string reconciled_donor_genotype;
    string reconciled_droplet_state;
    string original_demux_assignment;
    string proposed_donor_genotype;
    string reconciliation_nominated_swap;
    string candidate_b_fixed_identity;
    string source_reconciliation_event_id;
    string source_reconciliation_proposed_identity;
    string pair_construction_mode;
    Identity identity;
    bool scoreable = false;
    ScoreAccum accum;
    double log_likelihood = 0.0;
};

static double binom_log_kernel(double ref, double alt, double expected_alt){
    const double eps = 1e-12;
    double q = expected_alt;
    if (q < eps) q = eps;
    if (q > 1.0 - eps) q = 1.0 - eps;
    return alt * log(q) + ref * log(1.0 - q);
}

static map<string,int> header_index_simple(const vector<string>& header){
    map<string,int> out;
    for (int i=0; i<(int)header.size(); ++i) out[header[i]] = i;
    return out;
}

static int required_simple_col(const map<string,int>& idx, const string& name, const string& path){
    auto it = idx.find(name);
    if (it == idx.end()) die(path + ": required column missing: " + name);
    return it->second;
}

static int optional_simple_col(const map<string,int>& idx, const string& name){
    auto it = idx.find(name);
    return it == idx.end() ? -1 : it->second;
}

static vector<CandidateHypothesis> load_candidate_manifest(
    const string& path,
    const unordered_map<string,int>& sample2idx,
    unordered_map<unsigned long, vector<size_t>>& by_cell,
    bool candidate_axis_mode = false){
    gzFile gz = gzopen(path.c_str(), "rb");
    if (!gz) die("could not open candidate manifest: " + path);
    char buf[1<<20];
    if (!gzgets(gz, buf, sizeof(buf))){ gzclose(gz); die(path + ": empty candidate manifest"); }
    string hline(buf); hline.erase(remove(hline.begin(), hline.end(), '\n'), hline.end()); hline.erase(remove(hline.begin(), hline.end(), '\r'), hline.end());
    vector<string> header = split_tsv_strict(hline);
    map<string,int> idx = header_index_simple(header);
    const int lib_i = required_simple_col(idx, "library", path);
    const int bc_i = required_simple_col(idx, "barcode", path);
    const int hyp_i = required_simple_col(idx, "hypothesis_id", path);
    const int state_i = required_simple_col(idx, "state_notation", path);
    const int donor_i = required_simple_col(idx, "donor_genotype", path);
    const int current_i = required_simple_col(idx, "current_donor_genotype", path);
    const int score_pair_id_i = optional_simple_col(idx, "score_pair_id");
    const int score_pair_role_i = optional_simple_col(idx, "score_pair_role");
    const int score_scope_i = optional_simple_col(idx, "score_scope_contract");
    const int schema_i = optional_simple_col(idx, "schema_version");
    const int origin_i = optional_simple_col(idx, "candidate_origin");
    const int expected_i = optional_simple_col(idx, "expected_genotype_status");
    const int project_i = optional_simple_col(idx, "project_genotype_status");
    const int admissibility_i = optional_simple_col(idx, "biological_admissibility");
    const int pair_source_i = optional_simple_col(idx, "score_pair_source");
    const int population_i = optional_simple_col(idx, "score_population_scope");
    const int votes_i = optional_simple_col(idx, "population_votes_in_authoritative_event");
    const int event_key_i = optional_simple_col(idx, "supported_event_key");
    const int selected_event_i = optional_simple_col(idx, "selected_supported_event_id");
    const int selected_proposal_i = optional_simple_col(idx, "selected_supported_event_proposal");
    const int event_id_i = optional_simple_col(idx, "reconciliation_event_id");
    const int event_class_i = optional_simple_col(idx, "reconciliation_event_class");
    const int event_confidence_i = optional_simple_col(idx, "reconciliation_event_confidence");
    const int final_action_i = optional_simple_col(idx, "reconciliation_final_action");
    const int decision_confidence_i = optional_simple_col(idx, "reconciliation_decision_confidence");
    const int applied_i = optional_simple_col(idx, "reconciliation_reassignment_applied");
    const int current_refined_i = optional_simple_col(idx, "reconciliation_current_refined_assignment");
    const int reconciled_genotype_i = optional_simple_col(idx, "reconciled_donor_genotype");
    const int reconciled_droplet_i = optional_simple_col(idx, "reconciled_droplet_state");
    const int original_i = optional_simple_col(idx, "original_demux_assignment");
    const int proposed_i = optional_simple_col(idx, "proposed_donor_genotype");
    const int nominated_i = optional_simple_col(idx, "reconciliation_nominated_swap");
    const int fixed_b_i = optional_simple_col(idx, "candidate_b_fixed_identity");
    const int source_event_i = optional_simple_col(idx, "source_reconciliation_event_id");
    const int source_proposal_i = optional_simple_col(idx, "source_reconciliation_proposed_identity");
    const int construction_i = optional_simple_col(idx, "pair_construction_mode");
    if (candidate_axis_mode){
        const vector<string> required_axis = {
            "score_pair_id", "score_pair_role", "score_scope_contract",
            "schema_version", "candidate_origin", "score_pair_source",
            "score_population_scope", "population_votes_in_authoritative_event",
            "supported_event_key", "selected_supported_event_id",
            "selected_supported_event_proposal", "original_demux_assignment",
            "candidate_b_fixed_identity", "pair_construction_mode"
        };
        for (const string& name : required_axis) required_simple_col(idx, name, path);
    }
    vector<CandidateHypothesis> out;
    unordered_map<unsigned long,string> encoded_barcodes;
    long line_no = 1;
    while (gzgets(gz, buf, sizeof(buf))){
        ++line_no;
        string line(buf); line.erase(remove(line.begin(), line.end(), '\n'), line.end()); line.erase(remove(line.begin(), line.end(), '\r'), line.end());
        if (line.empty()) continue;
        vector<string> f = split_tsv_strict(line);
        if (f.size() != header.size()) die(path + ": malformed candidate row at line " + to_string(line_no));
        CandidateHypothesis c;
        c.library = f[lib_i]; c.barcode = f[bc_i]; c.hypothesis_id = f[hyp_i]; c.state_notation = f[state_i]; c.donor_genotype = trim(f[donor_i]); c.current_donor_genotype = trim(f[current_i]);
        if (score_pair_id_i >= 0) c.score_pair_id = trim(f[score_pair_id_i]);
        if (score_pair_role_i >= 0) c.score_pair_role = trim(f[score_pair_role_i]);
        if (score_scope_i >= 0) c.score_scope_contract = trim(f[score_scope_i]);
        if (schema_i >= 0) c.schema_version = trim(f[schema_i]);
        if (origin_i >= 0) c.candidate_origin = trim(f[origin_i]);
        if (expected_i >= 0) c.expected_genotype_status = trim(f[expected_i]);
        if (project_i >= 0) c.project_genotype_status = trim(f[project_i]);
        if (admissibility_i >= 0) c.biological_admissibility = trim(f[admissibility_i]);
        if (pair_source_i >= 0) c.score_pair_source = trim(f[pair_source_i]);
        if (population_i >= 0) c.score_population_scope = trim(f[population_i]);
        if (votes_i >= 0) c.population_votes_in_authoritative_event = trim(f[votes_i]);
        if (event_key_i >= 0) c.supported_event_key = trim(f[event_key_i]);
        if (selected_event_i >= 0) c.selected_supported_event_id = trim(f[selected_event_i]);
        if (selected_proposal_i >= 0) c.selected_supported_event_proposal = trim(f[selected_proposal_i]);
        if (event_id_i >= 0) c.reconciliation_event_id = trim(f[event_id_i]);
        if (event_class_i >= 0) c.reconciliation_event_class = trim(f[event_class_i]);
        if (event_confidence_i >= 0) c.reconciliation_event_confidence = trim(f[event_confidence_i]);
        if (final_action_i >= 0) c.reconciliation_final_action = trim(f[final_action_i]);
        if (decision_confidence_i >= 0) c.reconciliation_decision_confidence = trim(f[decision_confidence_i]);
        if (applied_i >= 0) c.reconciliation_reassignment_applied = trim(f[applied_i]);
        if (current_refined_i >= 0) c.reconciliation_current_refined_assignment = trim(f[current_refined_i]);
        if (reconciled_genotype_i >= 0) c.reconciled_donor_genotype = trim(f[reconciled_genotype_i]);
        if (reconciled_droplet_i >= 0) c.reconciled_droplet_state = trim(f[reconciled_droplet_i]);
        if (original_i >= 0) c.original_demux_assignment = trim(f[original_i]);
        if (proposed_i >= 0) c.proposed_donor_genotype = trim(f[proposed_i]);
        if (nominated_i >= 0) c.reconciliation_nominated_swap = trim(f[nominated_i]);
        if (fixed_b_i >= 0) c.candidate_b_fixed_identity = trim(f[fixed_b_i]);
        if (source_event_i >= 0) c.source_reconciliation_event_id = trim(f[source_event_i]);
        if (source_proposal_i >= 0) c.source_reconciliation_proposed_identity = trim(f[source_proposal_i]);
        if (construction_i >= 0) c.pair_construction_mode = trim(f[construction_i]);
        // Candidate scoring intentionally supports only one/two SNP-resolvable donors.
        const size_t first_plus = c.donor_genotype.find('+');
        const bool too_many = first_plus != string::npos && c.donor_genotype.find('+', first_plus + 1) != string::npos;
        if (!c.donor_genotype.empty() && !too_many){
            try { c.identity = parse_identity(c.donor_genotype, sample2idx); c.scoreable = true; }
            catch (...) { c.scoreable = false; }
        }
        const unsigned long encoded = bc_ul(c.barcode);
        if (candidate_axis_mode){
            auto seen = encoded_barcodes.find(encoded);
            if (seen != encoded_barcodes.end() && seen->second != c.barcode)
                die(path + ": barcode encoding collision between '" + seen->second +
                    "' and '" + c.barcode + "'");
            encoded_barcodes[encoded] = c.barcode;
        }
        const size_t index = out.size();
        out.push_back(c);
        by_cell[encoded].push_back(index);
    }
    if (gzclose(gz) != Z_OK) die("failed closing candidate manifest: " + path);
    return out;
}

static void score_candidate_counts(
    const string& path,
    vector<CandidateHypothesis>& candidates,
    const unordered_map<unsigned long, vector<size_t>>& by_cell,
    double e_ref,
    double e_alt){
    gzFile gz = gzopen(path.c_str(), "rb");
    if (!gz) die("could not open counts: " + path);
    char buf[1<<20];
    while (gzgets(gz, buf, sizeof(buf))){
        vector<string> f = split(string(buf), '\t');
        if (f.size() < 7) continue;
        const unsigned long cell = strtoul(f[0].c_str(), nullptr, 10);
        auto hit = by_cell.find(cell); if (hit == by_cell.end()) continue;
        const int indv1 = atoi(f[1].c_str()), type1 = atoi(f[2].c_str());
        const int indv2 = atoi(f[3].c_str()), type2 = atoi(f[4].c_str());
        const double ref = atof(f[5].c_str()), alt = atof(f[6].c_str());
        for (size_t ci : hit->second){
            CandidateHypothesis& c = candidates[ci]; if (!c.scoreable) continue;
            double expected = 0.0;
            if (!expected_for_record(c.identity, indv1, type1, indv2, type2, expected)) continue;
            update_score(c.accum, expected, ref, alt, e_ref, e_alt);
            c.log_likelihood += binom_log_kernel(ref, alt, adjust_expected(expected, e_ref, e_alt));
        }
    }
    if (gzclose(gz) != Z_OK) die("failed closing counts: " + path);
}

static void write_candidate_scores(
    const string& output,
    const vector<CandidateHypothesis>& candidates,
    const unordered_map<unsigned long, vector<size_t>>& by_cell,
    long min_evidence,
    const string& score_prefix){
    vector<int> ranks(candidates.size(), 0);
    vector<double> current_ll(candidates.size(), NAN);
    for (const auto& kv : by_cell){
        vector<size_t> scoreable;
        double cur = NAN;
        for (size_t ci : kv.second){
            const CandidateHypothesis& c = candidates[ci];
            if (!c.scoreable || c.accum.depth < min_evidence) continue;
            scoreable.push_back(ci);
            if (c.state_notation == "" || c.donor_genotype == "") continue;
        }
        sort(scoreable.begin(), scoreable.end(), [&](size_t a, size_t b){ return candidates[a].log_likelihood > candidates[b].log_likelihood; });
        for (size_t r=0; r<scoreable.size(); ++r) ranks[scoreable[r]] = (int)r + 1;
        // Compare every hypothesis to the frozen current donor genotype from the
        // candidate manifest.  Do not infer the reference from row order or IDs.
        size_t cur_idx = (size_t)-1;
        for (size_t ci : kv.second){
            if (candidates[ci].scoreable &&
                candidates[ci].donor_genotype == candidates[ci].current_donor_genotype){
                cur_idx = ci;
                break;
            }
        }
        if (cur_idx == (size_t)-1 && !kv.second.empty()) cur_idx = kv.second.front();
        if (cur_idx != (size_t)-1 && candidates[cur_idx].scoreable && candidates[cur_idx].accum.depth >= min_evidence) cur = candidates[cur_idx].log_likelihood;
        for (size_t ci : kv.second) current_ll[ci] = cur;
    }

    gzFile out = gzopen(output.c_str(), "wb"); if (!out) die("could not open candidate score output: " + output);
    const string pre = score_prefix.empty() ? "" : score_prefix + "_";
    gzprintf(out, "library\tbarcode\thypothesis_id\tstate_notation\tdonor_genotype\t%slog_likelihood\t%sdelta_ll_vs_current\t%srank_within_candidates\t%sinformative_depth\t%sinformative_units\t%sdosage_concordance\t%sresidual_mismatch\t%sdepth_normalized_delta\tcomparison_state\t%sscore_status\tschema_version\n",
        pre.c_str(), pre.c_str(), pre.c_str(), pre.c_str(), pre.c_str(), pre.c_str(), pre.c_str(), pre.c_str(), pre.c_str());
    for (size_t i=0; i<candidates.size(); ++i){
        const CandidateHypothesis& c = candidates[i];
        double conc = concordance(c.accum, min_evidence);
        double residual = isfinite(conc) ? 1.0 - conc : NAN;
        const double cur = current_ll[i];
        const double delta = c.scoreable && isfinite(cur) && c.accum.depth >= min_evidence ? c.log_likelihood - cur : NAN;
        const double norm = isfinite(delta) && c.accum.depth > 0 ? delta / c.accum.depth : NAN;
        string status = !c.scoreable ? "UNSUPPORTED_COMPONENT_COUNT_OR_PANEL_IDENTITY" : (c.accum.depth < min_evidence ? "LOW_EVIDENCE" : "PASS");
        gzprintf(out, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%ld\t%s\t%s\t%s\t%s\t%s\tidentity_hypothesis_scores_v1\n",
            c.library.c_str(), c.barcode.c_str(), c.hypothesis_id.c_str(), c.state_notation.c_str(), c.donor_genotype.c_str(),
            fmt(c.scoreable && c.accum.depth >= min_evidence ? c.log_likelihood : NAN).c_str(), fmt(delta).c_str(), ranks[i], fmt(c.accum.depth).c_str(), c.accum.bins,
            fmt(conc).c_str(), fmt(residual).c_str(), fmt(norm).c_str(), isfinite(delta) ? "present" : "unavailable", status.c_str());
    }
    if (gzclose(out) != Z_OK) die("failed closing candidate score output: " + output);
}

struct FoldSite { vector<int8_t> genotype; };
struct FoldAccum { double ll = 0.0; double depth = 0.0; long units = 0; };

static uint64_t site_key(int tid, int pos){ return (uint64_t)(uint32_t)tid << 32 | (uint32_t)pos; }
static int stable_fold_cpp(int tid, int pos, int folds){
    uint64_t x = site_key(tid,pos); x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
    return folds > 0 ? (int)(x % (uint64_t)folds) : 0;
}

static void score_site_folds(
    const string& sites_path,
    const string& obs_path,
    const string& output,
    vector<CandidateHypothesis>& candidates,
    const unordered_map<unsigned long, vector<size_t>>& by_cell,
    int n_samples,
    int folds,
    double e_ref,
    double e_alt){
    if (output.empty()) return;
    if (folds < 1) folds = 1;
    vector<vector<FoldAccum>> accum(candidates.size(), vector<FoldAccum>(folds));
    if (!file_exists(sites_path) || !file_exists(obs_path)){
        gzFile out = gzopen(output.c_str(), "wb"); if (!out) die("could not open fold output: " + output);
        gzprintf(out, "library\tbarcode\thypothesis_id\tfold\tfold_delta_ll\tfold_informative_depth\tfold_status\tschema_version\n");
        for (const auto& c : candidates) for (int f=0; f<folds; ++f) gzprintf(out, "%s\t%s\t%s\t%d\tNA\t0\tSITE_FOLD_UNAVAILABLE\tidentity_site_fold_scores_v1\n", c.library.c_str(), c.barcode.c_str(), c.hypothesis_id.c_str(), f);
        gzclose(out); return;
    }
    // The site sidecar contains the complete panel (often ~20 million sites).
    // Discover the much smaller set observed in candidate cells before loading
    // genotypes; loading the whole panel scales as panel_sites * roster_size and
    // can exceed memory for large-pool libraries.
    unordered_set<uint64_t> observed_sites;
    observed_sites.reserve(by_cell.size() * 32);
    char buf[1<<20];
    gzFile discovery_gz = gzopen(obs_path.c_str(), "rb");
    if (!discovery_gz) die("could not open pileup observations: " + obs_path);
    while (gzgets(discovery_gz, buf, sizeof(buf))){
        vector<string> f = split(string(buf), '\t');
        if (f.size() < 5) continue;
        const unsigned long bc = strtoul(f[0].c_str(), nullptr, 10);
        if (by_cell.find(bc) == by_cell.end()) continue;
        const double ref = atof(f[3].c_str());
        const double alt = atof(f[4].c_str());
        if (ref + alt <= 0.0) continue;
        observed_sites.insert(site_key(atoi(f[1].c_str()), atoi(f[2].c_str())));
    }
    if (gzclose(discovery_gz) != Z_OK)
        die("failed closing pileup observations: " + obs_path);

    unordered_map<uint64_t, FoldSite> sites;
    sites.reserve(observed_sites.size());
    gzFile sgz = gzopen(sites_path.c_str(), "rb"); if (!sgz) die("could not open pileup sites: " + sites_path);
    while (gzgets(sgz, buf, sizeof(buf))){
        vector<string> f = split(string(buf), '\t'); if ((int)f.size() < 5 + n_samples) continue;
        int tid = atoi(f[0].c_str()), pos = atoi(f[2].c_str());
        const uint64_t key = site_key(tid, pos);
        if (observed_sites.find(key) == observed_sites.end()) continue;
        FoldSite site; site.genotype.resize(n_samples, -1);
        for (int i=0;i<n_samples;++i)
            site.genotype[i]=(int8_t)atoi(f[5+i].c_str());
        sites[key] = site;
    }
    if (gzclose(sgz) != Z_OK) die("failed closing pileup sites: " + sites_path);
    observed_sites.clear();
    observed_sites.rehash(0);
    gzFile ogz = gzopen(obs_path.c_str(), "rb"); if (!ogz) die("could not open pileup observations: " + obs_path);
    while (gzgets(ogz, buf, sizeof(buf))){
        vector<string> f = split(string(buf), '\t'); if (f.size() < 5) continue;
        unsigned long bc = strtoul(f[0].c_str(), nullptr, 10); auto hit = by_cell.find(bc); if (hit == by_cell.end()) continue;
        int tid=atoi(f[1].c_str()), pos=atoi(f[2].c_str()); auto sit=sites.find(site_key(tid,pos)); if (sit==sites.end()) continue;
        double ref=atof(f[3].c_str()), alt=atof(f[4].c_str()); int fold=stable_fold_cpp(tid,pos,folds);
        for (size_t ci : hit->second){
            const CandidateHypothesis& c=candidates[ci]; if(!c.scoreable) continue;
            int ga = c.identity.a >= 0 && c.identity.a < n_samples ? sit->second.genotype[c.identity.a] : -1;
            int gb = c.identity.b >= 0 && c.identity.b < n_samples ? sit->second.genotype[c.identity.b] : -1;
            if (ga < 0 || ga > 2) continue;
            double p = 0.0;
            if (c.identity.b < 0 || c.identity.a == c.identity.b) p = ga / 2.0;
            else { if (gb < 0 || gb > 2) continue; p = (ga + gb) / 4.0; }
            FoldAccum& a=accum[ci][fold]; a.ll += binom_log_kernel(ref,alt,adjust_expected(p,e_ref,e_alt)); a.depth += ref+alt; a.units += 1;
        }
    }
    if (gzclose(ogz) != Z_OK) die("failed closing pileup observations: " + obs_path);
    vector<vector<double>> current(candidates.size(), vector<double>(folds,NAN));
    for (const auto& kv:by_cell){
        size_t cur=(size_t)-1; for(size_t ci:kv.second){if(candidates[ci].donor_genotype==candidates[ci].current_donor_genotype){cur=ci;break;}}
        if(cur==(size_t)-1&&!kv.second.empty()) cur=kv.second.front();
        if(cur!=(size_t)-1) for(int f=0;f<folds;++f) for(size_t ci:kv.second) current[ci][f]=accum[cur][f].depth>0?accum[cur][f].ll:NAN;
    }
    gzFile out=gzopen(output.c_str(),"wb"); if(!out)die("could not open fold output: "+output);
    gzprintf(out,"library\tbarcode\thypothesis_id\tfold\tfold_delta_ll\tfold_informative_depth\tfold_status\tschema_version\n");
    for(size_t i=0;i<candidates.size();++i) for(int f=0;f<folds;++f){
        double delta=(accum[i][f].depth>0&&isfinite(current[i][f]))?accum[i][f].ll-current[i][f]:NAN;
        gzprintf(out,"%s\t%s\t%s\t%d\t%s\t%s\t%s\tidentity_site_fold_scores_v1\n",candidates[i].library.c_str(),candidates[i].barcode.c_str(),candidates[i].hypothesis_id.c_str(),f,fmt(delta).c_str(),fmt(accum[i][f].depth).c_str(),accum[i][f].depth>0?"PASS":"NO_INFORMATIVE_SITES");
    }
    gzclose(out);
}

// -------------------------------------------------------------------------
// Common-evidence post-reconciliation identity probabilities
// -------------------------------------------------------------------------

struct PairSiteEvidence {
    uint64_t key = 0;
    int tid = -1;
    int pos = -1;
    double ref = 0.0;
    double alt = 0.0;
};

struct PairMoleculeEvidence {
    uint64_t molecule = 0;
    uint64_t key = 0;
    string basis;
    double ref = 0.0;
    double alt = 0.0;
};

struct PairSiteDefinition {
    int tid = -1;
    int pos = -1;
    vector<int8_t> genotype;
};

struct PairContribution {
    uint64_t key = 0;
    int tid = -1;
    double delta_a_minus_b = 0.0;
};

struct PairEvaluation {
    string status = "NO_COMMON_EVIDENCE";
    string probability_basis = "nuclear_site_likelihood_equal_priors";
    vector<string> warnings;
    double delta_a_minus_b = NAN;
    double probability_a = NAN;
    double probability_b = NAN;
    double site_delta_a_minus_b = NAN;
    double site_probability_a = NAN;
    double site_probability_b = NAN;
    double molecule_delta_a_minus_b = NAN;
    double molecule_probability_a = NAN;
    double molecule_probability_b = NAN;
    double residual_a = NAN;
    double residual_b = NAN;
    double genotype_similarity = NAN;
    long common_sites = 0;
    long discriminating_sites = 0;
    double common_depth = 0.0;
    double discriminating_depth = 0.0;
    long sites_favor_a = 0;
    long sites_favor_b = 0;
    long sites_neutral = 0;
    int chromosomes_covered = 0;
    double effective_snps = NAN;
    double maximum_site_fraction = NAN;
    double top_five_site_fraction = NAN;
    long independent_molecules = 0;
    double effective_molecules = NAN;
    double maximum_molecule_fraction = NAN;
    double umi_gene_molecule_fraction = NAN;
    string molecule_status = "MOLECULE_SIDECAR_UNAVAILABLE";
    double probability_without_top_site = NAN;
    double probability_without_top_five_sites = NAN;
    double probability_without_top_molecule = NAN;
    double minimum_leave_one_chromosome_out_probability = NAN;
    bool winner_changed_after_influence_removal = false;
    double site_bootstrap_win_fraction = NAN;
    double downsample_50pct_win_fraction = NAN;
    string downsample_basis = "UNAVAILABLE";
    double minimum_error_sensitivity_probability = NAN;
    bool error_sensitivity_stable = false;
    vector<PairContribution> site_contributions;
    unordered_map<uint64_t,double> molecule_contributions;
    unordered_map<uint64_t,long> molecule_site_units;
    unordered_map<uint64_t,string> molecule_basis;
    map<int,double> chromosome_contributions;
};

static uint64_t stable_text_hash(const string& value){
    uint64_t hash = 1469598103934665603ULL;
    for (unsigned char ch : value){
        hash ^= (uint64_t)ch;
        hash *= 1099511628211ULL;
    }
    return hash;
}

static double probability_from_delta(double delta){
    if (!isfinite(delta)) return NAN;
    if (delta >= 0.0){
        return 1.0 / (1.0 + exp(-std::min(delta, 745.0)));
    }
    const double e = exp(std::max(delta, -745.0));
    return e / (1.0 + e);
}

static bool expected_probability_at_site(
    const Identity& identity,
    const PairSiteDefinition& site,
    double& expected){
    if (identity.a < 0 || identity.a >= (int)site.genotype.size()) return false;
    const int ga = site.genotype[identity.a];
    if (ga < 0 || ga > 2) return false;
    if (identity.b < 0 || identity.a == identity.b){
        expected = ga / 2.0;
        return true;
    }
    if (identity.b >= (int)site.genotype.size()) return false;
    const int gb = site.genotype[identity.b];
    if (gb < 0 || gb > 2) return false;
    expected = (ga + gb) / 4.0;
    return true;
}

static void load_pair_observations(
    const string& path,
    const unordered_map<unsigned long, vector<size_t>>& candidate_by_cell,
    unordered_map<unsigned long, vector<PairSiteEvidence>>& by_cell,
    unordered_set<uint64_t>& observed_sites){
    gzFile gz = gzopen(path.c_str(), "rb");
    if (!gz) die("could not open common-evidence pileup observations: " + path);
    by_cell.reserve(candidate_by_cell.size());
    char buf[1<<20];
    while (gzgets(gz, buf, sizeof(buf))){
        vector<string> f = split(string(buf), '\t');
        if (f.size() < 5) continue;
        const unsigned long barcode = strtoul(f[0].c_str(), NULL, 10);
        if (candidate_by_cell.find(barcode) == candidate_by_cell.end()) continue;
        PairSiteEvidence observation;
        observation.tid = atoi(f[1].c_str());
        observation.pos = atoi(f[2].c_str());
        observation.key = site_key(observation.tid, observation.pos);
        observation.ref = atof(f[3].c_str());
        observation.alt = atof(f[4].c_str());
        if (observation.ref + observation.alt <= 0.0) continue;
        by_cell[barcode].push_back(observation);
        observed_sites.insert(observation.key);
    }
    if (gzclose(gz) != Z_OK) die("failed closing pileup observations: " + path);
    for (auto& cell : by_cell){
        vector<PairSiteEvidence>& target = cell.second;
        sort(target.begin(), target.end(),
             [](const PairSiteEvidence& a, const PairSiteEvidence& b){
                 return a.key < b.key;
             });
        size_t write = 0;
        for (size_t read = 0; read < target.size(); ++read){
            if (write > 0 && target[write - 1].key == target[read].key){
                target[write - 1].ref += target[read].ref;
                target[write - 1].alt += target[read].alt;
            } else {
                if (write != read) target[write] = target[read];
                ++write;
            }
        }
        target.resize(write);
    }
}

static bool load_pair_molecules(
    const string& path,
    const unordered_map<unsigned long, vector<size_t>>& candidate_by_cell,
    unordered_map<unsigned long, vector<PairMoleculeEvidence>>& by_cell,
    unordered_set<uint64_t>& observed_sites){
    if (path.empty() || !file_exists(path)) return false;
    gzFile gz = gzopen(path.c_str(), "rb");
    if (!gz) die("could not open molecule-aware pileup observations: " + path);
    by_cell.reserve(candidate_by_cell.size());
    char buf[1<<20];
    while (gzgets(gz, buf, sizeof(buf))){
        vector<string> f = split(string(buf), '\t');
        if (f.size() < 7) continue;
        const unsigned long barcode = strtoul(f[0].c_str(), NULL, 10);
        if (candidate_by_cell.find(barcode) == candidate_by_cell.end()) continue;
        PairMoleculeEvidence observation;
        observation.molecule = strtoull(f[1].c_str(), NULL, 10);
        observation.basis = trim(f[2]);
        const int tid = atoi(f[3].c_str());
        const int pos = atoi(f[4].c_str());
        observation.key = site_key(tid, pos);
        observation.ref = atof(f[5].c_str());
        observation.alt = atof(f[6].c_str());
        if (observation.ref + observation.alt <= 0.0) continue;
        by_cell[barcode].push_back(observation);
        observed_sites.insert(observation.key);
    }
    if (gzclose(gz) != Z_OK) die("failed closing molecule pileup: " + path);
    for (auto& cell : by_cell){
        vector<PairMoleculeEvidence>& target = cell.second;
        sort(target.begin(), target.end(),
             [](const PairMoleculeEvidence& a, const PairMoleculeEvidence& b){
                 if (a.molecule != b.molecule) return a.molecule < b.molecule;
                 return a.key < b.key;
             });
        size_t write = 0;
        for (size_t read = 0; read < target.size(); ++read){
            if (write > 0 &&
                    target[write - 1].molecule == target[read].molecule &&
                    target[write - 1].key == target[read].key){
                target[write - 1].ref += target[read].ref;
                target[write - 1].alt += target[read].alt;
                if (target[write - 1].basis == "QNAME_FALLBACK" &&
                        target[read].basis != "QNAME_FALLBACK"){
                    target[write - 1].basis = target[read].basis;
                }
            } else {
                if (write != read) target[write] = target[read];
                ++write;
            }
        }
        target.resize(write);
    }
    return true;
}

static unordered_map<uint64_t, PairSiteDefinition> load_pair_sites(
    const string& path,
    const unordered_set<uint64_t>& observed_sites,
    int n_samples){
    gzFile gz = gzopen(path.c_str(), "rb");
    if (!gz) die("could not open common-evidence pileup sites: " + path);
    unordered_map<uint64_t, PairSiteDefinition> sites;
    sites.reserve(observed_sites.size());
    char buf[1<<20];
    while (gzgets(gz, buf, sizeof(buf))){
        vector<string> f = split(string(buf), '\t');
        if ((int)f.size() < 5 + n_samples) continue;
        const int tid = atoi(f[0].c_str());
        const int pos = atoi(f[2].c_str());
        const uint64_t key = site_key(tid, pos);
        if (observed_sites.find(key) == observed_sites.end()) continue;
        PairSiteDefinition site;
        site.tid = tid;
        site.pos = pos;
        site.genotype.resize(n_samples, -1);
        for (int i = 0; i < n_samples; ++i){
            site.genotype[i] = (int8_t)atoi(f[5+i].c_str());
        }
        sites[key] = site;
    }
    if (gzclose(gz) != Z_OK) die("failed closing pileup sites: " + path);
    return sites;
}

static double effective_count_from_contributions(const vector<double>& values){
    double total = 0.0;
    for (double value : values) total += fabs(value);
    if (total <= 0.0) return NAN;
    double squares = 0.0;
    for (double value : values){
        const double weight = fabs(value) / total;
        squares += weight * weight;
    }
    return squares > 0.0 ? 1.0 / squares : NAN;
}

static double preferred_probability(double delta_a_minus_b, int preferred_sign){
    return 100.0 * probability_from_delta(preferred_sign * delta_a_minus_b);
}

static double pair_delta_only(
    const CandidateHypothesis& a,
    const CandidateHypothesis& b,
    const vector<PairSiteEvidence>& observations,
    const unordered_map<uint64_t, PairSiteDefinition>& sites,
    double e_ref,
    double e_alt,
    long* discriminating_sites = NULL){
    double delta = 0.0;
    long informative = 0;
    for (const PairSiteEvidence& observation : observations){
        auto found = sites.find(observation.key);
        if (found == sites.end()) continue;
        double pa = 0.0, pb = 0.0;
        if (!expected_probability_at_site(a.identity, found->second, pa) ||
                !expected_probability_at_site(b.identity, found->second, pb) ||
                fabs(pa - pb) <= 1e-12) continue;
        const double aa = adjust_expected(pa, e_ref, e_alt);
        const double ab = adjust_expected(pb, e_ref, e_alt);
        delta += binom_log_kernel(observation.ref, observation.alt, aa) -
            binom_log_kernel(observation.ref, observation.alt, ab);
        ++informative;
    }
    if (discriminating_sites) *discriminating_sites = informative;
    return informative > 0 ? delta : NAN;
}

static PairEvaluation evaluate_pair(
    const CandidateHypothesis& a,
    const CandidateHypothesis& b,
    const vector<PairSiteEvidence>& observations,
    const vector<PairMoleculeEvidence>& molecule_observations,
    const unordered_map<uint64_t, PairSiteDefinition>& sites,
    double e_ref,
    double e_alt,
    long min_evidence,
    int resamples,
    uint64_t seed,
    double poor_fit_residual){
    PairEvaluation result;
    double ll_a = 0.0, ll_b = 0.0;
    double residual_a_sum = 0.0, residual_b_sum = 0.0;
    double similarity_sum = 0.0;
    set<int> chromosomes;
    for (const PairSiteEvidence& observation : observations){
        auto found = sites.find(observation.key);
        if (found == sites.end()) continue;
        double pa = 0.0, pb = 0.0;
        if (!expected_probability_at_site(a.identity, found->second, pa) ||
                !expected_probability_at_site(b.identity, found->second, pb)) continue;
        const double depth = observation.ref + observation.alt;
        if (depth <= 0.0) continue;
        const double observed = observation.alt / depth;
        const double adjusted_a = adjust_expected(pa, e_ref, e_alt);
        const double adjusted_b = adjust_expected(pb, e_ref, e_alt);
        ++result.common_sites;
        result.common_depth += depth;
        similarity_sum += depth * (1.0 - fabs(pa - pb));
        if (fabs(pa - pb) <= 1e-12) continue;
        residual_a_sum += depth * fabs(observed - adjusted_a);
        residual_b_sum += depth * fabs(observed - adjusted_b);
        const double site_ll_a = binom_log_kernel(
            observation.ref, observation.alt, adjusted_a);
        const double site_ll_b = binom_log_kernel(
            observation.ref, observation.alt, adjusted_b);
        const double delta = site_ll_a - site_ll_b;
        ll_a += site_ll_a;
        ll_b += site_ll_b;
        ++result.discriminating_sites;
        result.discriminating_depth += depth;
        if (delta > 1e-12) ++result.sites_favor_a;
        else if (delta < -1e-12) ++result.sites_favor_b;
        else ++result.sites_neutral;
        chromosomes.insert(found->second.tid);
        PairContribution contribution;
        contribution.key = observation.key;
        contribution.tid = found->second.tid;
        contribution.delta_a_minus_b = delta;
        result.site_contributions.push_back(contribution);
        result.chromosome_contributions[found->second.tid] += delta;
    }
    if (result.discriminating_depth > 0.0){
        result.residual_a = residual_a_sum / result.discriminating_depth;
        result.residual_b = residual_b_sum / result.discriminating_depth;
    }
    if (result.common_depth > 0.0){
        result.genotype_similarity = 100.0 * similarity_sum / result.common_depth;
    }
    result.chromosomes_covered = (int)chromosomes.size();
    if (result.discriminating_sites == 0){
        result.status = result.common_sites > 0
            ? "PANEL_NONDISCRIMINATING" : "NO_COMMON_EVIDENCE";
        return result;
    }
    result.site_delta_a_minus_b = ll_a - ll_b;
    result.site_probability_a = 100.0 * probability_from_delta(
        result.site_delta_a_minus_b);
    result.site_probability_b = 100.0 - result.site_probability_a;
    result.delta_a_minus_b = result.site_delta_a_minus_b;
    result.probability_a = result.site_probability_a;
    result.probability_b = result.site_probability_b;

    vector<double> site_values;
    double site_abs_total = 0.0;
    for (const PairContribution& contribution : result.site_contributions){
        site_values.push_back(contribution.delta_a_minus_b);
        site_abs_total += fabs(contribution.delta_a_minus_b);
    }
    result.effective_snps = effective_count_from_contributions(site_values);
    vector<double> sorted_abs;
    for (double value : site_values) sorted_abs.push_back(fabs(value));
    sort(sorted_abs.begin(), sorted_abs.end(), greater<double>());
    if (site_abs_total > 0.0){
        result.maximum_site_fraction = sorted_abs.front() / site_abs_total;
        result.top_five_site_fraction = accumulate(
            sorted_abs.begin(), sorted_abs.begin() + min<size_t>(5, sorted_abs.size()),
            0.0) / site_abs_total;
    }

    // Molecule rows use the same alleles and genotype expectations as the site
    // comparison, but group all covered SNPs by corrected UMI+gene (or the
    // explicit QNAME fallback).  When molecule evidence is available it is
    // the primary probability basis, preventing long reads and PCR depth from
    // manufacturing independent support. Site likelihood remains explicit as
    // a separate comparison and is the fallback when the sidecar is absent.
    for (const PairMoleculeEvidence& observation : molecule_observations){
        auto found = sites.find(observation.key);
        if (found == sites.end()) continue;
        double pa = 0.0, pb = 0.0;
        if (!expected_probability_at_site(a.identity, found->second, pa) ||
                !expected_probability_at_site(b.identity, found->second, pb) ||
                fabs(pa - pb) <= 1e-12) continue;
        const double molecule_depth = observation.ref + observation.alt;
        if (molecule_depth <= 0.0) continue;
        // Collapse all reads for one molecule/site to one fractional allele
        // observation.  Averaging again across covered sites below gives each
        // corrected UMI+gene (or explicit QNAME fallback) one total evidence
        // unit, so long reads and PCR depth cannot manufacture independence.
        const double delta = binom_log_kernel(
            observation.ref / molecule_depth,
            observation.alt / molecule_depth,
            adjust_expected(pa, e_ref, e_alt)) - binom_log_kernel(
            observation.ref / molecule_depth,
            observation.alt / molecule_depth,
            adjust_expected(pb, e_ref, e_alt));
        result.molecule_contributions[observation.molecule] += delta;
        ++result.molecule_site_units[observation.molecule];
        auto prior = result.molecule_basis.find(observation.molecule);
        if (prior == result.molecule_basis.end() ||
                (prior->second == "QNAME_FALLBACK" &&
                 observation.basis != "QNAME_FALLBACK")){
            result.molecule_basis[observation.molecule] = observation.basis;
        }
    }
    if (!molecule_observations.empty()){
        result.molecule_status = result.molecule_contributions.empty()
            ? "MOLECULE_SIDECAR_NONDISCRIMINATING"
            : "MOLECULE_AWARE";
    }
    if (!result.molecule_contributions.empty()){
        result.molecule_delta_a_minus_b = 0.0;
        for (auto& item : result.molecule_contributions){
            const long units = result.molecule_site_units[item.first];
            if (units > 0) item.second /= (double)units;
            result.molecule_delta_a_minus_b += item.second;
        }
        result.molecule_probability_a = 100.0 * probability_from_delta(
            result.molecule_delta_a_minus_b);
        result.molecule_probability_b = 100.0 -
            result.molecule_probability_a;
        result.delta_a_minus_b = result.molecule_delta_a_minus_b;
        result.probability_a = result.molecule_probability_a;
        result.probability_b = result.molecule_probability_b;
        result.probability_basis =
            "nuclear_molecule_balanced_likelihood_equal_priors";
        result.independent_molecules = (long)result.molecule_contributions.size();
        vector<double> molecule_values;
        double molecule_abs_total = 0.0;
        long umi_gene = 0;
        for (const auto& item : result.molecule_contributions){
            molecule_values.push_back(item.second);
            molecule_abs_total += fabs(item.second);
            const string basis = result.molecule_basis[item.first];
            if (basis == "UB_GX" || basis == "UB_GN") ++umi_gene;
        }
        result.effective_molecules = effective_count_from_contributions(
            molecule_values);
        if (molecule_abs_total > 0.0){
            double maximum = 0.0;
            for (double value : molecule_values)
                maximum = max(maximum, fabs(value));
            result.maximum_molecule_fraction = maximum / molecule_abs_total;
        }
        result.umi_gene_molecule_fraction =
            (double)umi_gene / (double)result.independent_molecules;
    }

    const int preferred_sign = result.delta_a_minus_b >= 0.0 ? 1 : -1;
    if (!result.site_contributions.empty()){
        vector<PairContribution> influence = result.site_contributions;
        sort(influence.begin(), influence.end(), [](const PairContribution& x,
                                                     const PairContribution& y){
            return fabs(x.delta_a_minus_b) > fabs(y.delta_a_minus_b);
        });
        double removed = influence.front().delta_a_minus_b;
        result.probability_without_top_site = preferred_probability(
            result.site_delta_a_minus_b - removed, preferred_sign);
        removed = 0.0;
        for (size_t i = 0; i < min<size_t>(5, influence.size()); ++i)
            removed += influence[i].delta_a_minus_b;
        result.probability_without_top_five_sites = preferred_probability(
            result.site_delta_a_minus_b - removed, preferred_sign);
    }
    if (!result.molecule_contributions.empty()){
        double top_value = 0.0;
        for (const auto& item : result.molecule_contributions){
            if (fabs(item.second) > fabs(top_value)) top_value = item.second;
        }
        result.probability_without_top_molecule = preferred_probability(
            result.delta_a_minus_b - top_value, preferred_sign);
    }
    result.minimum_leave_one_chromosome_out_probability = 100.0;
    for (const auto& item : result.chromosome_contributions){
        result.minimum_leave_one_chromosome_out_probability = min(
            result.minimum_leave_one_chromosome_out_probability,
            preferred_probability(result.site_delta_a_minus_b - item.second,
                                  preferred_sign));
    }
    result.winner_changed_after_influence_removal =
        (isfinite(result.probability_without_top_site) &&
         result.probability_without_top_site <= 50.0) ||
        (isfinite(result.probability_without_top_five_sites) &&
         result.probability_without_top_five_sites <= 50.0) ||
        (isfinite(result.probability_without_top_molecule) &&
         result.probability_without_top_molecule <= 50.0) ||
        (isfinite(result.minimum_leave_one_chromosome_out_probability) &&
         result.minimum_leave_one_chromosome_out_probability <= 50.0);

    if (resamples > 0 && !site_values.empty()){
        mt19937_64 rng(seed);
        uniform_int_distribution<size_t> choose_site(0, site_values.size()-1);
        long wins = 0;
        for (int replicate = 0; replicate < resamples; ++replicate){
            double delta = 0.0;
            for (size_t i = 0; i < site_values.size(); ++i)
                delta += site_values[choose_site(rng)];
            if (preferred_sign * delta > 0.0) ++wins;
        }
        result.site_bootstrap_win_fraction =
            (double)wins / (double)resamples;

        vector<double> downsample_values = site_values;
        if (!result.molecule_contributions.empty()){
            downsample_values.clear();
            for (const auto& item : result.molecule_contributions)
                downsample_values.push_back(item.second);
            result.downsample_basis = "INDEPENDENT_MOLECULES";
        } else {
            result.downsample_basis = "SITE_PROXY_NO_MOLECULE_SIDECAR";
        }
        bernoulli_distribution keep(0.5);
        wins = 0;
        long evaluable = 0;
        for (int replicate = 0; replicate < resamples; ++replicate){
            double delta = 0.0;
            long kept = 0;
            for (double value : downsample_values){
                if (keep(rng)){ delta += value; ++kept; }
            }
            if (kept == 0) continue;
            ++evaluable;
            if (preferred_sign * delta > 0.0) ++wins;
        }
        result.downsample_50pct_win_fraction = evaluable > 0
            ? (double)wins / (double)evaluable : NAN;
    }

    // A small explicit sequencing-error grid detects calls that exist only at
    // one error assumption.  These are data-model probabilities, not an
    // empirical calibration transform.
    vector<double> error_grid = {e_ref, 0.005, 0.01, 0.02};
    result.minimum_error_sensitivity_probability = 100.0;
    for (double error : error_grid){
        long n_sites = 0;
        const double delta = pair_delta_only(
            a, b, observations, sites, error, error, &n_sites);
        if (n_sites > 0 && isfinite(delta)){
            result.minimum_error_sensitivity_probability = min(
                result.minimum_error_sensitivity_probability,
                preferred_probability(delta, preferred_sign));
        }
    }
    result.error_sensitivity_stable =
        result.minimum_error_sensitivity_probability > 50.0;

    const bool no_fit = isfinite(result.residual_a) &&
        isfinite(result.residual_b) && result.residual_a > poor_fit_residual &&
        result.residual_b > poor_fit_residual;
    if (no_fit){
        result.status = "NO_CANDIDATE_FITS";
    } else if (result.discriminating_depth < min_evidence ||
               result.discriminating_sites < 2){
        result.status = "LOW_EVIDENCE";
    } else {
        result.status = "PASS";
    }
    if (result.independent_molecules == 0)
        result.warnings.push_back("MOLECULE_INDEPENDENCE_UNAVAILABLE");
    if (isfinite(result.molecule_delta_a_minus_b) &&
            result.molecule_delta_a_minus_b * result.site_delta_a_minus_b < 0.0)
        result.warnings.push_back("SITE_MOLECULE_WINNER_DISAGREEMENT");
    if (result.independent_molecules > 0 &&
            isfinite(result.umi_gene_molecule_fraction) &&
            result.umi_gene_molecule_fraction < 0.50)
        result.warnings.push_back("QNAME_FALLBACK_DOMINANT");
    if (result.maximum_site_fraction > 0.50)
        result.warnings.push_back("SITE_DOMINATED");
    if (result.discriminating_sites >= 10 &&
            result.top_five_site_fraction > 0.80)
        result.warnings.push_back("TOP_SITES_DOMINATED");
    if (isfinite(result.maximum_molecule_fraction) &&
            result.maximum_molecule_fraction > 0.50)
        result.warnings.push_back("MOLECULE_DOMINATED");
    if (result.winner_changed_after_influence_removal)
        result.warnings.push_back("INFLUENCE_REMOVAL_UNSTABLE");
    if (isfinite(result.downsample_50pct_win_fraction) &&
            result.downsample_50pct_win_fraction < 0.80)
        result.warnings.push_back("DOWNSAMPLE_UNSTABLE");
    if (!result.error_sensitivity_stable)
        result.warnings.push_back("ERROR_MODEL_SENSITIVE");
    return result;
}

static string bool_text(bool value){ return value ? "TRUE" : "FALSE"; }

static void gzwrite_tsv_row(gzFile out, const vector<string>& fields){
    string line;
    for (size_t i = 0; i < fields.size(); ++i){
        if (i) line.push_back('\t');
        line += fields[i].empty() ? "NA" : fields[i];
    }
    line.push_back('\n');
    if (gzwrite(out, line.data(), (unsigned int)line.size()) == 0)
        die("failed writing common-evidence probability output");
}

static string shared_components(const string& a, const string& b){
    multiset<string> left;
    for (const string& item : split(a, '+')) if (!item.empty()) left.insert(item);
    vector<string> shared;
    for (const string& item : split(b, '+')){
        auto found = left.find(item);
        if (found != left.end()){
            shared.push_back(item);
            left.erase(found);
        }
    }
    string result;
    for (const string& item : shared){
        if (!result.empty()) result += ",";
        result += item;
    }
    return result.empty() ? "NONE" : result;
}

struct TargetedReconciliationPair {
    size_t original = (size_t)-1;
    size_t proposed = (size_t)-1;
    string pair_id;
};

static bool legacy_pair_manifest_schema_allowed(const string& schema){
    return schema == "identity_reconciliation_score_pair_manifest_v1" ||
        schema == "identity_reconciliation_score_pair_manifest_v2";
}

static TargetedReconciliationPair targeted_reconciliation_pair(
    const vector<CandidateHypothesis>& candidates,
    const vector<size_t>& indices,
    unsigned long barcode){
    if (indices.size() != 2){
        die("reconciliation probability manifest must contain exactly two "
            "hypotheses per barcode; barcode=" + to_string(barcode) +
            " rows=" + to_string(indices.size()));
    }
    TargetedReconciliationPair result;
    for (size_t index : indices){
        const CandidateHypothesis& candidate = candidates[index];
        if (candidate.score_pair_id.empty())
            die("reconciliation probability manifest is missing score_pair_id "
                "for barcode=" + candidate.barcode);
        if (result.pair_id.empty()) result.pair_id = candidate.score_pair_id;
        if (candidate.score_pair_id != result.pair_id)
            die("reconciliation probability manifest has inconsistent "
                "score_pair_id values for barcode=" + candidate.barcode);
        if (candidate.score_pair_role == "ORIGINAL_ALLOWED_DEMUX"){
            if (result.original != (size_t)-1)
                die("duplicate ORIGINAL_ALLOWED_DEMUX role for barcode=" +
                    candidate.barcode);
            result.original = index;
        } else if (candidate.score_pair_role ==
                   "RECONCILIATION_NOMINATED_SWAP"){
            if (result.proposed != (size_t)-1)
                die("duplicate RECONCILIATION_NOMINATED_SWAP role for barcode=" +
                    candidate.barcode);
            result.proposed = index;
        } else {
            die("probability scoring accepts only post-reconciliation pair roles; "
                "barcode=" + candidate.barcode + " role=" +
                (candidate.score_pair_role.empty() ? "MISSING" :
                 candidate.score_pair_role));
        }
    }
    if (result.original == (size_t)-1 || result.proposed == (size_t)-1)
        die("reconciliation probability manifest must contain one original and "
            "one nominated-swap hypothesis per barcode=" + to_string(barcode));
    const CandidateHypothesis& original = candidates[result.original];
    const CandidateHypothesis& proposed = candidates[result.proposed];
    const string required_contract =
        "ORIGINAL_ALLOWED_DEMUX_VS_RECONCILIATION_NOMINATED_SWAP_ONLY";
    for (const CandidateHypothesis* candidate : {&original, &proposed}){
        if (candidate->score_scope_contract != required_contract)
            die("reconciliation probability pair has an invalid score-scope "
                "contract for barcode=" + original.barcode);
        if (!legacy_pair_manifest_schema_allowed(candidate->schema_version))
            die("reconciliation probability pair has an invalid manifest "
                "schema for barcode=" + original.barcode);
        if (candidate->biological_admissibility !=
                "SINGLET_IDENTITY_CANDIDATE" &&
                candidate->biological_admissibility !=
                "BIOLOGICAL_SINGLE_CELL_ALLOWED")
            die("reconciliation probability pair contains a non-biological "
                "identity for barcode=" + original.barcode);
    }
    if (original.candidate_origin != "ORIGINAL_ALLOWED_DEMUX" ||
            proposed.candidate_origin != "RECONCILIATION_NOMINATED_SWAP")
        die("reconciliation probability pair has invalid candidate provenance "
            "for barcode=" + original.barcode);
    if (original.expected_genotype_status != "EXPECTED" ||
            proposed.expected_genotype_status == "EXPECTED")
        die("reconciliation probability pair violates the original-allowed "
            "versus library-unexpected contract for barcode=" +
            original.barcode);
    const bool original_project_valid =
        original.project_genotype_status ==
            "GLOBAL_REAL_DONOR_LIBRARY_EXPECTED" ||
        original.project_genotype_status ==
            "GLOBAL_REAL_LINE_LIBRARY_EXPECTED";
    const bool proposed_project_valid =
        proposed.project_genotype_status ==
            "GLOBAL_REAL_DONOR_LIBRARY_UNEXPECTED" ||
        proposed.project_genotype_status ==
            "GLOBAL_REAL_LINE_LIBRARY_UNEXPECTED";
    if (!original_project_valid || !proposed_project_valid)
        die("reconciliation probability pair contains a project identity that "
            "is not an allowed original or real unexpected biological line; "
            "barcode=" + original.barcode);
    if (!original.scoreable || !proposed.scoreable)
        die("reconciliation probability pair contains an identity absent from "
            "the nuclear sample roster for barcode=" + original.barcode);
    if (original.donor_genotype != original.current_donor_genotype ||
            proposed.current_donor_genotype != original.donor_genotype)
        die("reconciliation probability pair does not freeze candidate A as "
            "the original allowed demux assignment for barcode=" +
            original.barcode);
    if (original.donor_genotype == proposed.donor_genotype)
        die("reconciliation probability pair contains identical original and "
            "proposed identities for barcode=" + original.barcode);
    return result;
}

static void write_pairwise_probability_scores(
    const string& sites_path,
    const string& observations_path,
    const string& molecules_path,
    const string& output,
    const vector<CandidateHypothesis>& candidates,
    const unordered_map<unsigned long, vector<size_t>>& candidate_by_cell,
    int n_samples,
    double e_ref,
    double e_alt,
    long min_evidence,
    int resamples,
    uint64_t random_seed,
    double poor_fit_residual){
    if (output.empty()) return;
    if (!file_exists(sites_path) || !file_exists(observations_path))
        die("common-evidence probability scoring requires pileup sites and observations");
    unordered_map<unsigned long, TargetedReconciliationPair> targeted_pairs;
    targeted_pairs.reserve(candidate_by_cell.size());
    for (const auto& cell_entry : candidate_by_cell){
        targeted_pairs.emplace(
            cell_entry.first,
            targeted_reconciliation_pair(
                candidates, cell_entry.second, cell_entry.first));
    }
    unordered_map<unsigned long, vector<PairSiteEvidence>> observations;
    unordered_map<unsigned long, vector<PairMoleculeEvidence>> molecules;
    unordered_set<uint64_t> observed_sites;
    load_pair_observations(
        observations_path, candidate_by_cell, observations, observed_sites);
    const bool molecule_sidecar = load_pair_molecules(
        molecules_path, candidate_by_cell, molecules, observed_sites);
    const unordered_map<uint64_t, PairSiteDefinition> sites = load_pair_sites(
        sites_path, observed_sites, n_samples);
    observed_sites.clear();
    observed_sites.rehash(0);

    gzFile out = gzopen(output.c_str(), "wb");
    if (!out) die("could not open probability output: " + output);
    const vector<string> header = {
        "library","barcode","score_pair_id","comparison","comparison_status",
        "candidate_a","candidate_b","candidate_a_role","candidate_b_role",
        "preferred_assignment",
        "preferred_probability_pct","alternative_assignment",
        "alternative_probability_pct","probability_gap_pp",
        "alternative_closeness_pct","candidate_a_probability_pct",
        "candidate_b_probability_pct","delta_log_likelihood_a_minus_b",
        "site_delta_log_likelihood_a_minus_b",
        "site_candidate_a_probability_pct","site_candidate_b_probability_pct",
        "molecule_balanced_delta_log_likelihood_a_minus_b",
        "molecule_balanced_candidate_a_probability_pct",
        "molecule_balanced_candidate_b_probability_pct",
        "probability_basis","evidence_basis","shared_donor_components",
        "n_common_observed_snps","n_discriminating_snps",
        "n_nondiscriminating_snps","common_evidence_depth",
        "discriminating_evidence_depth","n_snps_favor_preferred",
        "n_snps_favor_alternative","n_snps_neutral",
        "genotype_model_similarity_pct","chromosomes_covered",
        "effective_independent_snps","maximum_site_contribution_fraction",
        "top_five_site_contribution_fraction","n_independent_molecules",
        "effective_independent_molecules","maximum_molecule_contribution_fraction",
        "molecule_umi_gene_fraction","molecule_evidence_status",
        "preferred_residual_mismatch","alternative_residual_mismatch",
        "absolute_fit_status","probability_without_top_site_pct",
        "probability_without_top_five_sites_pct",
        "probability_without_top_molecule_pct",
        "minimum_leave_one_chromosome_out_probability_pct",
        "winner_changed_after_influence_removal","site_bootstrap_win_fraction",
        "downsample_50pct_win_fraction","downsample_basis","resamples",
        "candidate_scope_size","candidate_scope_complete","candidate_cap_applied",
        "candidate_scope_stable","strongest_expansion_challenger",
        "preferred_probability_vs_strongest_expansion_challenger_pct",
        "minimum_error_sensitivity_probability_pct","error_sensitivity_stable",
        "ambient_sensitivity_status","error_ref","error_alt","warnings",
        "schema_version"
    };
    gzwrite_tsv_row(out, header);

    vector<unsigned long> barcode_order;
    barcode_order.reserve(candidate_by_cell.size());
    for (const auto& cell_entry : candidate_by_cell)
        barcode_order.push_back(cell_entry.first);
    sort(barcode_order.begin(), barcode_order.end());
    for (unsigned long barcode : barcode_order){
        const TargetedReconciliationPair& pair = targeted_pairs.at(barcode);
        const CandidateHypothesis& candidate_a = candidates[pair.original];
        const CandidateHypothesis& candidate_b = candidates[pair.proposed];
        const vector<PairSiteEvidence> empty_observations;
        const vector<PairMoleculeEvidence> empty_molecules;
        auto obs_it = observations.find(barcode);
        auto mol_it = molecules.find(barcode);
        const vector<PairSiteEvidence>& cell_observations = obs_it == observations.end()
            ? empty_observations : obs_it->second;
        const vector<PairMoleculeEvidence>& cell_molecules = mol_it == molecules.end()
            ? empty_molecules : mol_it->second;
        const uint64_t seed = random_seed ^ barcode ^
            stable_text_hash(pair.pair_id);
        PairEvaluation evaluation = evaluate_pair(
            candidate_a, candidate_b, cell_observations, cell_molecules,
            sites, e_ref, e_alt, min_evidence, resamples, seed,
            poor_fit_residual);
        const bool probability_available =
            isfinite(evaluation.delta_a_minus_b) &&
            isfinite(evaluation.probability_a) &&
            isfinite(evaluation.probability_b);
        const bool a_preferred = !probability_available ||
            evaluation.delta_a_minus_b >= 0.0;
        const CandidateHypothesis& preferred =
            a_preferred ? candidate_a : candidate_b;
        const CandidateHypothesis& alternative =
            a_preferred ? candidate_b : candidate_a;
        const double preferred_probability_value = probability_available
            ? (a_preferred ? evaluation.probability_a : evaluation.probability_b)
            : NAN;
        const double alternative_probability_value = probability_available
            ? (a_preferred ? evaluation.probability_b : evaluation.probability_a)
            : NAN;
        const double preferred_residual = a_preferred
            ? evaluation.residual_a : evaluation.residual_b;
        const double alternative_residual = a_preferred
            ? evaluation.residual_b : evaluation.residual_a;
        const long preferred_sites = a_preferred
            ? evaluation.sites_favor_a : evaluation.sites_favor_b;
        const long alternative_sites = a_preferred
            ? evaluation.sites_favor_b : evaluation.sites_favor_a;
        const string fit_status = !probability_available ? "UNAVAILABLE" :
            (evaluation.status == "NO_CANDIDATE_FITS" ?
             "NO_CANDIDATE_FITS" :
             (isfinite(preferred_residual) &&
              preferred_residual <= poor_fit_residual ? "PASS" : "POOR_FIT"));
        const string warnings = join_flags(evaluation.warnings);
        gzwrite_tsv_row(out, {
            candidate_a.library, candidate_a.barcode, pair.pair_id,
            "reconciliation_swap", evaluation.status,
            candidate_a.donor_genotype, candidate_b.donor_genotype,
            candidate_a.score_pair_role, candidate_b.score_pair_role,
            probability_available ? preferred.donor_genotype : "NA",
            fmt(preferred_probability_value),
            probability_available ? alternative.donor_genotype : "NA",
            fmt(alternative_probability_value),
            fmt(preferred_probability_value - alternative_probability_value),
            fmt(2.0 * alternative_probability_value),
            fmt(evaluation.probability_a), fmt(evaluation.probability_b),
            fmt(evaluation.delta_a_minus_b),
            fmt(evaluation.site_delta_a_minus_b),
            fmt(evaluation.site_probability_a),
            fmt(evaluation.site_probability_b),
            fmt(evaluation.molecule_delta_a_minus_b),
            fmt(evaluation.molecule_probability_a),
            fmt(evaluation.molecule_probability_b),
            evaluation.probability_basis,
            "identical_common_observed_sites",
            shared_components(
                candidate_a.donor_genotype, candidate_b.donor_genotype),
            to_string(evaluation.common_sites),
            to_string(evaluation.discriminating_sites),
            to_string(max<long>(0, evaluation.common_sites -
                evaluation.discriminating_sites)),
            fmt(evaluation.common_depth),fmt(evaluation.discriminating_depth),
            to_string(preferred_sites),to_string(alternative_sites),
            to_string(evaluation.sites_neutral),
            fmt(evaluation.genotype_similarity),
            to_string(evaluation.chromosomes_covered),
            fmt(evaluation.effective_snps),fmt(evaluation.maximum_site_fraction),
            fmt(evaluation.top_five_site_fraction),
            to_string(evaluation.independent_molecules),
            fmt(evaluation.effective_molecules),
            fmt(evaluation.maximum_molecule_fraction),
            fmt(evaluation.umi_gene_molecule_fraction),
            evaluation.molecule_status,fmt(preferred_residual),
            fmt(alternative_residual),fit_status,
            fmt(evaluation.probability_without_top_site),
            fmt(evaluation.probability_without_top_five_sites),
            fmt(evaluation.probability_without_top_molecule),
            fmt(evaluation.minimum_leave_one_chromosome_out_probability),
            bool_text(evaluation.winner_changed_after_influence_removal),
            fmt(evaluation.site_bootstrap_win_fraction),
            fmt(evaluation.downsample_50pct_win_fraction),
            evaluation.downsample_basis,to_string(resamples),
            "2","TRUE","FALSE","TRUE",
            "NOT_APPLICABLE_TARGETED_PAIR","NA",
            fmt(evaluation.minimum_error_sensitivity_probability),
            bool_text(evaluation.error_sensitivity_stable),
            "NOT_EVALUATED_NO_FROZEN_ALLELE_PROFILE",fmt(e_ref),fmt(e_alt),
            warnings.empty() ? "NONE" : warnings,
            "identity_pair_probability_v3_reconciliation_targeted"
        });
    }
    if (gzclose(out) != Z_OK)
        die("failed closing probability output: " + output);
    fprintf(stderr,
        "Wrote targeted original-vs-reconciliation-swap probabilities to %s (%s molecule sidecar)\n",
        output.c_str(), molecule_sidecar ? "with" : "without");
}

// -------------------------------------------------------------------------
// Fixed-pair candidate-axis pilot (standalone, bounded-memory mode)
// -------------------------------------------------------------------------

static const char* AXIS_SCHEMA =
    "identity_candidate_axis_pair_score_v2_sampling_adjusted_fit_diagnostic";
static const char* AXIS_MANIFEST_SCHEMA =
    "identity_candidate_axis_pair_manifest_v1";
static const char* AXIS_FORMULA =
    "WEIGHTED_BRIER_FIXED_PAIR_PROJECTION_UNCLIPPED_V1";
static const char* AXIS_PREDICTION_TRANSFORM =
    "HARD_GT_DOSAGE_EXISTING_ASYMMETRIC_ERROR_TRANSFORM_V1";
static const char* AXIS_BASIS = "NUCLEAR_SITE_BALANCED_FIXED_PRIMARY";
static const char* AXIS_BASIS_INTERPRETATION =
    "SITE_BALANCED_FALLBACK_PROTOTYPE_NOT_MOLECULE_INDEPENDENCE";
static const char* AXIS_FOLD_BASIS = "GENOMIC_SITE_GROUP";
static const char* AXIS_FOLD_VERSION =
    "CANDIDATE_AXIS_GREEDY_DESIGN_MASS_SITE_GROUPS_PROJECT_FNV1A64_COMPAT_V1";
static const char* AXIS_TOLERANCE_VERSION =
    "IEEE754_LONG_DOUBLE_SCALE64_V1";
static const char* AXIS_OPERATIONAL_CONTRACT =
    "ORIGINAL_ALLOWED_DEMUX_VS_RECONCILIATION_NOMINATED_SWAP_ONLY";
static const char* AXIS_RETAINED_CONTRACT =
    "ORIGINAL_ALLOWED_DEMUX_VS_FROZEN_SUPPORTED_EVENT_PROPOSAL_CONTRAST_ONLY";

static long long strict_ll(const string& raw, const string& context){
    const string value = trim(raw);
    errno = 0;
    char* end = NULL;
    const long long parsed = strtoll(value.c_str(), &end, 10);
    if (value.empty() || errno != 0 || end == value.c_str() || *end != '\0')
        throw runtime_error(context + ": expected an integer, saw '" + raw + "'");
    return parsed;
}

static unsigned long strict_barcode_number(
        const string& raw, const string& context){
    const string value = trim(raw);
    errno = 0;
    char* end = NULL;
    const unsigned long parsed = strtoul(value.c_str(), &end, 10);
    if (value.empty() || value[0] == '-' || errno != 0 ||
            end == value.c_str() || *end != '\0')
        throw runtime_error(context + ": expected a nonnegative barcode integer, saw '" + raw + "'");
    return parsed;
}

static long double strict_ld(const string& raw, const string& context){
    const string value = trim(raw);
    errno = 0;
    char* end = NULL;
    const long double parsed = strtold(value.c_str(), &end);
    if (value.empty() || errno != 0 || end == value.c_str() || *end != '\0' ||
            !isfinite(parsed))
        throw runtime_error(context + ": expected a finite number, saw '" + raw + "'");
    return parsed;
}

static string axis_fmt(long double value){
    if (!isfinite(value)) return "NA";
    ostringstream out;
    out.precision(numeric_limits<long double>::max_digits10);
    out << value;
    return out.str();
}

static string axis_fmt6(double value){
    if (!isfinite(value)) return "NA";
    char buffer[128];
    snprintf(buffer, sizeof(buffer), "%.6g", value);
    return string(buffer);
}

static bool axis_true(const string& raw){
    const string value = lowercase(trim(raw));
    return value == "true" || value == "1" || value == "yes" || value == "y";
}

static bool axis_false(const string& raw){
    const string value = lowercase(trim(raw));
    return value == "false" || value == "0" || value == "no" || value == "n";
}

struct AxisKahan {
    long double value = 0.0L;
    long double correction = 0.0L;
    void add(long double item){
        const long double adjusted = item - correction;
        const long double next = value + adjusted;
        correction = (next - value) - adjusted;
        value = next;
    }
};

struct AxisObservationRecord {
    unsigned long barcode = 0;
    int32_t tid = -1;
    int32_t pos = -1;
    double ref = 0.0;
    double alt = 0.0;
};

static bool axis_observation_less(
        const AxisObservationRecord& left,
        const AxisObservationRecord& right){
    if (left.barcode != right.barcode) return left.barcode < right.barcode;
    if (left.tid != right.tid) return left.tid < right.tid;
    return left.pos < right.pos;
}

struct AxisSiteDefinition {
    uint64_t key = 0;
    int32_t tid = -1;
    int32_t pos = -1;
    string contig;
    string ref_allele;
    string alt_allele;
    bool found = false;
    bool mitochondrial = false;
    vector<int8_t> genotype;
};

struct CandidateAxisPair {
    size_t original = (size_t)-1;
    size_t proposed = (size_t)-1;
    string pair_id;
};

static void axis_require_equal(
        const string& field,
        const string& left,
        const string& right,
        const string& barcode){
    if (left != right)
        throw runtime_error("candidate-axis pair metadata mismatch for barcode=" +
            barcode + " field=" + field + " left='" + left + "' right='" + right + "'");
}

static CandidateAxisPair candidate_axis_pair(
        const vector<CandidateHypothesis>& candidates,
        const vector<size_t>& indices,
        unsigned long encoded_barcode,
        const string& libname){
    if (indices.size() != 2)
        throw runtime_error("candidate-axis manifest must contain exactly two rows per barcode; encoded_barcode=" +
            to_string(encoded_barcode) + " rows=" + to_string(indices.size()));
    CandidateAxisPair result;
    for (size_t index : indices){
        const CandidateHypothesis& candidate = candidates[index];
        if (candidate.schema_version != AXIS_MANIFEST_SCHEMA)
            throw runtime_error("candidate-axis manifest has unsupported schema for barcode=" + candidate.barcode);
        if (candidate.score_pair_id.empty())
            throw runtime_error("candidate-axis manifest is missing score_pair_id for barcode=" + candidate.barcode);
        if (result.pair_id.empty()) result.pair_id = candidate.score_pair_id;
        if (candidate.score_pair_id != result.pair_id)
            throw runtime_error("candidate-axis manifest has inconsistent score_pair_id values for barcode=" + candidate.barcode);
        if (candidate.score_pair_role == "ORIGINAL_ALLOWED_DEMUX"){
            if (result.original != (size_t)-1)
                throw runtime_error("duplicate ORIGINAL_ALLOWED_DEMUX role for barcode=" + candidate.barcode);
            result.original = index;
        } else if (candidate.score_pair_role == "RECONCILIATION_NOMINATED_SWAP" ||
                   candidate.score_pair_role == "FROZEN_SUPPORTED_EVENT_PROPOSAL_CONTRAST"){
            if (result.proposed != (size_t)-1)
                throw runtime_error("duplicate candidate-B role for barcode=" + candidate.barcode);
            result.proposed = index;
        } else {
            throw runtime_error("unsupported candidate-axis role for barcode=" + candidate.barcode +
                " role=" + candidate.score_pair_role);
        }
    }
    if (result.original == (size_t)-1 || result.proposed == (size_t)-1)
        throw runtime_error("candidate-axis pair must contain one A and one B role; encoded_barcode=" +
            to_string(encoded_barcode));
    const CandidateHypothesis& a = candidates[result.original];
    const CandidateHypothesis& b = candidates[result.proposed];
    if (a.library != libname || b.library != libname)
        throw runtime_error("candidate-axis manifest library does not match --libname for barcode=" + a.barcode);
    if (a.barcode != b.barcode)
        throw runtime_error("candidate-axis pair rows have different barcode text");
    if (!a.scoreable || !b.scoreable)
        throw runtime_error("candidate-axis pair contains an identity absent from the nuclear sample vector for barcode=" + a.barcode);
    if (a.donor_genotype == b.donor_genotype)
        throw runtime_error("candidate-axis pair contains identical A/B identities for barcode=" + a.barcode);
    if (a.donor_genotype != a.original_demux_assignment ||
            b.original_demux_assignment != a.donor_genotype ||
            a.current_donor_genotype != a.donor_genotype ||
            b.current_donor_genotype != a.donor_genotype)
        throw runtime_error("candidate-axis A is not the frozen original allowed demux identity for barcode=" + a.barcode);
    if (a.candidate_origin != "ORIGINAL_ALLOWED_DEMUX")
        throw runtime_error("candidate-axis A has invalid origin for barcode=" + a.barcode);
    if (b.donor_genotype != b.candidate_b_fixed_identity ||
            b.donor_genotype != b.selected_supported_event_proposal)
        throw runtime_error("candidate-axis B does not equal the selected fixed proposal for barcode=" + a.barcode);
    const bool retained = b.score_pair_role ==
        "FROZEN_SUPPORTED_EVENT_PROPOSAL_CONTRAST";
    if (retained){
        if (a.score_scope_contract != AXIS_RETAINED_CONTRACT ||
                b.score_scope_contract != AXIS_RETAINED_CONTRACT ||
                a.pair_construction_mode != "SUPPORTED_EVENT_CHALLENGE" ||
                b.pair_construction_mode != "SUPPORTED_EVENT_CHALLENGE" ||
                a.score_population_scope != "RETAINED_ORIGINAL_CONTRAST_ONLY" ||
                b.score_population_scope != "RETAINED_ORIGINAL_CONTRAST_ONLY" ||
                a.candidate_origin != "ORIGINAL_ALLOWED_DEMUX" ||
                b.candidate_origin != "FROZEN_SUPPORTED_EVENT_PROPOSAL_CONTRAST" ||
                axis_true(a.population_votes_in_authoritative_event) ||
                axis_true(b.population_votes_in_authoritative_event))
            throw runtime_error("retained candidate-axis pair violates the retained-contrast contract for barcode=" + a.barcode);
    } else {
        if (a.score_scope_contract != AXIS_OPERATIONAL_CONTRACT ||
                b.score_scope_contract != AXIS_OPERATIONAL_CONTRACT ||
                a.pair_construction_mode != "RECONCILIATION_NOMINATED_SWAP" ||
                b.pair_construction_mode != "RECONCILIATION_NOMINATED_SWAP" ||
                b.candidate_origin != "RECONCILIATION_NOMINATED_SWAP" ||
                b.donor_genotype != b.reconciliation_nominated_swap)
            throw runtime_error("operational candidate-axis pair violates the nominated-swap contract for barcode=" + a.barcode);
        const bool reassignment = a.score_population_scope == "APPLIED_REASSIGNMENT" ||
            a.score_population_scope == "RECOMMENDED_NOT_APPLIED";
        if (reassignment != axis_true(a.population_votes_in_authoritative_event) ||
                reassignment != axis_true(b.population_votes_in_authoritative_event))
            throw runtime_error("candidate-axis voting annotation disagrees with population scope for barcode=" + a.barcode);
    }
    if (!axis_true(a.population_votes_in_authoritative_event) &&
            !axis_false(a.population_votes_in_authoritative_event))
        throw runtime_error("candidate-axis voting annotation must be explicit TRUE/FALSE for barcode=" + a.barcode);
    const vector<pair<string,pair<string,string>>> shared = {
        {"library", {a.library,b.library}},
        {"score_pair_id", {a.score_pair_id,b.score_pair_id}},
        {"score_population_scope", {a.score_population_scope,b.score_population_scope}},
        {"population_votes_in_authoritative_event", {a.population_votes_in_authoritative_event,b.population_votes_in_authoritative_event}},
        {"supported_event_key", {a.supported_event_key,b.supported_event_key}},
        {"selected_supported_event_id", {a.selected_supported_event_id,b.selected_supported_event_id}},
        {"selected_supported_event_proposal", {a.selected_supported_event_proposal,b.selected_supported_event_proposal}},
        {"reconciliation_event_id", {a.reconciliation_event_id,b.reconciliation_event_id}},
        {"reconciliation_event_class", {a.reconciliation_event_class,b.reconciliation_event_class}},
        {"reconciliation_event_confidence", {a.reconciliation_event_confidence,b.reconciliation_event_confidence}},
        {"reconciliation_final_action", {a.reconciliation_final_action,b.reconciliation_final_action}},
        {"reconciliation_decision_confidence", {a.reconciliation_decision_confidence,b.reconciliation_decision_confidence}},
        {"reconciliation_reassignment_applied", {a.reconciliation_reassignment_applied,b.reconciliation_reassignment_applied}},
        {"original_demux_assignment", {a.original_demux_assignment,b.original_demux_assignment}},
        {"pair_construction_mode", {a.pair_construction_mode,b.pair_construction_mode}},
        {"score_scope_contract", {a.score_scope_contract,b.score_scope_contract}},
        {"schema_version", {a.schema_version,b.schema_version}}
    };
    for (const auto& field : shared)
        axis_require_equal(field.first, field.second.first, field.second.second, a.barcode);
    return result;
}

static string axis_parent(const string& path){
    const size_t slash = path.find_last_of('/');
    return slash == string::npos ? string(".") :
        (slash == 0 ? string("/") : path.substr(0, slash));
}

class AxisTempGuard {
public:
    explicit AxisTempGuard(const string& root){
        struct stat info;
        if (root.empty() || root[0] != '/' || stat(root.c_str(), &info) != 0 ||
                !S_ISDIR(info.st_mode))
            throw runtime_error("--candidate-axis-temp-dir must be an existing absolute directory: " + root);
        string pattern = root;
        if (pattern.back() != '/') pattern.push_back('/');
        pattern += "tetra_candidate_axis_XXXXXX";
        vector<char> buffer(pattern.begin(), pattern.end());
        buffer.push_back('\0');
        char* made = mkdtemp(buffer.data());
        if (!made) throw runtime_error("mkdtemp failed beneath " + root + ": " + strerror(errno));
        path_ = made;
        const string prefix = root.back() == '/' ? root : root + "/";
        if (path_.compare(0, prefix.size(), prefix) != 0){
            rmdir(path_.c_str());
            path_.clear();
            throw runtime_error("mkdtemp returned a child outside the requested temp root");
        }
    }
    ~AxisTempGuard(){ cleanup(); }
    const string& path() const { return path_; }
private:
    string path_;
    void cleanup(){
        if (path_.empty()) return;
        DIR* directory = opendir(path_.c_str());
        if (directory){
            struct dirent* entry = NULL;
            while ((entry = readdir(directory)) != NULL){
                const string name = entry->d_name;
                if (name == "." || name == "..") continue;
                const string child = path_ + "/" + name;
                unlink(child.c_str());
            }
            closedir(directory);
        }
        rmdir(path_.c_str());
        path_.clear();
    }
};

struct AxisUnit {
    int32_t tid = -1;
    int32_t pos = -1;
    long double ref = 0.0L;
    long double alt = 0.0L;
    long double y = 0.0L;
    long double p_a = 0.0L;
    long double p_b = 0.0L;
    long double a = 0.0L;
    long double b = 0.0L;
    long double n = 0.0L;
    long double d = 0.0L;
    long double m = 0.0L;
    long double ll_a = 0.0L;
    long double ll_b = 0.0L;
    bool discriminating = false;
};

struct AxisSums {
    long double w = 0.0L;
    long double n = 0.0L;
    long double d = 0.0L;
    long double sum_abs_n = 0.0L;
    long double sum_abs_m = 0.0L;
};

struct AxisNumericResult {
    string status = "NO_COMMON_NUCLEAR_EVIDENCE";
    string direction = "UNAVAILABLE";
    string segment = "UNAVAILABLE";
    long double position = NAN;
    long double margin = NAN;
    long double tau_d = NAN;
    long double tau_m = NAN;
};

struct AxisSamplingBrier {
    long double observed = NAN;
    long double expected_sampling = NAN;
    long double excess = NAN;
};

static AxisSamplingBrier axis_sampling_brier(
        long double observed_fraction,
        long double expected_fraction,
        long double depth){
    AxisSamplingBrier result;
    if (depth <= 0.0L || observed_fraction < 0.0L ||
            observed_fraction > 1.0L || expected_fraction < 0.0L ||
            expected_fraction > 1.0L)
        return result;
    const long double difference = observed_fraction - expected_fraction;
    result.observed = difference * difference;
    result.expected_sampling =
        expected_fraction * (1.0L - expected_fraction) / depth;
    result.excess = result.observed - result.expected_sampling;
    return result;
}

static AxisNumericResult axis_numeric(const AxisSums& sums){
    AxisNumericResult result;
    const long double epsilon = numeric_limits<long double>::epsilon();
    result.tau_d = 64.0L * epsilon * max(sums.w, 1.0L);
    result.tau_m = 64.0L * epsilon *
        max(2.0L * sums.sum_abs_n + sums.d, 1.0L);
    if (sums.w == 0.0L) return result;
    if (sums.d <= result.tau_d){
        result.status = "INSUFFICIENT_CANDIDATE_SEPARATION";
        return result;
    }
    result.status = "AVAILABLE";
    result.position = 100.0L * sums.n / sums.d;
    result.margin = 2.0L * sums.n - sums.d;
    if (result.margin > result.tau_m) result.direction = "PROPOSAL_SIDE";
    else if (result.margin < -result.tau_m) result.direction = "ORIGINAL_SIDE";
    else result.direction = "TIE";
    if (result.position < 0.0L)
        result.segment = "BEYOND_ORIGINAL_AWAY_FROM_PROPOSAL";
    else if (result.position <= 100.0L)
        result.segment = "BETWEEN_FIXED_CANDIDATE_EXPECTATIONS";
    else result.segment = "BEYOND_PROPOSAL_AWAY_FROM_ORIGINAL";
    return result;
}

static AxisSums axis_sum_units(const vector<AxisUnit>& units){
    AxisKahan w, n, d, abs_n, abs_m;
    for (const AxisUnit& unit : units){
        w.add(1.0L); n.add(unit.n); d.add(unit.d);
        abs_n.add(fabsl(unit.n)); abs_m.add(fabsl(unit.m));
    }
    AxisSums result;
    result.w = w.value; result.n = n.value; result.d = d.value;
    result.sum_abs_n = abs_n.value; result.sum_abs_m = abs_m.value;
    return result;
}

static long double axis_probability_from_delta(long double delta){
    if (!isfinite(delta)) return NAN;
    if (delta >= 0.0L)
        return 100.0L / (1.0L + expl(-min(delta, 11350.0L)));
    const long double value = expl(max(delta, -11350.0L));
    return 100.0L * value / (1.0L + value);
}

static long double axis_binom(
        long double ref, long double alt, long double expected){
    const long double protection = 1e-12L;
    const long double q = min(max(expected, protection), 1.0L - protection);
    return alt * logl(q) + ref * logl(1.0L - q);
}

static string normalized_contig(string value){
    value = trim(value);
    if (value.size() >= 3 && lowercase(value.substr(0,3)) == "chr")
        value = value.substr(3);
    return lowercase(value);
}

static bool axis_expected(
        const Identity& identity,
        const AxisSiteDefinition& site,
        const unordered_map<int,size_t>& donor_slot,
        long double& expected){
    auto first = donor_slot.find(identity.a);
    if (first == donor_slot.end()) return false;
    const int ga = site.genotype[first->second];
    if (ga < 0 || ga > 2) return false;
    if (identity.b < 0 || identity.a == identity.b){
        expected = (long double)ga / 2.0L;
        return true;
    }
    auto second = donor_slot.find(identity.b);
    if (second == donor_slot.end()) return false;
    const int gb = site.genotype[second->second];
    if (gb < 0 || gb > 2) return false;
    expected = (long double)(ga + gb) / 4.0L;
    return true;
}

struct AxisCellResult {
    AxisSums sums;
    AxisNumericResult numeric;
    vector<AxisUnit> units;
    long n_unique_merged = 0;
    long n_duplicate_rows = 0;
    long excluded_mito = 0;
    long excluded_missing_a = 0;
    long excluded_missing_b = 0;
    long excluded_missing_both = 0;
    long excluded_missing_definition = 0;
    long excluded_nonpositive = 0;
    long common_sites = 0;
    long discriminating_sites = 0;
    long double total_common_depth = 0.0L;
    long double discriminating_depth = 0.0L;
    long double separation = NAN;
    long double similarity = NAN;
    long double ll_a = NAN;
    long double ll_b = NAN;
    long double ll_delta = NAN;
    long double legacy_probability_a = NAN;
    long double legacy_probability_b = NAN;
    long sites_favor_a = 0;
    long sites_favor_b = 0;
    long sites_tied = 0;
    long double residual_a = NAN;
    long double residual_b = NAN;
    long double candidate_a_observed_brier_mean = NAN;
    long double candidate_a_expected_sampling_brier_mean = NAN;
    long double candidate_a_excess_brier_mean = NAN;
    long double candidate_b_observed_brier_mean = NAN;
    long double candidate_b_expected_sampling_brier_mean = NAN;
    long double candidate_b_excess_brier_mean = NAN;
    string raw_residual_threshold_flag = "UNAVAILABLE";
    string comparison_status = "NO_COMMON_EVIDENCE";
    string absolute_fit_status = "UNAVAILABLE";
    long double legacy_without_top = NAN;
    long double legacy_without_top_five = NAN;
    long double minimum_error_probability = NAN;
    string error_stable = "FALSE";
    long double without_top_position = NAN;
    long double without_five_position = NAN;
    string without_top_direction = "UNAVAILABLE";
    string without_five_direction = "UNAVAILABLE";
    string preserve_top = "NA";
    string preserve_five = "NA";
    long removed_top = 0;
    long removed_five = 0;
    string removal_top_status = "FULL_SCORE_UNAVAILABLE";
    string removal_five_status = "FULL_SCORE_UNAVAILABLE";
    long double maximum_margin_fraction = NAN;
    long double top_five_margin_fraction = NAN;
    string concentration_status = "FULL_SCORE_UNAVAILABLE";
    string top_unit_id = "NA";
    string top_five_unit_ids = "NA";
    int fold_count = 0;
    int folds_evaluable = 0;
    long double fold_min = NAN;
    long double fold_median = NAN;
    long double fold_max = NAN;
    long double fold_proposal_fraction = NAN;
    string folds_preserved = "NA";
    string fold_status = "FULL_SCORE_UNAVAILABLE";
    vector<string> warnings;
};

static string axis_raw_residual_threshold_flag(
        long double residual_a,
        long double residual_b,
        long double legacy_threshold){
    if (!isfinite(residual_a) || !isfinite(residual_b))
        return "UNAVAILABLE";
    const bool a_above = residual_a > legacy_threshold;
    const bool b_above = residual_b > legacy_threshold;
    if (a_above && b_above) return "BOTH_ABOVE_LEGACY_THRESHOLD";
    if (a_above) return "CANDIDATE_A_ONLY_ABOVE_LEGACY_THRESHOLD";
    if (b_above) return "CANDIDATE_B_ONLY_ABOVE_LEGACY_THRESHOLD";
    return "NEITHER_ABOVE_LEGACY_THRESHOLD";
}

static string axis_unit_id(const AxisUnit& unit){
    return to_string(unit.tid) + ":" + to_string(unit.pos);
}

static long double axis_median(vector<long double> values){
    if (values.empty()) return NAN;
    sort(values.begin(), values.end());
    const size_t n = values.size();
    return n % 2 ? values[n/2] : (values[n/2-1] + values[n/2]) / 2.0L;
}

static AxisSums axis_remove(
        const AxisSums& full, const vector<AxisUnit>& units,
        const vector<size_t>& order, size_t count){
    AxisSums result = full;
    for (size_t i = 0; i < min(count, order.size()); ++i){
        const AxisUnit& unit = units[order[i]];
        result.w -= 1.0L;
        result.n -= unit.n;
        result.d -= unit.d;
        result.sum_abs_n -= fabsl(unit.n);
        result.sum_abs_m -= fabsl(unit.m);
    }
    result.w = max(result.w, 0.0L);
    result.sum_abs_n = max(result.sum_abs_n, 0.0L);
    result.sum_abs_m = max(result.sum_abs_m, 0.0L);
    return result;
}

static void axis_influence(AxisCellResult& result){
    if (result.numeric.status != "AVAILABLE") return;
    vector<size_t> order;
    for (size_t i = 0; i < result.units.size(); ++i)
        if (result.units[i].discriminating) order.push_back(i);
    sort(order.begin(), order.end(), [&](size_t left, size_t right){
        const long double a = fabsl(result.units[left].m);
        const long double b = fabsl(result.units[right].m);
        if (a != b) return a > b;
        if (result.units[left].tid != result.units[right].tid)
            return result.units[left].tid < result.units[right].tid;
        return result.units[left].pos < result.units[right].pos;
    });
    if (order.empty()){
        result.removal_top_status = result.removal_five_status =
            "NO_AXIS_DISCRIMINATING_UNITS";
        result.concentration_status = "NO_NONZERO_BRIER_MARGIN_CONTRIBUTIONS";
        return;
    }
    result.removed_top = 1;
    result.removed_five = (long)min<size_t>(5, order.size());
    result.top_unit_id = axis_unit_id(result.units[order[0]]);
    result.top_five_unit_ids.clear();
    for (size_t i = 0; i < (size_t)result.removed_five; ++i){
        if (!result.top_five_unit_ids.empty()) result.top_five_unit_ids += ",";
        result.top_five_unit_ids += axis_unit_id(result.units[order[i]]);
    }
    const AxisNumericResult top = axis_numeric(axis_remove(
        result.sums, result.units, order, 1));
    const AxisNumericResult five = axis_numeric(axis_remove(
        result.sums, result.units, order, 5));
    if (top.status == "AVAILABLE"){
        result.removal_top_status = "AVAILABLE";
        result.without_top_position = top.position;
        result.without_top_direction = top.direction;
        result.preserve_top = result.numeric.direction == "TIE" ? "NA" :
            bool_text(top.direction == result.numeric.direction);
    } else result.removal_top_status =
        "INSUFFICIENT_CANDIDATE_SEPARATION_AFTER_REMOVAL";
    if (five.status == "AVAILABLE"){
        result.removal_five_status = "AVAILABLE";
        result.without_five_position = five.position;
        result.without_five_direction = five.direction;
        result.preserve_five = result.numeric.direction == "TIE" ? "NA" :
            bool_text(five.direction == result.numeric.direction);
    } else result.removal_five_status =
        "INSUFFICIENT_CANDIDATE_SEPARATION_AFTER_REMOVAL";
    if (result.sums.sum_abs_m > 0.0L){
        result.concentration_status = "AVAILABLE";
        result.maximum_margin_fraction = fabsl(result.units[order[0]].m) /
            result.sums.sum_abs_m;
        long double top_five = 0.0L;
        for (size_t i = 0; i < (size_t)result.removed_five; ++i)
            top_five += fabsl(result.units[order[i]].m);
        result.top_five_margin_fraction = top_five / result.sums.sum_abs_m;
    } else result.concentration_status =
        "NO_NONZERO_BRIER_MARGIN_CONTRIBUTIONS";
}

static uint64_t axis_fold_hash(
        const string& library, const string& barcode,
        const AxisUnit& unit){
    return stable_text_hash(library + "|" + barcode + "|" +
        to_string(unit.tid) + "|" + to_string(unit.pos) + "|" +
        AXIS_FOLD_VERSION);
}

static bool axis_sum_nearly_equal(
        long double observed,
        long double expected,
        size_t n_terms){
    const long double scale = max(
        max(fabsl(observed), fabsl(expected)), 1.0L);
    const long double tolerance =
        64.0L * max<size_t>(n_terms, 1) *
        numeric_limits<long double>::epsilon() * scale;
    return fabsl(observed - expected) <= tolerance;
}

static void axis_folds(
        AxisCellResult& result,
        const string& library,
        const string& barcode){
    if (result.numeric.status != "AVAILABLE"){
        result.fold_status = "FULL_SCORE_UNAVAILABLE";
        return;
    }
    vector<size_t> positive, zero;
    for (size_t i = 0; i < result.units.size(); ++i){
        if (result.units[i].d > 0.0L) positive.push_back(i);
        else zero.push_back(i);
    }
    if (positive.size() < 2){
        result.fold_status = "INSUFFICIENT_POSITIVE_DESIGN_GROUPS";
        return;
    }
    result.fold_count = (int)min<size_t>(10, positive.size());
    sort(positive.begin(), positive.end(), [&](size_t left, size_t right){
        const AxisUnit& a = result.units[left];
        const AxisUnit& b = result.units[right];
        if (a.d != b.d) return a.d > b.d;
        const uint64_t ah = axis_fold_hash(library, barcode, a);
        const uint64_t bh = axis_fold_hash(library, barcode, b);
        if (ah != bh) return ah < bh;
        if (a.tid != b.tid) return a.tid < b.tid;
        return a.pos < b.pos;
    });
    sort(zero.begin(), zero.end(), [&](size_t left, size_t right){
        const uint64_t ah = axis_fold_hash(library, barcode, result.units[left]);
        const uint64_t bh = axis_fold_hash(library, barcode, result.units[right]);
        if (ah != bh) return ah < bh;
        if (result.units[left].tid != result.units[right].tid)
            return result.units[left].tid < result.units[right].tid;
        return result.units[left].pos < result.units[right].pos;
    });
    vector<AxisSums> fold(result.fold_count);
    vector<long> counts(result.fold_count, 0);
    auto add = [&](int target, const AxisUnit& unit){
        fold[target].w += 1.0L;
        fold[target].n += unit.n;
        fold[target].d += unit.d;
        fold[target].sum_abs_n += fabsl(unit.n);
        fold[target].sum_abs_m += fabsl(unit.m);
        ++counts[target];
    };
    for (size_t index : positive){
        int best = 0;
        for (int f = 1; f < result.fold_count; ++f){
            const tuple<long double,long,long,int> candidate(
                fold[f].d, counts[f], counts[f], f);
            const tuple<long double,long,long,int> incumbent(
                fold[best].d, counts[best], counts[best], best);
            if (candidate < incumbent) best = f;
        }
        add(best, result.units[index]);
    }
    for (size_t index : zero){
        int best = 0;
        for (int f = 1; f < result.fold_count; ++f){
            const tuple<long double,long,int> candidate(
                fold[f].w, counts[f], f);
            const tuple<long double,long,int> incumbent(
                fold[best].w, counts[best], best);
            if (candidate < incumbent) best = f;
        }
        add(best, result.units[index]);
    }
    AxisSums reconstructed;
    for (const AxisSums& item : fold){
        reconstructed.w += item.w; reconstructed.n += item.n;
        reconstructed.d += item.d;
        reconstructed.sum_abs_n += item.sum_abs_n;
        reconstructed.sum_abs_m += item.sum_abs_m;
    }
    const size_t n_terms = result.units.size();
    const bool reconstructed_ok =
        reconstructed.w == result.sums.w &&
        axis_sum_nearly_equal(reconstructed.n, result.sums.n, n_terms) &&
        axis_sum_nearly_equal(reconstructed.d, result.sums.d, n_terms) &&
        axis_sum_nearly_equal(
            reconstructed.sum_abs_n, result.sums.sum_abs_n, n_terms) &&
        axis_sum_nearly_equal(
            reconstructed.sum_abs_m, result.sums.sum_abs_m, n_terms);
    if (!reconstructed_ok){
        result.fold_status = "FOLD_RECONSTRUCTION_MISMATCH";
        result.folds_preserved = "NA";
        result.warnings.push_back("FOLD_RECONSTRUCTION_MISMATCH");
        return;
    }
    vector<long double> positions;
    long proposal = 0;
    bool unavailable = false;
    bool preserved = true;
    for (int f = 0; f < result.fold_count; ++f){
        AxisSums remaining = result.sums;
        remaining.w -= fold[f].w; remaining.n -= fold[f].n;
        remaining.d -= fold[f].d;
        remaining.sum_abs_n -= fold[f].sum_abs_n;
        remaining.sum_abs_m -= fold[f].sum_abs_m;
        remaining.w = max(remaining.w, 0.0L);
        remaining.sum_abs_n = max(remaining.sum_abs_n, 0.0L);
        remaining.sum_abs_m = max(remaining.sum_abs_m, 0.0L);
        const AxisNumericResult score = axis_numeric(remaining);
        if (score.status != "AVAILABLE"){
            unavailable = true;
            continue;
        }
        ++result.folds_evaluable;
        positions.push_back(score.position);
        if (score.direction == "PROPOSAL_SIDE") ++proposal;
        if (score.direction != result.numeric.direction) preserved = false;
    }
    if (unavailable || result.folds_evaluable != result.fold_count){
        result.fold_status = "FOLD_UNAVAILABLE";
        result.folds_preserved = "NA";
        return;
    }
    sort(positions.begin(), positions.end());
    result.fold_min = positions.front();
    result.fold_median = axis_median(positions);
    result.fold_max = positions.back();
    result.fold_proposal_fraction = (long double)proposal /
        (long double)result.folds_evaluable;
    if (result.numeric.direction == "TIE"){
        result.fold_status = "FULL_SCORE_TIE";
        result.folds_preserved = "NA";
    } else {
        result.folds_preserved = bool_text(preserved);
        result.fold_status = preserved ? "PRESERVED_ALL" :
            "DIRECTION_CHANGED_OR_TIED";
    }
}

static AxisCellResult evaluate_candidate_axis_cell(
        const CandidateHypothesis& candidate_a,
        const CandidateHypothesis& candidate_b,
        const vector<AxisObservationRecord>& merged,
        const vector<long>& duplicate_counts,
        const vector<AxisSiteDefinition>& sites,
        const unordered_map<int,size_t>& donor_slot,
        long double e_ref,
        long double e_alt,
        long min_evidence,
        long double poor_fit_residual){
    AxisCellResult result;
    AxisKahan residual_a, residual_b, ll_a, ll_b;
    AxisKahan observed_brier_a, expected_sampling_brier_a, excess_brier_a;
    AxisKahan observed_brier_b, expected_sampling_brier_b, excess_brier_b;
    for (size_t i = 0; i < merged.size(); ++i){
        const AxisObservationRecord& observation = merged[i];
        ++result.n_unique_merged;
        result.n_duplicate_rows += duplicate_counts[i];
        const uint64_t key = site_key(observation.tid, observation.pos);
        auto found = lower_bound(sites.begin(), sites.end(), key,
            [](const AxisSiteDefinition& site, uint64_t value){ return site.key < value; });
        if (found == sites.end() || found->key != key || !found->found){
            ++result.excluded_missing_definition;
            continue;
        }
        if (found->mitochondrial){ ++result.excluded_mito; continue; }
        const long double depth = (long double)observation.ref + observation.alt;
        if (depth <= 0.0L){ ++result.excluded_nonpositive; continue; }
        long double p_a = 0.0L, p_b = 0.0L;
        const bool has_a = axis_expected(candidate_a.identity, *found, donor_slot, p_a);
        const bool has_b = axis_expected(candidate_b.identity, *found, donor_slot, p_b);
        if (!has_a && !has_b){ ++result.excluded_missing_both; continue; }
        if (!has_a){ ++result.excluded_missing_a; continue; }
        if (!has_b){ ++result.excluded_missing_b; continue; }
        const long double observed = observation.alt / depth;
        const long double adjusted_a = p_a * (1.0L - e_alt) +
            (1.0L - p_a) * e_ref;
        const long double adjusted_b = p_b * (1.0L - e_alt) +
            (1.0L - p_b) * e_ref;
        if (observed < 0.0L || observed > 1.0L || adjusted_a < 0.0L ||
                adjusted_a > 1.0L || adjusted_b < 0.0L || adjusted_b > 1.0L)
            throw runtime_error("candidate-axis prediction/observation out of [0,1] for barcode=" + candidate_a.barcode);
        AxisUnit unit;
        unit.tid = observation.tid; unit.pos = observation.pos;
        unit.ref = observation.ref; unit.alt = observation.alt;
        unit.y = observed; unit.p_a = p_a; unit.p_b = p_b;
        unit.a = adjusted_a; unit.b = adjusted_b;
        const long double delta = adjusted_b - adjusted_a;
        unit.n = delta * (observed - adjusted_a);
        unit.d = delta * delta;
        unit.m = 2.0L * unit.n - unit.d;
        unit.discriminating = p_a != p_b;
        result.units.push_back(unit);
        ++result.common_sites;
        result.total_common_depth += depth;
        if (!unit.discriminating) continue;
        ++result.discriminating_sites;
        result.discriminating_depth += depth;
        const AxisSamplingBrier brier_a = axis_sampling_brier(
            observed, adjusted_a, depth);
        const AxisSamplingBrier brier_b = axis_sampling_brier(
            observed, adjusted_b, depth);
        observed_brier_a.add(brier_a.observed);
        expected_sampling_brier_a.add(brier_a.expected_sampling);
        excess_brier_a.add(brier_a.excess);
        observed_brier_b.add(brier_b.observed);
        expected_sampling_brier_b.add(brier_b.expected_sampling);
        excess_brier_b.add(brier_b.excess);
        residual_a.add(depth * fabsl(observed - adjusted_a));
        residual_b.add(depth * fabsl(observed - adjusted_b));
        unit.ll_a = axis_binom(observation.ref, observation.alt, adjusted_a);
        unit.ll_b = axis_binom(observation.ref, observation.alt, adjusted_b);
        result.units.back().ll_a = unit.ll_a;
        result.units.back().ll_b = unit.ll_b;
        ll_a.add(unit.ll_a); ll_b.add(unit.ll_b);
        const long double delta_ll = unit.ll_a - unit.ll_b;
        if (delta_ll > 1e-12L) ++result.sites_favor_a;
        else if (delta_ll < -1e-12L) ++result.sites_favor_b;
        else ++result.sites_tied;
    }
    const long classified = result.excluded_missing_definition +
        result.excluded_mito + result.excluded_nonpositive +
        result.excluded_missing_both + result.excluded_missing_a +
        result.excluded_missing_b + result.common_sites;
    if (classified != result.n_unique_merged)
        throw runtime_error("candidate-axis site-accounting identity failed for barcode=" + candidate_a.barcode);
    result.sums = axis_sum_units(result.units);
    result.numeric = axis_numeric(result.sums);
    if (result.sums.w > 0.0L){
        result.separation = 100.0L * sqrtl(max(result.sums.d, 0.0L) /
            result.sums.w);
        result.similarity = 100.0L - result.separation;
    }
    if (result.numeric.status == "AVAILABLE"){
        AxisKahan margin_units;
        for (const AxisUnit& unit : result.units) margin_units.add(unit.m);
        if (fabsl(margin_units.value - result.numeric.margin) > result.numeric.tau_m)
            throw runtime_error("candidate-axis Brier-margin identity failed for barcode=" + candidate_a.barcode);
    }
    if (result.discriminating_sites > 0){
        const long double denominator =
            (long double)result.discriminating_sites;
        result.candidate_a_observed_brier_mean =
            observed_brier_a.value / denominator;
        result.candidate_a_expected_sampling_brier_mean =
            expected_sampling_brier_a.value / denominator;
        result.candidate_a_excess_brier_mean =
            excess_brier_a.value / denominator;
        result.candidate_b_observed_brier_mean =
            observed_brier_b.value / denominator;
        result.candidate_b_expected_sampling_brier_mean =
            expected_sampling_brier_b.value / denominator;
        result.candidate_b_excess_brier_mean =
            excess_brier_b.value / denominator;
        result.ll_a = ll_a.value; result.ll_b = ll_b.value;
        result.ll_delta = result.ll_a - result.ll_b;
        result.legacy_probability_a = axis_probability_from_delta(result.ll_delta);
        result.legacy_probability_b = 100.0L - result.legacy_probability_a;
        if (result.discriminating_depth > 0.0L){
            result.residual_a = residual_a.value / result.discriminating_depth;
            result.residual_b = residual_b.value / result.discriminating_depth;
        }
    }
    if (result.common_sites == 0) result.comparison_status = "NO_COMMON_EVIDENCE";
    else if (result.discriminating_sites == 0)
        result.comparison_status = "PANEL_NONDISCRIMINATING";
    else if (isfinite(result.residual_a) && isfinite(result.residual_b) &&
            result.residual_a > poor_fit_residual &&
            result.residual_b > poor_fit_residual)
        result.comparison_status = "NO_CANDIDATE_FITS";
    else if (result.discriminating_depth < min_evidence ||
            result.discriminating_sites < 2)
        result.comparison_status = "LOW_EVIDENCE";
    else result.comparison_status = "PASS";
    if (result.discriminating_sites > 0){
        const long double preferred_residual = result.ll_delta >= 0.0L ?
            result.residual_a : result.residual_b;
        result.absolute_fit_status = result.comparison_status == "NO_CANDIDATE_FITS" ?
            "NO_CANDIDATE_FITS" :
            (isfinite(preferred_residual) && preferred_residual <= poor_fit_residual ?
             "PASS" : "POOR_FIT");
    }
    result.raw_residual_threshold_flag = axis_raw_residual_threshold_flag(
        result.residual_a, result.residual_b, poor_fit_residual);

    vector<size_t> legacy_order;
    for (size_t i = 0; i < result.units.size(); ++i)
        if (result.units[i].discriminating) legacy_order.push_back(i);
    sort(legacy_order.begin(), legacy_order.end(), [&](size_t left, size_t right){
        const long double a = fabsl(result.units[left].ll_a - result.units[left].ll_b);
        const long double b = fabsl(result.units[right].ll_a - result.units[right].ll_b);
        if (a != b) return a > b;
        if (result.units[left].tid != result.units[right].tid)
            return result.units[left].tid < result.units[right].tid;
        return result.units[left].pos < result.units[right].pos;
    });
    if (isfinite(result.ll_delta) && !legacy_order.empty()){
        const int preferred_sign = result.ll_delta >= 0.0L ? 1 : -1;
        long double removed = result.units[legacy_order[0]].ll_a -
            result.units[legacy_order[0]].ll_b;
        result.legacy_without_top = axis_probability_from_delta(
            preferred_sign * (result.ll_delta - removed));
        removed = 0.0L;
        for (size_t i = 0; i < min<size_t>(5, legacy_order.size()); ++i)
            removed += result.units[legacy_order[i]].ll_a - result.units[legacy_order[i]].ll_b;
        result.legacy_without_top_five = axis_probability_from_delta(
            preferred_sign * (result.ll_delta - removed));
        const long double error_values[] = {e_ref, 0.005L, 0.01L, 0.02L};
        result.minimum_error_probability = 100.0L;
        for (long double error : error_values){
            AxisKahan delta;
            for (size_t index : legacy_order){
                const AxisUnit& unit = result.units[index];
                const long double a = unit.p_a * (1.0L - error) +
                    (1.0L - unit.p_a) * error;
                const long double b = unit.p_b * (1.0L - error) +
                    (1.0L - unit.p_b) * error;
                delta.add(axis_binom(unit.ref, unit.alt, a) -
                    axis_binom(unit.ref, unit.alt, b));
            }
            result.minimum_error_probability = min(
                result.minimum_error_probability,
                axis_probability_from_delta(preferred_sign * delta.value));
        }
        result.error_stable = bool_text(result.minimum_error_probability > 50.0L);
    }
    axis_influence(result);
    axis_folds(result, candidate_a.library, candidate_a.barcode);
    return result;
}

struct AxisResourceAudit {
    string status = "PASS";
    unsigned long long target_rows = 0;
    unsigned long long target_barcodes = 0;
    unsigned long long unique_cell_sites = 0;
    unsigned long long duplicate_rows = 0;
    unsigned long long zero_depth_sites = 0;
    unsigned long long malformed_target_rows = 0;
    unsigned long long mitochondrial_sites = 0;
    unsigned long long missing_site_definitions = 0;
    unsigned long long missing_both_candidates = 0;
    unsigned long long missing_candidate_a_only = 0;
    unsigned long long missing_candidate_b_only = 0;
    unsigned long long common_nuclear_sites = 0;
    unsigned long long unique_site_keys = 0;
    unsigned long long spill_runs = 0;
    unsigned long long bucket_count = 0;
    unsigned long long largest_barcode_rows = 0;
    unsigned long long estimated_largest_bucket_rows = 0;
    unsigned long long observed_peak_bucket_rows = 0;
    unsigned long long selected_key_bytes = 0;
    unsigned long long compact_site_definition_bytes = 0;
};

static void axis_write_resource_audit(
        const string& path,
        const string& library,
        const string& temp_root,
        const AxisResourceAudit& audit){
    const string temporary = path + ".tmp." + to_string((long long)getpid());
    ofstream out(temporary.c_str());
    if (!out) throw runtime_error("could not write candidate-axis resource audit: " + temporary);
    out << "status\tcheck\tvalue\tdetail\tschema_version\n";
    const vector<pair<string,string>> values = {
        {"library", library},
        {"candidate_axis_temp_root", temp_root},
        {"candidate_axis_temp_policy", "MKDTEMP_UNIQUE_CHILD_CLEAN_EXACT_CHILD_ONLY"},
        {"n_target_observation_rows", to_string(audit.target_rows)},
        {"n_distinct_target_barcodes", to_string(audit.target_barcodes)},
        {"n_unique_merged_target_cell_sites", to_string(audit.unique_cell_sites)},
        {"n_duplicate_observation_rows_merged", to_string(audit.duplicate_rows)},
        {"n_sites_excluded_nonpositive_observation", to_string(audit.zero_depth_sites)},
        {"n_malformed_target_observation_rows", to_string(audit.malformed_target_rows)},
        {"n_sites_excluded_mitochondrial", to_string(audit.mitochondrial_sites)},
        {"n_sites_excluded_missing_site_definition", to_string(audit.missing_site_definitions)},
        {"n_sites_excluded_missing_both_candidates", to_string(audit.missing_both_candidates)},
        {"n_sites_excluded_missing_candidate_a_only", to_string(audit.missing_candidate_a_only)},
        {"n_sites_excluded_missing_candidate_b_only", to_string(audit.missing_candidate_b_only)},
        {"n_common_observed_nuclear_sites", to_string(audit.common_nuclear_sites)},
        {"candidate_axis_unique_selected_site_keys", to_string(audit.unique_site_keys)},
        {"candidate_axis_first_pass_spill_runs", to_string(audit.spill_runs)},
        {"candidate_axis_bucket_count", to_string(audit.bucket_count)},
        {"candidate_axis_largest_barcode_rows", to_string(audit.largest_barcode_rows)},
        {"candidate_axis_estimated_largest_bucket_rows", to_string(audit.estimated_largest_bucket_rows)},
        {"candidate_axis_observed_peak_bucket_rows", to_string(audit.observed_peak_bucket_rows)},
        {"selected_site_key_bytes", to_string(audit.selected_key_bytes)},
        {"compact_site_definition_bytes", to_string(audit.compact_site_definition_bytes)},
        {"observation_pass_contract", "EXACTLY_TWO_PILEUP_OBSERVATION_PASSES_ONE_PILEUP_SITE_PASS"}
    };
    for (const auto& item : values)
        out << audit.status << '\t' << item.first << '\t' << item.second <<
            "\tMEASURED_OR_DERIVED_BEFORE_BUCKET_SCORING\tidentity_candidate_axis_resource_evidence_audit_v1\n";
    out.close();
    if (!out) throw runtime_error("failed closing candidate-axis resource audit: " + temporary);
    if (rename(temporary.c_str(), path.c_str()) != 0){
        unlink(temporary.c_str());
        throw runtime_error("failed publishing candidate-axis resource audit: " + path);
    }
}

static void axis_write_binary_record(
        ofstream& out, const AxisObservationRecord& record,
        const string& path){
    out.write(reinterpret_cast<const char*>(&record), sizeof(record));
    if (!out) throw runtime_error("failed writing candidate-axis temporary record: " + path);
}

static bool axis_read_binary_record(
        ifstream& in, AxisObservationRecord& record,
        const string& path){
    in.read(reinterpret_cast<char*>(&record), sizeof(record));
    if (in.gcount() == 0 && in.eof()) return false;
    if (in.gcount() != (streamsize)sizeof(record))
        throw runtime_error("truncated candidate-axis temporary record file: " + path);
    return true;
}

static string axis_spill_run(
        vector<AxisObservationRecord>& records,
        const string& temp_path,
        size_t run_index){
    sort(records.begin(), records.end(), axis_observation_less);
    const string path = temp_path + "/first_pass_run_" +
        to_string(run_index) + ".bin";
    ofstream out(path.c_str(), ios::binary);
    if (!out) throw runtime_error("could not create candidate-axis spill run: " + path);
    for (const AxisObservationRecord& record : records)
        axis_write_binary_record(out, record, path);
    out.close();
    if (!out) throw runtime_error("failed closing candidate-axis spill run: " + path);
    records.clear();
    return path;
}

static string axis_spill_key_run(
        vector<uint64_t>& keys,
        const string& temp_path,
        size_t run_index){
    sort(keys.begin(), keys.end());
    keys.erase(unique(keys.begin(), keys.end()), keys.end());
    const string path = temp_path + "/selected_key_run_" +
        to_string(run_index) + ".bin";
    ofstream out(path.c_str(), ios::binary);
    if (!out) throw runtime_error(
        "could not create candidate-axis selected-key spill run: " + path);
    for (uint64_t key : keys)
        out.write(reinterpret_cast<const char*>(&key), sizeof(key));
    out.close();
    if (!out) throw runtime_error(
        "failed closing candidate-axis selected-key spill run: " + path);
    vector<uint64_t>().swap(keys);
    return path;
}

static bool axis_read_key(ifstream& input, uint64_t& key, const string& path){
    input.read(reinterpret_cast<char*>(&key), sizeof(key));
    if (input.gcount() == 0 && input.eof()) return false;
    if (input.gcount() != (streamsize)sizeof(key))
        throw runtime_error(
            "truncated candidate-axis selected-key spill run: " + path);
    return true;
}

struct AxisKeyCursor {
    uint64_t key = 0;
    size_t run = 0;
};

struct AxisKeyCursorGreater {
    bool operator()(const AxisKeyCursor& left, const AxisKeyCursor& right) const {
        if (left.key != right.key) return left.key > right.key;
        return left.run > right.run;
    }
};

static AxisObservationRecord axis_parse_observation(
        const vector<string>& fields,
        const string& path,
        unsigned long long line_no){
    if (fields.size() != 5)
        throw runtime_error(path + ": target observation row must have exactly five fields at line " + to_string(line_no));
    const string context = path + ": line " + to_string(line_no);
    AxisObservationRecord record;
    record.barcode = strict_barcode_number(fields[0], context + " barcode");
    const long long tid = strict_ll(fields[1], context + " tid");
    const long long pos = strict_ll(fields[2], context + " position");
    if (tid < 0 || tid > INT_MAX || pos < 0 || pos > INT_MAX)
        throw runtime_error(context + ": tid and position must be nonnegative 32-bit integers");
    record.tid = (int32_t)tid; record.pos = (int32_t)pos;
    const long double ref = strict_ld(fields[3], context + " REF");
    const long double alt = strict_ld(fields[4], context + " ALT");
    if (ref < 0.0L || alt < 0.0L)
        throw runtime_error(context + ": REF and ALT must be nonnegative");
    record.ref = (double)ref; record.alt = (double)alt;
    return record;
}

struct AxisRunCursor {
    AxisObservationRecord record;
    size_t run = 0;
};

struct AxisRunCursorGreater {
    bool operator()(const AxisRunCursor& left, const AxisRunCursor& right) const {
        return axis_observation_less(right.record, left.record);
    }
};

static void axis_first_observation_pass(
        const string& path,
        const unordered_map<unsigned long,CandidateAxisPair>& pairs,
        const string& temp_path,
        unsigned long long chunk_bytes,
        AxisResourceAudit& audit,
        unordered_map<unsigned long,unsigned long long>& rows_by_barcode,
        vector<string>& runs){
    gzFile input = gzopen(path.c_str(), "rb");
    if (!input) throw runtime_error("could not open candidate-axis observations: " + path);
    vector<AxisObservationRecord> chunk;
    const size_t capacity = max<size_t>(1, chunk_bytes /
        max<size_t>(sizeof(AxisObservationRecord), 1));
    chunk.reserve(min<size_t>(capacity, 10000000));
    char buffer[1<<20];
    unsigned long long line_no = 0;
    try {
        while (gzgets(input, buffer, sizeof(buffer))){
            ++line_no;
            string line(buffer);
            line.erase(remove(line.begin(), line.end(), '\n'), line.end());
            line.erase(remove(line.begin(), line.end(), '\r'), line.end());
            if (line.empty()) continue;
            vector<string> fields = split_tsv_strict(line);
            if (fields.empty()) continue;
            unsigned long barcode = 0;
            try { barcode = strict_barcode_number(fields[0], path + ": line " + to_string(line_no)); }
            catch (const exception&) { continue; }
            if (pairs.find(barcode) == pairs.end()) continue;
            AxisObservationRecord record = axis_parse_observation(
                fields, path, line_no);
            ++audit.target_rows;
            ++rows_by_barcode[barcode];
            chunk.push_back(record);
            if (chunk.size() >= capacity){
                runs.push_back(axis_spill_run(chunk, temp_path, runs.size()));
            }
        }
        if (gzclose(input) != Z_OK)
            throw runtime_error("failed closing candidate-axis observations after first pass: " + path);
        input = NULL;
    } catch (...) {
        if (input) gzclose(input);
        throw;
    }
    if (!chunk.empty() || runs.empty())
        runs.push_back(axis_spill_run(chunk, temp_path, runs.size()));
    audit.spill_runs = runs.size();
    audit.target_barcodes = rows_by_barcode.size();
    for (const auto& item : rows_by_barcode)
        audit.largest_barcode_rows = max(audit.largest_barcode_rows, item.second);
}

static void axis_merge_first_pass_runs(
        const vector<string>& runs,
        const string& merged_path,
        const string& temp_path,
        unsigned long long key_chunk_bytes,
        AxisResourceAudit& audit,
        vector<uint64_t>& selected_site_keys){
    vector<ifstream> owned(runs.size());
    priority_queue<AxisRunCursor,vector<AxisRunCursor>,AxisRunCursorGreater> heap;
    for (size_t i = 0; i < runs.size(); ++i){
        owned[i].open(runs[i].c_str(), ios::binary);
        if (!owned[i]) throw runtime_error("could not open candidate-axis spill run: " + runs[i]);
        AxisObservationRecord record;
        if (axis_read_binary_record(owned[i], record, runs[i])){
            AxisRunCursor cursor; cursor.record = record; cursor.run = i;
            heap.push(cursor);
        }
    }
    ofstream merged(merged_path.c_str(), ios::binary);
    if (!merged) throw runtime_error("could not create merged candidate-axis first-pass file: " + merged_path);
    bool have = false;
    AxisObservationRecord current;
    AxisKahan ref, alt;
    unsigned long long copies = 0;
    vector<uint64_t> key_chunk;
    const size_t key_capacity = max<size_t>(1,
        key_chunk_bytes / sizeof(uint64_t));
    key_chunk.reserve(key_capacity);
    vector<string> key_runs;
    auto flush = [&](){
        if (!have) return;
        current.ref = (double)ref.value;
        current.alt = (double)alt.value;
        axis_write_binary_record(merged, current, merged_path);
        ++audit.unique_cell_sites;
        audit.duplicate_rows += copies > 0 ? copies - 1 : 0;
        if (ref.value + alt.value <= 0.0L) ++audit.zero_depth_sites;
        key_chunk.push_back(site_key(current.tid, current.pos));
        if (key_chunk.size() >= key_capacity)
            key_runs.push_back(axis_spill_key_run(
                key_chunk, temp_path, key_runs.size()));
    };
    while (!heap.empty()){
        AxisRunCursor cursor = heap.top(); heap.pop();
        const AxisObservationRecord& record = cursor.record;
        if (!have || record.barcode != current.barcode ||
                record.tid != current.tid || record.pos != current.pos){
            flush();
            current = record; ref = AxisKahan(); alt = AxisKahan();
            copies = 0; have = true;
        }
        ref.add(record.ref); alt.add(record.alt); ++copies;
        AxisObservationRecord next;
        if (axis_read_binary_record(owned[cursor.run], next, runs[cursor.run])){
            AxisRunCursor following; following.record = next;
            following.run = cursor.run; heap.push(following);
        }
    }
    flush();
    merged.close();
    if (!merged) throw runtime_error("failed closing merged candidate-axis first-pass file: " + merged_path);
    if (!key_chunk.empty() || key_runs.empty())
        key_runs.push_back(axis_spill_key_run(
            key_chunk, temp_path, key_runs.size()));
    vector<ifstream> key_inputs(key_runs.size());
    priority_queue<AxisKeyCursor,vector<AxisKeyCursor>,AxisKeyCursorGreater>
        key_heap;
    for (size_t i = 0; i < key_runs.size(); ++i){
        key_inputs[i].open(key_runs[i].c_str(), ios::binary);
        if (!key_inputs[i]) throw runtime_error(
            "could not open candidate-axis selected-key spill run: " +
            key_runs[i]);
        uint64_t key = 0;
        if (axis_read_key(key_inputs[i], key, key_runs[i])){
            AxisKeyCursor cursor; cursor.key = key; cursor.run = i;
            key_heap.push(cursor);
        }
    }
    bool have_key = false;
    uint64_t prior_key = 0;
    unsigned long long unique_key_count = 0;
    const string unique_key_path = temp_path + "/selected_keys_unique.bin";
    ofstream unique_key_output(unique_key_path.c_str(), ios::binary);
    if (!unique_key_output) throw runtime_error(
        "could not create candidate-axis unique selected-key file: " +
        unique_key_path);
    while (!key_heap.empty()){
        const AxisKeyCursor cursor = key_heap.top(); key_heap.pop();
        if (!have_key || cursor.key != prior_key){
            unique_key_output.write(
                reinterpret_cast<const char*>(&cursor.key),
                sizeof(cursor.key));
            if (!unique_key_output) throw runtime_error(
                "failed writing candidate-axis unique selected-key file: " +
                unique_key_path);
            ++unique_key_count;
            prior_key = cursor.key;
            have_key = true;
        }
        uint64_t next = 0;
        if (axis_read_key(key_inputs[cursor.run], next, key_runs[cursor.run])){
            AxisKeyCursor following; following.key = next;
            following.run = cursor.run; key_heap.push(following);
        }
    }
    unique_key_output.close();
    if (!unique_key_output) throw runtime_error(
        "failed closing candidate-axis unique selected-key file: " +
        unique_key_path);
    for (size_t i = 0; i < key_inputs.size(); ++i){
        key_inputs[i].close();
        unlink(key_runs[i].c_str());
    }
    audit.unique_site_keys = unique_key_count;
    audit.selected_key_bytes = unique_key_count * sizeof(uint64_t);
    selected_site_keys.resize(unique_key_count);
    ifstream unique_key_input(unique_key_path.c_str(), ios::binary);
    if (!unique_key_input) throw runtime_error(
        "could not reopen candidate-axis unique selected-key file: " +
        unique_key_path);
    if (unique_key_count > 0)
        unique_key_input.read(
            reinterpret_cast<char*>(&selected_site_keys[0]),
            unique_key_count * sizeof(uint64_t));
    if (!unique_key_input && unique_key_count > 0)
        throw runtime_error(
            "failed reading candidate-axis unique selected-key file: " +
            unique_key_path);
    unique_key_input.close();
    unlink(unique_key_path.c_str());
}

static int8_t axis_parse_gt(
        const string& raw, const string& context){
    const string value = lowercase(trim(raw));
    if (value.empty() || value == "." || value == "na") return -1;
    const long long parsed = strict_ll(raw, context);
    if (parsed < -1 || parsed > 2)
        throw runtime_error(context + ": hard GT must be -1, 0, 1, or 2");
    return (int8_t)parsed;
}

static vector<AxisSiteDefinition> axis_load_site_definitions(
        const string& path,
        const vector<uint64_t>& selected_keys,
        int n_samples,
        const vector<int>& donors,
        AxisResourceAudit& audit){
    vector<AxisSiteDefinition> sites(selected_keys.size());
    for (size_t i = 0; i < selected_keys.size(); ++i) sites[i].key = selected_keys[i];
    gzFile input = gzopen(path.c_str(), "rb");
    if (!input) throw runtime_error("could not open candidate-axis site definitions: " + path);
    char buffer[1<<20];
    unsigned long long line_no = 0;
    try {
        while (gzgets(input, buffer, sizeof(buffer))){
            ++line_no;
            string line(buffer);
            line.erase(remove(line.begin(), line.end(), '\n'), line.end());
            line.erase(remove(line.begin(), line.end(), '\r'), line.end());
            if (line.empty()) continue;
            vector<string> fields = split_tsv_strict(line);
            if (fields.size() < 3)
                throw runtime_error(path + ": malformed site row at line " + to_string(line_no));
            const long long tid_ll = strict_ll(fields[0], path + ": line " + to_string(line_no) + " tid");
            const long long pos_ll = strict_ll(fields[2], path + ": line " + to_string(line_no) + " position");
            if (tid_ll < 0 || tid_ll > INT_MAX || pos_ll < 0 || pos_ll > INT_MAX)
                throw runtime_error(path + ": tid and position must be nonnegative 32-bit integers at line " + to_string(line_no));
            const uint64_t key = site_key((int)tid_ll, (int)pos_ll);
            auto selected = lower_bound(selected_keys.begin(), selected_keys.end(), key);
            if (selected == selected_keys.end() || *selected != key) continue;
            if ((int)fields.size() != 5 + n_samples)
                throw runtime_error(path + ": selected site row has wrong field count at line " + to_string(line_no));
            const size_t index = selected - selected_keys.begin();
            AxisSiteDefinition definition;
            definition.key = key; definition.tid = (int32_t)tid_ll;
            definition.pos = (int32_t)pos_ll; definition.contig = trim(fields[1]);
            definition.ref_allele = trim(fields[3]);
            definition.alt_allele = trim(fields[4]);
            const string normalized = normalized_contig(definition.contig);
            definition.mitochondrial = normalized == "m" || normalized == "mt";
            definition.found = true;
            definition.genotype.reserve(donors.size());
            for (int donor : donors){
                if (donor < 0 || donor >= n_samples)
                    throw runtime_error("candidate-axis donor index outside sample vector");
                definition.genotype.push_back(axis_parse_gt(
                    fields[5 + donor], path + ": line " + to_string(line_no) +
                    " sample_index=" + to_string(donor)));
            }
            if (sites[index].found){
                const AxisSiteDefinition& prior = sites[index];
                if (prior.tid != definition.tid || prior.pos != definition.pos ||
                        prior.contig != definition.contig ||
                        prior.ref_allele != definition.ref_allele ||
                        prior.alt_allele != definition.alt_allele ||
                        prior.genotype != definition.genotype)
                    throw runtime_error(path + ": inconsistent duplicate selected site definition for tid=" +
                        to_string(definition.tid) + " pos=" + to_string(definition.pos));
            } else sites[index] = definition;
        }
        if (gzclose(input) != Z_OK)
            throw runtime_error("failed closing candidate-axis site definitions: " + path);
        input = NULL;
    } catch (...) {
        if (input) gzclose(input);
        throw;
    }
    unsigned long long bytes = sites.capacity() * sizeof(AxisSiteDefinition);
    for (const AxisSiteDefinition& site : sites){
        bytes += site.contig.capacity() + site.ref_allele.capacity() +
            site.alt_allele.capacity() + site.genotype.capacity() * sizeof(int8_t);
    }
    audit.compact_site_definition_bytes = bytes;
    return sites;
}

static void axis_classify_first_pass(
        const string& merged_path,
        const vector<AxisSiteDefinition>& sites,
        const vector<CandidateHypothesis>& candidates,
        const unordered_map<unsigned long,CandidateAxisPair>& pairs,
        const unordered_map<int,size_t>& donor_slot,
        AxisResourceAudit& audit){
    ifstream input(merged_path.c_str(), ios::binary);
    if (!input) throw runtime_error("could not open merged candidate-axis first-pass file: " + merged_path);
    AxisObservationRecord record;
    while (axis_read_binary_record(input, record, merged_path)){
        const uint64_t key = site_key(record.tid, record.pos);
        auto site = lower_bound(sites.begin(), sites.end(), key,
            [](const AxisSiteDefinition& item, uint64_t value){ return item.key < value; });
        if (site == sites.end() || site->key != key || !site->found){
            ++audit.missing_site_definitions; continue;
        }
        if (site->mitochondrial){ ++audit.mitochondrial_sites; continue; }
        if ((long double)record.ref + record.alt <= 0.0L) continue;
        auto pair = pairs.find(record.barcode);
        if (pair == pairs.end())
            throw runtime_error("internal candidate-axis barcode lookup failure during resource classification");
        long double a = 0.0L, b = 0.0L;
        const bool has_a = axis_expected(candidates[pair->second.original].identity,
            *site, donor_slot, a);
        const bool has_b = axis_expected(candidates[pair->second.proposed].identity,
            *site, donor_slot, b);
        if (!has_a && !has_b) ++audit.missing_both_candidates;
        else if (!has_a) ++audit.missing_candidate_a_only;
        else if (!has_b) ++audit.missing_candidate_b_only;
        else ++audit.common_nuclear_sites;
    }
}

static vector<unsigned long long> axis_assign_buckets(
        const unordered_map<unsigned long,unsigned long long>& rows_by_barcode,
        unsigned long long bucket_target_bytes,
        unordered_map<unsigned long,size_t>& assignment,
        AxisResourceAudit& audit){
    vector<pair<unsigned long,unsigned long long>> ordered(rows_by_barcode.begin(),
        rows_by_barcode.end());
    sort(ordered.begin(), ordered.end(), [](const pair<unsigned long,unsigned long long>& a,
                                            const pair<unsigned long,unsigned long long>& b){
        if (a.second != b.second) return a.second > b.second;
        return a.first < b.first;
    });
    if (ordered.empty()){
        audit.bucket_count = 0;
        return vector<unsigned long long>();
    }
    const unsigned long long row_bytes = sizeof(AxisObservationRecord);
    unsigned long long total_rows = 0;
    for (const auto& item : ordered) total_rows += item.second;
    size_t bucket_count = max<unsigned long long>(1,
        (total_rows * row_bytes + bucket_target_bytes - 1) / bucket_target_bytes);
    bucket_count = min(bucket_count, ordered.size());
    vector<unsigned long long> loads;
    while (true){
        loads.assign(bucket_count, 0);
        assignment.clear();
        for (const auto& item : ordered){
            size_t best = 0;
            for (size_t i = 1; i < loads.size(); ++i)
                if (make_pair(loads[i],i) < make_pair(loads[best],best)) best = i;
            assignment[item.first] = best;
            loads[best] += item.second;
        }
        const unsigned long long largest = *max_element(loads.begin(), loads.end());
        if (largest * row_bytes <= bucket_target_bytes || bucket_count == ordered.size()) break;
        ++bucket_count;
    }
    audit.bucket_count = loads.size();
    audit.estimated_largest_bucket_rows = *max_element(loads.begin(), loads.end());
    return loads;
}

static vector<string> axis_second_observation_pass(
        const string& path,
        const unordered_map<unsigned long,size_t>& assignment,
        size_t bucket_count,
        const string& temp_path){
    vector<string> paths(bucket_count);
    vector<ofstream> outputs(bucket_count);
    for (size_t i = 0; i < bucket_count; ++i){
        paths[i] = temp_path + "/bucket_" + to_string(i) + ".bin";
        outputs[i].open(paths[i].c_str(), ios::binary);
        if (!outputs[i]) throw runtime_error("could not create candidate-axis bucket: " + paths[i]);
    }
    gzFile input = gzopen(path.c_str(), "rb");
    if (!input) throw runtime_error("could not reopen candidate-axis observations: " + path);
    char buffer[1<<20];
    unsigned long long line_no = 0;
    try {
        while (gzgets(input, buffer, sizeof(buffer))){
            ++line_no;
            string line(buffer);
            line.erase(remove(line.begin(), line.end(), '\n'), line.end());
            line.erase(remove(line.begin(), line.end(), '\r'), line.end());
            if (line.empty()) continue;
            vector<string> fields = split_tsv_strict(line);
            if (fields.empty()) continue;
            unsigned long barcode = 0;
            try { barcode = strict_barcode_number(fields[0], path + ": line " + to_string(line_no)); }
            catch (const exception&) { continue; }
            auto found = assignment.find(barcode);
            if (found == assignment.end()) continue;
            const AxisObservationRecord record = axis_parse_observation(fields, path, line_no);
            axis_write_binary_record(outputs[found->second], record, paths[found->second]);
        }
        if (gzclose(input) != Z_OK)
            throw runtime_error("failed closing candidate-axis observations after second pass: " + path);
        input = NULL;
    } catch (...) {
        if (input) gzclose(input);
        throw;
    }
    for (size_t i = 0; i < outputs.size(); ++i){
        outputs[i].close();
        if (!outputs[i]) throw runtime_error("failed closing candidate-axis bucket: " + paths[i]);
    }
    return paths;
}

static vector<string> axis_output_header(){
    return {
        "schema_version","library","barcode","score_pair_id",
        "supported_event_key","selected_supported_event_id",
        "selected_supported_event_proposal","source_reconciliation_event_id",
        "source_reconciliation_proposed_identity","score_population_scope",
        "population_votes_in_authoritative_event",
        "original_allowed_demux_assignment","reconciliation_nominated_swap",
        "candidate_b_fixed_identity","candidate_a","candidate_b",
        "candidate_a_role","candidate_b_role","candidate_a_origin",
        "candidate_b_origin","score_pair_source","pair_construction_mode",
        "score_scope_contract","candidate_axis_status","candidate_axis_direction",
        "candidate_axis_position_raw","candidate_axis_distance_from_midpoint_absolute",
        "candidate_axis_segment","candidate_axis_evidence_basis",
        "candidate_axis_evidence_basis_interpretation","candidate_axis_numerator",
        "candidate_axis_design_mass","candidate_axis_brier_margin_original_minus_proposal",
        "candidate_axis_common_weight_sum","candidate_axis_discriminating_weight_sum",
        "observed_design_candidate_separation_rms_out_of_100",
        "observed_design_candidate_similarity_complement_out_of_100",
        "n_common_observed_nuclear_sites","n_unique_merged_target_cell_sites",
        "n_candidate_axis_discriminating_sites","n_primary_evidence_units",
        "n_corrected_molecules","total_common_observed_evidence_depth",
        "discriminating_evidence_depth","n_duplicate_observation_rows_merged",
        "n_sites_excluded_mitochondrial","n_sites_excluded_missing_candidate_a_only",
        "n_sites_excluded_missing_candidate_b_only",
        "n_sites_excluded_missing_both_candidates",
        "n_sites_excluded_missing_site_definition",
        "n_sites_excluded_nonpositive_observation",
        "candidate_axis_position_without_top_primary_unit",
        "candidate_axis_position_without_top_five_primary_units",
        "candidate_axis_direction_without_top_primary_unit",
        "candidate_axis_direction_without_top_five_primary_units",
        "candidate_axis_direction_preserved_without_top_primary_unit",
        "candidate_axis_direction_preserved_without_top_five_primary_units",
        "n_primary_units_removed_top_one","n_primary_units_removed_top_five",
        "top_primary_unit_removal_status","top_five_primary_units_removal_status",
        "maximum_primary_unit_absolute_brier_margin_fraction",
        "top_five_primary_units_absolute_brier_margin_fraction",
        "primary_unit_brier_margin_concentration_status",
        "sum_absolute_primary_unit_axis_numerator_contributions",
        "sum_absolute_primary_unit_brier_margins","top_primary_unit_id",
        "top_five_primary_unit_ids","candidate_axis_fold_group_basis",
        "candidate_axis_fold_definition_version","candidate_axis_fold_count",
        "n_candidate_axis_folds_evaluable",
        "minimum_leave_one_fold_out_candidate_axis_position",
        "median_leave_one_fold_out_candidate_axis_position",
        "maximum_leave_one_fold_out_candidate_axis_position",
        "fraction_leave_one_fold_out_positions_proposal_side",
        "candidate_axis_direction_preserved_all_evaluable_folds",
        "candidate_axis_fold_direction_stability_status",
        "site_log_likelihood_candidate_a","site_log_likelihood_candidate_b",
        "site_delta_log_likelihood_a_minus_b",
        "v6_3_compatible_site_delta_log_likelihood_a_minus_b",
        "v6_3_compatible_discriminating_evidence_depth",
        "v6_3_compatible_site_candidate_a_probability_pct",
        "v6_3_compatible_site_candidate_b_probability_pct",
        "comparison_status_legacy","probability_basis_legacy",
        "molecule_evidence_status_legacy","n_independent_molecules_legacy",
        "site_candidate_a_probability_pct","site_candidate_b_probability_pct",
        "legacy_saturated_preferred_probability_pct",
        "legacy_saturated_probability_interpretation","n_sites_favor_candidate_a",
        "n_sites_favor_candidate_b","n_sites_tied_legacy_log_score",
        "candidate_a_residual_mismatch_legacy","candidate_b_residual_mismatch_legacy",
        "absolute_fit_status_legacy","probability_without_top_site_pct_legacy",
        "probability_without_top_five_sites_pct_legacy",
        "minimum_error_sensitivity_probability_pct_legacy",
        "error_sensitivity_stable_legacy","error_ref","error_alt",
        "min_evidence_legacy","poor_fit_residual_threshold_legacy","resamples",
        "candidate_axis_formula_version","candidate_prediction_transform_version",
        "numerical_tolerance_version","long_double_mantissa_digits",
        "long_double_epsilon",
        "candidate_axis_bucket_count","candidate_axis_target_observation_rows",
        "candidate_axis_unique_selected_site_keys","candidate_axis_peak_bucket_rows",
        "warnings","candidate_a_observed_brier_mean",
        "candidate_a_expected_sampling_brier_mean",
        "candidate_a_excess_brier_mean","candidate_b_observed_brier_mean",
        "candidate_b_expected_sampling_brier_mean",
        "candidate_b_excess_brier_mean","raw_residual_threshold_flag"
    };
}

static vector<string> axis_output_row(
        const CandidateHypothesis& a,
        const CandidateHypothesis& b,
        const AxisCellResult& result,
        const AxisResourceAudit& audit,
        long double e_ref,
        long double e_alt,
        long min_evidence,
        long double poor_fit_residual){
    const long double preferred_probability = !isfinite(result.ll_delta) ? NAN :
        max(result.legacy_probability_a, result.legacy_probability_b);
    vector<string> warnings = result.warnings;
    if (result.numeric.position < 0.0L || result.numeric.position > 100.0L)
        warnings.push_back("RAW_POSITION_OUTSIDE_FIXED_CANDIDATE_SEGMENT");
    if (result.fold_status == "DIRECTION_CHANGED_OR_TIED")
        warnings.push_back("FOLD_DIRECTION_UNSTABLE");
    if (result.preserve_top == "FALSE" || result.preserve_five == "FALSE")
        warnings.push_back("PRIMARY_UNIT_REMOVAL_DIRECTION_UNSTABLE");
    string warning_text = warnings.empty() ? "NONE" : join_flags(warnings);
    return {
        AXIS_SCHEMA,a.library,a.barcode,a.score_pair_id,a.supported_event_key,
        a.selected_supported_event_id,a.selected_supported_event_proposal,
        a.source_reconciliation_event_id,a.source_reconciliation_proposed_identity,
        a.score_population_scope,a.population_votes_in_authoritative_event,
        a.original_demux_assignment,a.reconciliation_nominated_swap,
        a.candidate_b_fixed_identity,a.donor_genotype,b.donor_genotype,
        a.score_pair_role,b.score_pair_role,a.candidate_origin,b.candidate_origin,
        a.score_pair_source,a.pair_construction_mode,a.score_scope_contract,
        result.numeric.status,result.numeric.direction,axis_fmt(result.numeric.position),
        axis_fmt(isfinite(result.numeric.position) ? fabsl(result.numeric.position - 50.0L) : NAN),
        result.numeric.segment,AXIS_BASIS,AXIS_BASIS_INTERPRETATION,
        axis_fmt(result.sums.n),axis_fmt(result.sums.d),axis_fmt(result.numeric.margin),
        axis_fmt(result.sums.w),to_string(result.discriminating_sites),
        axis_fmt(result.separation),axis_fmt(result.similarity),
        to_string(result.common_sites),to_string(result.n_unique_merged),
        to_string(result.discriminating_sites),to_string(result.common_sites),"NA",
        axis_fmt(result.total_common_depth),axis_fmt(result.discriminating_depth),
        to_string(result.n_duplicate_rows),to_string(result.excluded_mito),
        to_string(result.excluded_missing_a),to_string(result.excluded_missing_b),
        to_string(result.excluded_missing_both),to_string(result.excluded_missing_definition),
        to_string(result.excluded_nonpositive),axis_fmt(result.without_top_position),
        axis_fmt(result.without_five_position),result.without_top_direction,
        result.without_five_direction,result.preserve_top,result.preserve_five,
        to_string(result.removed_top),to_string(result.removed_five),
        result.removal_top_status,result.removal_five_status,
        axis_fmt(result.maximum_margin_fraction),axis_fmt(result.top_five_margin_fraction),
        result.concentration_status,axis_fmt(result.sums.sum_abs_n),
        axis_fmt(result.sums.sum_abs_m),result.top_unit_id,result.top_five_unit_ids,
        AXIS_FOLD_BASIS,AXIS_FOLD_VERSION,to_string(result.fold_count),
        to_string(result.folds_evaluable),axis_fmt(result.fold_min),
        axis_fmt(result.fold_median),axis_fmt(result.fold_max),
        axis_fmt(result.fold_proposal_fraction),result.folds_preserved,result.fold_status,
        axis_fmt(result.ll_a),axis_fmt(result.ll_b),axis_fmt(result.ll_delta),
        axis_fmt6((double)result.ll_delta),axis_fmt6((double)result.discriminating_depth),
        axis_fmt6((double)result.legacy_probability_a),
        axis_fmt6((double)result.legacy_probability_b),result.comparison_status,
        "nuclear_site_likelihood_equal_priors","MOLECULE_SIDECAR_UNAVAILABLE","NA",
        axis_fmt(result.legacy_probability_a),axis_fmt(result.legacy_probability_b),
        axis_fmt(preferred_probability),
        "LEGACY_AUDIT_ONLY_TOTAL_SITE_LIKELIHOOD_PERCENTAGE_NOT_CORRECTNESS",
        to_string(result.sites_favor_a),to_string(result.sites_favor_b),
        to_string(result.sites_tied),axis_fmt(result.residual_a),axis_fmt(result.residual_b),
        result.absolute_fit_status,axis_fmt(result.legacy_without_top),
        axis_fmt(result.legacy_without_top_five),axis_fmt(result.minimum_error_probability),
        result.error_stable,axis_fmt(e_ref),axis_fmt(e_alt),to_string(min_evidence),
        axis_fmt(poor_fit_residual),"0",AXIS_FORMULA,AXIS_PREDICTION_TRANSFORM,
        AXIS_TOLERANCE_VERSION,to_string(LDBL_MANT_DIG),
        axis_fmt(numeric_limits<long double>::epsilon()),
        to_string(audit.bucket_count),to_string(audit.target_rows),
        to_string(audit.unique_site_keys),to_string(audit.observed_peak_bucket_rows),
        warning_text,axis_fmt(result.candidate_a_observed_brier_mean),
        axis_fmt(result.candidate_a_expected_sampling_brier_mean),
        axis_fmt(result.candidate_a_excess_brier_mean),
        axis_fmt(result.candidate_b_observed_brier_mean),
        axis_fmt(result.candidate_b_expected_sampling_brier_mean),
        axis_fmt(result.candidate_b_excess_brier_mean),
        result.raw_residual_threshold_flag
    };
}

static void axis_gzwrite(gzFile output, const vector<string>& fields){
    string line;
    for (size_t i = 0; i < fields.size(); ++i){
        if (i) line.push_back('\t');
        const string& value = fields[i];
        if (value.empty()) line += "NA";
        else {
            const string lower = lowercase(value);
            if (lower == "nan" || lower == "inf" || lower == "+inf" ||
                    lower == "-inf")
                throw runtime_error("candidate-axis attempted to emit a nonfinite TSV token");
            line += value;
        }
    }
    line.push_back('\n');
    if (gzwrite(output, line.data(), (unsigned int)line.size()) == 0)
        throw runtime_error("failed writing candidate-axis output");
}

static unordered_map<unsigned long,AxisCellResult> axis_process_buckets(
        const vector<string>& bucket_paths,
        const vector<CandidateHypothesis>& candidates,
        const unordered_map<unsigned long,CandidateAxisPair>& pairs,
        const vector<AxisSiteDefinition>& sites,
        const unordered_map<int,size_t>& donor_slot,
        long double e_ref,
        long double e_alt,
        long min_evidence,
        long double poor_fit_residual,
        AxisResourceAudit& audit){
    unordered_map<unsigned long,AxisCellResult> results;
    results.reserve(pairs.size());
    for (const string& path : bucket_paths){
        struct stat info;
        if (stat(path.c_str(), &info) != 0)
            throw runtime_error("could not stat candidate-axis bucket: " + path);
        const size_t n = (size_t)info.st_size / sizeof(AxisObservationRecord);
        audit.observed_peak_bucket_rows = max<unsigned long long>(
            audit.observed_peak_bucket_rows, n);
        vector<AxisObservationRecord> records(n);
        ifstream input(path.c_str(), ios::binary);
        if (!input) throw runtime_error("could not open candidate-axis bucket: " + path);
        if (n){
            input.read(reinterpret_cast<char*>(records.data()),
                n * sizeof(AxisObservationRecord));
            if (!input) throw runtime_error("failed reading candidate-axis bucket: " + path);
        }
        sort(records.begin(), records.end(), axis_observation_less);
        size_t begin = 0;
        while (begin < records.size()){
            const unsigned long barcode = records[begin].barcode;
            size_t end = begin;
            while (end < records.size() && records[end].barcode == barcode) ++end;
            vector<AxisObservationRecord> merged;
            vector<long> duplicate_counts;
            size_t cursor = begin;
            while (cursor < end){
                const size_t first = cursor;
                AxisObservationRecord aggregate = records[cursor];
                AxisKahan ref, alt;
                while (cursor < end && records[cursor].tid == aggregate.tid &&
                        records[cursor].pos == aggregate.pos){
                    ref.add(records[cursor].ref); alt.add(records[cursor].alt);
                    ++cursor;
                }
                aggregate.ref = (double)ref.value;
                aggregate.alt = (double)alt.value;
                merged.push_back(aggregate);
                duplicate_counts.push_back((long)(cursor - first - 1));
            }
            auto pair = pairs.find(barcode);
            if (pair == pairs.end())
                throw runtime_error("candidate-axis bucket contains a nontarget barcode");
            results[barcode] = evaluate_candidate_axis_cell(
                candidates[pair->second.original], candidates[pair->second.proposed],
                merged, duplicate_counts, sites, donor_slot, e_ref, e_alt,
                min_evidence, poor_fit_residual);
            begin = end;
        }
    }
    for (const auto& item : pairs){
        if (results.count(item.first) == 0){
            const vector<AxisObservationRecord> empty;
            const vector<long> empty_counts;
            results[item.first] = evaluate_candidate_axis_cell(
                candidates[item.second.original], candidates[item.second.proposed],
                empty, empty_counts, sites, donor_slot, e_ref, e_alt,
                min_evidence, poor_fit_residual);
        }
    }
    return results;
}

static void run_candidate_axis(
        const string& samples_path,
        const string& manifest_path,
        const string& sites_path,
        const string& observations_path,
        const string& output_path,
        const string& temp_root,
        const string& library,
        long double e_ref,
        long double e_alt,
        long min_evidence,
        long double poor_fit_residual){
    if (e_ref < 0.0L || e_ref > 1.0L || e_alt < 0.0L || e_alt > 1.0L ||
            e_ref + e_alt >= 1.0L)
        throw runtime_error("candidate-axis errors must each be in [0,1] and sum to less than one");
    if (min_evidence < 0)
        throw runtime_error("--min_evidence must be nonnegative in candidate-axis mode");
    if (poor_fit_residual < 0.0L || poor_fit_residual > 1.0L)
        throw runtime_error("--poor-fit-residual must be within [0,1]");
    vector<string> samples = load_samples(samples_path);
    unordered_map<string,int> sample2idx;
    for (int i = 0; i < (int)samples.size(); ++i){
        if (sample2idx.count(samples[i]))
            throw runtime_error("duplicate sample name in candidate-axis sample vector: " + samples[i]);
        if (trim(samples[i]).empty())
            throw runtime_error("blank sample name in candidate-axis sample vector");
        sample2idx[samples[i]] = i;
    }
    unordered_map<unsigned long,vector<size_t>> candidate_by_cell;
    vector<CandidateHypothesis> candidates = load_candidate_manifest(
        manifest_path, sample2idx, candidate_by_cell, true);
    unordered_map<unsigned long,CandidateAxisPair> pairs;
    pairs.reserve(candidate_by_cell.size());
    for (const auto& item : candidate_by_cell)
        pairs[item.first] = candidate_axis_pair(
            candidates, item.second, item.first, library);

    const string audit_path = axis_parent(output_path) + "/" + library +
        ".candidate_axis_resource_evidence_audit.tsv";
    AxisResourceAudit audit;
    if (pairs.empty()){
        axis_write_resource_audit(audit_path, library, temp_root, audit);
        const string temporary = output_path + ".tmp." + to_string((long long)getpid());
        gzFile output = gzopen(temporary.c_str(), "wb");
        if (!output) throw runtime_error("could not create candidate-axis output: " + temporary);
        try {
            axis_gzwrite(output, axis_output_header());
            if (gzclose(output) != Z_OK)
                throw runtime_error("failed closing header-only candidate-axis output");
            output = NULL;
        } catch (...) {
            if (output) gzclose(output);
            unlink(temporary.c_str());
            throw;
        }
        if (rename(temporary.c_str(), output_path.c_str()) != 0){
            unlink(temporary.c_str());
            throw runtime_error("failed publishing header-only candidate-axis output: " + output_path);
        }
        return;
    }

    vector<int> donors;
    for (const auto& item : pairs){
        const Identity* identities[] = {
            &candidates[item.second.original].identity,
            &candidates[item.second.proposed].identity
        };
        for (const Identity* identity : identities){
            donors.push_back(identity->a);
            if (identity->b >= 0) donors.push_back(identity->b);
        }
    }
    sort(donors.begin(), donors.end());
    donors.erase(unique(donors.begin(), donors.end()), donors.end());
    unordered_map<int,size_t> donor_slot;
    for (size_t i = 0; i < donors.size(); ++i) donor_slot[donors[i]] = i;

    AxisTempGuard temporary(temp_root);
    const unsigned long long first_pass_chunk_bytes = 512ULL * 1024ULL * 1024ULL;
    const unsigned long long bucket_target_bytes = 1024ULL * 1024ULL * 1024ULL;
    vector<string> runs;
    unordered_map<unsigned long,unsigned long long> rows_by_barcode;
    axis_first_observation_pass(observations_path, pairs, temporary.path(),
        first_pass_chunk_bytes, audit, rows_by_barcode, runs);
    const string merged_path = temporary.path() + "/merged_first_pass.bin";
    vector<uint64_t> selected_keys;
    axis_merge_first_pass_runs(
        runs, merged_path, temporary.path(), first_pass_chunk_bytes,
        audit, selected_keys);
    for (const string& run : runs) unlink(run.c_str());
    vector<AxisSiteDefinition> sites = axis_load_site_definitions(
        sites_path, selected_keys, (int)samples.size(), donors, audit);
    axis_classify_first_pass(merged_path, sites, candidates, pairs, donor_slot, audit);
    const unsigned long long classified = audit.missing_site_definitions +
        audit.mitochondrial_sites + audit.zero_depth_sites +
        audit.missing_both_candidates + audit.missing_candidate_a_only +
        audit.missing_candidate_b_only + audit.common_nuclear_sites;
    if (classified != audit.unique_cell_sites)
        throw runtime_error("global candidate-axis first-pass site-accounting identity failed");

    unordered_map<unsigned long,size_t> bucket_assignment;
    vector<unsigned long long> bucket_loads = axis_assign_buckets(
        rows_by_barcode, bucket_target_bytes, bucket_assignment, audit);
    const unsigned long long largest_bucket_rows = bucket_loads.empty() ? 0 :
        *max_element(bucket_loads.begin(), bucket_loads.end());
    audit.observed_peak_bucket_rows = largest_bucket_rows;
    axis_write_resource_audit(audit_path, library, temp_root, audit);

    vector<string> buckets = axis_second_observation_pass(
        observations_path, bucket_assignment, bucket_loads.size(), temporary.path());
    unordered_map<unsigned long,AxisCellResult> results = axis_process_buckets(
        buckets, candidates, pairs, sites, donor_slot, e_ref, e_alt,
        min_evidence, poor_fit_residual, audit);
    for (const auto& item : pairs)
        if (results.find(item.first) == results.end())
            results[item.first] = AxisCellResult();
    if (audit.observed_peak_bucket_rows != largest_bucket_rows)
        throw runtime_error("candidate-axis observed bucket size did not match the first-pass exact assignment");

    const string output_temporary = output_path + ".tmp." +
        to_string((long long)getpid());
    gzFile output = gzopen(output_temporary.c_str(), "wb");
    if (!output) throw runtime_error("could not create candidate-axis output: " + output_temporary);
    try {
        axis_gzwrite(output, axis_output_header());
        vector<unsigned long> order;
        for (const auto& item : pairs) order.push_back(item.first);
        sort(order.begin(), order.end(), [&](unsigned long left, unsigned long right){
            const CandidateHypothesis& a = candidates[pairs.at(left).original];
            const CandidateHypothesis& b = candidates[pairs.at(right).original];
            if (a.library != b.library) return a.library < b.library;
            if (a.barcode != b.barcode) return a.barcode < b.barcode;
            return a.score_pair_id < b.score_pair_id;
        });
        const vector<string> header = axis_output_header();
        for (unsigned long barcode : order){
            const CandidateAxisPair& pair = pairs.at(barcode);
            const vector<string> row = axis_output_row(
                candidates[pair.original], candidates[pair.proposed],
                results.at(barcode), audit, e_ref, e_alt,
                min_evidence, poor_fit_residual);
            if (row.size() != header.size())
                throw runtime_error(
                    "candidate-axis internal output schema/value count mismatch");
            axis_gzwrite(output, row);
        }
        if (gzclose(output) != Z_OK)
            throw runtime_error("failed closing candidate-axis output: " + output_temporary);
        output = NULL;
    } catch (...) {
        if (output) gzclose(output);
        unlink(output_temporary.c_str());
        throw;
    }
    if (rename(output_temporary.c_str(), output_path.c_str()) != 0){
        unlink(output_temporary.c_str());
        throw runtime_error("failed publishing candidate-axis output: " + output_path);
    }
    fprintf(stderr, "Wrote candidate-axis scores for %lu fixed pairs to %s\n",
        (unsigned long)pairs.size(), output_path.c_str());
}

static int candidate_axis_self_test(){
    auto require = [](bool condition, const string& message){
        if (!condition) throw runtime_error("candidate-axis self-test failed: " + message);
    };
    auto one = [](long double a, long double b, long double y, int tid, int pos){
        AxisUnit unit; unit.a = a; unit.b = b; unit.y = y;
        unit.tid = tid; unit.pos = pos; unit.discriminating = a != b;
        const long double delta = b - a;
        unit.n = delta * (y - a); unit.d = delta * delta;
        unit.m = 2.0L * unit.n - unit.d;
        return unit;
    };
    vector<AxisUnit> units(1, one(0.2L,0.8L,0.2L,1,10));
    AxisNumericResult score = axis_numeric(axis_sum_units(units));
    require(fabsl(score.position) < 1e-15L, "y=a must give T=0");
    units[0] = one(0.2L,0.8L,0.8L,1,10);
    score = axis_numeric(axis_sum_units(units));
    require(fabsl(score.position - 100.0L) < 1e-15L, "y=b must give T=100");
    units[0] = one(0.2L,0.8L,0.5L,1,10);
    score = axis_numeric(axis_sum_units(units));
    require(fabsl(score.position - 50.0L) < 1e-15L && score.direction == "TIE",
        "midpoint must give T=50 and tie");
    const long double loss_a = (0.7L-0.2L)*(0.7L-0.2L);
    const long double loss_b = (0.7L-0.8L)*(0.7L-0.8L);
    units[0] = one(0.2L,0.8L,0.7L,1,10);
    score = axis_numeric(axis_sum_units(units));
    require((score.margin > 0.0L) == (loss_b < loss_a), "Brier margin direction");
    const AxisSamplingBrier half_zero = axis_sampling_brier(0.0L,0.5L,1.0L);
    const AxisSamplingBrier half_one = axis_sampling_brier(1.0L,0.5L,1.0L);
    require(half_zero.observed == 0.25L &&
        half_zero.expected_sampling == 0.25L && half_zero.excess == 0.0L &&
        half_one.observed == 0.25L &&
        half_one.expected_sampling == 0.25L && half_one.excess == 0.0L,
        "q=0.5 depth-one sampling-adjusted Brier identity");
    const AxisSamplingBrier quarter_zero =
        axis_sampling_brier(0.0L,0.25L,1.0L);
    const AxisSamplingBrier quarter_one =
        axis_sampling_brier(1.0L,0.25L,1.0L);
    const AxisSamplingBrier three_quarter_zero =
        axis_sampling_brier(0.0L,0.75L,1.0L);
    const AxisSamplingBrier three_quarter_one =
        axis_sampling_brier(1.0L,0.75L,1.0L);
    require(0.75L*quarter_zero.excess + 0.25L*quarter_one.excess == 0.0L,
        "q=0.25 depth-one expected excess Brier is zero");
    require(0.25L*three_quarter_zero.excess +
        0.75L*three_quarter_one.excess == 0.0L,
        "q=0.75 depth-one expected excess Brier is zero");
    const long double wrong_candidate_expected_excess =
        0.75L*three_quarter_zero.excess +
        0.25L*three_quarter_one.excess;
    require(wrong_candidate_expected_excess > 0.0L,
        "wrong candidate has positive expected excess Brier");
    const AxisNumericResult forward = score;
    vector<AxisUnit> swapped(1, one(0.8L,0.2L,0.7L,1,10));
    const AxisSums swapped_sums = axis_sum_units(swapped);
    const AxisNumericResult reverse_score = axis_numeric(swapped_sums);
    require(fabsl(reverse_score.position - (100.0L-forward.position)) < 1e-14L &&
        fabsl(reverse_score.margin + forward.margin) < 1e-14L &&
        fabsl(swapped_sums.d - axis_sum_units(units).d) < 1e-14L,
        "candidate swap symmetry");
    const AxisSamplingBrier forward_brier_a =
        axis_sampling_brier(0.7L,0.2L,5.0L);
    const AxisSamplingBrier forward_brier_b =
        axis_sampling_brier(0.7L,0.8L,5.0L);
    const AxisSamplingBrier reverse_brier_a =
        axis_sampling_brier(0.7L,0.8L,5.0L);
    const AxisSamplingBrier reverse_brier_b =
        axis_sampling_brier(0.7L,0.2L,5.0L);
    require(forward_brier_a.observed == reverse_brier_b.observed &&
        forward_brier_a.expected_sampling ==
            reverse_brier_b.expected_sampling &&
        forward_brier_a.excess == reverse_brier_b.excess &&
        forward_brier_b.observed == reverse_brier_a.observed &&
        forward_brier_b.expected_sampling ==
            reverse_brier_a.expected_sampling &&
        forward_brier_b.excess == reverse_brier_a.excess,
        "candidate swap exchanges sampling-adjusted Brier diagnostics");
    vector<AxisUnit> replicated(20, units[0]);
    AxisNumericResult repeated = axis_numeric(axis_sum_units(replicated));
    const AxisSums single_sums = axis_sum_units(units);
    const AxisSums repeated_sums = axis_sum_units(replicated);
    require(fabsl(repeated.position-forward.position) < 1e-14L &&
        repeated_sums.w == 20.0L * single_sums.w &&
        fabsl(repeated_sums.n - 20.0L*single_sums.n) < 1e-14L &&
        fabsl(repeated_sums.d - 20.0L*single_sums.d) < 1e-14L &&
        fabsl(repeated.margin - 20.0L*forward.margin) < 1e-14L,
        "replication must scale sufficient statistics and preserve geometry");
    const long double split_y = (3.0L+7.0L)/(10.0L+10.0L);
    const long double merged_y = 10.0L/20.0L;
    require(split_y == merged_y &&
        one(0.2L,0.8L,split_y,1,10).n == one(0.2L,0.8L,merged_y,1,10).n,
        "proportional duplicate observation rows must merge identically");
    AxisSums scaled = single_sums;
    scaled.w *= 7.0L; scaled.n *= 7.0L; scaled.d *= 7.0L;
    scaled.sum_abs_n *= 7.0L; scaled.sum_abs_m *= 7.0L;
    const AxisNumericResult scaled_score = axis_numeric(scaled);
    require(fabsl(scaled_score.position-forward.position) < 1e-14L &&
        fabsl(sqrtl(scaled.d/scaled.w)-sqrtl(single_sums.d/single_sums.w)) < 1e-14L,
        "uniform weight scaling must preserve position and separation");
    vector<AxisUnit> reordered;
    for (int i = 0; i < 100; ++i)
        reordered.push_back(one(0.1L,0.9L,(i%3)/2.0L,i,i+1));
    const AxisSums ordered_sums = axis_sum_units(reordered);
    reverse(reordered.begin(), reordered.end());
    const AxisSums reversed_sums = axis_sum_units(reordered);
    const long double reorder_tolerance = 128.0L *
        numeric_limits<long double>::epsilon() *
        max(max(fabsl(ordered_sums.n),ordered_sums.d),1.0L);
    require(fabsl(ordered_sums.n-reversed_sums.n) <= reorder_tolerance &&
        fabsl(ordered_sums.d-reversed_sums.d) <= reorder_tolerance,
        "row ordering must stay within the numerical tolerance");
    vector<AxisSamplingBrier> ordered_brier;
    for (int i = 0; i < 100; ++i){
        const long double expected = 0.1L + 0.008L*(long double)(i%100);
        const long double observed = (long double)(i%7) / 6.0L;
        ordered_brier.push_back(axis_sampling_brier(
            observed, expected, (long double)(i%11+1)));
    }
    auto excess_sum = [](const vector<AxisSamplingBrier>& values){
        AxisKahan sum;
        for (const AxisSamplingBrier& value : values) sum.add(value.excess);
        return sum.value;
    };
    const long double ordered_excess = excess_sum(ordered_brier);
    reverse(ordered_brier.begin(), ordered_brier.end());
    const long double reversed_excess = excess_sum(ordered_brier);
    require(fabsl(ordered_excess-reversed_excess) <=
        128.0L*numeric_limits<long double>::epsilon()*
        max(fabsl(ordered_excess),1.0L),
        "sampling-adjusted Brier accumulation is stable across site order");
    require(axis_raw_residual_threshold_flag(0.2L,0.4L,0.3L) ==
        "CANDIDATE_B_ONLY_ABOVE_LEGACY_THRESHOLD" &&
        axis_raw_residual_threshold_flag(0.4L,0.4L,0.3L) ==
        "BOTH_ABOVE_LEGACY_THRESHOLD" &&
        axis_raw_residual_threshold_flag(NAN,0.4L,0.3L) == "UNAVAILABLE",
        "legacy residual threshold flag remains neutral and explicit");
    const CandidateHypothesis output_a, output_b;
    const AxisCellResult output_result;
    const AxisResourceAudit output_audit;
    require(axis_output_header().size() == axis_output_row(
            output_a,output_b,output_result,output_audit,
            0.001L,0.001L,10,0.3L).size() &&
        string(AXIS_SCHEMA) ==
            "identity_candidate_axis_pair_score_v2_sampling_adjusted_fit_diagnostic",
        "sampling-adjusted raw schema header/value alignment");
    struct SyntheticAxisDesign {
        vector<long double> a, b;
        vector<int> depth;
    };
    vector<SyntheticAxisDesign> synthetic_designs;
    const vector<long double> low_a = {0.42L,0.37L,0.55L,0.48L,0.33L};
    const vector<long double> low_b = {0.58L,0.51L,0.69L,0.62L,0.47L};
    const vector<long double> high_a = {0.05L,0.20L,0.35L,0.10L,0.25L};
    const vector<long double> high_b = {0.95L,0.80L,0.65L,0.90L,0.75L};
    const vector<int> low_depth = {4,7,5,9,6};
    const vector<int> high_depth = {80,120,160,100,200};
    for (int separation = 0; separation < 2; ++separation){
        for (int depth_pattern = 0; depth_pattern < 2; ++depth_pattern){
            for (int orientation = 0; orientation < 2; ++orientation){
                SyntheticAxisDesign design;
                design.a = separation ? high_a : low_a;
                design.b = separation ? high_b : low_b;
                if (orientation) swap(design.a,design.b);
                design.depth = depth_pattern ? high_depth : low_depth;
                synthetic_designs.push_back(design);
            }
        }
    }
    require(synthetic_designs.size() == 8,
        "bounded synthetic design count");
    const vector<long double> generating_coordinates = {
        0.0L,25.0L,50.0L,75.0L,100.0L};
    for (const SyntheticAxisDesign& design : synthetic_designs){
        for (long double coordinate : generating_coordinates){
            vector<AxisUnit> noiseless;
            for (size_t i = 0; i < design.a.size(); ++i){
                const long double y = design.a[i] +
                    (coordinate/100.0L)*(design.b[i]-design.a[i]);
                noiseless.push_back(one(
                    design.a[i],design.b[i],y,(int)i,(int)i+1));
            }
            const AxisNumericResult recovered = axis_numeric(
                axis_sum_units(noiseless));
            const long double geometry_tolerance = 512.0L *
                max<size_t>(noiseless.size(),1) *
                numeric_limits<long double>::epsilon() * 100.0L;
            require(recovered.status == "AVAILABLE" &&
                fabsl(recovered.position-coordinate) <= geometry_tolerance,
                "known 0/25/50/75/100 coordinate recovery");
            vector<AxisUnit> reversed_candidates;
            for (size_t i = 0; i < noiseless.size(); ++i)
                reversed_candidates.push_back(one(
                    design.b[i],design.a[i],noiseless[i].y,
                    (int)i,(int)i+1));
            const AxisNumericResult reversed_candidates_score = axis_numeric(
                axis_sum_units(reversed_candidates));
            require(reversed_candidates_score.status == "AVAILABLE" &&
                fabsl(reversed_candidates_score.position-
                    (100.0L-coordinate)) <= geometry_tolerance,
                "multi-unit candidate reversal recovery");
            reverse(noiseless.begin(),noiseless.end());
            const AxisNumericResult reversed_input_score = axis_numeric(
                axis_sum_units(noiseless));
            require(reversed_input_score.status == "AVAILABLE" &&
                fabsl(reversed_input_score.position-coordinate) <=
                    geometry_tolerance &&
                (coordinate == 50.0L ||
                 reversed_input_score.direction == recovered.direction),
                "multi-unit input-order numerical and qualitative stability");
        }
    }
    auto sampling_summaries = [&](uint64_t seed){
        mt19937_64 rng(seed);
        vector<pair<long double,long double>> summaries;
        for (const SyntheticAxisDesign& design : synthetic_designs){
            vector<long double> design_biases, design_variances;
            for (long double coordinate : generating_coordinates){
                vector<long double> positions;
                for (int replicate_index = 0; replicate_index < 100;
                        ++replicate_index){
                    vector<AxisUnit> sampled;
                    for (size_t i = 0; i < design.a.size(); ++i){
                        const long double expected = design.a[i] +
                            (coordinate/100.0L)*
                            (design.b[i]-design.a[i]);
                        binomial_distribution<int> draw(
                            design.depth[i],(double)expected);
                        const long double observed =
                            (long double)draw(rng)/
                            (long double)design.depth[i];
                        sampled.push_back(one(
                            design.a[i],design.b[i],observed,
                            (int)i,(int)i+1));
                    }
                    const AxisNumericResult sampled_score = axis_numeric(
                        axis_sum_units(sampled));
                    require(sampled_score.status == "AVAILABLE",
                        "sampling design must remain informative");
                    positions.push_back(sampled_score.position);
                }
                const long double mean = accumulate(
                    positions.begin(),positions.end(),0.0L)/
                    (long double)positions.size();
                long double squared = 0.0L;
                for (long double position : positions)
                    squared += (position-mean)*(position-mean);
                const long double variance = squared/
                    (long double)(positions.size()-1);
                const long double standard_error = sqrtl(
                    variance/(long double)positions.size());
                const long double numerical_tolerance = 512.0L *
                    design.a.size() *
                    numeric_limits<long double>::epsilon() *
                    max(fabsl(mean),100.0L);
                const long double bias = mean-coordinate;
                require(fabsl(bias) <=
                    4.0L*standard_error+numerical_tolerance,
                    "bounded sampling bias");
                design_biases.push_back(bias);
                design_variances.push_back(
                    standard_error*standard_error);
                summaries.push_back(make_pair(mean,standard_error));
            }
            const long double aggregate_bias = accumulate(
                design_biases.begin(),design_biases.end(),0.0L)/
                (long double)design_biases.size();
            const long double aggregate_standard_error = sqrtl(accumulate(
                design_variances.begin(),design_variances.end(),0.0L))/
                (long double)design_variances.size();
            const long double numerical_tolerance = 512.0L *
                design.a.size() *
                numeric_limits<long double>::epsilon() * 100.0L;
            require(fabsl(aggregate_bias) <=
                4.0L*aggregate_standard_error+numerical_tolerance &&
                aggregate_bias <=
                4.0L*aggregate_standard_error+numerical_tolerance,
                "symmetric-coordinate aggregate bias has no Candidate-B drift");
        }
        return summaries;
    };
    const vector<pair<long double,long double>> sampling_first =
        sampling_summaries(0x8a5cd789635d2dffULL);
    const vector<pair<long double,long double>> sampling_second =
        sampling_summaries(0x8a5cd789635d2dffULL);
    require(sampling_first == sampling_second,
        "fixed-seed sampling summaries must be deterministic");
    AxisUnit agreeing = one(0.4L,0.4L,0.9L,2,20);
    replicated.push_back(agreeing);
    AxisSums with_agreeing = axis_sum_units(replicated);
    require(with_agreeing.w == 21.0L && agreeing.n == 0.0L &&
        agreeing.d == 0.0L && agreeing.m == 0.0L,
        "candidate-agreeing sites contribute only W");
    vector<AxisUnit> identical(1, agreeing);
    const AxisNumericResult identical_score = axis_numeric(
        axis_sum_units(identical));
    require(identical_score.status ==
        "INSUFFICIENT_CANDIDATE_SEPARATION" &&
        identical_score.direction == "UNAVAILABLE",
        "identical candidates cannot produce directional evidence");
    vector<AxisUnit> effectively_identical(4,
        one(0.5L,0.5L+1e-20L,0.9L,3,30));
    const AxisNumericResult effectively_identical_score = axis_numeric(
        axis_sum_units(effectively_identical));
    require(effectively_identical_score.status ==
        "INSUFFICIENT_CANDIDATE_SEPARATION" &&
        effectively_identical_score.direction == "UNAVAILABLE",
        "effectively zero-separation candidates cannot produce strong evidence");
    require(axis_numeric(axis_sum_units(vector<AxisUnit>(1,
        one(0.2L,0.8L,-0.1L,1,1)))).position < 0.0L, "negative raw position retained");
    require(axis_numeric(axis_sum_units(vector<AxisUnit>(1,
        one(0.2L,0.8L,1.1L,1,1)))).position > 100.0L, "raw position above 100 retained");
    require(axis_sum_units(replicated).w == (long double)replicated.size(),
        "every primary site has exactly unit weight");
    vector<size_t> removal_order = {0};
    const AxisSums removed = axis_remove(single_sums, units, removal_order, 1);
    require(removed.w == single_sums.w-1.0L &&
        removed.n == single_sums.n-units[0].n &&
        removed.d == single_sums.d-units[0].d,
        "top-unit removal must subtract W, N, and D");
    AxisCellResult influence;
    influence.units = {one(0.1L,0.9L,0.9L,2,2),
                       one(0.1L,0.9L,0.55L,1,1)};
    influence.sums = axis_sum_units(influence.units);
    influence.numeric = axis_numeric(influence.sums);
    axis_influence(influence);
    require(influence.top_unit_id == axis_unit_id(influence.units[0]),
        "top-unit ranking must use absolute Brier margin");
    AxisCellResult folded;
    folded.units = {
        one(0.1L,0.9L,0.2L,1,1), one(0.2L,0.8L,0.4L,1,2),
        one(0.3L,0.7L,0.6L,1,3), one(0.1L,0.6L,0.5L,2,1),
        one(0.4L,0.9L,0.8L,2,2), one(0.2L,0.7L,0.3L,2,3)};
    folded.sums = axis_sum_units(folded.units);
    folded.numeric = axis_numeric(folded.sums);
    axis_folds(folded,"lib19","BC");
    AxisCellResult folded_reordered;
    folded_reordered.units = folded.units;
    reverse(folded_reordered.units.begin(),folded_reordered.units.end());
    folded_reordered.sums = axis_sum_units(folded_reordered.units);
    folded_reordered.numeric = axis_numeric(folded_reordered.sums);
    axis_folds(folded_reordered,"lib19","BC");
    require(folded.fold_count == folded_reordered.fold_count &&
        fabsl(folded.fold_min-folded_reordered.fold_min) <= reorder_tolerance &&
        fabsl(folded.fold_median-folded_reordered.fold_median) <= reorder_tolerance &&
        fabsl(folded.fold_max-folded_reordered.fold_max) <= reorder_tolerance,
        "fold assignment and reconstruction must be deterministic across row order");
    AxisCellResult fold_mismatch = folded;
    fold_mismatch.sums.sum_abs_n += 0.01L;
    axis_folds(fold_mismatch,"lib19","BC_MISMATCH");
    require(fold_mismatch.numeric.status == folded.numeric.status &&
        fold_mismatch.fold_status == "FOLD_RECONSTRUCTION_MISMATCH" &&
        fold_mismatch.folds_preserved == "NA",
        "fold reconstruction mismatch must preserve the primary score and remain nonfatal");
    Identity missing_identity; missing_identity.a = 0;
    AxisSiteDefinition missing_site; missing_site.genotype.push_back(-1);
    unordered_map<int,size_t> missing_slot; missing_slot[0] = 0;
    long double missing_expected = 0.123L;
    require(!axis_expected(missing_identity,missing_site,missing_slot,missing_expected) &&
        missing_expected == 0.123L,
        "missing genotype must never be converted to REF");
    auto must_throw = [&](const function<void()>& operation, const string& message){
        try { operation(); } catch (const exception&) { return; }
        throw runtime_error("candidate-axis self-test failed: " + message);
    };
    must_throw([&](){ axis_parse_observation(
        vector<string>{"1","1","1","-1","2"},"synthetic",1); },
        "negative observations must fail validation");
    must_throw([&](){ strict_ld("nan","synthetic"); },
        "nonfinite values must fail validation");
    require(normalized_contig("chrM") == "m" && normalized_contig("MT") == "mt",
        "mitochondrial contig normalization");
    CandidateHypothesis mito_a, mito_b;
    mito_a.identity.a=0; mito_b.identity.a=1;
    AxisSiteDefinition mito_site; mito_site.key=site_key(1,1); mito_site.tid=1;
    mito_site.pos=1; mito_site.found=true; mito_site.mitochondrial=true;
    mito_site.genotype={0,2};
    AxisObservationRecord mito_observation; mito_observation.tid=1;
    mito_observation.pos=1; mito_observation.ref=5; mito_observation.alt=5;
    unordered_map<int,size_t> mito_slots; mito_slots[0]=0; mito_slots[1]=1;
    const AxisCellResult mito_result = evaluate_candidate_axis_cell(
        mito_a,mito_b,vector<AxisObservationRecord>{mito_observation},
        vector<long>{0},vector<AxisSiteDefinition>{mito_site},mito_slots,
        0.001L,0.001L,10,0.3L);
    require(mito_result.excluded_mito == 1 && mito_result.sums.w == 0.0L,
        "mitochondrial sites must never enter W, N, D, or M");
    require(legacy_pair_manifest_schema_allowed(
        "identity_reconciliation_score_pair_manifest_v1"), "legacy V1 accepted");
    require(legacy_pair_manifest_schema_allowed(
        "identity_reconciliation_score_pair_manifest_v2"), "legacy V2 accepted");
    require(!legacy_pair_manifest_schema_allowed(AXIS_MANIFEST_SCHEMA),
        "candidate-axis manifest rejected by legacy mode");
    require(!legacy_pair_manifest_schema_allowed("unsupported"),
        "unsupported legacy schema rejected");
    auto operational_pair = [](){
        vector<CandidateHypothesis> pair(2);
        for (CandidateHypothesis& candidate : pair){
            candidate.library="lib19"; candidate.barcode="BC";
            candidate.score_pair_id="lib19:BC:CANDIDATE_AXIS_PILOT";
            candidate.schema_version=AXIS_MANIFEST_SCHEMA;
            candidate.score_population_scope="APPLIED_REASSIGNMENT";
            candidate.population_votes_in_authoritative_event="TRUE";
            candidate.supported_event_key="lib19|EVENT|A+B";
            candidate.selected_supported_event_id="EVENT";
            candidate.selected_supported_event_proposal="A+B";
            candidate.reconciliation_event_id="EVENT";
            candidate.reconciliation_event_class="SUPPORTED_EXACT_SWAP";
            candidate.reconciliation_event_confidence="DECISIVE";
            candidate.reconciliation_final_action="REASSIGN_GENOTYPE";
            candidate.reconciliation_decision_confidence="DECISIVE";
            candidate.reconciliation_reassignment_applied="TRUE";
            candidate.original_demux_assignment="A";
            candidate.current_donor_genotype="A";
            candidate.pair_construction_mode="RECONCILIATION_NOMINATED_SWAP";
            candidate.score_scope_contract=AXIS_OPERATIONAL_CONTRACT;
            candidate.scoreable=true;
        }
        pair[0].donor_genotype="A";
        pair[0].score_pair_role="ORIGINAL_ALLOWED_DEMUX";
        pair[0].candidate_origin="ORIGINAL_ALLOWED_DEMUX";
        pair[0].identity.a=0;
        pair[1].donor_genotype="A+B";
        pair[1].score_pair_role="RECONCILIATION_NOMINATED_SWAP";
        pair[1].candidate_origin="RECONCILIATION_NOMINATED_SWAP";
        pair[1].candidate_b_fixed_identity="A+B";
        pair[1].reconciliation_nominated_swap="A+B";
        pair[1].identity.a=0; pair[1].identity.b=1;
        return pair;
    };
    vector<CandidateHypothesis> valid_pair = operational_pair();
    const CandidateAxisPair resolved = candidate_axis_pair(
        valid_pair,vector<size_t>{0,1},1,"lib19");
    require(resolved.original == 0 && resolved.proposed == 1,
        "operational pair roles and fixed A/B orientation");
    vector<CandidateHypothesis> duplicate_role = operational_pair();
    duplicate_role[1].score_pair_role="ORIGINAL_ALLOWED_DEMUX";
    must_throw([&](){ candidate_axis_pair(
        duplicate_role,vector<size_t>{0,1},1,"lib19"); },
        "duplicate candidate-axis role must fail");
    vector<CandidateHypothesis> interchanged = operational_pair();
    interchanged[0].score_scope_contract=AXIS_RETAINED_CONTRACT;
    interchanged[1].score_scope_contract=AXIS_RETAINED_CONTRACT;
    must_throw([&](){ candidate_axis_pair(
        interchanged,vector<size_t>{0,1},1,"lib19"); },
        "operational and retained contracts must not interchange");
    cout << "PASS: tetra_score_calls candidate-axis self-test" << endl;
    return 0;
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
        "  --strict | --best-effort       default best-effort\n"
        "  --candidate_manifest FILE       targeted identity-reconciliation hypotheses\n"
        "  --pileup-sites FILE             optional headerless pileup site table for fold scoring\n"
        "  --pileup-observations FILE      optional headerless pileup observation table\n"
        "  --pileup-molecules FILE         optional molecule-aware pileup sidecar\n"
        "  --site-folds N                  deterministic site folds (default 5)\n"
        "  --site-fold-output FILE         candidate-by-fold score output\n"
        "  --probability-output FILE       targeted original-vs-reconciliation-swap probability\n"
        "  --probability-resamples N       bootstrap/downsample replicates (default 100)\n"
        "  --probability-seed N            deterministic resampling seed (default 1729)\n"
        "  --poor-fit-residual X           flag both candidates as poor fits above X (default 0.30)\n"
        "  --score-prefix STR              prefix candidate score columns (e.g. atac)\n"
        "\nStandalone fixed-pair candidate-axis pilot (not a correctness probability):\n"
        "  --candidate-axis-output FILE    site-balanced raw candidate-axis rows\n"
        "  --candidate-axis-temp-dir DIR   existing absolute job-local temp root\n"
        "  --candidate-axis-self-test      deterministic math/parser/pair checks\n"
        "  Candidate-axis mode also requires explicit --samples, --candidate_manifest,\n"
        "  --pileup-sites, --pileup-observations, --libname, --error_ref,\n"
        "  --error_alt, --min_evidence, and --poor-fit-residual. It is standalone:\n"
        "  counts, assignments, molecule evidence, resamples, and legacy outputs are rejected.\n"
        "  candidate_axis_position_raw is uncalibrated and must never be called confidence,\n"
        "  certainty, accuracy, FDR, posterior probability, or correctness probability.\n");
}

int main(int argc, char** argv){
    string counts, samples_path, assignments, diagnostics, output, runnerups, panel_path, libname="NA";
    string species_counts, species_condf, species_samples_path, condf;
    string candidate_manifest, pileup_sites, pileup_observations;
    string pileup_molecules, site_fold_output, probability_output, score_prefix;
    string candidate_axis_output, candidate_axis_temp_dir;
    int site_folds = 5;
    int probability_resamples = 100;
    uint64_t probability_seed = 1729;
    long min_evidence = 10;
    double e_ref = 0.001, e_alt = 0.001;
    double poor_fit_residual = 0.30;
    double species_support_threshold = 0.70;
    bool strict = false;
    bool candidate_axis_self_test_requested = false;
    bool explicit_error_ref = false, explicit_error_alt = false;
    bool explicit_min_evidence = false, explicit_poor_fit = false;
    bool explicit_probability_resamples = false, explicit_probability_seed = false;
    bool explicit_site_folds = false, explicit_threads = false;
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
        else if (a == "--threads") { string tmp; need(tmp); explicit_threads = true; }
        else if (a == "--min_evidence") { string tmp; need(tmp); min_evidence = atol(tmp.c_str()); explicit_min_evidence = true; }
        else if (a == "--error_ref") { string tmp; need(tmp); e_ref = atof(tmp.c_str()); explicit_error_ref = true; }
        else if (a == "--error_alt") { string tmp; need(tmp); e_alt = atof(tmp.c_str()); explicit_error_alt = true; }
        else if (a == "--species_support_threshold") { string tmp; need(tmp); species_support_threshold = atof(tmp.c_str()); }
        else if (a == "--candidate_manifest") need(candidate_manifest);
        else if (a == "--pileup-sites") need(pileup_sites);
        else if (a == "--pileup-observations") need(pileup_observations);
        else if (a == "--pileup-molecules") need(pileup_molecules);
        else if (a == "--site-fold-output") need(site_fold_output);
        else if (a == "--probability-output") need(probability_output);
        else if (a == "--candidate-axis-output") need(candidate_axis_output);
        else if (a == "--candidate-axis-temp-dir") need(candidate_axis_temp_dir);
        else if (a == "--candidate-axis-self-test") candidate_axis_self_test_requested = true;
        else if (a == "--score-prefix") need(score_prefix);
        else if (a == "--site-folds") { string tmp; need(tmp); site_folds = atoi(tmp.c_str()); explicit_site_folds = true; }
        else if (a == "--probability-resamples") { string tmp; need(tmp); probability_resamples = atoi(tmp.c_str()); explicit_probability_resamples = true; }
        else if (a == "--probability-seed") { string tmp; need(tmp); probability_seed = strtoull(tmp.c_str(), NULL, 10); explicit_probability_seed = true; }
        else if (a == "--poor-fit-residual") { string tmp; need(tmp); poor_fit_residual = atof(tmp.c_str()); explicit_poor_fit = true; }
        else if (a == "--strict") strict = true;
        else if (a == "--best-effort") strict = false;
        else if (a == "--help" || a == "-h") { usage(); return 0; }
        else die("unknown argument: " + a);
    }
    if (candidate_axis_self_test_requested){
        try { return candidate_axis_self_test(); }
        catch (const exception& error){
            fprintf(stderr, "ERROR: %s\n", error.what());
            return 1;
        }
    }
    if (!candidate_axis_output.empty()){
        const bool incompatible = !counts.empty() || !assignments.empty() ||
            !diagnostics.empty() || !output.empty() || !runnerups.empty() ||
            !panel_path.empty() || !species_counts.empty() ||
            !species_condf.empty() || !species_samples_path.empty() ||
            !condf.empty() || !pileup_molecules.empty() ||
            !site_fold_output.empty() || !probability_output.empty() ||
            !score_prefix.empty() || explicit_probability_resamples ||
            explicit_probability_seed || explicit_site_folds || explicit_threads;
        if (incompatible)
            die("--candidate-axis-output is standalone and rejects legacy scoring/count/molecule/fold/resample options");
        if (samples_path.empty() || candidate_manifest.empty() ||
                pileup_sites.empty() || pileup_observations.empty() ||
                candidate_axis_temp_dir.empty() ||
                libname.empty() || libname == "NA" || !explicit_error_ref ||
                !explicit_error_alt || !explicit_min_evidence || !explicit_poor_fit){
            usage();
            return 1;
        }
        try {
            if (setlocale(LC_NUMERIC, "C") == NULL)
                throw runtime_error(
                    "candidate-axis mode could not establish the C numeric locale");
            run_candidate_axis(samples_path, candidate_manifest, pileup_sites,
                pileup_observations, candidate_axis_output,
                candidate_axis_temp_dir, libname,
                strict_ld(axis_fmt(e_ref), "--error_ref"),
                strict_ld(axis_fmt(e_alt), "--error_alt"), min_evidence,
                strict_ld(axis_fmt(poor_fit_residual), "--poor-fit-residual"));
            return 0;
        } catch (const exception& error){
            fprintf(stderr, "ERROR: %s\n", error.what());
            return 1;
        }
    }
    if (counts.empty() || samples_path.empty() || assignments.empty() || output.empty() || (candidate_manifest.empty() && diagnostics.empty())){
        usage();
        return 1;
    }
    if (probability_resamples < 0)
        die("--probability-resamples must be non-negative");
    if (poor_fit_residual < 0.0 || poor_fit_residual > 1.0)
        die("--poor-fit-residual must be within [0,1]");

    vector<string> samples = load_samples(samples_path);
    unordered_map<string,int> sample2idx;
    for (int i=0; i<(int)samples.size(); ++i) sample2idx[samples[i]] = i;
    if (!candidate_manifest.empty()){
        unordered_map<unsigned long, vector<size_t>> candidate_by_cell;
        vector<CandidateHypothesis> candidates = load_candidate_manifest(candidate_manifest, sample2idx, candidate_by_cell);
        if (candidates.empty()) die("candidate manifest contained no hypotheses");
        score_candidate_counts(counts, candidates, candidate_by_cell, e_ref, e_alt);
        write_candidate_scores(output, candidates, candidate_by_cell, min_evidence, score_prefix);
        if (!site_fold_output.empty()){
            score_site_folds(pileup_sites, pileup_observations, site_fold_output, candidates, candidate_by_cell, (int)samples.size(), site_folds, e_ref, e_alt);
        }
        if (!probability_output.empty()){
            write_pairwise_probability_scores(
                pileup_sites, pileup_observations, pileup_molecules,
                probability_output, candidates, candidate_by_cell,
                (int)samples.size(), e_ref, e_alt, min_evidence,
                probability_resamples, probability_seed,
                poor_fit_residual);
        }
        fprintf(stderr, "Wrote targeted identity hypothesis scores for %lu candidate rows to %s\n", (unsigned long)candidates.size(), output.c_str());
        return 0;
    }
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
