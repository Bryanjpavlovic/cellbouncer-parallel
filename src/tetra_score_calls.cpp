#include <algorithm>
#include <atomic>
#include <cctype>
#include <cerrno>
#include <chrono>
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
#include <iomanip>
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

static string uppercase(string s){
    transform(s.begin(),s.end(),s.begin(),[](unsigned char c){
        return static_cast<char>(toupper(c));
    });
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

struct JointStreamingDigest {
    uint64_t value = 1469598103934665603ULL;
    unsigned long long uncompressed_bytes = 0;
    void update(const char* data, size_t size){
        for (size_t i=0;i<size;++i){
            value ^= static_cast<unsigned char>(data[i]);
            value *= 1099511628211ULL;
        }
        uncompressed_bytes += size;
    }
    string text() const {
        ostringstream out;
        out << "fnv1a64:" << hex << setw(16) << setfill('0') << value;
        return out.str();
    }
};

static string joint_file_content_digest(const string& path){
    ifstream input(path.c_str(),ios::binary);
    if (!input) throw runtime_error("could not open cache component for digest: "+path);
    JointStreamingDigest digest;
    vector<char> buffer(1024*1024);
    while (input){
        input.read(buffer.data(),buffer.size());
        const streamsize count=input.gcount();
        if (count>0) digest.update(buffer.data(),static_cast<size_t>(count));
    }
    if (!input.eof())
        throw runtime_error("failed reading cache component for digest: "+path);
    return digest.text();
}

struct AxisObservationRecord {
    unsigned long barcode = 0;
    int32_t tid = -1;
    int32_t pos = -1;
    double ref = 0.0;
    double alt = 0.0;
};

// Fixed-width on-disk record used only in the bounded joint-doublet molecule
// spool.  Strings are represented by a validated basis code so no pointers are
// ever serialized.
struct JointMoleculeRecord {
    unsigned long barcode = 0;
    uint64_t molecule = 0;
    uint8_t basis = 0;
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

static void axis_write_binary_record(
        ofstream& out, const JointMoleculeRecord& record,
        const string& path){
    out.write(reinterpret_cast<const char*>(&record),sizeof(record));
    if (!out)
        throw runtime_error("failed writing joint molecule temporary record: "+path);
}

static bool axis_read_binary_record(
        ifstream& in, JointMoleculeRecord& record,
        const string& path){
    in.read(reinterpret_cast<char*>(&record),sizeof(record));
    if (in.gcount()==0 && in.eof()) return false;
    if (in.gcount()!=(streamsize)sizeof(record))
        throw runtime_error("truncated joint molecule temporary record: "+path);
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
        AxisResourceAudit& audit,
        JointStreamingDigest* streaming_digest = NULL,
        uint64_t memory_limit = UINT64_MAX){
    const uint64_t per_site=sizeof(AxisSiteDefinition)+donors.size()+64;
    if (selected_keys.size()>memory_limit/per_site)
        throw runtime_error("selected genotype dictionary exceeds extraction memory budget");
    uint64_t allocated=selected_keys.size()*per_site;
    vector<AxisSiteDefinition> sites(selected_keys.size());
    for (size_t i = 0; i < selected_keys.size(); ++i) sites[i].key = selected_keys[i];
    gzFile input = gzopen(path.c_str(), "rb");
    if (!input) throw runtime_error("could not open candidate-axis site definitions: " + path);
    char buffer[1<<20];
    unsigned long long line_no = 0;
    try {
        while (gzgets(input, buffer, sizeof(buffer))){
            ++line_no;
            if (streaming_digest)
                streaming_digest->update(buffer,strlen(buffer));
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
            } else {
                const uint64_t strings=definition.contig.capacity()+
                    definition.ref_allele.capacity()+definition.alt_allele.capacity();
                if (strings>memory_limit-allocated)
                    throw runtime_error("selected allele metadata exceeds extraction memory budget");
                allocated+=strings;
                sites[index] = move(definition);
            }
        }
        const int close_status=gzclose(input);
        input = NULL;
        if (close_status != Z_OK)
            throw runtime_error("failed closing candidate-axis site definitions: " + path);
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
        AxisResourceAudit& audit,
        size_t minimum_bucket_count = 1){
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
    bucket_count = max(bucket_count, minimum_bucket_count);
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

static string joint_stage_observations_one_pass(
        const string& path,
        const unordered_map<unsigned long,CandidateAxisPair>& targets,
        const string& temp_path,
        AxisResourceAudit& audit,
        unordered_map<unsigned long,unsigned long long>& rows_by_barcode,
        vector<uint64_t>& selected_site_keys,
        JointStreamingDigest* streaming_digest = NULL){
    const string spool_path = temp_path + "/joint_selected_observations.bin";
    ofstream spool(spool_path.c_str(), ios::binary);
    if (!spool)
        throw runtime_error("could not create joint-doublet observation spool: " +
            spool_path);
    gzFile input = gzopen(path.c_str(), "rb");
    if (!input)
        throw runtime_error("could not open joint-doublet observations: " + path);

    const size_t key_capacity = (128ULL * 1024ULL * 1024ULL) /
        sizeof(uint64_t);
    vector<uint64_t> key_chunk;
    key_chunk.reserve(min<size_t>(key_capacity, 1000000));
    vector<string> key_runs;
    char buffer[1<<20];
    unsigned long long line_no = 0;
    try {
        while (gzgets(input, buffer, sizeof(buffer))){
            ++line_no;
            if (streaming_digest)
                streaming_digest->update(buffer,strlen(buffer));
            string line(buffer);
            line.erase(remove(line.begin(), line.end(), '\n'), line.end());
            line.erase(remove(line.begin(), line.end(), '\r'), line.end());
            if (line.empty()) continue;
            vector<string> fields = split_tsv_strict(line);
            if (fields.empty()) continue;
            unsigned long barcode = 0;
            try {
                barcode = strict_barcode_number(
                    fields[0], path + ": line " + to_string(line_no));
            } catch (const exception&) {
                continue;
            }
            if (targets.find(barcode) == targets.end()) continue;
            const AxisObservationRecord record = axis_parse_observation(
                fields, path, line_no);
            axis_write_binary_record(spool, record, spool_path);
            ++audit.target_rows;
            ++rows_by_barcode[barcode];
            key_chunk.push_back(site_key(record.tid, record.pos));
            if (key_chunk.size() >= key_capacity){
                key_runs.push_back(axis_spill_key_run(
                    key_chunk, temp_path, key_runs.size()));
                key_chunk.reserve(min<size_t>(key_capacity, 1000000));
            }
        }
        if (gzclose(input) != Z_OK)
            throw runtime_error(
                "failed closing joint-doublet observations: " + path);
        input = NULL;
        spool.close();
        if (!spool)
            throw runtime_error(
                "failed closing joint-doublet observation spool: " +
                spool_path);
    } catch (...) {
        if (input) gzclose(input);
        spool.close();
        throw;
    }

    audit.target_barcodes = rows_by_barcode.size();
    for (const auto& item : rows_by_barcode)
        audit.largest_barcode_rows = max(
            audit.largest_barcode_rows, item.second);

    if (key_runs.empty()){
        sort(key_chunk.begin(), key_chunk.end());
        key_chunk.erase(unique(key_chunk.begin(), key_chunk.end()),
            key_chunk.end());
        selected_site_keys.swap(key_chunk);
    } else {
        if (!key_chunk.empty())
            key_runs.push_back(axis_spill_key_run(
                key_chunk, temp_path, key_runs.size()));
        vector<ifstream> inputs(key_runs.size());
        priority_queue<AxisKeyCursor,vector<AxisKeyCursor>,AxisKeyCursorGreater>
            heap;
        for (size_t i = 0; i < key_runs.size(); ++i){
            inputs[i].open(key_runs[i].c_str(), ios::binary);
            if (!inputs[i])
                throw runtime_error(
                    "could not open joint-doublet selected-site run: " +
                    key_runs[i]);
            uint64_t key = 0;
            if (axis_read_key(inputs[i], key, key_runs[i])){
                AxisKeyCursor cursor;
                cursor.key = key;
                cursor.run = i;
                heap.push(cursor);
            }
        }
        bool have_prior = false;
        uint64_t prior = 0;
        while (!heap.empty()){
            const AxisKeyCursor cursor = heap.top();
            heap.pop();
            if (!have_prior || cursor.key != prior){
                selected_site_keys.push_back(cursor.key);
                prior = cursor.key;
                have_prior = true;
            }
            uint64_t next = 0;
            if (axis_read_key(
                    inputs[cursor.run], next, key_runs[cursor.run])){
                AxisKeyCursor following;
                following.key = next;
                following.run = cursor.run;
                heap.push(following);
            }
        }
        for (size_t i = 0; i < inputs.size(); ++i){
            inputs[i].close();
            unlink(key_runs[i].c_str());
        }
    }
    audit.spill_runs = key_runs.size();
    audit.unique_site_keys = selected_site_keys.size();
    audit.selected_key_bytes =
        selected_site_keys.capacity() * sizeof(uint64_t);
    return spool_path;
}

static vector<string> joint_partition_staged_observations(
        const string& spool_path,
        const unordered_map<unsigned long,size_t>& assignment,
        size_t bucket_count,
        const string& temp_path){
    vector<string> paths(bucket_count);
    vector<ofstream> outputs(bucket_count);
    for (size_t i = 0; i < bucket_count; ++i){
        paths[i] = temp_path + "/joint_bucket_" + to_string(i) + ".bin";
        outputs[i].open(paths[i].c_str(), ios::binary);
        if (!outputs[i])
            throw runtime_error(
                "could not create joint-doublet bucket: " + paths[i]);
    }
    ifstream input(spool_path.c_str(), ios::binary);
    if (!input)
        throw runtime_error(
            "could not reopen joint-doublet observation spool: " +
            spool_path);
    AxisObservationRecord record;
    while (axis_read_binary_record(input, record, spool_path)){
        const auto found = assignment.find(record.barcode);
        if (found == assignment.end())
            throw runtime_error(
                "joint-doublet observation spool contains nontarget barcode");
        axis_write_binary_record(outputs[found->second], record,
            paths[found->second]);
    }
    input.close();
    for (size_t i = 0; i < outputs.size(); ++i){
        outputs[i].close();
        if (!outputs[i])
            throw runtime_error(
                "failed closing joint-doublet bucket: " + paths[i]);
    }
    unlink(spool_path.c_str());
    return paths;
}

static uint8_t joint_molecule_basis_code(const string& raw){
    const string basis=trim(raw);
    if (basis=="UB_GX") return 1;
    if (basis=="UB_GN") return 2;
    if (basis=="QNAME_FALLBACK") return 3;
    throw runtime_error("unsupported pileup molecule basis: " + basis);
}

static string joint_molecule_basis_name(uint8_t basis){
    if (basis==1) return "UB_GX";
    if (basis==2) return "UB_GN";
    if (basis==3) return "QNAME_FALLBACK";
    return "UNKNOWN";
}

static string joint_molecule_fold_basis_name(uint8_t basis){
    if (basis==1 || basis==2) return "UMI_GENE";
    if (basis==3) return "READ_NAME";
    return "UNKNOWN";
}

static bool joint_molecule_less(
        const JointMoleculeRecord& left, const JointMoleculeRecord& right){
    if (left.barcode!=right.barcode) return left.barcode<right.barcode;
    if (left.basis!=right.basis) return left.basis<right.basis;
    if (left.molecule!=right.molecule) return left.molecule<right.molecule;
    if (left.tid!=right.tid) return left.tid<right.tid;
    return left.pos<right.pos;
}

static string joint_stage_molecules_one_pass(
        const string& path,
        const unordered_map<unsigned long,CandidateAxisPair>& targets,
        const vector<uint64_t>& selected_site_keys,
        const string& temp_path,
        unordered_map<unsigned long,unsigned long long>& rows_by_barcode,
        unordered_map<unsigned long,long>& malformed_by_barcode,
        JointStreamingDigest* streaming_digest = NULL,
        set<uint64_t>* molecule_site_union = NULL){
    if (path.empty() || !file_exists(path)) return "";
    const string spool_path=temp_path+"/joint_selected_molecules.bin";
    ofstream spool(spool_path.c_str(),ios::binary);
    if (!spool)
        throw runtime_error("could not create joint-doublet molecule spool: "+spool_path);
    gzFile input=gzopen(path.c_str(),"rb");
    if (!input)
        throw runtime_error("could not open joint-doublet molecule sidecar: "+path);
    char buffer[1<<20];
    unsigned long long line_no=0, parsed_rows=0, malformed_rows=0;
    try {
        while (gzgets(input,buffer,sizeof(buffer))){
            ++line_no;
            if (streaming_digest)
                streaming_digest->update(buffer,strlen(buffer));
            string line(buffer);
            line.erase(remove(line.begin(),line.end(),'\n'),line.end());
            line.erase(remove(line.begin(),line.end(),'\r'),line.end());
            if (line.empty()) continue;
            const vector<string> fields=split_tsv_strict(line);
            unsigned long barcode=0;
            try {
                if (fields.empty()) throw runtime_error("missing barcode");
                barcode=strict_barcode_number(fields[0],path+": line "+to_string(line_no));
            } catch (const exception&){
                ++malformed_rows;
                continue;
            }
            if (targets.find(barcode)==targets.end()) continue;
            try {
                if (fields.size()!=7)
                    throw runtime_error("expected exactly seven columns");
                char* molecule_end=NULL;
                errno=0;
                const uint64_t molecule=strtoull(fields[1].c_str(),&molecule_end,10);
                if (errno!=0 || molecule_end==fields[1].c_str() ||
                        *molecule_end!='\0')
                    throw runtime_error("invalid molecule hash");
                char* tid_end=NULL; char* pos_end=NULL;
                errno=0;
                const long tid=strtol(fields[3].c_str(),&tid_end,10);
                const long pos=strtol(fields[4].c_str(),&pos_end,10);
                if (errno!=0 || tid_end==fields[3].c_str() || *tid_end!='\0' ||
                        pos_end==fields[4].c_str() || *pos_end!='\0' ||
                        tid<INT32_MIN || tid>INT32_MAX ||
                        pos<INT32_MIN || pos>INT32_MAX)
                    throw runtime_error("invalid molecule site coordinate");
                const long double ref=strict_ld(fields[5],path+": molecule ref");
                const long double alt=strict_ld(fields[6],path+": molecule alt");
                if (ref<0.0L || alt<0.0L || ref+alt<=0.0L)
                    throw runtime_error("nonpositive molecule allele depth");
                const uint64_t key=site_key((int32_t)tid,(int32_t)pos);
                if (molecule_site_union) molecule_site_union->insert(key);
                if (!molecule_site_union && !binary_search(selected_site_keys.begin(),selected_site_keys.end(),key))
                    continue;
                JointMoleculeRecord record;
                record.barcode=barcode;
                record.molecule=molecule;
                record.basis=joint_molecule_basis_code(fields[2]);
                record.tid=(int32_t)tid; record.pos=(int32_t)pos;
                record.ref=(double)ref; record.alt=(double)alt;
                axis_write_binary_record(spool,record,spool_path);
                if(!spool)throw runtime_error("MOLECULE_SPOOL_WRITE_FAILED");
                ++rows_by_barcode[barcode];
                ++parsed_rows;
            } catch (const exception& error){
                if(!spool)throw;
                ++malformed_rows;
                ++malformed_by_barcode[barcode];
            }
        }
        const int close_status=gzclose(input);
        input=NULL;
        if (close_status!=Z_OK)
            throw runtime_error("failed closing joint-doublet molecule sidecar: "+path);
        spool.close();
        if (!spool)
            throw runtime_error("failed closing joint-doublet molecule spool: "+spool_path);
    } catch (...){
        if (input) gzclose(input);
        spool.close();
        throw;
    }
    if (parsed_rows==0 && malformed_rows>0)
        throw runtime_error(path+": molecule sidecar has no valid seven-column target records");
    return spool_path;
}

static vector<string> joint_partition_staged_molecules(
        const string& spool_path,
        const unordered_map<unsigned long,size_t>& assignment,
        size_t bucket_count, const string& temp_path){
    if (spool_path.empty()) return vector<string>();
    vector<string> paths(bucket_count);
    vector<ofstream> outputs(bucket_count);
    for (size_t i=0;i<bucket_count;++i){
        paths[i]=temp_path+"/joint_molecule_bucket_"+to_string(i)+".bin";
        outputs[i].open(paths[i].c_str(),ios::binary);
        if (!outputs[i])
            throw runtime_error("could not create joint molecule bucket: "+paths[i]);
    }
    ifstream input(spool_path.c_str(),ios::binary);
    if (!input)
        throw runtime_error("could not reopen joint molecule spool: "+spool_path);
    JointMoleculeRecord record;
    while (axis_read_binary_record(input,record,spool_path)){
        const auto found=assignment.find(record.barcode);
        if (found==assignment.end())
            throw runtime_error("joint molecule spool contains nontarget barcode");
        axis_write_binary_record(outputs[found->second],record,paths[found->second]);
    }
    input.close();
    for (size_t i=0;i<outputs.size();++i){
        outputs[i].close();
        if (!outputs[i])
            throw runtime_error("failed closing joint molecule bucket: "+paths[i]);
    }
    unlink(spool_path.c_str());
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

// -------------------------------------------------------------------------
// Generalized K=1 versus K=2 + ambient scorer (standalone bounded mode)
// -------------------------------------------------------------------------

static const char* JOINT_DOUBLET_SCHEMA =
    "joint_doublet_site_and_molecule_candidate_score_v3";
static const char* JOINT_DOUBLET_MANIFEST_SCHEMA =
    "joint_doublet_candidate_manifest_v4";
static const char* JOINT_DOUBLET_FORMULA =
    "K1_LOCKED_STATE_VS_K2_ADDED_STATE_PLUS_FIXED_AMBIENT_V1";
static const char* JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION =
    "JOINT_DOUBLET_TARGETED_V4_20260920";
static const char* JOINT_DOUBLET_FOLD_VERSION =
    "LEGACY_PROJECT_FNV1A64_FOLD_V1";
static const char* JOINT_DOUBLET_MOLECULE_FORMULA =
    "LINKED_UNIT_EQUAL_WEIGHT_NORMALIZED_SITE_LL_V1";
static const char* JOINT_DOUBLET_MOLECULE_FOLD_VERSION =
    "JOINT_DOUBLET_LINKED_UNIT_SORTED_ROUND_ROBIN_FNV1A64_V1";
static const char* JOINT_DOUBLET_MOLECULE_SIDECAR_SCHEMA =
    "pileup_molecules_headerless_7col_v1";
static const char* JOINT_DOUBLET_TARGETED_CACHE_SCHEMA =
    "joint_doublet_targeted_evidence_cache_v4";

struct JointComposition {
    vector<pair<int,long double>> members;
    string text;
    string missing_donors = "NONE";
    bool valid = false;
};

struct JointHypothesis {
    string library;
    string barcode;
    unsigned long encoded_barcode = 0;
    string candidate_id;
    string locked_state;
    JointComposition locked;
    string second_state;
    JointComposition second;
    JointComposition ambient;
    string candidate_origin = "UNSPECIFIED";
    string exhaustive_fallback = "FALSE";
    string nomination_modalities = "NONE";
    string ambient_status = "UNAVAILABLE";
    long double rho_requested = 0.0L;
    long double rho_effective = 0.0L;
};

struct JointSiteUnit {
    int32_t tid = -1;
    int32_t pos = -1;
    // The targeted workflow must never treat a numeric tid as a portable
    // genomic identity.  These fields are copied from the authoritative site
    // definition and are used for folds and synthetic-control pooling.
    string contig;
    string ref_allele;
    string alt_allele;
    // Numeric identifiers and the linear probability representation are
    // populated once when a cell is compiled.  Replicate evaluators use
    // these fields directly and never rebuild string-keyed site maps.
    uint32_t compiled_site_id = numeric_limits<uint32_t>::max();
    long double ref = 0.0L;
    long double alt = 0.0L;
    long double q_locked = 0.0L;
    long double q_second = 0.0L;
    long double q_ambient = 0.0L;
    long double probability_intercept = NAN;
    long double probability_slope = NAN;
};

struct JointLinkedUnit {
    uint64_t molecule = 0;
    string molecule_id;
    uint8_t basis = 0;
    string parent_origin = "RECIPIENT";
    uint32_t compiled_unit_id = numeric_limits<uint32_t>::max();
    // Candidate applicability is compiled once.  Replicate fits consume this
    // byte directly instead of rediscovering informative units or allocating
    // a vector of string/unit keys for every candidate and replicate.
    uint8_t genotype_distinguishable = 0;
    vector<JointSiteUnit> sites;
    int fold = -1;
};

struct JointFit {
    long double alpha = NAN;
    long double balanced = NAN;
    long double raw = NAN;
};

// Enabled only by the normalized-cache analysis entry point.  These counters
// are intentionally process-local: a shard owns one process and writes one
// deterministic metrics row after all scientific rows are complete.
static bool joint_count_targeted_work = false;
static atomic<uint64_t> joint_targeted_optimizer_calls(0);
static atomic<uint64_t> joint_targeted_likelihood_evaluations(0);
static atomic<uint64_t> joint_targeted_derivative_evaluations(0);
static atomic<uint64_t> joint_targeted_optimized_fit_calls(0);
static atomic<uint64_t> joint_targeted_reference_fit_calls(0);
static atomic<uint64_t> joint_targeted_row_scans(0);

struct JointResult {
    string status = "NO_OBSERVATIONS";
    string genotype_visibility = "UNAVAILABLE";
    string warnings = "NONE";
    long double k1_balanced = NAN;
    long double k2_balanced = NAN;
    long double delta_balanced = NAN;
    long double k1_raw = NAN;
    long double k2_raw = NAN;
    long double delta_raw = NAN;
    long double contributor_only_balanced = NAN;
    long double contributor_only_raw = NAN;
    string preferred_site_model = "UNAVAILABLE";
    long double alpha = NAN;
    long double alpha_low = NAN;
    long double alpha_high = NAN;
    long double discriminating_depth = 0.0L;
    long double top_site_fraction = NAN;
    long double fold_support_fraction = NAN;
    long double fold_min_delta = NAN;
    long double fold_median_delta = NAN;
    int common_sites = 0;
    int discriminating_sites = 0;
    int folds_evaluable = 0;
    long excluded_missing_definition = 0;
    long excluded_mitochondrial = 0;
    long excluded_nonpositive = 0;
    long excluded_missing_genotype = 0;
    bool genotype_equivalent = false;
    bool replacement_like = false;
    string molecule_status = "MOLECULE_SIDECAR_UNAVAILABLE";
    string molecule_unusable_reason = "PILEUP_MOLECULES_NOT_PROVIDED";
    string molecule_warnings = "NONE";
    long double molecule_k1 = NAN;
    long double molecule_k2 = NAN;
    long double molecule_delta = NAN;
    long double molecule_contributor_only = NAN;
    string preferred_molecule_model = "UNAVAILABLE";
    long double molecule_alpha = NAN;
    long double molecule_alpha_low = NAN;
    long double molecule_alpha_high = NAN;
    long double molecule_effective_units = NAN;
    long double molecule_maximum_influence_fraction = NAN;
    long double molecule_heldout_delta = NAN;
    long double molecule_heldout_support_fraction = NAN;
    long double molecule_fold_min_delta = NAN;
    long double molecule_fold_median_delta = NAN;
    long double molecule_without_top_delta = NAN;
    long double molecule_without_top_alpha = NAN;
    long double molecule_umi_gene_fraction = NAN;
    long double molecule_query_name_fraction = NAN;
    int molecule_units = 0;
    int molecule_discriminating_units = 0;
    int molecule_multi_snp_units = 0;
    int molecule_total_snps = 0;
    int molecule_fold_count = 0;
    int molecule_folds_evaluable = 0;
    string molecule_snps_per_unit = "NA";
    string molecule_fold_units = "NA";
    string molecule_fold_fitted_fractions = "NA";
    long molecule_malformed_rows = 0;
    vector<JointSiteUnit> cache_site_units;
    vector<JointLinkedUnit> cache_molecule_units;
};

static map<string,int> joint_header_index(
        const vector<string>& header, const string& path){
    map<string,int> result;
    for (size_t i = 0; i < header.size(); ++i){
        const string key = lowercase(trim(header[i]));
        if (key.empty())
            throw runtime_error(path + ": blank manifest header field");
        if (result.count(key))
            throw runtime_error(path + ": duplicate manifest header field: " + key);
        result[key] = (int)i;
    }
    return result;
}

static string joint_field(
        const vector<string>& fields, const map<string,int>& index,
        const string& name, const string& path, long long line_no,
        bool required = true){
    auto found = index.find(lowercase(name));
    if (found == index.end()){
        if (!required) return "";
        throw runtime_error(path + ": missing required manifest column: " + name);
    }
    if (found->second < 0 || found->second >= (int)fields.size()){
        if (!required) return "";
        throw runtime_error(path + ": short manifest row at line " +
            to_string(line_no) + " for column " + name);
    }
    const string value = trim(fields[found->second]);
    if (required && (value.empty() || lowercase(value) == "na"))
        throw runtime_error(path + ": blank required value at line " +
            to_string(line_no) + " column " + name);
    return value;
}

static JointComposition joint_parse_composition(
        const string& raw, const unordered_map<string,int>& sample2idx,
        const string& context, bool required){
    JointComposition result;
    result.text = trim(raw);
    const string lower = lowercase(result.text);
    if (result.text.empty() || lower == "na" || lower == "none" ||
            result.text == "."){
        if (required) result.missing_donors = "MISSING_COMPOSITION";
        return result;
    }
    map<int,long double> aggregated;
    vector<string> missing;
    for (const string& raw_token : split(result.text, ',')){
        const string token = trim(raw_token);
        const size_t colon = token.find_last_of(':');
        if (colon == string::npos)
            throw runtime_error(context + ": composition token must be donor:copies: " + token);
        const string donor = trim(token.substr(0, colon));
        const long double copies = strict_ld(
            token.substr(colon + 1), context + " donor=" + donor);
        if (donor.empty() || copies <= 0.0L)
            throw runtime_error(context + ": donor names and copies must be positive");
        auto found = sample2idx.find(donor);
        if (found == sample2idx.end()) missing.push_back(donor);
        else aggregated[found->second] += copies;
    }
    if (!missing.empty()){
        sort(missing.begin(), missing.end());
        missing.erase(unique(missing.begin(), missing.end()), missing.end());
        result.missing_donors.clear();
        for (size_t i = 0; i < missing.size(); ++i){
            if (i) result.missing_donors += ",";
            result.missing_donors += missing[i];
        }
        return result;
    }
    for (const auto& item : aggregated)
        result.members.push_back(item);
    result.valid = !result.members.empty();
    return result;
}

static vector<JointHypothesis> joint_load_manifest(
        const string& path, const unordered_map<string,int>& sample2idx,
        const string& library,
        unordered_map<unsigned long,vector<size_t>>& by_cell,
        JointStreamingDigest* streaming_digest = NULL){
    gzFile input = gzopen(path.c_str(), "rb");
    if (!input) throw runtime_error("could not open joint-doublet manifest: " + path);
    char buffer[1<<20];
    if (!gzgets(input, buffer, sizeof(buffer))){
        gzclose(input);
        throw runtime_error("empty joint-doublet manifest: " + path);
    }
    if (streaming_digest) streaming_digest->update(buffer,strlen(buffer));
    string header_line(buffer);
    header_line.erase(remove(header_line.begin(), header_line.end(), '\n'), header_line.end());
    header_line.erase(remove(header_line.begin(), header_line.end(), '\r'), header_line.end());
    const map<string,int> index = joint_header_index(
        split_tsv_strict(header_line), path);
    vector<JointHypothesis> hypotheses;
    unordered_set<string> candidate_ids;
    long long line_no = 1;
    try {
        while (gzgets(input, buffer, sizeof(buffer))){
            if (streaming_digest)
                streaming_digest->update(buffer,strlen(buffer));
            ++line_no;
            string line(buffer);
            line.erase(remove(line.begin(), line.end(), '\n'), line.end());
            line.erase(remove(line.begin(), line.end(), '\r'), line.end());
            if (line.empty()) continue;
            const vector<string> fields = split_tsv_strict(line);
            const string schema = joint_field(
                fields,index,"schema_version",path,line_no);
            if (schema != JOINT_DOUBLET_MANIFEST_SCHEMA)
                throw runtime_error(path + ": unsupported schema at line " +
                    to_string(line_no) + ": " + schema);
            if (joint_field(fields,index,"scientific_method_version",path,line_no)!=
                    JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION)
                throw runtime_error(path+
                    ": scientific-method mismatch at line "+to_string(line_no));
            if (joint_field(fields,index,"calibration_library",path,line_no)!="25")
                throw runtime_error(path+
                    ": calibration-library mismatch at line "+to_string(line_no));
            JointHypothesis hypothesis;
            hypothesis.library = joint_field(
                fields,index,"library",path,line_no);
            if (hypothesis.library != library)
                throw runtime_error(path + ": manifest library does not match --libname at line " +
                    to_string(line_no));
            hypothesis.barcode = joint_field(
                fields,index,"barcode",path,line_no);
            if (hypothesis.barcode.size() != 16 ||
                    any_of(hypothesis.barcode.begin(),hypothesis.barcode.end(),
                        [](char base){
                            return base!='A' && base!='C' &&
                                   base!='G' && base!='T';
                        }))
                throw runtime_error(path +
                    ": joint-doublet barcodes must be canonical 16-bp A/C/G/T strings at line " +
                    to_string(line_no));
            hypothesis.encoded_barcode = bc_ul(hypothesis.barcode);
            hypothesis.candidate_id = joint_field(
                fields,index,"candidate_id",path,line_no);
            const string unique_id = hypothesis.barcode + "\t" + hypothesis.candidate_id;
            if (!candidate_ids.insert(unique_id).second)
                throw runtime_error(path + ": duplicate barcode/candidate_id at line " +
                    to_string(line_no));
            hypothesis.locked_state = joint_field(
                fields,index,"locked_state",path,line_no);
            const string locked_text = joint_field(
                fields,index,"locked_copy_vector",path,line_no);
            hypothesis.locked = joint_parse_composition(
                locked_text,sample2idx,path + ": line " + to_string(line_no) +
                " locked_copy_vector",true);
            hypothesis.second_state = joint_field(
                fields,index,"second_state",path,line_no);
            const string second_text = joint_field(
                fields,index,"second_copy_vector",path,line_no);
            hypothesis.second = joint_parse_composition(
                second_text,sample2idx,path + ": line " + to_string(line_no) +
                " second_copy_vector",true);
            const string ambient_text = joint_field(
                fields,index,"ambient_copy_vector",path,line_no,false);
            hypothesis.ambient = joint_parse_composition(
                ambient_text,sample2idx,path + ": line " + to_string(line_no) +
                " ambient_copy_vector",false);
            hypothesis.candidate_origin = joint_field(
                fields,index,"candidate_origin",path,line_no,false);
            if (hypothesis.candidate_origin.empty())
                hypothesis.candidate_origin = "UNSPECIFIED";
            hypothesis.exhaustive_fallback = joint_field(
                fields,index,"exhaustive_fallback",path,line_no,false);
            if (hypothesis.exhaustive_fallback.empty())
                hypothesis.exhaustive_fallback = "FALSE";
            hypothesis.nomination_modalities = joint_field(
                fields,index,"nomination_modalities",path,line_no,false);
            if (hypothesis.nomination_modalities.empty())
                hypothesis.nomination_modalities = "NONE";
            hypothesis.ambient_status = joint_field(
                fields,index,"ambient_status",path,line_no,false);
            if (hypothesis.ambient_status.empty())
                hypothesis.ambient_status = "UNAVAILABLE";
            const string rho_text = joint_field(
                fields,index,"rho",path,line_no,false);
            hypothesis.rho_requested = rho_text.empty() || lowercase(rho_text) == "na" ?
                0.0L : strict_ld(rho_text,path + ": line " + to_string(line_no) + " rho");
            if (hypothesis.rho_requested < 0.0L ||
                    hypothesis.rho_requested > 0.99L)
                throw runtime_error(path + ": rho must be in [0,0.99] at line " +
                    to_string(line_no));
            hypothesis.rho_effective = hypothesis.ambient.valid ?
                hypothesis.rho_requested : 0.0L;
            by_cell[hypothesis.encoded_barcode].push_back(hypotheses.size());
            hypotheses.push_back(hypothesis);
        }
        if (gzclose(input) != Z_OK)
            throw runtime_error("failed closing joint-doublet manifest: " + path);
        input = NULL;
    } catch (...) {
        if (input) gzclose(input);
        throw;
    }
    return hypotheses;
}

static bool joint_expected(
        const JointComposition& composition,
        const AxisSiteDefinition& site,
        const unordered_map<int,size_t>& donor_slot,
        long double& expected, bool allow_missing_members = false){
    if (!composition.valid) return false;
    AxisKahan numerator, denominator;
    for (const auto& member : composition.members){
        auto slot = donor_slot.find(member.first);
        if (slot == donor_slot.end() || slot->second >= site.genotype.size()){
            if (allow_missing_members) continue;
            return false;
        }
        const int8_t gt = site.genotype[slot->second];
        if (gt < 0){
            if (allow_missing_members) continue;
            return false;
        }
        numerator.add(member.second * ((long double)gt / 2.0L));
        denominator.add(member.second);
    }
    if (denominator.value <= 0.0L) return false;
    expected = numerator.value / denominator.value;
    return true;
}

static long double joint_adjust_probability(
        long double probability, long double e_ref, long double e_alt){
    const long double adjusted = probability * (1.0L - e_alt) +
        (1.0L - probability) * e_ref;
    return min(max(adjusted, 1e-15L), 1.0L - 1e-15L);
}

static long double joint_site_log_likelihood(
        const JointSiteUnit& unit, long double alpha, long double rho,
        long double e_ref, long double e_alt, bool balanced){
    if (joint_count_targeted_work){
        ++joint_targeted_likelihood_evaluations;
        ++joint_targeted_row_scans;
    }
    const long double biological =
        (1.0L-alpha)*unit.q_locked + alpha*unit.q_second;
    const long double mixture = (1.0L-rho)*biological +
        rho*unit.q_ambient;
    const long double probability = joint_adjust_probability(
        mixture,e_ref,e_alt);
    const long double raw = unit.alt*logl(probability) +
        unit.ref*logl(1.0L-probability);
    const long double depth = unit.ref + unit.alt;
    return balanced && depth > 0.0L ? raw/depth : raw;
}

static pair<long double,long double> joint_site_log_likelihood_derivatives(
        const JointSiteUnit& unit, long double alpha, long double rho,
        long double e_ref, long double e_alt, bool balanced){
    const long double error_scale=1.0L-e_ref-e_alt;
    const long double biological=
        (1.0L-alpha)*unit.q_locked+alpha*unit.q_second;
    const long double mixture=(1.0L-rho)*biological+rho*unit.q_ambient;
    const long double unclamped=e_ref+error_scale*mixture;
    const long double epsilon=1e-15L;
    const long double probability=min(max(unclamped,epsilon),1.0L-epsilon);
    long double slope=error_scale*(1.0L-rho)*
        (unit.q_second-unit.q_locked);
    if (unclamped<=epsilon || unclamped>=1.0L-epsilon) slope=0.0L;
    const long double depth=unit.ref+unit.alt;
    const long double weight=balanced && depth>0.0L ? 1.0L/depth : 1.0L;
    const long double first=weight*slope*
        (unit.alt/probability-unit.ref/(1.0L-probability));
    const long double second=-weight*slope*slope*
        (unit.alt/(probability*probability)+
         unit.ref/((1.0L-probability)*(1.0L-probability)));
    return make_pair(first,second);
}

template <typename DerivativeFunction,typename ScoreFunction>
static pair<long double,long double> joint_concave_maximum(
        long double max_alpha, const DerivativeFunction& derivatives,
        const ScoreFunction& score){
    if (joint_count_targeted_work) ++joint_targeted_optimizer_calls;
    const pair<long double,long double> at_zero=derivatives(0.0L);
    const pair<long double,long double> at_max=derivatives(max_alpha);
    long double alpha=0.0L;
    if (at_zero.first<=0.0L){
        alpha=0.0L;
    } else if (at_max.first>=0.0L){
        alpha=max_alpha;
    } else {
        long double left=0.0L,right=max_alpha;
        long double current=(left+right)/2.0L;
        for (int iteration=0;iteration<32;++iteration){
            const pair<long double,long double> value=derivatives(current);
            if (value.first>0.0L) left=current;
            else right=current;
            if (right-left<=1e-10L*max(1.0L,max_alpha)) break;
            long double proposed=NAN;
            if (isfinite(value.second) && value.second<0.0L)
                proposed=current-value.first/value.second;
            if (!isfinite(proposed) || proposed<=left || proposed>=right)
                proposed=(left+right)/2.0L;
            if (proposed==current){ left=current; right=current; break; }
            current=proposed;
        }
        alpha=(left+right)/2.0L;
    }
    long double best=score(alpha);
    const long double zero_score=score(0.0L);
    if (zero_score>best){ alpha=0.0L; best=zero_score; }
    const long double max_score=score(max_alpha);
    if (max_score>best){ alpha=max_alpha; best=max_score; }
    return make_pair(alpha,best);
}

template <typename ScoreFunction>
static pair<long double,long double> joint_profile_interval(
        long double max_alpha, long double maximum_alpha,
        long double cutoff, const ScoreFunction& score){
    (void)maximum_alpha;
    // The scientific contract is the inclusive passing interval on exactly
    // this 201-point grid.  Evaluating the literal grid avoids changing the
    // result through an optimization shortcut or interpolation.
    const int steps=200;
    int low=-1,high=-1;
    for (int index=0;index<=steps;++index){
        const long double alpha=max_alpha*(long double)index/(long double)steps;
        const long double value=score(alpha);
        if (isfinite(value) && value>=cutoff){
            if (low<0) low=index;
            high=index;
        }
    }
    if (low<0) return make_pair(NAN,NAN);
    return make_pair(max_alpha*(long double)low/(long double)steps,
                     max_alpha*(long double)high/(long double)steps);
}

static JointFit joint_fit(
        const vector<JointSiteUnit>& units, const vector<size_t>& indices,
        long double rho, long double e_ref, long double e_alt,
        long double max_alpha){
    JointFit result;
    if (indices.empty()) return result;
    auto score = [&](long double alpha, bool balanced){
        AxisKahan sum;
        for (size_t index : indices)
            sum.add(joint_site_log_likelihood(
                units[index],alpha,rho,e_ref,e_alt,balanced));
        return sum.value;
    };
    const auto derivatives=[&](long double alpha){
        AxisKahan first,second;
        for (size_t index : indices){
            const pair<long double,long double> value=
                joint_site_log_likelihood_derivatives(
                    units[index],alpha,rho,e_ref,e_alt,true);
            first.add(value.first); second.add(value.second);
        }
        return make_pair(first.value,second.value);
    };
    const auto balanced_score=[&](long double alpha){
        return score(alpha,true);
    };
    const pair<long double,long double> maximum=joint_concave_maximum(
        max_alpha,derivatives,balanced_score);
    result.alpha=maximum.first;
    result.balanced=maximum.second;
    result.balanced /= (long double)indices.size();
    result.raw = score(result.alpha,false);
    return result;
}

static long double joint_linked_unit_log_likelihood(
        const JointLinkedUnit& unit, long double alpha, long double rho,
        long double e_ref, long double e_alt){
    if (unit.sites.empty()) return NAN;
    AxisKahan total;
    for (const JointSiteUnit& site : unit.sites)
        total.add(joint_site_log_likelihood(
            site,alpha,rho,e_ref,e_alt,true));
    return total.value/(long double)unit.sites.size();
}

static JointFit joint_fit_linked_units(
        const vector<JointLinkedUnit>& units, const vector<size_t>& indices,
        long double rho, long double e_ref, long double e_alt,
        long double max_alpha){
    JointFit result;
    if (indices.empty()) return result;
    const auto score=[&](long double alpha){
        AxisKahan total;
        for (size_t index : indices)
            total.add(joint_linked_unit_log_likelihood(
                units[index],alpha,rho,e_ref,e_alt));
        return total.value/(long double)indices.size();
    };
    const auto derivatives=[&](long double alpha){
        AxisKahan first,second;
        for (size_t index : indices){
            AxisKahan unit_first,unit_second;
            for (const JointSiteUnit& site : units[index].sites){
                const pair<long double,long double> value=
                    joint_site_log_likelihood_derivatives(
                        site,alpha,rho,e_ref,e_alt,true);
                unit_first.add(value.first); unit_second.add(value.second);
            }
            const long double count=(long double)units[index].sites.size();
            first.add(unit_first.value/count);
            second.add(unit_second.value/count);
        }
        const long double count=(long double)indices.size();
        return make_pair(first.value/count,second.value/count);
    };
    const pair<long double,long double> maximum=joint_concave_maximum(
        max_alpha,derivatives,score);
    result.alpha=maximum.first;
    result.balanced=maximum.second;
    result.raw=result.balanced*(long double)indices.size();
    return result;
}

static uint64_t joint_linked_fold_hash(
        const JointHypothesis& hypothesis, const JointLinkedUnit& unit){
    return stable_text_hash(
        hypothesis.library+"|"+hypothesis.barcode+"|"+
        unit.parent_origin+"|"+
        joint_molecule_fold_basis_name(unit.basis)+"|"+
        (unit.molecule_id.empty() ? to_string(unit.molecule) :
         unit.molecule_id)+"|"+JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION+
        "|FOLD_V1");
}

static string joint_summary_integers(vector<int> values){
    if (values.empty()) return "NA";
    sort(values.begin(),values.end());
    const auto quantile=[&](long double p){
        const size_t index=(size_t)floorl(
            p*(long double)(values.size()-1));
        return values[index];
    };
    return string("min=")+to_string(values.front())+
        ";q1="+to_string(quantile(0.25L))+
        ";median="+to_string(quantile(0.5L))+
        ";q3="+to_string(quantile(0.75L))+
        ";max="+to_string(values.back());
}

static void joint_evaluate_molecules(
        const JointHypothesis& hypothesis,
        const vector<JointMoleculeRecord>& records,
        const vector<AxisSiteDefinition>& sites,
        const unordered_map<int,size_t>& donor_slot,
        long double e_ref, long double e_alt, long double max_alpha,
        long malformed_rows, bool capture_cache, JointResult& result,
        bool build_only=false){
    result.molecule_malformed_rows=malformed_rows;
    if (!hypothesis.locked.valid || !hypothesis.second.valid){
        result.molecule_status="UNAVAILABLE";
        result.molecule_unusable_reason="DONOR_NOT_IN_MODALITY_PANEL";
        return;
    }
    vector<JointLinkedUnit> units;
    size_t cursor=0;
    while (cursor<records.size()){
        const uint64_t molecule=records[cursor].molecule;
        const uint8_t basis=records[cursor].basis;
        JointLinkedUnit linked;
        linked.molecule=molecule; linked.molecule_id=to_string(molecule);
        linked.basis=basis;
        while (cursor<records.size() && records[cursor].molecule==molecule &&
                records[cursor].basis==basis){
            const JointMoleculeRecord& record=records[cursor++];
            const uint64_t key=site_key(record.tid,record.pos);
            auto found=lower_bound(sites.begin(),sites.end(),key,
                [](const AxisSiteDefinition& site,uint64_t value){
                    return site.key<value;
                });
            if (found==sites.end() || found->key!=key || !found->found ||
                    found->mitochondrial) continue;
            const long double depth=record.ref+record.alt;
            if (depth<=0.0L) continue;
            JointSiteUnit site;
            site.tid=record.tid; site.pos=record.pos;
            site.contig=found->contig;
            site.ref_allele=uppercase(found->ref_allele);
            site.alt_allele=uppercase(found->alt_allele);
            // Normalize within linked-unit/site so duplicated reads do not
            // increase weight.
            site.ref=record.ref/depth; site.alt=record.alt/depth;
            if (!joint_expected(hypothesis.locked,*found,donor_slot,site.q_locked) ||
                    !joint_expected(hypothesis.second,*found,donor_slot,site.q_second) ||
                    (hypothesis.rho_effective>0.0L &&
                     !joint_expected(hypothesis.ambient,*found,donor_slot,
                                     site.q_ambient,true))) continue;
            if (hypothesis.rho_effective==0.0L) site.q_ambient=site.q_locked;
            linked.sites.push_back(site);
        }
        if (!linked.sites.empty()) units.push_back(linked);
    }
    result.molecule_units=(int)units.size();
    if (capture_cache) result.cache_molecule_units=units;
    if (build_only){
        result.molecule_status=units.empty() ? "UNAVAILABLE" :
            "UNCOMPUTED_NONEMPTY";
        result.molecule_unusable_reason=units.empty() ?
            (malformed_rows>0 ?
             "NO_USABLE_LINKED_UNITS_AFTER_MALFORMED_ROWS" :
             "NO_USABLE_LINKED_UNITS") : "NONE";
        return;
    }
    if (units.empty()){
        result.molecule_status="UNAVAILABLE";
        result.molecule_unusable_reason=malformed_rows>0 ?
            "NO_USABLE_LINKED_UNITS_AFTER_MALFORMED_ROWS" :
            "NO_USABLE_LINKED_UNITS";
        return;
    }
    vector<int> snps_per_unit;
    int umi_gene=0,qname=0;
    vector<size_t> informative;
    for (size_t index=0;index<units.size();++index){
        JointLinkedUnit& unit=units[index];
        snps_per_unit.push_back((int)unit.sites.size());
        result.molecule_total_snps+=(int)unit.sites.size();
        if (unit.sites.size()>1) ++result.molecule_multi_snp_units;
        if (unit.basis==1 || unit.basis==2) ++umi_gene;
        if (unit.basis==3) ++qname;
        bool distinguishes=false;
        for (const JointSiteUnit& site : unit.sites)
            if (fabsl(site.q_locked-site.q_second)>1e-18L){
                distinguishes=true; break;
            }
        if (distinguishes) informative.push_back(index);
    }
    result.molecule_discriminating_units=(int)informative.size();
    result.molecule_snps_per_unit=joint_summary_integers(snps_per_unit);
    result.molecule_umi_gene_fraction=(long double)umi_gene/(long double)units.size();
    result.molecule_query_name_fraction=(long double)qname/(long double)units.size();
    if (informative.empty()){
        result.molecule_status="GENOTYPE_EQUIVALENT";
        result.molecule_unusable_reason="NO_GENOTYPE_DISTINGUISHABLE_LINKED_UNITS";
        result.molecule_alpha=0.0L;
        result.preferred_molecule_model=capture_cache ?
            "LOCKED_IDENTITY" : "TARGETED_ONLY_NOT_REQUESTED";
        return;
    }
    const JointFit fit=joint_fit_linked_units(
        units,informative,hypothesis.rho_effective,e_ref,e_alt,max_alpha);
    result.molecule_alpha=fit.alpha;
    result.molecule_k2=fit.balanced;
    AxisKahan frozen;
    for (size_t index : informative)
        frozen.add(joint_linked_unit_log_likelihood(
            units[index],0.0L,hypothesis.rho_effective,e_ref,e_alt));
    result.molecule_k1=frozen.value/(long double)informative.size();
    result.molecule_delta=result.molecule_k2-result.molecule_k1;

    const long double cutoff=result.molecule_k2-
        1.920729410347062L/(long double)informative.size();
    const auto profile_score=[&](long double alpha){
        AxisKahan score;
        for (size_t index : informative)
            score.add(joint_linked_unit_log_likelihood(
                units[index],alpha,hypothesis.rho_effective,e_ref,e_alt));
        return score.value/(long double)informative.size();
    };
    const pair<long double,long double> profile=joint_profile_interval(
        max_alpha,result.molecule_alpha,cutoff,profile_score);
    result.molecule_alpha_low=profile.first;
    result.molecule_alpha_high=profile.second;
    if (capture_cache){
        result.molecule_contributor_only=profile_score(1.0L);
        result.preferred_molecule_model="LOCKED_IDENTITY";
        long double preferred=result.molecule_k1;
        if (result.molecule_alpha>1e-12L && result.molecule_k2>preferred){
            preferred=result.molecule_k2;
            result.preferred_molecule_model="INTERIOR_LOCKED_PLUS_CONTRIBUTOR";
        }
        if (result.molecule_contributor_only>preferred){
            result.preferred_molecule_model="CONTRIBUTOR_ONLY_FRACTION_1";
        }
    } else {
        result.preferred_molecule_model="TARGETED_ONLY_NOT_REQUESTED";
    }

    vector<long double> influences(informative.size());
    AxisKahan absolute_sum,squared_sum;
    size_t top_position=0;
    for (size_t i=0;i<informative.size();++i){
        const JointLinkedUnit& unit=units[informative[i]];
        const long double signed_delta=joint_linked_unit_log_likelihood(
            unit,result.molecule_alpha,hypothesis.rho_effective,e_ref,e_alt)-
            joint_linked_unit_log_likelihood(
                unit,0.0L,hypothesis.rho_effective,e_ref,e_alt);
        influences[i]=fabsl(signed_delta);
        absolute_sum.add(influences[i]);
        squared_sum.add(influences[i]*influences[i]);
        if (influences[i]>influences[top_position]) top_position=i;
    }
    if (squared_sum.value>0.0L)
        result.molecule_effective_units=
            absolute_sum.value*absolute_sum.value/squared_sum.value;
    if (absolute_sum.value>0.0L)
        result.molecule_maximum_influence_fraction=
            influences[top_position]/absolute_sum.value;
    if (informative.size()>1){
        vector<size_t> without_top;
        for (size_t i=0;i<informative.size();++i)
            if (i!=top_position) without_top.push_back(informative[i]);
        const JointFit reduced=joint_fit_linked_units(
            units,without_top,hypothesis.rho_effective,e_ref,e_alt,max_alpha);
        AxisKahan reduced_frozen;
        for (size_t index : without_top)
            reduced_frozen.add(joint_linked_unit_log_likelihood(
                units[index],0.0L,hypothesis.rho_effective,e_ref,e_alt));
        result.molecule_without_top_alpha=reduced.alpha;
        result.molecule_without_top_delta=reduced.balanced-
            reduced_frozen.value/(long double)without_top.size();
    }

    const int fold_count=5;
    result.molecule_fold_count=fold_count;
    vector<size_t> fold_order(units.size());
    iota(fold_order.begin(),fold_order.end(),0);
    vector<int> fold_sizes(fold_count,0);
    for (size_t rank=0;rank<fold_order.size();++rank){
        const int fold=(int)(joint_linked_fold_hash(
            hypothesis,units[fold_order[rank]])%5ULL);
        units[fold_order[rank]].fold=fold;
        ++fold_sizes[fold];
    }
    result.molecule_fold_units=joint_summary_integers(fold_sizes);
    vector<long double> heldout_unit_deltas,fold_means,fold_alphas;
    int positive_folds=0;
    if (informative.size()>=2 && fold_count>=2){
        for (int fold=0;fold<fold_count;++fold){
            vector<size_t> training,heldout;
            for (size_t index : informative){
                if (units[index].fold==fold) heldout.push_back(index);
                else training.push_back(index);
            }
            if (training.empty() || heldout.empty()) continue;
            const JointFit fold_fit=joint_fit_linked_units(
                units,training,hypothesis.rho_effective,e_ref,e_alt,max_alpha);
            AxisKahan fold_total;
            for (size_t index : heldout){
                const long double delta=joint_linked_unit_log_likelihood(
                    units[index],fold_fit.alpha,hypothesis.rho_effective,
                    e_ref,e_alt)-joint_linked_unit_log_likelihood(
                    units[index],0.0L,hypothesis.rho_effective,e_ref,e_alt);
                heldout_unit_deltas.push_back(delta);
                fold_total.add(delta);
            }
            const long double fold_mean=fold_total.value/(long double)heldout.size();
            fold_means.push_back(fold_mean);
            fold_alphas.push_back(fold_fit.alpha);
            if (fold_mean>0.0L) ++positive_folds;
        }
    }
    result.molecule_folds_evaluable=(int)fold_means.size();
    if (heldout_unit_deltas.size()>=2 && fold_means.size()>=2){
        AxisKahan heldout_total;
        for (long double value : heldout_unit_deltas) heldout_total.add(value);
        result.molecule_heldout_delta=heldout_total.value/
            (long double)heldout_unit_deltas.size();
        result.molecule_heldout_support_fraction=(long double)positive_folds/
            (long double)fold_means.size();
        sort(fold_means.begin(),fold_means.end());
        result.molecule_fold_min_delta=fold_means.front();
        result.molecule_fold_median_delta=fold_means[fold_means.size()/2];
        string alpha_text;
        for (size_t i=0;i<fold_alphas.size();++i){
            if (i) alpha_text+=",";
            alpha_text+=axis_fmt(fold_alphas[i]);
        }
        result.molecule_fold_fitted_fractions=alpha_text;
    }
    result.molecule_status=informative.size()<2 ? "LIMITED_EVIDENCE" : "AVAILABLE";
    result.molecule_unusable_reason=informative.size()<2 ?
        "FEWER_THAN_TWO_DISCRIMINATING_LINKED_UNITS" : "NONE";
    vector<string> warnings;
    if (malformed_rows>0)
        warnings.push_back("MALFORMED_CELL_MOLECULE_ROWS_SKIPPED:"+
            to_string(malformed_rows));
    if (isfinite(result.molecule_maximum_influence_fraction) &&
            result.molecule_maximum_influence_fraction>0.5L)
        warnings.push_back("TOP_LINKED_UNIT_DOMINATES_ABSOLUTE_INFLUENCE");
    if (result.molecule_folds_evaluable<2)
        warnings.push_back("HELD_OUT_LINKED_UNIT_SCORE_UNAVAILABLE");
    result.molecule_warnings=warnings.empty() ? "NONE" : join_flags(warnings);
}

static uint64_t joint_fold_hash(const JointSiteUnit& unit){
    if (unit.contig.empty() || unit.ref_allele.empty() || unit.alt_allele.empty())
        throw runtime_error("site fold requires canonical contig/REF/ALT metadata");
    const string key=unit.contig+"|"+to_string(unit.pos)+"|"+
        uppercase(unit.ref_allele)+"|"+uppercase(unit.alt_allele);
    return stable_text_hash(key+"|"+
        JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION+"|FOLD_V1");
}

static JointResult joint_evaluate(
        const JointHypothesis& hypothesis,
        const vector<AxisObservationRecord>& observations,
        const vector<AxisSiteDefinition>& sites,
        const unordered_map<int,size_t>& donor_slot,
        long double e_ref, long double e_alt, long min_evidence,
        long double max_alpha, int requested_folds, bool capture_cache,
        bool build_only=false){
    JointResult result;
    vector<string> warnings;
    if (!hypothesis.locked.valid || !hypothesis.second.valid){
        result.status = "DONOR_NOT_IN_MODALITY_PANEL";
        result.genotype_visibility = "UNAVAILABLE_DONOR";
        if (!hypothesis.locked.valid)
            warnings.push_back("LOCKED_DONOR_MISSING:" + hypothesis.locked.missing_donors);
        if (!hypothesis.second.valid)
            warnings.push_back("SECOND_DONOR_MISSING:" + hypothesis.second.missing_donors);
        result.warnings = join_flags(warnings);
        return result;
    }
    if (hypothesis.rho_requested > 0.0L && !hypothesis.ambient.valid)
        warnings.push_back("AMBIENT_PROFILE_UNAVAILABLE_RHO_ZEROED");
    vector<JointSiteUnit> units;
    units.reserve(observations.size());
    for (const AxisObservationRecord& observation : observations){
        const uint64_t key = site_key(observation.tid,observation.pos);
        auto found = lower_bound(sites.begin(),sites.end(),key,
            [](const AxisSiteDefinition& site, uint64_t value){
                return site.key < value;
            });
        if (found == sites.end() || found->key != key || !found->found){
            ++result.excluded_missing_definition; continue;
        }
        if (found->mitochondrial){
            ++result.excluded_mitochondrial; continue;
        }
        const long double depth = observation.ref + observation.alt;
        if (depth <= 0.0L){
            ++result.excluded_nonpositive; continue;
        }
        JointSiteUnit unit;
        unit.tid=observation.tid; unit.pos=observation.pos;
        unit.contig=found->contig;
        unit.ref_allele=uppercase(found->ref_allele);
        unit.alt_allele=uppercase(found->alt_allele);
        unit.ref=observation.ref; unit.alt=observation.alt;
        if (!joint_expected(hypothesis.locked,*found,donor_slot,unit.q_locked) ||
                !joint_expected(hypothesis.second,*found,donor_slot,unit.q_second) ||
                (hypothesis.rho_effective > 0.0L &&
                 !joint_expected(hypothesis.ambient,*found,donor_slot,
                                 unit.q_ambient,true))){
            ++result.excluded_missing_genotype; continue;
        }
        if (hypothesis.rho_effective == 0.0L)
            unit.q_ambient = unit.q_locked;
        units.push_back(unit);
    }
    result.common_sites = (int)units.size();
    if (capture_cache) result.cache_site_units=units;
    if (build_only){
        result.status=units.empty() ?
            (observations.empty() ? "NO_OBSERVATIONS" :
             "NO_COMMON_NUCLEAR_GENOTYPES") : "UNCOMPUTED_NONEMPTY";
        result.genotype_visibility=units.empty() ? "UNAVAILABLE" :
            "NOT_EVALUATED";
        result.warnings=warnings.empty() ? "NONE" : join_flags(warnings);
        return result;
    }
    if (units.empty()){
        result.status = observations.empty() ? "NO_OBSERVATIONS" :
            "NO_COMMON_NUCLEAR_GENOTYPES";
        result.genotype_visibility = "UNAVAILABLE";
        result.warnings = warnings.empty() ? "NONE" : join_flags(warnings);
        return result;
    }
    vector<size_t> discriminating;
    discriminating.reserve(units.size());
    for (size_t i = 0; i < units.size(); ++i){
        if (fabsl(units[i].q_locked-units[i].q_second) > 1e-18L){
            discriminating.push_back(i);
            ++result.discriminating_sites;
            result.discriminating_depth += units[i].ref+units[i].alt;
        }
    }
    result.genotype_equivalent = result.discriminating_sites == 0;
    result.genotype_visibility = result.genotype_equivalent ?
        "GENOTYPE_EQUIVALENT_REQUIRES_OCCUPANCY" : "GENOTYPE_VISIBLE";
    if (result.genotype_equivalent){
        AxisKahan baseline_balanced, baseline_raw;
        for (const JointSiteUnit& unit : units){
            baseline_balanced.add(joint_site_log_likelihood(
                unit,0.0L,hypothesis.rho_effective,e_ref,e_alt,true));
            baseline_raw.add(joint_site_log_likelihood(
                unit,0.0L,hypothesis.rho_effective,e_ref,e_alt,false));
        }
        result.alpha=0.0L;
        result.k1_balanced=baseline_balanced.value/(long double)units.size();
        result.k2_balanced=result.k1_balanced;
        result.delta_balanced=0.0L;
        result.k1_raw=baseline_raw.value;
        result.k2_raw=result.k1_raw;
        result.delta_raw=0.0L;
        if (capture_cache){
            result.contributor_only_balanced=result.k1_balanced;
            result.contributor_only_raw=result.k1_raw;
            result.preferred_site_model="LOCKED_IDENTITY";
        } else {
            result.preferred_site_model="TARGETED_ONLY_NOT_REQUESTED";
        }
        result.status="GENOTYPE_EQUIVALENT";
        result.warnings=warnings.empty() ? "NONE" : join_flags(warnings);
        return result;
    }
    const JointFit fitted = joint_fit(units,discriminating,hypothesis.rho_effective,
        e_ref,e_alt,max_alpha);
    result.alpha=fitted.alpha;
    result.k2_balanced=fitted.balanced;
    result.k2_raw=fitted.raw;
    AxisKahan k1_balanced_sum, k1_raw_sum;
    for (size_t index : discriminating){
        const JointSiteUnit& unit=units[index];
        k1_balanced_sum.add(joint_site_log_likelihood(
            unit,0.0L,hypothesis.rho_effective,e_ref,e_alt,true));
        k1_raw_sum.add(joint_site_log_likelihood(
            unit,0.0L,hypothesis.rho_effective,e_ref,e_alt,false));
    }
    result.k1_balanced=k1_balanced_sum.value/(long double)discriminating.size();
    result.k1_raw=k1_raw_sum.value;
    result.delta_balanced=result.k2_balanced-result.k1_balanced;
    result.delta_raw=result.k2_raw-result.k1_raw;
    result.replacement_like=result.alpha > 0.5L;

    // A likelihood-profile interval on the site-balanced likelihood scale.
    const long double cutoff = result.k2_balanced-
        1.920729410347062L/(long double)discriminating.size();
    const auto profile_score=[&](long double alpha){
        AxisKahan value;
        for (size_t index : discriminating)
            value.add(joint_site_log_likelihood(
                units[index],alpha,hypothesis.rho_effective,e_ref,e_alt,true));
        return value.value/(long double)discriminating.size();
    };
    const pair<long double,long double> profile=joint_profile_interval(
        max_alpha,result.alpha,cutoff,profile_score);
    result.alpha_low=profile.first;
    result.alpha_high=profile.second;
    if (capture_cache){
        result.contributor_only_balanced=profile_score(1.0L);
        AxisKahan contributor_raw;
        for (size_t index : discriminating)
            contributor_raw.add(joint_site_log_likelihood(
                units[index],1.0L,hypothesis.rho_effective,e_ref,e_alt,false));
        result.contributor_only_raw=contributor_raw.value;
        result.preferred_site_model="LOCKED_IDENTITY";
        long double preferred=result.k1_balanced;
        if (result.alpha>1e-12L && result.k2_balanced>preferred){
            preferred=result.k2_balanced;
            result.preferred_site_model="INTERIOR_LOCKED_PLUS_CONTRIBUTOR";
        }
        if (result.contributor_only_balanced>preferred){
            result.preferred_site_model="CONTRIBUTOR_ONLY_FRACTION_1";
        }
    } else {
        result.preferred_site_model="TARGETED_ONLY_NOT_REQUESTED";
    }

    vector<long double> absolute_site_delta;
    AxisKahan absolute_total;
    for (const JointSiteUnit& unit : units){
        const long double delta = joint_site_log_likelihood(
            unit,result.alpha,hypothesis.rho_effective,e_ref,e_alt,true) -
            joint_site_log_likelihood(
                unit,0.0L,hypothesis.rho_effective,e_ref,e_alt,true);
        absolute_site_delta.push_back(fabsl(delta));
        absolute_total.add(fabsl(delta));
    }
    if (absolute_total.value > 0.0L)
        result.top_site_fraction=*max_element(
            absolute_site_delta.begin(),absolute_site_delta.end())/
            absolute_total.value;

    const int folds=5;
    vector<long double> fold_deltas;
    int supporting=0;
    if (discriminating.size() >= 2){
        for (int fold = 0; fold < folds; ++fold){
            vector<size_t> training;
            for (size_t index : discriminating)
                if ((int)(joint_fold_hash(units[index])%(uint64_t)folds) != fold)
                    training.push_back(index);
            if (training.empty() || training.size() == discriminating.size()) continue;
            JointFit fold_fit=joint_fit(units,training,hypothesis.rho_effective,
                e_ref,e_alt,max_alpha);
            AxisKahan fold_k1;
            for (size_t index : training)
                fold_k1.add(joint_site_log_likelihood(
                    units[index],0.0L,hypothesis.rho_effective,e_ref,e_alt,true));
            const long double delta=fold_fit.balanced-
                fold_k1.value/(long double)training.size();
            fold_deltas.push_back(delta);
            if (fold_fit.alpha > 0.01L && delta > 0.0L) ++supporting;
        }
    }
    result.folds_evaluable=(int)fold_deltas.size();
    if (!fold_deltas.empty()){
        sort(fold_deltas.begin(),fold_deltas.end());
        result.fold_min_delta=fold_deltas.front();
        result.fold_median_delta=fold_deltas[fold_deltas.size()/2];
        result.fold_support_fraction=(long double)supporting/
            (long double)fold_deltas.size();
    }

    if (result.discriminating_sites < 2 ||
            result.discriminating_depth < min_evidence)
        result.status="LOW_EVIDENCE";
    else result.status="AVAILABLE";
    if (result.replacement_like)
        warnings.push_back("FITTED_SECOND_FRACTION_ABOVE_ONE_HALF");
    if (isfinite(result.top_site_fraction) && result.top_site_fraction > 0.5L)
        warnings.push_back("TOP_SITE_DOMINATES_BALANCED_DELTA");
    if (result.folds_evaluable > 0 && result.fold_support_fraction < 0.8L)
        warnings.push_back("LEAVE_ONE_FOLD_OUT_SUPPORT_UNSTABLE");
    result.warnings=warnings.empty() ? "NONE" : join_flags(warnings);
    return result;
}

static vector<string> joint_output_header(){
    return {
        "schema_version","scientific_method_version","library","barcode","modality","candidate_id",
        "locked_state","locked_copy_vector","second_state",
        "second_copy_vector","candidate_origin","exhaustive_fallback",
        "nomination_modalities","rho_requested","rho_effective",
        "ambient_copy_vector","ambient_status","score_status",
        "genotype_visibility","genotype_equivalent","replacement_like",
        "fitted_second_fraction","fitted_second_fraction_profile_low",
        "fitted_second_fraction_profile_high",
        "k1_site_balanced_log_likelihood",
        "k2_site_balanced_log_likelihood",
        "delta_site_balanced_log_likelihood_k2_minus_k1",
        "contributor_only_fraction1_site_balanced_log_likelihood",
        "preferred_site_model_interpretation",
        "k1_raw_log_likelihood","k2_raw_log_likelihood",
        "delta_raw_log_likelihood_k2_minus_k1",
        "contributor_only_fraction1_raw_log_likelihood",
        "n_common_nuclear_sites","n_discriminating_sites",
        "discriminating_depth","n_leave_one_fold_out_evaluable",
        "leave_one_fold_out_support_fraction",
        "minimum_leave_one_fold_out_balanced_delta",
        "median_leave_one_fold_out_balanced_delta",
        "maximum_single_site_absolute_balanced_delta_fraction",
        "n_excluded_missing_site_definition","n_excluded_mitochondrial",
        "n_excluded_nonpositive","n_excluded_missing_genotype",
        "error_ref","error_alt","min_evidence","max_second_fraction",
        "site_fold_count_requested","formula_version","fold_version",
        "observation_bucket_count","target_observation_rows",
        "unique_selected_site_keys","peak_bucket_rows","warnings"
        ,"molecule_score_status","molecule_unusable_reason",
        "primary_evidence_basis","n_independent_linked_units",
        "n_discriminating_linked_units","effective_linked_unit_count",
        "n_multi_snp_linked_units","multi_snp_linked_unit_fraction",
        "total_snps_in_linked_units","snps_per_linked_unit_summary",
        "maximum_single_linked_unit_absolute_contribution_fraction",
        "molecule_balanced_k1_log_likelihood",
        "molecule_balanced_k2_log_likelihood",
        "molecule_balanced_delta_log_likelihood_k2_minus_k1",
        "molecule_balanced_contributor_only_fraction1_log_likelihood",
        "preferred_molecule_model_interpretation",
        "molecule_balanced_fitted_second_fraction",
        "molecule_balanced_fitted_second_fraction_profile_low",
        "molecule_balanced_fitted_second_fraction_profile_high",
        "molecule_fold_count","molecule_folds_evaluable",
        "molecule_heldout_equal_unit_mean_delta",
        "molecule_heldout_fold_support_fraction",
        "molecule_heldout_minimum_fold_mean_delta",
        "molecule_heldout_median_fold_mean_delta",
        "molecule_fold_unit_count_summary",
        "molecule_fold_fitted_second_fractions",
        "molecule_without_top_unit_delta",
        "molecule_without_top_unit_fitted_second_fraction",
        "rna_umi_gene_basis_fraction","query_name_fallback_basis_fraction",
        "molecule_malformed_rows","molecule_formula_version",
        "molecule_sidecar_schema_version","molecule_fold_version",
        "molecule_warnings"
    };
}

static vector<string> joint_output_row(
        const JointHypothesis& hypothesis, const JointResult& result,
        const string& modality, const AxisResourceAudit& audit,
        long double e_ref, long double e_alt, long min_evidence,
        long double max_alpha, int folds){
    return {
        JOINT_DOUBLET_SCHEMA,JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION,
        hypothesis.library,hypothesis.barcode,modality,
        hypothesis.candidate_id,hypothesis.locked_state,hypothesis.locked.text,
        hypothesis.second_state,hypothesis.second.text,
        hypothesis.candidate_origin,hypothesis.exhaustive_fallback,
        hypothesis.nomination_modalities,axis_fmt(hypothesis.rho_requested),
        axis_fmt(hypothesis.rho_effective),
        hypothesis.ambient.text.empty() ? "NA" : hypothesis.ambient.text,
        hypothesis.ambient_status,result.status,result.genotype_visibility,
        bool_text(result.genotype_equivalent),bool_text(result.replacement_like),
        axis_fmt(result.alpha),axis_fmt(result.alpha_low),axis_fmt(result.alpha_high),
        axis_fmt(result.k1_balanced),axis_fmt(result.k2_balanced),
        axis_fmt(result.delta_balanced),axis_fmt(result.contributor_only_balanced),
        result.preferred_site_model,axis_fmt(result.k1_raw),
        axis_fmt(result.k2_raw),axis_fmt(result.delta_raw),
        axis_fmt(result.contributor_only_raw),
        to_string(result.common_sites),to_string(result.discriminating_sites),
        axis_fmt(result.discriminating_depth),to_string(result.folds_evaluable),
        axis_fmt(result.fold_support_fraction),axis_fmt(result.fold_min_delta),
        axis_fmt(result.fold_median_delta),axis_fmt(result.top_site_fraction),
        to_string(result.excluded_missing_definition),
        to_string(result.excluded_mitochondrial),
        to_string(result.excluded_nonpositive),
        to_string(result.excluded_missing_genotype),axis_fmt(e_ref),
        axis_fmt(e_alt),to_string(min_evidence),axis_fmt(max_alpha),
        to_string(folds),JOINT_DOUBLET_FORMULA,JOINT_DOUBLET_FOLD_VERSION,
        to_string(audit.bucket_count),to_string(audit.target_rows),
        to_string(audit.unique_site_keys),
        to_string(audit.observed_peak_bucket_rows),result.warnings,
        result.molecule_status,result.molecule_unusable_reason,
        "SITE_AND_MOLECULE_SEPARATE_SENSITIVITIES",
        to_string(result.molecule_units),
        to_string(result.molecule_discriminating_units),
        axis_fmt(result.molecule_effective_units),
        to_string(result.molecule_multi_snp_units),
        result.molecule_units>0 ? axis_fmt(
            (long double)result.molecule_multi_snp_units/
            (long double)result.molecule_units) : "NA",
        to_string(result.molecule_total_snps),result.molecule_snps_per_unit,
        axis_fmt(result.molecule_maximum_influence_fraction),
        axis_fmt(result.molecule_k1),axis_fmt(result.molecule_k2),
        axis_fmt(result.molecule_delta),axis_fmt(result.molecule_contributor_only),
        result.preferred_molecule_model,axis_fmt(result.molecule_alpha),
        axis_fmt(result.molecule_alpha_low),axis_fmt(result.molecule_alpha_high),
        to_string(result.molecule_fold_count),
        to_string(result.molecule_folds_evaluable),
        axis_fmt(result.molecule_heldout_delta),
        axis_fmt(result.molecule_heldout_support_fraction),
        axis_fmt(result.molecule_fold_min_delta),
        axis_fmt(result.molecule_fold_median_delta),
        result.molecule_fold_units,result.molecule_fold_fitted_fractions,
        axis_fmt(result.molecule_without_top_delta),
        axis_fmt(result.molecule_without_top_alpha),
        axis_fmt(result.molecule_umi_gene_fraction),
        axis_fmt(result.molecule_query_name_fraction),
        to_string(result.molecule_malformed_rows),
        JOINT_DOUBLET_MOLECULE_FORMULA,
        JOINT_DOUBLET_MOLECULE_SIDECAR_SCHEMA,
        JOINT_DOUBLET_MOLECULE_FOLD_VERSION,result.molecule_warnings
    };
}

static void joint_write_targeted_cache(
        const string& output_path, const vector<JointHypothesis>& hypotheses,
        const vector<JointResult>& results, const string& modality,
        long double e_ref, long double e_alt){
    if (output_path.empty()) return;
    if (hypotheses.size()!=results.size())
        throw runtime_error("targeted cache hypothesis/result size mismatch");
    const string temporary_output=output_path+".tmp."+
        to_string((long long)getpid());
    gzFile output=gzopen(temporary_output.c_str(),"wb");
    if (!output)
        throw runtime_error("could not create targeted evidence cache: "+
            temporary_output);
    const vector<string> header={
        "schema_version","scientific_method_version","calibration_library",
        "library","barcode","modality","candidate_id",
        "locked_state","second_state","evidence_basis","linked_unit_id",
        "tid","pos","contig","site_ref_allele","site_alt_allele",
        "ref","alt","q_locked","q_second","q_ambient",
        "rho","error_ref","error_alt"
    };
    try {
        axis_gzwrite(output,header);
        vector<size_t> order(hypotheses.size());
        iota(order.begin(),order.end(),0);
        sort(order.begin(),order.end(),[&](size_t left,size_t right){
            return make_pair(hypotheses[left].barcode,hypotheses[left].candidate_id)<
                make_pair(hypotheses[right].barcode,hypotheses[right].candidate_id);
        });
        for (size_t index : order){
            const JointHypothesis& hypothesis=hypotheses[index];
            const JointResult& result=results[index];
            axis_gzwrite(output,{
                JOINT_DOUBLET_TARGETED_CACHE_SCHEMA,
                JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION,"25",hypothesis.library,
                hypothesis.barcode,modality,hypothesis.candidate_id,
                hypothesis.locked_state,hypothesis.second_state,
                "CANDIDATE_METADATA","NA","-1","-1","NA","NA","NA","NA","NA",
                "NA","NA","NA",axis_fmt(hypothesis.rho_effective),
                axis_fmt(e_ref),axis_fmt(e_alt)
            });
            for (const JointSiteUnit& site : result.cache_site_units){
                axis_gzwrite(output,{
                    JOINT_DOUBLET_TARGETED_CACHE_SCHEMA,
                    JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION,"25",hypothesis.library,
                    hypothesis.barcode,modality,hypothesis.candidate_id,
                    hypothesis.locked_state,hypothesis.second_state,"SITE",
                    string("SITE:")+to_string(site.tid)+":"+to_string(site.pos),
                    to_string(site.tid),to_string(site.pos),site.contig,
                    site.ref_allele,site.alt_allele,axis_fmt(site.ref),
                    axis_fmt(site.alt),axis_fmt(site.q_locked),
                    axis_fmt(site.q_second),axis_fmt(site.q_ambient),
                    axis_fmt(hypothesis.rho_effective),axis_fmt(e_ref),
                    axis_fmt(e_alt)
                });
            }
            for (const JointLinkedUnit& unit : result.cache_molecule_units){
                const string basis=string("MOLECULE_")+
                    joint_molecule_basis_name(unit.basis);
                const string unit_id=basis+":"+to_string(unit.molecule);
                for (const JointSiteUnit& site : unit.sites){
                    axis_gzwrite(output,{
                        JOINT_DOUBLET_TARGETED_CACHE_SCHEMA,
                        JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION,"25",hypothesis.library,
                        hypothesis.barcode,modality,hypothesis.candidate_id,
                        hypothesis.locked_state,hypothesis.second_state,basis,
                        unit_id,to_string(site.tid),to_string(site.pos),site.contig,
                        site.ref_allele,site.alt_allele,
                        axis_fmt(site.ref),axis_fmt(site.alt),
                        axis_fmt(site.q_locked),axis_fmt(site.q_second),
                        axis_fmt(site.q_ambient),
                        axis_fmt(hypothesis.rho_effective),axis_fmt(e_ref),
                        axis_fmt(e_alt)
                    });
                }
            }
        }
        if (gzclose(output)!=Z_OK)
            throw runtime_error("failed closing targeted evidence cache: "+
                temporary_output);
        output=NULL;
    } catch (...) {
        if (output) gzclose(output);
        unlink(temporary_output.c_str());
        throw;
    }
    if (rename(temporary_output.c_str(),output_path.c_str())!=0){
        unlink(temporary_output.c_str());
        throw runtime_error("failed publishing targeted evidence cache: "+
            output_path);
    }
}

// -------------------------------------------------------------------------
// Evidence-normalized indexed cache and compiled targeted validation (v4)
// -------------------------------------------------------------------------

static const char* JOINT_NORMALIZED_CACHE_SCHEMA =
    "joint_doublet_normalized_indexed_cache_v5";
static const char* JOINT_TARGETED_ANALYSIS_SCHEMA =
    "joint_doublet_compiled_targeted_analysis_v5";

struct JointNormalizedCellIndex {
    uint64_t observation_offset = 0;
    uint64_t observation_count = 0;
    uint64_t molecule_offset = 0;
    uint64_t molecule_count = 0;
    long malformed_molecule_rows = 0;
};

struct JointBucketPlan {
    unordered_map<unsigned long,size_t> assignment;
    vector<unsigned long long> measured_bytes;
    vector<unsigned long long> projected_peak_bytes;
    unsigned long long limit_bytes = 256ULL*1024ULL*1024ULL;
};

static JointBucketPlan joint_plan_normalized_buckets(
        const unordered_map<unsigned long,vector<size_t>>& by_cell,
        const unordered_map<unsigned long,unsigned long long>& observation_rows,
        const unordered_map<unsigned long,unsigned long long>& molecule_rows,
        int worker_threads,unsigned long long worker_memory_budget_bytes){
    // The plan is deliberately byte based.  It includes the two persistent
    // record streams plus pessimistic linked-unit boundaries, compiled menu
    // tables, thread-local workspaces, bounded null blocks, output staging and
    // a 50% safety margin.  Cells with only molecule evidence are therefore
    // first-class planning inputs.
    JointBucketPlan plan;
    if (worker_memory_budget_bytes<16ULL*1024ULL*1024ULL)
        throw runtime_error("normalized-cache worker memory budget is too small");
    plan.limit_bytes=worker_memory_budget_bytes;
    const unsigned long long threads=max(1,worker_threads);
    vector<pair<unsigned long,unsigned long long>> cells;
    for (const auto& item : by_cell){
        const unsigned long barcode=item.first;
        const unsigned long long observations=observation_rows.count(barcode) ?
            observation_rows.at(barcode) : 0ULL;
        const unsigned long long molecules=molecule_rows.count(barcode) ?
            molecule_rows.at(barcode) : 0ULL;
        // Extraction sorts fixed-width records in bounded runs. In particular,
        // a high-evidence cell may span runs without splitting linked units in
        // the published cache. Candidate fitting is not performed here.
        const unsigned long long measured=observations*sizeof(AxisObservationRecord)+
            molecules*sizeof(JointMoleculeRecord);
        const unsigned long long projected=min<unsigned long long>(
            measured*2+2ULL*1024*1024,256ULL*1024*1024);
        cells.push_back(make_pair(barcode,projected));
    }
    sort(cells.begin(),cells.end(),[](
            const pair<unsigned long,unsigned long long>& left,
            const pair<unsigned long,unsigned long long>& right){
        return left.second!=right.second ? left.second>right.second : left.first<right.first;
    });
    if (cells.empty()) throw runtime_error("normalized-cache workload has no selected cells");
    const size_t minimum=min(cells.size(),max<size_t>(1,threads*4));
    plan.projected_peak_bytes.assign(minimum,0);
    plan.measured_bytes.assign(minimum,0);
    for (const auto& item : cells){
        size_t best=0;
        for (size_t i=1;i<plan.projected_peak_bytes.size();++i)
            if (make_pair(plan.projected_peak_bytes[i],i)<
                    make_pair(plan.projected_peak_bytes[best],best)) best=i;
        if (plan.projected_peak_bytes[best]+item.second>plan.limit_bytes){
            plan.projected_peak_bytes.push_back(0);
            plan.measured_bytes.push_back(0);
            best=plan.projected_peak_bytes.size()-1;
        }
        plan.assignment[item.first]=best;
        plan.projected_peak_bytes[best]+=item.second;
        const unsigned long long observations=observation_rows.count(item.first) ?
            observation_rows.at(item.first) : 0ULL;
        const unsigned long long molecules=molecule_rows.count(item.first) ?
            molecule_rows.at(item.first) : 0ULL;
        plan.measured_bytes[best]+=
            observations*sizeof(AxisObservationRecord)+
            molecules*sizeof(JointMoleculeRecord);
    }
    if (plan.assignment.size()!=by_cell.size() ||
            *max_element(plan.projected_peak_bytes.begin(),
                         plan.projected_peak_bytes.end())>plan.limit_bytes)
        throw runtime_error("normalized-cache bucket plan violates bounded-memory contract");
    return plan;
}

static unsigned long long joint_file_bytes(const string& path){
    struct stat info;
    if (stat(path.c_str(),&info)!=0 || info.st_size<0)
        throw runtime_error("could not stat source/cache file: "+path);
    return static_cast<unsigned long long>(info.st_size);
}

static string joint_json_escape(const string& value){
    string result;
    for (unsigned char character : value){
        if (character=='"') result += "\\\"";
        else if (character=='\\') result += "\\\\";
        else if (character=='\n') result += "\\n";
        else if (character=='\r') result += "\\r";
        else if (character=='\t') result += "\\t";
        else if (character<0x20){
            char buffer[8];
            snprintf(buffer,sizeof(buffer),"\\u%04x",(unsigned int)character);
            result += buffer;
        } else result.push_back(static_cast<char>(character));
    }
    return result;
}

static vector<string> joint_load_samples_once(
        const string& path,JointStreamingDigest& digest){
    ifstream input(path.c_str(),ios::binary);
    if (!input) throw runtime_error("could not open samples source: "+path);
    string content;
    char buffer[1<<20];
    while (input){
        input.read(buffer,sizeof(buffer));
        const streamsize count=input.gcount();
        if (count>0){
            digest.update(buffer,static_cast<size_t>(count));
            content.append(buffer,static_cast<size_t>(count));
        }
    }
    if (!input.eof()) throw runtime_error("failed reading samples source: "+path);
    istringstream parsed(content);
    vector<string> samples;
    string sample;
    while (parsed>>sample) samples.push_back(sample);
    if (samples.empty()) throw runtime_error("empty samples source: "+path);
    return samples;
}

static void joint_publish_file(const string& temporary,const string& final_path){
    if (rename(temporary.c_str(),final_path.c_str())!=0){
        unlink(temporary.c_str());
        throw runtime_error("failed atomically publishing cache component: "+final_path);
    }
}

template <typename Value>
static void joint_write_binary_value(ofstream& output,const Value& value,
                                     const string& path){
    output.write(reinterpret_cast<const char*>(&value),sizeof(Value));
    if (!output) throw runtime_error("failed writing binary cache component: "+path);
}

template <typename Value>
static void joint_read_binary_value(ifstream& input,Value& value,
                                    const string& path){
    input.read(reinterpret_cast<char*>(&value),sizeof(Value));
    if (!input) throw runtime_error("truncated binary cache component: "+path);
}

static void joint_write_binary_string(
        ofstream& output,const string& value,const string& path){
    if (value.size()>numeric_limits<uint32_t>::max())
        throw runtime_error("binary cache string is too long: "+path);
    const uint32_t size=static_cast<uint32_t>(value.size());
    joint_write_binary_value(output,size,path);
    if (size) output.write(value.data(),size);
    if (!output) throw runtime_error("failed writing binary cache string: "+path);
}

static string joint_read_binary_string(ifstream& input,const string& path){
    uint32_t size=0;
    joint_read_binary_value(input,size,path);
    if (size>(1U<<20)) throw runtime_error("implausible binary cache string: "+path);
    string value(size,'\0');
    if (size) input.read(&value[0],size);
    if (!input) throw runtime_error("truncated binary cache string: "+path);
    return value;
}

static void joint_write_normalized_sites(
        const string& temporary_path,const vector<AxisSiteDefinition>& sites,
        const vector<int>& donors){
    ofstream output(temporary_path.c_str(),ios::binary);
    if (!output) throw runtime_error("could not create normalized site cache: "+temporary_path);
    char magic[24]={0};
    const string label="JDNORMCACHEV4SITES";
    memcpy(magic,label.data(),min(label.size(),sizeof(magic)));
    output.write(magic,sizeof(magic));
    const uint32_t version=4;
    const uint64_t n_sites=sites.size(),n_donors=donors.size();
    joint_write_binary_value(output,version,temporary_path);
    joint_write_binary_value(output,n_sites,temporary_path);
    joint_write_binary_value(output,n_donors,temporary_path);
    for (int donor : donors){
        const int32_t value=donor;
        joint_write_binary_value(output,value,temporary_path);
    }
    for (const AxisSiteDefinition& site : sites){
        joint_write_binary_value(output,site.tid,temporary_path);
        joint_write_binary_value(output,site.pos,temporary_path);
        const uint8_t found=site.found ? 1 : 0;
        const uint8_t mitochondrial=site.mitochondrial ? 1 : 0;
        joint_write_binary_value(output,found,temporary_path);
        joint_write_binary_value(output,mitochondrial,temporary_path);
        joint_write_binary_string(output,site.contig,temporary_path);
        joint_write_binary_string(output,uppercase(site.ref_allele),temporary_path);
        joint_write_binary_string(output,uppercase(site.alt_allele),temporary_path);
        if (site.genotype.size()!=donors.size() && site.found)
            throw runtime_error("normalized site genotype width mismatch");
        for (size_t index=0;index<donors.size();++index){
            const int8_t genotype=site.found ? site.genotype[index] : -1;
            joint_write_binary_value(output,genotype,temporary_path);
        }
    }
    output.close();
    if (!output) throw runtime_error("failed closing normalized site cache: "+temporary_path);
}

// Sort a numerical bucket in bounded runs with at most 32 open readers.
// The caller streams the final merge and preserves cell/unit identity across
// every run boundary. Source gzip streams are never rescanned here.
template<class Record,class Less>
static string joint_external_sort(const string& path,const string& temp_root,Less less){
    const size_t capacity=max<size_t>(1,(128ULL*1024*1024)/sizeof(Record));
    ifstream input(path.c_str(),ios::binary);
    if(!input)throw runtime_error("cannot open numerical sort input: "+path);
    vector<string> runs;size_t serial=0;
    while(input.peek()!=EOF){
        vector<Record> records(capacity);
        input.read(reinterpret_cast<char*>(records.data()),records.size()*sizeof(Record));
        const streamsize bytes=input.gcount();
        if(input.bad() || bytes%sizeof(Record))throw runtime_error("truncated numerical sort input");
        records.resize(static_cast<size_t>(bytes)/sizeof(Record));
        sort(records.begin(),records.end(),less);
        const string run=path+".sort."+to_string(serial++);
        ofstream out(run.c_str(),ios::binary);
        out.write(reinterpret_cast<const char*>(records.data()),records.size()*sizeof(Record));
        out.close();if(!out)throw runtime_error("failed writing numerical sort run");runs.push_back(run);
    }
    input.close();
    if(runs.empty())return path;
    while(runs.size()>1){
        vector<string> next;
        for(size_t begin=0;begin<runs.size();begin+=32){
            const size_t count=min<size_t>(32,runs.size()-begin);
            vector<ifstream> streams(count);vector<Record> heads(count);vector<bool> present(count,false);
            const string output=path+".sort."+to_string(serial++);
            ofstream out(output.c_str(),ios::binary);
            const auto advance=[&](size_t i){
                streams[i].read(reinterpret_cast<char*>(&heads[i]),sizeof(Record));
                if(streams[i].gcount()!=0 && streams[i].gcount()!=sizeof(Record))throw runtime_error("truncated sort merge run");
                if(streams[i].bad())throw runtime_error("failed sort merge read");
                present[i]=streams[i].gcount()==sizeof(Record);
            };
            for(size_t i=0;i<count;++i){streams[i].open(runs[begin+i].c_str(),ios::binary);if(!streams[i])throw runtime_error("cannot open sort run");advance(i);}
            for(;;){
                size_t best=count;
                for(size_t i=0;i<count;++i)if(present[i] && (best==count || less(heads[i],heads[best])))best=i;
                if(best==count)break;
                out.write(reinterpret_cast<const char*>(&heads[best]),sizeof(Record));advance(best);
            }
            out.close();if(!out)throw runtime_error("failed closing merged numerical run");
            for(size_t i=0;i<count;++i){streams[i].close();unlink(runs[begin+i].c_str());}
            next.push_back(output);
        }runs.swap(next);
    }
    unlink(path.c_str());return runs.front();
}

static void run_joint_doublet_normalized_extract(
        const string& samples_path,const string& manifest_path,
        const string& manifest_digest,const string& sites_path,
        const string& observations_path,const string& molecules_path,
        const string& cache_prefix,const string& workload_generation,
        const string& temp_root,const string& library,const string& modality,
        long double e_ref,long double e_alt,long min_evidence,
        long double max_alpha,int worker_threads,
        unsigned long long scheduler_memory_bytes,
        unsigned long long launcher_runtime_reserve_bytes,
        unsigned long long worker_memory_budget_bytes){
    if (cache_prefix.empty() || workload_generation.empty() ||
            manifest_digest.empty())
        throw runtime_error("normalized extraction requires cache prefix, generation, and manifest digest");
    if ((library!="lib7" && library!="lib9" && library!="lib12" && library!="lib17" && library!="lib20" && library!="lib25" && library!="lib29"))
        throw runtime_error("normalized extraction rejected non-target library: "+library);
    if (modality!="RNA" && modality!="ATAC")
        throw runtime_error("normalized extraction modality must be RNA or ATAC");
    if (worker_threads<1) throw runtime_error("normalized extraction threads must be positive");
    if (e_ref<0.0L || e_alt<0.0L || e_ref+e_alt>=1.0L)
        throw runtime_error("normalized extraction error rates are invalid");
    if (min_evidence<0 || max_alpha<=0.0L || max_alpha>1.0L)
        throw runtime_error("normalized extraction model parameters are invalid");
    const unsigned long long expected_reserve=max<unsigned long long>(
        16ULL*1024ULL*1024ULL,(scheduler_memory_bytes+9ULL)/10ULL);
    if (!scheduler_memory_bytes || scheduler_memory_bytes<=expected_reserve ||
            launcher_runtime_reserve_bytes!=expected_reserve ||
            worker_memory_budget_bytes!=scheduler_memory_bytes-
                launcher_runtime_reserve_bytes)
        throw runtime_error("normalized extraction memory-budget contract mismatch");
    JointStreamingDigest samples_digest;
    const vector<string> samples=joint_load_samples_once(
        samples_path,samples_digest);
    unordered_map<string,int> sample2idx;
    for (int index=0;index<(int)samples.size();++index){
        if (samples[index].empty() || sample2idx.count(samples[index]))
            throw runtime_error("normalized extraction samples must be unique and nonblank");
        sample2idx[samples[index]]=index;
    }
    unordered_map<unsigned long,vector<size_t>> by_cell;
    JointStreamingDigest manifest_content_digest;
    const vector<JointHypothesis> hypotheses=joint_load_manifest(
        manifest_path,sample2idx,library,by_cell,&manifest_content_digest);
    if (manifest_content_digest.text()!=manifest_digest)
        throw runtime_error("candidate manifest content digest mismatch");
    if (hypotheses.empty())
        throw runtime_error("normalized extraction manifest contains no candidates");
    unordered_map<unsigned long,CandidateAxisPair> targets;
    unordered_map<unsigned long,string> barcode_text;
    for (const JointHypothesis& hypothesis : hypotheses){
        targets[hypothesis.encoded_barcode]=CandidateAxisPair();
        barcode_text[hypothesis.encoded_barcode]=hypothesis.barcode;
    }
    vector<int> donors;
    for (const JointHypothesis& hypothesis : hypotheses){
        const JointComposition* compositions[]={
            &hypothesis.locked,&hypothesis.second,&hypothesis.ambient};
        for (const JointComposition* composition : compositions)
            for (const auto& member : composition->members)
                donors.push_back(member.first);
    }
    sort(donors.begin(),donors.end());
    donors.erase(unique(donors.begin(),donors.end()),donors.end());

    AxisTempGuard temporary(temp_root);
    AxisResourceAudit audit;
    JointStreamingDigest observations_digest,sites_digest,molecules_digest;
    unordered_map<unsigned long,unsigned long long> rows_by_barcode;
    vector<uint64_t> selected_keys;
    const string observation_spool=joint_stage_observations_one_pass(
        observations_path,targets,temporary.path(),audit,rows_by_barcode,
        selected_keys,&observations_digest);
    unordered_map<unsigned long,unsigned long long> molecule_rows_by_barcode;
    unordered_map<unsigned long,long> malformed_by_barcode;
    set<uint64_t> molecule_site_union;
    const string molecule_spool=joint_stage_molecules_one_pass(
        molecules_path,targets,selected_keys,temporary.path(),
        molecule_rows_by_barcode,malformed_by_barcode,&molecules_digest,&molecule_site_union);
    selected_keys.insert(selected_keys.end(),molecule_site_union.begin(),molecule_site_union.end());
    molecule_site_union.clear();
    sort(selected_keys.begin(),selected_keys.end());
    selected_keys.erase(unique(selected_keys.begin(),selected_keys.end()),selected_keys.end());
    vector<AxisSiteDefinition> sites=axis_load_site_definitions(
        sites_path,selected_keys,(int)samples.size(),donors,audit,&sites_digest,
        worker_memory_budget_bytes/2);
    unsigned long long dictionary_bytes=selected_keys.capacity()*sizeof(uint64_t)+sites.capacity()*sizeof(AxisSiteDefinition);
    for(const auto& site:sites)dictionary_bytes+=site.genotype.capacity()*sizeof(site.genotype[0])+
        site.contig.capacity()+site.ref_allele.capacity()+site.alt_allele.capacity();
    if(dictionary_bytes+512ULL*1024*1024>worker_memory_budget_bytes)
        throw runtime_error("selected genotype dictionary exceeds extraction memory budget");
    const JointBucketPlan bucket_plan=joint_plan_normalized_buckets(
        by_cell,rows_by_barcode,molecule_rows_by_barcode,worker_threads,
        worker_memory_budget_bytes);
    audit.bucket_count=bucket_plan.projected_peak_bytes.size();
    vector<string> observation_buckets=joint_partition_staged_observations(
        observation_spool,bucket_plan.assignment,
        bucket_plan.projected_peak_bytes.size(),temporary.path());
    vector<string> molecule_buckets=joint_partition_staged_molecules(
        molecule_spool,bucket_plan.assignment,
        bucket_plan.projected_peak_bytes.size(),temporary.path());

    const string suffix=".tmp."+to_string((long long)getpid());
    const string observations_tmp=cache_prefix+".observations.bin"+suffix;
    const string molecules_tmp=cache_prefix+".molecules.bin"+suffix;
    const string sites_tmp=cache_prefix+".sites.bin"+suffix;
    const string cells_tmp=cache_prefix+".cells.tsv"+suffix;
    const string samples_tmp=cache_prefix+".samples.tsv"+suffix;
    const string genotypes_tmp=cache_prefix+".genotypes.tsv.gz"+suffix;
    const string metadata_tmp=cache_prefix+".metadata.json"+suffix;
    ofstream observations_output(observations_tmp.c_str(),ios::binary);
    ofstream molecules_output(molecules_tmp.c_str(),ios::binary);
    if (!observations_output || !molecules_output)
        throw runtime_error("could not create normalized evidence record files");
    map<unsigned long,JointNormalizedCellIndex> index;
    uint64_t observation_offset=0,molecule_offset=0,linked_unit_count=0;
    for (const string& bucket : observation_buckets){
        const string path=joint_external_sort<AxisObservationRecord>(bucket,temporary.path(),axis_observation_less);
        ifstream input(path.c_str(),ios::binary);AxisObservationRecord row,merged;bool have=false;
        AxisKahan ref,alt;
        const auto flush=[&](){
            if(!have)return;
            merged.ref=static_cast<double>(ref.value);merged.alt=static_cast<double>(alt.value);
            JointNormalizedCellIndex& cell=index[merged.barcode];
            if(!cell.observation_count)cell.observation_offset=observation_offset;
            observations_output.write(reinterpret_cast<const char*>(&merged),sizeof(merged));
            if(!observations_output)throw runtime_error("failed writing normalized observation");
            ++cell.observation_count;++observation_offset;
        };
        while(input.read(reinterpret_cast<char*>(&row),sizeof(row))){
            if(!have || row.barcode!=merged.barcode || row.tid!=merged.tid || row.pos!=merged.pos){
                flush();merged=row;ref=AxisKahan();alt=AxisKahan();have=true;
            }ref.add(row.ref);alt.add(row.alt);
        }
        if(input.bad() || input.gcount())throw runtime_error("truncated sorted observation stream");
        flush();input.close();unlink(path.c_str());
    }
    for (const string& bucket : molecule_buckets){
        const string path=joint_external_sort<JointMoleculeRecord>(bucket,temporary.path(),joint_molecule_less);
        ifstream input(path.c_str(),ios::binary);JointMoleculeRecord row,merged,previous;bool have=false,have_previous=false;
        AxisKahan ref,alt;
        const auto flush=[&](){
            if(!have)return;
            merged.ref=static_cast<double>(ref.value);merged.alt=static_cast<double>(alt.value);
            JointNormalizedCellIndex& cell=index[merged.barcode];
            if(!cell.molecule_count)cell.molecule_offset=molecule_offset;
            if(!have_previous || previous.barcode!=merged.barcode || previous.basis!=merged.basis || previous.molecule!=merged.molecule)++linked_unit_count;
            previous=merged;have_previous=true;
            molecules_output.write(reinterpret_cast<const char*>(&merged),sizeof(merged));
            if(!molecules_output)throw runtime_error("failed writing normalized molecule");
            ++cell.molecule_count;++molecule_offset;
        };
        while(input.read(reinterpret_cast<char*>(&row),sizeof(row))){
            if(!have || row.barcode!=merged.barcode || row.basis!=merged.basis || row.molecule!=merged.molecule || row.tid!=merged.tid || row.pos!=merged.pos){
                flush();merged=row;ref=AxisKahan();alt=AxisKahan();have=true;
            }ref.add(row.ref);alt.add(row.alt);
        }
        if(input.bad() || input.gcount())throw runtime_error("truncated sorted molecule stream");
        flush();input.close();unlink(path.c_str());
    }
    observations_output.close(); molecules_output.close();
    if (!observations_output || !molecules_output)
        throw runtime_error("failed closing normalized evidence record files");
    for (const auto& item : by_cell){
        JointNormalizedCellIndex& cell=index[item.first];
        cell.malformed_molecule_rows=malformed_by_barcode.count(item.first) ?
            malformed_by_barcode.at(item.first) : 0;
    }
    joint_write_normalized_sites(sites_tmp,sites,donors);

    ofstream cell_output(cells_tmp.c_str());
    if (!cell_output) throw runtime_error("could not create normalized cell index");
    cell_output << "schema_version\tworkload_generation_id\tlibrary\tmodality\tbarcode\tencoded_barcode\tobservation_offset\tobservation_count\tmolecule_offset\tmolecule_count\tmalformed_molecule_rows\n";
    for (const auto& item : index){
        if (!barcode_text.count(item.first))
            throw runtime_error("normalized cache index contains an unknown barcode");
        const JointNormalizedCellIndex& cell=item.second;
        cell_output << JOINT_NORMALIZED_CACHE_SCHEMA << '\t' << workload_generation
            << '\t' << library << '\t' << modality << '\t'
            << barcode_text.at(item.first) << '\t' << item.first << '\t'
            << cell.observation_offset << '\t' << cell.observation_count << '\t'
            << cell.molecule_offset << '\t' << cell.molecule_count << '\t'
            << cell.malformed_molecule_rows << '\n';
    }
    cell_output.close();
    if (!cell_output) throw runtime_error("failed closing normalized cell index");

    ofstream sample_output(samples_tmp.c_str());
    if (!sample_output) throw runtime_error("could not create normalized sample table");
    sample_output << "cache_donor_slot\tsample_index\tsample\n";
    for (size_t slot=0;slot<donors.size();++slot)
        sample_output << slot << '\t' << donors[slot] << '\t'
                      << samples[donors[slot]] << '\n';
    sample_output.close();
    if (!sample_output) throw runtime_error("failed closing normalized sample table");

    gzFile genotype_output=gzopen(genotypes_tmp.c_str(),"wb");
    if (!genotype_output)
        throw runtime_error("could not create normalized genotype dictionary");
    axis_gzwrite(genotype_output,{"schema_version","library","modality",
        "tid","pos","contig","ref","alt","sample","genotype"});
    for (const AxisSiteDefinition& site : sites){
        if (!site.found || site.mitochondrial) continue;
        for (size_t slot=0;slot<donors.size();++slot)
            axis_gzwrite(genotype_output,{JOINT_NORMALIZED_CACHE_SCHEMA,library,
                modality,to_string(site.tid),to_string(site.pos),site.contig,
                site.ref_allele,site.alt_allele,samples[donors[slot]],
                to_string((int)site.genotype[slot])});
    }
    if (gzclose(genotype_output)!=Z_OK)
        throw runtime_error("failed closing normalized genotype dictionary");

    const vector<pair<string,string>> cache_components={
        {"observations_binary",observations_tmp},
        {"molecules_binary",molecules_tmp},
        {"sites_binary",sites_tmp},
        {"cells_index",cells_tmp},
        {"samples_dictionary",samples_tmp},
        {"genotypes_dictionary",genotypes_tmp}
    };
    map<string,string> cache_component_digests;
    for (const auto& component : cache_components)
        cache_component_digests[component.first]=
            joint_file_content_digest(component.second);
    uint64_t staged_record_bytes=0;
    for(const auto& item:rows_by_barcode)staged_record_bytes+=item.second*sizeof(AxisObservationRecord);
    for(const auto& item:molecule_rows_by_barcode)staged_record_bytes+=item.second*sizeof(JointMoleculeRecord);
    ofstream metadata(metadata_tmp.c_str());
    if (!metadata) throw runtime_error("could not create normalized cache metadata");
    metadata << "{\n"
        << "  \"schema_version\": \"" << JOINT_NORMALIZED_CACHE_SCHEMA << "\",\n"
        << "  \"scientific_method_version\": \""
        << JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION << "\",\n"
        << "  \"calibration_library\": 25,\n"
        << "  \"scheduler_memory_bytes\": " << scheduler_memory_bytes << ",\n"
        << "  \"launcher_runtime_reserve_bytes\": "
        << launcher_runtime_reserve_bytes << ",\n"
        << "  \"worker_memory_budget_bytes\": "
        << worker_memory_budget_bytes << ",\n"
        << "  \"workload_generation_id\": \"" << joint_json_escape(workload_generation) << "\",\n"
        << "  \"library\": \"" << joint_json_escape(library) << "\",\n"
        << "  \"modality\": \"" << joint_json_escape(modality) << "\",\n"
        << "  \"manifest_digest\": \"" << joint_json_escape(manifest_digest) << "\",\n"
        << "  \"manifest_path\": \"" << joint_json_escape(manifest_path) << "\",\n"
        << "  \"samples_path\": \"" << joint_json_escape(samples_path) << "\",\n"
        << "  \"sites_path\": \"" << joint_json_escape(sites_path) << "\",\n"
        << "  \"observations_path\": \"" << joint_json_escape(observations_path) << "\",\n"
        << "  \"molecules_path\": \"" << joint_json_escape(molecules_path) << "\",\n"
        << "  \"samples_source_bytes\": " << joint_file_bytes(samples_path) << ",\n"
        << "  \"sites_source_bytes\": " << joint_file_bytes(sites_path) << ",\n"
        << "  \"observations_source_bytes\": " << joint_file_bytes(observations_path) << ",\n"
        << "  \"molecules_source_bytes\": " << joint_file_bytes(molecules_path) << ",\n"
        << "  \"samples_content_digest\": \"" << samples_digest.text() << "\",\n"
        << "  \"sites_content_digest\": \"" << sites_digest.text() << "\",\n"
        << "  \"observations_content_digest\": \"" << observations_digest.text() << "\",\n"
        << "  \"molecules_content_digest\": \"" << molecules_digest.text() << "\",\n"
        << "  \"digest_definition\": \"FNV1A64 over decompressed bytes during the permitted source scan\",\n"
        << "  \"selected_cells\": " << by_cell.size() << ",\n"
        << "  \"candidate_rows\": " << hypotheses.size() << ",\n"
        << "  \"selected_sites\": " << sites.size() << ",\n"
        << "  \"donors\": " << donors.size() << ",\n"
        << "  \"observation_records\": " << observation_offset << ",\n"
        << "  \"molecule_records\": " << molecule_offset << ",\n"
        << "  \"linked_units\": " << linked_unit_count << ",\n"
        << "  \"staged_record_bytes_before_deduplication\": " << staged_record_bytes << ",\n"
        << "  \"temporary_numeric_record_bytes_upper_bound\": " << 3*staged_record_bytes << ",\n"
        << "  \"dictionary_memory_bytes\": " << dictionary_bytes << ",\n"
        << "  \"sort_working_memory_bound_bytes\": " << 256ULL*1024*1024 << ",\n"
        << "  \"observation_record_bytes\": " << sizeof(AxisObservationRecord) << ",\n"
        << "  \"molecule_record_bytes\": " << sizeof(JointMoleculeRecord) << ",\n"
        << "  \"observations_binary_bytes\": " << joint_file_bytes(observations_tmp) << ",\n"
        << "  \"molecules_binary_bytes\": " << joint_file_bytes(molecules_tmp) << ",\n"
        << "  \"sites_binary_bytes\": " << joint_file_bytes(sites_tmp) << ",\n"
        << "  \"cells_index_bytes\": " << joint_file_bytes(cells_tmp) << ",\n"
        << "  \"samples_dictionary_bytes\": " << joint_file_bytes(samples_tmp) << ",\n"
        << "  \"genotypes_dictionary_bytes\": " << joint_file_bytes(genotypes_tmp) << ",\n"
        << "  \"observations_binary_content_digest\": \""
        << cache_component_digests.at("observations_binary") << "\",\n"
        << "  \"molecules_binary_content_digest\": \""
        << cache_component_digests.at("molecules_binary") << "\",\n"
        << "  \"sites_binary_content_digest\": \""
        << cache_component_digests.at("sites_binary") << "\",\n"
        << "  \"cells_index_content_digest\": \""
        << cache_component_digests.at("cells_index") << "\",\n"
        << "  \"samples_dictionary_content_digest\": \""
        << cache_component_digests.at("samples_dictionary") << "\",\n"
        << "  \"genotypes_dictionary_content_digest\": \""
        << cache_component_digests.at("genotypes_dictionary") << "\",\n"
        << "  \"cache_component_digest_definition\": "
        << "\"LEGACY_PROJECT_FNV1A64 over exact published component bytes\",\n"
        << "  \"error_ref\": " << axis_fmt(e_ref) << ",\n"
        << "  \"error_alt\": " << axis_fmt(e_alt) << ",\n"
        << "  \"min_evidence\": " << min_evidence << ",\n"
        << "  \"max_second_fraction\": " << axis_fmt(max_alpha) << ",\n"
        << "  \"bucket_memory_limit_bytes\": " << bucket_plan.limit_bytes << ",\n"
        << "  \"bucket_memory_safety_margin\": 1.5,\n"
        << "  \"bucket_measured_bytes\": [";
    for (size_t i=0;i<bucket_plan.measured_bytes.size();++i){
        if (i) metadata << ",";
        metadata << bucket_plan.measured_bytes[i];
    }
    metadata << "],\n  \"bucket_projected_peak_bytes\": [";
    for (size_t i=0;i<bucket_plan.projected_peak_bytes.size();++i){
        if (i) metadata << ",";
        metadata << bucket_plan.projected_peak_bytes[i];
    }
    metadata << "],\n"
        << "  \"source_scans\": {\"sites\": 1, \"observations\": 1, \"molecules\": 1},\n"
        << "  \"atomic_publication\": \"metadata published last\"\n"
        << "}\n";
    metadata.close();
    if (!metadata) throw runtime_error("failed closing normalized cache metadata");
    if (joint_file_bytes(observations_tmp)!=observation_offset*sizeof(AxisObservationRecord) ||
            joint_file_bytes(molecules_tmp)!=molecule_offset*sizeof(JointMoleculeRecord))
        throw runtime_error("normalized cache record-size validation failed before publication");
    joint_publish_file(observations_tmp,cache_prefix+".observations.bin");
    joint_publish_file(molecules_tmp,cache_prefix+".molecules.bin");
    joint_publish_file(sites_tmp,cache_prefix+".sites.bin");
    joint_publish_file(cells_tmp,cache_prefix+".cells.tsv");
    joint_publish_file(samples_tmp,cache_prefix+".samples.tsv");
    joint_publish_file(genotypes_tmp,cache_prefix+".genotypes.tsv.gz");
    joint_publish_file(metadata_tmp,cache_prefix+".metadata.json");
}

struct JointCompiledFit {
    string status = "UNAVAILABLE";
    string unavailable_reason = "NONE";
    long double locked = NAN;
    long double interior = NAN;
    long double contributor = NAN;
    long double delta = NAN;
    long double alpha = NAN;
    long double alpha_low = NAN;
    long double alpha_high = NAN;
    string preferred = "UNAVAILABLE";
    size_t usable_units = 0;
    size_t units = 0;
    long double discriminating_depth = 0.0L;
    long double maximum_influence_fraction = NAN;
    int folds_requested = 5;
    int folds_evaluable = 0;
    int fold_support_numerator = 0;
    long double fold_support_fraction = NAN;
};

#include "joint_doublet_bounded.h"

struct JointCompiledCandidate {
    JointHypothesis hypothesis;
    JointStoredUnits site_units;
    JointStoredUnits molecule_units;
    JointCompiledFit site_reference;
    JointCompiledFit molecule_reference;
};

struct JointCompiledUniverse {
    vector<string> unit_keys;
    vector<uint64_t> unit_hashes;
    vector<string> site_keys;
    JointStoredUnits templates;
};

struct JointCachedCell {
    vector<JointCompiledCandidate> candidates;
    JointRecordView<AxisObservationRecord> raw_observations;
    JointRecordView<JointMoleculeRecord> raw_molecules;
    long malformed_molecule_rows = 0;
};

struct JointAnalysisTask {
    int task_index = -1;
    string generation;
    string cache_generation;
    string action;
    string library;
    string modality;
    string role;
    string barcode;
    string control_ids;
    string cache_prefix;
    string scientific_method_version;
    string calibration_library;
    uint64_t scheduler_memory_bytes = 0;
    uint64_t launcher_runtime_reserve_bytes = 0;
    uint64_t worker_memory_budget_bytes = 0;
    long double error_ref = NAN;
    long double error_alt = NAN;
    long min_evidence = -1;
    long double max_second_fraction = NAN;
};

struct JointAnalysisRow {
    map<string,string> value;
};

static vector<string> joint_analysis_header(){
    return {
        "schema_version","scientific_method_version","calibration_library",
        "workload_generation_id","task_index","result_key",
        "equivalence_key","result_class","action","role","library",
        "modality","barcode","control_id","control_class","evidence_channel",
        "candidate_id","second_state","engine","threads","status","fraction",
        "replicate_count","replicate_index","legal_menu_candidates","rank","winner",
        "locked_log_likelihood","interior_log_likelihood",
        "contributor_only_log_likelihood","delta_log_likelihood",
        "fitted_fraction","fitted_fraction_profile_low",
        "fitted_fraction_profile_high","preferred_model","observed_maximum_delta",
        "evidence_category","full_menu_category","category_retention_fraction",
        "full_menu_preferred_model","model_retention_fraction","model_counts",
        "fitted_fraction_p025","fitted_fraction_p975","category_counts",
        "null_median","null_maximum","null_p95","null_p99",
        "empirical_upper_tail_probability","full_menu_winner",
        "winner_retention_fraction","median_winner_fitted_fraction",
        "winner_counts","expected_contributor","expected_contributor_rank",
        "decoy_contributor","decoy_rank","expected_recovered","decoy_won",
        "recipient_units","source_units","realized_source_fraction",
        "planned_source_fraction","successful_source_units",
        "recipient_evidence_basis","source_evidence_basis",
        "requested_replicates","attempted_replicates","successful_replicates",
        "scientifically_unavailable_replicates","unavailable_replicates",
        "technical_failures","empirical_p_numerator","empirical_p_denominator",
        "decision_eligible","derived_seed_key","runner_up",
        "winner_margin","raw_delta_margin","error_ref","error_alt","min_evidence",
        "max_second_fraction","fitted_fraction_at_cap","complete_menu_unit_count",
        "candidate_evaluated_unit_count_min",
        "candidate_evaluated_unit_count_max",
        "candidate_usable_unit_count_min","candidate_usable_unit_count_max",
        "candidate_evaluated_unit_count","usable_units","evidence_units",
        "site_layout","fit_interior_eligible","profile_low_pass",
        "profile_high_pass","fraction_range_pass","fold_support_pass",
        "influence_pass","folds_requested","folds_evaluable",
        "fold_support_numerator","fold_support_fraction",
        "maximum_influence_fraction","calibrated_support_status",
        "scheduler_memory_bytes","launcher_runtime_reserve_bytes",
        "worker_memory_budget_bytes","optimized_fit_calls","reference_fit_calls",
        "optimizer_calls","likelihood_evaluations","derivative_evaluations",
        "row_scans","candidate_owned_raw_evidence_copies",
        "per_fraction_string_rehashes","per_fraction_string_sorts",
        "unchanged_80_pass_refits",
        "selection_statistic_definition","reference_definition"
    };
}

static vector<string> joint_analysis_values(
        const JointAnalysisRow& row,const vector<string>& header){
    vector<string> values;
    values.reserve(header.size());
    for (const string& field : header){
        auto found=row.value.find(field);
        values.push_back(found==row.value.end() ? "NA" : found->second);
    }
    return values;
}

static map<string,string> joint_read_named_row(
        const string& path,int requested_index,const string& index_field){
    gzFile input=gzopen(path.c_str(),"rb");
    if (!input) throw runtime_error("could not open TSV: "+path);
    char buffer[1<<20];
    if (!gzgets(input,buffer,sizeof(buffer))){
        gzclose(input); throw runtime_error("empty TSV: "+path);
    }
    string line(buffer);
    line.erase(remove(line.begin(),line.end(),'\n'),line.end());
    line.erase(remove(line.begin(),line.end(),'\r'),line.end());
    const vector<string> header=split_tsv_strict(line);
    const map<string,int> index=joint_header_index(header,path);
    auto index_column=index.find(lowercase(index_field));
    if (index_column==index.end()){
        gzclose(input); throw runtime_error(path+": missing index field "+index_field);
    }
    map<string,string> result;
    while (gzgets(input,buffer,sizeof(buffer))){
        line=buffer;
        line.erase(remove(line.begin(),line.end(),'\n'),line.end());
        line.erase(remove(line.begin(),line.end(),'\r'),line.end());
        if (line.empty()) continue;
        const vector<string> fields=split_tsv_strict(line);
        if (index_column->second>=(int)fields.size()) continue;
        const int value=stoi(trim(fields[index_column->second]));
        if (value!=requested_index) continue;
        for (size_t i=0;i<header.size();++i)
            result[lowercase(trim(header[i]))]=i<fields.size() ? trim(fields[i]) : "";
        break;
    }
    if (gzclose(input)!=Z_OK) throw runtime_error("failed closing TSV: "+path);
    if (result.empty()) throw runtime_error(path+": requested task index not found");
    return result;
}

static vector<map<string,string>> joint_read_tsv_maps(const string& path){
    gzFile input=gzopen(path.c_str(),"rb");
    if (!input) throw runtime_error("could not open TSV: "+path);
    char buffer[1<<20];
    if (!gzgets(input,buffer,sizeof(buffer))){
        gzclose(input); throw runtime_error("empty TSV: "+path);
    }
    string line(buffer);
    line.erase(remove(line.begin(),line.end(),'\n'),line.end());
    line.erase(remove(line.begin(),line.end(),'\r'),line.end());
    const vector<string> header=split_tsv_strict(line);
    vector<map<string,string>> rows;
    while (gzgets(input,buffer,sizeof(buffer))){
        line=buffer;
        line.erase(remove(line.begin(),line.end(),'\n'),line.end());
        line.erase(remove(line.begin(),line.end(),'\r'),line.end());
        if (line.empty()) continue;
        const vector<string> fields=split_tsv_strict(line);
        map<string,string> row;
        for (size_t i=0;i<header.size();++i)
            row[lowercase(trim(header[i]))]=i<fields.size() ? trim(fields[i]) : "";
        rows.push_back(row);
    }
    if (gzclose(input)!=Z_OK) throw runtime_error("failed closing TSV: "+path);
    return rows;
}

static string joint_map_value(
        const map<string,string>& row,const string& key,bool required=true){
    auto found=row.find(lowercase(key));
    if (found==row.end() || (required && found->second.empty())){
        if (required) throw runtime_error("missing required TSV value: "+key);
        return "";
    }
    return found->second;
}

static JointAnalysisTask joint_load_analysis_task(
        const string& path,int task_index){
    const map<string,string> row=joint_read_named_row(path,task_index,"task_index");
    JointAnalysisTask task;
    task.task_index=stoi(joint_map_value(row,"task_index"));
    task.generation=joint_map_value(row,"workload_generation_id");
    task.cache_generation=joint_map_value(row,"cache_generation_id");
    task.action=joint_map_value(row,"action");
    task.library=joint_map_value(row,"library");
    task.modality=joint_map_value(row,"modality");
    task.role=joint_map_value(row,"role");
    task.barcode=joint_map_value(row,"barcode",false);
    task.control_ids=joint_map_value(row,"control_ids",false);
    task.cache_prefix=joint_map_value(row,"cache_prefix");
    task.scientific_method_version=joint_map_value(
        row,"scientific_method_version");
    task.calibration_library=joint_map_value(row,"calibration_library");
    task.scheduler_memory_bytes=stoull(joint_map_value(
        row,"scheduler_memory_bytes"));
    task.launcher_runtime_reserve_bytes=stoull(joint_map_value(
        row,"launcher_runtime_reserve_bytes"));
    task.worker_memory_budget_bytes=stoull(joint_map_value(
        row,"worker_memory_budget_bytes"));
    task.error_ref=strict_ld(joint_map_value(row,"error_ref"),"task error_ref");
    task.error_alt=strict_ld(joint_map_value(row,"error_alt"),"task error_alt");
    task.min_evidence=stol(joint_map_value(row,"min_evidence"));
    task.max_second_fraction=strict_ld(
        joint_map_value(row,"max_second_fraction"),"task max_second_fraction");
    if (task.action!="CELL" && task.action!="CONTROL_BIN")
        throw runtime_error("unsupported targeted analysis action: "+task.action);
    if ((task.library!="lib7" && task.library!="lib9" && task.library!="lib12" && task.library!="lib17" && task.library!="lib20" && task.library!="lib25" && task.library!="lib29"))
        throw runtime_error("targeted analysis rejected non-target library: "+task.library);
    if (task.modality!="RNA" && task.modality!="ATAC")
        throw runtime_error("targeted analysis modality must be RNA or ATAC");
    if (task.scientific_method_version!=JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION)
        throw runtime_error("targeted analysis scientific-method mismatch");
    if (task.calibration_library!="25")
        throw runtime_error("targeted analysis requires calibration Library 25");
    const uint64_t expected_reserve=max<uint64_t>(16ULL*1024ULL*1024ULL,
        (task.scheduler_memory_bytes+9ULL)/10ULL);
    if (task.scheduler_memory_bytes==0 ||
            task.launcher_runtime_reserve_bytes!=expected_reserve ||
            task.worker_memory_budget_bytes!=task.scheduler_memory_bytes-
                task.launcher_runtime_reserve_bytes)
        throw runtime_error("targeted analysis memory-budget contract mismatch");
    if (task.error_ref<0.0L || task.error_alt<0.0L ||
            task.error_ref+task.error_alt>=1.0L || task.min_evidence<0 ||
            task.max_second_fraction<=0.0L || task.max_second_fraction>1.0L)
        throw runtime_error("targeted analysis task contains invalid model parameters");
    return task;
}

static void joint_load_normalized_samples(
        const string& prefix,unordered_map<string,int>& sample2idx,
        unordered_map<int,size_t>& donor_slot){
    const vector<map<string,string>> rows=joint_read_tsv_maps(prefix+".samples.tsv");
    for (const auto& row : rows){
        const size_t slot=static_cast<size_t>(stoull(joint_map_value(row,"cache_donor_slot")));
        const int sample_index=stoi(joint_map_value(row,"sample_index"));
        const string sample=joint_map_value(row,"sample");
        if (sample2idx.count(sample) || donor_slot.count(sample_index))
            throw runtime_error("duplicate sample in normalized cache");
        sample2idx[sample]=sample_index;
        donor_slot[sample_index]=slot;
    }
    if (sample2idx.empty()) throw runtime_error("normalized cache sample table is empty");
}

static vector<AxisSiteDefinition> joint_load_normalized_sites(
        const string& prefix,vector<int>& donors,uint64_t memory_limit){
    const string path=prefix+".sites.bin";
    ifstream input(path.c_str(),ios::binary);
    if (!input) throw runtime_error("could not open normalized site cache: "+path);
    char magic[24]={0}; input.read(magic,sizeof(magic));
    if (!input || string(magic).find("JDNORMCACHEV4SITES")!=0)
        throw runtime_error("normalized site cache magic mismatch: "+path);
    uint32_t version=0; uint64_t n_sites=0,n_donors=0;
    joint_read_binary_value(input,version,path);
    joint_read_binary_value(input,n_sites,path);
    joint_read_binary_value(input,n_donors,path);
    if (version!=4 || n_sites>100000000ULL || n_donors>1000000ULL)
        throw runtime_error("normalized site cache header is invalid: "+path);
    const uint64_t per_site=sizeof(AxisSiteDefinition)+n_donors+64;
    if (n_sites>memory_limit/per_site)
        throw runtime_error("normalized genotype dictionary exceeds task memory budget");
    uint64_t allocated=n_sites*per_site;
    donors.resize(static_cast<size_t>(n_donors));
    for (size_t i=0;i<donors.size();++i){
        int32_t donor=0; joint_read_binary_value(input,donor,path); donors[i]=donor;
    }
    vector<AxisSiteDefinition> sites(static_cast<size_t>(n_sites));
    for (AxisSiteDefinition& site : sites){
        uint8_t found=0,mitochondrial=0;
        joint_read_binary_value(input,site.tid,path);
        joint_read_binary_value(input,site.pos,path);
        joint_read_binary_value(input,found,path);
        joint_read_binary_value(input,mitochondrial,path);
        site.key=site_key(site.tid,site.pos);
        site.found=found!=0; site.mitochondrial=mitochondrial!=0;
        site.contig=joint_read_binary_string(input,path);
        site.ref_allele=uppercase(joint_read_binary_string(input,path));
        site.alt_allele=uppercase(joint_read_binary_string(input,path));
        const uint64_t strings=site.contig.capacity()+site.ref_allele.capacity()+site.alt_allele.capacity();
        if(strings>memory_limit-allocated)
            throw runtime_error("normalized allele metadata exceeds task memory budget");
        allocated+=strings;
        if (site.found && (site.contig.empty() || site.ref_allele.empty() ||
                site.alt_allele.empty()))
            throw runtime_error("normalized site cache lacks canonical allele metadata");
        site.genotype.resize(donors.size());
        for (int8_t& genotype : site.genotype)
            joint_read_binary_value(input,genotype,path);
    }
    char extra=0;
    if (input.read(&extra,1))
        throw runtime_error("normalized site cache has trailing bytes: "+path);
    if (!input.eof()) throw runtime_error("failed reading normalized site cache: "+path);
    sort(sites.begin(),sites.end(),[](const AxisSiteDefinition& left,
                                     const AxisSiteDefinition& right){
        return left.key<right.key;
    });
    return sites;
}

static map<unsigned long,JointNormalizedCellIndex> joint_load_normalized_index(
        const string& prefix,const JointAnalysisTask& task){
    const vector<map<string,string>> rows=joint_read_tsv_maps(prefix+".cells.tsv");
    map<unsigned long,JointNormalizedCellIndex> result;
    for (const auto& row : rows){
        if (joint_map_value(row,"schema_version")!=JOINT_NORMALIZED_CACHE_SCHEMA ||
                joint_map_value(row,"workload_generation_id")!=task.cache_generation ||
                joint_map_value(row,"library")!=task.library ||
                joint_map_value(row,"modality")!=task.modality)
            throw runtime_error("normalized cell index provenance mismatch");
        JointNormalizedCellIndex item;
        item.observation_offset=stoull(joint_map_value(row,"observation_offset"));
        item.observation_count=stoull(joint_map_value(row,"observation_count"));
        item.molecule_offset=stoull(joint_map_value(row,"molecule_offset"));
        item.molecule_count=stoull(joint_map_value(row,"molecule_count"));
        item.malformed_molecule_rows=stol(joint_map_value(row,"malformed_molecule_rows"));
        const unsigned long barcode=strtoul(
            joint_map_value(row,"encoded_barcode").c_str(),NULL,10);
        if (!result.emplace(barcode,item).second)
            throw runtime_error("duplicate barcode in normalized cell index");
    }
    return result;
}

template <typename Record>
static vector<Record> joint_read_normalized_records(
        const string& path,uint64_t offset,uint64_t count){
    if (count>100000000ULL)
        throw runtime_error("refusing an implausibly large single-cell cache slice");
    ifstream input(path.c_str(),ios::binary);
    if (!input) throw runtime_error("could not open normalized records: "+path);
    input.seekg(static_cast<streamoff>(offset*sizeof(Record)),ios::beg);
    if (!input) throw runtime_error("could not seek normalized records: "+path);
    vector<Record> records(static_cast<size_t>(count));
    if (count) input.read(reinterpret_cast<char*>(records.data()),
                          static_cast<streamsize>(count*sizeof(Record)));
    if (!input) throw runtime_error("truncated normalized record slice: "+path);
    return records;
}

static vector<JointLinkedUnit> joint_site_as_linked(
        const vector<JointSiteUnit>& sites){
    vector<JointLinkedUnit> result;
    result.reserve(sites.size());
    for (const JointSiteUnit& site : sites){
        JointLinkedUnit unit;
        unit.molecule=site_key(site.tid,site.pos); unit.basis=0;
        unit.sites.push_back(site); result.push_back(unit);
    }
    return result;
}

static string joint_unit_key(const JointLinkedUnit& unit){
    if (unit.basis==0 && unit.sites.size()==1){
        const JointSiteUnit& site=unit.sites.front();
        if (site.contig.empty() || site.ref_allele.empty() ||
                site.alt_allele.empty())
            throw runtime_error("SITE unit lacks canonical allele metadata");
        return "SITE:"+site.contig+"|"+to_string(site.pos)+"|"+
            uppercase(site.ref_allele)+"|"+uppercase(site.alt_allele);
    }
    return to_string((unsigned int)unit.basis)+":"+
        (unit.molecule_id.empty() ? to_string(unit.molecule) : unit.molecule_id);
}

static pair<long double,long double> joint_compiled_counts(
        const JointSiteUnit& site,
        const vector<pair<double,double>>* replacement_counts){
    if (replacement_counts && site.compiled_site_id<replacement_counts->size()){
        const pair<double,double>& value=(*replacement_counts)[site.compiled_site_id];
        if (isfinite(value.first) && isfinite(value.second))
            return make_pair((long double)value.first,(long double)value.second);
    }
    return make_pair(site.ref,site.alt);
}

static long double joint_compiled_site_log_likelihood(
        const JointSiteUnit& site,long double alpha,
        const vector<pair<double,double>>* replacement_counts){
    if (joint_count_targeted_work){
        ++joint_targeted_likelihood_evaluations;
        ++joint_targeted_row_scans;
    }
    const pair<long double,long double> counts=
        joint_compiled_counts(site,replacement_counts);
    const long double epsilon=1e-15L;
    const long double probability=min(max(
        site.probability_intercept+alpha*site.probability_slope,
        epsilon),1.0L-epsilon);
    const long double depth=counts.first+counts.second;
    if (depth<=0.0L) return 0.0L;
    return (counts.second*logl(probability)+
            counts.first*logl(1.0L-probability))/depth;
}

static pair<long double,long double> joint_compiled_site_derivatives(
        const JointSiteUnit& site,long double alpha,
        const vector<pair<double,double>>* replacement_counts){
    if (joint_count_targeted_work){
        ++joint_targeted_derivative_evaluations;
        ++joint_targeted_row_scans;
    }
    const pair<long double,long double> counts=
        joint_compiled_counts(site,replacement_counts);
    const long double depth=counts.first+counts.second;
    if (depth<=0.0L) return make_pair(0.0L,0.0L);
    const long double epsilon=1e-15L;
    const long double raw=site.probability_intercept+
        alpha*site.probability_slope;
    const long double probability=min(max(raw,epsilon),1.0L-epsilon);
    const long double slope=(raw<=epsilon || raw>=1.0L-epsilon) ?
        0.0L : site.probability_slope;
    const long double first=slope*(counts.second/probability-
        counts.first/(1.0L-probability))/depth;
    const long double second=-slope*slope*(
        counts.second/(probability*probability)+
        counts.first/((1.0L-probability)*(1.0L-probability)))/depth;
    return make_pair(first,second);
}

static long double joint_compiled_unit_score(
        const JointLinkedUnit& unit,long double alpha,
        const vector<pair<double,double>>* replacement_counts){
    AxisKahan total;
    for (const JointSiteUnit& site : unit.sites)
        total.add(joint_compiled_site_log_likelihood(
            site,alpha,replacement_counts));
    return unit.sites.empty() ? NAN :
        total.value/(long double)unit.sites.size();
}

static string joint_canonical_site_key(const JointSiteUnit& site){
    if (site.contig.empty() || site.ref_allele.empty() || site.alt_allele.empty())
        throw runtime_error("site is missing canonical contig/REF/ALT metadata");
    return site.contig+"|"+to_string(site.pos)+"|"+
        uppercase(site.ref_allele)+"|"+uppercase(site.alt_allele);
}

static int joint_compiled_fold(
        const JointLinkedUnit& unit,const JointHypothesis& hypothesis,
        const string& channel){
    string token;
    if (channel=="SITE"){
        if (unit.sites.size()!=1)
            throw runtime_error("SITE fold unit must contain exactly one site");
        token=joint_canonical_site_key(unit.sites.front());
    } else {
        token=hypothesis.library+"|"+hypothesis.barcode+"|"+
            unit.parent_origin+"|"+joint_molecule_fold_basis_name(unit.basis)+"|"+
            (unit.molecule_id.empty() ? to_string(unit.molecule) :
             unit.molecule_id);
    }
    token+="|"+string(JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION)+"|FOLD_V1";
    return static_cast<int>(stable_text_hash(token)%5ULL);
}

template<class ScoreFunction>
static pair<long double,long double> joint_bounded_profile_interval(
        long double cap,long double maximum,long double cutoff,const ScoreFunction& score){
    const auto passes=[&](int i){const long double value=score(cap*(long double)i/200.0L);return isfinite(value)&&value>=cutoff;};
    int peak=max(0,min(200,(int)floorl(maximum/cap*200.0L)));
    if(!passes(peak)){
        if(peak<200 && passes(peak+1))++peak;
        else return make_pair(NAN,NAN);
    }
    int lo=0,hi=peak;
    while(lo<hi){int mid=(lo+hi)/2;if(passes(mid))hi=mid;else lo=mid+1;}
    int first=lo;lo=peak;hi=200;
    while(lo<hi){int mid=(lo+hi+1)/2;if(passes(mid))lo=mid;else hi=mid-1;}
    return make_pair(cap*(long double)first/200.0L,cap*(long double)lo/200.0L);
}

static JointCompiledFit joint_compiled_fit(
        const JointStoredUnits& units,const JointHypothesis& hypothesis,
        const string& channel,long double e_ref,long double e_alt,
        long min_evidence,long double max_alpha,
        const vector<uint8_t>* retained_mask=NULL,
        const vector<pair<double,double>>* replacement_counts=NULL,
        bool evaluate_folds=true){
    (void)e_ref; (void)e_alt; // already represented by compiled coefficients
    if (joint_count_targeted_work) ++joint_targeted_optimized_fit_calls;
    JointCompiledFit result;
    if (!hypothesis.locked.valid || !hypothesis.second.valid){
        result.status=channel=="SITE" ?
            "DONOR_NOT_IN_MODALITY_PANEL" : "UNAVAILABLE";
        result.unavailable_reason="DONOR_NOT_IN_MODALITY_PANEL";
        return result;
    }
    size_t informative_count=0;
    size_t selected_units=0;
    for (const JointLinkedUnit& unit : units){
        if (retained_mask && (unit.compiled_unit_id>=retained_mask->size() ||
                !(*retained_mask)[unit.compiled_unit_id])) continue;
        ++selected_units;
        if (unit.genotype_distinguishable){
            ++informative_count;
            if (channel=="SITE")
                for (const JointSiteUnit& site : unit.sites){
                    const auto counts=joint_compiled_counts(
                        site,replacement_counts);
                    result.discriminating_depth+=counts.first+counts.second;
                }
        }
    }
    result.usable_units=selected_units;
    result.units=informative_count;
    if (informative_count==0){
        result.status=selected_units==0 ?
            (channel=="SITE" ? "NO_OBSERVATIONS" : "UNAVAILABLE") :
            "GENOTYPE_EQUIVALENT";
        result.unavailable_reason=selected_units==0 ?
            (channel=="SITE" ? "NO_OBSERVATIONS" : "NO_USABLE_LINKED_UNITS") :
            (channel=="SITE" ? "NONE" :
             "NO_GENOTYPE_DISTINGUISHABLE_LINKED_UNITS");
        if (!selected_units) return result;
        result.alpha=0.0L;
        if (selected_units){
            if (channel=="SITE"){
                AxisKahan baseline;
                for (const JointLinkedUnit& unit : units){
                    if (retained_mask && (unit.compiled_unit_id>=retained_mask->size() ||
                            !(*retained_mask)[unit.compiled_unit_id])) continue;
                    baseline.add(joint_compiled_unit_score(
                        unit,0.0L,replacement_counts));
                }
                result.locked=baseline.value/(long double)selected_units;
                result.interior=result.locked;
                result.contributor=result.locked;
                result.delta=0.0L;
            }
            result.preferred="LOCKED_IDENTITY";
        }
        return result;
    }
    const auto score=[&](long double alpha){
        AxisKahan total;
        for (const JointLinkedUnit& unit : units){
            if (retained_mask && (unit.compiled_unit_id>=retained_mask->size() ||
                    !(*retained_mask)[unit.compiled_unit_id])) continue;
            if (unit.genotype_distinguishable)
                total.add(joint_compiled_unit_score(
                    unit,alpha,replacement_counts));
        }
        return total.value/(long double)informative_count;
    };
    const auto derivatives=[&](long double alpha){
        AxisKahan first,second;
        for (const JointLinkedUnit& unit : units){
            if (retained_mask && (unit.compiled_unit_id>=retained_mask->size() ||
                    !(*retained_mask)[unit.compiled_unit_id])) continue;
            if (!unit.genotype_distinguishable) continue;
            AxisKahan unit_first,unit_second;
            for (const JointSiteUnit& site : unit.sites){
                const auto value=joint_compiled_site_derivatives(
                    site,alpha,replacement_counts);
                unit_first.add(value.first); unit_second.add(value.second);
            }
            const long double count=(long double)unit.sites.size();
            first.add(unit_first.value/count); second.add(unit_second.value/count);
        }
        const long double count=(long double)informative_count;
        return make_pair(first.value/count,second.value/count);
    };
    const pair<long double,long double> fitted=joint_concave_maximum(
        max_alpha,derivatives,score);
    AxisKahan locked,contributor;
    for (const JointLinkedUnit& unit : units){
        if (retained_mask && (unit.compiled_unit_id>=retained_mask->size() ||
                !(*retained_mask)[unit.compiled_unit_id])) continue;
        if (!unit.genotype_distinguishable) continue;
        locked.add(joint_compiled_unit_score(
            unit,0.0L,replacement_counts));
        contributor.add(joint_compiled_unit_score(
            unit,1.0L,replacement_counts));
    }
    result.status=channel=="SITE" ?
        (informative_count<2 || result.discriminating_depth<min_evidence ?
         "LOW_EVIDENCE" : "AVAILABLE") :
        (informative_count<2 ? "LIMITED_EVIDENCE" : "AVAILABLE");
    result.alpha=fitted.first;
    result.locked=locked.value/(long double)informative_count;
    result.interior=fitted.second;
    result.contributor=contributor.value/(long double)informative_count;
    result.delta=result.interior-result.locked;
    if (evaluate_folds){
    const long double cutoff=result.interior-
        1.920729410347062L/(long double)informative_count;
    const pair<long double,long double> profile=joint_bounded_profile_interval(
        max_alpha,result.alpha,cutoff,score);
    result.alpha_low=profile.first;
    result.alpha_high=profile.second;
    }
    result.preferred="LOCKED_IDENTITY";
    long double preferred=result.locked;
    if (result.alpha>1e-12L && result.interior>preferred){
        preferred=result.interior;
        result.preferred="INTERIOR_LOCKED_PLUS_CONTRIBUTOR";
    }
    if (result.contributor>preferred)
        result.preferred="CONTRIBUTOR_ONLY_FRACTION_1";
    if(!evaluate_folds)return result;
    AxisKahan absolute;
    long double largest=0.0L;
    for (const JointLinkedUnit& unit : units){
        if (retained_mask && (unit.compiled_unit_id>=retained_mask->size() ||
                !(*retained_mask)[unit.compiled_unit_id])) continue;
        if (!unit.genotype_distinguishable) continue;
        const long double influence=fabsl(joint_compiled_unit_score(
            unit,result.alpha,replacement_counts)-
            joint_compiled_unit_score(unit,0.0L,replacement_counts));
        absolute.add(influence); largest=max(largest,influence);
    }
    if (absolute.value>0.0L)
        result.maximum_influence_fraction=largest/absolute.value;
    if (evaluate_folds){
        vector<vector<size_t>> heldout(5);
        for (size_t index=0;index<units.size();++index){
            const JointLinkedUnit& unit=units[index];
            if (retained_mask && (unit.compiled_unit_id>=retained_mask->size() ||
                    !(*retained_mask)[unit.compiled_unit_id])) continue;
            if (!unit.genotype_distinguishable) continue;
            heldout[joint_compiled_fold(unit,hypothesis,channel)].push_back(index);
        }
        for (int fold=0;fold<5;++fold){
            if (heldout[fold].empty() || heldout[fold].size()==informative_count)
                continue;
            vector<uint8_t> training_mask;
            const vector<uint8_t>* mask_pointer=NULL;
            if (retained_mask) training_mask=*retained_mask;
            else {
                size_t mask_size=0;
                for (const JointLinkedUnit& unit : units)
                    mask_size=max(mask_size,(size_t)unit.compiled_unit_id+1);
                training_mask.assign(mask_size,1);
            }
            for (size_t index : heldout[fold])
                training_mask[units[index].compiled_unit_id]=0;
            mask_pointer=&training_mask;
            JointCompiledFit training=joint_compiled_fit(
                units,hypothesis,channel,e_ref,e_alt,min_evidence,max_alpha,
                mask_pointer,replacement_counts,false);
            if (!isfinite(training.alpha) || !isfinite(training.delta)) continue;
            ++result.folds_evaluable;
            bool supports=false;
            if (channel=="SITE"){
                supports=training.alpha>0.01L && training.delta>0.0L;
            } else {
                AxisKahan heldout_improvement;
                for (size_t index : heldout[fold])
                    heldout_improvement.add(joint_compiled_unit_score(
                        units[index],training.alpha,replacement_counts)-
                        joint_compiled_unit_score(units[index],0.0L,
                                                  replacement_counts));
                supports=heldout_improvement.value/
                    (long double)heldout[fold].size()>0.0L;
            }
            if (supports) ++result.fold_support_numerator;
        }
        if (result.folds_evaluable>0)
            result.fold_support_fraction=(long double)result.fold_support_numerator/
                (long double)result.folds_evaluable;
    }
    return result;
}

static JointCompiledFit joint_reference_fit(
        const JointStoredUnits& source,const JointHypothesis& hypothesis,
        const string& channel,long double e_ref,long double e_alt,
        long min_evidence,long double max_alpha,
        const vector<uint8_t>* retained_mask=NULL,
        const vector<pair<double,double>>* replacement_counts=NULL,
        bool evaluate_folds=true){
    (void)evaluate_folds;
    if (joint_count_targeted_work) ++joint_targeted_reference_fit_calls;
    // Independent scalar path: materialize only the bounded candidate slice,
    // then call the established likelihood and optimizer functions directly.
    // It intentionally never calls joint_compiled_fit().
    vector<JointLinkedUnit> units;
    for (const JointLinkedUnit& original : source){
        if (retained_mask && (original.compiled_unit_id>=retained_mask->size() ||
                !(*retained_mask)[original.compiled_unit_id])) continue;
        JointLinkedUnit unit=original;
        if (replacement_counts)
            for (JointSiteUnit& site : unit.sites){
                const auto counts=joint_compiled_counts(site,replacement_counts);
                site.ref=counts.first; site.alt=counts.second;
            }
        units.push_back(move(unit));
    }
    JointCompiledFit result;
    if (!hypothesis.locked.valid || !hypothesis.second.valid){
        result.status=channel=="SITE" ?
            "DONOR_NOT_IN_MODALITY_PANEL" : "UNAVAILABLE";
        result.unavailable_reason="DONOR_NOT_IN_MODALITY_PANEL";
        return result;
    }
    vector<size_t> informative;
    for (size_t i=0;i<units.size();++i){
        bool distinguishes=false;
        for (const JointSiteUnit& site : units[i].sites)
            if (fabsl(site.q_locked-site.q_second)>1e-18L){
                distinguishes=true; break;
            }
        if (distinguishes){
            informative.push_back(i);
            if (channel=="SITE")
                for (const JointSiteUnit& site : units[i].sites)
                    result.discriminating_depth+=site.ref+site.alt;
        }
    }
    result.usable_units=units.size();
    result.units=informative.size();
    if (informative.empty()){
        result.status=units.empty() ?
            (channel=="SITE" ? "NO_OBSERVATIONS" : "UNAVAILABLE") :
            "GENOTYPE_EQUIVALENT";
        result.unavailable_reason=units.empty() ?
            (channel=="SITE" ? "NO_OBSERVATIONS" : "NO_USABLE_LINKED_UNITS") :
            (channel=="SITE" ? "NONE" :
             "NO_GENOTYPE_DISTINGUISHABLE_LINKED_UNITS");
        if (units.empty()) return result;
        result.alpha=0.0L;
        if (channel=="SITE"){
            AxisKahan baseline;
            for (const JointLinkedUnit& unit : units)
                baseline.add(joint_linked_unit_log_likelihood(
                    unit,0.0L,hypothesis.rho_effective,e_ref,e_alt));
            result.locked=baseline.value/(long double)units.size();
            result.interior=result.contributor=result.locked;
            result.delta=0.0L;
        }
        result.preferred="LOCKED_IDENTITY";
        return result;
    }
    const JointFit fitted=joint_fit_linked_units(
        units,informative,hypothesis.rho_effective,e_ref,e_alt,max_alpha);
    AxisKahan locked,contributor;
    for (size_t index : informative){
        locked.add(joint_linked_unit_log_likelihood(
            units[index],0.0L,hypothesis.rho_effective,e_ref,e_alt));
        contributor.add(joint_linked_unit_log_likelihood(
            units[index],1.0L,hypothesis.rho_effective,e_ref,e_alt));
    }
    result.alpha=fitted.alpha; result.interior=fitted.balanced;
    result.locked=locked.value/(long double)informative.size();
    result.contributor=contributor.value/(long double)informative.size();
    result.delta=result.interior-result.locked;
    result.status=channel=="SITE" ?
        (informative.size()<2 || result.discriminating_depth<min_evidence ?
         "LOW_EVIDENCE" : "AVAILABLE") :
        (informative.size()<2 ? "LIMITED_EVIDENCE" : "AVAILABLE");
    const long double cutoff=result.interior-
        1.920729410347062L/(long double)informative.size();
    const auto score=[&](long double alpha){
        AxisKahan total;
        for (size_t index : informative)
            total.add(joint_linked_unit_log_likelihood(
                units[index],alpha,hypothesis.rho_effective,e_ref,e_alt));
        return total.value/(long double)informative.size();
    };
    const auto profile=joint_profile_interval(
        max_alpha,result.alpha,cutoff,score);
    result.alpha_low=profile.first; result.alpha_high=profile.second;
    result.preferred="LOCKED_IDENTITY";
    long double preferred=result.locked;
    if (result.alpha>1e-12L && result.interior>preferred){
        preferred=result.interior;
        result.preferred="INTERIOR_LOCKED_PLUS_CONTRIBUTOR";
    }
    if (result.contributor>preferred)
        result.preferred="CONTRIBUTOR_ONLY_FRACTION_1";
    AxisKahan absolute_influence;
    long double largest_influence=0.0L;
    for (size_t index : informative){
        const long double influence=fabsl(joint_linked_unit_log_likelihood(
            units[index],result.alpha,hypothesis.rho_effective,e_ref,e_alt)-
            joint_linked_unit_log_likelihood(
                units[index],0.0L,hypothesis.rho_effective,e_ref,e_alt));
        absolute_influence.add(influence);
        largest_influence=max(largest_influence,influence);
    }
    if (absolute_influence.value>0.0L)
        result.maximum_influence_fraction=
            largest_influence/absolute_influence.value;
    if (evaluate_folds){
        vector<vector<size_t>> heldout(5);
        for (size_t index : informative)
            heldout[joint_compiled_fold(
                units[index],hypothesis,channel)].push_back(index);
        for (int fold=0;fold<5;++fold){
            if (heldout[fold].empty() || heldout[fold].size()==informative.size())
                continue;
            vector<JointLinkedUnit> training;
            for (size_t index : informative)
                if (find(heldout[fold].begin(),heldout[fold].end(),index)==
                        heldout[fold].end())
                    training.push_back(units[index]);
            JointCompiledFit training_fit=joint_reference_fit(
                training,hypothesis,channel,e_ref,e_alt,min_evidence,max_alpha,
                NULL,NULL,false);
            if (!isfinite(training_fit.alpha) || !isfinite(training_fit.delta))
                continue;
            ++result.folds_evaluable;
            bool supports=training_fit.alpha>0.01L && training_fit.delta>0.0L;
            if (channel=="MOLECULE"){
                AxisKahan heldout_delta;
                for (size_t index : heldout[fold])
                    heldout_delta.add(joint_linked_unit_log_likelihood(
                        units[index],training_fit.alpha,
                        hypothesis.rho_effective,e_ref,e_alt)-
                        joint_linked_unit_log_likelihood(units[index],0.0L,
                            hypothesis.rho_effective,e_ref,e_alt));
                supports=heldout_delta.value/(long double)heldout[fold].size()>0.0L;
            }
            if (supports) ++result.fold_support_numerator;
        }
        if (result.folds_evaluable>0)
            result.fold_support_fraction=(long double)result.fold_support_numerator/
                (long double)result.folds_evaluable;
    }
    return result;
}

static const JointStoredUnits& joint_channel_units_ref(
        const JointCompiledCandidate& candidate,const string& channel);
static string joint_fit_category(
        const JointCompiledFit& fit,long double winner_margin,
        bool assay_or_evidence_conflict);

static JointCompiledCandidate joint_build_cached_candidate(
        const JointHypothesis& hypothesis,
        const JointRecordView<AxisObservationRecord>& observations,
        const JointRecordView<JointMoleculeRecord>& molecules,
        const vector<AxisSiteDefinition>& sites,const unordered_map<int,size_t>& donor_slot){
    JointCompiledCandidate candidate;candidate.hypothesis=hypothesis;
    const auto make_site=[&](int32_t tid,int32_t pos,double ref,double alt,
                              bool molecule,JointSiteUnit& site){
        if(!hypothesis.locked.valid || !hypothesis.second.valid)return false;
        const uint64_t key=site_key(tid,pos);
        auto found=lower_bound(sites.begin(),sites.end(),key,
            [](const AxisSiteDefinition& value,uint64_t k){return value.key<k;});
        const long double depth=ref+alt;
        if(found==sites.end() || found->key!=key || !found->found || found->mitochondrial || depth<=0)return false;
        site.tid=tid;site.pos=pos;site.contig=found->contig;
        site.ref_allele=uppercase(found->ref_allele);site.alt_allele=uppercase(found->alt_allele);
        site.ref=molecule?ref/depth:ref;site.alt=molecule?alt/depth:alt;
        if(!joint_expected(hypothesis.locked,*found,donor_slot,site.q_locked) ||
           !joint_expected(hypothesis.second,*found,donor_slot,site.q_second) ||
           (hypothesis.rho_effective>0 && !joint_expected(hypothesis.ambient,*found,donor_slot,site.q_ambient,true)))return false;
        if(hypothesis.rho_effective==0)site.q_ambient=site.q_locked;
        return true;
    };
    for(const AxisObservationRecord& record:observations){
        JointSiteUnit site;
        if(!make_site(record.tid,record.pos,record.ref,record.alt,false,site))continue;
        JointLinkedUnit unit;unit.molecule=site_key(site.tid,site.pos);unit.basis=0;
        unit.sites.push_back(move(site));candidate.site_units.push_back(unit);
    }
    JointLinkedUnit unit;bool started=false;
    for(const JointMoleculeRecord& record:molecules){
        if(started && (unit.molecule!=record.molecule || unit.basis!=record.basis)){
            if(!unit.sites.empty())candidate.molecule_units.push_back(unit);
            unit=JointLinkedUnit();started=false;
        }
        if(!started){unit.molecule=record.molecule;unit.molecule_id=to_string(record.molecule);unit.basis=record.basis;started=true;}
        JointSiteUnit site;
        if(make_site(record.tid,record.pos,record.ref,record.alt,true,site))unit.sites.push_back(move(site));
    }
    if(!unit.sites.empty())candidate.molecule_units.push_back(unit);
    candidate.site_reference.status=!hypothesis.locked.valid || !hypothesis.second.valid?
        "DONOR_NOT_IN_MODALITY_PANEL":candidate.site_units.empty()?
        (observations.empty()?"NO_OBSERVATIONS":"NO_COMMON_NUCLEAR_GENOTYPES"):"AVAILABLE";
    candidate.molecule_reference.status=candidate.molecule_units.empty()?"UNAVAILABLE":"AVAILABLE";
    return candidate;
}

static vector<JointCompiledCandidate> joint_evaluate_cached_cell(
        unsigned long encoded,const vector<JointHypothesis>& hypotheses,
        const unordered_map<unsigned long,vector<size_t>>& by_cell,
        const JointRecordView<AxisObservationRecord>& observations,
        const JointRecordView<JointMoleculeRecord>& molecules,
        const vector<AxisSiteDefinition>& sites,
        const unordered_map<int,size_t>& donor_slot,long malformed,
        long double e_ref,long double e_alt,long min_evidence,
        long double max_alpha,bool authoritative_reference=false){
    (void)malformed;(void)e_ref;(void)e_alt;(void)min_evidence;(void)max_alpha;(void)authoritative_reference;
    auto found=by_cell.find(encoded);
    if(found==by_cell.end())throw runtime_error("requested cell absent from complete menu");
    vector<JointCompiledCandidate> result;
    for(size_t index:found->second)result.push_back(joint_build_cached_candidate(
        hypotheses[index],observations,molecules,sites,donor_slot));
    sort(result.begin(),result.end(),[](const JointCompiledCandidate& a,const JointCompiledCandidate& b){
        return a.hypothesis.candidate_id<b.hypothesis.candidate_id;});
    return result;
}

static JointCompiledUniverse joint_compile_channel_universe(
        vector<JointCompiledCandidate>& candidates,const string& channel,
        long double e_ref,long double e_alt){
    set<string> unit_keys,site_keys;
    for (const JointCompiledCandidate& candidate : candidates)
        for (const JointLinkedUnit& unit :
                joint_channel_units_ref(candidate,channel)){
            const string unit_key=joint_unit_key(unit);
            unit_keys.insert(unit_key);
            for (const JointSiteUnit& site : unit.sites)
                site_keys.insert(unit_key+":"+joint_canonical_site_key(site));
        }
    JointCompiledUniverse universe;
    universe.unit_keys.assign(unit_keys.begin(),unit_keys.end());
    universe.unit_hashes.reserve(universe.unit_keys.size());
    for (const string& key : universe.unit_keys)
        universe.unit_hashes.push_back(stable_text_hash(key));
    universe.site_keys.assign(site_keys.begin(),site_keys.end());
    map<string,uint32_t> unit_id,site_id;
    for (size_t i=0;i<universe.unit_keys.size();++i)
        unit_id[universe.unit_keys[i]]=static_cast<uint32_t>(i);
    for (size_t i=0;i<universe.site_keys.size();++i)
        site_id[universe.site_keys[i]]=static_cast<uint32_t>(i);
    map<uint32_t,JointLinkedUnit> templates;
    for (JointCompiledCandidate& candidate : candidates){
        const JointStoredUnits& units=channel=="SITE" ?
            candidate.site_units : candidate.molecule_units;
        JointStoredUnits compiled;
        for (JointLinkedUnit unit : units){
            const string unit_key=joint_unit_key(unit);
            unit.compiled_unit_id=unit_id.at(unit_key);
            unit.genotype_distinguishable=0;
            for (JointSiteUnit& site : unit.sites){
                site.compiled_site_id=site_id.at(unit_key+":"+
                    joint_canonical_site_key(site));
                if (fabsl(site.q_locked-site.q_second)>1e-18L)
                    unit.genotype_distinguishable=1;
                const long double error_scale=1.0L-e_ref-e_alt;
                const long double rho=candidate.hypothesis.rho_effective;
                const long double base=(1.0L-rho)*site.q_locked+
                    rho*site.q_ambient;
                site.probability_intercept=e_ref+error_scale*base;
                site.probability_slope=error_scale*(1.0L-rho)*
                    (site.q_second-site.q_locked);
            }
            sort(unit.sites.begin(),unit.sites.end(),[](
                    const JointSiteUnit& left,const JointSiteUnit& right){
                return left.compiled_site_id<right.compiled_site_id;
            });
            if (!templates.count(unit.compiled_unit_id)){
                templates[unit.compiled_unit_id]=unit;
            } else {
                // Units can be non-nested across legal candidates.  The null
                // template must be the numeric union of every site in this
                // complete-menu unit, never the first/largest candidate's
                // vector.  Candidate-specific fits still use their own subset.
                JointLinkedUnit& combined=templates[unit.compiled_unit_id];
                set<uint32_t> present;
                for (const JointSiteUnit& site : combined.sites)
                    present.insert(site.compiled_site_id);
                for (const JointSiteUnit& site : unit.sites)
                    if (present.insert(site.compiled_site_id).second)
                        combined.sites.push_back(site);
                combined.genotype_distinguishable=
                    combined.genotype_distinguishable ||
                    unit.genotype_distinguishable;
                sort(combined.sites.begin(),combined.sites.end(),[](
                        const JointSiteUnit& left,const JointSiteUnit& right){
                    return left.compiled_site_id<right.compiled_site_id;
                });
            }
            compiled.push_back(unit);
        }
        if(channel=="SITE")candidate.site_units=move(compiled);
        else candidate.molecule_units=move(compiled);
    }
    for (size_t i=0;i<universe.unit_keys.size();++i){
        auto found=templates.find(static_cast<uint32_t>(i));
        if (found==templates.end())
            throw runtime_error("compiled complete-menu universe has a missing unit template");
        universe.templates.push_back(found->second);
    }
    return universe;
}

static void joint_assert_compiled_reference_equivalence(
        const vector<JointCompiledCandidate>& candidates,const string& channel,
        long double e_ref,long double e_alt,long min_evidence,
        long double max_alpha){
    const long double tolerance=1e-9L;
    for (const JointCompiledCandidate& candidate : candidates){
        JointCompiledFit optimized=joint_compiled_fit(
            joint_channel_units_ref(candidate,channel),candidate.hypothesis,
            channel,e_ref,e_alt,min_evidence,max_alpha);
        JointCompiledFit scalar=joint_reference_fit(
            joint_channel_units_ref(candidate,channel),candidate.hypothesis,
            channel,e_ref,e_alt,min_evidence,max_alpha);
        // Preserve the cache-construction reason for an empty candidate slice.
        // The independent scalar evaluator checks the objective on this same
        // support; the optional Python comparator checks retained scorer rows.
        if (joint_channel_units_ref(candidate,channel).empty()){
            const string status=channel=="SITE" ? candidate.site_reference.status :
                candidate.molecule_reference.status;
            optimized.status=status;
            scalar.status=status;
        }
        if (optimized.status!=scalar.status)
            throw runtime_error("optimized/reference eligibility status mismatch for "+
                channel+":"+candidate.hypothesis.candidate_id+":"+
                optimized.status+":"+scalar.status);
        if (optimized.preferred!=scalar.preferred ||
                optimized.usable_units!=scalar.usable_units ||
                optimized.units!=scalar.units ||
                optimized.folds_evaluable!=scalar.folds_evaluable ||
                optimized.fold_support_numerator!=scalar.fold_support_numerator ||
                joint_fit_category(optimized,NAN,false)!=
                    joint_fit_category(scalar,NAN,false))
            throw runtime_error("optimized/reference model or unit-count mismatch for "+
                channel+":"+candidate.hypothesis.candidate_id);
        const long double pairs[][2]={{optimized.locked,scalar.locked},
            {optimized.interior,scalar.interior},{optimized.contributor,scalar.contributor},
            {optimized.delta,scalar.delta},{optimized.alpha,scalar.alpha},
            {optimized.alpha_low,scalar.alpha_low},{optimized.alpha_high,scalar.alpha_high},
            {optimized.maximum_influence_fraction,scalar.maximum_influence_fraction},
            {optimized.fold_support_fraction,scalar.fold_support_fraction}};
        for (const auto& pair : pairs)
            if (isfinite(pair[0])!=isfinite(pair[1]) ||
                    (isfinite(pair[0]) && fabsl(pair[0]-pair[1])>
                     tolerance*max(1.0L,max(fabsl(pair[0]),fabsl(pair[1])))))
                throw runtime_error("optimized/reference numeric mismatch for "+
                    channel+":"+candidate.hypothesis.candidate_id);
    }
}

static const JointStoredUnits& joint_channel_units_ref(
        const JointCompiledCandidate& candidate,const string& channel){
    return channel=="SITE" ? candidate.site_units : candidate.molecule_units;
}

static vector<pair<size_t,JointCompiledFit>> joint_rank_candidates(
        const vector<JointCompiledCandidate>& candidates,const string& channel,
        long double e_ref,long double e_alt,long min_evidence,
        long double max_alpha,const vector<uint8_t>* retained_mask=NULL,
        const vector<pair<double,double>>* replacement_counts=NULL,
        bool scalar_reference=false,bool evaluate_folds=true){
    vector<pair<size_t,JointCompiledFit>> ranked;
    for (size_t index=0;index<candidates.size();++index){
        JointCompiledFit fit=scalar_reference ? joint_reference_fit(
            joint_channel_units_ref(candidates[index],channel),
            candidates[index].hypothesis,channel,e_ref,e_alt,min_evidence,
            max_alpha,retained_mask,replacement_counts,evaluate_folds) : joint_compiled_fit(
            joint_channel_units_ref(candidates[index],channel),
            candidates[index].hypothesis,channel,e_ref,e_alt,min_evidence,
            max_alpha,retained_mask,replacement_counts,evaluate_folds);
        if (!retained_mask && !replacement_counts &&
                joint_channel_units_ref(candidates[index],channel).empty())
            fit.status=(channel=="SITE" ? candidates[index].site_reference :
                        candidates[index].molecule_reference).status;
        if (fit.status=="AVAILABLE") ranked.push_back(make_pair(index,fit));
    }
    sort(ranked.begin(),ranked.end(),[&](
            const pair<size_t,JointCompiledFit>& left,
            const pair<size_t,JointCompiledFit>& right){
        if (left.second.delta!=right.second.delta)
            return left.second.delta>right.second.delta;
        return candidates[left.first].hypothesis.candidate_id>
            candidates[right.first].hypothesis.candidate_id;
    });
    return ranked;
}

static long double joint_quantile(vector<long double> values,long double p){
    if (values.empty()) return NAN;
    sort(values.begin(),values.end());
    if (p<=0.0L) return values.front();
    if (p>=1.0L) return values.back();
    const long double position=p*(long double)(values.size()-1);
    const size_t lower=(size_t)floorl(position),upper=(size_t)ceill(position);
    if (lower==upper) return values[lower];
    const long double weight=position-(long double)lower;
    return values[lower]*(1.0L-weight)+values[upper]*weight;
}

static string joint_count_json(const map<string,int>& values){
    string result="{";
    bool first=true;
    for (const auto& item : values){
        if (!first) result+=",";
        first=false;
        result+="\""+joint_json_escape(item.first)+"\":"+to_string(item.second);
    }
    return result+"}";
}

static bool joint_fit_interior_eligible(const JointCompiledFit& fit){
    return fit.status=="AVAILABLE" &&
        fit.preferred=="INTERIOR_LOCKED_PLUS_CONTRIBUTOR" &&
        isfinite(fit.alpha) && fit.alpha>0.01L && fit.alpha<=0.50L &&
        isfinite(fit.alpha_low) && fit.alpha_low>0.0L &&
        isfinite(fit.alpha_high) && fit.alpha_high<0.90L &&
        fit.folds_evaluable==5 && isfinite(fit.fold_support_fraction) &&
        fit.fold_support_fraction>=0.80L &&
        isfinite(fit.maximum_influence_fraction) &&
        fit.maximum_influence_fraction<=0.50L;
}

static string joint_fit_category(
        const JointCompiledFit& fit,long double winner_margin=NAN,
        bool assay_or_evidence_conflict=false){
    (void)winner_margin; // raw candidate separation is descriptive only
    if (assay_or_evidence_conflict)
        return "CONFLICTING_ASSAY_OR_EVIDENCE";
    if (fit.status=="LOW_EVIDENCE" || fit.status=="LIMITED_EVIDENCE")
        return "LOW_OR_LIMITED_EVIDENCE";
    if (fit.status!="AVAILABLE" && fit.status!="GENOTYPE_EQUIVALENT")
        return "UNAVAILABLE";
    if (fit.status=="GENOTYPE_EQUIVALENT")
        return "LOCKED_SOURCE_ONLY";
    if (fit.preferred=="CONTRIBUTOR_ONLY_FRACTION_1" ||
            (isfinite(fit.alpha) && fit.alpha>0.5L))
        return "REPLACEMENT_OR_CONTRIBUTOR_ONLY_BOUNDARY";
    if (fit.preferred=="INTERIOR_LOCKED_PLUS_CONTRIBUTOR" &&
            isfinite(fit.alpha) && fit.alpha>0.0L){
        if (joint_fit_interior_eligible(fit)) return "FIT_INTERIOR_ELIGIBLE";
        return "WEAK_ADDITION_COMPATIBLE";
    }
    return "LOCKED_SOURCE_ONLY";
}

static void joint_add_fit_contract_fields(
        JointAnalysisRow& row,const JointCompiledFit& fit){
    const long double cap=stold(row.value.at("max_second_fraction"));
    row.value["fitted_fraction_at_cap"]=bool_text(isfinite(fit.alpha)&&fabsl(fit.alpha-cap)<=1e-10L);
    row.value["usable_units"]=to_string(fit.usable_units);
    row.value["candidate_evaluated_unit_count"]=to_string(fit.units);
    row.value["evidence_units"]=to_string(fit.units);
    row.value["fit_interior_eligible"]=bool_text(
        joint_fit_interior_eligible(fit));
    row.value["profile_low_pass"]=bool_text(
        isfinite(fit.alpha_low) && fit.alpha_low>0.0L);
    row.value["profile_high_pass"]=bool_text(
        isfinite(fit.alpha_high) && fit.alpha_high<0.90L);
    row.value["fraction_range_pass"]=bool_text(
        isfinite(fit.alpha) && fit.alpha>0.01L && fit.alpha<=0.50L);
    row.value["fold_support_pass"]=bool_text(
        fit.folds_evaluable==5 && isfinite(fit.fold_support_fraction) &&
        fit.fold_support_fraction>=0.80L);
    row.value["influence_pass"]=bool_text(
        isfinite(fit.maximum_influence_fraction) &&
        fit.maximum_influence_fraction<=0.50L);
    row.value["folds_requested"]=to_string(fit.folds_requested);
    row.value["folds_evaluable"]=to_string(fit.folds_evaluable);
    row.value["fold_support_numerator"]=to_string(
        fit.fold_support_numerator);
    row.value["fold_support_fraction"]=axis_fmt(
        fit.fold_support_fraction);
    row.value["maximum_influence_fraction"]=axis_fmt(
        fit.maximum_influence_fraction);
}

static JointAnalysisRow joint_analysis_base(
        const JointAnalysisTask& task,const string& result_class,
        const string& equivalence_key,int threads){
    JointAnalysisRow row;
    const string scoped_key=task.library+":"+task.modality+":"+
        to_string(task.task_index)+":"+equivalence_key;
    row.value={
        {"schema_version",JOINT_TARGETED_ANALYSIS_SCHEMA},
        {"scientific_method_version",task.scientific_method_version},
        {"calibration_library",task.calibration_library},
        {"workload_generation_id",task.generation},
        {"task_index",to_string(task.task_index)},
        {"equivalence_key",scoped_key},
        {"result_class",result_class},{"action",task.action},
        {"role",task.role},{"library",task.library},
        {"modality",task.modality},{"threads",to_string(threads)},
        {"error_ref",axis_fmt(task.error_ref)},
        {"error_alt",axis_fmt(task.error_alt)},
        {"min_evidence",to_string(task.min_evidence)},
        {"max_second_fraction",axis_fmt(task.max_second_fraction)},
        {"selection_statistic_definition","maximum over all AVAILABLE legal candidates of fitted locked-plus-contributor minus locked-only mean balanced log likelihood on each candidate common support"},
        {"scheduler_memory_bytes",to_string(task.scheduler_memory_bytes)},
        {"launcher_runtime_reserve_bytes",to_string(
            task.launcher_runtime_reserve_bytes)},
        {"worker_memory_budget_bytes",to_string(task.worker_memory_budget_bytes)}
    };
    return row;
}

static vector<JointAnalysisRow> joint_observed_rows(
        const JointAnalysisTask& task,const vector<JointCompiledCandidate>& candidates,
        const string& barcode,long double e_ref,long double e_alt,int threads,
        bool scalar_reference=false){
    vector<JointAnalysisRow> output;
    for (const string& channel : {string("SITE"),string("MOLECULE")}){
        vector<pair<size_t,JointCompiledFit>> ranked=joint_rank_candidates(
            candidates,channel,e_ref,e_alt,task.min_evidence,
            task.max_second_fraction,NULL,NULL,scalar_reference);
        const long double winner_margin=ranked.empty() ? NAN :
            (ranked.size()>1 ? ranked[0].second.delta-ranked[1].second.delta :
             numeric_limits<long double>::infinity());
        map<size_t,size_t> rank;
        for (size_t i=0;i<ranked.size();++i) rank[ranked[i].first]=i+1;
        for (size_t index=0;index<candidates.size();++index){
            JointCompiledFit fit=scalar_reference ? joint_reference_fit(
                joint_channel_units_ref(candidates[index],channel),
                candidates[index].hypothesis,channel,e_ref,e_alt,
                task.min_evidence,task.max_second_fraction) : joint_compiled_fit(
                joint_channel_units_ref(candidates[index],channel),
                candidates[index].hypothesis,channel,e_ref,e_alt,
                task.min_evidence,task.max_second_fraction);
            if (joint_channel_units_ref(candidates[index],channel).empty())
                fit.status=(channel=="SITE" ? candidates[index].site_reference :
                            candidates[index].molecule_reference).status;
            const string key="observed:"+barcode+":"+channel+":"+
                candidates[index].hypothesis.candidate_id;
            JointAnalysisRow row=joint_analysis_base(
                task,"OBSERVED_FULL_MENU_CANDIDATE",key,threads);
            row.value.insert({
                {"barcode",barcode},{"evidence_channel",channel},
                {"candidate_id",candidates[index].hypothesis.candidate_id},
                {"second_state",candidates[index].hypothesis.second_state},
                {"status",fit.status},{"legal_menu_candidates",to_string(candidates.size())},
                {"rank",rank.count(index) ? to_string(rank[index]) : "NA"},
                {"winner",(!ranked.empty() && ranked.front().first==index) ? "TRUE" : "FALSE"},
                {"locked_log_likelihood",axis_fmt(fit.locked)},
                {"interior_log_likelihood",axis_fmt(fit.interior)},
                {"contributor_only_log_likelihood",axis_fmt(fit.contributor)},
                {"delta_log_likelihood",axis_fmt(fit.delta)},
                {"fitted_fraction",axis_fmt(fit.alpha)},
                {"fitted_fraction_profile_low",axis_fmt(fit.alpha_low)},
                {"fitted_fraction_profile_high",axis_fmt(fit.alpha_high)},
                {"preferred_model",fit.preferred},
                {"evidence_category",joint_fit_category(
                    fit,rank.count(index) && rank[index]==1 ? winner_margin : NAN)},
                {"runner_up",ranked.size()>1 ?
                    candidates[ranked[1].first].hypothesis.candidate_id : "NA"},
                {"winner_margin",axis_fmt(winner_margin)},
                {"raw_delta_margin",axis_fmt(winner_margin)},
                {"calibrated_support_status","PENDING_PRIMARY_EXACT_MASK_AND_COVERAGE"}
            });
            joint_add_fit_contract_fields(row,fit);
            if (scalar_reference) row.value["engine"]="REFERENCE_SCALAR";
            output.push_back(row);
        }
    }
    return output;
}

static size_t joint_python_round_nonnegative(long double value){
    if (value<=0.0L) return 0;
    const long double lower=floorl(value);
    const long double fraction=value-lower;
    if (fraction<0.5L) return static_cast<size_t>(lower);
    if (fraction>0.5L) return static_cast<size_t>(lower+1.0L);
    const unsigned long long integer=static_cast<unsigned long long>(lower);
    return static_cast<size_t>((integer%2ULL)==0ULL ? lower : lower+1.0L);
}

static uint64_t joint_mix_u64(uint64_t left,uint64_t right){
    // SplitMix64 finalizer over two already-canonical hashes.  Unit strings
    // are hashed once per channel, not rebuilt/sorted in every replicate.
    uint64_t value=left^(right+0x9e3779b97f4a7c15ULL+(left<<6)+(left>>2));
    value=(value^(value>>30))*0xbf58476d1ce4e5b9ULL;
    value=(value^(value>>27))*0x94d049bb133111ebULL;
    return value^(value>>31);
}

static vector<JointAnalysisRow> joint_downsample_rows(
        const JointAnalysisTask& task,const vector<JointCompiledCandidate>& candidates,
        const JointCompiledUniverse& site_universe,
        const JointCompiledUniverse& molecule_universe,
        const string& barcode,uint64_t master_seed,int replicates,
        long double e_ref,long double e_alt,int threads){
    vector<JointAnalysisRow> output;
    if (candidates.empty()) return output;
    for (const string& channel : {string("SITE"),string("MOLECULE")}){
        const JointCompiledUniverse& universe=channel=="SITE" ?
            site_universe : molecule_universe;
        const vector<pair<size_t,JointCompiledFit>> full=joint_rank_candidates(
            candidates,channel,e_ref,e_alt,task.min_evidence,
            task.max_second_fraction);
        const long double full_margin=full.empty() ? NAN :
            (full.size()>1 ? full[0].second.delta-full[1].second.delta :
             numeric_limits<long double>::infinity());
        const string full_winner=full.empty() ? "" :
            candidates[full.front().first].hypothesis.candidate_id;
        const string full_category=full.empty() ? "UNAVAILABLE" :
            joint_fit_category(full.front().second,full_margin);
        const string full_model=full.empty() ? "UNAVAILABLE" :
            full.front().second.preferred;
        for (long double fraction : {0.25L,0.50L,0.75L}){
            const size_t retain=universe.unit_keys.empty() ? 0 : max<size_t>(
                1,joint_python_round_nonnegative(
                    (long double)universe.unit_keys.size()*fraction));
            map<string,int> winners,categories,models;
            vector<long double> alphas;
            int full_retained=0,category_retained=0,model_retained=0;
            int successful=0,unavailable=0;
            vector<string> replicate_winners(replicates),
                replicate_categories(replicates),replicate_models(replicates);
            vector<long double> replicate_alphas(replicates,NAN);
            vector<size_t> replicate_candidate_units(replicates,0),
                replicate_usable_units(replicates,0);
            vector<JointAnalysisRow> replicate_rows(replicates);
            vector<string> replicate_errors(replicates);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
            for (int replicate=0;replicate<replicates;++replicate){
                try {
                    vector<pair<uint64_t,uint32_t>> ordered;
                    ordered.reserve(universe.unit_keys.size());
                    const string context=to_string(master_seed)+"|DOWNSAMPLE|"+
                        JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION+"|"+
                        task.library+"|"+task.modality+"|"+channel+"|"+
                        barcode+"|"+axis_fmt(fraction)+"|"+
                        to_string(replicate);
                    const uint64_t context_hash=stable_text_hash(context);
                    for (size_t id=0;id<universe.unit_keys.size();++id){
                        ordered.push_back(make_pair(
                            joint_mix_u64(context_hash,universe.unit_hashes[id]),
                            static_cast<uint32_t>(id)));
                    }
                    if (retain<ordered.size())
                        nth_element(ordered.begin(),ordered.begin()+retain,
                                    ordered.end());
                    vector<uint8_t> selected(universe.unit_keys.size(),0);
                    for (size_t i=0;i<retain && i<ordered.size();++i)
                        selected[ordered[i].second]=1;
                    const vector<pair<size_t,JointCompiledFit>> ranked=
                        joint_rank_candidates(candidates,channel,e_ref,e_alt,
                            task.min_evidence,task.max_second_fraction,&selected);
                    string winner,category,model;
                    long double margin=NAN;
                    if (!ranked.empty()){
                        winner=candidates[ranked.front().first].hypothesis.candidate_id;
                        margin=ranked.size()>1 ?
                            ranked[0].second.delta-ranked[1].second.delta :
                            numeric_limits<long double>::infinity();
                        category=joint_fit_category(ranked.front().second,margin);
                        model=ranked.front().second.preferred;
                        replicate_winners[replicate]=winner;
                        replicate_alphas[replicate]=ranked.front().second.alpha;
                        replicate_categories[replicate]=category;
                        replicate_models[replicate]=model;
                        replicate_candidate_units[replicate]=
                            ranked.front().second.units;
                        replicate_usable_units[replicate]=
                            ranked.front().second.usable_units;
                    }
                    const string replicate_key="downsample_replicate:"+barcode+":"+
                        channel+":"+axis_fmt(fraction)+":"+to_string(replicate);
                    JointAnalysisRow replicate_row=joint_analysis_base(
                        task,"FULL_MENU_DOWNSAMPLING_REPLICATE",replicate_key,threads);
                    replicate_row.value.insert({
                        {"barcode",barcode},{"evidence_channel",channel},
                        {"fraction",axis_fmt(fraction)},
                        {"replicate_index",to_string(replicate)},
                        {"legal_menu_candidates",to_string(candidates.size())},
                        {"status",ranked.empty() ? "UNAVAILABLE_EVIDENCE" : "AVAILABLE"},
                        {"full_menu_winner",winner.empty() ? "NA" : winner},
                        {"runner_up",ranked.size()>1 ?
                            candidates[ranked[1].first].hypothesis.candidate_id : "NA"},
                        {"winner_margin",axis_fmt(margin)},
                        {"delta_log_likelihood",ranked.empty() ? "NA" :
                            axis_fmt(ranked.front().second.delta)},
                        {"locked_log_likelihood",ranked.empty()?"NA":axis_fmt(ranked.front().second.locked)},
                        {"interior_log_likelihood",ranked.empty()?"NA":axis_fmt(ranked.front().second.interior)},
                        {"contributor_only_log_likelihood",ranked.empty()?"NA":axis_fmt(ranked.front().second.contributor)},
                        {"fitted_fraction",ranked.empty() ? "NA" :
                            axis_fmt(ranked.front().second.alpha)},
                        {"preferred_model",model.empty() ? "UNAVAILABLE" : model},
                        {"evidence_category",category.empty() ? "UNAVAILABLE" : category},
                        {"calibrated_support_status",
                         "CALIBRATED_SUPPORT_NOT_EVALUATED"},
                        {"complete_menu_unit_count",to_string(universe.unit_keys.size())},
                        {"derived_seed_key",to_string(master_seed)+
                         "|DOWNSAMPLE|"+JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION+
                         "|"+task.library+"|"+task.modality+"|"+channel+"|"+
                         barcode+"|"+axis_fmt(fraction)+"|"+to_string(replicate)}
                    });
                    if (!ranked.empty())
                        joint_add_fit_contract_fields(
                            replicate_row,ranked.front().second);
                    replicate_rows[replicate]=move(replicate_row);
                } catch (const exception& error){
                    replicate_errors[replicate]=error.what();
                }
            }
            for (int replicate=0;replicate<replicates;++replicate){
                if (!replicate_errors[replicate].empty())
                    throw runtime_error("downsample replicate failed: "+
                                        replicate_errors[replicate]);
                if (!replicate_winners[replicate].empty()){
                    ++successful;
                    ++winners[replicate_winners[replicate]];
                    if (replicate_winners[replicate]==full_winner) ++full_retained;
                    alphas.push_back(replicate_alphas[replicate]);
                    ++categories[replicate_categories[replicate]];
                    if (replicate_categories[replicate]==full_category)
                        ++category_retained;
                    const string& model=replicate_models[replicate];
                    if (model.empty())
                        throw runtime_error("successful downsample replicate has empty model");
                    ++models[model];
                    if (model==full_model) ++model_retained;
                } else ++unavailable;
                output.push_back(move(replicate_rows[replicate]));
            }
            int model_total=0;
            for (const auto& item : models){
                if (item.first.empty())
                    throw runtime_error("empty downsample model-retention key");
                model_total+=item.second;
            }
            if (successful+unavailable!=replicates || model_total!=successful)
                throw runtime_error("downsample replicate/model accounting invariant failed");
            vector<size_t> successful_candidate_units,successful_usable_units;
            for (int replicate=0;replicate<replicates;++replicate)
                if (!replicate_winners[replicate].empty()){
                    successful_candidate_units.push_back(
                        replicate_candidate_units[replicate]);
                    successful_usable_units.push_back(
                        replicate_usable_units[replicate]);
                }
            const string key="downsample:"+barcode+":"+channel+":"+
                axis_fmt(fraction);
            JointAnalysisRow row=joint_analysis_base(
                task,"FULL_MENU_DOWNSAMPLING",key,threads);
            row.value.insert({
                {"barcode",barcode},{"evidence_channel",channel},
                {"fraction",axis_fmt(fraction)},
                {"replicate_count",to_string(replicates)},
                {"requested_replicates",to_string(replicates)},
                {"attempted_replicates",to_string(successful+unavailable)},
                {"successful_replicates",to_string(successful)},
                {"scientifically_unavailable_replicates",to_string(unavailable)},
                {"unavailable_replicates",to_string(unavailable)},
                {"technical_failures","0"},
                {"decision_eligible",bool_text(
                    replicates==100 && successful==100)},
                {"legal_menu_candidates",to_string(candidates.size())},
                {"status",universe.unit_keys.empty() ? "EVIDENCE_UNAVAILABLE_NO_UNITS" :
                    full.empty() || successful==0 ?
                    "EVIDENCE_UNAVAILABLE_NO_SUCCESSFUL_REPLICATE" :
                    successful==replicates ? "AVAILABLE" : "DESCRIPTIVE_PARTIAL"},
                {"full_menu_winner",full_winner.empty() ? "NA" : full_winner},
                {"runner_up",full.size()>1 ?
                    candidates[full[1].first].hypothesis.candidate_id : "NA"},
                {"winner_margin",axis_fmt(full_margin)},
                {"full_menu_category",full_category},
                {"full_menu_preferred_model",full_model},
                {"winner_retention_fraction",successful>0 ? axis_fmt(
                    (long double)full_retained/(long double)successful) : "NA"},
                {"median_winner_fitted_fraction",axis_fmt(joint_quantile(alphas,0.5L))},
                {"fitted_fraction_p025",axis_fmt(joint_quantile(alphas,0.025L))},
                {"fitted_fraction_p975",axis_fmt(joint_quantile(alphas,0.975L))},
                {"category_retention_fraction",successful>0 ? axis_fmt(
                    (long double)category_retained/(long double)successful) : "NA"},
                {"model_retention_fraction",successful>0 ? axis_fmt(
                    (long double)model_retained/(long double)successful) : "NA"},
                {"winner_counts",joint_count_json(winners)},
                {"category_counts",joint_count_json(categories)},
                {"model_counts",joint_count_json(models)},
                {"candidate_evaluated_unit_count_min",
                 successful_candidate_units.empty() ? "NA" : to_string(
                    *min_element(successful_candidate_units.begin(),
                                 successful_candidate_units.end()))},
                {"candidate_evaluated_unit_count_max",
                 successful_candidate_units.empty() ? "NA" : to_string(
                    *max_element(successful_candidate_units.begin(),
                                 successful_candidate_units.end()))},
                {"candidate_usable_unit_count_min",
                 successful_usable_units.empty() ? "NA" : to_string(
                    *min_element(successful_usable_units.begin(),
                                 successful_usable_units.end()))},
                {"candidate_usable_unit_count_max",
                 successful_usable_units.empty() ? "NA" : to_string(
                    *max_element(successful_usable_units.begin(),
                                 successful_usable_units.end()))},
                {"complete_menu_unit_count",to_string(universe.unit_keys.size())},
                {"calibrated_support_status",
                 "CALIBRATED_SUPPORT_NOT_EVALUATED"}
            });
            output.push_back(row);
        }
    }
    return output;
}

static vector<JointAnalysisRow> joint_null_rows(
        const JointAnalysisTask& task,const vector<JointCompiledCandidate>& candidates,
        const JointCompiledUniverse& site_universe,
        const JointCompiledUniverse& molecule_universe,
        const string& barcode,uint64_t master_seed,int replicates,
        long double e_ref,long double e_alt,int threads){
    vector<JointAnalysisRow> output;
    if (candidates.empty()) return output;
    for (const string& channel : {string("SITE"),string("MOLECULE")}){
        const JointCompiledUniverse& universe=channel=="SITE" ?
            site_universe : molecule_universe;
        const JointStoredUnits& baseline=universe.templates;
        const JointHypothesis& null_hypothesis=candidates.front().hypothesis;
        for (const JointCompiledCandidate& candidate : candidates){
            if (candidate.hypothesis.locked_state!=null_hypothesis.locked_state ||
                    fabsl(candidate.hypothesis.rho_effective-
                          null_hypothesis.rho_effective)>1e-18L)
                throw runtime_error(
                    "TECHNICAL_INVALID: locked-null hypothesis invariant mismatch");
        }
        map<uint32_t,tuple<string,long double,long double,long double,long double>>
            invariant;
        for (const JointCompiledCandidate& candidate : candidates)
            for (const JointLinkedUnit& unit :
                    joint_channel_units_ref(candidate,channel))
                for (const JointSiteUnit& site : unit.sites){
                    const auto counts=joint_compiled_counts(site,NULL);
                    const auto value=make_tuple(joint_canonical_site_key(site),
                        site.q_locked,site.q_ambient,counts.first,counts.second);
                    auto found=invariant.find(site.compiled_site_id);
                    if (found==invariant.end()) invariant[site.compiled_site_id]=value;
                    else if (get<0>(found->second)!=get<0>(value) ||
                            fabsl(get<1>(found->second)-get<1>(value))>1e-18L ||
                            fabsl(get<2>(found->second)-get<2>(value))>1e-18L ||
                            fabsl(get<3>(found->second)-get<3>(value))>1e-12L ||
                            fabsl(get<4>(found->second)-get<4>(value))>1e-12L)
                        throw runtime_error(
                            "TECHNICAL_INVALID: locked-null site/count invariant mismatch");
                }
        const vector<pair<size_t,JointCompiledFit>> observed=joint_rank_candidates(
            candidates,channel,e_ref,e_alt,task.min_evidence,
            task.max_second_fraction);
        const long double observed_maximum=observed.empty() ? NAN :
            observed.front().second.delta;
        vector<long double> maxima;
        vector<size_t> successful_candidate_units,successful_usable_units;
        map<string,int> winners;
        size_t layout=0;
        for (const JointLinkedUnit& unit : baseline) layout+=unit.sites.size();
        const int block_size=max(1,min(64,threads));
        for (int block_begin=0;block_begin<replicates;block_begin+=block_size){
            const int block_end=min(replicates,block_begin+block_size);
            const int count=block_end-block_begin;
            vector<vector<pair<double,double>>> replacements(
                count,vector<pair<double,double>>(universe.site_keys.size(),
                                                  make_pair(NAN,NAN)));
            vector<long double> block_maxima(count,NAN);
            vector<size_t> block_candidate_units(count,0),
                block_usable_units(count,0);
            vector<string> block_winners(count),block_errors(count);
            vector<JointAnalysisRow> block_rows(count);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
            for (int local=0;local<count;++local){
                const int replicate=block_begin+local;
                try {
                    const string seed_key=to_string(master_seed)+"|NULL|"+
                        JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION+"|"+
                        task.library+"|"+task.modality+"|"+channel+"|"+
                        barcode+"|"+to_string(replicate);
                    const uint64_t seed=stable_text_hash(seed_key);
                    mt19937_64 random(seed);
                    for (const JointLinkedUnit& unit : baseline)
                        for (const JointSiteUnit& site : unit.sites){
                            const long double probability=joint_adjust_probability(
                                (1.0L-null_hypothesis.rho_effective)*
                                site.q_locked+
                                null_hypothesis.rho_effective*
                                site.q_ambient,e_ref,e_alt);
                            const int depth=static_cast<int>(
                                joint_python_round_nonnegative(site.ref+site.alt));
                            binomial_distribution<int> draw(depth,(double)probability);
                            const int alt=draw(random);
                            replacements[local][site.compiled_site_id]=make_pair(
                                (double)(depth-alt),(double)alt);
                        }
                    const vector<pair<size_t,JointCompiledFit>> ranked=
                        joint_rank_candidates(candidates,channel,e_ref,e_alt,
                            task.min_evidence,task.max_second_fraction,NULL,
                            &replacements[local],false,false);
                    string winner;
                    long double margin=NAN;
                    if (!ranked.empty()){
                        block_maxima[local]=ranked.front().second.delta;
                        winner=candidates[ranked.front().first].hypothesis.candidate_id;
                        block_winners[local]=winner;
                        block_candidate_units[local]=ranked.front().second.units;
                        block_usable_units[local]=
                            ranked.front().second.usable_units;
                        margin=ranked.size()>1 ?
                            ranked[0].second.delta-ranked[1].second.delta :
                            numeric_limits<long double>::infinity();
                    }
                    const string replicate_key="null_replicate:"+barcode+":"+
                        channel+":"+to_string(replicate);
                    JointAnalysisRow replicate_row=joint_analysis_base(
                        task,"CELL_CONDITIONAL_NULL_REPLICATE",replicate_key,threads);
                    replicate_row.value.insert({
                        {"barcode",barcode},{"evidence_channel",channel},
                        {"replicate_index",to_string(replicate)},
                        {"legal_menu_candidates",to_string(candidates.size())},
                        {"status",ranked.empty() ? "UNAVAILABLE_EVIDENCE" : "AVAILABLE"},
                        {"full_menu_winner",winner.empty() ? "NA" : winner},
                        {"runner_up",ranked.size()>1 ?
                            candidates[ranked[1].first].hypothesis.candidate_id : "NA"},
                        {"winner_margin",axis_fmt(margin)},
                        {"delta_log_likelihood",ranked.empty() ? "NA" :
                            axis_fmt(ranked.front().second.delta)},
                        {"locked_log_likelihood",ranked.empty()?"NA":axis_fmt(ranked.front().second.locked)},
                        {"interior_log_likelihood",ranked.empty()?"NA":axis_fmt(ranked.front().second.interior)},
                        {"contributor_only_log_likelihood",ranked.empty()?"NA":axis_fmt(ranked.front().second.contributor)},
                        {"fitted_fraction",ranked.empty() ? "NA" :
                            axis_fmt(ranked.front().second.alpha)},
                        {"preferred_model",ranked.empty() ? "UNAVAILABLE" :
                            ranked.front().second.preferred},
                        {"evidence_category","SUPPORT_CATEGORY_NOT_EVALUATED"},
                        {"calibrated_support_status",
                         "CALIBRATED_SUPPORT_NOT_EVALUATED"},
                        {"candidate_evaluated_unit_count",ranked.empty() ? "0" :
                            to_string(ranked.front().second.units)},
                        {"complete_menu_unit_count",to_string(
                            universe.unit_keys.size())},
                        {"site_layout",to_string(layout)},
                        {"derived_seed_key",seed_key}
                    });
                    block_rows[local]=move(replicate_row);
                } catch (const exception& error){
                    block_errors[local]=error.what();
                }
            }
            for (int local=0;local<count;++local){
                if (!block_errors[local].empty())
                    throw runtime_error("null replicate failed: "+block_errors[local]);
                output.push_back(move(block_rows[local]));
                if (isfinite(block_maxima[local])){
                    maxima.push_back(block_maxima[local]);
                    successful_candidate_units.push_back(
                        block_candidate_units[local]);
                    successful_usable_units.push_back(
                        block_usable_units[local]);
                    ++winners[block_winners[local]];
                }
            }
        }
        int exceed=0;
        if (isfinite(observed_maximum))
            for (long double value : maxima) if (value>=observed_maximum) ++exceed;
        const string key="null:"+barcode+":"+channel;
        JointAnalysisRow row=joint_analysis_base(
            task,"CELL_CONDITIONAL_LOCKED_MODEL_NULL",key,threads);
        row.value.insert({
            {"barcode",barcode},{"evidence_channel",channel},
            {"replicate_count",to_string(replicates)},
            {"requested_replicates",to_string(replicates)},
            {"attempted_replicates",to_string(replicates)},
            {"successful_replicates",to_string(maxima.size())},
            {"scientifically_unavailable_replicates",to_string(
                replicates-(int)maxima.size())},
            {"unavailable_replicates",to_string(
                replicates-(int)maxima.size())},
            {"technical_failures","0"},
            {"legal_menu_candidates",to_string(candidates.size())},
            {"status",baseline.empty() || observed.empty() || maxima.empty() ?
                "UNAVAILABLE_EVIDENCE" :
                maxima.size()==(size_t)replicates ? "AVAILABLE" :
                "DESCRIPTIVE_PARTIAL"},
            {"observed_maximum_delta",axis_fmt(observed_maximum)},
            {"null_median",axis_fmt(joint_quantile(maxima,0.5L))},
            {"null_maximum",maxima.empty() ? "NA" : axis_fmt(
                *max_element(maxima.begin(),maxima.end()))},
            {"null_p95",axis_fmt(joint_quantile(maxima,0.95L))},
            {"null_p99",axis_fmt(joint_quantile(maxima,0.99L))},
            {"empirical_upper_tail_probability",
             isfinite(observed_maximum) && !maxima.empty() ? axis_fmt(
                (long double)(exceed+1)/(long double)(maxima.size()+1)) : "NA"},
            {"empirical_p_numerator",isfinite(observed_maximum) && !maxima.empty() ?
                to_string(exceed+1) : "NA"},
            {"empirical_p_denominator",isfinite(observed_maximum) && !maxima.empty() ?
                to_string(maxima.size()+1) : "NA"},
            {"decision_eligible",bool_text(
                replicates==1000 && maxima.size()==1000)},
            {"winner_counts",joint_count_json(winners)},
            {"candidate_evaluated_unit_count_min",
             successful_candidate_units.empty() ? "NA" : to_string(
                *min_element(successful_candidate_units.begin(),
                             successful_candidate_units.end()))},
            {"candidate_evaluated_unit_count_max",
             successful_candidate_units.empty() ? "NA" : to_string(
                *max_element(successful_candidate_units.begin(),
                             successful_candidate_units.end()))},
            {"candidate_usable_unit_count_min",
             successful_usable_units.empty() ? "NA" : to_string(
                *min_element(successful_usable_units.begin(),
                             successful_usable_units.end()))},
            {"candidate_usable_unit_count_max",
             successful_usable_units.empty() ? "NA" : to_string(
                *max_element(successful_usable_units.begin(),
                             successful_usable_units.end()))},
            {"complete_menu_unit_count",to_string(universe.unit_keys.size())},
            {"site_layout",to_string(layout)},
            {"evidence_category","SUPPORT_CATEGORY_NOT_EVALUATED"},
            {"calibrated_support_status","CALIBRATED_SUPPORT_NOT_EVALUATED"}
        });
        output.push_back(row);
    }
    return output;
}

static vector<long double> joint_parse_fractions(const string& text){
    vector<long double> result;
    for (const string& token : split(text,',')){
        if (trim(token).empty()) continue;
        const long double value=strict_ld(token,"control requested_fractions");
        if (value<=0.0L || value>=1.0L)
            throw runtime_error("control fractions must be strictly between zero and one");
        result.push_back(value);
    }
    if (result.empty()) throw runtime_error("control has no requested fractions");
    return result;
}

static vector<JointCompiledCandidate> joint_score_external_counts_for_menu(
        const vector<JointCompiledCandidate>& menu,
        const JointRecordView<AxisObservationRecord>& observations,
        const JointRecordView<JointMoleculeRecord>& molecules,
        const vector<AxisSiteDefinition>& sites,
        const unordered_map<int,size_t>& donor_slot,long malformed,
        long double e_ref,long double e_alt,long min_evidence,
        long double max_alpha){
    vector<JointCompiledCandidate> result;
    result.reserve(menu.size());
    for (const JointCompiledCandidate& prototype : menu){
        const JointHypothesis& hypothesis=prototype.hypothesis;
        JointCompiledCandidate candidate=joint_build_cached_candidate(
            hypothesis,observations,molecules,sites,donor_slot);
        result.push_back(move(candidate));
    }
    return result;
}

static set<uint8_t> joint_raw_molecule_bases(
        const JointRecordView<JointMoleculeRecord>& rows){
    set<uint8_t> result;
    for (const JointMoleculeRecord& row : rows) result.insert(row.basis);
    return result;
}

static string joint_basis_set_text(const set<uint8_t>& values){
    if (values.empty()) return "UNAVAILABLE";
    string result;
    for (uint8_t value : values){
        if (!result.empty()) result+=",";
        result+=joint_molecule_basis_name(value);
    }
    return result;
}

static vector<uint8_t> joint_deterministic_unit_mask(
        const JointCompiledUniverse& universe,size_t requested,
        const string& seed){
    vector<pair<uint64_t,uint32_t>> order;
    const uint64_t context_hash=stable_text_hash(seed);
    for (size_t i=0;i<universe.unit_keys.size();++i)
        order.push_back(make_pair(joint_mix_u64(
            context_hash,universe.unit_hashes[i]),static_cast<uint32_t>(i)));
    if (requested<order.size())
        nth_element(order.begin(),order.begin()+requested,order.end());
    vector<uint8_t> mask(universe.unit_keys.size(),0);
    for (size_t i=0;i<requested && i<order.size();++i)
        mask[order[i].second]=1;
    return mask;
}

static JointStoredUnits joint_control_candidate_units(
        const JointStoredUnits& recipient,
        const JointStoredUnits& source,
        const vector<uint8_t>& recipient_mask,const vector<uint8_t>& source_mask,
        size_t& recipient_count,size_t& source_count,const string& channel){
    JointStoredUnits mixed;
    recipient_count=source_count=0;
    long double expected_ref=0.0L,expected_alt=0.0L;
    const auto selected=[&](const JointLinkedUnit& unit,
                            const vector<uint8_t>& mask){
        return unit.compiled_unit_id<mask.size() && mask[unit.compiled_unit_id];
    };
    if (channel=="SITE"){
        map<string,JointLinkedUnit> pooled;
        const auto append=[&](const JointLinkedUnit& original,bool is_source){
            if (original.sites.size()!=1)
                throw runtime_error("SITE control input unit is not one canonical site");
            const JointSiteUnit& site=original.sites.front();
            const string key=joint_canonical_site_key(site);
            expected_ref+=site.ref; expected_alt+=site.alt;
            auto found=pooled.find(key);
            if (found==pooled.end()){
                JointLinkedUnit unit=original;
                unit.parent_origin=is_source ? "SOURCE" : "RECIPIENT";
                pooled.emplace(key,move(unit));
            } else {
                JointSiteUnit& combined=found->second.sites.front();
                if (combined.tid!=site.tid ||
                        fabsl(combined.q_locked-site.q_locked)>1e-18L ||
                        fabsl(combined.q_second-site.q_second)>1e-18L ||
                        fabsl(combined.q_ambient-site.q_ambient)>1e-18L)
                    throw runtime_error(
                        "SITE control canonical-key definition mismatch");
                combined.ref+=site.ref;
                combined.alt+=site.alt;
                found->second.genotype_distinguishable=
                    found->second.genotype_distinguishable ||
                    original.genotype_distinguishable;
                found->second.parent_origin="POOLED_RECIPIENT_SOURCE";
            }
        };
        for (const JointLinkedUnit& unit : recipient)
            if (selected(unit,recipient_mask)){ append(unit,false); ++recipient_count; }
        for (const JointLinkedUnit& unit : source)
            if (selected(unit,source_mask)){ append(unit,true); ++source_count; }
        long double actual_ref=0.0L,actual_alt=0.0L;
        for (auto& item : pooled){
            item.second.compiled_unit_id=static_cast<uint32_t>(mixed.size());
            item.second.sites.front().compiled_site_id=
                static_cast<uint32_t>(mixed.size());
            actual_ref+=item.second.sites.front().ref;
            actual_alt+=item.second.sites.front().alt;
            mixed.push_back(move(item.second));
        }
        if (fabsl(actual_ref-expected_ref)>1e-12L ||
                fabsl(actual_alt-expected_alt)>1e-12L ||
                fabsl((actual_ref+actual_alt)-(expected_ref+expected_alt))>1e-12L)
            throw runtime_error("SITE control count/depth conservation failed");
    } else {
        const auto append=[&](const JointLinkedUnit& original,
                              const string& origin){
            JointLinkedUnit unit=original;
            unit.parent_origin=origin;
            unit.compiled_unit_id=static_cast<uint32_t>(mixed.size());
            for (JointSiteUnit& site : unit.sites)
                site.compiled_site_id=static_cast<uint32_t>(mixed.size());
            mixed.push_back(move(unit));
        };
        for (const JointLinkedUnit& unit : recipient)
            if (selected(unit,recipient_mask)){ append(unit,"RECIPIENT"); ++recipient_count; }
        for (const JointLinkedUnit& unit : source)
            if (selected(unit,source_mask)){ append(unit,"SOURCE"); ++source_count; }
    }
    return mixed;
}

static vector<JointAnalysisRow> joint_control_rows(
        const JointAnalysisTask& task,const map<string,string>& control,
        const JointCachedCell& recipient_cell,const JointCachedCell& source_cell,
        const vector<AxisSiteDefinition>& sites,
        const unordered_map<int,size_t>& donor_slot,
        long double e_ref,long double e_alt,int threads,
        bool scalar_reference=false){
    vector<JointAnalysisRow> output;
    const string control_id=joint_map_value(control,"control_id");
    const string control_class=joint_map_value(control,"control_class");
    const string recipient_barcode=joint_map_value(control,"recipient_barcode");
    const string source_identity=joint_map_value(control,"source_identity");
    const string expected=joint_map_value(control,"expected_contributor");
    const string decoy=joint_map_value(control,"genotype_distance_matched_decoy",false);
    const string seed=joint_map_value(control,"seed");
    if (joint_map_value(control,"library")!=task.library)
        throw runtime_error("control/task library mismatch for "+control_id);
    if (recipient_cell.candidates.empty() || source_cell.candidates.empty())
        throw runtime_error("control has empty recipient/source candidate menu");
    const bool source_identity_valid=any_of(
        source_cell.candidates.begin(),source_cell.candidates.end(),
        [&](const JointCompiledCandidate& candidate){
            return candidate.hypothesis.locked_state==source_identity;
        });
    if (!source_identity_valid)
        throw runtime_error("control source identity does not match source-cell locked state: "+
                            control_id);
    vector<JointCompiledCandidate> recipient_candidates=recipient_cell.candidates;
    vector<JointCompiledCandidate> source_candidates=
        joint_score_external_counts_for_menu(
            recipient_candidates,source_cell.raw_observations,
            source_cell.raw_molecules,sites,donor_slot,
            source_cell.malformed_molecule_rows,e_ref,e_alt,
            task.min_evidence,task.max_second_fraction);
    const JointCompiledUniverse recipient_site_universe=
        joint_compile_channel_universe(recipient_candidates,"SITE",e_ref,e_alt);
    const JointCompiledUniverse recipient_molecule_universe=
        joint_compile_channel_universe(recipient_candidates,"MOLECULE",e_ref,e_alt);
    const JointCompiledUniverse source_site_universe=
        joint_compile_channel_universe(source_candidates,"SITE",e_ref,e_alt);
    const JointCompiledUniverse source_molecule_universe=
        joint_compile_channel_universe(source_candidates,"MOLECULE",e_ref,e_alt);
    const set<uint8_t> recipient_bases=joint_raw_molecule_bases(
        recipient_cell.raw_molecules);
    const set<uint8_t> source_bases=joint_raw_molecule_bases(
        source_cell.raw_molecules);
    for (const string& channel : {string("SITE"),string("MOLECULE")}){
        const JointCompiledUniverse& recipient_universe=channel=="SITE" ?
            recipient_site_universe : recipient_molecule_universe;
        const JointCompiledUniverse& source_universe=channel=="SITE" ?
            source_site_universe : source_molecule_universe;
        const bool basis_compatible=channel=="SITE" ||
            (!recipient_bases.empty() && recipient_bases==source_bases);
        size_t control_union_count=0;
        if (channel=="SITE"){
            set<string> union_keys(recipient_universe.unit_keys.begin(),
                                   recipient_universe.unit_keys.end());
            union_keys.insert(source_universe.unit_keys.begin(),
                              source_universe.unit_keys.end());
            control_union_count=union_keys.size();
        } else control_union_count=recipient_universe.unit_keys.size()+
            source_universe.unit_keys.size();
        for (long double fraction : joint_parse_fractions(
                joint_map_value(control,"requested_fractions"))){
            const size_t total=min(recipient_universe.unit_keys.size(),
                                   source_universe.unit_keys.size());
            const size_t requested_source=joint_python_round_nonnegative(
                (long double)total*fraction);
            const size_t requested_recipient=total-requested_source;
            const vector<uint8_t> selected_recipient=joint_deterministic_unit_mask(
                recipient_universe,requested_recipient,
                seed+"|"+control_class+"|"+
                JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION+"|"+task.library+"|"+
                task.modality+"|"+channel+"|"+control_id+"|"+
                axis_fmt(fraction)+"|RECIPIENT");
            const vector<uint8_t> selected_source=joint_deterministic_unit_mask(
                source_universe,requested_source,
                seed+"|"+control_class+"|"+
                JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION+"|"+task.library+"|"+
                task.modality+"|"+channel+"|"+control_id+"|"+
                axis_fmt(fraction)+"|SOURCE");
            vector<pair<size_t,JointCompiledFit>> ranked;
            vector<pair<size_t,size_t>> realized(recipient_candidates.size());
            vector<JointCompiledFit> candidate_fits(recipient_candidates.size());
            if (basis_compatible){
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
              for (long long signed_index=0;
                    signed_index<(long long)recipient_candidates.size();
                    ++signed_index){
                const size_t index=static_cast<size_t>(signed_index);
                size_t recipient_count=0,source_count=0;
                const JointStoredUnits mixed=joint_control_candidate_units(
                    joint_channel_units_ref(recipient_candidates[index],channel),
                    joint_channel_units_ref(source_candidates[index],channel),
                    selected_recipient,selected_source,recipient_count,source_count,
                    channel);
                realized[index]=make_pair(recipient_count,source_count);
                JointCompiledFit fit=scalar_reference ? joint_reference_fit(
                    mixed,recipient_candidates[index].hypothesis,channel,
                    e_ref,e_alt,task.min_evidence,task.max_second_fraction) :
                    joint_compiled_fit(mixed,recipient_candidates[index].hypothesis,
                    channel,e_ref,e_alt,task.min_evidence,
                    task.max_second_fraction);
                candidate_fits[index]=fit;
              }
              for (size_t index=0;index<candidate_fits.size();++index)
                if (candidate_fits[index].status=="AVAILABLE")
                    ranked.push_back(make_pair(index,candidate_fits[index]));
            }
            sort(ranked.begin(),ranked.end(),[&](
                    const pair<size_t,JointCompiledFit>& left,
                    const pair<size_t,JointCompiledFit>& right){
                if (left.second.delta!=right.second.delta)
                    return left.second.delta>right.second.delta;
                return recipient_candidates[left.first].hypothesis.candidate_id>
                    recipient_candidates[right.first].hypothesis.candidate_id;
            });
            const long double winner_margin=ranked.empty() ? NAN :
                (ranked.size()>1 ? ranked[0].second.delta-ranked[1].second.delta :
                 numeric_limits<long double>::infinity());
            map<size_t,size_t> rank;
            for (size_t i=0;i<ranked.size();++i) rank[ranked[i].first]=i+1;
            size_t expected_index=recipient_candidates.size();
            size_t decoy_index=recipient_candidates.size();
            for (size_t i=0;i<recipient_candidates.size();++i){
                if (recipient_candidates[i].hypothesis.second_state==expected)
                    expected_index=i;
                if (!decoy.empty() && recipient_candidates[i].hypothesis.second_state==decoy)
                    decoy_index=i;
            }
            // Emit every legal candidate, including low/limited/unavailable
            // candidates.  Equivalence therefore compares complete menus,
            // not merely the candidates accepted by both engines.
            for (size_t index=0;index<recipient_candidates.size();++index){
                const JointCompiledFit& fit=candidate_fits[index];
                const string key="control_candidate:"+control_id+":"+channel+":"+
                    axis_fmt(fraction)+":"+
                    recipient_candidates[index].hypothesis.candidate_id;
                JointAnalysisRow row=joint_analysis_base(
                    task,"SOURCE_DISJOINT_CONTROL_CANDIDATE",key,threads);
                row.value.insert({
                    {"barcode",recipient_barcode},{"control_id",control_id},
                    {"control_class",control_class},{"evidence_channel",channel},
                    {"candidate_id",recipient_candidates[index].hypothesis.candidate_id},
                    {"second_state",recipient_candidates[index].hypothesis.second_state},
                    {"fraction",axis_fmt(fraction)},
                    {"status",basis_compatible ? fit.status :
                        "UNAVAILABLE_INCOMPATIBLE_MOLECULE_BASIS"},
                    {"legal_menu_candidates",to_string(recipient_candidates.size())},
                    {"rank",rank.count(index) ? to_string(rank[index]) : "NA"},
                    {"winner",(!ranked.empty() && ranked.front().first==index) ? "TRUE" : "FALSE"},
                    {"locked_log_likelihood",axis_fmt(fit.locked)},
                    {"interior_log_likelihood",axis_fmt(fit.interior)},
                    {"contributor_only_log_likelihood",axis_fmt(fit.contributor)},
                    {"delta_log_likelihood",axis_fmt(fit.delta)},
                    {"fitted_fraction",axis_fmt(fit.alpha)},
                    {"fitted_fraction_profile_low",axis_fmt(fit.alpha_low)},
                    {"fitted_fraction_profile_high",axis_fmt(fit.alpha_high)},
                    {"preferred_model",fit.preferred},
                    {"evidence_category",joint_fit_category(
                        fit,rank.count(index) && rank[index]==1 ?
                        winner_margin : NAN)},
                    {"calibrated_support_status",
                     "CALIBRATED_SUPPORT_NOT_EVALUATED"},
                    {"recipient_units",to_string(realized[index].first)},
                    {"source_units",to_string(realized[index].second)},
                    {"successful_source_units",to_string(realized[index].second)},
                    {"planned_source_fraction",axis_fmt(fraction)},
                    {"recipient_evidence_basis",channel=="SITE" ? "SITE_COUNTS" :
                        joint_basis_set_text(recipient_bases)},
                    {"source_evidence_basis",channel=="SITE" ? "SITE_COUNTS" :
                        joint_basis_set_text(source_bases)},
                    {"winner_margin",axis_fmt(winner_margin)},
                    {"raw_delta_margin",axis_fmt(winner_margin)},
                    {"realized_source_fraction",
                     realized[index].first+realized[index].second ? axis_fmt(
                        (long double)realized[index].second/
                        (long double)(realized[index].first+realized[index].second)) : "NA"}
                });
                joint_add_fit_contract_fields(row,fit);
                if (scalar_reference) row.value["engine"]="REFERENCE_SCALAR";
                output.push_back(row);
            }
            const string key="control_summary:"+control_id+":"+channel+":"+
                axis_fmt(fraction);
            JointAnalysisRow summary=joint_analysis_base(
                task,"SOURCE_DISJOINT_CONTROL_SUMMARY",key,threads);
            const size_t winner=ranked.empty() ? recipient_candidates.size() :
                ranked.front().first;
            const bool addition_supported=winner<recipient_candidates.size() &&
                joint_fit_category(ranked.front().second,winner_margin)==
                    "FIT_INTERIOR_ELIGIBLE";
            const size_t winner_recipient=winner<realized.size() ?
                realized[winner].first : 0;
            const size_t winner_source=winner<realized.size() ?
                realized[winner].second : 0;
            summary.value.insert({
                {"barcode",recipient_barcode},{"control_id",control_id},
                {"control_class",control_class},{"evidence_channel",channel},
                {"fraction",axis_fmt(fraction)},
                {"status",!basis_compatible ? "UNAVAILABLE_INCOMPATIBLE_MOLECULE_BASIS" :
                    ranked.empty() ? "UNAVAILABLE_EVIDENCE" : "AVAILABLE"},
                {"legal_menu_candidates",to_string(recipient_candidates.size())},
                {"full_menu_winner",winner<recipient_candidates.size() ?
                    recipient_candidates[winner].hypothesis.second_state : "NA"},
                {"expected_contributor",expected},
                {"expected_contributor_rank",rank.count(expected_index) ?
                    to_string(rank[expected_index]) : "NA"},
                {"decoy_contributor",decoy.empty() ? "NA" : decoy},
                {"decoy_rank",rank.count(decoy_index) ?
                    to_string(rank[decoy_index]) : "NA"},
                {"expected_recovered",winner==expected_index && addition_supported ?
                    "TRUE" : "FALSE"},
                {"decoy_won",winner==decoy_index &&
                    decoy_index<recipient_candidates.size() && addition_supported ?
                    "TRUE" : "FALSE"},
                {"preferred_model",ranked.empty() ? "UNAVAILABLE" :
                    ranked.front().second.preferred},
                {"evidence_category",ranked.empty() ? "UNAVAILABLE" :
                    joint_fit_category(ranked.front().second,winner_margin)},
                {"calibrated_support_status",
                 "CALIBRATED_SUPPORT_NOT_EVALUATED"},
                {"fitted_fraction",ranked.empty() ? "NA" :
                    axis_fmt(ranked.front().second.alpha)},
                {"planned_source_fraction",axis_fmt(fraction)},
                {"recipient_units",to_string(winner_recipient)},
                {"source_units",to_string(winner_source)},
                {"successful_source_units",to_string(winner_source)},
                {"realized_source_fraction",winner_recipient+winner_source ?
                    axis_fmt((long double)winner_source/
                             (long double)(winner_recipient+winner_source)) : "NA"},
                {"recipient_evidence_basis",channel=="SITE" ? "SITE_COUNTS" :
                    joint_basis_set_text(recipient_bases)},
                {"source_evidence_basis",channel=="SITE" ? "SITE_COUNTS" :
                    joint_basis_set_text(source_bases)},
                {"winner_margin",axis_fmt(winner_margin)},
                {"raw_delta_margin",axis_fmt(winner_margin)},
                {"complete_menu_unit_count",to_string(
                    control_union_count)}
            });
            if (!ranked.empty())
                joint_add_fit_contract_fields(summary,ranked.front().second);
            if (scalar_reference) summary.value["engine"]="REFERENCE_SCALAR";
            output.push_back(summary);
        }
    }
    return output;
}

static void joint_require_cache_metadata(
        const JointAnalysisTask& task,const string& manifest_digest){
    const string path=task.cache_prefix+".metadata.json";
    ifstream input(path.c_str());
    if (!input) throw runtime_error("missing normalized cache metadata: "+path);
    string content((istreambuf_iterator<char>(input)),istreambuf_iterator<char>());
    const vector<pair<string,string>> required={
        {"schema_version",JOINT_NORMALIZED_CACHE_SCHEMA},
        {"scientific_method_version",JOINT_DOUBLET_SCIENTIFIC_METHOD_VERSION},
        {"workload_generation_id",task.cache_generation},
        {"library",task.library},{"modality",task.modality},
        {"manifest_digest",manifest_digest}
    };
    for (const auto& item : required){
        const string token="\""+item.first+"\": \""+
            joint_json_escape(item.second)+"\"";
        if (content.find(token)==string::npos)
            throw runtime_error("normalized cache metadata mismatch for "+item.first);
    }
    const vector<pair<string,string>> numeric_required={
        {"error_ref",axis_fmt(task.error_ref)},
        {"error_alt",axis_fmt(task.error_alt)},
        {"min_evidence",to_string(task.min_evidence)},
        {"max_second_fraction",axis_fmt(task.max_second_fraction)}
    };
    for (const auto& item : numeric_required)
        if (content.find("\""+item.first+"\": "+item.second)==string::npos)
            throw runtime_error("normalized cache model-parameter mismatch for "+
                                item.first);
    if (content.find("\"calibration_library\": 25")==string::npos)
        throw runtime_error("normalized cache calibration-library mismatch");
    if (content.find("\"source_scans\": {\"sites\": 1, \"observations\": 1, \"molecules\": 1}")==string::npos)
        throw runtime_error("normalized cache metadata lacks one-scan provenance");
    const vector<pair<string,string>> components={
        {"observations_binary",task.cache_prefix+".observations.bin"},
        {"molecules_binary",task.cache_prefix+".molecules.bin"},
        {"sites_binary",task.cache_prefix+".sites.bin"},
        {"cells_index",task.cache_prefix+".cells.tsv"},
        {"samples_dictionary",task.cache_prefix+".samples.tsv"},
        {"genotypes_dictionary",task.cache_prefix+".genotypes.tsv.gz"}
    };
    for (const auto& component : components){
        const string digest=joint_file_content_digest(component.second);
        const string token="\""+component.first+
            "_content_digest\": \""+digest+"\"";
        if (content.find(token)==string::npos)
            throw runtime_error("normalized cache component digest mismatch for "+
                                component.first);
    }
}

static int joint_doublet_v4_self_test(){
    const long double tolerance=1e-9L;
    const auto require=[](bool condition,const string& message){
        if (!condition) throw runtime_error(
            "joint-doublet-v4 self-test failed: "+message);
    };
    const auto close=[&](long double observed,long double expected){
        return isfinite(observed) &&
            fabsl(observed-expected)<=tolerance*max(
                1.0L,max(fabsl(observed),fabsl(expected)));
    };
    JointHypothesis hypothesis;
    hypothesis.library="lib12"; hypothesis.barcode="ORACLE";
    hypothesis.candidate_id="ORACLE";
    hypothesis.locked.valid=true; hypothesis.second.valid=true;
    hypothesis.rho_effective=0.0L;
    const auto site_unit=[&](int position,uint32_t id){
        JointSiteUnit site;
        site.tid=0; site.pos=position; site.contig="chr1";
        site.ref_allele="A"; site.alt_allele="G";
        site.ref=9.0L; site.alt=1.0L;
        site.q_locked=0.0L; site.q_second=1.0L; site.q_ambient=0.0L;
        site.probability_intercept=0.001L;
        site.probability_slope=0.998L;
        site.compiled_site_id=id;
        JointLinkedUnit unit;
        unit.molecule=position; unit.compiled_unit_id=id;
        unit.genotype_distinguishable=1; unit.sites.push_back(site);
        return unit;
    };
    const auto molecule_unit=[&](const string& molecule,uint32_t id){
        JointLinkedUnit unit=site_unit(100+(int)id,id);
        unit.molecule=id+1; unit.molecule_id=molecule;
        unit.basis=1; unit.parent_origin="RECIPIENT";
        return unit;
    };
    const auto check_numeric=[&](const JointCompiledFit& fit,
                                 const string& label){
        require(close(fit.alpha,0.09919839679358718L),label+" alpha");
        require(close(fit.locked,-0.6916759781984388L),label+" locked");
        require(close(fit.interior,-0.3250829733914482L),label+" interior");
        require(close(fit.contributor,-6.217079801117281L),
                label+" contributor-at-one");
        require(close(fit.delta,0.36659300480699064L),label+" delta");
        require(fit.preferred=="INTERIOR_LOCKED_PLUS_CONTRIBUTOR",
                label+" preferred model");
    };
    vector<JointLinkedUnit> one_site{site_unit(100,0)};
    JointCompiledFit one_site_fit=joint_compiled_fit(
        one_site,hypothesis,"SITE",0.001L,0.001L,10,0.95L,NULL,NULL,false);
    require(one_site_fit.status=="LOW_EVIDENCE" &&
            joint_fit_category(one_site_fit,NAN,false)==
                "LOW_OR_LIMITED_EVIDENCE","one-site edge state");
    check_numeric(one_site_fit,"one-site");
    vector<JointLinkedUnit> two_site{site_unit(100,0),site_unit(101,1)};
    JointCompiledFit two_site_fit=joint_compiled_fit(
        two_site,hypothesis,"SITE",0.001L,0.001L,10,0.95L,NULL,NULL,false);
    JointCompiledFit two_site_reference=joint_reference_fit(
        two_site,hypothesis,"SITE",0.001L,0.001L,10,0.95L,NULL,NULL,false);
    require(two_site_fit.status=="AVAILABLE" &&
            joint_fit_category(two_site_fit,NAN,false)==
                "WEAK_ADDITION_COMPATIBLE","two-site edge state");
    check_numeric(two_site_fit,"two-site optimized");
    check_numeric(two_site_reference,"two-site reference");
    vector<JointLinkedUnit> one_molecule{molecule_unit("M001",0)};
    JointCompiledFit one_molecule_fit=joint_compiled_fit(
        one_molecule,hypothesis,"MOLECULE",0.001L,0.001L,10,0.95L,
        NULL,NULL,false);
    require(one_molecule_fit.status=="LIMITED_EVIDENCE" &&
            joint_fit_category(one_molecule_fit,NAN,false)==
                "LOW_OR_LIMITED_EVIDENCE","one-molecule edge state");
    check_numeric(one_molecule_fit,"one-molecule");
    vector<JointLinkedUnit> two_molecule{
        molecule_unit("M001",0),molecule_unit("M005",1)};
    JointCompiledFit two_molecule_fit=joint_compiled_fit(
        two_molecule,hypothesis,"MOLECULE",0.001L,0.001L,10,0.95L,
        NULL,NULL,false);
    JointCompiledFit two_molecule_reference=joint_reference_fit(
        two_molecule,hypothesis,"MOLECULE",0.001L,0.001L,10,0.95L,
        NULL,NULL,false);
    require(two_molecule_fit.status=="AVAILABLE" &&
            joint_fit_category(two_molecule_fit,NAN,false)==
                "WEAK_ADDITION_COMPATIBLE","two-molecule edge state");
    check_numeric(two_molecule_fit,"two-molecule optimized");
    check_numeric(two_molecule_reference,"two-molecule reference");
    vector<JointLinkedUnit> equivalent={site_unit(100,0),site_unit(101,1)};
    for (JointLinkedUnit& unit : equivalent){
        unit.genotype_distinguishable=0;
        unit.sites[0].q_second=unit.sites[0].q_locked;
        unit.sites[0].probability_slope=0.0L;
    }
    JointCompiledFit site_equivalent=joint_compiled_fit(
        equivalent,hypothesis,"SITE",0.001L,0.001L,10,0.95L);
    JointCompiledFit molecule_equivalent=joint_compiled_fit(
        equivalent,hypothesis,"MOLECULE",0.001L,0.001L,10,0.95L);
    require(site_equivalent.status=="GENOTYPE_EQUIVALENT" &&
            isfinite(site_equivalent.locked) && site_equivalent.delta==0.0L &&
            joint_fit_category(site_equivalent,NAN,false)=="LOCKED_SOURCE_ONLY",
            "site genotype-equivalent contract");
    require(molecule_equivalent.status=="GENOTYPE_EQUIVALENT" &&
            !isfinite(molecule_equivalent.locked) &&
            joint_fit_category(molecule_equivalent,NAN,false)=="LOCKED_SOURCE_ONLY",
            "molecule genotype-equivalent contract");
    vector<JointLinkedUnit> empty;
    JointCompiledFit empty_site=joint_compiled_fit(
        empty,hypothesis,"SITE",0.001L,0.001L,10,0.95L);
    JointCompiledFit empty_molecule=joint_compiled_fit(
        empty,hypothesis,"MOLECULE",0.001L,0.001L,10,0.95L);
    require(empty_site.status=="NO_OBSERVATIONS" &&
            !isfinite(empty_site.alpha),"empty SITE contract");
    require(empty_molecule.status=="UNAVAILABLE" &&
            empty_molecule.unavailable_reason=="NO_USABLE_LINKED_UNITS" &&
            !isfinite(empty_molecule.alpha),"empty MOLECULE contract");
    vector<int> positions={102,104,100,103,109,111,101,105,106,110};
    vector<JointLinkedUnit> ten_site;
    for (size_t i=0;i<positions.size();++i)
        ten_site.push_back(site_unit(positions[i],(uint32_t)i));
    JointCompiledFit ten_site_fit=joint_compiled_fit(
        ten_site,hypothesis,"SITE",0.001L,0.001L,10,0.95L);
    require(ten_site_fit.folds_evaluable==5 &&
            ten_site_fit.fold_support_numerator==5 &&
            close(ten_site_fit.fold_support_fraction,1.0L) &&
            close(ten_site_fit.alpha_low,0.0095L) &&
            close(ten_site_fit.alpha_high,0.3705L) &&
            close(ten_site_fit.maximum_influence_fraction,0.1L) &&
            joint_fit_category(ten_site_fit,NAN,false)=="FIT_INTERIOR_ELIGIBLE",
            "ten-site fold/profile/influence oracle");
    const vector<string> molecule_ids={
        "M001","M005","M004","M012","M003",
        "M015","M002","M010","M007","M008"};
    vector<JointLinkedUnit> ten_molecule;
    for (size_t i=0;i<molecule_ids.size();++i)
        ten_molecule.push_back(molecule_unit(molecule_ids[i],(uint32_t)i));
    JointCompiledFit ten_molecule_fit=joint_compiled_fit(
        ten_molecule,hypothesis,"MOLECULE",0.001L,0.001L,10,0.95L);
    require(ten_molecule_fit.folds_evaluable==5 &&
            ten_molecule_fit.fold_support_numerator==5 &&
            close(ten_molecule_fit.fold_support_fraction,1.0L) &&
            close(ten_molecule_fit.alpha_low,0.0095L) &&
            close(ten_molecule_fit.alpha_high,0.3705L) &&
            close(ten_molecule_fit.maximum_influence_fraction,0.1L) &&
            joint_fit_category(ten_molecule_fit,NAN,false)==
                "FIT_INTERIOR_ELIGIBLE",
            "ten-molecule fold/profile/influence oracle: evaluable="+
            to_string(ten_molecule_fit.folds_evaluable)+" numerator="+
            to_string(ten_molecule_fit.fold_support_numerator)+" profile="+
            axis_fmt(ten_molecule_fit.alpha_low)+","+
            axis_fmt(ten_molecule_fit.alpha_high)+" influence="+
            axis_fmt(ten_molecule_fit.maximum_influence_fraction)+" category="+
            joint_fit_category(ten_molecule_fit,NAN,false));
    require(joint_python_round_nonnegative(0.5L)==0 &&
            joint_python_round_nonnegative(1.5L)==2 &&
            joint_python_round_nonnegative(2.5L)==2,
            "Python half-even rounding");
    cout << "JOINT_DOUBLET_V4_SELF_TEST_PASS "
         << "alpha=" << axis_fmt(ten_site_fit.alpha)
         << " locked=" << axis_fmt(ten_site_fit.locked)
         << " interior=" << axis_fmt(ten_site_fit.interior)
         << " contributor=" << axis_fmt(ten_site_fit.contributor)
         << " delta=" << axis_fmt(ten_site_fit.delta)
         << " profile=[" << axis_fmt(ten_site_fit.alpha_low) << ","
         << axis_fmt(ten_site_fit.alpha_high) << "] influence="
         << axis_fmt(ten_site_fit.maximum_influence_fraction)
         << " folds=" << ten_site_fit.fold_support_numerator << "/"
         << ten_site_fit.folds_evaluable << endl;
    return 0;
}

static void run_joint_doublet_targeted_analysis(
        const string& task_manifest_path,int task_index,
        const string& cache_prefix_argument,const string& candidate_manifest,
        const string& control_manifest,const string& output_path,
        uint64_t master_seed,int downsample_replicates,int null_replicates,
        long double cli_e_ref,long double cli_e_alt,long cli_min_evidence,
        long double cli_max_alpha,bool reference_mode,int worker_threads){
    if (downsample_replicates!=100 || null_replicates!=1000)
        throw runtime_error("production targeted analysis requires exactly 100 downsample and 1000 null replicates");
    if (worker_threads<1) throw runtime_error("targeted analysis threads must be positive");
    JointAnalysisTask task=joint_load_analysis_task(task_manifest_path,task_index);
    joint_spill_root=output_path+".work";
    const long double parameter_tolerance=1e-15L;
    if (fabsl(task.error_ref-cli_e_ref)>parameter_tolerance ||
            fabsl(task.error_alt-cli_e_alt)>parameter_tolerance ||
            task.min_evidence!=cli_min_evidence ||
            fabsl(task.max_second_fraction-cli_max_alpha)>parameter_tolerance)
        throw runtime_error("CLI model parameters do not match analysis task contract");
    if (task.cache_prefix!=cache_prefix_argument)
        throw runtime_error("analysis cache prefix does not match task manifest");
    vector<map<string,string>> control_rows_manifest;
    if (task.action=="CONTROL_BIN"){
        control_rows_manifest=joint_read_tsv_maps(control_manifest);
        for (const auto& row : control_rows_manifest){
            const string library=joint_map_value(row,"library");
            if ((library!="lib7" && library!="lib9" && library!="lib12" && library!="lib17" && library!="lib20" && library!="lib25" && library!="lib29"))
                throw runtime_error(
                    "control manifest rejected non-target library: "+library);
        }
    }
    const vector<map<string,string>> extraction_rows=
        joint_read_tsv_maps(task_manifest_path.substr(0,task_manifest_path.find_last_of('/'))+
                            "/extraction_tasks.tsv");
    string manifest_digest;
    for (const auto& row : extraction_rows)
        if (joint_map_value(row,"library")==task.library &&
                joint_map_value(row,"modality")==task.modality){
            manifest_digest=joint_map_value(row,"manifest_digest"); break;
        }
    if (manifest_digest.empty())
        throw runtime_error("no matching extraction task for targeted analysis");
    joint_require_cache_metadata(task,manifest_digest);

    unordered_map<string,int> sample2idx;
    unordered_map<int,size_t> donor_slot;
    joint_load_normalized_samples(task.cache_prefix,sample2idx,donor_slot);
    vector<int> donors;
    const vector<AxisSiteDefinition> sites=joint_load_normalized_sites(
        task.cache_prefix,donors,task.worker_memory_budget_bytes/2);
    if (donors.size()!=donor_slot.size())
        throw runtime_error("normalized cache donor dictionaries disagree");
    for (size_t slot=0;slot<donors.size();++slot)
        if (!donor_slot.count(donors[slot]) || donor_slot.at(donors[slot])!=slot)
            throw runtime_error("normalized cache donor slot mismatch");
    unordered_map<unsigned long,vector<size_t>> by_cell;
    JointStreamingDigest current_manifest_digest;
    const vector<JointHypothesis> hypotheses=joint_load_manifest(
        candidate_manifest,sample2idx,task.library,by_cell,
        &current_manifest_digest);
    if (current_manifest_digest.text()!=manifest_digest)
        throw runtime_error("candidate manifest changed after cache extraction");
    const map<unsigned long,JointNormalizedCellIndex> cell_index=
        joint_load_normalized_index(task.cache_prefix,task);
    const long double e_ref=task.error_ref;
    const long double e_alt=task.error_alt;
    size_t largest_menu=1;
    for(const auto& item:by_cell)largest_menu=max(largest_menu,item.second.size());
    joint_spill_limit=max<size_t>(65536,min<uint64_t>(16ULL*1024*1024,
        task.worker_memory_budget_bytes/(32*largest_menu)));
    uint64_t dictionary_bytes=sites.capacity()*sizeof(AxisSiteDefinition);
    for(const auto& site:sites)dictionary_bytes+=site.genotype.capacity()*sizeof(site.genotype[0])+
        site.contig.capacity()+site.ref_allele.capacity()+site.alt_allele.capacity();
    const auto memory_bound=[&](const string& barcode){
        string text=barcode;auto found=cell_index.find(bc_ul(text));
        if(found==cell_index.end())throw runtime_error("cell absent from indexed evidence: "+barcode);
        const uint64_t records=found->second.observation_count+found->second.molecule_count;
        // Pessimistic strings/maps, immutable offset indexes, masks/worker
        // counters, two channels and temporary rebuilt coefficient stores.
        return records*(2048ULL+32ULL*largest_menu+128ULL*worker_threads)+
            8ULL*largest_menu*joint_spill_limit+128ULL*1024*1024;
    };
    uint64_t largest_task=0;
    if(task.action=="CELL")largest_task=memory_bound(task.barcode);
    else for(const auto& control:control_rows_manifest){
        const string id=joint_map_value(control,"control_id");
        const vector<string> ids=split(task.control_ids,',');
        if(find(ids.begin(),ids.end(),id)!=ids.end())largest_task=max<uint64_t>(largest_task,
            memory_bound(joint_map_value(control,"recipient_barcode"))+
            memory_bound(joint_map_value(control,"source_barcode")));
    }
    if(largest_task+dictionary_bytes>task.worker_memory_budget_bytes)
        throw runtime_error("bounded targeted metadata/workers exceed the configured memory budget; reduce worker count or increase the reviewed budget");
    auto load_cell=[&](const string& barcode)->JointCachedCell{
        joint_current_observed_pool.reset(new JointObservedPool);
        if (barcode.size()!=16)
            throw runtime_error("targeted analysis barcode must be canonical 16 bp");
        string mutable_barcode=barcode;
        const unsigned long encoded=bc_ul(mutable_barcode);
        auto position=cell_index.find(encoded);
        if (position==cell_index.end())
            throw runtime_error("requested barcode is missing from normalized cache index: "+barcode);
        const JointNormalizedCellIndex& item=position->second;
        const JointRecordView<AxisObservationRecord> observations(
                task.cache_prefix+".observations.bin",item.observation_offset,
                item.observation_count);
        const JointRecordView<JointMoleculeRecord> molecules(
                task.cache_prefix+".molecules.bin",item.molecule_offset,
                item.molecule_count);
        JointCachedCell cell;
        cell.raw_observations=observations;
        cell.raw_molecules=molecules;
        cell.malformed_molecule_rows=item.malformed_molecule_rows;
        cell.candidates=joint_evaluate_cached_cell(
            encoded,hypotheses,by_cell,observations,molecules,sites,donor_slot,
            item.malformed_molecule_rows,e_ref,e_alt,task.min_evidence,
            task.max_second_fraction,reference_mode);
        return cell;
    };

    joint_targeted_optimizer_calls=0;
    joint_targeted_likelihood_evaluations=0;
    joint_targeted_derivative_evaluations=0;
    joint_targeted_optimized_fit_calls=0;
    joint_targeted_reference_fit_calls=0;
    joint_targeted_row_scans=0;
    joint_count_targeted_work=true;
    vector<JointAnalysisRow> rows;
    if (task.action=="CELL"){
        JointCachedCell cell=load_cell(task.barcode);
        vector<JointCompiledCandidate>& candidates=cell.candidates;
        const JointCompiledUniverse site_universe=joint_compile_channel_universe(
            candidates,"SITE",e_ref,e_alt);
        const JointCompiledUniverse molecule_universe=joint_compile_channel_universe(
            candidates,"MOLECULE",e_ref,e_alt);
        vector<JointAnalysisRow> observed=joint_observed_rows(
            task,candidates,task.barcode,e_ref,e_alt,worker_threads);
        if (reference_mode){
            joint_assert_compiled_reference_equivalence(
                candidates,"SITE",e_ref,e_alt,task.min_evidence,
                task.max_second_fraction);
            joint_assert_compiled_reference_equivalence(
                candidates,"MOLECULE",e_ref,e_alt,task.min_evidence,
                task.max_second_fraction);
            vector<JointAnalysisRow> reference=joint_observed_rows(
                task,candidates,task.barcode,e_ref,e_alt,worker_threads,true);
            observed.insert(observed.end(),reference.begin(),reference.end());
        }
        vector<JointAnalysisRow> downsample=joint_downsample_rows(
            task,candidates,site_universe,molecule_universe,task.barcode,
            master_seed,downsample_replicates,
            e_ref,e_alt,worker_threads);
        vector<JointAnalysisRow> nulls=joint_null_rows(
            task,candidates,site_universe,molecule_universe,task.barcode,
            master_seed,null_replicates,
            e_ref,e_alt,worker_threads);
        rows.insert(rows.end(),make_move_iterator(observed.begin()),make_move_iterator(observed.end()));
        rows.insert(rows.end(),make_move_iterator(downsample.begin()),make_move_iterator(downsample.end()));
        rows.insert(rows.end(),make_move_iterator(nulls.begin()),make_move_iterator(nulls.end()));
    } else {
        map<string,map<string,string>> by_id;
        for (const auto& row : control_rows_manifest)
            by_id[joint_map_value(row,"control_id")]=row;
        for (const string& raw_id : split(task.control_ids,',')){
            const string control_id=trim(raw_id);
            if (control_id.empty()) continue;
            auto found=by_id.find(control_id);
            if (found==by_id.end())
                throw runtime_error("control id absent from finalized manifest: "+control_id);
            const string recipient=joint_map_value(found->second,"recipient_barcode");
            const string source=joint_map_value(found->second,"source_barcode");
            JointCachedCell recipient_cell=load_cell(recipient);
            JointCachedCell source_cell=load_cell(source);
            vector<JointAnalysisRow> control_rows=joint_control_rows(
                task,found->second,recipient_cell,source_cell,sites,donor_slot,
                e_ref,e_alt,worker_threads);
            rows.insert(rows.end(),control_rows.begin(),control_rows.end());
            if (reference_mode){
                vector<JointAnalysisRow> reference_rows=joint_control_rows(
                    task,found->second,recipient_cell,source_cell,sites,donor_slot,
                    e_ref,e_alt,worker_threads,true);
                rows.insert(rows.end(),reference_rows.begin(),reference_rows.end());
            }
        }
        if (rows.empty()) throw runtime_error("control analysis task selected no controls");
    }
    joint_count_targeted_work=false;
    JointAnalysisRow metrics=joint_analysis_base(
        task,"COMPILED_SHARD_ACCOUNTING","metrics:"+to_string(task.task_index),
        worker_threads);
    metrics.value.insert({
        {"status","COMPLETE"},
        {"optimized_fit_calls",to_string(
            joint_targeted_optimized_fit_calls.load())},
        {"reference_fit_calls",to_string(
            joint_targeted_reference_fit_calls.load())},
        {"optimizer_calls",to_string(joint_targeted_optimizer_calls.load())},
        {"likelihood_evaluations",to_string(
            joint_targeted_likelihood_evaluations.load())},
        {"derivative_evaluations",to_string(
            joint_targeted_derivative_evaluations.load())},
        {"row_scans",to_string(joint_targeted_row_scans.load())},
        {"candidate_owned_raw_evidence_copies","0"},
        {"per_fraction_string_rehashes","0"},
        {"per_fraction_string_sorts","0"},
        {"unchanged_80_pass_refits","0"},
        {"reference_definition","existing scalar likelihood and derivative optimizer; normalized cache and deterministic unit masks"}
    });
    rows.push_back(metrics);

    vector<JointAnalysisRow> emitted;
    for (JointAnalysisRow& row : rows){
        if (!row.value.count("engine")) row.value["engine"]="OPTIMIZED_BATCHED";
        row.value["result_key"]=row.value["equivalence_key"]+":"+
            row.value["engine"];
        row.value["reference_definition"]=row.value["engine"]=="REFERENCE_SCALAR" ?
            "completed joint-doublet scalar scorer using identical cached observations" :
            "normalized-cache compiled batch with shared deterministic unit masks";
        emitted.push_back(move(row));
    }
    sort(emitted.begin(),emitted.end(),[](const JointAnalysisRow& left,
                                          const JointAnalysisRow& right){
        return left.value.at("result_key")<right.value.at("result_key");
    });
    const string temporary=output_path+".tmp."+to_string((long long)getpid());
    gzFile output=gzopen(temporary.c_str(),"wb");
    if (!output) throw runtime_error("could not create targeted analysis output: "+temporary);
    const vector<string> header=joint_analysis_header();
    try {
        axis_gzwrite(output,header);
        for (const JointAnalysisRow& row : emitted)
            axis_gzwrite(output,joint_analysis_values(row,header));
        const int closed=gzclose(output);output=NULL;
        if(closed!=Z_OK)throw runtime_error("failed closing targeted analysis output");
    } catch (...) {
        joint_count_targeted_work=false;
        if (output) gzclose(output);
        unlink(temporary.c_str());
        throw;
    }
    if (rename(temporary.c_str(),output_path.c_str())!=0){
        unlink(temporary.c_str());
        throw runtime_error("failed publishing targeted analysis output: "+output_path);
    }
}

static void run_joint_doublet(
        const string& samples_path, const string& manifest_path,
        const string& sites_path, const string& observations_path,
        const string& molecules_path,
        const string& output_path, const string& targeted_cache_path,
        const string& temp_root,
        const string& library, const string& modality,
        long double e_ref, long double e_alt, long min_evidence,
        long double max_alpha, int folds, int worker_threads){
    if (e_ref < 0.0L || e_ref > 1.0L || e_alt < 0.0L || e_alt > 1.0L ||
            e_ref+e_alt >= 1.0L)
        throw runtime_error("joint-doublet errors must be in [0,1] and sum to less than one");
    if (min_evidence < 0)
        throw runtime_error("--min_evidence must be nonnegative in joint-doublet mode");
    if (max_alpha <= 0.0L || max_alpha > 1.0L)
        throw runtime_error("--max-second-fraction must be in (0,1]");
    if (folds < 2)
        throw runtime_error("--joint-folds must be at least 2");
    if (worker_threads < 1)
        throw runtime_error("--threads must be positive in joint-doublet mode");
    if (modality != "RNA" && modality != "ATAC")
        throw runtime_error("--modality must be RNA or ATAC");
    const bool capture_cache=!targeted_cache_path.empty();

    const vector<string> samples=load_samples(samples_path);
    unordered_map<string,int> sample2idx;
    for (int i=0; i<(int)samples.size(); ++i){
        if (trim(samples[i]).empty() || sample2idx.count(samples[i]))
            throw runtime_error("joint-doublet samples must be nonblank and unique");
        sample2idx[samples[i]]=i;
    }
    unordered_map<unsigned long,vector<size_t>> by_cell;
    vector<JointHypothesis> hypotheses=joint_load_manifest(
        manifest_path,sample2idx,library,by_cell);

    unordered_map<unsigned long,CandidateAxisPair> targets;
    for (const auto& item : by_cell) targets[item.first]=CandidateAxisPair();
    vector<int> donors;
    for (const JointHypothesis& hypothesis : hypotheses){
        const JointComposition* compositions[] = {
            &hypothesis.locked,&hypothesis.second,&hypothesis.ambient};
        for (const JointComposition* composition : compositions)
            for (const auto& member : composition->members)
                donors.push_back(member.first);
    }
    sort(donors.begin(),donors.end());
    donors.erase(unique(donors.begin(),donors.end()),donors.end());
    unordered_map<int,size_t> donor_slot;
    for (size_t i=0; i<donors.size(); ++i) donor_slot[donors[i]]=i;

    AxisResourceAudit audit;
    vector<JointResult> results(hypotheses.size());
    if (!hypotheses.empty()){
        AxisTempGuard temporary(temp_root);
        const unsigned long long bucket_bytes=256ULL*1024ULL*1024ULL;
        typedef chrono::steady_clock JointClock;
        JointClock::time_point phase_start=JointClock::now();
        const auto report_phase = [&](const string& phase){
            const JointClock::time_point now=JointClock::now();
            cerr << "JOINT_DOUBLET_PHASE\t" << phase << "\tseconds="
                 << chrono::duration<double>(now-phase_start).count() << "\n";
            phase_start=now;
        };
        unordered_map<unsigned long,unsigned long long> rows_by_barcode;
        vector<uint64_t> selected_keys;
        const string spool_path=joint_stage_observations_one_pass(
            observations_path,targets,temporary.path(),audit,
            rows_by_barcode,selected_keys);
        report_phase("read_and_stage_observations");
        vector<AxisSiteDefinition> sites=axis_load_site_definitions(
            sites_path,selected_keys,(int)samples.size(),donors,audit);
        report_phase("load_selected_site_definitions");
        unordered_map<unsigned long,size_t> bucket_assignment;
        const vector<unsigned long long> bucket_loads=axis_assign_buckets(
            rows_by_barcode,bucket_bytes,bucket_assignment,audit,
            max<size_t>(1,static_cast<size_t>(worker_threads)*4));
        vector<string> bucket_paths=joint_partition_staged_observations(
            spool_path,bucket_assignment,bucket_loads.size(),temporary.path());
        unordered_map<unsigned long,unsigned long long> molecule_rows_by_barcode;
        unordered_map<unsigned long,long> malformed_molecule_rows;
        const string molecule_spool=joint_stage_molecules_one_pass(
            molecules_path,targets,selected_keys,temporary.path(),
            molecule_rows_by_barcode,malformed_molecule_rows);
        vector<string> molecule_bucket_paths=joint_partition_staged_molecules(
            molecule_spool,bucket_assignment,bucket_loads.size(),temporary.path());
        vector<size_t> bucket_rows(bucket_paths.size(),0);
        for (size_t bucket_index=0; bucket_index<bucket_paths.size();
                ++bucket_index){
            const string& path=bucket_paths[bucket_index];
            struct stat info;
            if (stat(path.c_str(),&info) != 0)
                throw runtime_error("could not stat joint-doublet bucket: " + path);
            if (info.st_size < 0 ||
                    static_cast<unsigned long long>(info.st_size) %
                        sizeof(AxisObservationRecord) != 0)
                throw runtime_error(
                    "joint-doublet bucket has a truncated binary record: " + path);
            const size_t n=static_cast<size_t>(info.st_size) /
                sizeof(AxisObservationRecord);
            bucket_rows[bucket_index]=n;
            audit.observed_peak_bucket_rows=max<unsigned long long>(
                audit.observed_peak_bucket_rows,n);
        }
        report_phase("partition_binary_buckets");

        vector<string> bucket_errors(bucket_paths.size());
#pragma omp parallel for schedule(dynamic,1) num_threads(worker_threads)
        for (long long raw_bucket_index=0;
                raw_bucket_index<static_cast<long long>(bucket_paths.size());
                ++raw_bucket_index){
            const size_t bucket_index=static_cast<size_t>(raw_bucket_index);
            const string& path=bucket_paths[bucket_index];
            try {
                const size_t n=bucket_rows[bucket_index];
                vector<AxisObservationRecord> records(n);
                ifstream input(path.c_str(),ios::binary);
                if (!input)
                    throw runtime_error(
                        "could not open joint-doublet bucket: " + path);
                if (n){
                    input.read(reinterpret_cast<char*>(records.data()),
                        n*sizeof(AxisObservationRecord));
                    if (!input)
                        throw runtime_error(
                            "failed reading joint-doublet bucket: " + path);
                }
                sort(records.begin(),records.end(),axis_observation_less);
                size_t begin=0;
                while (begin<records.size()){
                    const unsigned long barcode=records[begin].barcode;
                    size_t end=begin;
                    while (end<records.size() &&
                            records[end].barcode==barcode) ++end;
                    vector<AxisObservationRecord> merged;
                    size_t cursor=begin;
                    while (cursor<end){
                        AxisObservationRecord aggregate=records[cursor];
                        AxisKahan ref,alt;
                        while (cursor<end &&
                                records[cursor].tid==aggregate.tid &&
                                records[cursor].pos==aggregate.pos){
                            ref.add(records[cursor].ref);
                            alt.add(records[cursor].alt);
                            ++cursor;
                        }
                        aggregate.ref=(double)ref.value;
                        aggregate.alt=(double)alt.value;
                        merged.push_back(aggregate);
                    }
                    auto found=by_cell.find(barcode);
                    if (found==by_cell.end())
                        throw runtime_error(
                            "joint-doublet bucket contains nontarget barcode");
                    for (size_t index : found->second)
                        results[index]=joint_evaluate(
                            hypotheses[index],merged,sites,donor_slot,e_ref,
                            e_alt,min_evidence,max_alpha,folds,capture_cache);
                    begin=end;
                }
                unlink(path.c_str());
            } catch (const exception& error){
                bucket_errors[bucket_index]=error.what();
            }
        }
        for (size_t i=0; i<bucket_errors.size(); ++i)
            if (!bucket_errors[i].empty())
                throw runtime_error(
                    "joint-doublet bucket worker " + to_string(i) +
                    " failed: " + bucket_errors[i]);
        report_phase("parallel_bucket_scoring");

        const vector<AxisObservationRecord> empty;
        vector<size_t> empty_indexes;
        for (const auto& item : by_cell)
            if (rows_by_barcode.count(item.first)==0)
                for (size_t index : item.second)
                    empty_indexes.push_back(index);
        vector<string> empty_errors(empty_indexes.size());
#pragma omp parallel for schedule(dynamic,16) num_threads(worker_threads)
        for (long long raw_empty_index=0;
                raw_empty_index<static_cast<long long>(empty_indexes.size());
                ++raw_empty_index){
            const size_t work_index=static_cast<size_t>(raw_empty_index);
            const size_t result_index=empty_indexes[work_index];
            try {
                results[result_index]=joint_evaluate(
                    hypotheses[result_index],empty,sites,donor_slot,e_ref,e_alt,
                    min_evidence,max_alpha,folds,capture_cache);
            } catch (const exception& error){
                empty_errors[work_index]=error.what();
            }
        }
        for (size_t i=0; i<empty_errors.size(); ++i)
            if (!empty_errors[i].empty())
                throw runtime_error(
                    "joint-doublet empty-observation worker failed for result " +
                    to_string(empty_indexes[i]) + ": " + empty_errors[i]);
        report_phase("empty_observation_scoring");

        if (!molecule_bucket_paths.empty()){
            vector<string> molecule_bucket_errors(molecule_bucket_paths.size());
#pragma omp parallel for schedule(dynamic,1) num_threads(worker_threads)
            for (long long raw_bucket_index=0;
                    raw_bucket_index<static_cast<long long>(molecule_bucket_paths.size());
                    ++raw_bucket_index){
                const size_t bucket_index=static_cast<size_t>(raw_bucket_index);
                const string& path=molecule_bucket_paths[bucket_index];
                try {
                    struct stat info;
                    if (stat(path.c_str(),&info)!=0 || info.st_size<0 ||
                            static_cast<unsigned long long>(info.st_size)%
                                sizeof(JointMoleculeRecord)!=0)
                        throw runtime_error(
                            "joint molecule bucket is missing or truncated: "+path);
                    const size_t n=(size_t)info.st_size/
                        sizeof(JointMoleculeRecord);
                    vector<JointMoleculeRecord> records(n);
                    ifstream input(path.c_str(),ios::binary);
                    if (!input)
                        throw runtime_error("could not open joint molecule bucket: "+path);
                    if (n){
                        input.read(reinterpret_cast<char*>(records.data()),
                            n*sizeof(JointMoleculeRecord));
                        if (!input)
                            throw runtime_error("failed reading joint molecule bucket: "+path);
                    }
                    sort(records.begin(),records.end(),joint_molecule_less);
                    vector<JointMoleculeRecord> merged;
                    size_t cursor=0;
                    while (cursor<records.size()){
                        JointMoleculeRecord aggregate=records[cursor];
                        AxisKahan ref,alt;
                        while (cursor<records.size() &&
                                records[cursor].barcode==aggregate.barcode &&
                                records[cursor].basis==aggregate.basis &&
                                records[cursor].molecule==aggregate.molecule &&
                                records[cursor].tid==aggregate.tid &&
                                records[cursor].pos==aggregate.pos){
                            ref.add(records[cursor].ref);
                            alt.add(records[cursor].alt);
                            ++cursor;
                        }
                        aggregate.ref=(double)ref.value;
                        aggregate.alt=(double)alt.value;
                        merged.push_back(aggregate);
                    }
                    size_t begin=0;
                    while (begin<merged.size()){
                        const unsigned long barcode=merged[begin].barcode;
                        size_t end=begin+1;
                        while (end<merged.size() && merged[end].barcode==barcode) ++end;
                        auto found=by_cell.find(barcode);
                        if (found==by_cell.end())
                            throw runtime_error(
                                "joint molecule bucket contains nontarget barcode");
                        const vector<JointMoleculeRecord> cell_records(
                            merged.begin()+begin,merged.begin()+end);
                        const long malformed=malformed_molecule_rows.count(barcode) ?
                            malformed_molecule_rows.at(barcode) : 0;
                        for (size_t index : found->second)
                            joint_evaluate_molecules(
                                hypotheses[index],cell_records,sites,donor_slot,
                                e_ref,e_alt,max_alpha,malformed,capture_cache,
                                results[index]);
                        begin=end;
                    }
                    unlink(path.c_str());
                } catch (const exception& error){
                    molecule_bucket_errors[bucket_index]=error.what();
                }
            }
            for (size_t i=0;i<molecule_bucket_errors.size();++i)
                if (!molecule_bucket_errors[i].empty())
                    throw runtime_error(
                        "joint molecule bucket worker "+to_string(i)+
                        " failed: "+molecule_bucket_errors[i]);
            for (const auto& item : by_cell){
                if (molecule_rows_by_barcode.count(item.first)>0) continue;
                const long malformed=malformed_molecule_rows.count(item.first) ?
                    malformed_molecule_rows.at(item.first) : 0;
                const vector<JointMoleculeRecord> empty_molecules;
                for (size_t index : item.second)
                    joint_evaluate_molecules(
                        hypotheses[index],empty_molecules,sites,donor_slot,
                        e_ref,e_alt,max_alpha,malformed,capture_cache,
                        results[index]);
            }
            report_phase("parallel_molecule_bucket_scoring");
        }
    }

    const string temporary_output=output_path+".tmp."+
        to_string((long long)getpid());
    gzFile output=gzopen(temporary_output.c_str(),"wb");
    if (!output) throw runtime_error("could not create joint-doublet output: " + temporary_output);
    try {
        const vector<string> header=joint_output_header();
        axis_gzwrite(output,header);
        vector<size_t> order(hypotheses.size());
        iota(order.begin(),order.end(),0);
        sort(order.begin(),order.end(),[&](size_t left,size_t right){
            return make_pair(hypotheses[left].barcode,hypotheses[left].candidate_id) <
                make_pair(hypotheses[right].barcode,hypotheses[right].candidate_id);
        });
        for (size_t index : order){
            const vector<string> row=joint_output_row(
                hypotheses[index],results[index],modality,audit,e_ref,e_alt,
                min_evidence,max_alpha,folds);
            if (row.size()!=header.size())
                throw runtime_error("joint-doublet internal output schema/value count mismatch");
            axis_gzwrite(output,row);
        }
        if (gzclose(output)!=Z_OK)
            throw runtime_error("failed closing joint-doublet output: " + temporary_output);
        output=NULL;
    } catch (...) {
        if (output) gzclose(output);
        unlink(temporary_output.c_str());
        throw;
    }
    if (rename(temporary_output.c_str(),output_path.c_str())!=0){
        unlink(temporary_output.c_str());
        throw runtime_error("failed publishing joint-doublet output: " + output_path);
    }
    joint_write_targeted_cache(
        targeted_cache_path,hypotheses,results,modality,e_ref,e_alt);
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
        "  --threads N                    joint-doublet barcode-bucket workers (default 1)\n"
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
        "  certainty, accuracy, FDR, posterior probability, or correctness probability.\n"
        "\nStandalone generalized joint-doublet scoring:\n"
        "  --joint-doublet-output FILE    K=1 versus K=2 + ambient candidate rows\n"
        "  --joint-doublet-manifest FILE  joint_doublet_candidate_manifest_v4 TSV[.gz]\n"
        "  --joint-doublet-temp-dir DIR   existing absolute job-local temp root\n"
        "  --joint-doublet-targeted-evidence-cache FILE\n"
        "                                 optional reusable site/linked-unit cache\n"
        "  --modality RNA|ATAC            preserve assay-specific evidence\n"
        "  --pileup-molecules FILE         optional seven-column linked-unit sidecar\n"
        "  --max-second-fraction FLOAT    fitted range upper bound (default 0.95)\n"
        "  --joint-folds N                genomic leave-one-fold-out groups (default 5)\n"
        "  Joint-doublet mode also requires --samples, --pileup-sites,\n"
        "  --pileup-observations, and --libname. It never changes identity assignments.\n");
    fprintf(stderr,
        "\nEvidence-normalized targeted validation v4:\n"
        "  --joint-doublet-normalized-cache-prefix PREFIX\n"
        "  --joint-doublet-manifest-digest DIGEST\n"
        "  --joint-doublet-workload-generation ID\n"
        "      Extract one library/modality cache in one source scan; also requires\n"
        "      --joint-doublet-manifest, --samples, all three pileup inputs,\n"
        "      --joint-doublet-temp-dir, --libname, and --modality.\n"
        "  --joint-doublet-targeted-analysis-output FILE\n"
        "  --joint-doublet-targeted-analysis-manifest FILE\n"
        "  --joint-doublet-targeted-analysis-index N\n"
        "  --joint-doublet-control-manifest FILE\n"
        "  --joint-doublet-seed N\n"
        "  --joint-doublet-downsample-replicates 100\n"
        "  --joint-doublet-null-replicates 1000\n"
        "  --joint-doublet-reference-mode\n"
        "  --joint-doublet-v4-self-test  literal edge/numeric/fold oracle\n"
        "      Run one bounded cache-backed cell/control shard. No raw pileup is read.\n");
}

int main(int argc, char** argv){
    if(argc==3 && string(argv[1])=="--joint-doublet-cache-digest"){
        try{cout<<joint_file_content_digest(argv[2])<<endl;return 0;}
        catch(const exception& error){cerr<<error.what()<<endl;return 1;}
    }
    string counts, samples_path, assignments, diagnostics, output, runnerups, panel_path, libname="NA";
    string species_counts, species_condf, species_samples_path, condf;
    string candidate_manifest, pileup_sites, pileup_observations;
    string pileup_molecules, site_fold_output, probability_output, score_prefix;
    string candidate_axis_output, candidate_axis_temp_dir;
    string joint_doublet_output, joint_doublet_manifest;
    string joint_doublet_targeted_evidence_cache;
    string joint_doublet_temp_dir, modality;
    string joint_doublet_normalized_cache_prefix;
    string joint_doublet_manifest_digest, joint_doublet_workload_generation;
    string joint_doublet_targeted_analysis_output;
    string joint_doublet_targeted_analysis_manifest;
    string joint_doublet_control_manifest;
    int joint_doublet_targeted_analysis_index = -1;
    uint64_t joint_doublet_seed = 20260920ULL;
    int joint_doublet_downsample_replicates = 100;
    int joint_doublet_null_replicates = 1000;
    uint64_t joint_doublet_scheduler_memory_bytes = 0;
    uint64_t joint_doublet_launcher_runtime_reserve_bytes = 0;
    uint64_t joint_doublet_worker_memory_budget_bytes = 0;
    bool joint_doublet_reference_mode = false;
    int site_folds = 5;
    int probability_resamples = 100;
    uint64_t probability_seed = 1729;
    long min_evidence = 10;
    double e_ref = 0.001, e_alt = 0.001;
    double poor_fit_residual = 0.30;
    double species_support_threshold = 0.70;
    double max_second_fraction = 0.95;
    int joint_folds = 5;
    int threads = 1;
    bool strict = false;
    bool candidate_axis_self_test_requested = false;
    bool joint_doublet_v4_self_test_requested = false;
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
        else if (a == "--threads") { string tmp; need(tmp); threads = atoi(tmp.c_str()); explicit_threads = true; }
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
        else if (a == "--joint-doublet-v4-self-test")
            joint_doublet_v4_self_test_requested = true;
        else if (a == "--joint-doublet-output") need(joint_doublet_output);
        else if (a == "--joint-doublet-manifest") need(joint_doublet_manifest);
        else if (a == "--joint-doublet-temp-dir") need(joint_doublet_temp_dir);
        else if (a == "--joint-doublet-targeted-evidence-cache")
            need(joint_doublet_targeted_evidence_cache);
        else if (a == "--joint-doublet-normalized-cache-prefix")
            need(joint_doublet_normalized_cache_prefix);
        else if (a == "--joint-doublet-manifest-digest")
            need(joint_doublet_manifest_digest);
        else if (a == "--joint-doublet-workload-generation")
            need(joint_doublet_workload_generation);
        else if (a == "--joint-doublet-targeted-analysis-output")
            need(joint_doublet_targeted_analysis_output);
        else if (a == "--joint-doublet-targeted-analysis-manifest")
            need(joint_doublet_targeted_analysis_manifest);
        else if (a == "--joint-doublet-targeted-analysis-index") {
            string tmp; need(tmp); joint_doublet_targeted_analysis_index=atoi(tmp.c_str());
        }
        else if (a == "--joint-doublet-control-manifest")
            need(joint_doublet_control_manifest);
        else if (a == "--joint-doublet-seed") {
            string tmp; need(tmp); joint_doublet_seed=strtoull(tmp.c_str(),NULL,10);
        }
        else if (a == "--joint-doublet-downsample-replicates") {
            string tmp; need(tmp); joint_doublet_downsample_replicates=atoi(tmp.c_str());
        }
        else if (a == "--joint-doublet-null-replicates") {
            string tmp; need(tmp); joint_doublet_null_replicates=atoi(tmp.c_str());
        }
        else if (a == "--joint-doublet-scheduler-memory-bytes") {
            string tmp; need(tmp); joint_doublet_scheduler_memory_bytes=
                strtoull(tmp.c_str(),NULL,10);
        }
        else if (a == "--joint-doublet-launcher-runtime-reserve-bytes") {
            string tmp; need(tmp); joint_doublet_launcher_runtime_reserve_bytes=
                strtoull(tmp.c_str(),NULL,10);
        }
        else if (a == "--joint-doublet-worker-memory-budget-bytes") {
            string tmp; need(tmp); joint_doublet_worker_memory_budget_bytes=
                strtoull(tmp.c_str(),NULL,10);
        }
        else if (a == "--joint-doublet-reference-mode")
            joint_doublet_reference_mode=true;
        else if (a == "--modality") need(modality);
        else if (a == "--max-second-fraction") { string tmp; need(tmp); max_second_fraction = atof(tmp.c_str()); }
        else if (a == "--joint-folds") { string tmp; need(tmp); joint_folds = atoi(tmp.c_str()); }
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
    if (joint_doublet_v4_self_test_requested){
        try { return joint_doublet_v4_self_test(); }
        catch (const exception& error){
            fprintf(stderr, "ERROR: %s\n", error.what());
            return 1;
        }
    }
    if (!joint_doublet_targeted_analysis_output.empty()){
        const bool incompatible=!joint_doublet_output.empty() || !counts.empty() ||
            !assignments.empty() || !diagnostics.empty() || !output.empty() ||
            !runnerups.empty() || !candidate_manifest.empty() ||
            !pileup_sites.empty() || !pileup_observations.empty() ||
            !pileup_molecules.empty() || !joint_doublet_temp_dir.empty() ||
            !joint_doublet_manifest_digest.empty() ||
            !joint_doublet_workload_generation.empty();
        if (incompatible)
            die("targeted normalized-cache analysis rejects raw-input and legacy scoring options");
        if (joint_doublet_targeted_analysis_manifest.empty() ||
                joint_doublet_targeted_analysis_index<0 ||
                joint_doublet_normalized_cache_prefix.empty() ||
                joint_doublet_manifest.empty() ||
                joint_doublet_control_manifest.empty()){
            usage(); return 1;
        }
        try {
            if (setlocale(LC_NUMERIC,"C")==NULL)
                throw runtime_error("targeted analysis could not establish C numeric locale");
            run_joint_doublet_targeted_analysis(
                joint_doublet_targeted_analysis_manifest,
                joint_doublet_targeted_analysis_index,
                joint_doublet_normalized_cache_prefix,joint_doublet_manifest,
                joint_doublet_control_manifest,joint_doublet_targeted_analysis_output,
                joint_doublet_seed,joint_doublet_downsample_replicates,
                joint_doublet_null_replicates,
                strict_ld(axis_fmt(e_ref),"--error_ref"),
                strict_ld(axis_fmt(e_alt),"--error_alt"),min_evidence,
                max_second_fraction,joint_doublet_reference_mode,threads);
            return 0;
        } catch (const exception& error){
            fprintf(stderr,"ERROR: %s\n",error.what()); return 1;
        }
    }
    if (!joint_doublet_normalized_cache_prefix.empty()){
        const bool incompatible=!joint_doublet_output.empty() || !counts.empty() ||
            !assignments.empty() || !diagnostics.empty() || !output.empty() ||
            !runnerups.empty() || !candidate_manifest.empty() ||
            !site_fold_output.empty() || !probability_output.empty() ||
            !candidate_axis_output.empty() ||
            !joint_doublet_targeted_evidence_cache.empty() ||
            joint_doublet_reference_mode;
        if (incompatible)
            die("normalized cache extraction is standalone and rejects legacy scoring options");
        if (samples_path.empty() || joint_doublet_manifest.empty() ||
                joint_doublet_manifest_digest.empty() ||
                joint_doublet_workload_generation.empty() ||
                pileup_sites.empty() || pileup_observations.empty() ||
                pileup_molecules.empty() || joint_doublet_temp_dir.empty() ||
                libname.empty() || libname=="NA" || modality.empty()){
            usage(); return 1;
        }
        try {
            if (setlocale(LC_NUMERIC,"C")==NULL)
                throw runtime_error("normalized extraction could not establish C numeric locale");
            run_joint_doublet_normalized_extract(
                samples_path,joint_doublet_manifest,joint_doublet_manifest_digest,
                pileup_sites,pileup_observations,pileup_molecules,
                joint_doublet_normalized_cache_prefix,
                joint_doublet_workload_generation,joint_doublet_temp_dir,
                libname,modality,strict_ld(axis_fmt(e_ref),"--error_ref"),
                strict_ld(axis_fmt(e_alt),"--error_alt"),min_evidence,
                max_second_fraction,threads,
                joint_doublet_scheduler_memory_bytes,
                joint_doublet_launcher_runtime_reserve_bytes,
                joint_doublet_worker_memory_budget_bytes);
            return 0;
        } catch (const exception& error){
            fprintf(stderr,"ERROR: %s\n",error.what()); return 1;
        }
    }
    if (!joint_doublet_output.empty()){
        const bool incompatible = !counts.empty() || !assignments.empty() ||
            !diagnostics.empty() || !output.empty() || !runnerups.empty() ||
            !panel_path.empty() || !species_counts.empty() ||
            !species_condf.empty() || !species_samples_path.empty() ||
            !condf.empty() || !candidate_manifest.empty() ||
            !site_fold_output.empty() ||
            !probability_output.empty() || !score_prefix.empty() ||
            !candidate_axis_output.empty() || !candidate_axis_temp_dir.empty() ||
            explicit_probability_resamples || explicit_probability_seed ||
            explicit_site_folds || explicit_poor_fit;
        if (incompatible)
            die("--joint-doublet-output is standalone and rejects legacy/candidate-axis options");
        if (samples_path.empty() || joint_doublet_manifest.empty() ||
                pileup_sites.empty() || pileup_observations.empty() ||
                joint_doublet_temp_dir.empty() || libname.empty() ||
                libname == "NA" || modality.empty()){
            usage();
            return 1;
        }
        if (!joint_doublet_targeted_evidence_cache.empty() &&
                joint_doublet_targeted_evidence_cache==joint_doublet_output)
            die("targeted evidence cache and score output must be distinct files");
        try {
            if (setlocale(LC_NUMERIC,"C") == NULL)
                throw runtime_error(
                    "joint-doublet mode could not establish the C numeric locale");
            run_joint_doublet(samples_path,joint_doublet_manifest,pileup_sites,
                pileup_observations,pileup_molecules,joint_doublet_output,
                joint_doublet_targeted_evidence_cache,joint_doublet_temp_dir,
                libname,modality,
                strict_ld(axis_fmt(e_ref),"--error_ref"),
                strict_ld(axis_fmt(e_alt),"--error_alt"),min_evidence,
                strict_ld(axis_fmt(max_second_fraction),"--max-second-fraction"),
                joint_folds,threads);
            return 0;
        } catch (const exception& error){
            fprintf(stderr,"ERROR: %s\n",error.what());
            return 1;
        }
    }
    if (!candidate_axis_output.empty()){
        const bool incompatible = !counts.empty() || !assignments.empty() ||
            !diagnostics.empty() || !output.empty() || !runnerups.empty() ||
            !panel_path.empty() || !species_counts.empty() ||
            !species_condf.empty() || !species_samples_path.empty() ||
            !condf.empty() || !pileup_molecules.empty() ||
            !joint_doublet_targeted_evidence_cache.empty() ||
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
    if (!joint_doublet_targeted_evidence_cache.empty())
        die("--joint-doublet-targeted-evidence-cache requires --joint-doublet-output");
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
