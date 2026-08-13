#include <getopt.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <tuple>
#include <utility>
#include <vector>

#include <htswrapper/robin_hood/robin_hood.h>

#include "ambient_rna.h"
#include "common.h"
#include "demux_vcf_io.h"

using namespace std;

namespace {

const char* TOOL_NAME = "legacy2c_contam_estimate";
const char* TOOL_VERSION = "2.1";
const char* RUN_CONTRACT_VERSION = "legacy2c_contam_estimate_run_contract_V1_R3";

struct Options {
    string variant;
    string output_prefix;
    string receiver_lines;
    string ambient_candidates;
    string run_contract;
    string run_class = "production";
    string assignments_basis = "unspecified";
    string expected_lines_basis = "unspecified";
    string ambient_candidates_basis = "unspecified";
    double min_signal_gap = 0.10;
    int num_threads = 1;
    int mixprop_trials = 0;
    bool strict_condf = false;
};

struct CoverageSummary {
    size_t required_lookups = 0;
    size_t missing_lookups = 0;
    double required_weight = 0.0;
    double missing_weight = 0.0;
};

struct CandidateRoster {
    vector<string> lines;
    // combined_sources is the bounded row domain needed to load receiver-keyed
    // evidence. ambient_sources is the scientific mixture domain and must not
    // be widened merely because receiver identities are present in the bundle.
    set<int> combined_sources;
    set<int> ambient_sources;
};

string trim(const string& input) {
    const string whitespace = " \t\r\n";
    const size_t first = input.find_first_not_of(whitespace);
    if (first == string::npos) return "";
    const size_t last = input.find_last_not_of(whitespace);
    return input.substr(first, last - first + 1);
}

bool nonempty_file(const string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode) && st.st_size > 0;
}

long long file_size_or_minus_one(const string& path) {
    struct stat st;
    if (path.empty() || stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) return -1;
    return static_cast<long long>(st.st_size);
}

void require_nonempty(const string& path, const string& label) {
    if (!nonempty_file(path)) {
        throw runtime_error("required " + label + " missing or empty: " + path);
    }
}

string json_escape(const string& value) {
    ostringstream out;
    for (string::const_iterator it = value.begin(); it != value.end(); ++it) {
        const unsigned char c = static_cast<unsigned char>(*it);
        switch (c) {
            case '\\': out << "\\\\"; break;
            case '"': out << "\\\""; break;
            case '\b': out << "\\b"; break;
            case '\f': out << "\\f"; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (c < 0x20) {
                    out << "\\u" << hex << setw(4) << setfill('0') << static_cast<int>(c)
                        << dec << setfill(' ');
                } else {
                    out << static_cast<char>(c);
                }
        }
    }
    return out.str();
}

string json_string(const string& value) {
    return string("\"") + json_escape(value) + "\"";
}

string json_path_record(const string& path) {
    ostringstream out;
    out << "{\"path\":" << json_string(path)
        << ",\"size_bytes\":" << file_size_or_minus_one(path) << "}";
    return out.str();
}

void atomic_write_text(const string& path, const string& content) {
    ostringstream tmp_builder;
    tmp_builder << path << ".tmp." << static_cast<long long>(getpid());
    const string tmp = tmp_builder.str();
    {
        ofstream out(tmp.c_str(), ios::out | ios::trunc);
        if (!out) throw runtime_error("cannot open temporary output: " + tmp);
        out << content;
        out.flush();
        if (!out) throw runtime_error("failed writing temporary output: " + tmp);
    }
    if (rename(tmp.c_str(), path.c_str()) != 0) {
        const string why = strerror(errno);
        remove(tmp.c_str());
        throw runtime_error("cannot publish output " + path + ": " + why);
    }
}

vector<string> read_tokens(const string& path) {
    ifstream in(path.c_str());
    if (!in) throw runtime_error("cannot open roster file: " + path);
    vector<string> tokens;
    string line;
    while (getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        istringstream fields(line);
        string token;
        while (fields >> token) tokens.push_back(token);
    }
    return tokens;
}

vector<string> identity_components(const string& token) {
    vector<string> parts;
    size_t start = 0;
    while (true) {
        const size_t plus = token.find('+', start);
        const string part = trim(token.substr(start, plus == string::npos ? string::npos : plus - start));
        if (part.empty()) throw runtime_error("invalid empty identity component in: " + token);
        parts.push_back(part);
        if (plus == string::npos) break;
        start = plus + 1;
    }
    if (parts.size() > 2) {
        throw runtime_error("unsupported identity with more than two components: " + token);
    }
    return parts;
}

CandidateRoster build_candidate_roster(
    const string& receiver_path,
    const string& ambient_path,
    const vector<string>& samples,
    const string& output_path
) {
    const vector<string> receiver = read_tokens(receiver_path);
    const vector<string> ambient = read_tokens(ambient_path);
    if (receiver.empty()) throw runtime_error("receiver-lines file contains no identities: " + receiver_path);
    if (ambient.empty()) throw runtime_error("ambient-candidates file contains no identities: " + ambient_path);

    map<string, int> sample_to_index;
    for (size_t i = 0; i < samples.size(); ++i) {
        if (!sample_to_index.insert(make_pair(samples[i], static_cast<int>(i))).second) {
            throw runtime_error("samples file contains duplicate sample name: " + samples[i]);
        }
    }

    CandidateRoster result;
    set<string> seen;
    vector<string> all = receiver;
    all.insert(all.end(), ambient.begin(), ambient.end());
    for (vector<string>::const_iterator it = all.begin(); it != all.end(); ++it) {
        const vector<string> parts = identity_components(*it);
        for (vector<string>::const_iterator part = parts.begin(); part != parts.end(); ++part) {
            map<string, int>::const_iterator found = sample_to_index.find(*part);
            if (found == sample_to_index.end()) {
                throw runtime_error("identity " + *it + " references sample absent from .samples: " + *part);
            }
            result.combined_sources.insert(found->second);
        }
        if (seen.insert(*it).second) result.lines.push_back(*it);
    }

    for (vector<string>::const_iterator it = ambient.begin(); it != ambient.end(); ++it) {
        const vector<string> parts = identity_components(*it);
        for (vector<string>::const_iterator part = parts.begin(); part != parts.end(); ++part) {
            result.ambient_sources.insert(sample_to_index.at(*part));
        }
    }
    if (result.ambient_sources.empty()) {
        throw runtime_error("ambient-candidates file contains no usable source components: " + ambient_path);
    }

    ostringstream content;
    for (vector<string>::const_iterator it = result.lines.begin(); it != result.lines.end(); ++it) {
        content << *it << '\n';
    }
    atomic_write_text(output_path, content.str());
    return result;
}

void parse_tetraploid_receiver_roster(
    const string& expected_lines,
    const vector<string>& samples,
    set<int>& locked_identities,
    set<int>& safe_singlets
) {
    map<string, int> name_to_idx;
    for (size_t i = 0; i < samples.size(); ++i) name_to_idx[samples[i]] = static_cast<int>(i);

    const vector<string> entries = read_tokens(expected_lines);
    set<string> combo_individuals;
    set<string> bare_singlets;
    for (vector<string>::const_iterator it = entries.begin(); it != entries.end(); ++it) {
        const vector<string> parts = identity_components(*it);
        if (parts.size() == 2 && parts[0] != parts[1]) {
            combo_individuals.insert(parts[0]);
            combo_individuals.insert(parts[1]);
            const int a = name_to_idx.at(parts[0]);
            const int b = name_to_idx.at(parts[1]);
            locked_identities.insert(hap_comb_to_idx(min(a, b), max(a, b), samples.size()));
        } else {
            bare_singlets.insert(parts[0]);
        }
    }

    for (set<string>::const_iterator it = bare_singlets.begin(); it != bare_singlets.end(); ++it) {
        const int idx = name_to_idx.at(*it);
        if (combo_individuals.count(*it)) locked_identities.insert(idx);
        else safe_singlets.insert(idx);
    }
    for (set<string>::const_iterator it = combo_individuals.begin(); it != combo_individuals.end(); ++it) {
        if (!bare_singlets.count(*it)) locked_identities.insert(name_to_idx.at(*it));
    }
}

CoverageSummary audit_condf_coverage(
    const robin_hood::unordered_map<unsigned long,
        map<pair<int, int>, map<pair<int, int>, pair<float, float> > > >& counts,
    const map<pair<int, int>, map<int, float> >& condf,
    const vector<string>& samples,
    const set<int>& active_sources,
    const string& report_path
) {
    typedef tuple<int, int, int> LookupKey;
    map<LookupKey, size_t> rows;
    map<LookupKey, double> weights;

    for (robin_hood::unordered_map<unsigned long,
             map<pair<int, int>, map<pair<int, int>, pair<float, float> > > >::const_iterator cell = counts.begin();
         cell != counts.end(); ++cell) {
        for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::const_iterator first = cell->second.begin();
             first != cell->second.end(); ++first) {
            const int indv1 = first->first.first;
            const int type1 = first->first.second;
            for (map<pair<int, int>, pair<float, float> >::const_iterator second = first->second.begin();
                 second != first->second.end(); ++second) {
                const int indv2 = second->first.first;
                const int type2 = second->first.second;
                const double weight = static_cast<double>(second->second.first) + static_cast<double>(second->second.second);
                if (!isfinite(weight) || weight < 0.0) {
                    throw runtime_error("non-finite or negative allele-count weight encountered during .condf audit");
                }
                if (weight <= 0.0) continue;
                for (set<int>::const_iterator source = active_sources.begin(); source != active_sources.end(); ++source) {
                    vector<LookupKey> required;
                    if (indv2 == -1) {
                        required.push_back(make_tuple(indv1, type1, *source));
                    } else if (*source != indv1 && *source != indv2) {
                        required.push_back(make_tuple(indv1, type1, *source));
                        required.push_back(make_tuple(indv2, type2, *source));
                    }
                    for (vector<LookupKey>::const_iterator key = required.begin(); key != required.end(); ++key) {
                        rows[*key] += 1;
                        weights[*key] += weight;
                    }
                }
            }
        }
    }

    CoverageSummary summary;
    summary.required_lookups = rows.size();
    for (map<LookupKey, double>::const_iterator it = weights.begin(); it != weights.end(); ++it) {
        summary.required_weight += it->second;
    }

    ostringstream report;
    report << "condition_sample_index\tcondition_sample\tnalt\tsource_index\tsource_sample\t"
              "lookup_rows\tlookup_weight\tstatus\n";
    report << setprecision(12);
    for (map<LookupKey, size_t>::const_iterator it = rows.begin(); it != rows.end(); ++it) {
        const int condition_index = get<0>(it->first);
        const int nalt = get<1>(it->first);
        const int source_index = get<2>(it->first);
        bool present = false;
        map<pair<int, int>, map<int, float> >::const_iterator outer = condf.find(make_pair(condition_index, nalt));
        if (outer != condf.end()) present = outer->second.count(source_index) > 0;
        if (!present) {
            ++summary.missing_lookups;
            summary.missing_weight += weights[it->first];
        }
        report << condition_index << '\t' << samples.at(condition_index) << '\t' << nalt << '\t'
               << source_index << '\t' << samples.at(source_index) << '\t' << it->second << '\t'
               << weights[it->first] << '\t' << (present ? "present" : "missing") << '\n';
    }
    atomic_write_text(report_path, report.str());
    return summary;
}

void normalize_singlet_null_keys(
    robin_hood::unordered_map<unsigned long,
        map<pair<int, int>, map<pair<int, int>, pair<float, float> > > >& counts
) {
    const pair<int, int> canonical_null = make_pair(-1, -1);
    for (robin_hood::unordered_map<unsigned long,
             map<pair<int, int>, map<pair<int, int>, pair<float, float> > > >::iterator cell = counts.begin();
         cell != counts.end(); ++cell) {
        for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator first = cell->second.begin();
             first != cell->second.end(); ++first) {
            vector<pair<int, int> > aliases;
            pair<float, float> total(0.0f, 0.0f);
            bool found = false;
            for (map<pair<int, int>, pair<float, float> >::const_iterator second = first->second.begin();
                 second != first->second.end(); ++second) {
                if (second->first.first == -1) {
                    aliases.push_back(second->first);
                    total.first += second->second.first;
                    total.second += second->second.second;
                    found = true;
                }
            }
            if (!found) continue;
            for (vector<pair<int, int> >::const_iterator alias = aliases.begin(); alias != aliases.end(); ++alias) {
                first->second.erase(*alias);
            }
            first->second[canonical_null] = total;
        }
    }
}

size_t assignment_changes(
    const robin_hood::unordered_map<unsigned long, int>& before,
    const robin_hood::unordered_map<unsigned long, int>& after
) {
    size_t changed = 0;
    for (robin_hood::unordered_map<unsigned long, int>::const_iterator it = before.begin(); it != before.end(); ++it) {
        robin_hood::unordered_map<unsigned long, int>::const_iterator found = after.find(it->first);
        if (found == after.end() || found->second != it->second) ++changed;
    }
    for (robin_hood::unordered_map<unsigned long, int>::const_iterator it = after.begin(); it != after.end(); ++it) {
        if (before.count(it->first) == 0) ++changed;
    }
    return changed;
}

void validate_rates(const robin_hood::unordered_map<unsigned long, double>& rates) {
    if (rates.empty()) throw runtime_error("estimator produced no contamination-rate estimates");
    for (robin_hood::unordered_map<unsigned long, double>::const_iterator it = rates.begin(); it != rates.end(); ++it) {
        if (!isfinite(it->second) || it->second < 0.0 || it->second > 1.0) {
            throw runtime_error("estimator produced a non-finite or out-of-range contamination rate");
        }
    }
}

double validate_profile(const map<int, double>& profile) {
    if (profile.empty()) throw runtime_error("estimator produced an empty contamination profile");
    double total = 0.0;
    for (map<int, double>::const_iterator it = profile.begin(); it != profile.end(); ++it) {
        if (!isfinite(it->second) || it->second < 0.0) {
            throw runtime_error("estimator produced a non-finite or negative contamination-profile component");
        }
        total += it->second;
    }
    if (!isfinite(total) || total <= 0.0 || fabs(total - 1.0) > 1e-4) {
        ostringstream msg;
        msg << setprecision(12) << "contamination profile does not sum to one: " << total;
        throw runtime_error(msg.str());
    }
    return total;
}

void write_model_outputs(
    const string& prefix,
    contamFinder& finder,
    vector<string>& samples
) {
    map<int, double> empty_concentration;
    string blank_libname;
    const bool cellranger = false;
    const bool seurat = false;
    const bool underscore = false;

    const string profile_path = prefix + ".contam_prof";
    const string rate_path = prefix + ".contam_rate";
    const string assignments_path = prefix + ".decontam.assignments";
    const string token = ".tmp.legacy2c." + std::to_string((long long)getpid());
    const string profile_tmp = profile_path + token;
    const string rate_tmp = rate_path + token;
    const string assignments_tmp = assignments_path + token;
    vector<string> installed;

    auto cleanup_temporary = [&]() {
        remove(profile_tmp.c_str());
        remove(rate_tmp.c_str());
        remove(assignments_tmp.c_str());
    };
    auto rollback_installed = [&]() {
        for (vector<string>::const_iterator it = installed.begin(); it != installed.end(); ++it) {
            remove(it->c_str());
        }
    };
    auto publish_one = [&](const string& temporary, const string& final_path) {
        if (!nonempty_file(temporary)) {
            throw runtime_error("temporary model output missing or empty: " + temporary);
        }
        if (rename(temporary.c_str(), final_path.c_str()) != 0) {
            throw runtime_error("cannot publish model output " + final_path + ": " + strerror(errno));
        }
        installed.push_back(final_path);
    };

    try {
        FILE* profile = fopen(profile_tmp.c_str(), "w");
        if (!profile) throw runtime_error("cannot open temporary contamination-profile output: " + profile_tmp);
        dump_contam_prof(profile, finder.contam_prof, empty_concentration, samples);
        if (fclose(profile) != 0) throw runtime_error("failed closing temporary contamination-profile output: " + profile_tmp);

        FILE* rate = fopen(rate_tmp.c_str(), "w");
        if (!rate) throw runtime_error("cannot open temporary contamination-rate output: " + rate_tmp);
        dump_contam_rates(rate, finder.contam_rate, finder.contam_rate_se, samples,
                          blank_libname, cellranger, seurat, underscore);
        if (fclose(rate) != 0) throw runtime_error("failed closing temporary contamination-rate output: " + rate_tmp);

        FILE* assignments = fopen(assignments_tmp.c_str(), "w");
        if (!assignments) throw runtime_error("cannot open temporary assignment-audit output: " + assignments_tmp);
        dump_assignments(assignments, finder.assn, finder.assn_llr, samples,
                         blank_libname, cellranger, seurat, underscore);
        if (fclose(assignments) != 0) throw runtime_error("failed closing temporary assignment-audit output: " + assignments_tmp);

        publish_one(profile_tmp, profile_path);
        publish_one(rate_tmp, rate_path);
        publish_one(assignments_tmp, assignments_path);
    } catch (...) {
        cleanup_temporary();
        rollback_installed();
        throw;
    }
}

string double_json(double value) {
    if (!isfinite(value)) return "null";
    ostringstream out;
    out << setprecision(15) << value;
    return out.str();
}

void write_profile_diagnostics(const string& path, int mixprop_trials, size_t n_components, double ll) {
    const size_t configured = 2 + n_components * static_cast<size_t>(mixprop_trials);
    ostringstream out;
    out << "configured_starts\tsuccessful_starts\tnear_optimal_starts\tbest_log_likelihood\t"
           "second_best_log_likelihood\tnear_optimal_l1_spread\tprofile_nonunique_flag\tdiagnostic_status\n";
    out << configured << "\t" << configured << "\t1\t" << setprecision(15) << ll
        << "\tnan\tnan\t0\tlegacy2c_solver_exports_selected_fit_only\n";
    atomic_write_text(path, out.str());
}

void write_legacy_diagnostics(
    const string& path,
    const Options& opt,
    size_t assigned_cells,
    size_t changed,
    size_t estimated_cells,
    size_t roster_lines,
    size_t active_sources,
    size_t profile_components,
    double profile_sum,
    const CoverageSummary& coverage,
    double ll
) {
    ostringstream out;
    out << "tool\ttool_version\tvariant\timplementation\tassignments_frozen\tassigned_cells\t"
           "changed_assignments\testimated_cells\tcandidate_roster_lines\tactive_profile_sources\t"
           "profile_components\tprofile_sum\tmixprop_trials\tdeterministic_profile_starts\tstrict_condf\t"
           "condf_required_lookups\tcondf_missing_lookups\tcondf_required_weight\tcondf_missing_weight\t"
           "tetraploid_aware\tmin_signal_gap\tselected_log_likelihood\n";
    out << TOOL_NAME << '\t' << TOOL_VERSION << '\t' << opt.variant
        << "\tcompiled_current_contamFinder\t1\t" << assigned_cells << '\t' << changed << '\t'
        << estimated_cells << '\t' << roster_lines << '\t' << active_sources << '\t'
        << profile_components << '\t' << setprecision(12) << profile_sum << '\t'
        << opt.mixprop_trials << '\t' << (opt.mixprop_trials == 0 ? 1 : 0) << "\t1\t"
        << coverage.required_lookups << '\t' << coverage.missing_lookups << '\t'
        << coverage.required_weight << '\t' << coverage.missing_weight << '\t'
        << (opt.variant == "tet_aware" ? 1 : 0) << '\t';
    if (opt.variant == "tet_aware") out << opt.min_signal_gap;
    else out << "NA";
    out << '\t' << ll << '\n';
    atomic_write_text(path, out.str());
}

void write_run_contract(
    const string& path,
    const Options& opt,
    const string& samples_path,
    const string& assignments_path,
    const string& counts_path,
    const string& condf_path,
    const string& candidate_roster,
    const vector<string>& argv_values
) {
    ostringstream out;
    out << "{\n";
    out << "  \"contract_version\": " << json_string(RUN_CONTRACT_VERSION) << ",\n";
    out << "  \"tool\": " << json_string(TOOL_NAME) << ",\n";
    out << "  \"tool_version\": " << json_string(TOOL_VERSION) << ",\n";
    out << "  \"implementation\": \"compiled_current_contamFinder\",\n";
    out << "  \"variant\": " << json_string(opt.variant) << ",\n";
    out << "  \"run_class\": " << json_string(opt.run_class) << ",\n";
    out << "  \"production_contract_pass\": true,\n";
    out << "  \"production_contract_reason\": \"ok\",\n";
    out << "  \"panel_mode\": \"interindividual\",\n";
    out << "  \"output_prefix\": " << json_string(opt.output_prefix) << ",\n";
    out << "  \"strict_condf\": true,\n";
    out << "  \"assignments_basis\": " << json_string(opt.assignments_basis) << ",\n";
    out << "  \"expected_lines_basis\": " << json_string(opt.expected_lines_basis) << ",\n";
    out << "  \"ambient_candidates_basis\": " << json_string(opt.ambient_candidates_basis) << ",\n";
    out << "  \"warm_start_basis\": \"none\",\n";
    out << "  \"fixed_ambient_basis\": \"none\",\n";
    out << "  \"fixed_r_basis\": \"none\",\n";
    out << "  \"fix_r_basis\": \"none\",\n";
    out << "  \"assignments_frozen\": true,\n";
    out << "  \"tetraploid_aware\": " << (opt.variant == "tet_aware" ? "true" : "false") << ",\n";
    out << "  \"min_signal_gap\": "
        << (opt.variant == "tet_aware" ? double_json(opt.min_signal_gap) : "null") << ",\n";
    out << "  \"mixprop_trials\": " << opt.mixprop_trials << ",\n";
    out << "  \"counts\": " << json_path_record(counts_path) << ",\n";
    out << "  \"condf\": " << json_path_record(condf_path) << ",\n";
    out << "  \"assignments\": " << json_path_record(assignments_path) << ",\n";
    out << "  \"samples\": " << json_path_record(samples_path) << ",\n";
    out << "  \"expected_lines\": " << json_path_record(opt.receiver_lines) << ",\n";
    out << "  \"ambient_candidates\": " << json_path_record(opt.ambient_candidates) << ",\n";
    out << "  \"legacy2c_candidate_roster\": " << json_path_record(candidate_roster) << ",\n";
    out << "  \"warm_start\": {\"path\":\"\",\"size_bytes\":-1},\n";
    out << "  \"fixed_ambient\": {\"path\":\"\",\"size_bytes\":-1},\n";
    out << "  \"fix_r\": {\"path\":\"\",\"size_bytes\":-1},\n";
    out << "  \"command\": [";
    for (size_t i = 0; i < argv_values.size(); ++i) {
        if (i) out << ',';
        out << json_string(argv_values[i]);
    }
    out << "]\n";
    out << "}\n";
    atomic_write_text(path, out.str());
}

void print_help(FILE* stream) {
    fprintf(stream,
        "legacy2c_contam_estimate [OPTIONS]\n"
        "Compiled two-component compatibility estimator using the current hardened\n"
        "CellBouncer contamFinder implementation. It never invokes quant_contam.\n\n"
        "Required:\n"
        "  --variant classic|tet_aware\n"
        "  --output_prefix, -o PREFIX\n"
        "  --receiver_lines, -X FILE\n"
        "  --ambient_candidates FILE\n"
        "  --strict_condf\n"
        "  --run_class production\n"
        "  --assignments_basis library\n"
        "  --expected_lines_basis library\n"
        "  --ambient_candidates_basis library\n\n"
        "Optional:\n"
        "  --min_signal_gap FLOAT   Tet-aware signal gap (default 0.10)\n"
        "  --num_threads, -T INT    Optimization worker count (default 1)\n"
        "  --mixprop_trials INT     Random starts per profile component (default 0)\n"
        "  --run_contract FILE      Output run-contract path\n"
        "  --version\n");
}

Options parse_options(int argc, char** argv) {
    enum {
        OPT_AMBIENT_CANDIDATES = 1000,
        OPT_MIN_SIGNAL_GAP,
        OPT_MIXPROP_TRIALS,
        OPT_RUN_CLASS,
        OPT_STRICT_CONDF,
        OPT_RUN_CONTRACT,
        OPT_ASSIGNMENTS_BASIS,
        OPT_EXPECTED_LINES_BASIS,
        OPT_AMBIENT_CANDIDATES_BASIS,
        OPT_VERSION
    };
    static struct option long_options[] = {
        {"variant", required_argument, 0, 'v'},
        {"output_prefix", required_argument, 0, 'o'},
        {"receiver_lines", required_argument, 0, 'X'},
        {"ambient_candidates", required_argument, 0, OPT_AMBIENT_CANDIDATES},
        {"min_signal_gap", required_argument, 0, OPT_MIN_SIGNAL_GAP},
        {"num_threads", required_argument, 0, 'T'},
        {"mixprop_trials", required_argument, 0, OPT_MIXPROP_TRIALS},
        {"run_class", required_argument, 0, OPT_RUN_CLASS},
        {"strict_condf", no_argument, 0, OPT_STRICT_CONDF},
        {"run_contract", required_argument, 0, OPT_RUN_CONTRACT},
        {"assignments_basis", required_argument, 0, OPT_ASSIGNMENTS_BASIS},
        {"expected_lines_basis", required_argument, 0, OPT_EXPECTED_LINES_BASIS},
        {"ambient_candidates_basis", required_argument, 0, OPT_AMBIENT_CANDIDATES_BASIS},
        {"version", no_argument, 0, OPT_VERSION},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    Options opt;
    int option_index = 0;
    int ch;
    while ((ch = getopt_long(argc, argv, "v:o:X:T:h", long_options, &option_index)) != -1) {
        switch (ch) {
            case 'v': opt.variant = optarg; break;
            case 'o': opt.output_prefix = optarg; break;
            case 'X': opt.receiver_lines = optarg; break;
            case 'T': opt.num_threads = atoi(optarg); break;
            case OPT_AMBIENT_CANDIDATES: opt.ambient_candidates = optarg; break;
            case OPT_MIN_SIGNAL_GAP: opt.min_signal_gap = atof(optarg); break;
            case OPT_MIXPROP_TRIALS: opt.mixprop_trials = atoi(optarg); break;
            case OPT_RUN_CLASS: opt.run_class = optarg; break;
            case OPT_STRICT_CONDF: opt.strict_condf = true; break;
            case OPT_RUN_CONTRACT: opt.run_contract = optarg; break;
            case OPT_ASSIGNMENTS_BASIS: opt.assignments_basis = optarg; break;
            case OPT_EXPECTED_LINES_BASIS: opt.expected_lines_basis = optarg; break;
            case OPT_AMBIENT_CANDIDATES_BASIS: opt.ambient_candidates_basis = optarg; break;
            case OPT_VERSION:
                cout << TOOL_NAME << ' ' << TOOL_VERSION << '\n';
                exit(0);
            case 'h':
                print_help(stdout);
                exit(0);
            default:
                print_help(stderr);
                throw runtime_error("invalid command-line arguments");
        }
    }

    if (opt.variant != "classic" && opt.variant != "tet_aware") {
        throw runtime_error("--variant must be classic or tet_aware");
    }
    if (opt.output_prefix.empty()) throw runtime_error("--output_prefix is required");
    if (opt.receiver_lines.empty()) throw runtime_error("--receiver_lines is required");
    if (opt.ambient_candidates.empty()) throw runtime_error("--ambient_candidates is required");
    if (opt.num_threads < 1) throw runtime_error("--num_threads must be >= 1");
    if (opt.mixprop_trials < 0) throw runtime_error("--mixprop_trials must be >= 0");
    if (!isfinite(opt.min_signal_gap) || opt.min_signal_gap < 0.0 || opt.min_signal_gap > 1.0) {
        throw runtime_error("--min_signal_gap must be between 0 and 1");
    }
    if (!opt.strict_condf) throw runtime_error("--strict_condf is mandatory");
    if (opt.run_class != "production" || opt.assignments_basis != "library" ||
        opt.expected_lines_basis != "library" || opt.ambient_candidates_basis != "library") {
        throw runtime_error(
            "legacy2c benchmark runs require run_class=production and library basis for "
            "assignments, expected lines, and ambient candidates"
        );
    }
    if (opt.run_contract.empty()) opt.run_contract = opt.output_prefix + ".run_contract.json";
    return opt;
}

int run(int argc, char** argv) {
    const Options opt = parse_options(argc, argv);
    vector<string> argv_values;
    for (int i = 0; i < argc; ++i) argv_values.push_back(argv[i]);

    const string samples_path = opt.output_prefix + ".samples";
    const string assignments_path = opt.output_prefix + ".assignments";
    const string counts_path = opt.output_prefix + ".counts";
    const string condf_path = opt.output_prefix + ".condf";
    const string candidate_roster_path = opt.output_prefix + ".legacy2c_candidates.txt";
    const string coverage_path = opt.output_prefix + ".condf_coverage.tsv";
    const string profile_diag_path = opt.output_prefix + ".profile_fit_diagnostics.tsv";
    const string legacy_diag_path = opt.output_prefix + ".legacy2c_diagnostics.tsv";
    const string rate_path = opt.output_prefix + ".contam_rate";
    const string profile_path = opt.output_prefix + ".contam_prof";
    const string decontam_path = opt.output_prefix + ".decontam.assignments";

    require_nonempty(samples_path, "samples input");
    require_nonempty(assignments_path, "assignments input");
    require_nonempty(counts_path, "counts input");
    require_nonempty(condf_path, "conditional-frequency input");
    require_nonempty(opt.receiver_lines, "receiver-lines input");
    require_nonempty(opt.ambient_candidates, "ambient-candidates input");
    if (nonempty_file(rate_path) || nonempty_file(profile_path) || nonempty_file(decontam_path)) {
        throw runtime_error("refusing to reuse or overwrite an existing estimator output for prefix: " + opt.output_prefix);
    }

    vector<string> samples;
    string samples_path_mut = samples_path;
    load_samples(samples_path_mut, samples);
    if (samples.empty()) throw runtime_error("samples input contains no sample names: " + samples_path);

    const CandidateRoster roster = build_candidate_roster(
        opt.receiver_lines, opt.ambient_candidates, samples, candidate_roster_path
    );

    set<int> allowed_ids;
    set<int> allowed_ids2;
    string roster_mut = candidate_roster_path;
    parse_idfile(roster_mut, samples, allowed_ids, allowed_ids2, false);
    if (allowed_ids.empty()) throw runtime_error("candidate roster contains no valid identities");

    map<pair<int, int>, map<int, float> > exp_match_fracs;
    string condf_path_mut = condf_path;
    load_exp_fracs(condf_path_mut, exp_match_fracs);
    if (exp_match_fracs.empty()) throw runtime_error("conditional-frequency input contains no entries: " + condf_path);

    robin_hood::unordered_map<unsigned long,
        map<pair<int, int>, map<pair<int, int>, pair<float, float> > > > allele_counts;
    string counts_path_mut = counts_path;
    load_counts_from_file(allele_counts, samples, counts_path_mut, allowed_ids);
    if (allele_counts.empty()) throw runtime_error("count input contains no rows accepted by the candidate roster");
    // Current count bundles encode the unused singlet second-component nalt as
    // zero, while the retained two-component contamFinder uses (-1,-1) as its
    // null sentinel.  Canonicalize only that representation; no source or row
    // selection is changed.
    normalize_singlet_null_keys(allele_counts);

    const CoverageSummary coverage = audit_condf_coverage(
        allele_counts, exp_match_fracs, samples, roster.ambient_sources, coverage_path
    );
    if (coverage.missing_lookups > 0) {
        ostringstream message;
        message << setprecision(12)
                << "strict conditional-frequency coverage failed: " << coverage.missing_lookups
                << " missing lookup(s), weighted observations=" << coverage.missing_weight
                << "; see " << coverage_path;
        throw runtime_error(message.str());
    }

    robin_hood::unordered_map<unsigned long, int> assignments;
    robin_hood::unordered_map<unsigned long, double> assignment_llr;
    string assignments_path_mut = assignments_path;
    load_assignments_from_file(assignments_path_mut, assignments, assignment_llr, samples);
    if (assignments.empty()) throw runtime_error("assignments input contains no assignments");
    const robin_hood::unordered_map<unsigned long, int> original_assignments = assignments;

    set<int> profile_ids = roster.ambient_sources;
    contamFinder finder(
        allele_counts, assignments, assignment_llr, exp_match_fracs,
        static_cast<int>(samples.size()), profile_ids, allowed_ids2
    );
    map<int, double> initial_profile;
    const double initial_mass = 1.0 / static_cast<double>(profile_ids.size());
    for (set<int>::const_iterator it = profile_ids.begin(); it != profile_ids.end(); ++it) {
        initial_profile[*it] = initial_mass;
    }
    finder.set_init_contam_prof(initial_profile);
    finder.set_doublet_rate(-1.0);
    finder.set_num_threads(opt.num_threads <= 1 ? 0 : opt.num_threads);
    finder.no_reassign();
    finder.set_error_rates(0.001, 0.001);
    finder.set_mixprop_trials(opt.mixprop_trials);
    finder.use_weights();

    set<int> locked_identities;
    set<int> safe_singlets;
    if (opt.variant == "tet_aware") {
        parse_tetraploid_receiver_roster(
            opt.receiver_lines, samples, locked_identities, safe_singlets
        );
        finder.set_tetraploid_aware(true);
        finder.set_locked_identities(locked_identities);
        finder.set_safe_singlets(safe_singlets);
        finder.set_min_signal_gap(opt.min_signal_gap);
        finder.set_ids_restricted(true);
    }

    cerr << TOOL_NAME << " v" << TOOL_VERSION << ": variant=" << opt.variant
         << ", implementation=compiled_current_contamFinder" << endl;
    finder.fit();
    const double selected_ll = finder.compute_ll();

    const size_t changed = assignment_changes(original_assignments, finder.assn);
    if (changed != 0) {
        ostringstream message;
        message << "assignment-freeze contract violated: " << changed
                << " barcode assignment(s) changed or disappeared";
        throw runtime_error(message.str());
    }
    validate_rates(finder.contam_rate);
    const double profile_sum = validate_profile(finder.contam_prof);
    write_model_outputs(opt.output_prefix, finder, samples);

    require_nonempty(rate_path, "contamination-rate output");
    require_nonempty(profile_path, "contamination-profile output");
    require_nonempty(decontam_path, "assignment-audit output");

    write_profile_diagnostics(
        profile_diag_path, opt.mixprop_trials, finder.contam_prof.size(), selected_ll
    );
    write_legacy_diagnostics(
        legacy_diag_path, opt, original_assignments.size(), changed,
        finder.contam_rate.size(), roster.lines.size(), roster.ambient_sources.size(),
        finder.contam_prof.size(), profile_sum, coverage, selected_ll
    );
    write_run_contract(
        opt.run_contract, opt, samples_path, assignments_path, counts_path, condf_path,
        candidate_roster_path, argv_values
    );
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(argc, argv);
    } catch (const exception& exc) {
        cerr << "ERROR: " << exc.what() << endl;
        return 1;
    }
}
