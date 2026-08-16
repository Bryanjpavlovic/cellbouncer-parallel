// =============================================================================
// genotype_llr.cpp
// Unified genotype likelihood, assignment, and comparison implementation.
// =============================================================================

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
#include <math.h>
#include <limits>
#include <omp.h>
#include <mixtureDist/functions.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "genotype_llr.h"

using std::cout;
using std::endl;
using namespace std;

// Global verbose flag (set from main)
bool g_verbose = false;

// Global debug flag (set from main)
bool g_debug = false;

/**
 * ===== Contains functions relating to identifying individuals of  ===== 
 * ===== origin, used by demux_vcf.                                 =====
 */

// ============================================================================
// LLR TABLE CLASS IMPLEMENTATION
// ============================================================================

const char* comparison_state_name(ComparisonState state) {
    switch (state) {
        case ComparisonState::PRESENT_NONZERO: return "present_nonzero";
        case ComparisonState::PRESENT_ZERO: return "present_zero";
        case ComparisonState::UNAVAILABLE: return "unavailable";
        case ComparisonState::PARTIAL_SUPPORT: return "partial_support";
        case ComparisonState::NOT_APPLICABLE: return "not_applicable";
    }
    return "unavailable";
}

llr_table::llr_table(int x) : n_indvs(0) {
    const int n_elt = std::max(0, x);
    included.reserve(n_elt);
    maxllr.reserve(n_elt);
    minllr.reserve(n_elt);
}

llr_table::~llr_table() {
    lookup_llr.clear();
    pairwise_llr.clear();
    pairwise_partial_support.clear();
    included.clear();
    maxllr.clear();
    minllr.clear();
}

void llr_table::print(string& bc_str, vector<string>& samples) {
    for (auto x = lookup_llr.begin(); x != lookup_llr.end(); ++x) {
        for (size_t i = 0; i < x->second.size(); ++i) {
            const int low = x->second[i].first;
            const int high = x->second[i].second;
            if (low < 0 || high < 0 || low >= (int)included.size() ||
                high >= (int)included.size() || !included[low] || !included[high]) {
                continue;
            }
            const string n1 = idx2name(low, samples);
            const string n2 = idx2name(high, samples);
            fprintf(stdout, "%s\t%s\t%s\t%f\n", bc_str.c_str(), n1.c_str(), n2.c_str(), x->first);
            fprintf(stdout, "%s\t%s\t%s\t%f\n", bc_str.c_str(), n2.c_str(), n1.c_str(), -x->first);
        }
    }
}

void llr_table::print_ranges(string& barcode, vector<string>& samples) {
    for (size_t i = 0; i < included.size(); ++i) {
        if (!included[i]) continue;
        fprintf(stdout, "%s\t%s\t%f\t%f\n", barcode.c_str(),
                idx2name((int)i, samples).c_str(), minllr[i], maxllr[i]);
    }
}

void llr_table::insert(int i1, int i2, double llr) {
    if (i1 < 0 || i2 < 0 || i1 == i2 || !std::isfinite(llr)) {
        return;
    }

    const int max_index = std::max(i1, i2);
    while ((int)included.size() <= max_index) {
        included.push_back(false);
        maxllr.push_back(0.0);
        minllr.push_back(0.0);
    }

    if (!included[i1]) ++n_indvs;
    if (!included[i2]) ++n_indvs;

    if (!included[i1] || maxllr[i1] < llr) maxllr[i1] = llr;
    if (!included[i1] || minllr[i1] > llr) minllr[i1] = llr;
    if (!included[i2] || maxllr[i2] < -llr) maxllr[i2] = -llr;
    if (!included[i2] || minllr[i2] > -llr) minllr[i2] = -llr;
    included[i1] = true;
    included[i2] = true;

    const std::pair<int, int> key = (i1 < i2) ? std::make_pair(i1, i2)
                                               : std::make_pair(i2, i1);
    pairwise_partial_support.erase(key);
    pairwise_llr[key] = (i1 < i2) ? llr : -llr;

    int low = i1;
    int high = i2;
    double stored = llr;
    if (llr > 0.0) {
        low = i2;
        high = i1;
        stored = -llr;
    }
    lookup_llr[stored].push_back(std::make_pair(low, high));
}

void llr_table::mark_partial_support(int i1, int i2) {
    if (i1 < 0 || i2 < 0 || i1 == i2) return;
    const std::pair<int, int> key = (i1 < i2) ? std::make_pair(i1, i2)
                                               : std::make_pair(i2, i1);
    if (pairwise_llr.count(key) == 0) pairwise_partial_support.insert(key);
}

void llr_table::disallow(int i) {
    if (i < 0 || i >= (int)included.size()) return;
    if (included[i]) {
        maxllr[i] = 0.0;
        minllr[i] = 0.0;
        --n_indvs;
    }
    included[i] = false;
}

void llr_table::recalculate_minmax() {
    const double inf = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < included.size(); ++i) {
        if (included[i]) {
            minllr[i] = inf;
            maxllr[i] = -inf;
        }
    }

    for (const auto& entry : pairwise_llr) {
        const int first = entry.first.first;
        const int second = entry.first.second;
        if (first < 0 || second < 0 || first >= (int)included.size() ||
            second >= (int)included.size() || !included[first] || !included[second]) {
            continue;
        }
        const double first_vs_second = entry.second;
        minllr[first] = std::min(minllr[first], first_vs_second);
        maxllr[first] = std::max(maxllr[first], first_vs_second);
        minllr[second] = std::min(minllr[second], -first_vs_second);
        maxllr[second] = std::max(maxllr[second], -first_vs_second);
    }

    for (size_t i = 0; i < included.size(); ++i) {
        if (included[i] && minllr[i] == inf) {
            minllr[i] = -inf;
            maxllr[i] = -inf;
        }
    }
}

bool llr_table::del(int n_keep) {
    if (n_indvs < n_keep) return false;

    it = lookup_llr.begin();
    while (n_indvs > n_keep && it != lookup_llr.end()) {
        bool del_all = it->second.size() <= static_cast<size_t>(n_indvs - n_keep);
        map<double, int> maxllr_counts;
        int indv_tot = 0;
        for (auto x = it->second.begin(); x != it->second.end();) {
            if (x->first < 0 || x->second < 0 ||
                x->first >= (int)included.size() || x->second >= (int)included.size() ||
                !included[x->first] || !included[x->second]) {
                x = it->second.erase(x);
            } else if (del_all) {
                if (included[x->first]) {
                    --n_indvs;
                    included[x->first] = false;
                    minllr[x->first] = 0.0;
                    maxllr[x->first] = 0.0;
                }
                x = it->second.erase(x);
            } else {
                ++maxllr_counts[maxllr[x->first]];
                ++indv_tot;
                ++x;
            }
        }

        if (indv_tot <= n_indvs - n_keep) del_all = true;

        if (!del_all) {
            double cutoff = 0.0;
            int runtot = 0;
            for (const auto& m : maxllr_counts) {
                runtot += m.second;
                cutoff = m.first;
                if (runtot >= n_indvs - n_keep) break;
            }
            for (auto x = it->second.begin(); x != it->second.end();) {
                if (maxllr[x->first] <= cutoff) {
                    if (included[x->first]) {
                        included[x->first] = false;
                        maxllr[x->first] = 0.0;
                        minllr[x->first] = 0.0;
                        --n_indvs;
                    }
                    x = it->second.erase(x);
                } else {
                    ++x;
                }
            }
        } else {
            for (auto x = it->second.begin(); x != it->second.end();) {
                if (included[x->first]) {
                    included[x->first] = false;
                    --n_indvs;
                    maxllr[x->first] = 0.0;
                    minllr[x->first] = 0.0;
                }
                x = it->second.erase(x);
            }
            lookup_llr.erase(it++);
        }
    }
    return true;
}

void llr_table::get_max(int& best_idx, double& best_llr) const {
    best_idx = -1;
    best_llr = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < included.size(); ++i) {
        if (included[i] && (best_idx == -1 || minllr[i] > best_llr)) {
            best_idx = static_cast<int>(i);
            best_llr = minllr[i];
        }
    }
}

double llr_table::get_min_margin(int identity) const {
    if (identity < 0 || identity >= (int)minllr.size() || !included[identity]) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return minllr[identity];
}

void llr_table::get_max_by_max_llr_comparator(int& best_idx, double& best_maxllr) const {
    best_idx = -1;
    best_maxllr = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < included.size(); ++i) {
        if (included[i] && (best_idx == -1 || maxllr[i] > best_maxllr)) {
            best_idx = static_cast<int>(i);
            best_maxllr = maxllr[i];
        }
    }
}

PairwiseComparison llr_table::get_pairwise(int lhs, int rhs) const {
    if (lhs < 0 || rhs < 0 || lhs == rhs) return PairwiseComparison();
    const std::pair<int, int> key = (lhs < rhs) ? std::make_pair(lhs, rhs)
                                                : std::make_pair(rhs, lhs);
    const auto found = pairwise_llr.find(key);
    if (found != pairwise_llr.end()) {
        return PairwiseComparison(true, lhs < rhs ? found->second : -found->second);
    }
    if (pairwise_partial_support.count(key) > 0) {
        return PairwiseComparison(false,
            std::numeric_limits<double>::quiet_NaN(), true);
    }
    return PairwiseComparison();
}

vector<int> llr_table::retained_identities() const {
    vector<int> result;
    result.reserve(n_indvs > 0 ? static_cast<size_t>(n_indvs) : 0U);
    for (size_t i = 0; i < included.size(); ++i) {
        if (included[i]) result.push_back(static_cast<int>(i));
    }
    return result;
}

bool llr_table::winner_has_complete_comparisons(
    int winner,
    vector<int>& missing_alternatives) const {
    missing_alternatives.clear();
    if (winner < 0 || winner >= (int)included.size() || !included[winner]) return false;
    for (size_t i = 0; i < included.size(); ++i) {
        if (!included[i] || static_cast<int>(i) == winner) continue;
        if (!get_pairwise(winner, static_cast<int>(i)).present) {
            missing_alternatives.push_back(static_cast<int>(i));
        }
    }
    return missing_alternatives.empty();
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

double adjust_p_err(double p, double e_r, double e_a){
    return p - p*e_a + (1-p)*e_r;
}

// ============================================================================
// DOUBLET COMPARISON COMPUTATION (from original demux_vcf_llr.cpp)
// ============================================================================

namespace {

bool get_directed_llr(
    const map<int, map<int, double> >& llrs,
    int lhs,
    int rhs,
    double& value) {

    if (lhs == rhs) {
        value = 0.0;
        return true;
    }
    const int low = std::min(lhs, rhs);
    const int high = std::max(lhs, rhs);
    const auto low_it = llrs.find(low);
    if (low_it == llrs.end()) return false;
    const auto high_it = low_it->second.find(high);
    if (high_it == low_it->second.end() || !std::isfinite(high_it->second)) return false;
    value = lhs < rhs ? high_it->second : -high_it->second;
    return true;
}

bool get_conditional_llr(
    const map<int, map<int, map<int, double> > >& conditional,
    int anchor,
    int lhs,
    int rhs,
    double& value) {

    if (lhs == rhs) {
        value = 0.0;
        return true;
    }
    const auto anchor_it = conditional.find(anchor);
    if (anchor_it == conditional.end()) return false;
    const int low = std::min(lhs, rhs);
    const int high = std::max(lhs, rhs);
    const auto low_it = anchor_it->second.find(low);
    if (low_it == anchor_it->second.end()) return false;
    const auto high_it = low_it->second.find(high);
    if (high_it == low_it->second.end() || !std::isfinite(high_it->second)) return false;
    value = lhs < rhs ? high_it->second : -high_it->second;
    return true;
}

}  // namespace

void get_kcomps_cond(
    map<int, map<int, double> >& kcomps,
    const map<int, map<int, double> >& llrs,
    int n_samples) {

    for (const auto& llr : llrs) {
        vector<int> component_ids;
        vector<double> component_scores;
        for (const auto& llr2 : llr.second) {
            if (llr2.first < n_samples || !std::isfinite(llr2.second)) continue;
            const pair<int, int> comb = idx_to_hap_comb(llr2.first, n_samples);
            const int other = comb.first == llr.first ? comb.second : comb.first;
            component_ids.push_back(other);
            component_scores.push_back(-llr2.second);
        }

        for (size_t i = 0; i + 1 < component_scores.size(); ++i) {
            for (size_t j = i + 1; j < component_scores.size(); ++j) {
                const int first = component_ids[i];
                const int second = component_ids[j];
                if (first == second) continue;
                const int low = std::min(first, second);
                const int high = std::max(first, second);
                const double directed = component_scores[i] - component_scores[j];
                kcomps[low][high] += first < second ? directed : -directed;
            }
        }
    }
}

void compute_k_comps(
    const map<int, map<int, double> >& llrs,
    llr_table& tab,
    const vector<int>& ks,
    int n_samples,
    set<int>& allowed_assignments,
    double doublet_rate,
    map<int, double>* prior_weights) {

    // Conditional comparisons are stored only when their supporting direct
    // comparisons exist. Exact zero remains present; absent paths remain absent.
    map<int, map<int, map<int, double> > > kcomp_cond;
    for (const auto& sample : llrs) {
        if (sample.first < 0 || sample.first >= (int)tab.included.size() ||
            !tab.included[sample.first]) {
            continue;
        }
        map<int, map<int, double> > conditional_for_anchor;
        map<int, map<int, double> > one_anchor;
        one_anchor.emplace(sample.first, sample.second);
        get_kcomps_cond(conditional_for_anchor, one_anchor, n_samples);
        if (!conditional_for_anchor.empty()) {
            kcomp_cond.emplace(sample.first, std::move(conditional_for_anchor));
        }
    }

    // Compute singlet-versus-doublet comparisons. The established complete-data
    // equation is retained exactly: a related singlet uses its one direct path;
    // an unrelated singlet uses the mean of both defined component paths. A
    // proper nonempty subset is recorded as partial support and is never inserted
    // as a production comparison.
    for (int i = 0; i < n_samples; ++i) {
        if ((!allowed_assignments.empty() && allowed_assignments.count(i) == 0) ||
            i >= (int)tab.included.size() || !tab.included[i]) {
            continue;
        }

        for (int k : ks) {
            const pair<int, int> comb = idx_to_hap_comb(k, n_samples);
            if (comb.first == i || comb.second == i) {
                double direct = 0.0;
                if (get_directed_llr(llrs, i, k, direct)) {
                    if (prior_weights != NULL &&
                        prior_weights->count(i) > 0 && prior_weights->count(k) > 0) {
                        direct += (*prior_weights)[i] - (*prior_weights)[k];
                    } else if (doublet_rate != 0.5 && doublet_rate < 1.0) {
                        direct += log2(1.0 - doublet_rate) - log2(doublet_rate);
                    }
                    tab.insert(i, k, direct);
                }
                continue;
            }

            vector<double> complete_paths;
            size_t atomic_support = 0;
            const int components[2] = {comb.first, comb.second};
            for (int component : components) {
                double component_vs_doublet = 0.0;
                double singlet_vs_component = 0.0;
                const bool have_component_vs_doublet =
                    get_directed_llr(llrs, component, k, component_vs_doublet);
                const bool have_singlet_vs_component =
                    get_directed_llr(llrs, i, component, singlet_vs_component);
                atomic_support += have_component_vs_doublet ? 1U : 0U;
                atomic_support += have_singlet_vs_component ? 1U : 0U;
                if (have_component_vs_doublet && have_singlet_vs_component) {
                    complete_paths.push_back(
                        component_vs_doublet + singlet_vs_component);
                }
            }

            if (complete_paths.size() != 2U) {
                if (atomic_support > 0U) tab.mark_partial_support(i, k);
                continue;
            }

            double llr = (complete_paths[0] + complete_paths[1]) / 2.0;
            if (prior_weights != NULL &&
                prior_weights->count(i) > 0 && prior_weights->count(k) > 0) {
                llr += (*prior_weights)[i] - (*prior_weights)[k];
            } else if (doublet_rate != 0.5 && doublet_rate < 1.0) {
                llr += log2(1.0 - doublet_rate) - log2(doublet_rate);
            }
            tab.insert(i, k, llr);
        }
    }

    // Compute doublet-versus-doublet comparisons. The existing complete-data
    // equation uses one conditional path for pairs sharing a component and all
    // eight two-anchor support paths for disjoint pairs. Every required path must
    // be present. A proper subset is represented explicitly as partial support.
    for (size_t ki = 0; ki + 1 < ks.size(); ++ki) {
        const int k1 = ks[ki];
        const pair<int, int> comb1 = idx_to_hap_comb(k1, n_samples);
        for (size_t kj = ki + 1; kj < ks.size(); ++kj) {
            const int k2 = ks[kj];
            const pair<int, int> comb2 = idx_to_hap_comb(k2, n_samples);
            const int a = comb1.first;
            const int b = comb1.second;
            const int c = comb2.first;
            const int d = comb2.second;

            int common = -1;
            int other1 = -1;
            int other2 = -1;
            if (a == c) { common = a; other1 = b; other2 = d; }
            else if (a == d) { common = a; other1 = b; other2 = c; }
            else if (b == c) { common = b; other1 = a; other2 = d; }
            else if (b == d) { common = b; other1 = a; other2 = c; }

            if (common >= 0) {
                double part = 0.0;
                if (get_conditional_llr(kcomp_cond, common, other1, other2, part)) {
                    if (prior_weights != NULL &&
                        prior_weights->count(k1) > 0 && prior_weights->count(k2) > 0) {
                        part += (*prior_weights)[k1] - (*prior_weights)[k2];
                    }
                    tab.insert(k1, k2, part);
                }
                continue;
            }

            vector<double> llr_parts;
            size_t atomic_support = 0;
            auto append_two_path = [&](int anchor1, int lhs1, int rhs1,
                                       int anchor2, int lhs2, int rhs2,
                                       double multiplier) {
                double first = 0.0;
                double second = 0.0;
                const bool have_first =
                    get_conditional_llr(kcomp_cond, anchor1, lhs1, rhs1, first);
                const bool have_second =
                    get_conditional_llr(kcomp_cond, anchor2, lhs2, rhs2, second);
                atomic_support += have_first ? 1U : 0U;
                atomic_support += have_second ? 1U : 0U;
                if (have_first && have_second) {
                    llr_parts.push_back(multiplier * (first + second));
                }
            };
            append_two_path(a, b, c, c, a, d, 1.0);
            append_two_path(b, a, c, c, b, d, 1.0);
            append_two_path(a, b, d, d, a, c, 1.0);
            append_two_path(b, a, d, d, b, c, 1.0);
            append_two_path(c, d, a, a, c, b, -1.0);
            append_two_path(d, c, a, a, d, b, -1.0);
            append_two_path(c, d, b, b, c, a, -1.0);
            append_two_path(d, c, b, b, d, a, -1.0);

            if (llr_parts.size() != 8U) {
                if (atomic_support > 0U) tab.mark_partial_support(k1, k2);
                continue;
            }

            double llr = 0.0;
            for (double part : llr_parts) llr += part;
            llr /= 8.0;
            if (prior_weights != NULL &&
                prior_weights->count(k1) > 0 && prior_weights->count(k2) > 0) {
                llr += (*prior_weights)[k1] - (*prior_weights)[k2];
            }
            tab.insert(k1, k2, llr);
        }
    }
}

// ============================================================================
// ORIGINAL LLR POPULATION (for legacy nested map structure)
// ============================================================================

bool populate_llr_table(map<pair<int, int>,
            map<pair<int, int>, 
                pair<float, float> > >& counts,
    map<int, map<int, double> >& llrs,
    llr_table& tab,
    int n_samples,
    set<int>& allowed_assignments,
    set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    map<int, double>* prior_weights,
    bool incl_contam,
    double contam_rate,
    double contam_rate_var,
    map<pair<int, int>, map<pair<int, int>, double> >* amb_fracs,
    int n_target){
    
    pair<int, int> nullkey = make_pair(-1, -1);
    
    for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator y = 
        counts.begin(); y != counts.end(); ++y){
        
        for (map<pair<int, int>, pair<float, float> >::iterator z = y->second.begin();
            z != y->second.end(); ++z){
            
            if (z->first.first == -1){
                continue;
            }
            
            // Expected alt allele fractions
            double exp1 = (double)y->first.second / 2.0;
            double exp2 = (double)z->first.second / 2.0;
            double exp3 = (double)(y->first.second + z->first.second) / 4.0;
            
            exp1 = adjust_p_err(exp1, error_rate_ref, error_rate_alt);
            exp2 = adjust_p_err(exp2, error_rate_ref, error_rate_alt);
            exp3 = adjust_p_err(exp3, error_rate_ref, error_rate_alt);
            
            double var1 = 0.0, var2 = 0.0, var3 = 0.0;
            
            if (incl_contam && amb_fracs != NULL){
                exp1 = (1.0-contam_rate)*(double)y->first.second/2.0 + 
                    contam_rate*((*amb_fracs)[y->first][nullkey]);
                exp1 = adjust_p_err(exp1, error_rate_ref, error_rate_alt);
                var1 = ((*amb_fracs)[y->first][nullkey] - (double)y->first.second/2.0);
                
                exp2 = (1.0-contam_rate)*(double)z->first.second/2.0 + 
                    contam_rate*((*amb_fracs)[z->first][nullkey]);
                exp2 = adjust_p_err(exp2, error_rate_ref, error_rate_alt);
                var2 = ((*amb_fracs)[z->first][nullkey] - (double)z->first.second/2.0);
                
                exp3 = (1.0-contam_rate)*(double)(y->first.second + z->first.second)/4.0 + 
                    contam_rate*((*amb_fracs)[y->first][z->first]);
                exp3 = adjust_p_err(exp3, error_rate_ref, error_rate_alt);
                var3 = ((*amb_fracs)[y->first][z->first] - (double)(y->first.second + z->first.second)/4.0);
            }

            int i = y->first.first;
            int j = z->first.first;
            int k = hap_comb_to_idx(i, j, n_samples);

            int ref = (int)round(z->second.first);
            int alt = (int)round(z->second.second);
            
            double ll1 = dbinom(ref+alt, alt, exp1);   
            double ll2 = dbinom(ref+alt, alt, exp2);
            double ll3 = dbinom(ref+alt, alt, exp3);
            
            if (incl_contam && contam_rate_var > 0){
                var1 *= (1.0 - error_rate_ref - error_rate_alt);
                var2 *= (1.0 - error_rate_ref - error_rate_alt);
                var3 *= (1.0 - error_rate_ref - error_rate_alt);
                var1 = var1*var1 * contam_rate_var;
                var2 = var2*var2 * contam_rate_var;
                var3 = var3*var3 * contam_rate_var;
                
                if (var1 > 0 && var2 > 0 && var3 > 0){
                    double fac1 = (exp1*(1.0-exp1))/var1 - 1.0;
                    double fac2 = (exp2*(1.0-exp2))/var2 - 1.0;
                    double fac3 = (exp3*(1.0-exp3))/var3 - 1.0;
                    if (fac1 > 0 && fac2 > 0 && fac3 > 0){
                        double a1 = fac1*exp1;
                        double b1 = fac1*(1.0-exp1);
                        double a2 = fac2*exp2;
                        double b2 = fac2*(1.0-exp2);
                        double a3 = fac3*exp3;
                        double b3 = fac3*(1.0-exp3);
                        ll1 = dbetabin(alt, ref+alt, a1, b1);
                        ll2 = dbetabin(alt, ref+alt, a2, b2);
                        ll3 = dbetabin(alt, ref+alt, a3, b3);
                    }
                }
            }

            map<int, double> m;
            if (llrs.count(i) == 0){
                llrs.insert(make_pair(i, m));
            }
            if (llrs.count(j) == 0){
                llrs.insert(make_pair(j, m));
            }
            if (llrs[i].count(j) == 0){
                llrs[i].insert(make_pair(j, 0.0));
            }
            llrs[i][j] += (ll1-ll2);
            
            if (doublet_rate > 0.0){
                if (llrs[i].count(k) == 0){
                    llrs[i].insert(make_pair(k, 0.0));
                }
                if (llrs[j].count(k) == 0){
                    llrs[j].insert(make_pair(k, 0.0));
                }
                
                llrs[i][k] += (ll1-ll3);
                llrs[j][k] += (ll2-ll3);
            }
        }
    }
    
    // Populate LLR table with singlet/singlet comparisons
    for (map<int, map<int, double> >::iterator x = llrs.begin(); x != llrs.end(); ++x){
        for (map<int, double>::iterator y = x->second.begin(); y != x->second.end(); ++y){
            if (y->first < n_samples){
                if (allowed_assignments.size() == 0 || 
                    (allowed_assignments.find(x->first) != allowed_assignments.end() &&
                     allowed_assignments.find(y->first) != allowed_assignments.end())){
                     
                     if (prior_weights != NULL && 
                        prior_weights->count(x->first) > 0 && prior_weights->count(y->first) > 0){
                        y->second += (*prior_weights)[x->first] - (*prior_weights)[y->first];
                     }
                     tab.insert(x->first, y->first, y->second);
                }
            }
        }
    }
    
    if (doublet_rate > 0.0){ 
        // n_target controls singlet pruning:
        //   -1 = auto: skip pruning if -i/-I provided, else use 10
        //    0 = never prune
        //   >0 = use this value
        int effective_n_target;
        if (n_target == 0){
            effective_n_target = 9999;  // effectively no pruning
        }
        else if (n_target > 0){
            effective_n_target = n_target;
        }
        else{
            // n_target == -1 (auto)
            if (allowed_assignments.size() > 0){
                effective_n_target = 9999;  // don't prune when -i/-I given
            }
            else{
                effective_n_target = 10;  // default
            }
        }
        
        if (g_verbose){
            fprintf(stderr, "[VERBOSE] populate_llr_table: n_indvs=%d, effective_n_target=%d\n",
                tab.n_indvs, effective_n_target);
        }
        
        if (tab.n_indvs > effective_n_target){
            tab.del(effective_n_target);
        }
        
        if (tab.n_indvs < 2 && tab.n_indvs < n_samples){
            return false;
        }
        
        vector<int> ks;
        for (int i = 0; i < n_samples-1; ++i){
            if (allowed_assignments.size() != 0 && 
                allowed_assignments.find(i) == allowed_assignments.end()){
                continue;
            }
            else if (!tab.included[i]){
                continue;
            }

            for (int j = i + 1; j < n_samples; ++j){
                if (allowed_assignments.size() != 0 && 
                    allowed_assignments.find(j) == allowed_assignments.end()){
                    continue;
                }
                else if (!tab.included[j]){
                    continue;
                }
                int k = hap_comb_to_idx(i, j, n_samples);
                if (allowed_assignments.size() == 0 || allowed_assignments.find(k) != 
                    allowed_assignments.end()){
                    ks.push_back(k);
                }
            }
        }
        
        if (ks.size() > 0){
            sort(ks.begin(), ks.end());
            compute_k_comps(llrs, tab, ks, n_samples, allowed_assignments, doublet_rate, prior_weights);
        }
    }
    
    // Disallow impossible combinations
    if (allowed_assignments.size() > 0 || doublet_rate == 1.0){
        for (int i = 0; i < n_samples; ++i){
            if (doublet_rate == 1 || (allowed_assignments2.size() > 0 && 
                allowed_assignments2.find(i) == allowed_assignments2.end())){
                tab.disallow(i);
            }
        }
    }

    return true;
}

// ============================================================================
// NEW OPTIMIZED LLR POPULATION (for CellCounts structure)
// ============================================================================

bool populate_llr_table_optimized(
    const CellCounts& counts,
    map<int, map<int, double> >& llrs,
    llr_table& tab,
    int n_samples,
    set<int>& allowed_assignments,
    set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    map<int, double>* prior_weights,
    bool incl_contam,
    double contam_rate,
    double contam_rate_var,
    map<pair<int, int>, map<pair<int, int>, double> >* amb_fracs,
    int n_target){
    
    // Iterate over all individual pairs and genotype combinations
    for (int i = 0; i < n_samples; ++i){
        for (int nalt_i = 0; nalt_i < 3; ++nalt_i){
            // Get total counts for individual i with this genotype
            auto total_i = counts.get_total(i, nalt_i);
            if (total_i.first + total_i.second == 0) continue;
            
            for (int j = i + 1; j < n_samples; ++j){
                for (int nalt_j = 0; nalt_j < 3; ++nalt_j){
                    // Get pairwise counts
                    auto pair_counts = counts.get(i, nalt_i, j, nalt_j);
                    
                    float ref_count = pair_counts.first;
                    float alt_count = pair_counts.second;
                    
                    if (ref_count + alt_count == 0) continue;
                    
                    // If same genotype, we can't distinguish between the
                    // two individuals from this piece of information
                    if (nalt_i == nalt_j) continue;
                    
                    // Expected alt allele fractions for each hypothesis
                    double exp1 = (double)nalt_i / 2.0;  // Cell is individual i
                    double exp2 = (double)nalt_j / 2.0;  // Cell is individual j
                    double exp3 = (double)(nalt_i + nalt_j) / 4.0;  // Cell is doublet i+j
                    
                    // Apply error rate correction
                    exp1 = adjust_p_err(exp1, error_rate_ref, error_rate_alt);
                    exp2 = adjust_p_err(exp2, error_rate_ref, error_rate_alt);
                    exp3 = adjust_p_err(exp3, error_rate_ref, error_rate_alt);
                    
                    // Handle contamination if requested
                    double var1 = 0.0, var2 = 0.0, var3 = 0.0;
                    if (incl_contam && amb_fracs != NULL){
                        pair<int, int> key_i = make_pair(i, nalt_i);
                        pair<int, int> key_j = make_pair(j, nalt_j);
                        pair<int, int> nullkey = make_pair(-1, -1);
                        
                        if (amb_fracs->count(key_i) > 0 && (*amb_fracs)[key_i].count(nullkey) > 0){
                            double amb_i = (*amb_fracs)[key_i].at(nullkey);
                            exp1 = (1.0-contam_rate)*(double)nalt_i/2.0 + contam_rate*amb_i;
                            exp1 = adjust_p_err(exp1, error_rate_ref, error_rate_alt);
                            var1 = (amb_i - (double)nalt_i/2.0);
                        }
                        if (amb_fracs->count(key_j) > 0 && (*amb_fracs)[key_j].count(nullkey) > 0){
                            double amb_j = (*amb_fracs)[key_j].at(nullkey);
                            exp2 = (1.0-contam_rate)*(double)nalt_j/2.0 + contam_rate*amb_j;
                            exp2 = adjust_p_err(exp2, error_rate_ref, error_rate_alt);
                            var2 = (amb_j - (double)nalt_j/2.0);
                        }
                    }
                    
                    // Compute log likelihoods
                    int ref = (int)round(ref_count);
                    int alt = (int)round(alt_count);
                    int n = ref + alt;
                    
                    double ll1 = dbinom(n, alt, exp1);
                    double ll2 = dbinom(n, alt, exp2);
                    double ll3 = dbinom(n, alt, exp3);
                    
                    // Apply beta-binomial if variance is provided
                    if (incl_contam && contam_rate_var > 0){
                        var1 *= (1.0 - error_rate_ref - error_rate_alt);
                        var2 *= (1.0 - error_rate_ref - error_rate_alt);
                        var3 *= (1.0 - error_rate_ref - error_rate_alt);
                        var1 = var1*var1 * contam_rate_var;
                        var2 = var2*var2 * contam_rate_var;
                        var3 = var3*var3 * contam_rate_var;
                        
                        if (var1 > 0 && var2 > 0 && var3 > 0){
                            double fac1 = (exp1*(1.0-exp1))/var1 - 1.0;
                            double fac2 = (exp2*(1.0-exp2))/var2 - 1.0;
                            double fac3 = (exp3*(1.0-exp3))/var3 - 1.0;
                            if (fac1 > 0 && fac2 > 0 && fac3 > 0){
                                double a1 = fac1*exp1, b1 = fac1*(1.0-exp1);
                                double a2 = fac2*exp2, b2 = fac2*(1.0-exp2);
                                double a3 = fac3*exp3, b3 = fac3*(1.0-exp3);
                                ll1 = dbetabin(alt, n, a1, b1);
                                ll2 = dbetabin(alt, n, a2, b2);
                                ll3 = dbetabin(alt, n, a3, b3);
                            }
                        }
                    }
                    
                    // Accumulate LLRs
                    int k = hap_comb_to_idx(i, j, n_samples);
                    
                    if (llrs.count(i) == 0){
                        llrs[i] = map<int, double>();
                    }
                    if (llrs.count(j) == 0){
                        llrs[j] = map<int, double>();
                    }
                    if (llrs[i].count(j) == 0){
                        llrs[i][j] = 0.0;
                    }
                    llrs[i][j] += (ll1 - ll2);
                    
                    if (doublet_rate > 0.0){
                        if (llrs[i].count(k) == 0){
                            llrs[i][k] = 0.0;
                        }
                        if (llrs[j].count(k) == 0){
                            llrs[j][k] = 0.0;
                        }
                        llrs[i][k] += (ll1 - ll3);
                        llrs[j][k] += (ll2 - ll3);
                    }
                }
            }
        }
    }
    
    // Populate LLR table with singlet/singlet comparisons
    for (auto& x : llrs){
        for (auto& y : x.second){
            if (y.first < n_samples){
                if (allowed_assignments.size() == 0 || 
                    (allowed_assignments.find(x.first) != allowed_assignments.end() &&
                     allowed_assignments.find(y.first) != allowed_assignments.end())){
                     
                    double llr = y.second;
                    if (prior_weights != NULL){
                        if (prior_weights->count(x.first) > 0){
                            llr += (*prior_weights)[x.first];
                        }
                        if (prior_weights->count(y.first) > 0){
                            llr -= (*prior_weights)[y.first];
                        }
                    }
                    tab.insert(x.first, y.first, llr);
                }
            }
        }
    }
    
    // Handle doublet identities
    if (doublet_rate > 0.0){
        // n_target controls singlet pruning:
        //   -1 = auto: skip pruning if -i/-I provided, else use 10
        //    0 = never prune
        //   >0 = use this value
        int effective_n_target;
        if (n_target == 0){
            effective_n_target = 9999;  // effectively no pruning
        }
        else if (n_target > 0){
            effective_n_target = n_target;
        }
        else{
            // n_target == -1 (auto)
            if (allowed_assignments.size() > 0){
                effective_n_target = 9999;  // don't prune when -i/-I given
            }
            else{
                effective_n_target = 10;  // default
            }
        }
        
        if (g_verbose){
            fprintf(stderr, "[VERBOSE] populate_llr_table_optimized: n_indvs=%d, effective_n_target=%d\n",
                tab.n_indvs, effective_n_target);
        }
        
        if (tab.n_indvs > effective_n_target){
            tab.del(effective_n_target);
        }
        
        if (tab.n_indvs < 2 && tab.n_indvs < n_samples){
            return false;
        }
        
        vector<int> ks;
        for (int i = 0; i < n_samples-1; ++i){
            if (allowed_assignments.size() != 0 && 
                allowed_assignments.find(i) == allowed_assignments.end()){
                continue;
            }
            if (!tab.included[i]) continue;

            for (int j = i + 1; j < n_samples; ++j){
                if (allowed_assignments.size() != 0 && 
                    allowed_assignments.find(j) == allowed_assignments.end()){
                    continue;
                }
                if (!tab.included[j]) continue;
                
                int k = hap_comb_to_idx(i, j, n_samples);
                if (allowed_assignments.size() == 0 || 
                    allowed_assignments.find(k) != allowed_assignments.end()){
                    ks.push_back(k);
                }
            }
        }
        
        if (ks.size() > 0){
            sort(ks.begin(), ks.end());
            compute_k_comps(llrs, tab, ks, n_samples, allowed_assignments, doublet_rate, prior_weights);
        }
    }
    
    // Disallow impossible combinations
    if (allowed_assignments.size() > 0 || doublet_rate == 1.0){
        for (int i = 0; i < n_samples; ++i){
            if (doublet_rate == 1 || (allowed_assignments2.size() > 0 && 
                allowed_assignments2.find(i) == allowed_assignments2.end())){
                tab.disallow(i);
            }
        }
        // Recalculate minllr/maxllr now that some identities are disallowed
        // This is critical: minllr values may reflect comparisons against
        // now-disallowed singlets, causing valid doublets to have negative LLRs
        tab.recalculate_minmax();
    }

    return true;
}

// ============================================================================
// DIAGNOSTIC EXTRACTION FUNCTIONS (NEW)
// ============================================================================

/**
 * Extract diagnostic information from LLR table before destruction
 */
void get_diagnostics_from_llrs(
    const map<int, map<int, double> >& llrs,
    const llr_table& tab,
    int winner,
    double winner_llr,
    int n_samples,
    int n_runner_ups,
    double close_threshold,
    CellDiagnostics& diag,
    vector<RunnerUp>& runners) {

    (void)llrs;
    (void)n_samples;
    runners.clear();

    diag.maximin_candidate = winner;
    diag.maximin_score = winner_llr;
    diag.min_margin = tab.get_min_margin(winner);
    tab.winner_has_complete_comparisons(
        winner, diag.missing_comparison_alternatives);
    diag.selection_resolved = winner >= 0 && winner_llr > 0.0 &&
                              diag.missing_comparison_alternatives.empty();

    const vector<int> retained = tab.retained_identities();
    const vector<double>& all_minllr = tab.get_minllr();
    vector<pair<double, int> > identity_scores;
    identity_scores.reserve(retained.size());

    double worst_value = std::numeric_limits<double>::infinity();
    bool saw_present_competitor = false;
    bool saw_partial_competitor = false;
    bool saw_any_competitor = false;

    for (int other : retained) {
        if (other == winner) continue;
        saw_any_competitor = true;

        const PairwiseComparison direct = tab.get_pairwise(winner, other);
        if (direct.state() == ComparisonState::PARTIAL_SUPPORT) {
            saw_partial_competitor = true;
        }
        if (direct.present && direct.value < worst_value) {
            worst_value = direct.value;
            diag.worst_competitor = other;
            diag.worst_comparison_state = direct.state();
            saw_present_competitor = true;
        }

        const double score = (other >= 0 && other < (int)all_minllr.size())
            ? all_minllr[other]
            : -std::numeric_limits<double>::infinity();
        identity_scores.push_back(make_pair(score, other));
        if (std::isfinite(diag.min_margin) && std::isfinite(score) &&
            (diag.min_margin - score) < close_threshold) {
            ++diag.n_close;
        }
    }

    if (!saw_any_competitor) {
        diag.worst_competitor = -1;
        diag.worst_comparison_state = ComparisonState::NOT_APPLICABLE;
    } else if (!saw_present_competitor) {
        diag.worst_competitor = -1;
        diag.worst_comparison_state = saw_partial_competitor
            ? ComparisonState::PARTIAL_SUPPORT
            : ComparisonState::UNAVAILABLE;
    }

    sort(identity_scores.begin(), identity_scores.end(),
         [](const pair<double, int>& a, const pair<double, int>& b) {
             if (a.first != b.first) return a.first > b.first;
             return a.second < b.second;
         });

    const int requested = std::max(0, n_runner_ups);
    const int n_to_extract = std::min((int)identity_scores.size(), requested);
    for (int i = 0; i < n_to_extract; ++i) {
        const int runner_id = identity_scores[i].second;
        const double runner_margin = identity_scores[i].first;
        const PairwiseComparison winner_vs_runner = tab.get_pairwise(winner, runner_id);
        const double runner_vs_winner = winner_vs_runner.present
            ? -winner_vs_runner.value
            : std::numeric_limits<double>::quiet_NaN();
        runners.push_back(RunnerUp(
            runner_id,
            runner_vs_winner,
            winner_vs_runner.state(),
            runner_margin));
    }

    if (runners.empty()) {
        diag.llr_vs_runnerup = std::numeric_limits<double>::quiet_NaN();
        diag.runnerup_comparison_state = ComparisonState::NOT_APPLICABLE;
    } else {
        const PairwiseComparison winner_vs_runner =
            tab.get_pairwise(winner, runners.front().identity);
        diag.llr_vs_runnerup = winner_vs_runner.present
            ? winner_vs_runner.value
            : std::numeric_limits<double>::quiet_NaN();
        diag.runnerup_comparison_state = winner_vs_runner.state();
    }
}

/**
 * Compute het site balance statistics for ploidy detection
 */
void compute_het_balance_stats(
    const CellCounts& het_counts,
    int assigned_id,
    int n_samples,
    CellDiagnostics& diag){
    
    vector<double> alt_fracs;
    double total_depth = 0.0;
    
    bool is_doublet = (assigned_id >= n_samples);
    
    if (is_doublet){
        // Doublet: use sites where BOTH component individuals are het
        pair<int, int> combo = idx_to_hap_comb(assigned_id, n_samples);
        int indv1 = combo.first;
        int indv2 = combo.second;
        
        // Get counts at sites where both are het (genotype = 1 for both)
        auto both_het = het_counts.get(indv1, 1, indv2, 1);
        float ref = both_het.first;
        float alt = both_het.second;
        
        if (ref + alt > 0){
            alt_fracs.push_back(alt / (ref + alt));
            total_depth += (ref + alt);
        }
        
        // Also consider sites where one is het and the other has known genotype
        // This gives us more data points
        for (int gt1 = 0; gt1 < 3; ++gt1){
            for (int gt2 = 0; gt2 < 3; ++gt2){
                // Skip if neither is het
                if (gt1 != 1 && gt2 != 1) continue;
                // Skip the both-het case (already handled)
                if (gt1 == 1 && gt2 == 1) continue;
                
                auto counts = het_counts.get(indv1, gt1, indv2, gt2);
                float r = counts.first;
                float a = counts.second;
                if (r + a > 0){
                    alt_fracs.push_back(a / (r + a));
                    total_depth += (r + a);
                }
            }
        }
    }
    else{
        // Singlet: use sites where the assigned individual is het
        int indv = assigned_id;
        
        // Get all het sites for this individual
        for (int other = 0; other < n_samples; ++other){
            if (other == indv) continue;
            
            for (int gt_other = 0; gt_other < 3; ++gt_other){
                // indv is het (genotype = 1)
                pair<float, float> counts;
                if (indv < other){
                    counts = het_counts.get(indv, 1, other, gt_other);
                }
                else{
                    counts = het_counts.get(other, gt_other, indv, 1);
                }
                
                float ref = counts.first;
                float alt = counts.second;
                if (ref + alt > 0){
                    alt_fracs.push_back(alt / (ref + alt));
                    total_depth += (ref + alt);
                }
            }
        }
        
        // Also use totals for het genotype
        auto total = het_counts.get_total(indv, 1);
        if (total.first + total.second > 0 && alt_fracs.empty()){
            // Fallback: use total counts if no pairwise data
            alt_fracs.push_back(total.second / (total.first + total.second));
            total_depth = total.first + total.second;
        }
    }
    
    diag.n_het_sites = (int)alt_fracs.size();
    diag.het_total_depth = total_depth;
    
    // Compute variance of alt_frac
    // Minimum number of sites needed for reliable variance estimate
    const int MIN_HET_SITES = 10;
    
    if (alt_fracs.size() < MIN_HET_SITES){
        diag.het_balance_var = std::numeric_limits<double>::quiet_NaN();
        diag.het_diagnostic_available = false;
        return;
    }
    
    // Compute mean
    double sum = 0.0;
    for (double af : alt_fracs){
        sum += af;
    }
    double mean = sum / alt_fracs.size();
    
    // Compute variance
    double var_sum = 0.0;
    for (double af : alt_fracs){
        double diff = af - mean;
        var_sum += diff * diff;
    }
    diag.het_balance_var = var_sum / (alt_fracs.size() - 1);  // Sample variance
    diag.het_diagnostic_available = true;
}

/**
 * Compute total depth from main demux counts
 */
double compute_total_depth(const CellCounts& counts, int n_samples){
    double total = 0.0;
    for (int indv = 0; indv < n_samples; ++indv){
        for (int nalt = 0; nalt < 3; ++nalt){
            auto tot = counts.get_total(indv, nalt);
            total += tot.first + tot.second;
        }
    }
    // Divide by n_samples since each read contributes to multiple individuals
    return total / n_samples;
}

// ============================================================================
// MARGIN SOFTMAX SCORE AND ENTROPY COMPUTATION
// ============================================================================

void compute_margin_softmax_scores(
    const map<int, map<int, double> >& llrs,
    const llr_table& tab,
    int winner,
    int n_samples,
    CellDiagnostics& diag) {

    (void)llrs;
    (void)n_samples;
    if (winner < 0 || winner >= (int)tab.included.size() ||
        !tab.included[winner]) {
        diag.margin_softmax_score = std::numeric_limits<double>::quiet_NaN();
        diag.margin_entropy = std::numeric_limits<double>::quiet_NaN();
        return;
    }

    const vector<double>& margins = tab.get_minllr();
    const double winner_margin = margins[winner];
    if (!std::isfinite(winner_margin)) {
        diag.margin_softmax_score = std::numeric_limits<double>::quiet_NaN();
        diag.margin_entropy = std::numeric_limits<double>::quiet_NaN();
        return;
    }

    vector<double> relative_margins;
    relative_margins.reserve(tab.n_indvs > 0 ? (size_t)tab.n_indvs : 0U);
    for (size_t i = 0; i < tab.included.size(); ++i) {
        if (!tab.included[i] || !std::isfinite(margins[i])) continue;
        relative_margins.push_back(margins[i] - winner_margin);
    }

    if (relative_margins.empty()) {
        diag.margin_softmax_score = std::numeric_limits<double>::quiet_NaN();
        diag.margin_entropy = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    if (relative_margins.size() == 1U) {
        diag.margin_softmax_score = 1.0;
        diag.margin_entropy = 0.0;
        return;
    }

    double sum_exp = 0.0;
    for (double value : relative_margins) sum_exp += exp(value);
    if (!(sum_exp > 0.0) || !std::isfinite(sum_exp)) {
        diag.margin_softmax_score = std::numeric_limits<double>::quiet_NaN();
        diag.margin_entropy = std::numeric_limits<double>::quiet_NaN();
        return;
    }

    diag.margin_softmax_score = 1.0 / sum_exp;
    double entropy = 0.0;
    for (double value : relative_margins) {
        const double p = exp(value) / sum_exp;
        if (p > 0.0) entropy -= p * log2(p);
    }
    diag.margin_entropy = entropy;
}

/**
 * Compute het balance using per-site data (PERSITE method)
 */
void compute_het_balance_persite(
    const HetSiteData& persite_data,
    const robin_hood::unordered_map<int, ChromSNPs>& het_snpdat,
    const vector<pair<int, int>>& idx_to_site,
    int assigned_id,
    int n_samples,
    double min_depth,
    int min_sites,
    CellDiagnostics& diag) {
    
    vector<double> alt_fracs;
    double total_depth = 0.0;
    
    bool is_doublet = (assigned_id >= n_samples);
    int indv1 = -1, indv2 = -1;
    
    if (is_doublet) {
        auto combo = idx_to_hap_comb(assigned_id, n_samples);
        indv1 = combo.first;
        indv2 = combo.second;
    } else {
        indv1 = assigned_id;
    }
    
    for (const auto& site : persite_data.sites) {
        if (site.site_idx < 0 || site.site_idx >= (int)idx_to_site.size()) continue;
        
        int tid = idx_to_site[site.site_idx].first;
        int pos = idx_to_site[site.site_idx].second;
        
        auto chrom_it = het_snpdat.find(tid);
        if (chrom_it == het_snpdat.end()) continue;
        
        const ChromSNPs& chrom_snps = chrom_it->second;
        SNPData target;
        target.pos = pos;
        auto snp_it = lower_bound(chrom_snps.snps.begin(), chrom_snps.snps.end(), target);
        if (snp_it == chrom_snps.snps.end() || snp_it->pos != pos) continue;
        
        const var& v = snp_it->data;
        
        bool use_site = false;
        if (is_doublet) {
            use_site = v.is_het(indv1) || v.is_het(indv2);
        } else {
            use_site = v.is_het(indv1);
        }
        
        if (!use_site) continue;
        
        double depth = site.ref + site.alt;
        if (depth < min_depth) continue;
        
        double alt_frac = site.alt / depth;
        alt_fracs.push_back(alt_frac);
        total_depth += depth;
    }
    
    diag.n_het_sites = (int)alt_fracs.size();
    diag.het_total_depth = total_depth;
    diag.het_method = HetBalanceMethod::PERSITE;
    
    if ((int)alt_fracs.size() < min_sites) {
        diag.het_balance_var = std::numeric_limits<double>::quiet_NaN();
        diag.het_diagnostic_available = false;
        return;
    }
    
    double sum = 0.0;
    for (double af : alt_fracs) sum += af;
    double mean = sum / alt_fracs.size();
    
    double var_sum = 0.0;
    for (double af : alt_fracs) {
        double diff = af - mean;
        var_sum += diff * diff;
    }
    diag.het_balance_var = var_sum / (alt_fracs.size() - 1);
    diag.het_diagnostic_available = true;
}

/**
 * Compute het balance using Welford stats (WELFORD method)
 */
void compute_het_balance_welford(
    const CellWelfordStats& welford_stats,
    int assigned_id,
    int n_samples,
    int min_sites,
    CellDiagnostics& diag) {

    diag.het_method = HetBalanceMethod::WELFORD;
    diag.het_balance_var = std::numeric_limits<double>::quiet_NaN();
    diag.het_diagnostic_available = false;

    const bool is_doublet = assigned_id >= n_samples;
    if (is_doublet) {
        const pair<int, int> combo = idx_to_hap_comb(assigned_id, n_samples);
        const WelfordStats& ws1 = welford_stats.get(combo.first);
        const WelfordStats& ws2 = welford_stats.get(combo.second);
        const double var1 = ws1.variance(min_sites);
        const double var2 = ws2.variance(min_sites);
        const bool valid1 = std::isfinite(var1) && var1 >= 0.0;
        const bool valid2 = std::isfinite(var2) && var2 >= 0.0;

        if (!valid1 && !valid2) {
            diag.n_het_sites = 0;
            diag.het_total_depth = 0.0;
            return;
        }
        if (!valid1) {
            diag.het_balance_var = var2;
            diag.n_het_sites = static_cast<int>(ws2.n);
            diag.het_total_depth = ws2.total_depth;
        } else if (!valid2) {
            diag.het_balance_var = var1;
            diag.n_het_sites = static_cast<int>(ws1.n);
            diag.het_total_depth = ws1.total_depth;
        } else {
            const double n_total = ws1.n + ws2.n;
            if (!(n_total > 0.0)) return;
            diag.het_balance_var = (var1 * ws1.n + var2 * ws2.n) / n_total;
            diag.n_het_sites = static_cast<int>(n_total);
            diag.het_total_depth = ws1.total_depth + ws2.total_depth;
        }
        diag.het_diagnostic_available = true;
        return;
    }

    const WelfordStats& ws = welford_stats.get(assigned_id);
    const double variance = ws.variance(min_sites);
    diag.n_het_sites = static_cast<int>(ws.n);
    diag.het_total_depth = ws.total_depth;
    if (std::isfinite(variance) && variance >= 0.0) {
        diag.het_balance_var = variance;
        diag.het_diagnostic_available = true;
    }
}

// ============================================================================
// PARALLEL IDENTITY ASSIGNMENT
// ============================================================================

static bool accepted_maximin_candidate(
    const llr_table& tab,
    int candidate,
    double score,
    vector<int>& missing_alternatives) {
    if (candidate < 0 || !(score > 0.0) || !std::isfinite(score)) {
        missing_alternatives.clear();
        return false;
    }
    return tab.winner_has_complete_comparisons(candidate, missing_alternatives);
}

static bool validate_assignment_dimensions(int n_samples) {
    string error;
    if (!validate_identity_and_allocation_request(n_samples, NULL, NULL, &error)) {
        fprintf(stderr, "ERROR: cannot score requested identity universe: %s\n", error.c_str());
        return false;
    }
    return true;
}

static bool reject_unsupported_identity_prior(const map<int, double>* identity_prior) {
    if (identity_prior == nullptr) return false;
    fprintf(stderr,
        "ERROR: --identity_prior is NOT IMPLEMENTED / NOT USED IN IDENTITY SCORING. "
        "Refusing to continue rather than silently ignoring it.\n");
    return true;
}

static void report_auxiliary_evidence_scope(
    const robin_hood::unordered_map<unsigned long, CellCounts>* atac_cell_counts,
    const robin_hood::unordered_map<unsigned long, CellCounts>* species_counts_rna,
    const robin_hood::unordered_map<unsigned long, CellCounts>* species_counts_atac) {
    if (atac_cell_counts != nullptr) {
        fprintf(stderr,
            "STATUS: ATAC counts are diagnostic/count-only and do not affect individual identity scoring.\n");
    }
    if (species_counts_rna != nullptr || species_counts_atac != nullptr) {
        fprintf(stderr,
            "STATUS: species counts are diagnostic/count-only and do not affect individual identity scoring.\n");
    }
}

bool assign_ids_parallel(
    robin_hood::unordered_map<unsigned long, CellCounts>& cell_counts,
    vector<string>& samples,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, double>& assignments_llr,
    set<int>& allowed_assignments,
    set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    bool use_prior_weights,
    map<int, double>& prior_weights,
    int n_threads,
    int n_target) {

    assignments.clear();
    assignments_llr.clear();
    const int n_samples = static_cast<int>(samples.size());
    if (!validate_assignment_dimensions(n_samples)) return false;

    if (g_verbose) {
        fprintf(stderr,
            "[VERBOSE] assign_ids_parallel: n_samples=%d, n_target=%d, doublet_rate=%.3f\n",
            n_samples, n_target, doublet_rate);
        fprintf(stderr, "[VERBOSE]   allowed_assignments: %lu, allowed_assignments2: %lu\n",
            allowed_assignments.size(), allowed_assignments2.size());
    }

    vector<unsigned long> barcodes;
    barcodes.reserve(cell_counts.size());
    for (const auto& kv : cell_counts) barcodes.push_back(kv.first);

    vector<int> result_assn(barcodes.size(), -1);
    vector<double> result_llr(barcodes.size(), std::numeric_limits<double>::quiet_NaN());

    fprintf(stderr, "Assigning identities to %lu cells using %d threads...\n",
        barcodes.size(), n_threads);

    #pragma omp parallel for num_threads(n_threads) schedule(dynamic, 100)
    for (size_t idx = 0; idx < barcodes.size(); ++idx) {
        const unsigned long bc = barcodes[idx];
        const CellCounts& counts = cell_counts[bc];
        if (counts.is_empty()) continue;

        map<int, map<int, double> > llrs;
        llr_table tab(n_samples);
        map<int, double>* pw_ptr = use_prior_weights ? &prior_weights : NULL;
        const bool success = populate_llr_table_optimized(
            counts, llrs, tab, n_samples,
            allowed_assignments, allowed_assignments2,
            doublet_rate, error_rate_ref, error_rate_alt,
            pw_ptr, false, 0.0, 0.0, NULL, n_target);
        if (!success) continue;

        int candidate = -1;
        double score = -std::numeric_limits<double>::infinity();
        tab.get_max(candidate, score);
        vector<int> missing;
        if (accepted_maximin_candidate(tab, candidate, score, missing)) {
            result_assn[idx] = candidate;
            result_llr[idx] = score;
        }
    }

    for (size_t idx = 0; idx < barcodes.size(); ++idx) {
        if (result_assn[idx] < 0) continue;
        assignments.emplace(barcodes[idx], result_assn[idx]);
        assignments_llr.emplace(barcodes[idx], result_llr[idx]);
    }

    fprintf(stderr, "Assigned %lu cells\n", assignments.size());
    if (g_verbose) {
        int n_singlets = 0;
        int n_doublets = 0;
        map<int, int> assignment_counts;
        for (const auto& assignment : assignments) {
            ++assignment_counts[assignment.second];
            if (assignment.second < n_samples) ++n_singlets;
            else ++n_doublets;
        }
        fprintf(stderr, "[VERBOSE] Assignment summary: %d singlets, %d doublets\n",
            n_singlets, n_doublets);
        for (const auto& entry : assignment_counts) {
            if (entry.first < n_samples) {
                fprintf(stderr, "[VERBOSE]   %s: %d cells\n",
                    samples[entry.first].c_str(), entry.second);
            } else {
                const pair<int, int> combo = idx_to_hap_comb(entry.first, n_samples);
                fprintf(stderr, "[VERBOSE]   %s+%s: %d cells\n",
                    samples[combo.first].c_str(), samples[combo.second].c_str(), entry.second);
            }
        }
    }
    return true;
}

bool assign_ids_parallel_with_diagnostics(
    robin_hood::unordered_map<unsigned long, CellCounts>& cell_counts,
    vector<string>& samples,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, double>& assignments_llr,
    set<int>& allowed_assignments,
    set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    bool use_prior_weights,
    map<int, double>& prior_weights,
    int n_threads,
    int n_target,
    bool compute_diagnostics,
    int n_runner_ups,
    double close_threshold,
    robin_hood::unordered_map<unsigned long, CellCounts>* het_counts,
    robin_hood::unordered_map<unsigned long, CellDiagnostics>& diagnostics,
    robin_hood::unordered_map<unsigned long, vector<RunnerUp> >& runner_ups,
    robin_hood::unordered_map<unsigned long, CellCounts>* atac_cell_counts,
    const map<int, double>* identity_prior,
    double z_doublet,
    robin_hood::unordered_map<unsigned long, CellCounts>* species_counts_rna,
    robin_hood::unordered_map<unsigned long, CellCounts>* species_counts_atac) {

    (void)z_doublet;
    assignments.clear();
    assignments_llr.clear();
    diagnostics.clear();
    runner_ups.clear();

    const int n_samples = static_cast<int>(samples.size());
    if (!validate_assignment_dimensions(n_samples) ||
        reject_unsupported_identity_prior(identity_prior)) return false;
    report_auxiliary_evidence_scope(atac_cell_counts, species_counts_rna, species_counts_atac);

    vector<unsigned long> barcodes;
    barcodes.reserve(cell_counts.size());
    for (const auto& kv : cell_counts) barcodes.push_back(kv.first);

    vector<int> result_assn(barcodes.size(), -1);
    vector<double> result_llr(barcodes.size(), std::numeric_limits<double>::quiet_NaN());
    vector<bool> result_has_diag(barcodes.size(), false);
    vector<CellDiagnostics> result_diag(barcodes.size());
    vector<vector<RunnerUp> > result_runners(barcodes.size());

    fprintf(stderr, "Assigning identities to %lu cells using %d threads%s...\n",
        barcodes.size(), n_threads, compute_diagnostics ? " (with diagnostics)" : "");

    #pragma omp parallel for num_threads(n_threads) schedule(dynamic, 100)
    for (size_t idx = 0; idx < barcodes.size(); ++idx) {
        const unsigned long bc = barcodes[idx];
        const CellCounts& counts = cell_counts[bc];
        if (counts.is_empty()) continue;

        map<int, map<int, double> > llrs;
        llr_table tab(n_samples);
        map<int, double>* pw_ptr = use_prior_weights ? &prior_weights : NULL;
        const bool success = populate_llr_table_optimized(
            counts, llrs, tab, n_samples,
            allowed_assignments, allowed_assignments2,
            doublet_rate, error_rate_ref, error_rate_alt,
            pw_ptr, false, 0.0, 0.0, NULL, n_target);
        if (!success) continue;

        int candidate = -1;
        double score = -std::numeric_limits<double>::infinity();
        tab.get_max(candidate, score);
        vector<int> missing;
        const bool accepted = accepted_maximin_candidate(tab, candidate, score, missing);
        if (accepted) {
            result_assn[idx] = candidate;
            result_llr[idx] = score;
        }

        if (compute_diagnostics && candidate >= 0) {
            CellDiagnostics diag;
            vector<RunnerUp> runners;
            get_diagnostics_from_llrs(llrs, tab, candidate, score, n_samples,
                n_runner_ups, close_threshold, diag, runners);
            compute_margin_softmax_scores(llrs, tab, candidate, n_samples, diag);
            tab.get_max_by_max_llr_comparator(
                diag.max_llr_comparator_winner,
                diag.max_llr_comparator_score);
            diag.total_depth = compute_total_depth(counts, n_samples);

            // Het/ploidy is read-only with respect to the frozen accepted identity.
            if (accepted && het_counts != NULL) {
                const auto het_it = het_counts->find(bc);
                if (het_it != het_counts->end()) {
                    compute_het_balance_stats(het_it->second, candidate, n_samples, diag);
                }
            }
            result_has_diag[idx] = true;
            result_diag[idx] = diag;
            result_runners[idx] = runners;
        }
    }

    for (size_t idx = 0; idx < barcodes.size(); ++idx) {
        const unsigned long bc = barcodes[idx];
        if (result_assn[idx] >= 0) {
            assignments.emplace(bc, result_assn[idx]);
            assignments_llr.emplace(bc, result_llr[idx]);
        }
        if (compute_diagnostics && result_has_diag[idx]) {
            diagnostics.emplace(bc, result_diag[idx]);
            runner_ups.emplace(bc, result_runners[idx]);
        }
    }

    fprintf(stderr, "Assigned %lu cells\n", assignments.size());
    if (compute_diagnostics) {
        fprintf(stderr, "Collected diagnostics for %lu evaluated cells\n", diagnostics.size());
    }
    return true;
}

bool assign_ids_parallel_with_diagnostics_extended(
    robin_hood::unordered_map<unsigned long, CellCounts>& cell_counts,
    vector<string>& samples,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, double>& assignments_llr,
    set<int>& allowed_assignments,
    set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    bool use_prior_weights,
    map<int, double>& prior_weights,
    int n_threads,
    int n_target,
    bool compute_diagnostics,
    int n_runner_ups,
    double close_threshold,
    robin_hood::unordered_map<unsigned long, CellHetData>* het_data,
    const robin_hood::unordered_map<int, ChromSNPs>* het_snpdat,
    const vector<pair<int, int> >* idx_to_site,
    HetBalanceMethod het_method,
    int min_het_sites,
    double min_het_depth,
    robin_hood::unordered_map<unsigned long, CellDiagnostics>& diagnostics,
    robin_hood::unordered_map<unsigned long, vector<RunnerUp> >& runner_ups,
    robin_hood::unordered_map<unsigned long, CellCounts>* atac_cell_counts,
    const map<int, double>* identity_prior,
    double z_doublet,
    robin_hood::unordered_map<unsigned long, CellCounts>* species_counts_rna,
    robin_hood::unordered_map<unsigned long, CellCounts>* species_counts_atac) {

    (void)z_doublet;
    assignments.clear();
    assignments_llr.clear();
    diagnostics.clear();
    runner_ups.clear();

    const int n_samples = static_cast<int>(samples.size());
    if (!validate_assignment_dimensions(n_samples) ||
        reject_unsupported_identity_prior(identity_prior)) return false;
    report_auxiliary_evidence_scope(atac_cell_counts, species_counts_rna, species_counts_atac);

    const char* method_name = het_method == HetBalanceMethod::PERSITE ? "per-site" : "Welford";
    vector<unsigned long> barcodes;
    barcodes.reserve(cell_counts.size());
    for (const auto& kv : cell_counts) barcodes.push_back(kv.first);

    vector<int> result_assn(barcodes.size(), -1);
    vector<double> result_llr(barcodes.size(), std::numeric_limits<double>::quiet_NaN());
    vector<bool> result_has_diag(barcodes.size(), false);
    vector<CellDiagnostics> result_diag(barcodes.size());
    vector<vector<RunnerUp> > result_runners(barcodes.size());
    atomic<int> cells_done(0);
    const int total_cells = static_cast<int>(barcodes.size());

    #pragma omp parallel num_threads(n_threads)
    {
        map<int, map<int, double> > llrs;
        #pragma omp for schedule(dynamic, 100)
        for (size_t idx = 0; idx < barcodes.size(); ++idx) {
            const unsigned long bc = barcodes[idx];
            const CellCounts& counts = cell_counts[bc];
            if (counts.is_empty()) continue;

            llr_table tab(n_samples);
            llrs.clear();
            const bool success = populate_llr_table_optimized(
                counts, llrs, tab, n_samples,
                allowed_assignments, allowed_assignments2,
                doublet_rate, error_rate_ref, error_rate_alt,
                use_prior_weights ? &prior_weights : NULL,
                false, 0.0, 0.0, NULL, n_target);
            if (!success) continue;

            int candidate = -1;
            double score = -std::numeric_limits<double>::infinity();
            tab.get_max(candidate, score);
            vector<int> missing;
            const bool accepted = accepted_maximin_candidate(tab, candidate, score, missing);
            if (accepted) {
                result_assn[idx] = candidate;
                result_llr[idx] = score;
            }

            if (compute_diagnostics && candidate >= 0) {
                CellDiagnostics diag;
                vector<RunnerUp> runners;
                get_diagnostics_from_llrs(llrs, tab, candidate, score,
                    n_samples, n_runner_ups, close_threshold, diag, runners);
                compute_margin_softmax_scores(llrs, tab, candidate, n_samples, diag);
                tab.get_max_by_max_llr_comparator(
                    diag.max_llr_comparator_winner,
                    diag.max_llr_comparator_score);
                diag.total_depth = compute_total_depth(counts, n_samples);

                // Freeze the accepted identity before optional het/ploidy diagnostics.
                if (accepted && het_data != NULL) {
                    const auto het_it = het_data->find(bc);
                    if (het_it != het_data->end()) {
                        const CellHetData& cell_het = het_it->second;
                        if (het_method == HetBalanceMethod::PERSITE &&
                            het_snpdat != NULL && idx_to_site != NULL) {
                            compute_het_balance_persite(
                                cell_het.persite_data, *het_snpdat, *idx_to_site,
                                candidate, n_samples, min_het_depth, min_het_sites, diag);
                        } else {
                            compute_het_balance_welford(
                                cell_het.welford_stats, candidate, n_samples,
                                min_het_sites, diag);
                        }
                    }
                }
                result_has_diag[idx] = true;
                result_diag[idx] = diag;
                result_runners[idx] = runners;
            }

            const int done = ++cells_done;
            if (done % 1000 == 0 || done == total_cells) {
                fprintf(stderr, "\rAssigning: %d/%d cells    ", done, total_cells);
            }
        }
    }

    for (size_t idx = 0; idx < barcodes.size(); ++idx) {
        const unsigned long bc = barcodes[idx];
        if (result_assn[idx] >= 0) {
            assignments.emplace(bc, result_assn[idx]);
            assignments_llr.emplace(bc, result_llr[idx]);
        }
        if (compute_diagnostics && result_has_diag[idx]) {
            diagnostics.emplace(bc, result_diag[idx]);
            runner_ups.emplace(bc, result_runners[idx]);
        }
    }

    fprintf(stderr, "\nAssigned %lu cells (%s het method)\n",
        assignments.size(), method_name);
    if (compute_diagnostics) {
        fprintf(stderr, "Collected diagnostics for %lu evaluated cells\n", diagnostics.size());
    }
    return true;
}


// Historical spelling retained as a source-compatibility wrapper only.
void llr_table::get_max_by_maxllr(int& best_idx, double& best_maxllr) const {
    get_max_by_max_llr_comparator(best_idx, best_maxllr);
}


bool populate_llr_table_peridentity(
    map<pair<int, int>, map<pair<int, int>, pair<float, float> > >& counts,
    map<int, map<int, double> >& llrs,
    llr_table& tab,
    int n_samples,
    set<int>& allowed_assignments,
    set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    map<int, double>* prior_weights,
    bool incl_contam,
    double contam_rate,
    double contam_rate_var,
    map<pair<int, int>, map<pair<int, int>, double> >* amb_fracs,
    int n_target) {
    return populate_llr_table(counts, llrs, tab, n_samples,
        allowed_assignments, allowed_assignments2, doublet_rate,
        error_rate_ref, error_rate_alt, prior_weights, incl_contam,
        contam_rate, contam_rate_var, amb_fracs, n_target);
}


// Pairwise ambient likelihood compatibility helpers.
double lbinom_antider_c(double n, double k, double c, double p_0, double p_c){
    double x = p_c;
    double y = p_0;

    double term1 = -c*(x-y)*(n-binom_coef_log(n,k));
    double term2 = -(-k-n)*(c*(x-y) + y - 1)*log2(-c*x + c*y - y + 1);
    double term3 = k*(c*(x-y) + y)*log2(c*(x-y) + y);
    return (1.0/(x-y))*(term1 + term2 + term3);
}

double lbinom_antider_c2(double n, 
    double k, 
    double c, 
    double p_0, 
    double p_c, 
    double e_r, 
    double e_a){
    
    double x = p_c;
    double y = p_0;
    double r = e_r;
    double a = e_a;

    double term1 = (k*(a-1)*y + r*(y-1))*log2(-c*(a+r-1)*(x-y) - a*y + r*(-y) + r + y);
    term1 /= ((a+r-1)*(x-y));
    double term2 = -k*(a*(y-1) + (r-1)*y+1)*log2(a*(c*(x-y) + y - 1) + c*(r-1)*(x-y) + r*y - y + 1);
    term2 /= ((a+r-1)*(x-y));
    double term3 = -c*k*log2((a+r-1)*(c*(x-y) + y - 1) + r);
    double term4 = c*k*log2(r - (a+r-1)*(c*(x-y) + y));
    double term5 = n*(a*(y-1) + (r-1)*y + 1)*log2(a*(c*(x-y) + y - 1) + c*(r-1)*(x-y) + r*y - y + 1);
    term5 /= ((a+r-1)*(x-y));
    double term6 = c*n*log2((a+r-1)*(c*(x-y) + y - 1) + r);
    double term7 = c*binom_coef_log(n,k) - c*n;
    return term1+term2+term3+term4+term5+term6+term7;
}

bool populate_llr_table_pairwise(map<pair<int, int>, 
    map<pair<int, int>, pair<float, float> > >& counts,
    map<int, map<int, double> >& llrs,
    llr_table& tab,
    int n_samples,
    set<int>& allowed_assignments,
    set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    map<int, double>* prior_weights,
    bool incl_contam,
    double contam_rate,
    double contam_rate_var,
    map<pair<int, int>, map<pair<int, int>, double> >* amb_fracs){
    
    for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator y = 
        counts.begin(); y != counts.end(); ++y){
        
        if (allowed_assignments.size() > 0 && allowed_assignments.find(y->first.first) == 
            allowed_assignments.end()){
            continue;
        }
        
        // Set default expectation for indv1
        // 0 = homozygous ref (~0% alt allele)
        // 1 = heterozygous (~50% alt allele)
        // 2 = homozygous alt (~100% alt allele)
        double exp1 = adjust_p_err((double)y->first.second / 2.0, error_rate_ref, error_rate_alt);
        double var1; 
        double exp1b = (double)y->first.second / 2.0;
        /*
        float exp1 = error_rate_ref;
        if (y->first.second == 1){
            //exp1 = 0.5;
            exp1 = 0.5*(1.0 - error_rate_alt + error_rate_ref);
        }
        else if (y->first.second == 2){
            exp1 = 1.0-error_rate_alt;
        }
        */

        for (map<pair<int, int>, pair<float, float> >::iterator z = 
            y->second.begin(); z != y->second.end(); ++z){
            
            if (allowed_assignments.size() > 0 && allowed_assignments.find(z->first.first) ==
                allowed_assignments.end()){
                continue;
            } 

            // If same site type, we can't distinguish between the
            // two individuals from this piece of information
            if (z->first.first != -1 && y->first.second != z->first.second){

                if (incl_contam){
                    exp1 = (1.0-contam_rate)*((double)y->first.second/2.0) + 
                        contam_rate*((*amb_fracs)[y->first][z->first]);
                    exp1 = adjust_p_err(exp1, error_rate_ref, error_rate_alt);
                    var1 = ((*amb_fracs)[y->first][z->first] - (double)y->first.second/2.0);
                }

                // Set default expectation for indv2
                double exp2 = adjust_p_err((double)z->first.second/2.0, error_rate_ref, error_rate_alt);
                double var2;
                double exp2b = (double)z->first.second/2.0;
                /*
                float exp2 = error_rate_ref;
                if (z->first.second == 1){
                    //exp2 = 0.5;
                    exp2 = 0.5*(1.0 - error_rate_alt + error_rate_ref);
                }
                else if (z->first.second == 2){
                    exp2 = 1.0-error_rate_alt;
                }
                */
                if (incl_contam){
                    exp2 = (1.0-contam_rate)*((double)z->first.second/2.0) + 
                        contam_rate*((*amb_fracs)[y->first][z->first]);
                    exp2 = adjust_p_err(exp2, error_rate_ref, error_rate_alt);
                    var2 = ((*amb_fracs)[y->first][z->first] - (double)z->first.second/2.0);
                }
                
                double exp3 = adjust_p_err((double)(y->first.second + z->first.second)/4.0, 
                    error_rate_ref, error_rate_alt);
                double var3;
                double exp3b = (double)(y->first.second + z->first.second)/4.0;
                /*
                float exp3;
                if (y->first.second == 0 && z->first.second == 0){
                    exp3 = error_rate_ref;
                }
                else if (y->first.second == 2 && z->first.second == 2){
                    exp3 = 1.0 - error_rate_alt;
                }
                else if ((y->first.second == 1 && z->first.second == 1) ||
                        (y->first.second == 0 && z->first.second == 2) ||
                        (y->first.second == 2 && z->first.second == 0)){
                    exp3 = 0.5*( 1.0 - error_rate_alt + error_rate_ref);
                }
                else if ((y->first.second == 0 && z->first.second == 1) ||
                    (y->first.second == 1 && z->first.second == 0)){
                    exp3 = 0.25*(1 - error_rate_alt + 3*error_rate_ref);
                }
                else if ((y->first.second == 1 && z->first.second == 2) ||
                    (y->first.second == 2 && z->first.second == 1)){
                    exp3 = 0.25*(3.0 - 3*error_rate_alt + error_rate_ref);
                }
                */
                if (incl_contam){
                    exp3 = (1.0-contam_rate)*(double)(y->first.second + z->first.second)/4.0 + 
                        contam_rate*((*amb_fracs)[y->first][z->first]);
                    exp3 = adjust_p_err(exp3, error_rate_ref, error_rate_alt);
                    var3 = ((*amb_fracs)[y->first][z->first] - (double)(y->first.second + z->first.second)/4.0);
                }

                int i = y->first.first;
                int j = z->first.first;
                int k = hap_comb_to_idx(i, j, n_samples);

                int ref = (int)round(z->second.first);
                int alt = (int)round(z->second.second);
                
                double ll1 = dbinom(ref+alt, alt, exp1);   
                double ll2 = dbinom(ref+alt, alt, exp2);
                double ll3 = dbinom(ref+alt, alt, exp3);
                
                if (incl_contam && contam_rate_var > 0){
                    /*
                    double p_c = (*amb_fracs)[y->first][z->first];
                    double delta = 0.05; 
                    ll1 = lbinom_antider_c(alt, ref+alt, contam_rate+delta,
                        exp1b, p_c) - 
                        lbinom_antider_c(alt, ref+alt, contam_rate,
                        exp1b, p_c);
                    ll2 = lbinom_antider_c(alt, ref+alt, contam_rate+delta,
                        exp2b, p_c) - 
                        lbinom_antider_c(alt, ref+alt, contam_rate,
                        exp2b, p_c);
                    ll3 = lbinom_antider_c(alt, ref+alt, contam_rate+delta,
                        exp3b, p_c) - 
                        lbinom_antider_c(alt, ref+alt, contam_rate,
                        exp3b, p_c);
                    */
                    /* 
                    ll1 = lbinom_antider_c2(alt, ref+alt, contam_rate+0.001,
                        exp1b, p_c, error_rate_ref, error_rate_alt) - 
                        lbinom_antider_c2(alt, ref+alt, contam_rate-0.001,
                            exp1b, p_c, error_rate_ref, error_rate_alt);
                    ll2 = lbinom_antider_c2(alt, ref+alt, contam_rate+0.001,
                        exp2b, p_c, error_rate_ref, error_rate_alt) - 
                        lbinom_antider_c2(alt, ref+alt, contam_rate-0.001,
                            exp2b, p_c, error_rate_ref, error_rate_alt);
                    ll3 = lbinom_antider_c2(alt, ref+alt, contam_rate+0.001,
                        exp3b, p_c, error_rate_ref, error_rate_alt) - 
                        lbinom_antider_c2(alt, ref+alt, contam_rate-0.001,
                            exp3b, p_c, error_rate_ref, error_rate_alt);
                    */

                    
                    var1 *= (1.0 - error_rate_ref - error_rate_alt);
                    var2 *= (1.0 - error_rate_ref - error_rate_alt);
                    var3 *= (1.0 - error_rate_ref - error_rate_alt);
                    var1 = var1*var1;
                    var2 = var2*var2;
                    var3 = var3*var3;
                    var1 *= contam_rate_var;
                    var2 *= contam_rate_var;
                    var3 *= contam_rate_var;
                    //var1 = var2 = var3 = contam_rate_var;
                    
                    //double varmu = (var1+ var2+var3)/3.0;
                    //var1 = var2 = var3 = varmu;

                    double fac1 = (exp1*(1.0-exp1))/var1 - 1.0;
                    double fac2 = (exp2*(1.0-exp2))/var2 - 1.0;
                    double fac3 = (exp3*(1.0-exp3))/var3 - 1.0;
                    double a1 = fac1*exp1;
                    double b1 = fac1*(1.0-exp1);
                    double a2 = fac2*exp2;
                    double b2 = fac2*(1.0-exp2);
                    double a3 = fac3*exp3;
                    double b3 = fac3*(1.0-exp3);
                    ll1 = dbetabin(alt, ref+alt, a1, b1);
                    ll2 = dbetabin(alt, ref+alt, a2, b2);
                    ll3 = dbetabin(alt, ref+alt, a3, b3);
                    
                }

                map<int, double> m;
                if (llrs.count(i) == 0){
                    llrs.insert(make_pair(i, m));
                }
                if (llrs.count(j) == 0){
                    llrs.insert(make_pair(j, m));
                }
                if (llrs[i].count(j) == 0){
                    llrs[i].insert(make_pair(j, 0.0));
                }
                llrs[i][j] += (ll1-ll2);
                if (doublet_rate > 0.0){
                    // Store comparisons between i and (i,j) combo and 
                    // between j and (i,j) combo
                    if (llrs[i].count(k) == 0){
                        llrs[i].insert(make_pair(k, 0.0));
                    }
                    if (llrs[j].count(k) == 0){
                        llrs[j].insert(make_pair(k, 0.0));
                    }
                    
                    llrs[i][k] += (ll1-ll3);
                    llrs[j][k] += (ll2-ll3);
                }
            }
        }
    }
    // Populate LLR table with singlet/singlet comparisons
    for (map<int, map<int, double> >::iterator x = llrs.begin(); x != llrs.end(); ++x){
        for (map<int, double>::iterator y = x->second.begin(); y != x->second.end(); ++y){
            if (y->first < n_samples){
                if (allowed_assignments.size() == 0 || 
                    (allowed_assignments.find(x->first) != allowed_assignments.end() &&
                     allowed_assignments.find(y->first) != allowed_assignments.end())){
                     
                     if (prior_weights != NULL && 
                        prior_weights->count(x->first) > 0 && prior_weights->count(y->first) > 0){
                        //y->second += log2((*prior_weights)[x->first]) - log2((*prior_weights)[y->first]);
                        y->second += (*prior_weights)[x->first] - (*prior_weights)[y->first];
                     }
                     tab.insert(x->first, y->first, y->second);
                }
            }
        }
    }
    if (doublet_rate > 0.0){ 
        
        // Toss out unlikely individuals (based on losing end of largest LLR
        // in the table, iteratively), until we are left with 10 individuals
        int n_target = 10;
        if (tab.n_indvs > n_target){
            bool success = tab.del(n_target);
            if (!success){
                //return false;
            }
        }
        
        // If we started with more than 3 individuals and threw out enough to get down
        // below 3, give up trying to make an assignment. This will only happen if there
        // are lots of ties, which would be the result of very sparse data.

        if (tab.n_indvs < 2 && tab.n_indvs < n_samples){
            return false;
        }
        
        // Get a list of all possible double identities to consider.
        vector<int> ks;
        for (int i = 0; i < n_samples-1; ++i){
            if (allowed_assignments.size() != 0 && 
                allowed_assignments.find(i) == allowed_assignments.end()){
                continue;
            }
            else if (!tab.included[i]){
                continue;
            }

            for (int j = i + 1; j < n_samples; ++j){
                if (allowed_assignments.size() != 0 && 
                    allowed_assignments.find(j) == allowed_assignments.end()){
                    continue;
                }
                else if (!tab.included[j]){
                    continue;
                }
                int k = hap_comb_to_idx(i, j, n_samples);
                if (allowed_assignments.size() == 0 || allowed_assignments.find(k) != 
                    allowed_assignments.end()){
                    
                    ks.push_back(k);
                }
            }
        }
        
        if (ks.size() > 0){
            // put k indices in increasing order so we know what order to store comparisons
            // in the data structure
            
            sort(ks.begin(), ks.end());
            // With all possible values of k to consider, we already have all member component vs 
            // double model comparisons computed. We now need to compute all other possible
            // single vs double model comparisons, as well as double model vs double model comparisons.
            
            compute_k_comps(llrs, tab, ks, n_samples, allowed_assignments, doublet_rate, prior_weights);
            
        }
    }
    
    // Disallow impossible combinations. We should already have excluded disallowed k combinations
    // so only need to exclude single individuals.
    if (allowed_assignments.size() > 0 || doublet_rate == 1.0){
        for (int i = 0; i < n_samples; ++i){
            if (doublet_rate == 1 || (allowed_assignments2.size() > 0 && 
                allowed_assignments2.find(i) == allowed_assignments2.end())){
                tab.disallow(i);
            }
        }
    }

    return true;
}
