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
#include "demux_parallel_llr.h"

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

llr_table::llr_table(int x){
    n_indvs = 0;
    int n_elt = x + (int)round(pow(2, binom_coef_log(x, 2)));
    included.reserve(n_elt);
    maxllr.reserve(n_elt);
}

llr_table::~llr_table(){
    lookup_llr.clear();
    included.clear();
    maxllr.clear();
    minllr.clear();
}

void llr_table::print(string& bc_str, vector<string>& samples){
    for (map<double, vector<pair<short, short> > >::iterator x = lookup_llr.begin();
        x != lookup_llr.end(); ++x){
        for (int i = 0; i < x->second.size(); ++i){
            if (included[x->second[i].first] && included[x->second[i].second]){
                string n1 = idx2name(x->second[i].first, samples);
                string n2 = idx2name(x->second[i].second, samples);
                fprintf(stdout, "%s\t%s\t%s\t%f\n", bc_str.c_str(),
                    n1.c_str(), n2.c_str(), -x->first);
                fprintf(stdout, "%s\t%s\t%s\t%f\n", bc_str.c_str(),
                    n2.c_str(), n1.c_str(), x->first);
            }
        }
    }
}

void llr_table::print_ranges(string& barcode, vector<string>& samples){
    int n_samples = samples.size();
    for (int i = 0; i < n_samples; ++i){
        if (included[i]){
            fprintf(stdout, "%s\t%s\t%f\t%f\n", barcode.c_str(), idx2name(i, samples).c_str(), 
                minllr[i], maxllr[i]);
        }
        for (int j = i + 1; j < n_samples; ++j){
            int k = hap_comb_to_idx(i, j, n_samples);
            fprintf(stdout, "%s\t%s\t%f\t%f\n", barcode.c_str(), idx2name(k, samples).c_str(),
                minllr[k], maxllr[k]);
        }
    }
}

void llr_table::insert(short i1, short i2, double llr){
    while(included.size() < i1+1){
        included.push_back(false);
        maxllr.push_back(0.0);
        minllr.push_back(0.0);
    }
    while(included.size() < i2+1){
        included.push_back(false);
        maxllr.push_back(0.0);
        minllr.push_back(0.0);
    }
    if (!included[i1]){
        ++n_indvs;
    }
    if (!included[i2]){
        ++n_indvs;
    }
    if (!included[i1] || maxllr[i1] < llr){
        maxllr[i1] = llr;
    }
    if (!included[i1] || minllr[i1] > llr){
        minllr[i1] = llr;
    }
    if (!included[i2] || maxllr[i2] < -llr){
        maxllr[i2] = -llr;
    }
    if (!included[i2] || minllr[i2] > -llr){
        minllr[i2] = -llr;
    }
    included[i1] = true;
    included[i2] = true;
    
    short low;
    short high;
    if (llr > 0){
        low = i2;
        high = i1;
        llr = -llr;
    }
    else{
        low = i1;
        high = i2;
    }
    if (lookup_llr.count(llr) == 0){
        vector<pair<short, short> > v;
        lookup_llr.emplace(llr, v);
    }
    lookup_llr[llr].push_back(make_pair(low, high));
}

void llr_table::disallow(short i){
    if (i < included.size()){
        if (included[i]){
            maxllr[i] = 0.0;
            minllr[i] = 0.0;
            --n_indvs;
        }
        included[i] = false;
    }
}

// Recalculate minllr/maxllr values considering only still-included identities
// This is needed after disallowing identities, since minllr may reflect
// comparisons against now-disallowed identities
void llr_table::recalculate_minmax(){
    // Reset all minllr/maxllr for included identities
    for (int i = 0; i < included.size(); ++i){
        if (included[i]){
            minllr[i] = std::numeric_limits<double>::max();
            maxllr[i] = std::numeric_limits<double>::lowest();
        }
    }
    
    // Recalculate from lookup_llr, considering only included pairs
    for (auto& entry : lookup_llr){
        double llr_neg = entry.first;  // Stored as negative
        for (auto& p : entry.second){
            short low = p.first;
            short high = p.second;
            
            if (!included[low] || !included[high]){
                continue;  // Skip pairs involving disallowed identities
            }
            
            // low has negative LLR vs high, high has positive LLR vs low
            // llr_neg is stored as negative, so:
            // LLR(low vs high) = llr_neg (negative)
            // LLR(high vs low) = -llr_neg (positive)
            
            if (llr_neg < minllr[low]){
                minllr[low] = llr_neg;
            }
            if (llr_neg > maxllr[low]){
                maxllr[low] = llr_neg;
            }
            if (-llr_neg < minllr[high]){
                minllr[high] = -llr_neg;
            }
            if (-llr_neg > maxllr[high]){
                maxllr[high] = -llr_neg;
            }
        }
    }
    
    // Handle identities that have no remaining comparisons
    for (int i = 0; i < included.size(); ++i){
        if (included[i] && minllr[i] == std::numeric_limits<double>::max()){
            minllr[i] = 0.0;
            maxllr[i] = 0.0;
        }
    }
}

bool llr_table::del(int n_keep){
    if (n_indvs < n_keep){
        return false;
    }
    
    it = lookup_llr.begin();
    while (n_indvs > n_keep && it != lookup_llr.end()){
        bool del_all = it->second.size() <= (n_indvs - n_keep);
        map<double, int> maxllr_counts;
        int indv_tot = 0;
        for (vector<pair<short, short> >::iterator x = it->second.begin();
            x != it->second.end();){
            if (!included[x->first] || !included[x->second]){
                it->second.erase(x);
            }
            else if (del_all){
                if (included[x->first]){
                    --n_indvs;
                    included[x->first] = false;
                    minllr[x->first] = 0.0;
                    maxllr[x->first] = 0.0;
                }
                it->second.erase(x);
            }
            else{
                double mllr = maxllr[x->first];
                if (maxllr_counts.count(mllr) == 0){
                    maxllr_counts.insert(make_pair(mllr, 1));
                }
                else{
                    maxllr_counts[mllr]++;
                }
                indv_tot++;
                ++x;
            }
        }
        
        if (indv_tot <= n_indvs - n_keep){
            del_all = true;
        }

        if (!del_all){
            double cutoff = 0.0;
            int runtot = 0;
            for (map<double, int>::iterator m = maxllr_counts.begin(); m != maxllr_counts.end();
                ++m){
                runtot += m->second;
                cutoff = m->first;
                if (runtot >= n_indvs - n_keep){
                    break;
                }
            }
            for (vector<pair<short, short> >::iterator x = it->second.begin(); x != it->second.end(); ){
                if (maxllr[x->first] <= cutoff){
                    if (included[x->first]){
                        included[x->first] = false;
                        maxllr[x->first] = 0.0;
                        minllr[x->first] = 0.0;
                        n_indvs--;
                    }
                    it->second.erase(x);
                }
                else{
                    ++x;
                }
            }
        }
        else{
            for (vector<pair<short, short> >::iterator x = it->second.begin(); x != 
                it->second.end(); ){
                if (included[x->first]){
                    included[x->first] = false;
                    n_indvs--;
                    maxllr[x->first] = 0.0;
                    minllr[x->first] = 0.0;
                }
                it->second.erase(x);
            }
            lookup_llr.erase(it++);
        }
    }
    return true;
}

void llr_table::get_max(int& best_idx, double& best_llr){
    best_idx = -1;
    best_llr = 0.0;
    
    for (int i = 0; i < included.size(); ++i){
        if (included[i]){
            if (best_idx == -1 || minllr[i] > best_llr){
                best_idx = i;
                best_llr = minllr[i];
            }
        }
    }
}

// NEW: Get min_margin for a specific identity
double llr_table::get_min_margin(int identity) const {
    if (identity < 0 || identity >= (int)minllr.size() || !included[identity]){
        return 0.0;
    }
    return minllr[identity];
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

void get_kcomps_cond(map<int, map<int, double> >& kcomps,
    map<int, map<int, double> >& llrs,
    llr_table& tab,
    int n_samples,
    set<int>& allowed_assignments){

    for (map<int, map<int, double> >::iterator llr = llrs.begin(); llr != llrs.end();
        ++llr){

        vector<int> js;
        vector<int> comps;
        for (map<int, double>::iterator llr2 = llr->second.begin(); llr2 != llr->second.end();
            ++llr2){
            if (llr2->first >= n_samples){
                if (!tab.included[llr2->first]){
                }

                if (true){
                    pair<int, int> comb = idx_to_hap_comb(llr2->first, n_samples);
                    int j = -1;
                    if (comb.first == llr->first){
                        j = comb.second;
                    }
                    else{
                        j = comb.first;
                    }
                    js.push_back(j);
                    comps.push_back(-llr2->second);
                }
            }
        }
        
        if (comps.size() > 1){
            for (int i = 0; i < comps.size()-1; ++i){
                int idx1 = js[i];
                for (int j = i +1; j < comps.size(); ++j){
                    int idx2 = js[j];
                    double llr = comps[i] - comps[j];
                    if (idx1 < idx2){
                        if (kcomps.count(idx1) == 0){
                            map<int, double> m;
                            kcomps.insert(make_pair(idx1, m));
                        }
                        if (kcomps[idx1].count(idx2) == 0){
                            kcomps[idx1].insert(make_pair(idx2, 0.0));
                        }
                        kcomps[idx1][idx2] += llr;
                    }
                    else{
                        if (kcomps.count(idx2) == 0){
                            map<int, double> m;
                            kcomps.insert(make_pair(idx2, m));
                        }
                        if (kcomps[idx2].count(idx1) == 0){
                            kcomps[idx2].insert(make_pair(idx1, 0.0));
                        }
                        kcomps[idx2][idx1] += llr;
                    }
                }
            }
        }
    }
}

void compute_k_comps(map<int, map<int, double> >& llrs, 
    llr_table& tab,
    vector<int>& ks,
    int n_samples,
    set<int>& allowed_assignments,
    double doublet_rate,
    map<int, double>* prior_weights){
    
    // Get a table of LLRs between different combination members, conditional on a cell
    // being half another member. In other words, keys/values here are
    // ID1 -> ID2A -> ID2B -> LLR
    // LLR is the log likelihood ratio of ID2A vs ID2B being included in a double model
    // with ID1. 
    map<int, map<int, map<int, double> > > kcomp_cond;
    for (map<int, map<int, double> >::iterator samp = llrs.begin(); samp != llrs.end(); ++samp){
        if (!tab.included[samp->first]){
            continue;
        }
        map<int, map<int, double> > llrstmp;
        llrstmp.insert(make_pair(samp->first, samp->second));
        map<int, map<int, double> > kcomp; 
        get_kcomps_cond(kcomp, llrstmp, tab, n_samples, allowed_assignments);
        kcomp_cond.insert(make_pair(samp->first, kcomp));
    }
    // compute all single vs double comparisons
    for (int i = 0; i < n_samples; ++i){
        if (allowed_assignments.size() != 0 && 
            allowed_assignments.find(i) == allowed_assignments.end()){
            continue;
        }
        else if (!tab.included[i]){
            continue;
        }

        for (int ki = 0; ki < ks.size(); ++ki){
            // NOTE: ks is already a filtered list (allowable)
            int k = ks[ki];
            pair<int, int> comb = idx_to_hap_comb(k, n_samples);
            if (comb.first != i && comb.second != i){
                double ll1 = llrs[comb.first][k];
                double ll2 = llrs[comb.second][k];
                if (comb.first < i){
                    ll1 += -llrs[comb.first][i];
                }
                else{
                    ll1 += llrs[i][comb.first];
                }
                if (comb.second < i){
                    ll2 += -llrs[comb.second][i];
                }
                else{
                    ll2 += llrs[i][comb.second];
                }
                double llr = 0.5*ll1 + 0.5*ll2;
                if (prior_weights != NULL && 
                    prior_weights->count(i) > 0 && prior_weights->count(k) > 0){
                    llr += (*prior_weights)[i] - (*prior_weights)[k];
                }
                else if (doublet_rate != 0.5 && doublet_rate < 1.0){
                    // i is singlet, k is doublet
                    llr += log2(1.0 - doublet_rate) - log2(doublet_rate);
                }
                tab.insert(i, k, llr);
            }
            else{
                // It's already been computed.
                double llr = llrs[i][k];
                if (prior_weights != NULL && 
                    prior_weights->count(i) > 0 && prior_weights->count(k) > 0){
                    llr += (*prior_weights)[i] - (*prior_weights)[k];
                }
                else if (doublet_rate != 0.5 && doublet_rate < 1.0){
                    llr += log2(1.0 - doublet_rate) - log2(doublet_rate);
                }
                tab.insert(i, k, llr);
            }
        }
    }

    // compute all double vs double comparisons
    for (int ki = 0; ki < ks.size()-1; ++ki){
        int k1 = ks[ki];
        
        pair<int, int> comb1 = idx_to_hap_comb(k1, n_samples);
        for (int kj = ki+1; kj < ks.size(); ++kj){
            int k2 = ks[kj];
            pair<int, int> comb2 = idx_to_hap_comb(k2, n_samples);
            vector<double> llr_parts;
            int a = comb1.first;
            int b = comb1.second;
            int c = comb2.first;
            int d = comb2.second;

            if (a == c){
                if (b < d){
                    llr_parts.push_back(kcomp_cond[a][b][d]);
                }
                else{
                    llr_parts.push_back(-kcomp_cond[a][d][b]);
                }
            }
            else if (a == d){
                if (b < c){
                    llr_parts.push_back(kcomp_cond[a][b][c]);
                }
                else{
                    llr_parts.push_back(-kcomp_cond[a][c][b]);
                }
            }
            else if (b == c){
                if (a < d){
                    llr_parts.push_back(kcomp_cond[b][a][d]);
                }         
                else{
                    llr_parts.push_back(-kcomp_cond[b][d][a]);
                }
            }
            else if (b == d){
                if (a < c){
                    llr_parts.push_back(kcomp_cond[b][a][c]);
                }
                else{
                    llr_parts.push_back(-kcomp_cond[b][c][a]);
                }
            }
            else{
                double llr1 = 0;
                double llr2 = 0;
                double llr3 = 0;
                double llr4 = 0;
                double llr5 = 0;
                double llr6 = 0;
                double llr7 = 0;
                double llr8 = 0;
                // Calculate (A+B)/(C+D) 
                if (b < c){
                    llr1 += kcomp_cond[a][b][c];
                }
                else{
                    llr1 += -kcomp_cond[a][c][b];
                }
                if (a < d){
                    llr1 += kcomp_cond[c][a][d];
                }
                else{
                    llr1 += -kcomp_cond[c][d][a];
                }
                
                if (a < c){
                    llr2 += kcomp_cond[b][a][c];
                }
                else{
                    llr2 += -kcomp_cond[b][c][a];
                }
                if (b < d){
                    llr2 += kcomp_cond[c][b][d];
                }
                else{
                    llr2 += -kcomp_cond[c][d][b];
                }

                if (b < d){
                    llr3 += kcomp_cond[a][b][d];
                }
                else{
                    llr3 += -kcomp_cond[a][d][b];
                }
                if (a < c){
                    llr3 += kcomp_cond[d][a][c];
                }
                else{
                    llr3 += -kcomp_cond[d][c][a];
                }

                if (a < d){
                    llr4 += kcomp_cond[b][a][d];
                }
                else{
                    llr4 += -kcomp_cond[b][d][a];
                }
                if (b < c){
                    llr4 += kcomp_cond[d][b][c];
                }
                else{
                    llr4 += -kcomp_cond[d][c][b];
                }

                if (d < a){
                    llr5 += kcomp_cond[c][d][a];
                }
                else{
                    llr5 += -kcomp_cond[c][a][d];
                }
                if (c < b){
                    llr5 += kcomp_cond[a][c][b];
                }
                else{
                    llr5 += -kcomp_cond[a][b][c];
                }
                llr5 = -llr5;

                if (c < a){
                    llr6 += kcomp_cond[d][c][a];
                }
                else{
                    llr6 += -kcomp_cond[d][a][c];
                }
                if (d < b){
                    llr6 += kcomp_cond[a][d][b];
                }
                else{
                    llr6 += -kcomp_cond[a][b][d];
                }
                llr6 = -llr6;

                if (d < b){
                    llr7 += kcomp_cond[c][d][b];
                }
                else{
                    llr7 += -kcomp_cond[c][b][d];
                }
                if (c < a){
                    llr7 += kcomp_cond[b][c][a];
                }
                else{
                    llr7 += -kcomp_cond[b][a][c];
                }
                llr7 = -llr7;

                if (c < b){
                    llr8 += kcomp_cond[d][c][b];
                }
                else{
                    llr8 += -kcomp_cond[d][b][c];
                }
                if (d < a){
                    llr8 += kcomp_cond[b][d][a];
                }
                else{
                    llr8 += -kcomp_cond[b][a][d];
                }
                llr8 = -llr8;
                llr_parts.push_back(llr1);
                llr_parts.push_back(llr2);
                llr_parts.push_back(llr3);
                llr_parts.push_back(llr4);
                llr_parts.push_back(llr5);
                llr_parts.push_back(llr6);
                llr_parts.push_back(llr7);
                llr_parts.push_back(llr8);
            }
            
            float llrsum = 0.0;
            float llrcount = 0;
            for (int idx = 0; idx < llr_parts.size(); ++idx){
                llrsum += llr_parts[idx];
                llrcount++;
            }
            double llr = llrsum / llrcount;
            if (prior_weights != NULL && 
                prior_weights->count(k1) > 0 && prior_weights->count(k2) > 0){
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
    vector<RunnerUp>& runners){
    
    runners.clear();
    
    // Get min_margin for the winner (already computed in tab)
    diag.min_margin = tab.get_min_margin(winner);
    
    // Find the worst competitor (identity that gave min_margin to winner)
    // We need to look through all LLRs involving the winner
    double min_llr_found = std::numeric_limits<double>::max();
    int worst_comp = -1;
    
    // Check all pairwise comparisons involving the winner
    if (llrs.count(winner) > 0){
        for (const auto& kv : llrs.at(winner)){
            int other = kv.first;
            if (!tab.included[other]) continue;
            
            // LLR(winner vs other) = llrs[winner][other] if winner < other in storage
            // Need to get the direct LLR from winner's perspective
            double direct_llr = kv.second;
            if (direct_llr < min_llr_found){
                min_llr_found = direct_llr;
                worst_comp = other;
            }
        }
    }
    
    // Also check cases where winner is in the "other" position
    for (const auto& outer : llrs){
        if (outer.first == winner) continue;
        if (!tab.included[outer.first]) continue;
        
        if (outer.second.count(winner) > 0){
            // llrs[outer.first][winner] is from outer.first's perspective
            // So from winner's perspective, it's -llrs[outer.first][winner]
            double llr_from_winner = -outer.second.at(winner);
            if (llr_from_winner < min_llr_found){
                min_llr_found = llr_from_winner;
                worst_comp = outer.first;
            }
        }
    }
    
    diag.worst_competitor = worst_comp;
    
    // Count n_close: identities within close_threshold of winner
    diag.n_close = 0;
    
    // Build list of all identity scores for runner-up extraction
    vector<pair<double, int> > identity_scores;  // (min_margin, identity)
    
    const vector<double>& all_minllr = tab.get_minllr();
    
    for (size_t i = 0; i < tab.included.size(); ++i){
        if (!tab.included[i]) continue;
        if ((int)i == winner) continue;
        
        double score = all_minllr[i];
        identity_scores.push_back(make_pair(score, (int)i));
        
        // Check if within close_threshold of winner
        // A competitor is "close" if winner's margin over them is small
        // winner_llr is the winner's min_margin
        // If (winner_llr - score) < close_threshold, they're close
        if ((diag.min_margin - score) < close_threshold){
            diag.n_close++;
        }
    }
    
    // Sort by score descending to get best runner-ups
    sort(identity_scores.begin(), identity_scores.end(), 
         [](const pair<double, int>& a, const pair<double, int>& b){
             return a.first > b.first;  // Descending
         });
    
    // Extract top n_runner_ups
    int n_to_extract = min((int)identity_scores.size(), n_runner_ups);
    for (int i = 0; i < n_to_extract; ++i){
        int runner_id = identity_scores[i].second;
        double runner_margin = identity_scores[i].first;
        
        // Get direct LLR vs winner
        double llr_vs_winner = 0.0;
        
        // Look up in llrs map
        if (winner < runner_id){
            if (llrs.count(winner) > 0 && llrs.at(winner).count(runner_id) > 0){
                llr_vs_winner = -llrs.at(winner).at(runner_id);  // Negative because runner lost
            }
        }
        else{
            if (llrs.count(runner_id) > 0 && llrs.at(runner_id).count(winner) > 0){
                llr_vs_winner = llrs.at(runner_id).at(winner);  // Already from runner's perspective
            }
        }
        
        runners.push_back(RunnerUp(runner_id, llr_vs_winner, runner_margin));
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
        diag.het_balance_var = -1.0;  // Insufficient data
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
        diag.het_balance_var = -1.0;
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
    
    bool is_doublet = (assigned_id >= n_samples);
    
    if (is_doublet) {
        auto combo = idx_to_hap_comb(assigned_id, n_samples);
        int indv1 = combo.first;
        int indv2 = combo.second;
        
        const WelfordStats& ws1 = welford_stats.get(indv1);
        const WelfordStats& ws2 = welford_stats.get(indv2);
        
        double var1 = ws1.variance(min_sites);
        double var2 = ws2.variance(min_sites);
        
        if (var1 < 0 && var2 < 0) {
            diag.het_balance_var = -1.0;
            diag.n_het_sites = 0;
            diag.het_total_depth = 0;
            return;
        }
        
        if (var1 < 0) {
            diag.het_balance_var = var2;
            diag.n_het_sites = (int)ws2.n;
            diag.het_total_depth = ws2.total_depth;
        } else if (var2 < 0) {
            diag.het_balance_var = var1;
            diag.n_het_sites = (int)ws1.n;
            diag.het_total_depth = ws1.total_depth;
        } else {
            double n_total = ws1.n + ws2.n;
            diag.het_balance_var = (var1 * ws1.n + var2 * ws2.n) / n_total;
            diag.n_het_sites = (int)n_total;
            diag.het_total_depth = ws1.total_depth + ws2.total_depth;
        }
    } else {
        const WelfordStats& ws = welford_stats.get(assigned_id);
        diag.het_balance_var = ws.variance(min_sites);
        diag.n_het_sites = (int)ws.n;
        diag.het_total_depth = ws.total_depth;
    }
}

// ============================================================================
// PARALLEL IDENTITY ASSIGNMENT
// ============================================================================

void assign_ids_parallel(
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
    int n_target){
    
    assignments.clear();
    assignments_llr.clear();
    
    if (g_verbose){
        fprintf(stderr, "[VERBOSE] assign_ids_parallel: n_samples=%d, n_target=%d, doublet_rate=%.3f\n",
            (int)samples.size(), n_target, doublet_rate);
        fprintf(stderr, "[VERBOSE]   allowed_assignments: %lu, allowed_assignments2: %lu\n",
            allowed_assignments.size(), allowed_assignments2.size());
    }
    
    // Collect barcodes into vector for parallel processing
    vector<unsigned long> barcodes;
    barcodes.reserve(cell_counts.size());
    for (auto& kv : cell_counts){
        barcodes.push_back(kv.first);
    }
    
    // Pre-allocate results
    vector<int> result_assn(barcodes.size(), -1);
    vector<double> result_llr(barcodes.size(), 0.0);
    
    fprintf(stderr, "Assigning identities to %lu cells using %d threads...\n",
        barcodes.size(), n_threads);
    
    int n_samples = samples.size();
    
    #pragma omp parallel for num_threads(n_threads) schedule(dynamic, 100)
    for (size_t idx = 0; idx < barcodes.size(); ++idx){
        unsigned long bc = barcodes[idx];
        const CellCounts& counts = cell_counts[bc];
        
        // Skip empty cells
        if (counts.is_empty()) continue;
        
        map<int, map<int, double> > llrs;
        llr_table tab(n_samples);
        
        map<int, double>* pw_ptr = use_prior_weights ? &prior_weights : NULL;
        
        bool success = populate_llr_table_optimized(
            counts, llrs, tab, n_samples,
            allowed_assignments, allowed_assignments2,
            doublet_rate, error_rate_ref, error_rate_alt,
            pw_ptr, false, 0.0, 0.0, NULL, n_target);
        
        if (success){
            int assn = -1;
            double llr_final = 0.0;
            tab.get_max(assn, llr_final);
            
            if (llr_final > 0.0){
                result_assn[idx] = assn;
                result_llr[idx] = llr_final;
            }
        }
    }
    
    // Collect results
    for (size_t idx = 0; idx < barcodes.size(); ++idx){
        if (result_assn[idx] >= 0){
            assignments.emplace(barcodes[idx], result_assn[idx]);
            assignments_llr.emplace(barcodes[idx], result_llr[idx]);
        }
    }
    
    fprintf(stderr, "Assigned %lu cells\n", assignments.size());
    
    if (g_verbose){
        // Count singlets vs doublets
        int n_samples = samples.size();
        int n_singlets = 0, n_doublets = 0;
        map<int, int> assignment_counts;
        for (auto& a : assignments){
            assignment_counts[a.second]++;
            if (a.second < n_samples) n_singlets++;
            else n_doublets++;
        }
        fprintf(stderr, "[VERBOSE] Assignment summary: %d singlets, %d doublets\n", 
            n_singlets, n_doublets);
        for (auto& ac : assignment_counts){
            if (ac.first < n_samples){
                fprintf(stderr, "[VERBOSE]   %s: %d cells\n", samples[ac.first].c_str(), ac.second);
            }
            else{
                pair<int, int> combo = idx_to_hap_comb(ac.first, n_samples);
                fprintf(stderr, "[VERBOSE]   %s+%s: %d cells\n", 
                    samples[combo.first].c_str(), samples[combo.second].c_str(), ac.second);
            }
        }
    }
}

// ============================================================================
// NEW: PARALLEL IDENTITY ASSIGNMENT WITH DIAGNOSTICS
// ============================================================================

void assign_ids_parallel_with_diagnostics(
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
    // Diagnostic options
    bool compute_diagnostics,
    int n_runner_ups,
    double close_threshold,
    robin_hood::unordered_map<unsigned long, CellCounts>* het_counts,
    // Diagnostic outputs
    robin_hood::unordered_map<unsigned long, CellDiagnostics>& diagnostics,
    robin_hood::unordered_map<unsigned long, vector<RunnerUp> >& runner_ups){
    
    assignments.clear();
    assignments_llr.clear();
    diagnostics.clear();
    runner_ups.clear();
    
    if (g_verbose){
        fprintf(stderr, "[VERBOSE] assign_ids_parallel_with_diagnostics: n_samples=%d, n_target=%d\n",
            (int)samples.size(), n_target);
        fprintf(stderr, "[VERBOSE]   compute_diagnostics=%d, n_runner_ups=%d, close_threshold=%.1f\n",
            compute_diagnostics, n_runner_ups, close_threshold);
        fprintf(stderr, "[VERBOSE]   het_counts provided: %s\n", het_counts ? "yes" : "no");
    }
    
    // Collect barcodes into vector for parallel processing
    vector<unsigned long> barcodes;
    barcodes.reserve(cell_counts.size());
    for (auto& kv : cell_counts){
        barcodes.push_back(kv.first);
    }
    
    // Pre-allocate results
    vector<int> result_assn(barcodes.size(), -1);
    vector<double> result_llr(barcodes.size(), 0.0);
    vector<CellDiagnostics> result_diag(barcodes.size());
    vector<vector<RunnerUp> > result_runners(barcodes.size());
    
    fprintf(stderr, "Assigning identities to %lu cells using %d threads%s...\n",
        barcodes.size(), n_threads, 
        compute_diagnostics ? " (with diagnostics)" : "");
    
    int n_samples = samples.size();
    
    #pragma omp parallel for num_threads(n_threads) schedule(dynamic, 100)
    for (size_t idx = 0; idx < barcodes.size(); ++idx){
        unsigned long bc = barcodes[idx];
        const CellCounts& counts = cell_counts[bc];
        
        // Skip empty cells
        if (counts.is_empty()) continue;
        
        map<int, map<int, double> > llrs;
        llr_table tab(n_samples);
        
        map<int, double>* pw_ptr = use_prior_weights ? &prior_weights : NULL;
        
        bool success = populate_llr_table_optimized(
            counts, llrs, tab, n_samples,
            allowed_assignments, allowed_assignments2,
            doublet_rate, error_rate_ref, error_rate_alt,
            pw_ptr, false, 0.0, 0.0, NULL, n_target);
        
        if (success){
            int assn = -1;
            double llr_final = 0.0;
            tab.get_max(assn, llr_final);
            
            if (llr_final > 0.0){
                result_assn[idx] = assn;
                result_llr[idx] = llr_final;
                
                // Extract diagnostics BEFORE table destruction
                if (compute_diagnostics){
                    CellDiagnostics diag;
                    vector<RunnerUp> runners;
                    
                    get_diagnostics_from_llrs(llrs, tab, assn, llr_final, n_samples,
                        n_runner_ups, close_threshold, diag, runners);
                    
                    // Compute total depth from main counts
                    diag.total_depth = compute_total_depth(counts, n_samples);
                    
                    // Compute het balance stats if het counts provided
                    if (het_counts != NULL){
                        auto het_it = het_counts->find(bc);
                        if (het_it != het_counts->end()){
                            compute_het_balance_stats(het_it->second, assn, n_samples, diag);
                        }
                    }
                    
                    result_diag[idx] = diag;
                    result_runners[idx] = runners;
                }
            }
        }
    }
    
    // Collect results
    for (size_t idx = 0; idx < barcodes.size(); ++idx){
        if (result_assn[idx] >= 0){
            unsigned long bc = barcodes[idx];
            assignments.emplace(bc, result_assn[idx]);
            assignments_llr.emplace(bc, result_llr[idx]);
            
            if (compute_diagnostics){
                diagnostics.emplace(bc, result_diag[idx]);
                runner_ups.emplace(bc, result_runners[idx]);
            }
        }
    }
    
    fprintf(stderr, "Assigned %lu cells\n", assignments.size());
    if (compute_diagnostics){
        fprintf(stderr, "Collected diagnostics for %lu cells\n", diagnostics.size());
    }
}

/**
 * Extended assignment with Welford het balance method
 */
void assign_ids_parallel_with_diagnostics_extended(
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
    const vector<pair<int, int>>* idx_to_site,
    HetBalanceMethod het_method,
    int min_het_sites,
    double min_het_depth,
    robin_hood::unordered_map<unsigned long, CellDiagnostics>& diagnostics,
    robin_hood::unordered_map<unsigned long, vector<RunnerUp> >& runner_ups) {
    
    int n_samples = samples.size();
    const char* method_name = (het_method == HetBalanceMethod::PERSITE) ? "per-site" : "Welford";
    
    if (g_verbose) {
        fprintf(stderr, "[VERBOSE] assign_ids_parallel_with_diagnostics_extended:\n");
        fprintf(stderr, "[VERBOSE]   cells: %lu, method: %s, min_het_sites: %d\n", 
                cell_counts.size(), method_name, min_het_sites);
    }
    
    // Convert to vector for parallel processing
    vector<unsigned long> barcodes;
    barcodes.reserve(cell_counts.size());
    for (const auto& kv : cell_counts) {
        barcodes.push_back(kv.first);
    }
    
    // Results storage
    vector<int> result_assn(barcodes.size(), -1);
    vector<double> result_llr(barcodes.size(), 0.0);
    vector<CellDiagnostics> result_diag(barcodes.size());
    vector<vector<RunnerUp> > result_runners(barcodes.size());
    
    atomic<int> cells_done(0);
    int total_cells = barcodes.size();
    
    #pragma omp parallel num_threads(n_threads)
    {
        map<int, map<int, double> > llrs;
        
        #pragma omp for schedule(dynamic, 100)
        for (size_t idx = 0; idx < barcodes.size(); idx++) {
            unsigned long bc = barcodes[idx];
            const CellCounts& counts = cell_counts[bc];
            
            if (counts.is_empty()) continue;
            
            int n_identities = n_samples + n_samples * (n_samples - 1) / 2;
            llr_table tab(n_identities);
            llrs.clear();
            
            bool success = populate_llr_table_optimized(
                counts, llrs, tab, n_samples,
                allowed_assignments, allowed_assignments2,
                doublet_rate, error_rate_ref, error_rate_alt,
                use_prior_weights ? &prior_weights : NULL,
                false, 0.0, 0.0, NULL, n_target);
            
            if (!success) continue;
            
            int best_idx;
            double best_llr;
            tab.get_max(best_idx, best_llr);
            
            result_assn[idx] = best_idx;
            result_llr[idx] = best_llr;
            
            if (compute_diagnostics) {
                CellDiagnostics diag;
                vector<RunnerUp> runners;
                
                get_diagnostics_from_llrs(llrs, tab, best_idx, best_llr,
                    n_samples, n_runner_ups, close_threshold, diag, runners);
                
                diag.total_depth = compute_total_depth(counts, n_samples);
                
                // Compute het balance using selected method
                if (het_data != NULL) {
                    auto het_it = het_data->find(bc);
                    if (het_it != het_data->end()) {
                        const CellHetData& cell_het = het_it->second;
                        
                        if (het_method == HetBalanceMethod::PERSITE && 
                            het_snpdat != NULL && idx_to_site != NULL) {
                            compute_het_balance_persite(
                                cell_het.persite_data,
                                *het_snpdat,
                                *idx_to_site,
                                best_idx,
                                n_samples,
                                min_het_depth,
                                min_het_sites,
                                diag);
                        } else {
                            compute_het_balance_welford(
                                cell_het.welford_stats,
                                best_idx,
                                n_samples,
                                min_het_sites,
                                diag);
                        }
                    }
                }
                
                result_diag[idx] = diag;
                result_runners[idx] = runners;
            }
            
            int done = ++cells_done;
            if (done % 1000 == 0 || done == total_cells) {
                fprintf(stderr, "\rAssigning: %d/%d cells    ", done, total_cells);
            }
        }
    }
    
    // Collect results
    for (size_t idx = 0; idx < barcodes.size(); ++idx) {
        if (result_assn[idx] >= 0) {
            unsigned long bc = barcodes[idx];
            assignments.emplace(bc, result_assn[idx]);
            assignments_llr.emplace(bc, result_llr[idx]);
            
            if (compute_diagnostics) {
                diagnostics.emplace(bc, result_diag[idx]);
                runner_ups.emplace(bc, result_runners[idx]);
            }
        }
    }
    
    fprintf(stderr, "\nAssigned %lu cells (%s het method)\n", assignments.size(), method_name);
    if (compute_diagnostics) {
        fprintf(stderr, "Collected diagnostics for %lu cells\n", diagnostics.size());
    }
}
