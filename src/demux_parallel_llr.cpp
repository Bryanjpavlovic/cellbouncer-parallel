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
