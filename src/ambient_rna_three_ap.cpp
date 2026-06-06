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
#include <cfloat>
#include <random>
#include <mixtureDist/functions.h>
#include <optimML/brent.h>
#include <optimML/multivar_ml.h>
#include <htswrapper/robin_hood/robin_hood.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "common.h"
#include "demux_vcf_llr.h"
#include "ambient_rna_three_ap.h"

using std::cout;
using std::endl;
using namespace std;

// ============================================================================
// ambient_rna_three_ap.cpp
// Three-component contamination model for tetraploid cells
// with r-feedback and adaptive prior extensions (ambient profile refinements)
//
// ============================================================================

const string AMBIENT_RNA_VERSION = "2.1-ns";
const string AMBIENT_RNA_VERSION_MSG = "Three-component model with r-feedback and adaptive prior";

/**
 * Three-component contamination model for single-cell RNA data with tetraploid
 * cell fusions.
 *
 * DIPLOID SINGLETS (standard two-component):
 *   p = (1-c)*p_e + c*p_c
 *   where p_e = nalt/2, optimized via Brent on c
 *
 * HETEROTYPIC TETRAPLOIDS (three-component):
 *   p = (1-c) * [r * p_A + (1-r) * p_B] + c * p_c
 *   where p_A = nalt_A/2, p_B = nalt_B/2
 *   Jointly optimized via BFGS on (r, c)
 *
 * HOMOTYPIC TETRAPLOIDS (A+A):
 *   p_A = p_B, so r drops out. Treated as diploid (two-component).
 *
 * Error rate adjustment:
 *   p_e_adjusted = p_e - p_e*e_a + (1-p_e)*e_r
 *   Applied identically to p_A, p_B, and p_e for consistency.
 */

// Define a value to bump 0 and 1 binomial expectations away from these
// numbers
double binom_0 = 1e-6;

// ===== Utility functions ===== 

/**
 * Integrate binomial log likelihood wrt p_c (expected alt allele match rate from contamination)
 * 
 * * NOT USED *
 */
double integral_ll_pc(double n, double k, double c, double p_e, double p_c){
    double x = c*p_c*(binom_coef_log(n,k)/log2(exp(1)) - n);
    double y = (k-n)*((c-1)*p_e - c*p_c + 1);
    double z = log((c-1)*p_e - c*p_c + 1);
    double zz = k*(-c*p_e + c*p_c + p_e)*log(-c*p_e + c*p_c + p_e);
    return (1/c)*(x + y*z + zz);
}

/**
 * Derivative wrt contam rate (c) of the integral of binomial log likelihood wrt p_c
 * (expected alt allele match rate from contamination)
 *
 * This allows us to find the likeliest c, over the range of all possible p_c.
 * 
 * * NOT USED *
 */
double d_integral_ll_pc_dc(double n, double k, double c, double p_e, double p_c){
    double x = (p_e-1)*(k-n)*log((c-1)*p_e - c*p_c + 1);
    double y = k*p_e*log(-c*p_e + c*p_c + p_e);
    double z = c*n*(p_c-p_e);
    return (x - y + z)/(c*c);
}

/**
 * * NOT USED *
 */
double integral_ll_c(double n, double k, double c, double p_e, double p_c){
    double denom = p_e - p_c;
    double x = -c*(p_e - p_c)*(n - binom_coef_log(n,k)/log2(exp(1)));
    double y = (k-n)*((c-1)*p_e - c*p_c + 1)*log((c-1)*p_e - c*p_c + 1);
    double z = k*((c-1)*p_e - c*p_c)*log(-c*p_e + c*p_c + p_e);
    
    return (x - y + z)/denom;
}

/**
 * * NOT USED *
 */
double d_integral_ll_c_dpc(double n, double k, double c, double p_e, double p_c){
    double denom = pow(p_e-p_c, 2);
    double x = (k-n)*(c*(p_e-p_c) + (p_e-1)*log((c-1)*p_e - c*p_c + 1));
    double y = -(k*(c*(p_e - p_c) + p_e*log(-c*p_e + c*p_c + p_e)));
    
    return (x+y)/denom;
}

// ===== Functions to use for Brent's root finding method - for maximizing
// ===== univariate log likelihood functions

/**
 * Log likelihood when estimating c
 */
double ll_c(double c, const map<string, double >& data_d, 
    const map<string, int >& data_i){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    double p_c = data_d.at("p_c");
    double binom_p = (1.0 - c)*p_e + c*p_c;
    return logbinom(n, k, binom_p);
}

/**
 * Derivative of log likelihood wrt c
 */
double dll_dc(double c, const map<string, double >& data_d, 
    const map<string, int >& data_i){

    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    double p_c = data_d.at("p_c");
    double binom_p = (1.0 - c)*p_e + c*p_c;
    double dy_dp = (k-n*binom_p)/(binom_p - binom_p * binom_p);
    double dp_dc = p_c - p_e;
    return dy_dp * dp_dc;
}

/**
 * Second derivative of log likelihood wrt c
 */
double d2ll_dc2(double c, const map<string, double >& data_d, 
    const map<string, int>& data_i){

    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    double p_c = data_d.at("p_c");
    double binom_p = (1.0 - c)*p_e + c*p_c;
    double dp_dc = p_c - p_e;
    double d2y_dp2 = (k*(2*binom_p - 1) - n*binom_p*binom_p)/(pow(binom_p-1,2)*binom_p*binom_p);
    return d2y_dp2 * dp_dc * dp_dc;
}

// ===== Three-component model functions for BFGS (r, c) optimization =====
// params[0] = r (genome A fraction of endogenous RNA)
// params[1] = c (contamination fraction)
// data: n, k, p_A, p_B, p_c

/**
 * Log likelihood for three-component model.
 * p = (1-c) * [r * p_A + (1-r) * p_B] + c * p_c
 */
double ll_three(const vector<double>& params, const map<string, double>& data_d,
    const map<string, int>& data_i){

    double r = params[0];
    double c = params[1];
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_A = data_d.at("p_A");
    double p_B = data_d.at("p_B");
    double p_c = data_d.at("p_c");

    double p_endo = r * p_A + (1.0 - r) * p_B;
    double binom_p = (1.0 - c) * p_endo + c * p_c;

    // Clamp away from 0 and 1
    if (binom_p < 1e-6) binom_p = 1e-6;
    if (binom_p > 1.0 - 1e-6) binom_p = 1.0 - 1e-6;

    return logbinom(n, k, binom_p);
}

/**
 * Gradient of three-component log likelihood wrt (r, c).
 * 
 * dp/dr = (1-c) * (p_A - p_B)
 * dp/dc = p_c - [r * p_A + (1-r) * p_B]  = p_c - p_endo
 * dll/dp = (k - n*p) / (p - p^2)
 * dll/dr = dll/dp * dp/dr
 * dll/dc = dll/dp * dp/dc
 */
void dll_three(const vector<double>& params, const map<string, double>& data_d,
    const map<string, int>& data_i, vector<double>& results){

    double r = params[0];
    double c = params[1];
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_A = data_d.at("p_A");
    double p_B = data_d.at("p_B");
    double p_c = data_d.at("p_c");

    double p_endo = r * p_A + (1.0 - r) * p_B;
    double binom_p = (1.0 - c) * p_endo + c * p_c;

    if (binom_p < 1e-6) binom_p = 1e-6;
    if (binom_p > 1.0 - 1e-6) binom_p = 1.0 - 1e-6;

    double dll_dp = (k - n * binom_p) / (binom_p - binom_p * binom_p);

    double dp_dr = (1.0 - c) * (p_A - p_B);
    double dp_dc = p_c - p_endo;

    results[0] += dll_dp * dp_dr;
    results[1] += dll_dp * dp_dc;
}

// ===== Functions to use for multivar_ml_solver -- for optimizing
// ===== multivariate log likelihood functions using BFGS

/**
 * Log likelihood for estimating best set of p_c parameters
 * in ambient RNA
 */
double ll_ambmu(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    int amb_idx = data_i.at("ef_idx");
    double p_c = params[amb_idx];
    double c = data_d.at("c");
    double binom_p = (1.0-c) * p_e + c*p_c;
    if (binom_p < DBL_MIN*1e6){
        binom_p = DBL_MIN*1e6;
    }
    else if (binom_p > 1.0 - DBL_MIN*1e6){
        binom_p = 1.0 - DBL_MIN*1e6;
    }
    return logbinom(n, k, binom_p); 
}

/**
 * First derivative of log likelihood for estimating best set
 * of p_c parameters in ambient RNA
 */
void dll_ambmu(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i, vector<double>& results){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    int amb_idx = data_i.at("ef_idx");
    double p_c = params[amb_idx];
    double c = data_d.at("c");
    double binom_p = (1.0-c) * p_e + c*p_c;

    // Clamp binom_p away from 0 and 1 to prevent division by zero in gradient
    if (binom_p < DBL_MIN*1e6){
        binom_p = DBL_MIN*1e6;
    }
    else if (binom_p > 1.0 - DBL_MIN*1e6){
        binom_p = 1.0 - DBL_MIN*1e6;
    }

    double dy_dp = (k-n*binom_p)/(binom_p - binom_p*binom_p);
    if (!std::isfinite(dy_dp)){
        return;
    }
    results[amb_idx] = dy_dp * c;
}

/**
 * Log likelihood function for finding best e_r & e_a
 */
double ll_err_rates(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i){
    double e_r = params[0];
    double e_a = params[1];
    double p_c = data_d.at("p_c");
    double p_e = data_d.at("p_e");
    double n = data_d.at("n");
    double k = data_d.at("k");
    double c = data_d.at("c");

    double p_e_a = adjust_p_err(p_e, e_r, e_a);
    double binom_p = (1-c)*p_e_a + c*p_c;
    if (binom_p < DBL_MIN*1e6){
        binom_p = DBL_MIN*1e6;
    }
    else if (binom_p > 1.0-DBL_MIN*1e6){
        binom_p = 1.0-DBL_MIN*1e6;
    }
    
    return logbinom(n, k, binom_p);
}

/**
 * First derivative log likelihood function for finding best e_r & e_a
 */
void dll_err_rates(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i, vector<double>& results){
    
    double e_r = params[0];
    double e_a = params[1];
    double p_c = data_d.at("p_c");
    double p_e = data_d.at("p_e");
    double n = data_d.at("n");
    double k = data_d.at("k");
    double c = data_d.at("c");

    double p_e_a = adjust_p_err(p_e, e_r, e_a);
    double binom_p = (1-c)*p_e_a + c*p_c;
    if (binom_p < DBL_MIN*1e6){
        binom_p = DBL_MIN*1e6;
    }
    else if (binom_p > 1.0-DBL_MIN*1e6){
        binom_p = 1.0-DBL_MIN*1e6;
    }
    double dy_dp = (k-n*binom_p)/(binom_p - binom_p*binom_p);
    
    results[0] += dy_dp * (1 - p_e + c*(p_e - p_c));
    results[1] += dy_dp * (c*(p_e - p_c) - p_e);
    
}

/**
 * optimML-format log likelihood function for updating the 
 * ambient RNA profile mixture proportions.
 */
double ll_amb_prof_mixture(const vector<double>& params, 
    const map<string, double>& data_d, 
    const map<string, int>& data_i){
    
    double c;
    if (data_d.count("c") > 0){
        c = data_d.at("c");
    }    
    else{
        c = params[0];
    }
    double p_c = params[params.size()-1];
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    double binom_p = (1.0 - c) * p_e + c * p_c;
    
    // Clamp binom_p to prevent NaN/Inf from logbinom
    if (binom_p < DBL_MIN*1e6){
        binom_p = DBL_MIN*1e6;
    }
    else if (binom_p > 1.0 - DBL_MIN*1e6){
        binom_p = 1.0 - DBL_MIN*1e6;
    }
    if (isnan(binom_p)){
        binom_p = 0.5;
    }
    return logbinom(n, k, binom_p);
}

/**
 * Derivative of log likelihood function for updating
 * ambient RNA profile mixture proportions
 * (optimML-format)
 */
void dll_amb_prof_mixture(const vector<double>& params, 
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    double c;
    bool c_in_data = false;
    if (data_d.count("c") > 0){
        c = data_d.at("c");
        c_in_data = true;
    }
    else{
        c = params[0];
    }
    double p_c = params[params.size()-1];
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    
    double binom_p = (1.0 - c) * p_e + c * p_c;

    // Clamp binom_p away from 0 and 1 to prevent division by zero in gradient
    if (binom_p < DBL_MIN*1e6){
        binom_p = DBL_MIN*1e6;
    }
    else if (binom_p > 1.0 - DBL_MIN*1e6){
        binom_p = 1.0 - DBL_MIN*1e6;
    }
    if (isnan(binom_p)){
        binom_p = 0.5;
    }

    double dy_dp = (k-n*binom_p)/(binom_p - binom_p*binom_p);

    // Guard against non-finite gradient (e.g., from degenerate data)
    if (!std::isfinite(dy_dp)){
        return;
    }

    results[results.size()-1] += dy_dp * c;
    if (c_in_data){
        results[0] += dy_dp * (p_c - p_e);
    }
}

// ===== contamFinder functions =====

/**
 * Class constructor
 *
 * Stores copies of data structures/values that will be used by
 * many functions.
 *
 * Compiles data in the form that will be needed by later functions.
 */
contamFinder3::contamFinder3(robin_hood::unordered_map<unsigned long, 
        map<pair<int, int>, map<pair<int, int>, pair<float, float> > > >& indv_allelecounts,
    robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    map<pair<int, int>, map<int, float> >& exp_match_fracs,
    int n_samples,
    set<int>& allowed_ids,
    set<int>& allowed_ids2){
    
    this->contam_prof_initialized = false;
    this->c_init = -1;

    fprintf(stderr, "ambient_rna v%s: %s\n", AMBIENT_RNA_VERSION.c_str(), 
        AMBIENT_RNA_VERSION_MSG.c_str());

    // Copy external data structures that we will need in the future 
    this->assn = assn;
    this->assn_llr = assn_llr;
    this->indv_allelecounts = indv_allelecounts;
    this->n_samples = n_samples;
    this->n_mixprop_trials = 10;
    this->expfracs = exp_match_fracs;
    this->allowed_ids = allowed_ids;
    this->allowed_ids2 = allowed_ids2;
    
    // Reassign cells by default
    this->skip_reassign = false;

    // Don't worry about doublet rate by default
    this->doublet_rate = -1;
    
    this->ef_all_avg = true;

    llrtot = 0.0;

    bool caller_restricted_ids = (allowed_ids.size() > 0);

    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != 
        assn.end(); ++a){
        // If the caller provided an allowed set (from -i/-I or --expected_lines),
        // treat it as authoritative. Do not re-add out-of-pool assignments from
        // the original demux call, because a few bad assignments can otherwise
        // make absent species legal ambient sources.
        if (!caller_restricted_ids){
            this->allowed_ids.insert(a->second);
            // Make sure we also allow sub-IDs of combinations
            if (a->second >= n_samples){
                pair<int, int> combo = idx_to_hap_comb(a->second, n_samples);
                this->allowed_ids.insert(combo.first);
                this->allowed_ids.insert(combo.second);
            }
        }
        if (id_llrsum.count(a->second) == 0){
            id_llrsum.insert(make_pair(a->second, 0.0));
            id_count.insert(make_pair(a->second, 0.0));
        }
        id_llrsum[a->second] += assn_llr[a->first];
        id_count[a->second]++;
        llrtot += assn_llr[a->first];
        if (a->second >= n_samples){
            pair<int, int> comb = idx_to_hap_comb(a->second, n_samples);
            if (id_llrsum2.count(comb.first) == 0){
                id_llrsum2.insert(make_pair(comb.first, 0.0));
            }
            if (id_llrsum2.count(comb.second) == 0){
                id_llrsum2.insert(make_pair(comb.second, 0.0));
            }
            id_llrsum2[comb.first] += 0.5*assn_llr[a->first];
            id_llrsum2[comb.second] += 0.5*assn_llr[a->first];
        }
        else{
            if (id_llrsum2.count(a->second) == 0){
                id_llrsum2.insert(make_pair(a->second, 0.0));
            }
            id_llrsum2[a->second] += assn_llr[a->first];
        }
    }

    // For the ambient RNA mixture, we only want to model single (not doublet) 
    // individual genotypes that were allowed in this run (if no filter was given,
    // allow all single individuals from the VCF).

    if (this->allowed_ids.size() == 0){
        // Default to all possible individuals
        for (int i = 0; i < n_samples; ++i){
            idx2samp.push_back(i);
        }
    }
    else{
        for (set<int>::iterator a = this->allowed_ids.begin(); a != this->allowed_ids.end(); ++a){
            if (*a < n_samples){
                idx2samp.push_back(*a);
            }
        }
    }

    // Set some default parameter values.
    this->e_r = 1e-3;
    this->e_a = 1e-3;
    
    this->inter_species = false;

    this->contam_cell_prior = -1;
    this->contam_cell_prior_var = -1;

    this->delta_thresh = 0.1;
    this->maxits = 100;
    
    this->weighted = false;
    
    this->num_threads = 1;
    this->thorough_multistart = false;
    this->adaptive_multistart = true;

    // Tetraploid-aware defaults (v1.4)
    this->tetraploid_aware = false;
    this->min_signal_gap = 0.10;
    this->amb_mu_available = false;
    this->ids_restricted = false;

    // LOO ambient profile default
    this->use_loo = false;

    // FI-weight default
    this->use_fi_weight = false;

    // User-specified prior default (disabled)
    this->user_prior_set = false;
    this->user_prior_mean = -1.0;
    this->user_prior_var = -1.0;

    // R-feedback default (disabled)
    this->r_feedback_enabled = false;

    // Adaptive prior defaults (disabled)
    this->adaptive_prior_enabled = false;
    this->adaptive_prior_mean = 0.05;
    this->adaptive_prior_init_var = 0.04;
    this->adaptive_prior_min_var = 0.001;
    this->adaptive_prior_boundary_thresh = 0.20;
    this->adaptive_prior_boundary_low = 0.02;
    this->adaptive_prior_boundary_high = 0.90;
    this->adaptive_prior_max_shrink_steps = 8;

    // Fixed ambient profile default (disabled)
    this->fixed_amb_prof = false;

    // Species mode default (disabled)
    this->species_mode = false;
    this->species_init_used = false;
    this->has_species_counts = false;
    this->primary_species_counts_enabled = false;

    // Bulk mode default (disabled)
    this->bulk_mode = false;

    // Compile data in the format needed by other functions
    this->compile_data(assn, indv_allelecounts);
}

void contamFinder3::set_init_contam_prof(map<int, double>& cp){
    // Filter to only include individuals present in idx2samp (active individuals).
    // Warm-start files may contain entries for all VCF individuals (e.g. 28), but
    // idx2samp only has the subset present in cell assignments. Without filtering,
    // update_amb_prof_mixture builds startprops from contam_prof (28 entries) but
    // mixfracs columns from idx2samp (fewer entries), causing an optimML dimension
    // mismatch error. The -1 key (other_species) is preserved if present.
    // Missing active individuals are filled with epsilon to prevent contam_prof
    // from being shorter than idx2samp.
    std::set<int> active(idx2samp.begin(), idx2samp.end());
    this->contam_prof.clear();
    double kept_sum = 0.0;
    int n_dropped = 0;
    for (auto& kv : cp){
        if (kv.first == -1){
            // Always keep the other_species entry
            this->contam_prof[-1] = kv.second;
            kept_sum += kv.second;
        } else if (active.count(kv.first) > 0){
            this->contam_prof[kv.first] = kv.second;
            kept_sum += kv.second;
        } else {
            n_dropped++;
        }
    }
    // Fill any active idx2samp individuals missing from the input profile
    int n_filled = 0;
    double eps = 1e-4;
    for (int i = 0; i < (int)idx2samp.size(); i++){
        int samp = idx2samp[i];
        if (this->contam_prof.count(samp) == 0){
            this->contam_prof[samp] = eps;
            kept_sum += eps;
            n_filled++;
        }
    }
    // Re-normalize so proportions sum to 1.0
    if (kept_sum > 0.0){
        for (auto& kv : this->contam_prof){
            kv.second /= kept_sum;
        }
    }
    if (n_dropped > 0 || n_filled > 0){
        fprintf(stderr, "  set_init_contam_prof: filtered %d inactive, filled %d missing "
            "(%lu final, %lu active individuals)\n",
            n_dropped, n_filled, this->contam_prof.size(), idx2samp.size());
    }
    contam_prof_initialized = true;
}

/**
 * Functions to set parameter values
 */

void contamFinder3::set_error_rates(double error_ref, double error_alt){
    this->e_r = error_ref;
    this->e_a = error_alt;
}

void contamFinder3::set_doublet_rate(double d){
    this->doublet_rate = d;
}

void contamFinder3::model_other_species(){
    if (this->species_mode){
        fprintf(stderr, "ERROR: --other_species cannot be used with --species_mode. "
            "Species mode already captures all cross-species contributions explicitly.\n");
        exit(1);
    }
    this->inter_species = true;
}

void contamFinder3::model_single_species(){
    this->inter_species = false;
}

void contamFinder3::set_mixprop_trials(int nt){
    this->n_mixprop_trials = nt;
}

void contamFinder3::set_delta(double delta){
    this->delta_thresh = delta;
}

void contamFinder3::set_maxiter(int its){
    this->maxits = its;
}

void contamFinder3::use_weights(){
    this->weighted = true;
}

void contamFinder3::no_weights(){
    this->weighted = false;
}

void contamFinder3::no_reassign(){
    this->skip_reassign = true;
}

void contamFinder3::set_init_c(double c){
    c_init = c;
}

void contamFinder3::set_num_threads(int nt){
    num_threads = nt;
}

void contamFinder3::set_thorough_multistart(bool enabled){
    thorough_multistart = enabled;
}

void contamFinder3::set_adaptive_multistart(bool enabled){
    adaptive_multistart = enabled;
}

void contamFinder3::set_locked_identities(const set<int>& locked){
    this->locked_identities = locked;
}

void contamFinder3::set_safe_singlets(const set<int>& safe){
    this->safe_singlets = safe;
}

void contamFinder3::set_tetraploid_aware(bool enabled){
    this->tetraploid_aware = enabled;
}

void contamFinder3::set_min_signal_gap(double gap){
    this->min_signal_gap = gap;
}

void contamFinder3::set_ids_restricted(bool restricted){
    this->ids_restricted = restricted;
}

void contamFinder3::set_use_loo(bool enabled){
    this->use_loo = enabled;
}

void contamFinder3::set_fi_weight(bool enabled){
    this->use_fi_weight = enabled;
}

void contamFinder3::set_user_prior(double mean, double var){
    this->user_prior_set = true;
    this->user_prior_mean = mean;
    this->user_prior_var = var;
    // Apply immediately so the first est_contam_cells pass uses it
    this->contam_cell_prior = mean;
    this->contam_cell_prior_var = var;
}

void contamFinder3::set_r_feedback(bool enabled){
    this->r_feedback_enabled = enabled;
}

void contamFinder3::set_adaptive_prior(bool enabled,
                                        double mean,
                                        double init_var,
                                        double min_var,
                                        double boundary_thresh){
    this->adaptive_prior_enabled = enabled;
    this->adaptive_prior_mean = mean;
    this->adaptive_prior_init_var = init_var;
    this->adaptive_prior_min_var = min_var;
    this->adaptive_prior_boundary_thresh = boundary_thresh;
}

// =============================================================================
// Fixed ambient profile (Step 0a)
// =============================================================================

void contamFinder3::set_fixed_amb_prof(bool enabled){
    this->fixed_amb_prof = enabled;
}

void contamFinder3::rebuild_amb_mu_from_contam_prof(){
    // Extracted from the "Update stored ambient RNA profile" block at the end of
    // update_amb_prof_mixture. Rebuilds amb_mu from contam_prof and expfracs.
    for (int x = 0; x < (int)idx2samp.size(); ++x){
        int i = idx2samp[x];
        for (int nalt = 0; nalt <= 2; ++nalt){
            pair<int, int> key = make_pair(i, nalt);
            if (amb_mu.count(key) == 0){
                map<pair<int, int>, double> m;
                amb_mu.insert(make_pair(key, m));
            }
            pair<int, int> nullkey = make_pair(-1,-1);
            if (amb_mu[key].count(nullkey) == 0){
                amb_mu[key].insert(make_pair(nullkey, 0));
            }
            double val = 0.0;
            for (map<int, double>::iterator cp = contam_prof.begin();
                cp != contam_prof.end(); ++cp){
                if (cp->first == -1){
                    val += cp->second * 0.0;
                }
                else if (ef_all_avg && cp->first == key.first){
                    val += cp->second * ((double)key.second/2.0);
                }
                else{
                    double ef_val = (double)key.second / 2.0;
                    if (expfracs.count(key) > 0 &&
                        expfracs[key].count(cp->first) > 0){
                        ef_val = expfracs[key][cp->first];
                    }
                    val += cp->second * ef_val;
                }
            }
            amb_mu[key][nullkey] = val;
            for (int y = x + 1; y < (int)idx2samp.size(); ++y){
                int j = idx2samp[y];
                for (int nalt2 = 0; nalt2 <= 2; ++nalt2){
                    pair<int, int> key2 = make_pair(j, nalt2);
                    if (amb_mu[key].count(key2) == 0){
                        amb_mu[key].insert(make_pair(key2, 0));
                    }
                    double val = 0.0;
                    for (map<int, double>::iterator cp = contam_prof.begin();
                        cp != contam_prof.end(); ++cp){
                        if (cp->first == -1){
                            val += cp->second * adjust_p_err(0.0, e_r, e_a);
                        }
                        else{
                            if (ef_all_avg && cp->first == key.first){
                                val += cp->second * adjust_p_err(
                                    (double)key.second / 2.0, e_r, e_a);
                            }
                            else if (ef_all_avg && cp->first == key2.first){
                                val += cp->second * adjust_p_err(
                                    (double)key2.second / 2.0, e_r, e_a);
                            }
                            else{
                                double ef1 = (double)key.second / 2.0;
                                double ef2 = (double)key2.second / 2.0;
                                if (expfracs.count(key) > 0 &&
                                    expfracs[key].count(cp->first) > 0){
                                    ef1 = expfracs[key][cp->first];
                                }
                                if (expfracs.count(key2) > 0 &&
                                    expfracs[key2].count(cp->first) > 0){
                                    ef2 = expfracs[key2][cp->first];
                                }
                                val += cp->second * (0.5*ef1 + 0.5*ef2);
                            }
                        }
                    }
                    amb_mu[key][key2] = val;
                }
            }
        }
    }
}

// =============================================================================
// Species mode
// =============================================================================

void contamFinder3::set_species_mode(const PanelMetadata& panel){
    this->species_mode = true;
    this->panel_meta = panel;
    // Two-way protection: species mode forces inter_species off
    if (this->inter_species){
        fprintf(stderr, "WARNING: --species_mode disables --other_species "
            "(species-level pi already captures cross-species contributions)\n");
    }
    this->inter_species = false;
}

void contamFinder3::set_species_prior(const map<string, double>& species_prof,
                                       const map<string, double>& species_prof_conc){
    this->species_prior_prof = species_prof;
    this->species_prior_conc = species_prof_conc;
}

void contamFinder3::set_species_init(const map<string, double>& init_prof){
    this->species_init_prof = init_prof;
    this->species_init_used = false;
}

void contamFinder3::set_indiv_init(const map<int, double>& init_prof){
    // Warm-start the individual-level ambient profile from an external file.
    // Filter to idx2samp, fill missing, normalize.
    set<int> active(idx2samp.begin(), idx2samp.end());
    contam_prof.clear();
    double sum = 0.0;
    for (int i = 0; i < (int)idx2samp.size(); i++){
        int samp = idx2samp[i];
        double val = 1e-6;
        if (init_prof.count(samp) > 0){
            val = init_prof.at(samp);
            if (val < 1e-6) val = 1e-6;
        }
        contam_prof[samp] = val;
        sum += val;
    }
    if (inter_species){
        double os_val = 1e-6;
        if (init_prof.count(-1) > 0){
            os_val = init_prof.at(-1);
            if (os_val < 1e-6) os_val = 1e-6;
        }
        contam_prof[-1] = os_val;
        sum += os_val;
    }
    // Normalize
    for (auto& cp : contam_prof){
        cp.second /= sum;
    }
    contam_prof_initialized = true;
}

void contamFinder3::compute_loading_weights(){
    indiv_loading_weights.clear();

    // Count cells assigned to each singlet individual.
    // Doublets contribute 0.5 to each constituent.
    for (auto& a : assn){
        int idx = a.second;
        if (idx < n_samples){
            indiv_loading_weights[idx] += 1.0;
        } else {
            std::pair<short, short> combo = idx_to_hap_comb(idx, n_samples);
            indiv_loading_weights[combo.first] += 0.5;
            indiv_loading_weights[combo.second] += 0.5;
        }
    }

    // If no assignments exist but contam_prof is initialized (warm start),
    // use the contam_prof values as proxy weights.
    if (indiv_loading_weights.empty() && contam_prof_initialized){
        for (auto& cp : contam_prof){
            if (cp.first >= 0){
                indiv_loading_weights[cp.first] = cp.second;
            }
        }
    }

    // Normalize within each species so weights sum to 1.0 per species.
    for (const auto& sp : panel_meta.species_list){
        double sp_sum = 0.0;
        if (panel_meta.species_to_sample_indices.count(sp) > 0){
            for (int idx : panel_meta.species_to_sample_indices.at(sp)){
                if (indiv_loading_weights.count(idx) > 0){
                    sp_sum += indiv_loading_weights[idx];
                }
            }
        }
        if (sp_sum > 0.0){
            if (panel_meta.species_to_sample_indices.count(sp) > 0){
                for (int idx : panel_meta.species_to_sample_indices.at(sp)){
                    if (indiv_loading_weights.count(idx) > 0){
                        indiv_loading_weights[idx] /= sp_sum;
                    }
                }
            }
        }
    }
}

void contamFinder3::set_primary_species_counts_enabled(bool enabled){
    this->primary_species_counts_enabled = enabled;
}

void contamFinder3::set_species_counts(
    robin_hood::unordered_map<unsigned long,
        map<pair<int, int>,
            map<pair<int, int>,
                pair<float, float> > > >& counts,
    map<pair<int, int>, map<int, float> >& condf){
    this->species_allelecounts = counts;
    this->species_expfracs = condf;
    this->has_species_counts = true;
    fprintf(stderr, "  Species-diagnostic counts: %lu cells\n", counts.size());
}

void contamFinder3::expand_species_prior_to_indiv(){
    // Expand species-level pi to individual-level via weighted split.
    // Each individual's contribution to a species is weighted by
    // panel_meta.get_weight(species, idx) -- normally 1.0, but 0.5 for
    // hybrid individuals folded into two parent species.
    // An individual appearing in multiple species (e.g. Hy folded into C+B)
    // accumulates contributions from each species.
    // Only populate contam_prof for individuals in idx2samp. When a species
    // has zero active weight, redistribute its mass proportionally.
    contam_prof.clear();
    species_contam_prof.clear();
    species_contam_prof_conc.clear();

    set<int> active_samples(idx2samp.begin(), idx2samp.end());

    // First pass: compute weighted active count per species, accumulate orphan mass
    double orphan_mass = 0.0;
    double active_mass = 0.0;
    map<string, double> sp_w_active;  // sum of weights for active individuals
    for (const auto& sp : panel_meta.species_list){
        double pi_s = 0.0;
        if (species_prior_prof.count(sp) > 0){
            pi_s = species_prior_prof[sp];
        }
        species_contam_prof[sp] = pi_s;
        if (species_prior_conc.count(sp) > 0){
            species_contam_prof_conc[sp] = species_prior_conc[sp];
        }

        double w_sum = 0.0;
        const auto& indices = panel_meta.species_to_sample_indices;
        if (indices.count(sp) > 0){
            for (int idx : indices.at(sp)){
                if (active_samples.count(idx) > 0){
                    w_sum += panel_meta.get_weight(sp, idx);
                }
            }
        }
        sp_w_active[sp] = w_sum;
        if (w_sum <= 0.0){
            orphan_mass += pi_s;
        } else {
            active_mass += pi_s;
        }
    }

    // Second pass: assign individual-level proportions with redistributed mass.
    //
    // Loading weight semantics: when indiv_loading_weights is non-empty, a
    // missing entry means zero assigned cells, so effective weight is 0.0 (not
    // implicit 1.0). If all individuals in a species have zero loading weight,
    // fall back to base panel weights for that species to avoid dropping its
    // entire mass.
    double scale = (active_mass > 0) ? (active_mass + orphan_mass) / active_mass : 1.0;
    for (const auto& sp : panel_meta.species_list){
        if (sp_w_active[sp] <= 0.0) continue;
        double pi_s = species_contam_prof[sp] * scale;
        const auto& indices = panel_meta.species_to_sample_indices;

        // Compute effective weight sum for this species.
        // Try loading-weighted first; fall back to base if the species has
        // zero total effective loading weight.
        double eff_w_sum = 0.0;
        bool use_loading = !indiv_loading_weights.empty();
        if (use_loading){
            for (int idx : indices.at(sp)){
                if (active_samples.count(idx) > 0){
                    auto it = indiv_loading_weights.find(idx);
                    if (it != indiv_loading_weights.end()){
                        eff_w_sum += panel_meta.get_weight(sp, idx) * it->second;
                    }
                    // missing entry => 0.0, contributes nothing
                }
            }
            if (eff_w_sum <= 0.0){
                // Entire species has zero loading weight (no cells assigned to
                // any individual in this species). Fall back to base panel
                // weights so the species mass is not silently dropped.
                use_loading = false;
            }
        }
        if (!use_loading){
            eff_w_sum = 0.0;
            for (int idx : indices.at(sp)){
                if (active_samples.count(idx) > 0){
                    eff_w_sum += panel_meta.get_weight(sp, idx);
                }
            }
        }
        if (eff_w_sum <= 0.0) continue;

        for (int idx : indices.at(sp)){
            if (active_samples.count(idx) > 0){
                double w;
                if (use_loading){
                    auto it = indiv_loading_weights.find(idx);
                    w = (it != indiv_loading_weights.end())
                      ? panel_meta.get_weight(sp, idx) * it->second
                      : 0.0;
                } else {
                    w = panel_meta.get_weight(sp, idx);
                }
                if (w > 0.0){
                    contam_prof[idx] += (w / eff_w_sum) * pi_s;
                }
            }
        }
    }
}

double contamFinder3::solve_species_level_pi(){
    // Build species-aggregated mixture components for the solver.
    // Each species gets one component whose mixfracs vector is the
    // average of its constituent individuals' expfracs entries.

    // Ensure loading weights are populated before any species-to-individual
    // expansion. On the first call (from init_params), loading weights may not
    // yet exist. Without them, expand_species_prior_to_indiv uses only base
    // panel weights, which is less accurate when individuals within a species
    // have unequal cell counts.
    // Skip in bulk mode: the synthetic placeholder assignment would produce
    // meaningless loading weights (one cell assigned to index 0).
    if (!bulk_mode && indiv_loading_weights.empty()){
        compute_loading_weights();
    }

    vector<vector<double> > mixfracs;
    vector<double> weights;
    vector<double> n;
    vector<double> k;
    vector<double> p_e;
    vector<double> c;

    // Choose how to get per-observation c values:
    // If per-cell contam_rate has been estimated, use per-cell c (use_global_c=false).
    // If contam_rate is empty (first call from init_params), use the global prior
    // or a default. Using per-cell=0 would zero out all mixture signal.
    bool use_global = contam_rate.empty();
    if (use_global && contam_cell_prior <= 0){
        // No global prior set yet either. Use a conservative default.
        contam_cell_prior = 0.05;
    }

    // When species-diagnostic counts are available, swap them in for
    // compile_amb_prof_dat. The species-diagnostic SNPs have near-fixed
    // allele frequencies within each species, giving the mixture solver
    // much cleaner signal than the individual-diagnostic SNPs.
    if (has_species_counts){
        std::swap(indv_allelecounts, species_allelecounts);
        std::swap(expfracs, species_expfracs);
    }

    if (bulk_mode){
        compile_bulk_amb_prof_dat(mixfracs, weights, n, k, p_e, c);
    } else {
        compile_amb_prof_dat(false, use_global, mixfracs, weights, n, k, p_e, c);
    }

    // Swap back so individual-level estimation uses the original data
    if (has_species_counts){
        std::swap(indv_allelecounts, species_allelecounts);
        std::swap(expfracs, species_expfracs);
    }

    if (mixfracs.empty()){
        fprintf(stderr, "WARNING: no data for species-level pi estimation\n");
        return 0.0;
    }

    // Remap per-individual mixture columns to per-species using weights.
    // An individual may contribute to multiple species (e.g. Hy folded into C+B),
    // so we use species_to_sample_indices + get_weight rather than a 1:1 map.
    int n_obs = (int)mixfracs.size();
    int n_species = (int)panel_meta.species_list.size();

    // Build idx2samp position lookup: sample_index -> position in idx2samp
    map<int, int> samp_to_pos;
    for (int i = 0; i < (int)idx2samp.size(); i++){
        samp_to_pos[idx2samp[i]] = i;
    }

    // Build samp_to_species_idx for the samp_to_species_idx lookups elsewhere
    // (used later for startprops). A sample can map to multiple species now,
    // but this map only stores the last one seen; that's OK because we don't
    // use it for the aggregation anymore.
    map<int, int> samp_to_species_idx;
    for (int s = 0; s < n_species; s++){
        const string& sp = panel_meta.species_list[s];
        if (panel_meta.species_to_sample_indices.count(sp) > 0){
            for (int idx : panel_meta.species_to_sample_indices.at(sp)){
                samp_to_species_idx[idx] = s;
            }
        }
    }

    // Compute weighted active count per species
    vector<double> species_weight_sum(n_species, 0.0);
    set<int> active_samples_set(idx2samp.begin(), idx2samp.end());
    for (int s = 0; s < n_species; s++){
        const string& sp = panel_meta.species_list[s];
        if (panel_meta.species_to_sample_indices.count(sp) > 0){
            for (int idx : panel_meta.species_to_sample_indices.at(sp)){
                if (active_samples_set.count(idx) > 0){
                    species_weight_sum[s] += panel_meta.get_weight(sp, idx);
                }
            }
        }
    }

    // Identify active species (those with at least one individual in idx2samp).
    // Species with zero active weight produce all-zero mixture columns, making
    // them invisible to the likelihood. The solver can assign arbitrary mass to
    // such columns without penalty, causing spurious nonzero estimates for
    // species not present in the pool. Exclude them from the solver entirely.
    vector<int> active_sp_idx;      // full-species-index for each active slot
    map<int, int> full_to_active;   // full species index -> active slot
    for (int s = 0; s < n_species; s++){
        if (species_weight_sum[s] > 0.0){
            full_to_active[s] = (int)active_sp_idx.size();
            active_sp_idx.push_back(s);
        }
    }
    int n_active_species = (int)active_sp_idx.size();

    if (n_active_species == 0){
        fprintf(stderr, "WARNING: no active species for species-level pi estimation\n");
        return 0.0;
    }

    // Log which species are active vs excluded
    {
        vector<string> active_names, excluded_names;
        for (int s = 0; s < n_species; s++){
            if (full_to_active.count(s) > 0){
                active_names.push_back(panel_meta.species_list[s]);
            } else {
                excluded_names.push_back(panel_meta.species_list[s]);
            }
        }
        fprintf(stderr, "  Species solver: %d active", n_active_species);
        for (const auto& nm : active_names) fprintf(stderr, " %s", nm.c_str());
        if (!excluded_names.empty()){
            fprintf(stderr, " | excluded (no individuals in pool):");
            for (const auto& nm : excluded_names) fprintf(stderr, " %s", nm.c_str());
        }
        fprintf(stderr, "\n");
    }

    // Aggregate: for each observation, build species-level mixfrac as
    // weighted average of individual mixfracs within each species.
    // Only active species get columns in the solver matrix.
    vector<vector<double> > species_mixfracs(n_obs, vector<double>(n_active_species, 0.0));

    for (int obs = 0; obs < n_obs; obs++){
        for (int a = 0; a < n_active_species; a++){
            int s = active_sp_idx[a];
            const string& sp = panel_meta.species_list[s];
            double wsum = 0.0;
            if (panel_meta.species_to_sample_indices.count(sp) > 0){
                for (int idx : panel_meta.species_to_sample_indices.at(sp)){
                    if (samp_to_pos.count(idx) > 0){
                        double w = panel_meta.get_weight(sp, idx);
                        wsum += w * mixfracs[obs][samp_to_pos[idx]];
                    }
                }
            }
            species_mixfracs[obs][a] = wsum / species_weight_sum[s];
        }
    }

    // Build starting proportions (active species only).
    // On the first call, use species_init_prof if provided (warm start).
    // On subsequent calls, use the previous solution stored in species_contam_prof.
    vector<double> startprops(n_active_species, 0.0);
    bool have_init = false;

    if (!species_init_used && !species_init_prof.empty()){
        // First call with warm-start from --species_init
        // Collect only active species values; mass from inactive species
        // is redistributed proportionally among active ones.
        double sp_sum = 0.0;
        double orphan = 0.0;
        for (int s = 0; s < n_species; s++){
            const string& sp = panel_meta.species_list[s];
            double val = 0.0;
            if (species_init_prof.count(sp) > 0){
                val = species_init_prof[sp];
            }
            if (full_to_active.count(s) > 0){
                startprops[full_to_active[s]] = val;
                sp_sum += val;
            } else {
                orphan += val;
            }
        }
        // Redistribute orphan mass proportionally
        if (orphan > 0.0 && sp_sum > 0.0){
            double scale = (sp_sum + orphan) / sp_sum;
            for (int a = 0; a < n_active_species; a++){
                startprops[a] *= scale;
            }
            sp_sum += orphan;
        }
        // Normalize
        if (sp_sum > 0.0){
            for (int a = 0; a < n_active_species; a++){
                startprops[a] /= sp_sum;
            }
            have_init = true;
        }
        species_init_used = true;
        fprintf(stderr, "  Species warm start from --species_init:");
        for (int a = 0; a < n_active_species; a++){
            fprintf(stderr, " %s=%.4f",
                panel_meta.species_list[active_sp_idx[a]].c_str(), startprops[a]);
        }
        fprintf(stderr, "\n");
    } else if (!species_contam_prof.empty()){
        // Subsequent calls: use previous solution (active species only)
        double sp_sum = 0.0;
        for (int a = 0; a < n_active_species; a++){
            const string& sp = panel_meta.species_list[active_sp_idx[a]];
            if (species_contam_prof.count(sp) > 0){
                startprops[a] = species_contam_prof[sp];
            }
            sp_sum += startprops[a];
        }
        if (sp_sum > 0.0){
            for (int a = 0; a < n_active_species; a++) startprops[a] /= sp_sum;
            have_init = true;
        }
    }

    if (!have_init){
        // Fall back to uniform across active species only
        for (int a = 0; a < n_active_species; a++){
            startprops[a] = 1.0 / (double)n_active_species;
        }
    }

    // Solve with active species only
    vector<double> params;
    optimML::multivar_ml_solver solver(params, ll_amb_prof_mixture, dll_amb_prof_mixture);
    if (num_threads > 1){
        solver.set_bfgs_threads(num_threads);
    }
    solver.add_mixcomp(species_mixfracs);
    solver.add_mixcomp_fracs(startprops);
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("p_e", p_e);
    solver.add_data("c", c);
    solver.add_weights(weights);

    bool solve_ok = false;
    double maxll = -1e30;
    vector<double> best_mixcomp;

    try {
        solver.solve();
        if (std::isfinite(solver.log_likelihood)){
            solve_ok = true;
            maxll = solver.log_likelihood;
            best_mixcomp = solver.results_mixcomp;
        }
    } catch (...) {
        fprintf(stderr, "WARNING: species-level mixture solve failed on initial attempt\n");
    }

    // Try randomized starting conditions
    if (!solve_ok || n_mixprop_trials > 0){
        int ntrials = n_active_species * n_mixprop_trials;
        for (int trial = 0; trial < ntrials; trial++){
            try {
                solver.randomize_mixcomps();
                solver.solve();
                if (std::isfinite(solver.log_likelihood) && solver.log_likelihood > maxll){
                    maxll = solver.log_likelihood;
                    best_mixcomp = solver.results_mixcomp;
                    solve_ok = true;
                }
            } catch (...) {
                // ignore, try next
            }
        }
    }

    if (!solve_ok){
        fprintf(stderr, "WARNING: all species-level mixture solves failed; "
            "falling back to uniform proportions across active species\n");
        // CRITICAL: contam_prof may contain stale entries from est_min_c()
        // with wrong dimensions (all VCF individuals vs idx2samp subset).
        // Initialize to uniform across ACTIVE species only to prevent
        // downstream dimension mismatch and spurious inactive species mass.
        species_contam_prof.clear();
        species_prior_prof.clear();
        for (int s = 0; s < n_species; s++){
            if (full_to_active.count(s) > 0){
                double unif = 1.0 / (double)n_active_species;
                species_contam_prof[panel_meta.species_list[s]] = unif;
                species_prior_prof[panel_meta.species_list[s]] = unif;
            } else {
                species_contam_prof[panel_meta.species_list[s]] = 0.0;
                species_prior_prof[panel_meta.species_list[s]] = 0.0;
            }
        }
        expand_species_prior_to_indiv();
        rebuild_amb_mu_from_contam_prof();
        return 0.0;
    }

    // Map active-species solver results back to full species list.
    // Inactive species get 0.0.
    species_contam_prof.clear();
    for (int s = 0; s < n_species; s++){
        if (full_to_active.count(s) > 0){
            species_contam_prof[panel_meta.species_list[s]] = best_mixcomp[full_to_active[s]];
        } else {
            species_contam_prof[panel_meta.species_list[s]] = 0.0;
        }
    }

    // Expand to individual-level via uniform-within-species split
    species_prior_prof = species_contam_prof;
    expand_species_prior_to_indiv();

    // Rebuild amb_mu from individual-level contam_prof
    rebuild_amb_mu_from_contam_prof();

    return maxll;
}

// =============================================================================
// Bulk mode (for quant3_contam_empty_drops)
// =============================================================================

void contamFinder3::set_bulk_mode(bool enabled){
    this->bulk_mode = enabled;
    if (enabled){
        this->skip_reassign = true;
    }
}

void contamFinder3::report_r_feedback_stats(){
    if (allele_ratio.empty()){
        fprintf(stderr, "R-feedback: no allele_ratio estimates available yet\n");
        return;
    }
    double r_sum = 0.0, r_min = 1.0, r_max = 0.0;
    int r_count = 0;
    for (auto& ar : allele_ratio){
        r_sum += ar.second;
        if (ar.second < r_min) r_min = ar.second;
        if (ar.second > r_max) r_max = ar.second;
        r_count++;
    }
    fprintf(stderr, "R-feedback: %d heterotypic cells, mean_r=%.4f, range=[%.4f, %.4f]\n",
        r_count, r_sum / r_count, r_min, r_max);
}

contamFinder3::BoundaryDiag contamFinder3::diagnose_contam_distribution(){
    BoundaryDiag diag = {0, 0, 0, 0, 0};
    vector<double> all_c;

    for (auto& cr : contam_rate){
        all_c.push_back(cr.second);
    }
    if (all_c.empty()) return diag;

    diag.n_cells = all_c.size();
    int n_low = 0, n_high = 0;
    for (double c : all_c){
        if (c < adaptive_prior_boundary_low) n_low++;
        if (c > adaptive_prior_boundary_high) n_high++;
    }
    diag.frac_low = (double)n_low / diag.n_cells;
    diag.frac_high = (double)n_high / diag.n_cells;
    diag.frac_boundary = diag.frac_low + diag.frac_high;

    sort(all_c.begin(), all_c.end());
    diag.median_c = all_c[all_c.size() / 2];
    return diag;
}

void contamFinder3::run_adaptive_prior(){
    BoundaryDiag bd = diagnose_contam_distribution();
    fprintf(stderr, "Adaptive prior check: boundary_frac=%.3f (low=%.3f, high=%.3f), "
            "median_c=%.4f, n_cells=%d\n",
            bd.frac_boundary, bd.frac_low, bd.frac_high, bd.median_c, bd.n_cells);

    if (bd.frac_boundary <= adaptive_prior_boundary_thresh){
        fprintf(stderr, "Distribution looks acceptable. Keeping current estimates.\n");
        return;
    }

    fprintf(stderr, "Distribution is pathological (boundary_frac %.3f > %.3f). "
            "Switching to adaptive fixed prior.\n",
            bd.frac_boundary, adaptive_prior_boundary_thresh);

    double current_var = adaptive_prior_init_var;

    for (int step = 0; step < adaptive_prior_max_shrink_steps; step++){
        if (current_var < adaptive_prior_min_var) break;

        contam_cell_prior = adaptive_prior_mean;
        contam_cell_prior_var = current_var;
        user_prior_set = true;

        est_contam_cells();

        BoundaryDiag bd2 = diagnose_contam_distribution();
        fprintf(stderr, "  Step %d: var=%.6f -> boundary_frac=%.3f "
                "(low=%.3f, high=%.3f), median_c=%.4f\n",
                step + 1, current_var, bd2.frac_boundary,
                bd2.frac_low, bd2.frac_high, bd2.median_c);

        if (bd2.frac_boundary <= adaptive_prior_boundary_thresh){
            fprintf(stderr, "  Distribution looks acceptable at var=%.6f. Stopping.\n",
                    current_var);
            return;
        }

        current_var *= 0.5;
    }

    fprintf(stderr, "  WARNING: Reached shrinkage limit (min_var=%.6f, max_steps=%d). "
            "Distribution may still be pathological.\n",
            adaptive_prior_min_var, adaptive_prior_max_shrink_steps);
}

// ============================================================================
// Helper: check if a SNP category is informative for contamination estimation
// ============================================================================

// Hard filter: for combo identities, skip categories where (nalt1 + nalt2) is 1 or 3
// These produce p_e = 0.25 or 0.75 which are close to typical pool averages
static bool category_passes_hard_filter(bool is_combo, int nalt1, int nalt2){
    if (!is_combo) return true;  // diploid categories always pass
    int sum = nalt1 + nalt2;
    return (sum != 1 && sum != 3);
}

// Adaptive filter: skip categories where |p_e - p_c| < min_signal_gap
static bool category_passes_adaptive_filter(double p_e, double p_c, double min_gap){
    return fabs(p_e - p_c) >= min_gap;
}

// Fisher Information weight: continuous weight based on how informative a category is
// for distinguishing contamination. The FI of a binomial observation for parameter c is:
//   FI(c) = n * (p_c - p_e)^2 / (p*(1-p))
// where p = (1-c)*p_e + c*p_c. We use (p_c - p_e)^2 as the signal strength since
// the denominator varies per-cell. Categories with p_c near p_e carry near-zero
// information and get downweighted continuously instead of being hard-filtered.
static double compute_fi_weight(double p_e, double p_c){
    double gap = p_c - p_e;
    return gap * gap;
}


bool contamFinder3::has_composition_override(unsigned long barcode) const{
    auto it = cell_composition_overrides.find(barcode);
    return it != cell_composition_overrides.end() && !it->second.empty();
}

double contamFinder3::composition_expected_from_row(
    const map<int, double>& comp,
    const pair<int, int>& key1,
    const pair<int, int>& key2) const{

    double out = 0.0;
    double wsum = 0.0;
    for (auto it = comp.begin(); it != comp.end(); ++it){
        int idx = it->first;
        double w = it->second;
        if (w <= 0.0) continue;

        double val = 0.5;
        bool have = false;

        if (key1.first == idx){
            val = (double)key1.second / 2.0;
            have = true;
        } else if (key2.first == idx){
            val = (double)key2.second / 2.0;
            have = true;
        } else {
            // For weighted identities with >2 native species components, the
            // counts file is still pairwise.  Estimate missing component dosage
            // from the conditional matching fractions for the observed row.  If
            // both pair axes can inform the missing species, average them.
            double acc = 0.0;
            int nhave = 0;
            auto e1 = expfracs.find(key1);
            if (e1 != expfracs.end()){
                auto v = e1->second.find(idx);
                if (v != e1->second.end()){
                    acc += v->second;
                    nhave++;
                }
            }
            if (key2.first != -1){
                auto e2 = expfracs.find(key2);
                if (e2 != expfracs.end()){
                    auto v = e2->second.find(idx);
                    if (v != e2->second.end()){
                        acc += v->second;
                        nhave++;
                    }
                }
            }
            if (nhave > 0){
                val = acc / (double)nhave;
                have = true;
            }
        }

        // If no conditional information is available, keep neutral 0.5 rather
        // than inventing a species-specific dosage.
        (void)have;
        out += w * val;
        wsum += w;
    }

    if (wsum > 0.0) out /= wsum;
    return out;
}

bool contamFinder3::composition_row_is_relevant(
    const map<int, double>& comp,
    const pair<int, int>& key1,
    const pair<int, int>& key2) const{

    if (comp.empty()) return false;
    if (comp.size() == 1){
        int only = comp.begin()->first;
        return key1.first == only && key2.first == -1;
    }

    // Pairwise counts are the native artifact currently available.  For a
    // weighted composition with k>1 species, use pair rows where both axes are
    // non-null and both species are among the nonzero composition components.
    if (key2.first < 0) return false;
    return comp.count(key1.first) > 0 && comp.count(key2.first) > 0;
}

void contamFinder3::set_cell_composition_overrides(
    const map<unsigned long, map<int, double> >& overrides){
    cell_composition_overrides = overrides;
    clear_data();
    compile_data(this->assn, this->indv_allelecounts);
}

/**
 * Populates internal data structures with data and expected values.
 */
void contamFinder3::compile_data(robin_hood::unordered_map<unsigned long, int>& assn,
     robin_hood::unordered_map<unsigned long, 
        map<pair<int, int>, map<pair<int, int>, pair<float, float> > > >& indv_allelecounts){

    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != 
        assn.end(); ++a){
        
        vector<double> n;
        vector<double> k;
        vector<double> p_e;
        vector<double> p_A;
        vector<double> p_B;
        vector<pair<int, int> > type1;
        vector<pair<int, int> > type2; 

        this->get_reads_expectations(a->first, a->second, indv_allelecounts[a->first],
            n, k, p_e, p_A, p_B, type1, type2);
        
        for (int i = 0; i < n.size(); ++i){
            int idx_all = n_all.size();

            n_all.push_back(n[i]);
            k_all.push_back(k[i]);
            p_e_all.push_back(p_e[i]);
            p_A_all.push_back(p_A[i]);
            p_B_all.push_back(p_B[i]);
            type1_all.push_back(type1[i]);
            type2_all.push_back(type2[i]);

            if (expfrac_to_idx.count(type1[i]) == 0){
                map<pair<int, int>, vector<int> > m;
                expfrac_to_idx.insert(make_pair(type1[i], m));
            }
            if (expfrac_to_idx[type1[i]].count(type2[i]) == 0){
                vector<int> v;
                expfrac_to_idx[type1[i]].insert(make_pair(type2[i], v));
            }
            expfrac_to_idx[type1[i]][type2[i]].push_back(idx_all);
            
            idx_to_cell.insert(make_pair(idx_all, a->first));
            
            if (cell_to_idx.count(a->first) == 0){
                vector<int> v;
                cell_to_idx.insert(make_pair(a->first, v));
            }
            cell_to_idx[a->first].push_back(idx_all);
            
            pair<int, int> sitekey;
            if (type2[i].second == -1){
                sitekey = make_pair(type1[i].second, -1);
            }
            else{
                int thistype1 = type1[i].second;
                int thistype2 = type2[i].second;
                if (thistype2 < thistype1){
                    int tmp = thistype2;
                    thistype2 = thistype1;
                    thistype1 = tmp;
                }
                sitekey = make_pair(thistype1, thistype2);
            }
            if (sitecomb_type_to_idx.count(sitekey) == 0){
                vector<int> v;
                sitecomb_type_to_idx.insert(make_pair(sitekey, v));
            }
            sitecomb_type_to_idx[sitekey].push_back(idx_all);
        }
    }
}

void contamFinder3::clear_data(){
    n_all.clear();
    k_all.clear();
    p_e_all.clear();
    p_A_all.clear();
    p_B_all.clear();
    type1_all.clear();
    type2_all.clear();
    expfrac_to_idx.clear();
    idx_to_cell.clear();
    cell_to_idx.clear();
    sitecomb_type_to_idx.clear();
}

/**
 * Helper function for compile_data.
 *
 * Given an assigned identity for a cell, collects all reads that should contain
 * each different type of allelic state and organize data in a way that is amenable
 * to solving for the parameters.
 */
void contamFinder3::get_reads_expectations(unsigned long barcode,
    int ident,
    map<pair<int, int>, map<pair<int, int>, pair<float, float> > >& allelecounts,
    vector<double>& n,
    vector<double>& k,
    vector<double>& p_e,
    vector<double>& p_A,
    vector<double>& p_B,
    vector<pair<int, int> >& type1,
    vector<pair<int, int> >& type2){
    
    static pair<int, int> nullkey = make_pair(-1, -1);

    if (has_composition_override(barcode)){
        const map<int, double>& comp = cell_composition_overrides[barcode];
        for (auto ac = allelecounts.begin(); ac != allelecounts.end(); ++ac){
            for (auto ac2 = ac->second.begin(); ac2 != ac->second.end(); ++ac2){
                if (!composition_row_is_relevant(comp, ac->first, ac2->first)){
                    continue;
                }
                double ref = ac2->second.first;
                double alt = ac2->second.second;
                if (ref + alt <= 0) continue;

                double expected_raw = composition_expected_from_row(comp, ac->first, ac2->first);
                // Weighted species-composition identities do not currently have
                // a generalized multi-r expression model.  Use the biologically
                // specified dosage composition as p_e and estimate only c.
                n.push_back(ref + alt);
                k.push_back(alt);
                p_e.push_back(expected_raw);
                p_A.push_back(-1.0);
                p_B.push_back(-1.0);
                type1.push_back(ac->first);
                type2.push_back(ac2->first);
            }
        }
        return;
    }

    bool is_combo = false;
    pair<int, int> combo;
    if (ident >= n_samples){
        combo = idx_to_hap_comb(ident, n_samples);
        is_combo = true;
    }
    for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator ac = 
        allelecounts.begin(); ac != allelecounts.end(); ++ac){
        if (is_combo){
            if (ac->first.first == combo.first){
                for (map<pair<int, int>, pair<float, float> >::iterator ac2 = ac->second.begin();
                    ac2 != ac->second.end(); ++ac2){
                    if (ac2->first.first == combo.second){
                        
                        double ref = ac2->second.first;
                        double alt = ac2->second.second;
                        if (ref + alt > 0){

                            // Two-component p_e (still needed for ambient profile estimation)
                            double expected = (double)(ac->first.second + ac2->first.second) / 4.0;
                            // Three-component per-genome expectations
                            double this_p_A = (double)ac->first.second / 2.0;
                            double this_p_B = (double)ac2->first.second / 2.0;
                            
                            n.push_back(ref+alt);
                            k.push_back(alt);
                            p_e.push_back(expected);
                            p_A.push_back(this_p_A);
                            p_B.push_back(this_p_B);
                            type1.push_back(ac->first);
                            type2.push_back(ac2->first);
                        }
                    }
                    else if (ac2->first.first > combo.second){
                        break;
                    }
                    else{
                        continue;
                    }
                }
            }
            else if (ac->first.first > combo.first){
                break;
            }
            else{
                continue;
            }
        }
        else{
            if (ac->first.first == ident){
                 // Store expectations without error incorporated (for now)
                double expected = (double)ac->first.second / 2.0;
                double ref = ac->second[nullkey].first;
                double alt = ac->second[nullkey].second;
                if (ref+alt > 0){
                    n.push_back(ref+alt);
                    k.push_back(alt);
                    p_e.push_back(expected);
                    // Diploid singlet: p_A/p_B not applicable, mark as -1
                    p_A.push_back(-1.0);
                    p_B.push_back(-1.0);
                    type1.push_back(ac->first);
                    type2.push_back(nullkey);
                }
            }
            else if (ac->first.first > ident){
                break;
            }
            else{
                continue;
            }
        }
    }
}

/**
 * When starting out, we don't know what the initial guess of contamination rate
 * should be. This is a hard problem, since we also don't know what the contamination
 * profile looks like yet.
 *
 * For each type of SNP, we can estimate the minimum bound on contamination rate
 *   in a given cell by assuming that ambient RNA has an alternate allele frequency
 *   as far as possible in the direction needed to "fix" the measured frequency.
 *
 *   In other words, if a cell belongs to individual 1 and we measure an alt allele
 *     frequency of 0.9 in sites where individual 1 is heterozygous alt, a minimum bound
 *     on the contamination rate could be found by assuming ambient RNA has 100% reference
 *     alleles at these sites (and c = 0.1). 
 *
 *   If a cell has alt allele frequency 0.67 at sites where the individual is heterozygous,
 *     then we can imagine ambient RNA is 100% alt alleles at these sites.
 *
 *     We then compute c = (measured-expected)/(ambient-expected), giving c = 0.34.
 *
 *     This gives us many measurements of c, and we return a weighted average.
 */
double contamFinder3::est_min_c(){
    map<pair<int, int>, map<pair<int, int>, double> > mean_f;
    map<pair<int, int>, map<pair<int, int>, double> > counts;
    
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != 
        assn.end(); ++a){
        
        if (a->second >= n_samples){
            pair<int, int> combo = idx_to_hap_comb(a->second, n_samples);
            for (int nalt1 = 0; nalt1 <= 2; ++nalt1){
                pair<int, int> key1 = make_pair(combo.first, nalt1);
                for (int nalt2 = 0; nalt2 <= 2; ++nalt2){
                    // Tetraploid-aware: skip weak categories (p_e = 0.25 or 0.75)
                    if (tetraploid_aware && !category_passes_hard_filter(true, nalt1, nalt2)){
                        continue;
                    }
                    pair<int, int> key2 = make_pair(combo.second, nalt2);
                    double ref = indv_allelecounts[a->first][key1][key2].first;
                    double alt = indv_allelecounts[a->first][key1][key2].second;
                    if (ref + alt == 0){
                        continue;
                    }
                    double f = alt/(ref+alt);
                    double p0 = (double)(nalt1 + nalt2)/4.0;
                    if (mean_f.count(key1) == 0){
                        map<pair<int, int>, double> m;
                        mean_f.insert(make_pair(key1, m));
                        counts.insert(make_pair(key1, m));
                    }
                    if (mean_f[key1].count(key2) == 0){
                        mean_f[key1].insert(make_pair(key2, 0));
                        counts[key1].insert(make_pair(key2, 0));
                    }
                    if (!weighted){
                        mean_f[key1][key2] += f;    
                        counts[key1][key2]++;
                    }
                    else{
                        mean_f[key1][key2] += f * assn_llr[a->first];
                        counts[key1][key2] += assn_llr[a->first];
                    }
                }
            }
        }
        else{
            pair<int, int> nullkey = make_pair(-1, -1);
            for (int nalt = 0; nalt <= 2; ++nalt){
                pair<int, int> key = make_pair(a->second, nalt);
                double ref = indv_allelecounts[a->first][key][nullkey].first;
                double alt = indv_allelecounts[a->first][key][nullkey].second;
                if (ref + alt == 0){
                    continue;
                }
                double f = alt/(ref+alt);
                double p0 = (double)nalt/2.0;
                if (mean_f.count(key) == 0){
                    map<pair<int, int>, double> m;
                    mean_f.insert(make_pair(key, m));
                    counts.insert(make_pair(key, m));
                }
                if (mean_f[key].count(nullkey) == 0){
                    mean_f[key].insert(make_pair(nullkey, 0));
                    counts[key].insert(make_pair(nullkey, 0));
                }
                if (!weighted){
                    mean_f[key][nullkey] += f;
                    counts[key][nullkey]++;
                }
                else{
                    mean_f[key][nullkey] += f*assn_llr[a->first];
                    counts[key][nullkey] += assn_llr[a->first];
                }
            }
        }
    }

    double minc_mean = 0.0;
    double minc_weightsum = 0.0;
    
    double minc_min = -1;
    double minc_max = -1;
    
    map<int, double> minc_by_id;
    map<int, double> minc_by_id_count;

    for (map<pair<int, int>, map<pair<int, int>, double> >::iterator d = mean_f.begin();
        d != mean_f.end(); ++d){
        for (map<pair<int, int>, double>::iterator d2 = d->second.begin(); d2 != d->second.end();
            ++d2){
            double expec;
            if (d2->first.first == -1){
                expec = d->first.second / 2.0;
            }
            else{
                expec = (d->first.second + d2->first.second)/4.0;
            }
            expec = adjust_p_err(expec, e_r, e_a);
            if (counts[d->first][d2->first] > 0){
               
                double inf = d2->second / counts[d->first][d2->first];
                double minc;
                if (inf > expec){
                    // treat p_c = 1
                    minc = (inf - expec)/(1.0 - expec);
                }   
                else{
                    // treat p_c = 0
                    minc = (inf - expec)/(0.0 - expec);
                }
                if (minc_min == -1 || minc < minc_min){
                    minc_min = minc;
                }
                if (minc_max == -1 || minc > minc_max){
                    minc_max = minc;
                }
                minc_mean += minc * counts[d->first][d2->first];
                minc_weightsum += counts[d->first][d2->first];
                
                if (minc_by_id.count(d->first.first) == 0){
                    minc_by_id.insert(make_pair(d->first.first, 0));
                    minc_by_id_count.insert(make_pair(d->first.first, 0));
                }
                
                if (d2->first.first != -1){
                    double w = 0.5*(double)counts[d->first][d2->first];
                    minc_by_id[d->first.first] += w*minc;
                    minc_by_id_count[d->first.first] += w;
                    if (minc_by_id.count(d2->first.first) == 0){
                        minc_by_id.insert(make_pair(d2->first.first, 0));
                        minc_by_id_count.insert(make_pair(d2->first.first, 0));
                    }
                    minc_by_id[d2->first.first] += w*minc;
                    minc_by_id_count[d2->first.first] += w;
                }
                else{
                    double w = (double)counts[d->first][d2->first];
                    minc_by_id[d->first.first] += w*minc;
                    minc_by_id_count[d->first.first] += w;
                }
            }
        }
    }
    
    double vsum = 0.0;
    double vcount = 0.0;
    for (map<int, double>::iterator mi = minc_by_id.begin(); mi != minc_by_id.end(); ++mi){
        vsum += mi->second/minc_by_id_count[mi->first];
        vcount++;
    }
    // Guard against divide-by-zero when only 0 or 1 individuals present.
    // The (vcount-1) denominator is meant to exclude "self" contribution,
    // but with <= 1 individual we fall back to a simple average.
    double c_est;
    if (vcount <= 1.0){
        c_est = (vcount > 0) ? vsum / vcount : 0.01;
    }
    else{
        c_est = vsum / (vcount - 1.0);
    }
    
    contam_prof.clear();
    double minval = 0.01;
    double denom = 0.0;

    // Build set of active individuals for filtering
    set<int> active_set(idx2samp.begin(), idx2samp.end());

    // Guard: if c_est is zero or negative, all fracs become degenerate.
    // Fall back to equal proportions.
    if (c_est <= 0.0 || isnan(c_est) || isinf(c_est)){
        for (map<int, double>::iterator mi = minc_by_id.begin(); mi != minc_by_id.end(); ++mi){
            if (active_set.count(mi->first) == 0) continue;
            double frac = 1.0 / (double)idx2samp.size();
            denom += frac;
            contam_prof.insert(make_pair(mi->first, frac));
        }
    }
    else{
        for (map<int, double>::iterator mi = minc_by_id.begin(); mi != minc_by_id.end(); ++mi){
            if (active_set.count(mi->first) == 0) continue;
            double val = mi->second/minc_by_id_count[mi->first];
            double frac = 1.0 - val/c_est;
            if (frac < minval){
                frac = minval;
            }
            denom += frac;
            contam_prof.insert(make_pair(mi->first, frac));
        }
    }
    // Fill any idx2samp individuals not in minc_by_id with minimum fraction
    for (int i = 0; i < (int)idx2samp.size(); i++){
        int samp = idx2samp[i];
        if (contam_prof.count(samp) == 0){
            contam_prof.insert(make_pair(samp, minval));
            denom += minval;
        }
    }
    if (inter_species){
        contam_prof.insert(make_pair(-1, 1.0/((double)idx2samp.size() + 1.0)));
        denom += contam_prof[-1];
    }
    for (map<int, double>::iterator cp = contam_prof.begin(); cp != contam_prof.end(); ++cp){
        cp->second /= denom;
    }
    
    return c_est;

    /*
    minc_mean /= minc_weightsum;
    fprintf(stderr, "mean %f\n", minc_mean);
    exit(0);

    return minc_mean;
    */
}

/**
 * When running the first time, create a starting guess of contamination profile
 * and infer the likeliest mixture of individuals.
 */
double contamFinder3::init_params(double& init_c){
    
    if (!contam_prof_initialized){
        // Init contam prof
        //contam_prof.clear();
        
        // This will compute an initial contamination profile based on the counts of different
        // types of cells.

        // The commented-out block will compute an initial contamination profile where
        // all proportions are equal.
        /*
        map<int, int> acounts;
        int atots = 0;
        for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
            a != assn.end(); ++a){
            if (a->second >= n_samples){
                pair<int, int> comb = idx_to_hap_comb(a->second, n_samples);
                if (acounts.count(comb.first) == 0){
                    acounts.insert(make_pair(comb.first, 0));
                }
                if (acounts.count(comb.second) == 0){
                    acounts.insert(make_pair(comb.second, 0));
                }
                acounts[comb.first]++;
                acounts[comb.second]++;
                atots += 2;
            }
            else{
                if (acounts.count(a->second) == 0){
                    acounts.insert(make_pair(a->second, 0));
                }
                acounts[a->second] += 2;
                atots += 2;
            }
        }
        double minfrac = 0.001;
        double proptot = 0.0;
        for (int i = 0; i < idx2samp.size(); ++i){
            int samp = idx2samp[i];
            double prop = (double)acounts[samp]/(double)atots;
            if (prop < minfrac){
                prop = minfrac;
            }
            contam_prof.insert(make_pair(samp, prop));
            proptot += prop;
        }
        if (inter_species){
            double prop = 1.0/(double)(idx2samp.size() + 1);
            contam_prof.insert(make_pair(-1, prop));
            proptot += prop;
        }
        for (map<int, double>::iterator cp = contam_prof.begin(); cp != contam_prof.end();
            ++cp){
            cp->second /= proptot;
        }
        */
        
        /*
        for (int i = 0; i < idx2samp.size(); ++i){
            int samp = idx2samp[i];
            if (inter_species){
                contam_prof.insert(make_pair(samp, 1.0/(double)(idx2samp.size() + 1)));
            }
            else{
                contam_prof.insert(make_pair(samp, 1.0/(double)idx2samp.size()));
            }
        }
        if (inter_species){
            contam_prof.insert(make_pair(-1, 1.0/(double)(idx2samp.size() + 1)));
        }
        */

    }
    return update_amb_prof_mixture(true, init_c, true);
}

/**
 * Produces a single, global contamination rate estimate.
 * Stores this as contam_cell_prior.
 */
void contamFinder3::est_contam_cells_global(){
    
    contam_rate.clear();
    contam_rate_se.clear();
    
    vector<double> n;
    vector<double> k;
    vector<double> p_e;
    vector<double> p_c;
    
    for (map<unsigned long, vector<int> >::iterator ci = cell_to_idx.begin(); 
        ci != cell_to_idx.end(); ++ci){
        
        for (vector<int>::iterator i = ci->second.begin(); i != ci->second.end(); ++i){
            n.push_back(n_all[*i]);
            k.push_back(k_all[*i]);
            p_e.push_back(adjust_p_err(p_e_all[*i], e_r, e_a));
            p_c.push_back(amb_mu[type1_all[*i]][type2_all[*i]]);
        }
    }
    optimML::brent_solver c_global(ll_c, dll_dc, d2ll_dc2);
    if (num_threads > 1){
        c_global.set_threads(num_threads);
    }
    c_global.add_data("n", n);
    c_global.add_data("k", k);
    c_global.add_data("p_e", p_e);
    c_global.add_data("p_c", p_c);
    c_global.constrain_01();
    c_global.set_maxiter(-1);
    try {
        contam_cell_prior = c_global.solve(0,1);
    } catch (...) {
        fprintf(stderr, "WARNING: global contamination rate solver failed; using initial estimate\n");
        contam_cell_prior = c_init > 0 ? c_init : 0.15;
    }
    contam_cell_prior_var = -1;
    fprintf(stderr, "MLE global contamination rate: %f\n", contam_cell_prior);
}

// ============================================================================
// LOO ambient profile: for each singlet identity, compute an ambient profile
// that excludes cells assigned to that identity. This prevents a cell's own
// signal from appearing in the ambient profile used to estimate its
// contamination rate, reducing circular bias.
// ============================================================================
void contamFinder3::compute_loo_profiles(){
    amb_mu_loo.clear();

    // Collect the set of singlet identities present in assignments
    set<int> singlet_ids;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){
        if (a->second < n_samples){
            singlet_ids.insert(a->second);
        }
    }

    // For each singlet identity, compute a leave-one-out ambient profile by
    // re-weighting the global contam_prof to exclude that identity's contribution.
    // amb_mu is a weighted sum over individuals:
    //   amb_mu[type1][type2] = sum_j { contam_prof[j] * expfracs[type1][j] }
    // The LOO version for identity i removes contam_prof[i] and re-normalizes.
    for (set<int>::iterator sid = singlet_ids.begin(); sid != singlet_ids.end(); ++sid){
        int excl = *sid;
        // Compute the renormalization factor: sum of contam_prof excluding excl
        double sum_other = 0.0;
        for (map<int, double>::iterator cp = contam_prof.begin();
            cp != contam_prof.end(); ++cp){
            if (cp->first != excl && cp->first >= 0){
                sum_other += cp->second;
            }
        }
        // If removing this identity leaves nothing, fall back to global profile
        if (sum_other < 1e-12){
            amb_mu_loo[excl] = amb_mu;
            continue;
        }
        double renorm = 1.0 / sum_other;

        // Build LOO ambient profile for all (type1, type2) categories
        map<pair<int,int>, map<pair<int,int>, double>> loo_prof;
        for (map<pair<int,int>, map<pair<int,int>, double>>::iterator t1 = amb_mu.begin();
            t1 != amb_mu.end(); ++t1){
            for (map<pair<int,int>, double>::iterator t2 = t1->second.begin();
                t2 != t1->second.end(); ++t2){
                // Recompute from contam_prof * expfracs, excluding excl
                double val = 0.0;
                if (expfracs.count(t1->first) > 0){
                    for (map<int, double>::iterator cp = contam_prof.begin();
                        cp != contam_prof.end(); ++cp){
                        if (cp->first == excl || cp->first < 0) continue;
                        if (expfracs[t1->first].count(cp->first) > 0){
                            val += (cp->second * renorm) * expfracs[t1->first][cp->first];
                        }
                    }
                }
                // If we got a valid value, use it; otherwise fall back to global
                if (val > 1e-15 && val < 1.0 - 1e-15){
                    loo_prof[t1->first][t2->first] = val;
                } else {
                    loo_prof[t1->first][t2->first] = t2->second;
                }
            }
        }
        amb_mu_loo[excl] = loo_prof;
    }
    fprintf(stderr, "  Computed LOO ambient profiles for %lu singlet identities\n",
        singlet_ids.size());
}

/**
 * Once an ambient RNA profile exists (we have estimates of p_c parameters for
 * every category of SNP, we can re-estimate the likeliest contamination rate per
 * cell given these parameters.
 *
 * Uses an Emprical Bayes method -- a prior distribution on contamination rate
 * per cell is based on the data set-wide distribution from the last round of 
 * estimate.
 *
 */

contamFinder3::CellContamFitResult contamFinder3::fit_one_contam_cell(
    unsigned long barcode,
    const vector<int>& obs_idx){

    CellContamFitResult out;
    out.barcode = barcode;
    out.c = 1.0;
    out.c_se = 0.0;
    out.ll = -INFINITY;
    out.r = 0.5;
    out.r_se = 0.0;
    out.is_heterotypic = false;
    out.bfgs_fallback = false;
    out.has_allele_ratio = false;

    // Determine this cell's identity for LOO lookup.
    int cell_ident = -1;
    if (use_loo && assn.count(barcode) > 0){
        cell_ident = assn[barcode];
        // For doublets, LOO is not applied (no single identity to exclude).
        if (cell_ident >= n_samples){
            cell_ident = -1;
        }
    }
    bool have_loo = (cell_ident >= 0 && amb_mu_loo.count(cell_ident) > 0);

    // Determine if this cell is a heterotypic combo eligible for the
    // three-component model.
    bool is_heterotypic = false;
    if (assn.count(barcode) > 0 && assn[barcode] >= n_samples){
        pair<int, int> combo = idx_to_hap_comb(assn[barcode], n_samples);
        if (combo.first != combo.second){
            is_heterotypic = true;
        }
    }
    out.is_heterotypic = is_heterotypic;

    vector<double> n;
    vector<double> k;
    vector<double> p_e;
    vector<double> p_c;
    vector<double> cell_p_A;
    vector<double> cell_p_B;
    vector<double> fi_weights;

    auto get_pc = [&](int obs_index) -> double {
        const pair<int,int>& t1 = type1_all[obs_index];
        const pair<int,int>& t2 = type2_all[obs_index];
        if (have_loo){
            auto loo_id = amb_mu_loo.find(cell_ident);
            if (loo_id != amb_mu_loo.end()){
                auto a1 = loo_id->second.find(t1);
                if (a1 != loo_id->second.end()){
                    auto a2 = a1->second.find(t2);
                    if (a2 != a1->second.end()){
                        return a2->second;
                    }
                }
            }
        }
        auto a1 = amb_mu.find(t1);
        if (a1 != amb_mu.end()){
            auto a2 = a1->second.find(t2);
            if (a2 != a1->second.end()){
                return a2->second;
            }
        }
        return 0.5;
    };

    // Compile data for this cell.
    for (vector<int>::const_iterator i = obs_idx.begin(); i != obs_idx.end(); ++i){
        double this_p_e = adjust_p_err(p_e_all[*i], e_r, e_a);
        double this_p_c = get_pc(*i);

        double this_p_A = -1.0;
        double this_p_B = -1.0;
        if (is_heterotypic && p_A_all[*i] >= 0){
            this_p_A = adjust_p_err(p_A_all[*i], e_r, e_a);
            this_p_B = adjust_p_err(p_B_all[*i], e_r, e_a);
        }

        if (use_fi_weight){
            double w = compute_fi_weight(this_p_e, this_p_c);
            if (w < 1e-12) continue;
            n.push_back(n_all[*i]);
            k.push_back(k_all[*i]);
            p_e.push_back(this_p_e);
            p_c.push_back(this_p_c);
            cell_p_A.push_back(this_p_A);
            cell_p_B.push_back(this_p_B);
            fi_weights.push_back(w);
        } else {
            if (!is_heterotypic && tetraploid_aware){
                bool is_combo_cat = (type2_all[*i].first != -1);
                if (is_combo_cat){
                    if (amb_mu_available){
                        if (!category_passes_adaptive_filter(this_p_e, this_p_c, min_signal_gap)){
                            continue;
                        }
                    } else {
                        if (!category_passes_hard_filter(true, type1_all[*i].second, type2_all[*i].second)){
                            continue;
                        }
                    }
                }
            }
            n.push_back(n_all[*i]);
            k.push_back(k_all[*i]);
            p_e.push_back(this_p_e);
            p_c.push_back(this_p_c);
            cell_p_A.push_back(this_p_A);
            cell_p_B.push_back(this_p_B);
        }
    }

    // If filtering removed all data for this cell, fall back to unfiltered.
    if (n.empty() && (tetraploid_aware || use_fi_weight)){
        fi_weights.clear();
        cell_p_A.clear();
        cell_p_B.clear();
        for (vector<int>::const_iterator i = obs_idx.begin(); i != obs_idx.end(); ++i){
            double this_p_e = adjust_p_err(p_e_all[*i], e_r, e_a);
            double this_p_c = get_pc(*i);
            n.push_back(n_all[*i]);
            k.push_back(k_all[*i]);
            p_e.push_back(this_p_e);
            p_c.push_back(this_p_c);
            double this_p_A = -1.0;
            double this_p_B = -1.0;
            if (is_heterotypic && p_A_all[*i] >= 0){
                this_p_A = adjust_p_err(p_A_all[*i], e_r, e_a);
                this_p_B = adjust_p_err(p_B_all[*i], e_r, e_a);
            }
            cell_p_A.push_back(this_p_A);
            cell_p_B.push_back(this_p_B);
            if (use_fi_weight){
                fi_weights.push_back(compute_fi_weight(this_p_e, this_p_c));
            }
        }
    }

    double c_cell_map = 1.0;
    double se = 0.0;
    double ll = -INFINITY;
    double r_cell_map = 0.5;
    double r_se = 0.0;

    if (is_heterotypic && !cell_p_A.empty() && cell_p_A[0] >= 0){
        double c_init_cell = (contam_cell_prior > 0) ? contam_cell_prior : 0.15;

        vector<double> params_init = {0.5, c_init_cell};
        optimML::multivar_ml_solver solver(params_init, ll_three, dll_three);
        solver.add_data("n", n);
        solver.add_data("k", k);
        solver.add_data("p_A", cell_p_A);
        solver.add_data("p_B", cell_p_B);
        solver.add_data("p_c", p_c);
        solver.constrain_01(0);
        solver.constrain_01(1);
        if (use_fi_weight && !fi_weights.empty()){
            solver.add_weights(fi_weights);
        }
        if (contam_cell_prior > 0 && contam_cell_prior_var > 0){
            pair<double, double> bm = beta_moments(contam_cell_prior, contam_cell_prior_var);
            solver.add_beta_prior(1, bm.first, bm.second);
        }

        try {
            bool success = solver.solve();
            if (success){
                r_cell_map = solver.results[0];
                c_cell_map = solver.results[1];
                ll = solver.log_likelihood;
                se = 0.0;
                r_se = 0.0;
            }
        } catch (...){
            out.bfgs_fallback = true;
            r_cell_map = 0.5;
            optimML::brent_solver c_fallback(ll_c, dll_dc, d2ll_dc2);
            c_fallback.add_data("n", n);
            c_fallback.add_data("k", k);
            c_fallback.add_data("p_e", p_e);
            c_fallback.add_data("p_c", p_c);
            c_fallback.constrain_01();
            if (contam_cell_prior > 0 && contam_cell_prior_var > 0){
                pair<double, double> bm = beta_moments(contam_cell_prior, contam_cell_prior_var);
                c_fallback.add_beta_prior(bm.first, bm.second);
            }
            c_fallback.set_maxiter(-1);
            try {
                c_cell_map = c_fallback.solve(0, 1);
                if (c_fallback.root_found){
                    ll = c_fallback.log_likelihood;
                    if (c_fallback.se_found) se = c_fallback.se;
                }
            } catch (...) {
                c_cell_map = 1.0;
            }
        }

        bool try_extra_starts = thorough_multistart;
        if (adaptive_multistart && !try_extra_starts){
            bool failed_or_bad = !std::isfinite(ll);
            bool boundary_r = (r_cell_map < 0.05 || r_cell_map > 0.95);
            bool boundary_c = (c_cell_map < 0.01 || c_cell_map > 0.95);
            try_extra_starts = failed_or_bad || boundary_r || boundary_c;
        }

        if (try_extra_starts){
            for (double r_alt : {0.2, 0.8}){
                vector<double> alt_init = {r_alt, c_init_cell};
                optimML::multivar_ml_solver solver_alt(alt_init, ll_three, dll_three);
                solver_alt.add_data("n", n);
                solver_alt.add_data("k", k);
                solver_alt.add_data("p_A", cell_p_A);
                solver_alt.add_data("p_B", cell_p_B);
                solver_alt.add_data("p_c", p_c);
                solver_alt.constrain_01(0);
                solver_alt.constrain_01(1);
                if (use_fi_weight && !fi_weights.empty()){
                    solver_alt.add_weights(fi_weights);
                }
                if (contam_cell_prior > 0 && contam_cell_prior_var > 0){
                    pair<double, double> bm = beta_moments(contam_cell_prior, contam_cell_prior_var);
                    solver_alt.add_beta_prior(1, bm.first, bm.second);
                }
                try {
                    bool success = solver_alt.solve();
                    if (success && solver_alt.log_likelihood > ll){
                        r_cell_map = solver_alt.results[0];
                        c_cell_map = solver_alt.results[1];
                        ll = solver_alt.log_likelihood;
                    }
                } catch (...) {
                    // ignore this alternate start
                }
            }
        }

        out.has_allele_ratio = true;
        out.r = r_cell_map;
        out.r_se = r_se;
    } else {
        optimML::brent_solver c_cell(ll_c, dll_dc, d2ll_dc2);
        c_cell.add_data("n", n);
        c_cell.add_data("k", k);
        c_cell.add_data("p_e", p_e);
        c_cell.add_data("p_c", p_c);
        c_cell.constrain_01();
        if (use_fi_weight && !fi_weights.empty()){
            c_cell.add_weights(fi_weights);
        }
        if (contam_cell_prior > 0 && contam_cell_prior_var > 0){
            pair<double, double> bm = beta_moments(contam_cell_prior, contam_cell_prior_var);
            c_cell.add_beta_prior(bm.first, bm.second);
        }
        c_cell.set_maxiter(-1);
        try{
            c_cell_map = c_cell.solve(0,1);
            if (c_cell.root_found){
                ll = c_cell.log_likelihood;
                if (c_cell.se_found){
                    se = c_cell.se;
                }
            }
        } catch (...){
            // keep default c=1.0 and ll=-inf
        }
    }

    out.c = c_cell_map;
    out.c_se = se;
    out.ll = ll;
    return out;
}

/**
 * Once an ambient RNA profile exists (we have estimates of p_c parameters for
 * every category of SNP), re-estimate the likeliest contamination rate per cell.
 * The expensive per-cell optimizer work is embarrassingly parallel, so worker
 * threads write CellContamFitResult records and the class-owned maps are merged
 * serially after the OpenMP loop.
 */
void contamFinder3::est_contam_cells(){

    contam_rate.clear();
    contam_rate_se.clear();
    contam_rate_ll.clear();
    allele_ratio.clear();
    allele_ratio_se.clear();

    vector<pair<unsigned long, vector<int> > > cells;
    cells.reserve(cell_to_idx.size());
    for (map<unsigned long, vector<int> >::iterator ci = cell_to_idx.begin();
        ci != cell_to_idx.end(); ++ci){
        cells.push_back(*ci);
    }

    vector<CellContamFitResult> results(cells.size());
    int nt = (num_threads > 1) ? num_threads : 1;

    #pragma omp parallel for num_threads(nt) schedule(dynamic, 16)
    for (int i = 0; i < (int)cells.size(); ++i){
        results[i] = fit_one_contam_cell(cells[i].first, cells[i].second);
    }

    vector<double> cell_c_maps;
    vector<double> cell_c_llr;
    vector<double> twocomp_c_maps;
    vector<double> twocomp_c_llr;
    vector<double> cell_r_maps;
    cell_c_maps.reserve(results.size());
    twocomp_c_maps.reserve(results.size());

    int n_three_comp = 0;
    int n_two_comp = 0;
    int n_bfgs_fallback = 0;

    for (vector<CellContamFitResult>::const_iterator it = results.begin();
        it != results.end(); ++it){
        contam_rate.emplace(it->barcode, it->c);
        contam_rate_se.emplace(it->barcode, it->c_se);
        contam_rate_ll.emplace(it->barcode, it->ll);

        if (it->has_allele_ratio){
            allele_ratio.emplace(it->barcode, it->r);
            allele_ratio_se.emplace(it->barcode, it->r_se);
            cell_r_maps.push_back(it->r);
        }
        if (it->is_heterotypic){
            n_three_comp++;
        } else {
            n_two_comp++;
        }
        if (it->bfgs_fallback){
            n_bfgs_fallback++;
        }

        // Keep cell_c_maps and cell_c_llr (and their twocomp counterparts)
        // in lockstep: in weighted mode, only push when a valid weight exists;
        // in unweighted mode, push unconditionally.
        if (weighted){
            if (assn.count(it->barcode) > 0){
                int aid = assn[it->barcode];
                if (id_llrsum.count(aid) > 0 && id_llrsum[aid] != 0.0){
                    double weight = assn_llr[it->barcode] / id_llrsum[aid];
                    cell_c_maps.push_back(it->c);
                    cell_c_llr.push_back(weight);
                    if (!it->is_heterotypic){
                        twocomp_c_maps.push_back(it->c);
                        twocomp_c_llr.push_back(weight);
                    }
                }
            }
        } else {
            cell_c_maps.push_back(it->c);
            if (!it->is_heterotypic){
                twocomp_c_maps.push_back(it->c);
            }
        }
    }

    // Re-compute data set-wide distribution (all cells, for reporting).
    pair<double, double> mu_var;
    if (weighted && !cell_c_llr.empty()){
        mu_var = welford_weights(cell_c_maps, cell_c_llr, false);
    }
    else{
        mu_var = welford(cell_c_maps);
    }

    // Compute prior from two-component cells only (singlets + homotypic),
    // to avoid heterotypic ridge-variance inflating the prior.
    pair<double, double> prior_mu_var;
    bool use_twocomp_prior = (twocomp_c_maps.size() >= 20);
    if (use_twocomp_prior){
        if (weighted && !twocomp_c_llr.empty()){
            prior_mu_var = welford_weights(twocomp_c_maps, twocomp_c_llr, false);
            if (prior_mu_var.second < 1e-3){
                prior_mu_var = welford(twocomp_c_maps);
            }
        } else {
            prior_mu_var = welford(twocomp_c_maps);
        }
    } else {
        prior_mu_var = mu_var;
    }

    if (weighted && mu_var.second < 1e-3){
        mu_var = welford(cell_c_maps);
    }
    if (mu_var.second > 1e-6){
        if (contam_cell_prior > 0 && contam_cell_prior_var > 0){
            fprintf(stderr, "Shrunken per-cell contamination rates:\n");
        }
        else{
            fprintf(stderr, "Per-cell contamination rates:\n");
        }
        fprintf(stderr, "  Mean: %f Std dev: %f\n", mu_var.first, sqrt(mu_var.second));
        fprintf(stderr, "  Three-component cells: %d  Two-component cells: %d\n",
            n_three_comp, n_two_comp);
        if (n_bfgs_fallback > 0){
            fprintf(stderr, "  BFGS fallbacks to Brent: %d / %d heterotypic cells\n",
                n_bfgs_fallback, n_three_comp);
        }
        if (!cell_r_maps.empty()){
            pair<double, double> r_stats = welford(cell_r_maps);
            fprintf(stderr, "  Allele ratio (r) mean: %f Std dev: %f\n",
                r_stats.first, sqrt(r_stats.second));
        }
        if (!user_prior_set){
            contam_cell_prior = prior_mu_var.first;
            contam_cell_prior_var = prior_mu_var.second;
            if (use_twocomp_prior && n_three_comp > 0){
                fprintf(stderr, "  Prior from two-component cells only (%lu cells): "
                    "mean=%f var=%f\n", twocomp_c_maps.size(),
                    prior_mu_var.first, prior_mu_var.second);
            }
        }
    }

    map<int, double> idcsum;
    map<int, double> idccount;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){
        if (contam_rate.count(a->first) == 0) continue;
        if (a->second >= n_samples){
            pair<int, int> comb = idx_to_hap_comb(a->second, n_samples);
            if (idcsum.count(comb.first) == 0){
                idcsum.insert(make_pair(comb.first, 0));
                idccount.insert(make_pair(comb.first,0));
            }
            if (idcsum.count(comb.second) == 0){
                idcsum.insert(make_pair(comb.second, 0));
                idccount.insert(make_pair(comb.second,0));
            }
            idcsum[comb.first] += 0.5*contam_rate[a->first];
            idcsum[comb.second] += 0.5*contam_rate[a->first];
            idccount[comb.first] += 0.5;
            idccount[comb.second] += 0.5;
        }
        else{
            if (idcsum.count(a->second) == 0){
                idcsum.insert(make_pair(a->second, 0));
                idccount.insert(make_pair(a->second, 0));
            }
            idcsum[a->second] += contam_rate[a->first];
            idccount[a->second] += 1.0;
        }
    }
}

/**
 * Updates ambient RNA profile based on current cell-contamination estimates,
 * without consideration of the fraction of individuals making up ambient
 * RNA.
 *
 * Currently not used.
 */
double contamFinder3::update_ambient_profile(bool global_c){
    
    vector<pair<int, int> > idx2expfrac1;
    vector<pair<int, int> > idx2expfrac2;
    
    vector<double> n;
    vector<double> k;
    vector<double> p_e;
    vector<double> p_c;
    vector<int> efp_idx;
    vector<double> c;
    vector<double> weights;

    map<int, int> idx2idx;
    vector<double> efparams;
    map<int, pair<int, int> > efp2ef1;
    map<int, pair<int, int> > efp2ef2;

    double floor = 1e-3;
    double ceil = 1-1e-3;

    for (map<pair<int, int>, map<pair<int, int>, vector<int> > >::iterator ei1 = 
        expfrac_to_idx.begin(); ei1 != expfrac_to_idx.end(); ++ei1){
        for (map<pair<int, int>, vector<int> >::iterator ei2 = ei1->second.begin(); ei2 != 
            ei1->second.end(); ++ei2){
            
            int ef_param_idx = efparams.size();
            double p_c = amb_mu[ei1->first][ei2->first];
            if (p_c < floor){
                p_c = floor;
            }
            else if (p_c > ceil){
                p_c = ceil;
            }
            efparams.push_back(p_c);
            efp2ef1.insert(make_pair(ef_param_idx, ei1->first));
            efp2ef2.insert(make_pair(ef_param_idx, ei2->first));
            
            for (vector<int>::iterator i = ei2->second.begin(); i != ei2->second.end(); ++i){
                if (global_c || (idx_to_cell.count(*i) > 0 && 
                    contam_rate.count(idx_to_cell[*i]) > 0)){ 
                    
                    idx2idx.insert(make_pair(*i, n.size()));    
                    n.push_back(n_all[*i]);
                    k.push_back(k_all[*i]);
                    p_e.push_back(adjust_p_err(p_e_all[*i], e_r, e_a));
                    efp_idx.push_back(ef_param_idx);
                    if (global_c){
                        c.push_back(contam_cell_prior);
                    }
                    else{
                        c.push_back(contam_rate[idx_to_cell[*i]]);
                    }
                    double weight = 1.0;
                    if (weighted){
                        double llr = assn_llr[idx_to_cell[*i]];
                        double llrtot = id_llrsum[assn[idx_to_cell[*i]]];
                        weight = llr / llrtot;
                    }
                    weights.push_back(weight);
                }
            }        
        }
    }
        
    optimML::multivar_ml_solver solver(efparams, ll_ambmu, dll_ambmu);
    if (num_threads > 1){
        solver.set_threads(num_threads);
        // Many parameters - use threads
        solver.set_bfgs_threads(num_threads);
    }
    for (int i = 0; i < efparams.size(); ++i){
        solver.constrain_01(i);
    }
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("p_e", p_e);
    solver.add_data("ef_idx", efp_idx);
    solver.add_data("c", c);
    solver.add_weights(weights);
    
    try {
        solver.solve();
    } catch (...) {
        fprintf(stderr, "WARNING: ambient profile direct solver failed; keeping current profile\n");
        return -1e30;
    }
    
    for (int i = 0; i < (int)efparams.size(); ++i){
        double updated = solver.results[i];
        if (amb_mu.count(efp2ef1[i]) == 0){
            map<pair<int, int>, double> m;
            amb_mu.insert(make_pair(efp2ef1[i], m));
        }
        if (amb_mu[efp2ef1[i]].count(efp2ef2[i]) == 0){
            amb_mu[efp2ef1[i]].insert(make_pair(efp2ef2[i], 0.0));
        }
        amb_mu[efp2ef1[i]][efp2ef2[i]] = updated;
    }
    return solver.log_likelihood;
}

/**
 * Compiles data for update_amb_prof_mixture(). Puts necessary values
 * in vectors, which can then be bootstrapped.
 */
void contamFinder3::compile_amb_prof_dat(bool solve_for_c, 
    bool use_global_c,
    vector<vector<double> >& mixfracs,
    vector<double>& weights,
    vector<double>& n,
    vector<double>& k,
    vector<double>& p_e,
    vector<double>& c){
    
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != 
        assn.end(); ++a){
        
        bool is_comb = false;
        pair<int, int> comb;
        if (a->second >= n_samples){
            is_comb = true;
            comb = idx_to_hap_comb(a->second, n_samples);
        }
        
        // get weight
        double weight = 1.0;
        if (weighted){
            weight = assn_llr[a->first] / id_llrsum[a->second];
        }

        double this_c;
        if (!solve_for_c){
            if (use_global_c){
                this_c = contam_cell_prior;
            }
            else{
                this_c = contam_rate[a->first];
            }
        }

        if (has_composition_override(a->first)){
            const map<int, double>& comp = cell_composition_overrides[a->first];
            for (auto ac1 = indv_allelecounts[a->first].begin();
                ac1 != indv_allelecounts[a->first].end(); ++ac1){

                for (auto ac2 = ac1->second.begin(); ac2 != ac1->second.end(); ++ac2){
                    if (!composition_row_is_relevant(comp, ac1->first, ac2->first)){
                        continue;
                    }

                    double ref = ac2->second.first;
                    double alt = ac2->second.second;
                    if (ref + alt <= 0) continue;

                    double expected = adjust_p_err(
                        composition_expected_from_row(comp, ac1->first, ac2->first),
                        e_r, e_a);

                    if (tetraploid_aware && comp.size() > 1){
                        if (amb_mu_available){
                            double this_p_c = 0.0;
                            if (amb_mu.count(ac1->first) > 0 &&
                                amb_mu[ac1->first].count(ac2->first) > 0){
                                this_p_c = amb_mu[ac1->first][ac2->first];
                            }
                            if (!category_passes_adaptive_filter(expected, this_p_c, min_signal_gap)){
                                continue;
                            }
                        } else {
                            // Use the original hard filter on pairwise rows.  For
                            // weighted multi-species identities these rows are an
                            // approximation until native triple-count artifacts exist.
                            if (!category_passes_hard_filter(true, ac1->first.second, ac2->first.second)){
                                continue;
                            }
                        }
                    }

                    n.push_back(ref + alt);
                    k.push_back(alt);
                    weights.push_back(weight);
                    if (!solve_for_c){
                        c.push_back(this_c);
                    }
                    p_e.push_back(expected);

                    vector<double> mixfrac_row;
                    for (int i = 0; i < idx2samp.size(); ++i){
                        int samp = idx2samp[i];
                        if (ac2->first.first == -1){
                            double ef_val = 0.5;
                            if (expfracs.count(ac1->first) > 0 &&
                                expfracs[ac1->first].count(samp) > 0){
                                ef_val = expfracs[ac1->first][samp];
                            }
                            mixfrac_row.push_back(ef_val);
                        } else {
                            if (ef_all_avg && ac1->first.first == samp){
                                mixfrac_row.push_back(adjust_p_err(
                                    ac1->first.second / 2.0, e_r, e_a));
                            }
                            else if (ef_all_avg && ac2->first.first == samp){
                                mixfrac_row.push_back(adjust_p_err(
                                    ac2->first.second / 2.0, e_r, e_a));
                            }
                            else{
                                double ef1 = 0.5, ef2 = 0.5;
                                if (expfracs.count(ac1->first) > 0 &&
                                    expfracs[ac1->first].count(samp) > 0){
                                    ef1 = expfracs[ac1->first][samp];
                                }
                                if (expfracs.count(ac2->first) > 0 &&
                                    expfracs[ac2->first].count(samp) > 0){
                                    ef2 = expfracs[ac2->first][samp];
                                }
                                mixfrac_row.push_back(0.5 * ef1 + 0.5 * ef2);
                            }
                        }
                    }
                    if (inter_species){
                        mixfrac_row.push_back(adjust_p_err(0.0, e_r, e_a));
                    }
                    mixfracs.push_back(mixfrac_row);
                }
            }
            continue;
        }

        for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator ac1 = 
            indv_allelecounts[a->first].begin(); ac1 != indv_allelecounts[a->first].end(); 
            ++ac1){
            
            if ((!is_comb && ac1->first.first == a->second) || 
                (is_comb && ac1->first.first == comb.first)){
                for (map<pair<int, int>, pair<float, float> >::iterator ac2 = 
                    ac1->second.begin(); ac2 != ac1->second.end(); ++ac2){
                    
                    if ((!is_comb && ac2->first.first == -1) || 
                        (is_comb && ac2->first.first == comb.second)){
                        
                        // For heterotypic combo rows, current r-feedback should
                        // use the same parent-specific mixing ratio for the
                        // endogenous expected fraction p_e and for the ambient
                        // candidate mixfrac columns.  When r-feedback is off or
                        // no allele_ratio has been estimated yet, cell_r remains
                        // 0.5, preserving the historical fixed-50/50 behavior.
                        double cell_r = 0.5;
                        if (r_feedback_enabled && is_comb && comb.first != comb.second &&
                            allele_ratio.count(a->first) > 0){
                            cell_r = allele_ratio[a->first];
                            if (cell_r < 0.01) cell_r = 0.01;
                            if (cell_r > 0.99) cell_r = 0.99;
                        }

                        double expected;
                        if (!is_comb){
                            expected = adjust_p_err((double)ac1->first.second / 2.0, 
                                e_r, e_a);
                        }
                        else{
                            double p1 = (double)ac1->first.second / 2.0;
                            double p2 = (double)ac2->first.second / 2.0;
                            expected = adjust_p_err(cell_r * p1 + (1.0 - cell_r) * p2,
                                e_r, e_a);
                        }
                        
                        // Tetraploid-aware filtering
                        if (tetraploid_aware && is_comb){
                            if (amb_mu_available){
                                // Adaptive filter using estimated ambient profile
                                double this_p_c = 0.0;
                                if (amb_mu.count(ac1->first) > 0 && 
                                    amb_mu[ac1->first].count(ac2->first) > 0){
                                    this_p_c = amb_mu[ac1->first][ac2->first];
                                }
                                if (!category_passes_adaptive_filter(expected, this_p_c, min_signal_gap)){
                                    continue;
                                }
                            } else {
                                // Hard filter fallback
                                if (!category_passes_hard_filter(true, ac1->first.second, ac2->first.second)){
                                    continue;
                                }
                            }
                        }
                        
                        double ref = ac2->second.first;
                        double alt = ac2->second.second;
                        
                        n.push_back(ref+alt);
                        k.push_back(alt);
                        weights.push_back(weight);
                        if (!solve_for_c){
                            c.push_back(this_c);
                        }
                        p_e.push_back(expected);
                        vector<double> mixfrac_row;
                        for (int i = 0; i < idx2samp.size(); ++i){
                            int samp = idx2samp[i];
                            if (ac2->first.first == -1){
                                // Safe lookup: use 0.5 default if key missing from .condf
                                double ef_val = 0.5;
                                if (expfracs.count(ac1->first) > 0 && 
                                    expfracs[ac1->first].count(samp) > 0){
                                    ef_val = expfracs[ac1->first][samp];
                                }
                                mixfrac_row.push_back(ef_val);
                            }
                            else{
                                if (ef_all_avg && ac1->first.first == samp){
                                    mixfrac_row.push_back(adjust_p_err(
                                        ac1->first.second / 2.0, e_r, e_a));
                                }
                                else if (ef_all_avg && ac2->first.first == samp){
                                    mixfrac_row.push_back(adjust_p_err(
                                        ac2->first.second / 2.0, e_r, e_a));
                                }
                                else{
                                    // Safe lookup for doublet conditional fracs
                                    double ef1 = 0.5, ef2 = 0.5;
                                    if (expfracs.count(ac1->first) > 0 &&
                                        expfracs[ac1->first].count(samp) > 0){
                                        ef1 = expfracs[ac1->first][samp];
                                    }
                                    if (expfracs.count(ac2->first) > 0 &&
                                        expfracs[ac2->first].count(samp) > 0){
                                        ef2 = expfracs[ac2->first][samp];
                                    }
                                    // R-feedback: use the same cell_r used above for p_e.
                                    mixfrac_row.push_back(cell_r * ef1 + (1.0 - cell_r) * ef2);
                                }
                            }
                        }
                        if (inter_species){
                            // Reference alleles
                            mixfrac_row.push_back(adjust_p_err(0.0, e_r, e_a));
                        }
                        mixfracs.push_back(mixfrac_row);
                    }
                }
            }
        }
    }
}

/**
 * Bulk-mode data compiler for ambient profile estimation from empty droplets.
 *
 * Unlike compile_amb_prof_dat, this does not filter count rows by cell identity.
 * Every count row with ref+alt > 0 is included, because c=1.0 means the entire
 * signal is ambient. The mixture column for each row uses the condf-derived
 * expected alt fraction per individual (same as the singlet branch of
 * compile_amb_prof_dat), giving the solver full signal across all SNP categories.
 */
void contamFinder3::compile_bulk_amb_prof_dat(
    vector<vector<double> >& mixfracs,
    vector<double>& weights,
    vector<double>& n,
    vector<double>& k,
    vector<double>& p_e,
    vector<double>& c){

    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){

        for (auto ac1 = indv_allelecounts[a->first].begin();
            ac1 != indv_allelecounts[a->first].end(); ++ac1){

            for (auto ac2 = ac1->second.begin(); ac2 != ac1->second.end(); ++ac2){

                double ref = ac2->second.first;
                double alt = ac2->second.second;
                if (ref + alt <= 0) continue;

                // For bulk mode, p_e is not used in the mixture model because
                // c=1.0 eliminates the endogenous term entirely.
                // Set to 0.5 as a harmless default.
                double expected = 0.5;

                n.push_back(ref + alt);
                k.push_back(alt);
                weights.push_back(1.0);
                c.push_back(1.0);
                p_e.push_back(expected);

                // Build mixfrac row: condf-derived expected alt frac per individual.
                // Mirror the non-bulk logic for paired categories: when
                // ac2->first.first != -1, the count row represents a paired
                // genotype category and needs a 50/50 blend of both sides.
                bool is_paired = (ac2->first.first != -1);
                vector<double> mixfrac_row;
                for (int i = 0; i < (int)idx2samp.size(); ++i){
                    int samp = idx2samp[i];
                    double ef1 = 0.5;
                    if (expfracs.count(ac1->first) > 0 &&
                        expfracs[ac1->first].count(samp) > 0){
                        ef1 = expfracs[ac1->first][samp];
                    }

                    if (!is_paired){
                        mixfrac_row.push_back(ef1);
                    } else {
                        double ef2 = 0.5;
                        if (expfracs.count(ac2->first) > 0 &&
                            expfracs[ac2->first].count(samp) > 0){
                            ef2 = expfracs[ac2->first][samp];
                        }
                        mixfrac_row.push_back(0.5 * ef1 + 0.5 * ef2);
                    }
                }
                if (inter_species){
                    mixfrac_row.push_back(adjust_p_err(0.0, e_r, e_a));
                }
                mixfracs.push_back(mixfrac_row);
            }
        }
    }
}

/**
 * Models ambient RNA as a mixture of individuals. Updates the 
 * ambient RNA alt allele fractions so est_contam_cells() can find
 * contam levels in cells.
 *
 * Optionally can solve for global contam rate, and can choose whether to
 * use global contam rate (stored) or per-cell estimates in calculations.
 *
 * Returns log likelihood.
 * Edits init_c to updated value (if solve_for_c == true)
 */
double contamFinder3::update_amb_prof_mixture(bool solve_for_c, double& init_c, bool use_global_c){
    
    // Species mode bypass (section 4.5)
    if (this->species_mode){
        if (this->fixed_amb_prof){
            // Species-level pi was loaded; expand to individual-level via split rule
            expand_species_prior_to_indiv();
            rebuild_amb_mu_from_contam_prof();
            return 0.0;
        } else {
            // Estimate species-level pi via species-aggregated mixture components
            return solve_species_level_pi();
        }
    }

    // Fixed ambient profile bypass (Step 0a)
    if (this->fixed_amb_prof){
        rebuild_amb_mu_from_contam_prof();
        return 0.0;
    }

    vector<double> params;
    if (solve_for_c){
        params.push_back(init_c);
    }
    optimML::multivar_ml_solver solver(params, ll_amb_prof_mixture, dll_amb_prof_mixture);
    if (num_threads > 1){
        // Avoid multi-threading for evaluation for mixture proportion problems
        //solver.set_threads(num_threads);
        solver.set_bfgs_threads(num_threads);
    }

    vector<vector<double> > mixfracs;
    vector<double> weights;
    vector<double> n;
    vector<double> k;
    vector<double> p_e;
    vector<double> c;
    
    // Get data
    if (bulk_mode){
        compile_bulk_amb_prof_dat(mixfracs, weights, n, k, p_e, c);
    } else {
        compile_amb_prof_dat(solve_for_c, use_global_c, mixfracs, weights,
            n, k, p_e, c);
    }
    
    // Starting proportions: build in idx2samp order to match mixfracs columns
    // produced by compile_amb_prof_dat (which iterates idx2samp).
    vector<double> startprops;
    for (int i = 0; i < (int)idx2samp.size(); i++){
        int samp = idx2samp[i];
        if (contam_prof.count(samp) > 0){
            startprops.push_back(contam_prof[samp]);
        } else {
            startprops.push_back(1.0 / (double)idx2samp.size());
        }
    }
    if (inter_species && contam_prof.count(-1) > 0){
        startprops.push_back(contam_prof[-1]);
    } else if (inter_species){
        startprops.push_back(1.0 / (double)(idx2samp.size() + 1));
    }

    // Safety check: if sizes still don't match, fall back to uniform
    if (!mixfracs.empty() && startprops.size() != mixfracs[0].size()){
        fprintf(stderr, "WARNING: startprops size (%lu) != mixfracs columns (%lu); "
            "falling back to uniform\n",
            startprops.size(), mixfracs[0].size());
        startprops.clear();
        size_t n_cols = mixfracs[0].size();
        for (size_t i = 0; i < n_cols; i++){
            startprops.push_back(1.0 / (double)n_cols);
        }
    }
    
    // Set up ML solver 
    solver.add_mixcomp(mixfracs);
    solver.add_mixcomp_fracs(startprops);
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("p_e", p_e);
    solver.add_weights(weights);
    if (solve_for_c){
        solver.constrain_01(0);
    }
    else{
        solver.add_data("c", c);
    }
    //solver.set_delta(1e-6);

    if (solve_for_c){
        // First time. Try a few starting conditions and get the maximum LL.
        bool any_solve_succeeded = false;
        
        vector<double> lls;
        vector<vector<double> > mcs;
        vector<double> cs;
        int maxidx = 0;
        double maxll = -1e30;

        try {
            solver.solve();
            if (std::isfinite(solver.log_likelihood) && std::isfinite(solver.results[0])){
                any_solve_succeeded = true;
                maxll = solver.log_likelihood;
                lls.push_back(solver.log_likelihood);
                mcs.push_back(solver.results_mixcomp);
                cs.push_back(solver.results[0]);
            } else {
                fprintf(stderr, "WARNING: ambient profile mixture solver returned non-finite results on initial attempt\n");
            }
        } catch (...) {
            fprintf(stderr, "WARNING: ambient profile mixture solver failed on initial attempt\n");
        }
        
        try {
            vector<double> mptest;
            for (int i = 0; i < (int)startprops.size(); ++i){
                mptest.push_back(1.0/startprops.size());
            }
            solver.add_mixcomp_fracs(mptest);
            solver.solve();
            if (std::isfinite(solver.log_likelihood) && std::isfinite(solver.results[0])){
                any_solve_succeeded = true;
                if (solver.log_likelihood > maxll){
                    maxll = solver.log_likelihood;
                    maxidx = lls.size();
                }
                lls.push_back(solver.log_likelihood);
                mcs.push_back(solver.results_mixcomp);
                cs.push_back(solver.results[0]);
            }
        } catch (...) {
            // ignore, try other starting conditions
        }

        int ntrials = contam_prof.size() * n_mixprop_trials;
        for (int n = 0; n < ntrials; ++n){
            try {
                solver.randomize_mixcomps();
                solver.solve();
                if (std::isfinite(solver.log_likelihood) && std::isfinite(solver.results[0])){
                    any_solve_succeeded = true;
                    if (solver.log_likelihood > maxll){
                        maxll = solver.log_likelihood;
                        maxidx = lls.size();
                        lls.push_back(solver.log_likelihood);
                        mcs.push_back(solver.results_mixcomp);
                        cs.push_back(solver.results[0]);
                    }
                }
            } catch (...) {
                // ignore, try next starting condition
            }
        }

        if (any_solve_succeeded){
            solver.log_likelihood = maxll;
            solver.results[0] = cs[maxidx];
            solver.results_mixcomp = mcs[maxidx];
            // Update stored contamination profile (mixture of individuals)
            contam_prof.clear();
            for (int i = 0; i < (int)idx2samp.size(); ++i){
                int samp = idx2samp[i];
                contam_prof.insert(make_pair(samp, solver.results_mixcomp[i]));
            }
            if (inter_species){
                contam_prof.insert(make_pair(-1, 
                    solver.results_mixcomp[solver.results_mixcomp.size()-1]));
            }
        } else {
            fprintf(stderr, "WARNING: all ambient profile mixture solves failed; keeping initial profile\n");
        }
    }
    else{
        vector<double> trialprops;
        bool solve_ok = false;
        vector<double> maxres;
        double maxll = -1e30;
        try {
            solver.solve();
            if (std::isfinite(solver.log_likelihood)){
                maxres = solver.results_mixcomp;
                maxll = solver.log_likelihood;
                solve_ok = true;
            } else {
                fprintf(stderr, "WARNING: ambient profile mixture solver returned non-finite LL (no-c path); keeping initial profile\n");
            }
        } catch (...) {
            fprintf(stderr, "WARNING: ambient profile mixture solver failed (no-c path); keeping initial profile\n");
        }
        /* 
        for (int i = 0; i < contam_prof.size() * n_mixprop_trials; ++i){
            //rdirichlet(startprops, trialprops);
            //solver.add_mixcomp_fracs(trialprops);
            solver.randomize_mixcomps();
            solver.solve();
            if (maxll == 0 || solver.log_likelihood > maxll){
                maxres = solver.results_mixcomp;
                maxll = solver.log_likelihood;
            }
        }
        */
        // Update stored contamination profile (mixture of individuals)
        if (solve_ok){
            contam_prof.clear();
            for (int i = 0; i < (int)idx2samp.size(); ++i){
                int samp = idx2samp[i];
                contam_prof.insert(make_pair(samp, maxres[i]));
            }
            if (inter_species){
                contam_prof.insert(make_pair(-1, 
                    maxres[maxres.size()-1]));
            }
        }
    }
    
    // Update stored ambient RNA profile (allele matching fractions)
    for (int x = 0; x < idx2samp.size(); ++x){
        int i = idx2samp[x];
        for (int nalt = 0; nalt <= 2; ++nalt){
            pair<int, int> key = make_pair(i, nalt);
            if (amb_mu.count(key) == 0){
                map<pair<int, int>, double> m;
                amb_mu.insert(make_pair(key, m));
            }
            pair<int, int> nullkey = make_pair(-1,-1);
            if (amb_mu[key].count(nullkey) == 0){
                amb_mu[key].insert(make_pair(nullkey, 0));
            }
            double val = 0.0;
            for (map<int, double>::iterator cp = contam_prof.begin(); 
                cp != contam_prof.end(); ++cp){
                
                if (cp->first == -1){
                    val += cp->second * 0.0;
                }
                else if (ef_all_avg && cp->first == key.first){
                    val += cp->second * ((double)key.second/2.0);
                }
                else{
                    // Safe lookup: default to nalt/2.0 if key missing
                    double ef_val = (double)key.second / 2.0;
                    if (expfracs.count(key) > 0 &&
                        expfracs[key].count(cp->first) > 0){
                        ef_val = expfracs[key][cp->first];
                    }
                    val += cp->second * ef_val;
                }
            }
            amb_mu[key][nullkey] = val;
            for (int y = x + 1; y < idx2samp.size(); ++y){
                int j = idx2samp[y];
                for (int nalt2 = 0; nalt2 <= 2; ++nalt2){
                    pair<int, int> key2 = make_pair(j, nalt2);
                    if (amb_mu[key].count(key2) == 0){
                        amb_mu[key].insert(make_pair(key2, 0));
                    }
                    double val = 0.0;
                    for (map<int, double>::iterator cp = contam_prof.begin(); 
                        cp != contam_prof.end(); ++cp){
                        
                        if (cp->first == -1){
                            val += cp->second * adjust_p_err(0.0, e_r, e_a);
                        }
                        else{
                            if (ef_all_avg && cp->first == key.first){
                                val += cp->second * adjust_p_err(
                                    (double)key.second / 2.0, e_r, e_a);
                            }
                            else if (ef_all_avg && cp->first == key2.first){
                                val += cp->second * adjust_p_err(
                                    (double)key2.second / 2.0, e_r, e_a);
                            }
                            else{
                                // Safe lookup for doublet conditional fracs
                                double ef1 = (double)key.second / 2.0;
                                double ef2 = (double)key2.second / 2.0;
                                if (expfracs.count(key) > 0 &&
                                    expfracs[key].count(cp->first) > 0){
                                    ef1 = expfracs[key][cp->first];
                                }
                                if (expfracs.count(key2) > 0 &&
                                    expfracs[key2].count(cp->first) > 0){
                                    ef2 = expfracs[key2][cp->first];
                                }
                                val += cp->second * (0.5*ef1 + 0.5*ef2);
                            }
                        }
                    }
                    amb_mu[key][key2] = val;
                }
            }
        }
    }
    if (solve_for_c){
        init_c = solver.results[0];
    }
    return solver.log_likelihood;
}

/**
 * Gets variances on contamination profile proportions by bootstrapping.
 * Fits a Dirichlet distribution to all bootstrap samples; Dirichlet
 * concentration parameters will then be reported in output files.
 */
void contamFinder3::bootstrap_amb_prof(int n_boots, map<int, double>& dirichlet_params){
    // Assumes we have already solved everything and that this is at the end.
    if (n_boots <= 0){
        fprintf(stderr, "Bootstrap disabled (n_boots <= 0); skipping concentration fit\n");
        return;
    }

    // Compile individual-level data (used by both species and non-species paths)
    vector<vector<double> > mixfracs;
    vector<double> weights;
    vector<double> n;
    vector<double> k;
    vector<double> p_e;
    vector<double> c;
    
    if (bulk_mode){
        compile_bulk_amb_prof_dat(mixfracs, weights, n, k, p_e, c);
    } else {
        compile_amb_prof_dat(false, false, mixfracs, weights, n, k, p_e, c);
    }

    if (mixfracs.empty() || n.empty()){
        fprintf(stderr, "WARNING: no data for bootstrap; skipping\n");
        return;
    }

    // ---- Species-mode bootstrap: aggregate to species level ----
    if (species_mode){
        // Keep bootstrap alpha expansion consistent with solve_species_level_pi()
        // and expand_species_prior_to_indiv(). Skip in bulk mode because the
        // synthetic placeholder assignment would create meaningless loading weights.
        if (!bulk_mode && indiv_loading_weights.empty()){
            compute_loading_weights();
        }

        int n_obs = (int)mixfracs.size();
        int n_sp = (int)panel_meta.species_list.size();

        // Build sample position lookup
        map<int, int> samp_to_pos;
        for (int i = 0; i < (int)idx2samp.size(); i++){
            samp_to_pos[idx2samp[i]] = i;
        }

        // Compute weighted active count per species
        set<int> active_set(idx2samp.begin(), idx2samp.end());
        vector<double> sp_weight_sum(n_sp, 0.0);
        for (int s = 0; s < n_sp; s++){
            const string& sp = panel_meta.species_list[s];
            if (panel_meta.species_to_sample_indices.count(sp) > 0){
                for (int idx : panel_meta.species_to_sample_indices.at(sp)){
                    if (active_set.count(idx) > 0){
                        sp_weight_sum[s] += panel_meta.get_weight(sp, idx);
                    }
                }
            }
        }

        // Identify active species (exclude zero-weight species from solver
        // to prevent degenerate mass assignment to absent species)
        vector<int> active_sp_idx;
        map<int, int> full_to_active;
        for (int s = 0; s < n_sp; s++){
            if (sp_weight_sum[s] > 0.0){
                full_to_active[s] = (int)active_sp_idx.size();
                active_sp_idx.push_back(s);
            }
        }
        int n_active_sp = (int)active_sp_idx.size();

        if (n_active_sp == 0){
            fprintf(stderr, "WARNING: no active species for bootstrap; skipping\n");
            return;
        }

        // Remap mixfracs to active-species level using weighted average
        vector<vector<double> > sp_mixfracs(n_obs, vector<double>(n_active_sp, 0.0));
        for (int obs = 0; obs < n_obs; obs++){
            for (int a = 0; a < n_active_sp; a++){
                int s = active_sp_idx[a];
                const string& sp = panel_meta.species_list[s];
                double wsum = 0.0;
                if (panel_meta.species_to_sample_indices.count(sp) > 0){
                    for (int idx : panel_meta.species_to_sample_indices.at(sp)){
                        if (samp_to_pos.count(idx) > 0){
                            double w = panel_meta.get_weight(sp, idx);
                            wsum += w * mixfracs[obs][samp_to_pos[idx]];
                        }
                    }
                }
                sp_mixfracs[obs][a] = wsum / sp_weight_sum[active_sp_idx[a]];
            }
        }

        // Active-species startprops and MLE fracs from species_contam_prof
        vector<double> startprops(n_active_sp, 1.0 / (double)n_active_sp);
        vector<double> mle_fracs(n_active_sp, 0.0);
        vector<vector<double> > dirprops(n_active_sp);
        for (int a = 0; a < n_active_sp; a++){
            const string& sp = panel_meta.species_list[active_sp_idx[a]];
            if (species_contam_prof.count(sp) > 0){
                startprops[a] = species_contam_prof[sp];
                mle_fracs[a] = species_contam_prof[sp];
            }
        }
        // Normalize startprops to sum to 1 (inactive species mass excluded)
        {
            double sp_sum = 0.0;
            for (int a = 0; a < n_active_sp; a++) sp_sum += startprops[a];
            if (sp_sum > 0.0){
                for (int a = 0; a < n_active_sp; a++) startprops[a] /= sp_sum;
            }
        }

        // Bootstrap at species level. Parallelize over independent replicates.
        // Each worker stores one result vector; dirprops is populated serially
        // after the OpenMP region to avoid concurrent vector push_back.
        vector<vector<double> > boot_results(n_boots);
        int nt_boot = (num_threads > 1) ? std::min(num_threads, n_boots) : 1;

        #pragma omp parallel for num_threads(nt_boot) schedule(dynamic, 1)
        for (int b = 0; b < n_boots; ++b){
            mt19937 rand_gen(1337 + b);
            uniform_int_distribution<int> uni_dist(0, n_obs - 1);

            vector<vector<double> > mf_boot;
            vector<double> w_boot, n_boot, k_boot, pe_boot, c_boot;
            mf_boot.reserve(n_obs);
            w_boot.reserve(n_obs);
            n_boot.reserve(n_obs);
            k_boot.reserve(n_obs);
            pe_boot.reserve(n_obs);
            c_boot.reserve(n_obs);

            for (int x = 0; x < n_obs; ++x){
                int r = uni_dist(rand_gen);
                mf_boot.push_back(sp_mixfracs[r]);
                w_boot.push_back(weights[r]);
                n_boot.push_back(n[r]);
                k_boot.push_back(k[r]);
                pe_boot.push_back(p_e[r]);
                c_boot.push_back(c[r]);
            }

            vector<double> params;
            optimML::multivar_ml_solver solver(params, ll_amb_prof_mixture,
                dll_amb_prof_mixture);
            solver.add_mixcomp(mf_boot);
            solver.add_mixcomp_fracs(startprops);
            solver.add_data("n", n_boot);
            solver.add_data("k", k_boot);
            solver.add_data("p_e", pe_boot);
            solver.add_data("c", c_boot);

            try {
                solver.solve();
                if (std::isfinite(solver.log_likelihood)){
                    for (int x = 0; x < (int)solver.results_mixcomp.size() &&
                         x < n_active_sp; ++x){
                        boot_results[b].push_back(solver.results_mixcomp[x]);
                    }
                }
            } catch (...) {
                // Skip failed bootstrap sample
            }
        }

        for (int b = 0; b < n_boots; ++b){
            if ((int)boot_results[b].size() != n_active_sp) continue;
            for (int x = 0; x < n_active_sp; ++x){
                dirprops[x].push_back(boot_results[b][x]);
            }
        }
        fprintf(stderr, "Bootstrap samples complete: %d requested\n", n_boots);

        // Fit Dirichlet at active-species level
        vector<double> dirichlet_soln;
        fit_dirichlet(mle_fracs, dirprops, dirichlet_soln);

        // Expand active-species Dirichlet params to individual level using the
        // same effective split as expand_species_prior_to_indiv().
        //
        // Semantics:
        //   - base panel weight comes from panel_meta.get_weight(sp, idx)
        //   - if indiv_loading_weights is non-empty, missing idx means zero
        //   - if an entire species has zero loading-weighted mass, fall back to
        //     base panel weights for that species only
        dirichlet_params.clear();

        for (int a = 0; a < n_active_sp && a < (int)dirichlet_soln.size(); a++){
            int s = active_sp_idx[a];
            const string& sp = panel_meta.species_list[s];

            if (panel_meta.species_to_sample_indices.count(sp) == 0) continue;

            const auto& sp_indices = panel_meta.species_to_sample_indices.at(sp);

            double eff_w_sum = 0.0;
            bool use_loading = !indiv_loading_weights.empty();

            if (use_loading){
                for (int idx : sp_indices){
                    if (active_set.count(idx) > 0){
                        auto it = indiv_loading_weights.find(idx);
                        if (it != indiv_loading_weights.end()){
                            eff_w_sum += panel_meta.get_weight(sp, idx) * it->second;
                        }
                    }
                }
                if (eff_w_sum <= 0.0){
                    use_loading = false;
                }
            }

            if (!use_loading){
                eff_w_sum = 0.0;
                for (int idx : sp_indices){
                    if (active_set.count(idx) > 0){
                        eff_w_sum += panel_meta.get_weight(sp, idx);
                    }
                }
            }

            if (eff_w_sum <= 0.0) continue;

            for (int idx : sp_indices){
                if (active_set.count(idx) > 0){
                    double w = 0.0;
                    if (use_loading){
                        auto it = indiv_loading_weights.find(idx);
                        w = (it != indiv_loading_weights.end())
                          ? panel_meta.get_weight(sp, idx) * it->second
                          : 0.0;
                    } else {
                        w = panel_meta.get_weight(sp, idx);
                    }
                    if (w > 0.0){
                        dirichlet_params[idx] += (w / eff_w_sum) * dirichlet_soln[a];
                    }
                }
            }
        }
        // Store species-level Dirichlet concentrations for output.
        // Active species get fitted values; inactive species get 0.0.
        species_contam_prof_conc.clear();
        for (int s = 0; s < n_sp; s++){
            if (full_to_active.count(s) > 0 &&
                full_to_active[s] < (int)dirichlet_soln.size()){
                species_contam_prof_conc[panel_meta.species_list[s]] =
                    dirichlet_soln[full_to_active[s]];
            } else {
                species_contam_prof_conc[panel_meta.species_list[s]] = 0.0;
            }
        }

        return;
    }

    // ---- Individual-level bootstrap (non-species mode) ----

    // Build startprops/mle_fracs/dirprops in idx2samp order to match
    // mixfracs columns from compile_amb_prof_dat.
    vector<vector<double> > dirprops;
    vector<double> mle_fracs;
    vector<double> startprops;
    for (int i = 0; i < (int)idx2samp.size(); i++){
        int samp = idx2samp[i];
        double val = 1.0 / (double)idx2samp.size();
        if (contam_prof.count(samp) > 0){
            val = contam_prof[samp];
        }
        startprops.push_back(val);
        dirprops.push_back(vector<double>());
        mle_fracs.push_back(val);
    }
    if (inter_species && contam_prof.count(-1) > 0){
        startprops.push_back(contam_prof[-1]);
        dirprops.push_back(vector<double>());
        mle_fracs.push_back(contam_prof[-1]);
    }

    // Safety check: if sizes still don't match, fall back to uniform
    if (!mixfracs.empty() && startprops.size() != mixfracs[0].size()){
        fprintf(stderr, "WARNING: bootstrap startprops size (%lu) != mixfracs columns (%lu); "
            "falling back to uniform\n", startprops.size(), mixfracs[0].size());
        size_t n_cols = mixfracs[0].size();
        startprops.clear();
        dirprops.clear();
        mle_fracs.clear();
        for (size_t i = 0; i < n_cols; i++){
            startprops.push_back(1.0 / (double)n_cols);
            dirprops.push_back(vector<double>());
            mle_fracs.push_back(1.0 / (double)n_cols);
        }
    }
    
    // Parallel bootstrap over independent replicates.  Each replicate gets a
    // deterministic seed and writes to boot_results[b]; dirprops is filled
    // serially afterward.
    int n_obs = (int)n.size();
    int n_props = (int)startprops.size();
    vector<vector<double> > boot_results(n_boots);
    int nt_boot = (num_threads > 1) ? std::min(num_threads, n_boots) : 1;

    #pragma omp parallel for num_threads(nt_boot) schedule(dynamic, 1)
    for (int b = 0; b < n_boots; ++b){
        mt19937 rand_gen(7331 + b);
        uniform_int_distribution<int> uni_dist(0, n_obs - 1);

        vector<vector<double> > mixfracs_boot;
        vector<double> weights_boot;
        vector<double> n_boot;
        vector<double> k_boot;
        vector<double> p_e_boot;
        vector<double> c_boot;
        mixfracs_boot.reserve(n_obs);
        weights_boot.reserve(n_obs);
        n_boot.reserve(n_obs);
        k_boot.reserve(n_obs);
        p_e_boot.reserve(n_obs);
        c_boot.reserve(n_obs);

        for (int x = 0; x < n_obs; ++x){
            int r = uni_dist(rand_gen);
            mixfracs_boot.push_back(mixfracs[r]);
            weights_boot.push_back(weights[r]);
            n_boot.push_back(n[r]);
            k_boot.push_back(k[r]);
            p_e_boot.push_back(p_e[r]);
            c_boot.push_back(c[r]);
        }

        vector<double> params;
        optimML::multivar_ml_solver solver(params, ll_amb_prof_mixture,
            dll_amb_prof_mixture);
        solver.add_mixcomp(mixfracs_boot);
        solver.add_mixcomp_fracs(startprops);
        solver.add_data("n", n_boot);
        solver.add_data("k", k_boot);
        solver.add_data("p_e", p_e_boot);
        solver.add_data("c", c_boot);

        try {
            solver.solve();
            if (std::isfinite(solver.log_likelihood)){
                for (int x = 0; x < (int)solver.results_mixcomp.size() && x < n_props; ++x){
                    boot_results[b].push_back(solver.results_mixcomp[x]);
                }
            }
        } catch (...) {
            // Skip failed bootstrap sample
        }
    }

    for (int b = 0; b < n_boots; ++b){
        if ((int)boot_results[b].size() != n_props) continue;
        for (int x = 0; x < n_props; ++x){
            dirprops[x].push_back(boot_results[b][x]);
        }
    }
    fprintf(stderr, "Bootstrap samples complete: %d requested\n", n_boots);
    
    vector<double> dirichlet_soln;
    fit_dirichlet(mle_fracs, dirprops, dirichlet_soln);
    
    dirichlet_params.clear();
    for (int i = 0; i < (int)dirichlet_soln.size() && i < (int)idx2samp.size(); ++i){
        int samp = idx2samp[i];
        dirichlet_params.insert(make_pair(samp, dirichlet_soln[i]));
    }
    if (inter_species){
        dirichlet_params.insert(make_pair(-1, 
           dirichlet_soln[dirichlet_soln.size()-1]));
    }
}

/**
 * Returns whether anything was changed.
 */
bool contamFinder3::reclassify_cells(){
    // Assume uniform global contam rate
    double c = contam_cell_prior;
    bool changed = false;
    
    int n_reassigned = 0;
    
    bool reweight_doublets = (doublet_rate > 0 && doublet_rate < 1);
    
    pair<double, double> betaparams;
    if (contam_cell_prior > 0 && contam_cell_prior_var > 0){
        betaparams = beta_moments(contam_cell_prior, contam_cell_prior_var);
    }

    map<int, double> priorweights;
    map<int, double>* priorweights_ptr = NULL;
    if (reweight_doublets){
        
        // Notify the populate_llr_table function that we have LLR elements
        // to add into the table
        priorweights_ptr = &priorweights;
        
        // Compute total fraction of RNA per identity (as if bulk)
        // Store total # cells mapped to each identity, counting halves of doublets
        // once and singlets twice
        map<int, int> acounts;

        // Compute total number of cells with each identity
        // Store total # cells mapped to each identity, including doublets
        map<int, int> acounts_incld;
        int atot = 0;
        int atot_incld = 0;
        
        for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); 
            a != assn.end(); ++a){
            if (a->second >= n_samples){
                pair<int, int> comb = idx_to_hap_comb(a->second, n_samples);
                if (acounts.count(comb.first) == 0){
                    acounts.insert(make_pair(comb.first, 0));
                }
                if (acounts.count(comb.second) == 0){
                    acounts.insert(make_pair(comb.second, 0));
                }
                acounts[comb.first]++;
                acounts[comb.second]++;
                atot += 2;
                if (acounts_incld.count(a->second) == 0){
                    acounts_incld.insert(make_pair(a->second, 0));
                }
                acounts_incld[a->second]++;
                atot_incld++;
            }

            else{
                if (acounts.count(a->second) == 0){
                    acounts.insert(make_pair(a->second, 0));
                }
                acounts[a->second] += 2;
                atot += 2;   
                if (acounts_incld.count(a->second) == 0){
                   acounts_incld.insert(make_pair(a->second, 0));
                }
                acounts_incld[a->second]++;
                atot_incld++; 
            }
        }
        vector<int> samps;
        for (int i = 0; i < idx2samp.size(); ++i){
            int samp = idx2samp[i];
            samps.push_back(samp);
        }
        sort(samps.begin(), samps.end());
        
        for (int i = 0; i < samps.size(); ++i){
            int si = samps[i];
            
            // Compute expected fraction of this singlet identity
            // if d is doublet rate and x is the bulk fraction of this
            // identity, then p = (1-d)*s + d*s*s
            double si_p = (double)acounts[si]/(double)atot;
            double prob = (1.0 - doublet_rate)*si_p + doublet_rate*si_p*si_p;
            
            // Add in log likelihood of present number of counts under the model
            double ll = dbinom(atot_incld, acounts_incld[si], prob);
            priorweights.insert(make_pair(si, ll));
            
            for (int j = i + 1; j < samps.size(); ++j){
                int sj = samps[j];
                int sk = hap_comb_to_idx(si, sj, n_samples);

                // Compute expected fraction of this doublet identity
                // if d is doublet rate, x is bulk fraction of ID1 and y is
                // bulk fraction of ID2, then p = 2*d*x*y
                double sj_p = (double)acounts[sj]/(double)atot;
                double prob = 2.0*doublet_rate*si_p*sj_p;
                
                double ll = dbinom(atot_incld, acounts_incld[sk], prob);
                priorweights.insert(make_pair(sk, ll));
            }
            
        }
    }
    else if (false){
        priorweights_ptr = &priorweights;
        pair<double, double> betaparams = beta_moments(contam_cell_prior, contam_cell_prior_var);
        map<int, double> crmean;
        map<int, int> crtot;
        map<int, vector<double> > crmed;
        for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
            a != assn.end(); ++a){
            if (contam_rate.count(a->first) > 0){
                if (crmean.count(a->second) > 0){
                    crmean.insert(make_pair(a->second, 0.0));
                    crtot.insert(make_pair(a->second, 0));
                    vector<double> v;
                    crmed.insert(make_pair(a->second, v));
                }
                crmean[a->second] += contam_rate[a->first];
                crtot[a->second]++;
                crmed[a->second].push_back(contam_rate[a->first]);
            }
        }
        for (map<int, double>::iterator crm = crmean.begin(); crm != crmean.end(); ++crm){
            crm->second /= (double)crtot[crm->first];
            /*
            double p = 1.0 - pbeta(crm->second, betaparams.first, betaparams.second);
            if (p == 0){
                p += 1e-6;
            }
            else if (p == 1){
                p -= 1e-6;
            }
            p = log2(p);
            */
            /*
            double med;
            if (crmed[crm->first].size() % 2 == 0){
                double m1 = crmed[crm->first][crmed[crm->first].size()/2];
                double m2 = crmed[crm->first][crmed[crm->first].size()/2-1];
                med = (m1+m2)/2.0;
            }
            else{
                med = crmed[crm->first][(crmed[crm->first].size()-1)/2];
            }
            */
            double p = dbeta(crm->second, betaparams.first, betaparams.second);
            priorweights.insert(make_pair(crm->first, p));
        }
    }
    
    // If the user has set doublet rate to 0 or 1, then instead of doing
    // complicated stuff above to try to enforce correct proportions of 
    // all individuals, we'll just pass that value to populate_llr_table,
    // which will exclude doublets or singlets accordingly.
    
    // Otherwise, no need to include doublet rate since we will already be
    // including it with the "prior_weights" parameter.
    double doub_rate_table = 0.5;
    if (doublet_rate == 0 || doublet_rate == 1){
        doub_rate_table = doublet_rate;
    }

    set<unsigned long> cell_rm;
   
    /*
    vector<pair<double, int> > cpsort;
    for (map<int, double>::iterator cp = contam_prof.begin(); cp != contam_prof.end(); ++cp){
        if (cp->first != -1){
            cpsort.push_back(make_pair(-cp->second, cp->first));
        }
    }
    sort(cpsort.begin(), cpsort.end());
    
    double newcsum = 0.0;
    double newccount = 0.0;

    int nproc = 0;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){
        
        vector<int> possibilities;
        vector<double> possibility_lls;
        vector<double> possibility_cs;
    
        string bcstr = bc2str(a->first);
        fprintf(stderr, "%s\n", bcstr.c_str());

        if (a->second >= n_samples){
            pair<int, int> combo = idx_to_hap_comb(a->second, n_samples);
            fprintf(stderr, "orig (%d %d) %f %f\n", combo.first, combo.second, contam_rate[a->first],
                contam_rate_ll[a->first]);
            if (allowed_ids2.size() == 0 || allowed_ids2.find(combo.first) != allowed_ids2.end()){
                possibilities.push_back(combo.first);
            }
            if (allowed_ids2.size() == 0 || allowed_ids2.find(combo.second) != allowed_ids2.end()){
                possibilities.push_back(combo.second);
            }
            for (int x = 0; x < cpsort.size(); ++x){
                int j = cpsort[x].second;
                if (j != combo.first && j != combo.second){
                    int k1;
                    int k2;
                    if (j < combo.first){
                        k1 = hap_comb_to_idx(j, combo.first, n_samples);
                    }
                    else{
                        k1 = hap_comb_to_idx(combo.first, j, n_samples);
                    }
                    if (j < combo.second){
                        k2 = hap_comb_to_idx(j, combo.second, n_samples);
                    }
                    else{
                        k2 = hap_comb_to_idx(combo.second, j, n_samples);
                    }
                    if (allowed_ids2.size() == 0 || allowed_ids2.find(k1) != allowed_ids2.end()){
                        possibilities.push_back(k1);
                    }
                    if (allowed_ids2.size() == 0 || allowed_ids2.find(k2) != allowed_ids2.end()){
                        possibilities.push_back(k2);
                    }
                }
            }
        }
        else{
            fprintf(stderr, "orig (%d) %f %f\n", a->second, contam_rate[a->first], contam_rate_ll[a->first]);
            //for (int x = 0; x < idx2samp.size(); ++x){
            for (int x = 0; x < cpsort.size(); ++x){
                int j = cpsort[x].second;    
            //int j = idx2samp[x];
                if (j != a->second){
                    int k;
                    if (j < a->second){
                        k = hap_comb_to_idx(j, a->second, n_samples);
                    }
                    else{
                        k = hap_comb_to_idx(a->second, k, n_samples);
                    }
                    if (allowed_ids2.size() == 0 || allowed_ids2.find(k) != allowed_ids2.end()){
                        possibilities.push_back(k);
                    }
                }
            }
        }
        int maxidx = -1;
        double maxll = contam_rate_ll[a->first];
        double maxc = contam_rate[a->first];
        double maxcse = contam_rate_se[a->first];
        for (int x = 0; x < possibilities.size(); ++x){
            int j = possibilities[x];
            bool is_comb = false;
            pair<int, int> comb;
            if (j >= n_samples){
                is_comb = true;
                comb = idx_to_hap_comb(j, n_samples);
            }
            vector<double> n;
            vector<double> k;
            vector<double> p_c;
            vector<double> p_e;
            for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator ac = 
                indv_allelecounts[a->first].begin(); ac != indv_allelecounts[a->first].end();
                ++ac){
                if ((!is_comb && ac->first.first == j) || (is_comb && ac->first.first == comb.first)){
                    for (map<pair<int, int>, pair<float, float> >::iterator ac2 = 
                        ac->second.begin(); ac2 != ac->second.end(); ++ac2){
                        bool breakout = false;
                        double p_e_this;
                        bool keep = false;
                        if ((!is_comb && ac2->first.first == -1)){
                            p_e_this = (double)ac->first.second/2.0;
                            keep = true;
                            breakout = true;
                        }
                        else if ((is_comb && ac2->first.first == comb.second)){
                            p_e_this = (double)(ac->first.second + ac2->first.second)/4.0;
                            keep = true;
                        }
                        if (keep){
                            double p_c_this = amb_mu[ac->first][ac2->first];
                            n.push_back(ac2->second.first + ac2->second.second);
                            k.push_back(ac2->second.second);
                            p_c.push_back(p_c_this);
                            p_e.push_back(p_e_this);
                            
                        }
                        if (breakout){
                            break;
                        }
                    }
                }
            }
            if (n.size() == 0){
                continue;
            }
            optimML::brent_solver c_cell(ll_c, dll_dc, d2ll_dc2);
            if (num_threads > 1){
                c_cell.set_threads(num_threads);
            }
            c_cell.constrain_01();
            c_cell.add_data("n", n);
            c_cell.add_data("k", k);
            c_cell.add_data("p_c", p_c);
            c_cell.add_data("p_e", p_e);
            if (contam_cell_prior > 0 && contam_cell_prior_var > 0){
                c_cell.add_beta_prior(betaparams.first, betaparams.second);
            }
            c_cell.set_maxiter(-1);
            double c_cell_map;
            bool solve_ok = false;
            try {
                c_cell_map = c_cell.solve(contam_rate[a->first],1);
                solve_ok = true;
            } catch (...) {
                // Solver failed for this alternative identity - skip it
                continue;
            }
            //double c_cell_map = c_cell.solve(0,1);
            if (is_comb){
                fprintf(stderr, "  (%d %d) %f %f\n", comb.first, comb.second, c_cell_map,
                    c_cell.log_likelihood);
            }
            else{
                fprintf(stderr, "  (%d) %f %f\n", j, c_cell_map, c_cell.log_likelihood);
            }
            if (c_cell.root_found){
                if (maxll == 0 || c_cell.log_likelihood > maxll){
                    maxll = c_cell.log_likelihood;
                    maxidx = j;
                    maxc = c_cell_map;
                    maxcse = c_cell.se;
                }
                else if (a->second < n_samples || x > 1){
                    //break;
                }    
            }
        }
        if (maxll != 0 && maxidx != -1){
            ++n_reassigned;
            fprintf(stderr, "reassign %d / %d\n", n_reassigned, nproc);
            // Update
            contam_rate_ll[a->first] = maxll;
            contam_rate[a->first] = maxc;
            contam_rate_se[a->first] = maxcse;
            assn[a->first] = maxidx;
            
            newcsum += maxc;
            newccount++;
            
            fprintf(stderr, " c = %f\n", newcsum / newccount);
        }
        nproc++;
    }
    fprintf(stderr, " %d cells reassigned\n", n_reassigned);
    return n_reassigned > 0;
    */ 

    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); 
        a != assn.end(); ++a){
        
        // Tetraploid-aware: skip locked identities
        if (tetraploid_aware && locked_identities.count(a->second) > 0){
            continue;
        }
        
        // Tetraploid-aware: for safe singlets, restrict reassignment candidates
        // to other safe singlets only. Build a restricted allowed_ids2 for this cell.
        set<int> cell_allowed_ids2;
        set<int>* cell_allowed_ids2_ptr = &allowed_ids2;
        if (tetraploid_aware && safe_singlets.count(a->second) > 0){
            // Only allow reassignment to other safe singlets
            for (set<int>::iterator ss = safe_singlets.begin(); ss != safe_singlets.end(); ++ss){
                cell_allowed_ids2.insert(*ss);
            }
            cell_allowed_ids2_ptr = &cell_allowed_ids2;
        }
        
        // Get a table of log likelihood ratios between every possible
        // pair of identities
        map<int, map<int, double> > llrs;
        llr_table tab(n_samples);
        
        //c = contam_rate[a->first]; 
        bool success = populate_llr_table(indv_allelecounts[a->first], llrs, tab, n_samples, 
            allowed_ids, *cell_allowed_ids2_ptr, doub_rate_table, e_r, e_a, priorweights_ptr,
            true, contam_cell_prior, 0, &amb_mu);

        if (success){
            
            // If we reassigned an identity and it's different than the old one,
            // also infer a new contamination rate. Accept the change if the
            // overall log likelihood of the new identity + new contamination rate
            // beats the log likelihood of the old identity + old contamination rate.

            int a_new;
            double llr_new; 
            tab.get_max(a_new, llr_new);
            if (a_new != -1 && llr_new > 0){
                if (a_new != a->second){
                    
                    // Get new contam rate conditional on candidate new identity
                    vector<double> n;
                    vector<double> k;
                    vector<double> p_e;
                    vector<double> p_c;

                    if (a_new >= n_samples){
                        pair<int, int> comb = idx_to_hap_comb(a_new, n_samples);
                        for (int x = 0; x <= 2; ++x){
                            pair<int, int> k1 = make_pair(comb.first, x);
                            for (int y = 0; y <= 2; ++y){
                                pair<int, int> k2 = make_pair(comb.second, y);
                                if (true){
                                //if (x != y){
                                    double ref = indv_allelecounts[a->first][k1][k2].first;
                                    double alt = indv_allelecounts[a->first][k1][k2].second;
                                    n.push_back(ref+alt);
                                    k.push_back(alt);
                                    p_e.push_back(adjust_p_err((double)(x+y)/4.0, e_r, e_a));
                                    p_c.push_back(amb_mu[k1][k2]);
                                }
                            }
                        }
                    }
                    else{
                        pair<int, int> nullkey = make_pair(-1, -1);
                        for (int x = 0; x <= 2; ++x){
                            pair<int, int> key = make_pair(a_new, x);
                            double ref = indv_allelecounts[a->first][key][nullkey].first;
                            double alt = indv_allelecounts[a->first][key][nullkey].second;
                            n.push_back(ref+alt);
                            k.push_back(alt);
                            p_e.push_back(adjust_p_err((double)x/2.0, e_r, e_a));
                            p_c.push_back(amb_mu[key][nullkey]);
                        }
                    }
                    optimML::brent_solver c_cell(ll_c, dll_dc, d2ll_dc2);
                    
                    if (num_threads > 1){
                        c_cell.set_threads(num_threads);
                    }
                    c_cell.add_data("n", n);
                    c_cell.add_data("k", k);
                    c_cell.add_data("p_e", p_e);
                    c_cell.add_data("p_c", p_c);
                    c_cell.constrain_01();
                    if (contam_cell_prior > 0 && contam_cell_prior_var > 0){
                        c_cell.add_beta_prior(betaparams.first, betaparams.second);
                    }
                    c_cell.set_maxiter(-1);
                    
                    double c_cell_map = 1.0;
                    bool rc_solve_ok = false;
                    try {
                        c_cell_map = c_cell.solve(0,1);
                        rc_solve_ok = true;
                    } catch (...) {
                        // Solver failed for candidate reassignment - skip
                    }
                    
                    if (!rc_solve_ok) {
                        continue;
                    }
                    
                    // Only accept the change if the new assignment + contam rate inference has a higher
                    // log likelihood than the older assignment + contam rate inference
                    //if (c_cell.root_found && c_cell_map < contam_rate[a->first] && 
                    //    (contam_rate_ll[a->first] == 0 || 
                    //    c_cell.log_likelihood > contam_rate_ll[a->first])){
                    
                    //if (c_cell.root_found &&
                    //    (contam_rate_ll[a->first] == 0 || c_cell.log_likelihood < contam_rate_ll[a->first])){
  
                    if (c_cell.root_found &&
//                        c_cell_map >= contam_rate[a->first] &&  
                        (contam_rate_ll[a->first] == 0 || c_cell.log_likelihood >= contam_rate_ll[a->first])){
                        n_reassigned++;
                        
                        //fprintf(stdout, "%s\t%f\t%f\t%f\t%f\n", bc2str(a->first).c_str(),
                        //    contam_rate[a->first], c_cell_map, contam_rate_ll[a->first],
                        //    c_cell.log_likelihood);


                        changed = true;
                        a->second = a_new;
                        assn_llr[a->first] = llr_new;
                        contam_rate[a->first] = c_cell_map;
                        contam_rate_se[a->first] = c_cell.se;
                        contam_rate_ll[a->first] = c_cell.log_likelihood;
                        // Clean up allele_ratio if new identity is not heterotypic
                        if (a_new < n_samples){
                            // Reassigned to singlet: r is meaningless
                            allele_ratio.erase(a->first);
                            allele_ratio_se.erase(a->first);
                        } else {
                            pair<int, int> new_combo = idx_to_hap_comb(a_new, n_samples);
                            if (new_combo.first == new_combo.second){
                                // Reassigned to homotypic combo: r is meaningless
                                allele_ratio.erase(a->first);
                                allele_ratio_se.erase(a->first);
                            }
                        }
                    }
                }
            }
            else{
                cell_rm.insert(a->first);        
            }
        }
        else{
            // Keep original? Delete original?
            cell_rm.insert(a->first);
        }
    }
    
    fprintf(stderr, "  %d cells reassigned\n", n_reassigned);
    fprintf(stderr, "  %ld cells removed\n", cell_rm.size());

    if (cell_rm.size() > 0){
        changed = true;
    }
    for (set<unsigned long>::iterator rm = cell_rm.begin(); rm != cell_rm.end(); ++rm){
        assn.erase(*rm);
        assn_llr.erase(*rm);
        contam_rate.erase(*rm);
        contam_rate_se.erase(*rm);
        allele_ratio.erase(*rm);
        allele_ratio_se.erase(*rm);
    }
    if (species_mode){
        compute_loading_weights();
    }
    return changed;
}

/**
 * Get best estimates of (residual) reference & alt allele misreading
 * error rates. This should be done after estimating all other parameters,
 * to see how close to zero these values are. High values indicate 
 * discordance between the data and variant calls not reflected in the
 * contamination model.
 */
pair<double, double> contamFinder3::est_error_rates(bool init){
    
    vector<double> n;
    vector<double> k;
    vector<double> c;
    vector<double> p_e;
    vector<double> p_c; 
    vector<double> weights;

    for (map<unsigned long, vector<int> >::iterator ci = cell_to_idx.begin();
        ci != cell_to_idx.end(); ++ci){
        
        if (assn_llr.count(ci->first) > 0 && (init || contam_rate.count(ci->first) > 0)){
            for (vector<int>::iterator i = ci->second.begin(); i != 
                ci->second.end(); ++i){
                if (init){
                   c.push_back(0.0);
                }
                else{
                    c.push_back(contam_rate[ci->first]);
                }
                n.push_back(n_all[*i]);
                k.push_back(k_all[*i]);
                p_e.push_back(p_e_all[*i]);
                double weight = 1.0;
                if (weighted){
                    weight = assn_llr[ci->first] / id_llrsum[assn[ci->first]];
                }
                weights.push_back(weight);
                if (init){
                    p_c.push_back(0.0);
                }
                else{
                    p_c.push_back(amb_mu[type1_all[*i]][type2_all[*i]]);
                }
            }
        }
    }
    
    optimML::multivar_ml_solver solver({e_r, e_a}, ll_err_rates, dll_err_rates);
    if (num_threads > 1){
        solver.set_threads(num_threads);
    }
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("c", c);
    solver.add_data("p_e", p_e);
    solver.add_data("p_c", p_c);
    solver.constrain_01(0);
    solver.constrain_01(1);
    solver.add_weights(weights);

    double this_e_r = e_r;
    double this_e_a = e_a;

    try{
        solver.solve();
    
        this_e_r = solver.results[0];
        this_e_a = solver.results[1];
    
    }
    catch (...){
        fprintf(stderr, "Error inferring error rates; keeping initial values.\n");
    }
    return make_pair(this_e_r, this_e_a);
}

/**
 * Compute log likelihood of data set with current parameters
 * Assumes pre-compiled data. This will need to be re-generated after
 * assignments are changed via reclassify_cells().
 */
double contamFinder3::compute_ll(){
    // Compute log likelihood of data set
    double loglik = 0.0;
    pair<double, double> betaparams = beta_moments(contam_cell_prior, contam_cell_prior_var);
    for (int i = 0; i < (int)n_all.size(); ++i){
        unsigned long cell = idx_to_cell[i];
        if (contam_rate.count(cell) > 0){
            double c = contam_rate[cell];
            double p_c = amb_mu[type1_all[i]][type2_all[i]];
            double binom_p;

            // Three-component for heterotypic cells that have allele_ratio
            if (allele_ratio.count(cell) > 0 && p_A_all[i] >= 0){
                double r = allele_ratio[cell];
                double p_A = adjust_p_err(p_A_all[i], e_r, e_a);
                double p_B = adjust_p_err(p_B_all[i], e_r, e_a);
                double p_endo = r * p_A + (1.0 - r) * p_B;
                binom_p = (1.0 - c) * p_endo + c * p_c;
            } else {
                double p_e = adjust_p_err(p_e_all[i], e_r, e_a);
                binom_p = (1.0 - c) * p_e + c * p_c;
            }

            if (binom_p < 1e-6) binom_p = 1e-6;
            if (binom_p > 1.0 - 1e-6) binom_p = 1.0 - 1e-6;

            loglik += logbinom(n_all[i], k_all[i], binom_p);
            loglik += (dbeta(c, betaparams.first, betaparams.second)/log2(exp(1.0)));
        }
    }
    return loglik;
}

/**
 * Find maximum likelihood estimates of all parameters of interest.
 * Return overall log likelihood.
 */
void contamFinder3::fit(){
    // Report active mechanisms
    fprintf(stderr, "Three-component model active for heterotypic tetraploid cells\n");

    // Bulk mode: simplified flow for empty-droplet profiling
    if (bulk_mode){
        fprintf(stderr, "Bulk mode: solving for ambient profile only (c fixed at 1.0)\n");

        // Skip init_params() entirely in bulk mode. init_params() runs the mixture
        // solver at a freely-optimized c, producing starting proportions estimated
        // at the wrong contamination rate. It also runs 250 BFGS restarts that are
        // completely wasted since we immediately overwrite all c values to 1.0.
        //
        // Instead, set uniform starting proportions directly. If set_init_contam_prof()
        // was called before fit() (warm-start case), those proportions are already in
        // contam_prof and contam_prof_initialized is true, so we skip this block.
        if (!contam_prof_initialized){
            for (int i = 0; i < (int)idx2samp.size(); ++i){
                int samp = idx2samp[i];
                if (inter_species){
                    contam_prof[samp] = 1.0 / (double)(idx2samp.size() + 1);
                } else {
                    contam_prof[samp] = 1.0 / (double)idx2samp.size();
                }
            }
            if (inter_species){
                contam_prof[-1] = 1.0 / (double)(idx2samp.size() + 1);
            }
            contam_prof_initialized = true;
        }

        // Set all barcodes with assignments to c=1.0 (everything is ambient).
        // Use assn (public member) to get the full barcode set.
        for (auto& a : assn){
            contam_rate[a.first] = 1.0;
        }

        // Solve for pi (mixture proportions) with c fixed at 1.0
        double dummy = -1.0;
        double loglik = this->update_amb_prof_mixture(false, dummy, false);
        fprintf(stderr, "Bulk mode log likelihood: %f\n", loglik);
        this->amb_mu_available = true;
        return;
    }

    if (use_fi_weight){
        fprintf(stderr, "FI-weight mode: continuous Fisher Information weighting active\n");
    }
    if (use_loo){
        fprintf(stderr, "LOO mode: leave-one-out ambient profiles will be computed\n");
    }
    if (user_prior_set){
        fprintf(stderr, "Fixed prior: mean=%.4f var=%.6f (overrides empirical Bayes)\n",
            contam_cell_prior, contam_cell_prior_var);
    }
    if (r_feedback_enabled){
        fprintf(stderr, "R-feedback: per-cell allele ratios will be fed back into "
            "ambient profile estimation\n");
    }
    if (adaptive_prior_enabled){
        fprintf(stderr, "Adaptive prior: will check distribution and apply fixed prior "
            "fallback if pathological (thresh=%.2f, mean=%.4f)\n",
            adaptive_prior_boundary_thresh, adaptive_prior_mean);
    }

    if (species_mode && !has_species_counts && !primary_species_counts_enabled && !fixed_amb_prof){
        fprintf(stderr, "WARNING: species-level estimation has no species-diagnostic data source.\n"
            "  Native --interspecies runs should call set_primary_species_counts_enabled(true)\n"
            "  after loading .species_counts/.species_condf as the primary inputs.\n");
    } else if (species_mode && primary_species_counts_enabled && !has_species_counts){
        fprintf(stderr, "Species-level estimation uses native .species_counts/.species_condf primary inputs.\n");
    }

    // -- Phase 1: Initialization (unchanged) --
    if (c_init <= 0){
        c_init = this->est_min_c();
    }
    double ll_init = init_params(c_init);
    fprintf(stderr, "Initial global contamination rate = %f\n", c_init);
    
    double dummy = -1.0;
    this->est_contam_cells_global();
    double loglik = this->update_amb_prof_mixture(false, dummy, true);
    
    this->amb_mu_available = true;

    if (use_loo){
        this->compute_loo_profiles();
    }

    // -- Phase 2: Empirical Bayes estimation (unchanged) --
    // First call: sets prior params from data
    this->est_contam_cells();
    // Second call: applies empirical Bayes shrinkage (or fixed prior if user_prior_set)
    this->est_contam_cells();

    if (species_mode){
        compute_loading_weights();
    }

    // -- Phase 3: Adaptive prior (if enabled) --
    // Checks the post-EB distribution. If pathological, overrides with a
    // fixed prior at adaptive_prior_mean and iteratively halves variance
    // until the distribution normalizes.
    if (adaptive_prior_enabled){
        this->run_adaptive_prior();
    }

    // -- Phase 4: R-feedback (if enabled) --
    // Uses per-cell allele ratios from est_contam_cells to re-weight
    // heterotypic cell contributions in the ambient profile mixture.
    //
    // V1_R12 / R3: make this update acceptance-gated. The species free
    // estimator can otherwise move into a coupled bad basin where a proposed
    // profile and the newly re-fit per-cell c/r reinforce each other. Save the
    // currently accepted state, propose the r-informed profile, re-fit c/r
    // against that profile, and keep it only if the comparable full objective
    // improves and remains finite.
    if (r_feedback_enabled && !allele_ratio.empty()){
        this->report_r_feedback_stats();

        map<int, double> old_contam_prof = this->contam_prof;
        map<string, double> old_species_contam_prof = this->species_contam_prof;
        map<string, double> old_species_contam_prof_conc = this->species_contam_prof_conc;
        map<string, double> old_species_prior_prof = this->species_prior_prof;
        map<string, double> old_species_prior_conc = this->species_prior_conc;
        map<pair<int, int>, map<pair<int, int>, double> > old_amb_mu = this->amb_mu;
        map<int, map<pair<int,int>, map<pair<int,int>, double> > > old_amb_mu_loo = this->amb_mu_loo;
        robin_hood::unordered_map<unsigned long, double> old_contam_rate = this->contam_rate;
        robin_hood::unordered_map<unsigned long, double> old_contam_rate_se = this->contam_rate_se;
        robin_hood::unordered_map<unsigned long, double> old_contam_rate_ll = this->contam_rate_ll;
        robin_hood::unordered_map<unsigned long, double> old_allele_ratio = this->allele_ratio;
        robin_hood::unordered_map<unsigned long, double> old_allele_ratio_se = this->allele_ratio_se;
        double old_contam_cell_prior = this->contam_cell_prior;
        double old_contam_cell_prior_var = this->contam_cell_prior_var;

        double old_obj = this->compute_ll();
        fprintf(stderr, "R-feedback acceptance gate: current objective before update = %f\n",
            old_obj);

        fprintf(stderr, "Re-estimating ambient profile with r-feedback...\n");
        double dummy2 = -1.0;
        double proposed_profile_ll = this->update_amb_prof_mixture(false, dummy2, false);

        fprintf(stderr, "Re-estimating per-cell contamination against r-informed profile...\n");
        this->est_contam_cells();

        // If adaptive prior triggered earlier, re-check after r-feedback. This
        // is part of the proposed state and will also be rolled back if the
        // objective gate rejects the proposal.
        if (adaptive_prior_enabled){
            fprintf(stderr, "Post-r-feedback adaptive prior re-check...\n");
            this->run_adaptive_prior();
        }

        double new_obj = this->compute_ll();
        double min_improvement = 1e-6;
        bool accept_r_feedback = std::isfinite(old_obj) && std::isfinite(new_obj) &&
            (new_obj > old_obj + min_improvement);

        if (accept_r_feedback){
            fprintf(stderr, "R-feedback acceptance gate: accepted update "
                "(old=%f, new=%f, delta=%f; profile_solver_ll=%f)\n",
                old_obj, new_obj, new_obj - old_obj, proposed_profile_ll);
            if (use_loo){
                fprintf(stderr, "Recomputing LOO profiles against accepted r-informed ambient...\n");
                this->compute_loo_profiles();
            }
        } else {
            fprintf(stderr, "R-feedback acceptance gate: rejected update "
                "(old=%f, new=%f, delta=%f; profile_solver_ll=%f). Restoring previous state.\n",
                old_obj, new_obj, new_obj - old_obj, proposed_profile_ll);

            this->contam_prof = old_contam_prof;
            this->species_contam_prof = old_species_contam_prof;
            this->species_contam_prof_conc = old_species_contam_prof_conc;
            this->species_prior_prof = old_species_prior_prof;
            this->species_prior_conc = old_species_prior_conc;
            this->amb_mu = old_amb_mu;
            this->amb_mu_loo = old_amb_mu_loo;
            this->contam_rate = old_contam_rate;
            this->contam_rate_se = old_contam_rate_se;
            this->contam_rate_ll = old_contam_rate_ll;
            this->allele_ratio = old_allele_ratio;
            this->allele_ratio_se = old_allele_ratio_se;
            this->contam_cell_prior = old_contam_cell_prior;
            this->contam_cell_prior_var = old_contam_cell_prior_var;
            this->amb_mu_available = true;
        }
    }

    // -- Phase 5: Reclassification (unchanged) --
    if (!skip_reassign){ 
        if (tetraploid_aware){
            fprintf(stderr, "Reclassifying cells (ploidy-aware: %lu locked, %lu safe singlets)...\n",
                locked_identities.size(), safe_singlets.size());
        } else {
            fprintf(stderr, "Reclassifying cells...\n");
        }
        bool reclassified = this->reclassify_cells();
        
        if (reclassified){
            clear_data();
            compile_data(assn, indv_allelecounts);
        }
    }
}

// ============================================================================
// Revision History
// V1_R1: Forked from ambient_rna_three.cpp V1_R3. Added r-feedback into
//        compile_amb_prof_dat (per-cell allele_ratio replaces hardcoded 0.5
//        for heterotypic cells when enabled). Added adaptive prior with
//        boundary fraction diagnostic and iterative variance shrinking.
//        New fit() flow: Phase 1 init, Phase 2 EB, Phase 3 adaptive prior,
//        Phase 4 r-feedback + ambient re-estimation, Phase 5 reclassification.
//        Both features off by default; controlled via set_r_feedback() and
//        set_adaptive_prior().
// V1_R2: Added fixed_amb_prof support (set_fixed_amb_prof, rebuild_amb_mu_from_contam_prof).
//        Added species mode (set_species_mode, set_species_prior,
//        expand_species_prior_to_indiv, solve_species_level_pi). Added bulk mode
//        (set_bulk_mode). Two-way protection between species_mode and inter_species.
//        update_amb_prof_mixture bypass for fixed_amb_prof and species mode paths.
// V1_R3: Fixed dimension mismatch crash in --species_mode estimate.
//        (1) expand_species_prior_to_indiv now filters to idx2samp-active
//            individuals only, preventing contam_prof from containing entries
//            for individuals not in the model. This fixes the bootstrap and
//            update_amb_prof_mixture dimension mismatch (contam_prof.size() vs
//            idx2samp.size() / mixfracs column count).
//        (2) solve_species_level_pi now detects empty contam_rate (first call
//            from init_params before any per-cell estimation) and falls back to
//            global c via use_global_c=true in compile_amb_prof_dat, preventing
//            degenerate c=0 data that eliminates mixture signal.
// V1_R4: Fixed remaining dimension mismatch in --species_mode estimate.
//        Root cause: est_min_c() populates contam_prof from minc_by_id which
//        contains ALL individuals observed in allele counts (up to 28), not
//        just the idx2samp subset (13 for lib13). When solve_species_level_pi
//        fails (all solver attempts throw), expand_species_prior_to_indiv is
//        never called, leaving contam_prof with 28 entries. bootstrap_amb_prof
//        then builds startprops (28) vs mixfracs (13 cols) -> mismatch crash.
//        Fixes: (a) est_min_c now filters contam_prof to idx2samp-active
//        individuals and fills missing ones with minimum fraction. (b)
//        solve_species_level_pi failure path now calls expand_species_prior_to
//        _indiv with uniform species proportions instead of leaving contam_prof
//        stale. (c) bootstrap_amb_prof dirichlet_soln indexing is bounds-checked.
// V1_R5: Species-level bootstrap in bootstrap_amb_prof. When species_mode is
//        true, the bootstrap now aggregates mixfracs to species level (same
//        remapping as solve_species_level_pi), runs the mixture solver with
//        n_species components, fits the Dirichlet at species level, then
//        expands the concentration parameters back to individual level via
//        the uniform-within-species split. Also stores species-level
//        concentrations in species_contam_prof_conc for output. Fixes the
//        dimension mismatch where individual-level bootstrap (28 components)
//        produced a degenerate Dirichlet that collapsed to 5 effective params,
//        then the indexing loop over idx2samp (28) overran the 5-element result.
//        Non-species path is unchanged.
// V1_R6: Fixed orphan mass bug in expand_species_prior_to_indiv. When a species
//        has nonzero estimated proportion but zero active individuals in idx2samp
//        (e.g. Hy at 22% with no Hy in expected_lines for lib13), the old code
//        silently dropped that mass via `continue`, causing contam_prof to sum to
//        ~78% instead of 100%. Now redistributes orphan mass proportionally across
//        species that have active individuals. species_contam_prof still records
//        the original solver-estimated proportions (including zero-active species).
// V1_R7: Hybrid fold support. PanelMetadata now carries per-(species, individual)
//        weights (default 1.0). load_panel_metadata folds Hy individuals into C
//        and B with weight 0.5 each, eliminating the collinear 5th species.
//        solve_species_level_pi, bootstrap_amb_prof, and expand_species_prior_to
//        _indiv all use weighted aggregation/expansion: an individual appearing
//        in multiple species contributes proportionally to each, and accumulates
//        proportion from each when expanding back. common.h adds get_weight()
//        convenience method and species_sample_weight map to PanelMetadata.
//        load_panel_metadata gains fold_hybrid param (default true).
// V1_R8: Added --species_init and --indiv_init warm-start features.
//        species_init_prof (map<string,double>) and species_init_used flag
//        on contamFinder3. set_species_init() stores the profile; on the first
//        solve_species_level_pi call, startprops are built from the init file
//        (with orphan mass redistribution for inactive species). Subsequent
//        calls use the previous solution. set_indiv_init() filters to idx2samp,
//        fills missing entries with 1e-6, normalizes, and sets
//        contam_prof_initialized=true so init_params uses these values.
//        quant3_contam_ap.cpp: CLI parsing for both flags, validation
//        (species_init requires species_mode estimate), file loading via
//        existing load_species_prior / load_contam_prof helpers.
// V1_R9: Added --species_counts support. New member variables:
//        has_species_counts, species_allelecounts, species_expfracs.
//        set_species_counts() stores species-diagnostic allele counts and condf.
//        solve_species_level_pi() swaps in species counts/condf before calling
//        compile_amb_prof_dat when available, swaps back after. Individual-level
//        estimation (per-cell c, bootstrap, est_contam_cells) is unchanged.
//        quant3_contam_ap.cpp: CLI --species_counts (long opt 1015), loads
//        .species_counts and companion .species_condf.
//        quant3_contam_empty_drops.cpp: CLI --species_counts (long opt 1002),
//        same loading pattern, requires --aggregate_to_species.
// V1_R10: Fixed degenerate zero-column bug in species-level solver and bootstrap.
//        When a species has no individuals in the pool (species_weight_sum == 0),
//        its mixture column is all zeros. The solver can assign arbitrary mass to
//        zero columns without changing the likelihood, causing spurious nonzero
//        estimates for absent species (e.g. Bonobo/Orangutan in human+chimp pools).
//        Fix: both solve_species_level_pi() and bootstrap_amb_prof() now build an
//        active-species index (active_sp_idx / full_to_active) and exclude inactive
//        species from the solver matrix, startprops, and MLE fracs. After solving,
//        results are mapped back to the full species list with 0.0 for inactive
//        species. The failure path in solve_species_level_pi also uses uniform
//        across active species only. Diagnostic logging reports which species are
//        active vs excluded per solver invocation.
// V1_R11: WP0 shared library changes (three changes):
//        (1) Bulk mode init fix in fit(): skip init_params() entirely in
//            bulk_mode. Set uniform starting proportions directly unless
//            contam_prof was pre-initialized via set_init_contam_prof().
//            Use assn (not contam_rate) to populate c=1.0 for all barcodes.
//        (2) Weighted within-species split: new compute_loading_weights()
//            method counts per-individual cell assignments (doublets
//            contribute 0.5 each). expand_species_prior_to_indiv() now
//            multiplies panel weights by loading weights for proportional
//            within-species distribution. Called after Phase 2 EB estimation
//            and at the end of reclassify_cells().
//        (3) Validation warning in fit(): warns when species_mode is active
//            without species-diagnostic counts and without fixed ambient
//            profile (interindividual SNPs produce poorly separated species
//            columns).
// ============================================================================
