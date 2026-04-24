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
#include "common.h"
#include "demux_vcf_llr.h"
#include "ambient_rna_three.h"

using std::cout;
using std::endl;
using namespace std;

// ============================================================================
// ambient_rna_three.cpp
// Three-component contamination model for tetraploid cells
//
// Version: V1_R3
// ============================================================================

const string AMBIENT_RNA_VERSION = "2.0-three";
const string AMBIENT_RNA_VERSION_MSG = "Three-component model: joint (r,c) estimation for heterotypic tetraploids";

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

    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != 
        assn.end(); ++a){
        this->allowed_ids.insert(a->second);
        // Make sure we also allow sub-IDs of combinations
        if (a->second >= n_samples){
            pair<int, int> combo = idx_to_hap_comb(a->second, n_samples);
            this->allowed_ids.insert(combo.first);
            this->allowed_ids.insert(combo.second);
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

    if (allowed_ids.size() == 0){
        // Default to all possible individuals
        for (int i = 0; i < n_samples; ++i){
            idx2samp.push_back(i);
        }
    }
    else{
        for (set<int>::iterator a = allowed_ids.begin(); a != allowed_ids.end(); ++a){
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

    // Compile data in the format needed by other functions
    this->compile_data(assn, indv_allelecounts);
}

void contamFinder3::set_init_contam_prof(map<int, double>& cp){
    this->contam_prof = cp;
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

        this->get_reads_expectations(a->second, indv_allelecounts[a->first],
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
void contamFinder3::get_reads_expectations(int ident,
    map<pair<int, int>, map<pair<int, int>, pair<float, float> > >& allelecounts,
    vector<double>& n,
    vector<double>& k,
    vector<double>& p_e,
    vector<double>& p_A,
    vector<double>& p_B,
    vector<pair<int, int> >& type1,
    vector<pair<int, int> >& type2){
    
    static pair<int, int> nullkey = make_pair(-1, -1);

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
    // Guard: if c_est is zero or negative, all fracs become degenerate.
    // Fall back to equal proportions.
    if (c_est <= 0.0 || isnan(c_est) || isinf(c_est)){
        for (map<int, double>::iterator mi = minc_by_id.begin(); mi != minc_by_id.end(); ++mi){
            double frac = 1.0 / (double)minc_by_id.size();
            denom += frac;
            contam_prof.insert(make_pair(mi->first, frac));
        }
    }
    else{
        for (map<int, double>::iterator mi = minc_by_id.begin(); mi != minc_by_id.end(); ++mi){
            double val = mi->second/minc_by_id_count[mi->first];
            double frac = 1.0 - val/c_est;
            if (frac < minval){
                frac = minval;
            }
            denom += frac;
            contam_prof.insert(make_pair(mi->first, frac));
        }
    }
    if (inter_species){
        contam_prof.insert(make_pair(-1, 1.0/((double)minc_by_id.size() + 1.0)));
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
                if (val > 1e-15 || val < 1.0 - 1e-15){
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
void contamFinder3::est_contam_cells(){
    
    contam_rate.clear();
    contam_rate_se.clear();
    contam_rate_ll.clear();
    allele_ratio.clear();
    allele_ratio_se.clear();

    // Store all successful estimates of contamination rate for computing the
    // mean and variance across the data set 
    vector<double> cell_c_maps;
    // Store cell weights (log likelihood ratio of most likely ID assignment)
    // to use in calculating this mean & variance
    vector<double> cell_c_llr;

    // Separate tracking for two-component (singlet/homotypic) c values only,
    // used for the Empirical Bayes prior. Heterotypic BFGS estimates on the
    // r/c likelihood ridge have inflated variance that would weaken the prior.
    vector<double> twocomp_c_maps;
    vector<double> twocomp_c_llr;

    // Track three-component stats for reporting
    int n_three_comp = 0;
    int n_two_comp = 0;
    int n_bfgs_fallback = 0;
    vector<double> cell_r_maps;

    for (map<unsigned long, vector<int> >::iterator ci = cell_to_idx.begin(); 
        ci != cell_to_idx.end(); ++ci){
        
        // Determine this cell's identity for LOO lookup
        int cell_ident = -1;
        if (use_loo && assn.count(ci->first) > 0){
            cell_ident = assn[ci->first];
            // For doublets, LOO is not applied (no single identity to exclude)
            if (cell_ident >= n_samples){
                cell_ident = -1;
            }
        }

        // Reference to the ambient profile for this cell (global or LOO)
        bool have_loo = (cell_ident >= 0 && amb_mu_loo.count(cell_ident) > 0);

        // Determine if this cell is a heterotypic combo (eligible for three-component)
        // A heterotypic combo has p_A_all >= 0 for its observations AND p_A != p_B
        // for at least some categories (i.e., the two genomes are distinguishable)
        bool is_heterotypic = false;
        bool has_any_data = false;
        if (assn.count(ci->first) > 0 && assn[ci->first] >= n_samples){
            // It's a combo assignment. Check if heterotypic (A != B)
            pair<int, int> combo = idx_to_hap_comb(assn[ci->first], n_samples);
            if (combo.first != combo.second){
                is_heterotypic = true;
            }
        }

        // Compile data for this cell
        vector<double> n;
        vector<double> k;
        vector<double> p_e;
        vector<double> p_c;
        vector<double> cell_p_A;
        vector<double> cell_p_B;
        vector<double> fi_weights;

        for (vector<int>::iterator i = ci->second.begin(); i != ci->second.end(); ++i){
            double this_p_e = adjust_p_err(p_e_all[*i], e_r, e_a);
            
            // LOO: use per-identity ambient profile if available
            double this_p_c;
            if (have_loo &&
                amb_mu_loo[cell_ident].count(type1_all[*i]) > 0 &&
                amb_mu_loo[cell_ident][type1_all[*i]].count(type2_all[*i]) > 0){
                this_p_c = amb_mu_loo[cell_ident][type1_all[*i]][type2_all[*i]];
            } else {
                this_p_c = amb_mu[type1_all[*i]][type2_all[*i]];
            }
            
            // For three-component model: get per-genome expectations
            double this_p_A = -1.0;
            double this_p_B = -1.0;
            if (is_heterotypic && p_A_all[*i] >= 0){
                this_p_A = adjust_p_err(p_A_all[*i], e_r, e_a);
                this_p_B = adjust_p_err(p_B_all[*i], e_r, e_a);
            }

            // FI-weight mode: compute continuous weight
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
                // Tetraploid-aware filtering (for two-component path)
                // Note: for three-component cells, we keep ALL categories since
                // the model can handle p_e near p_c by decomposing into p_A and p_B
                if (!is_heterotypic && tetraploid_aware){
                    bool is_combo_cat = (type2_all[*i].first != -1);
                    // Only filter combo-structured categories (homotypic tetraploids).
                    // Singlet cells have is_combo_cat == false for all categories,
                    // and their p_e = nalt/2 carries valid signal even when close to p_c.
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
        
        // If filtering removed all data for this cell, fall back to unfiltered
        if (n.empty() && (tetraploid_aware || use_fi_weight)){
            fi_weights.clear();
            cell_p_A.clear();
            cell_p_B.clear();
            for (vector<int>::iterator i = ci->second.begin(); i != ci->second.end(); ++i){
                double this_p_e = adjust_p_err(p_e_all[*i], e_r, e_a);
                double this_p_c;
                if (have_loo &&
                    amb_mu_loo[cell_ident].count(type1_all[*i]) > 0 &&
                    amb_mu_loo[cell_ident][type1_all[*i]].count(type2_all[*i]) > 0){
                    this_p_c = amb_mu_loo[cell_ident][type1_all[*i]][type2_all[*i]];
                } else {
                    this_p_c = amb_mu[type1_all[*i]][type2_all[*i]];
                }
                n.push_back(n_all[*i]);
                k.push_back(k_all[*i]);
                p_e.push_back(this_p_e);
                p_c.push_back(this_p_c);
                double this_p_A = -1.0, this_p_B = -1.0;
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

        // ================================================================
        // THREE-COMPONENT PATH: heterotypic combo cells
        // Jointly estimate (r, c) via BFGS
        // ================================================================
        if (is_heterotypic && !cell_p_A.empty() && cell_p_A[0] >= 0){
            // Initial guesses: r=0.5, c=global prior or 0.15
            double r_init = 0.5;
            double c_init_cell = (contam_cell_prior > 0) ? contam_cell_prior : 0.15;
            
            vector<double> params_init = {r_init, c_init_cell};
            optimML::multivar_ml_solver solver(params_init, ll_three, dll_three);
            
            solver.add_data("n", n);
            solver.add_data("k", k);
            solver.add_data("p_A", cell_p_A);
            solver.add_data("p_B", cell_p_B);
            solver.add_data("p_c", p_c);
            
            // Constrain both r and c to [0, 1]
            solver.constrain_01(0);  // r
            solver.constrain_01(1);  // c
            
            // FI-weight: add per-observation weights
            if (use_fi_weight && !fi_weights.empty()){
                solver.add_weights(fi_weights);
            }

            // Add Beta prior on c (param index 1) if available
            if (contam_cell_prior > 0 && contam_cell_prior_var > 0){
                pair<double, double> bm = beta_moments(contam_cell_prior, contam_cell_prior_var);
                solver.add_beta_prior(bm.first, bm.second, 1);
            }

            try {
                bool success = solver.solve();
                if (success){
                    r_cell_map = solver.results[0];
                    c_cell_map = solver.results[1];
                    ll = solver.log_likelihood;
                    // SE from BFGS: use diagonal of inverse Hessian if available
                    // The multivar solver stores se per parameter when available
                    // For now, set SE = 0 (BFGS does not always provide clean SEs)
                    se = 0.0;
                    r_se = 0.0;
                }
            }
            catch (...){
                // BFGS failed. Fall back to two-component estimate for this cell
                n_bfgs_fallback++;
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

            // Also try a few alternative starting conditions to avoid local optima
            // Try r=0.2 and r=0.8 in addition to r=0.5
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
                    solver_alt.add_beta_prior(bm.first, bm.second, 1);
                }
                try {
                    bool success = solver_alt.solve();
                    if (success && solver_alt.log_likelihood > ll){
                        r_cell_map = solver_alt.results[0];
                        c_cell_map = solver_alt.results[1];
                        ll = solver_alt.log_likelihood;
                    }
                } catch (...) {
                    // ignore
                }
            }

            n_three_comp++;
            cell_r_maps.push_back(r_cell_map);
            allele_ratio.emplace(ci->first, r_cell_map);
            allele_ratio_se.emplace(ci->first, r_se);
        }
        // ================================================================
        // TWO-COMPONENT PATH: diploid singlets and homotypic tetraploids
        // Estimate c via Brent's method (standard path)
        // ================================================================
        else {
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
            }
            catch (...){
                // pass
            }
            n_two_comp++;
        }

        contam_rate.emplace(ci->first, c_cell_map);
        cell_c_maps.push_back(c_cell_map);
        if (!is_heterotypic){
            twocomp_c_maps.push_back(c_cell_map);
        }
        if (weighted){
            double weight = assn_llr[ci->first] / id_llrsum[assn[ci->first]];
            cell_c_llr.push_back(weight);
            if (!is_heterotypic){
                twocomp_c_llr.push_back(weight);
            }
        }
        contam_rate_se.emplace(ci->first, se);
        contam_rate_ll.emplace(ci->first, ll);
    }
    
    // Re-compute data set-wide distribution (all cells, for reporting)
    pair<double, double> mu_var;
    if (weighted){
        mu_var = welford_weights(cell_c_maps, cell_c_llr, false);
    }
    else{
        mu_var = welford(cell_c_maps);
    }

    // Compute prior from two-component cells only (singlets + homotypic),
    // to avoid heterotypic ridge-variance inflating the prior.
    // Fall back to all-cell stats if there are too few two-component cells.
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
                n_bfgs_fallback, n_three_comp + n_bfgs_fallback);
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
    for (map<int, double>::iterator c = idcsum.begin(); c != idcsum.end(); ++c){
        //fprintf(stderr, "contam mean %d) %f\n", c->first, c->second/idccount[c->first]);
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
        for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator ac1 = 
            indv_allelecounts[a->first].begin(); ac1 != indv_allelecounts[a->first].end(); 
            ++ac1){
            
            if ((!is_comb && ac1->first.first == a->second) || 
                (is_comb && ac1->first.first == comb.first)){
                for (map<pair<int, int>, pair<float, float> >::iterator ac2 = 
                    ac1->second.begin(); ac2 != ac1->second.end(); ++ac2){
                    
                    if ((!is_comb && ac2->first.first == -1) || 
                        (is_comb && ac2->first.first == comb.second)){
                        
                        double expected;
                        if (!is_comb){
                            expected = adjust_p_err((double)ac1->first.second / 2.0, 
                                e_r, e_a);
                        }
                        else{
                            expected = adjust_p_err((double)(ac1->first.second + 
                                ac2->first.second)/4.0, e_r, e_a);
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
                                    mixfrac_row.push_back(0.5 * ef1 + 0.5 * ef2);
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
    compile_amb_prof_dat(solve_for_c, use_global_c, mixfracs, weights,
        n, k, p_e, c); 
    
    // Starting proportions should be those previously set
    vector<double> startprops;
    for (map<int, double>::iterator cp = contam_prof.begin(); 
        cp != contam_prof.end(); ++cp){
        if (cp->first != -1){
            startprops.push_back(cp->second);
        }
    }
    if (contam_prof.count(-1) > 0){
        startprops.push_back(contam_prof[-1]);
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

    // Compile everything we need
    vector<vector<double> > mixfracs;
    vector<double> weights;
    vector<double> n;
    vector<double> k;
    vector<double> p_e;
    vector<double> c;
    
    compile_amb_prof_dat(false, false, mixfracs, weights, n, k, p_e, c);
   
    // Store MLEs from bootstrap samples, which will serve as samples from
    // a Dirichlet distribution 
    vector<vector<double> > dirprops;
    
    // Store MLE solution from contam_prof
    vector<double> mle_fracs;

    // Starting proportions should be those previously set
    vector<double> startprops;
    for (map<int, double>::iterator cp = contam_prof.begin(); cp != contam_prof.end(); ++cp){
        if (cp->first != -1){
            startprops.push_back(cp->second);
            vector<double> v;
            dirprops.push_back(v);
            mle_fracs.push_back(cp->second);
        }
    }
    if (contam_prof.count(-1) > 0){
        startprops.push_back(contam_prof[-1]);
        vector<double> v;
        dirprops.push_back(v);
        mle_fracs.push_back(contam_prof[-1]);
    }
    
    // Initialize random stuff
    static random_device dev;
    static mt19937 rand_gen = mt19937(dev());
    static uniform_int_distribution<int> uni_dist(0, n.size()-1);
    
    for (int b = 0; b < n_boots; ++b){
        fprintf(stderr, "Bootstrap sample %d...\r", b+1);
        
        // Re-sample
        vector<vector<double> > mixfracs_boot;
        vector<double> weights_boot;
        vector<double> n_boot;
        vector<double> k_boot;
        vector<double> p_e_boot;
        vector<double> c_boot;
        
        for (int x = 0; x < n.size(); ++x){
            int r = uni_dist(rand_gen);
            mixfracs_boot.push_back(mixfracs[r]);
            weights_boot.push_back(weights[r]);
            n_boot.push_back(n[r]);
            k_boot.push_back(k[r]);
            p_e_boot.push_back(p_e[r]);
            c_boot.push_back(c[r]);
        }    
        
        // Solve
        vector<double> params;
        optimML::multivar_ml_solver solver(params, ll_amb_prof_mixture, dll_amb_prof_mixture);
        if (num_threads > 1){
            // Avoid multithreading for evaluation for mixture proportion problems
            //solver.set_threads(num_threads);
            solver.set_bfgs_threads(num_threads);
        }
        solver.add_mixcomp(mixfracs_boot);
        solver.add_mixcomp_fracs(startprops);
        solver.add_data("n", n_boot);
        solver.add_data("k", k_boot);
        solver.add_data("p_e", p_e_boot);
        solver.add_data("c", c_boot);
        
        try {
            solver.solve();
            if (std::isfinite(solver.log_likelihood)){
                for (int x = 0; x < (int)solver.results_mixcomp.size(); ++x){
                    dirprops[x].push_back(solver.results_mixcomp[x]);
                }
            }
            // else: skip this bootstrap sample (non-finite LL)
        } catch (...) {
            // Skip this bootstrap sample - solver failed on resampled data
        }
        
    }
    fprintf(stderr, "\n");
    
    // Now fit Dirichlet MLE
    vector<double> dirichlet_soln;
    fit_dirichlet(mle_fracs, dirprops, dirichlet_soln);
    
    dirichlet_params.clear();
    for (int i = 0; i < idx2samp.size(); ++i){
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
    // Begin with minimum estimate of c (global contamination rate)
    // if not provided from a previous run
    if (c_init <= 0){
        // Get (estimated/averaged) min bound on c
        c_init = this->est_min_c();
    }
    double ll_init = init_params(c_init);
    fprintf(stderr, "Initial global contamination rate = %f\n", c_init);
    
    double dummy = -1.0;
    this->est_contam_cells_global();
    double loglik = this->update_amb_prof_mixture(false, dummy, true);
    
    // Ambient profile now available for adaptive filtering
    this->amb_mu_available = true;

    // LOO: compute per-identity leave-one-out ambient profiles
    if (use_loo){
        this->compute_loo_profiles();
    }

    // Once to set prior params
    this->est_contam_cells();
    // Again to use empirical Bayes shrinkage
    this->est_contam_cells();
    
    // Allow ambient prof to update without considering individuals of origin
    //this->update_ambient_profile(false);
    
    // Update mixture components using per-cell contam estimates
    //this->update_amb_prof_mixture(false, dummy, false);
    
    if (!skip_reassign){ 
        // Update cell identities, considering contamination profile
        if (tetraploid_aware){
            fprintf(stderr, "Reclassifying cells (ploidy-aware: %lu locked, %lu safe singlets)...\n",
                locked_identities.size(), safe_singlets.size());
        } else {
            fprintf(stderr, "Reclassifying cells...\n");
        }
        bool reclassified = this->reclassify_cells();
        
        if (reclassified){
            // Allows log likelihood computation
            clear_data();
            compile_data(assn, indv_allelecounts);
        }
    }
    /* 
    // Check how low the error rates have dropped after modeling contamination
    pair<double, double> err_final = this->est_error_rates(false);
    fprintf(stderr, "Residual error rates:\n");
    fprintf(stderr, "  Reference alleles: %f\n", err_final.first);
    fprintf(stderr, "  Alt alleles: %f\n", err_final.second);
    */
}

// ============================================================================
// Revision History
// V1_R1: Initial three-component model implementation
// V1_R2: Fix ll init to -INFINITY, allele_ratio cleanup on reclassify,
//        num_threads default, split Welford prior, BFGS fallback counter
// V1_R3: Fix BFGS line search assertion failure on pathological libraries.
//        Clamp binom_p in dll_amb_prof_mixture, dll_ambmu, ll_ambmu to
//        prevent NaN gradients from division by zero. Add isfinite guard
//        on dy_dp in all gradient functions. Widen catch(int) to catch(...)
//        on all solver.solve() calls. Add post-solve isfinite checks on
//        solver.log_likelihood and solver.results to detect silent NaN
//        propagation in NDEBUG builds.
// ============================================================================
