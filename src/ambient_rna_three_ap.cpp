#include <limits>
// ============================================================================
// ambient_rna_three_ap.cpp
//
// Provenance: CellBouncer tet_contam_estimate ambient-RNA estimator
//   (contamFinder3). Part of the Tet2025 synthetic ambient-RNA benchmark.
//   Current version: V1_R23 (see Revision History block at end of file for the
//   full chronological log; the newest entries are summarized below).
//
// V1_R23 (this revision): keeps failed per-cell contamination fits unavailable
//   during profile compilation and optional reclassification instead of
//   inserting c=0, and preserves accepted state when reassignment is unresolved.
//
// V1_R22: adds endpoint-compatible continuous source exclusion,
//   explicit raw/effective mass accounting, parent-axis and conditional-information
//   diagnostics, and selector-limited r-C likelihood surfaces. The historical
//   global and hard-exclusion endpoints remain exact when the new option is absent
//   or set to 0/1 respectively.
//
// V1_R21: repairs ambient-profile bootstrap weighting while
//   preserving the complete V1_R20 model-objective instrumentation. Every
//   bootstrap replicate now receives its resampled observation-weight vector;
//   malformed and failed replicates are excluded and counted; invalid or empty
//   concentration fits are omitted; and candidate-keyed row-wise uncertainty
//   is suppressed until clustered cell-level resampling is implemented.
//
// V1_R20: reverts the withdrawn marginal self-column experiment and retains
//   the full-model likelihood instrument. The self-column hypothesis was
//   refuted by direct measurement, while final_loglik/final_loglik_valid remain
//   available so callers can compare fitted and supplied profiles on identical
//   data without changing point-estimation behavior.
//
// V1_R18: makes the deterministic free-(r,c) profile fit the production model
//   rather than a validator admitted only by BFGS success. Strict unpenalized
//   boundary, ridge, fixed-r-surrogate, regularized-boundary,
//   weighted-composition, and diagnostic safeguards are retained.
// ============================================================================
// ============================================================================
#include <algorithm>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
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
#include "genotype_llr.h"
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

const string AMBIENT_RNA_VERSION = "3.2-ck-phase1-cbr008-009";
const string AMBIENT_RNA_VERSION_MSG = "Prevents failed per-cell contamination fits from becoming implicit c=0 state and freezes unresolved reclassification without deleting supported cells";

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

// Named production safeguards for constrained contamination estimates. These
// are deliberately conservative policy thresholds. They are emitted through
// per-cell profile diagnostics and covered at/around the thresholds by the
// mathematical regression suite.
static const double C_BOUNDARY_TOLERANCE = 1e-5;
static const double MIN_UNIVARIATE_C_FISHER_INFORMATION = 1.0;
static const double MIN_BOUNDARY_VS_OPPOSITE_LL = 1.0;
static const double MIN_BOUNDARY_VS_INTERIOR_LL = 0.1;
static const double MAX_BOUNDARY_PROFILE_SUPPORT_WIDTH = 0.25;
// Raw objective gradients scale linearly with cell depth/weight and therefore
// are diagnostic only. Production acceptance is based on the independently
// profiled likelihood: evidence is allowed to grow with depth, but identical
// optima are not rejected merely because their raw gradient scale is larger.
static const double MAX_TWO_DIMENSIONAL_PROFILE_NEIGHBOR_GAIN = 0.05;
static const double MIN_INTERIOR_C_PROFILE_LL_SPAN = 0.1;
static const int INTERIOR_C_PROFILE_GRID_STEPS = 20;
static const double PROFILE_BOUNDARY_STEP = 0.01;
static const double PROFILE_SUPPORT_LL_DROP = 1.920729410347062;
static const int PROFILE_INTERVAL_BISECTION_STEPS = 48;
static const double PROFILE_INTERVAL_DOMAIN_TOLERANCE = 1e-8;

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
    
    results[0] += dy_dp * (1.0 - c) * (1.0 - p_e);
    results[1] += dy_dp * (-(1.0 - c) * p_e);
    
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
    if (!c_in_data){
        results[0] += dy_dp * (p_c - p_e);
    }
}


// ===== Tetraploid per-cell optimization helpers =====

static inline double clamp_unit_interval(double x, double eps = 1e-8){
    if (x < eps) return eps;
    if (x > 1.0 - eps) return 1.0 - eps;
    return x;
}

static inline double clamp_closed_unit_interval(double x){
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

UnivariateCInformationAssessment::UnivariateCInformationAssessment():
    identifiable(false), fisher_information(0.0), max_abs_contrast(0.0),
    status("invalid_univariate_c_inputs"){}

UnivariateCBoundaryAssessment::UnivariateCBoundaryAssessment():
    accepted(false), boundary(0), fisher_information(0.0),
    ll_boundary(-INFINITY), ll_opposite(-INFINITY),
    ll_best_interior(-INFINITY), derivative_near_boundary(
        std::numeric_limits<double>::quiet_NaN()),
    status("invalid_boundary_candidate"){}

HeterotypicCBoundaryAssessment::HeterotypicCBoundaryAssessment():
    accepted(false), boundary(0), ll_boundary(-INFINITY),
    ll_opposite(-INFINITY), ll_best_interior(-INFINITY),
    ll_candidate(-INFINITY), candidate_profile_gap(
        std::numeric_limits<double>::quiet_NaN()),
    derivative_near_boundary(std::numeric_limits<double>::quiet_NaN()),
    support_low(std::numeric_limits<double>::quiet_NaN()),
    support_high(std::numeric_limits<double>::quiet_NaN()),
    support_width(std::numeric_limits<double>::quiet_NaN()),
    nuisance_r_at_boundary(std::numeric_limits<double>::quiet_NaN()),
    status("invalid_two_dimensional_boundary_candidate"){}

ContinuousProfileIntervalAssessment::ContinuousProfileIntervalAssessment():
    success(false),
    lower(std::numeric_limits<double>::quiet_NaN()),
    upper(std::numeric_limits<double>::quiet_NaN()),
    width(std::numeric_limits<double>::quiet_NaN()),
    status("profile_interval_not_attempted"){}

TwoDimensionalFitAssessment::TwoDimensionalFitAssessment():
    accepted(false), c_boundary_candidate(false),
    profile_optimization_attempted(false),
    profile_optimization_succeeded(false), raw_candidate_available(false),
    finite_in_domain_raw_candidates(0),
    raw_candidate_r(std::numeric_limits<double>::quiet_NaN()),
    raw_candidate_c(std::numeric_limits<double>::quiet_NaN()),
    raw_candidate_objective(-INFINITY),
    projected_gradient_norm(std::numeric_limits<double>::quiet_NaN()),
    candidate_profile_gap(std::numeric_limits<double>::quiet_NaN()),
    profile_neighbor_gain(std::numeric_limits<double>::quiet_NaN()),
    profile_global_gain(std::numeric_limits<double>::quiet_NaN()),
    profile_likelihood_span(std::numeric_limits<double>::quiet_NaN()),
    profile_support_low(std::numeric_limits<double>::quiet_NaN()),
    profile_support_high(std::numeric_limits<double>::quiet_NaN()),
    profile_support_width(std::numeric_limits<double>::quiet_NaN()),
    r_profile_optimization_succeeded(false),
    r_profile_identifiable(false),
    r_profile_support_low(std::numeric_limits<double>::quiet_NaN()),
    r_profile_support_high(std::numeric_limits<double>::quiet_NaN()),
    r_profile_support_width(std::numeric_limits<double>::quiet_NaN()),
    r_validation_status("r_profile_not_attempted"),
    nuisance_r_argmax(std::numeric_limits<double>::quiet_NaN()),
    validated_r(std::numeric_limits<double>::quiet_NaN()),
    validated_c(std::numeric_limits<double>::quiet_NaN()),
    validated_objective(-INFINITY),
    status("invalid_two_dimensional_candidate"){}

static double unpenalized_c_dataset_ll(double c,
    const vector<double>& n, const vector<double>& k,
    const vector<double>& p_e, const vector<double>& p_c,
    const vector<double>& weights){

    if (n.size() != k.size() || n.size() != p_e.size() ||
        n.size() != p_c.size() ||
        (!weights.empty() && weights.size() != n.size())) return -INFINITY;

    double value = 0.0;
    for (size_t i = 0; i < n.size(); ++i){
        const double w = weights.empty() ? 1.0 : weights[i];
        if (!(w > 0.0) || !(n[i] > 0.0)) continue;
        double p = (1.0-c)*p_e[i] + c*p_c[i];
        p = clamp_unit_interval(p, 1e-9);
        const double term = w*logbinom(n[i], k[i], p);
        if (!std::isfinite(term)) return -INFINITY;
        value += term;
    }
    return value;
}

static double unpenalized_c_dataset_gradient(double c,
    const vector<double>& n, const vector<double>& k,
    const vector<double>& p_e, const vector<double>& p_c,
    const vector<double>& weights){

    if (n.size() != k.size() || n.size() != p_e.size() ||
        n.size() != p_c.size() ||
        (!weights.empty() && weights.size() != n.size())){
        return std::numeric_limits<double>::quiet_NaN();
    }

    double value = 0.0;
    map<string,int> data_i;
    for (size_t i = 0; i < n.size(); ++i){
        const double w = weights.empty() ? 1.0 : weights[i];
        if (!(w > 0.0) || !(n[i] > 0.0)) continue;
        map<string,double> data_d;
        data_d["n"] = n[i];
        data_d["k"] = k[i];
        data_d["p_e"] = p_e[i];
        data_d["p_c"] = p_c[i];
        const double term = w*dll_dc(c, data_d, data_i);
        if (!std::isfinite(term)) return std::numeric_limits<double>::quiet_NaN();
        value += term;
    }
    return value;
}

UnivariateCInformationAssessment assess_univariate_c_information(
    const vector<double>& n, const vector<double>& p_e,
    const vector<double>& p_c, const vector<double>& weights){

    UnivariateCInformationAssessment out;
    if (n.empty() || n.size() != p_e.size() || n.size() != p_c.size() ||
        (!weights.empty() && weights.size() != n.size())) return out;

    bool saw_observation = false;
    for (size_t i = 0; i < n.size(); ++i){
        const double w = weights.empty() ? 1.0 : weights[i];
        if (!(w > 0.0) || !(n[i] > 0.0) || !std::isfinite(w) ||
            !std::isfinite(n[i]) || !std::isfinite(p_e[i]) ||
            !std::isfinite(p_c[i])) continue;
        saw_observation = true;
        const double contrast = p_c[i] - p_e[i];
        out.max_abs_contrast = std::max(out.max_abs_contrast, fabs(contrast));
        const double p_mid = clamp_unit_interval(0.5*(p_e[i]+p_c[i]), 1e-6);
        out.fisher_information +=
            w*n[i]*contrast*contrast/(p_mid*(1.0-p_mid));
    }

    if (!saw_observation){
        out.status = "no_weighted_univariate_observations";
        return out;
    }
    if (out.max_abs_contrast <= 1e-10){
        out.status = "unidentifiable_flat_c_likelihood";
        return out;
    }
    // One unit of approximate Fisher information is a deliberately conservative
    // minimum for producing any scalar c estimate.  Smaller values represent a
    // nearly flat likelihood whose endpoint choice is numerically unstable.
    if (!std::isfinite(out.fisher_information) ||
        out.fisher_information < MIN_UNIVARIATE_C_FISHER_INFORMATION){
        out.status = "unidentifiable_weak_c_likelihood";
        return out;
    }
    out.identifiable = true;
    out.status = "identifiable_c_likelihood";
    return out;
}

UnivariateCBoundaryAssessment assess_unpenalized_c_boundary(
    double estimate, const vector<double>& n,
    const vector<double>& k, const vector<double>& p_e,
    const vector<double>& p_c, const vector<double>& weights,
    bool regularized_candidate){

    UnivariateCBoundaryAssessment out;
    if (regularized_candidate){
        out.status = "regularized_boundary_not_selectable";
        return out;
    }
    if (!std::isfinite(estimate) || n.size() != k.size()) return out;

    const double boundary_tolerance = C_BOUNDARY_TOLERANCE;
    if (estimate <= boundary_tolerance) out.boundary = -1;
    else if (estimate >= 1.0-boundary_tolerance) out.boundary = 1;
    else {
        out.status = "not_a_boundary_candidate";
        return out;
    }

    UnivariateCInformationAssessment info =
        assess_univariate_c_information(n, p_e, p_c, weights);
    out.fisher_information = info.fisher_information;
    if (!info.identifiable){
        out.status = info.status;
        return out;
    }

    const double lower_ll = unpenalized_c_dataset_ll(0.0, n, k, p_e, p_c, weights);
    const double upper_ll = unpenalized_c_dataset_ll(1.0, n, k, p_e, p_c, weights);
    out.ll_boundary = out.boundary < 0 ? lower_ll : upper_ll;
    out.ll_opposite = out.boundary < 0 ? upper_ll : lower_ll;

    static const double interior_grid[] = {
        0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95
    };
    for (double c : interior_grid){
        out.ll_best_interior = std::max(out.ll_best_interior,
            unpenalized_c_dataset_ll(c, n, k, p_e, p_c, weights));
    }
    if (!std::isfinite(out.ll_boundary) || !std::isfinite(out.ll_opposite) ||
        !std::isfinite(out.ll_best_interior)){
        out.status = "nonfinite_boundary_likelihood";
        return out;
    }

    // Require likelihood-ratio support, not merely a finite endpoint selected by
    // the vendored Brent tie-breaking rule.  The boundary must beat the opposite
    // endpoint by at least one log-likelihood unit and the nearest tested
    // interior alternatives by at least 0.1 log-likelihood units.
    if (out.ll_boundary-out.ll_opposite < MIN_BOUNDARY_VS_OPPOSITE_LL ||
        out.ll_boundary-out.ll_best_interior < MIN_BOUNDARY_VS_INTERIOR_LL){
        out.status = "boundary_likelihood_support_too_weak";
        return out;
    }

    const double epsilon = 1e-5;
    const double near_boundary = out.boundary < 0 ? epsilon : 1.0-epsilon;
    out.derivative_near_boundary = unpenalized_c_dataset_gradient(
        near_boundary, n, k, p_e, p_c, weights);
    const double mid_derivative = unpenalized_c_dataset_gradient(
        0.5, n, k, p_e, p_c, weights);
    if (!std::isfinite(out.derivative_near_boundary) ||
        !std::isfinite(mid_derivative)){
        out.status = "nonfinite_boundary_derivative";
        return out;
    }

    const double derivative_tolerance = 1e-8;
    const bool kkt_supported = out.boundary < 0 ?
        (out.derivative_near_boundary < -derivative_tolerance &&
         mid_derivative < -derivative_tolerance) :
        (out.derivative_near_boundary > derivative_tolerance &&
         mid_derivative > derivative_tolerance);
    if (!kkt_supported){
        out.status = "boundary_kkt_direction_not_supported";
        return out;
    }

    out.accepted = true;
    out.status = out.boundary < 0 ?
        "supported_unpenalized_lower_boundary" :
        "supported_unpenalized_upper_boundary";
    return out;
}

HeterotypicFitSource choose_heterotypic_fit_source(
    int successful_two_dimensional_starts,
    bool diagnostic_fixed_r_candidate_available){
    if (successful_two_dimensional_starts > 0){
        return HET_FIT_TWO_DIMENSIONAL;
    }
    if (diagnostic_fixed_r_candidate_available){
        return HET_FIT_DIAGNOSTIC_FIXED_R;
    }
    return HET_FIT_NONE;
}

bool heterotypic_fit_source_is_production_valid(HeterotypicFitSource source){
    return source == HET_FIT_TWO_DIMENSIONAL;
}

bool regularized_c_candidate_is_selectable(bool fit_success, double c){
    const double boundary_tolerance = C_BOUNDARY_TOLERANCE;
    return fit_success && std::isfinite(c) &&
        c > boundary_tolerance && c < 1.0-boundary_tolerance;
}

static bool beta_params_from_moments_safe(double mean, double var,
    double& alpha, double& beta){
    if (!std::isfinite(mean) || !std::isfinite(var) || mean <= 0.0 || mean >= 1.0 || var <= 0.0){
        return false;
    }
    const double vmax = mean * (1.0 - mean);
    if (var >= vmax){
        var = vmax * 0.95;
    }
    if (var <= 1e-12){
        var = 1e-12;
    }
    const double common = mean * (1.0 - mean) / var - 1.0;
    alpha = mean * common;
    beta = (1.0 - mean) * common;
    return std::isfinite(alpha) && std::isfinite(beta) && alpha > 0.0 && beta > 0.0;
}

bool contamination_beta_prior_is_log_concave(
    double mean, double var, double& alpha, double& beta){
    if (!beta_params_from_moments_safe(mean, var, alpha, beta)) return false;
    return alpha >= 1.0 && beta >= 1.0;
}

static double beta_log_density_no_const(double x, double mean, double var){
    double alpha = 0.0;
    double beta = 0.0;
    if (!beta_params_from_moments_safe(mean, var, alpha, beta)){
        return 0.0;
    }
    x = clamp_unit_interval(x);
    return (alpha - 1.0) * log(x) + (beta - 1.0) * log(1.0 - x);
}

static double beta_log_density_full(double x, double mean, double var){
    double alpha = 0.0;
    double beta = 0.0;
    if (!beta_params_from_moments_safe(mean, var, alpha, beta)){
        return 0.0;
    }
    x = clamp_unit_interval(x);
    return (alpha - 1.0) * log(x) + (beta - 1.0) * log(1.0 - x)
        - (lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta));
}

static double beta_log_gradient(double x, double mean, double var){
    double alpha = 0.0;
    double beta = 0.0;
    if (!beta_params_from_moments_safe(mean, var, alpha, beta)){
        return 0.0;
    }
    x = clamp_unit_interval(x);
    return (alpha - 1.0) / x - (beta - 1.0) / (1.0 - x);
}

static double eval_three_dataset_objective(double r, double c,
    const vector<double>& n, const vector<double>& k,
    const vector<double>& p_A, const vector<double>& p_B,
    const vector<double>& p_c, const vector<double>& obs_weights,
    bool apply_prior, double prior_mean, double prior_var){

    r = clamp_closed_unit_interval(r);
    c = clamp_closed_unit_interval(c);
    double value = 0.0;
    for (size_t i = 0; i < n.size(); ++i){
        const double p_endo = r * p_A[i] + (1.0 - r) * p_B[i];
        double p = (1.0 - c) * p_endo + c * p_c[i];
        p = clamp_unit_interval(p, 1e-6);
        const double w = obs_weights.empty() ? 1.0 : obs_weights[i];
        value += w * logbinom(n[i], k[i], p);
    }
    if (apply_prior){
        value += beta_log_density_no_const(c, prior_mean, prior_var);
    }
    return value;
}

static void eval_three_dataset_gradient(double r, double c,
    const vector<double>& n, const vector<double>& k,
    const vector<double>& p_A, const vector<double>& p_B,
    const vector<double>& p_c, const vector<double>& obs_weights,
    bool apply_prior, double prior_mean, double prior_var,
    double& grad_r, double& grad_c){

    r = clamp_closed_unit_interval(r);
    c = clamp_closed_unit_interval(c);
    grad_r = 0.0;
    grad_c = 0.0;
    for (size_t i = 0; i < n.size(); ++i){
        const double p_endo = r * p_A[i] + (1.0 - r) * p_B[i];
        double p = (1.0 - c) * p_endo + c * p_c[i];
        p = clamp_unit_interval(p, 1e-6);
        const double dy_dp = (k[i] - n[i] * p) / (p - p * p);
        const double w = obs_weights.empty() ? 1.0 : obs_weights[i];
        grad_r += w * dy_dp * (1.0 - c) * (p_A[i] - p_B[i]);
        grad_c += w * dy_dp * (p_c[i] - p_endo);
    }
    if (apply_prior){
        grad_c += beta_log_gradient(c, prior_mean, prior_var);
    }
}

static double maximize_profile_coordinate(bool fix_c, double fixed_value,
    const vector<double>& n, const vector<double>& k,
    const vector<double>& p_A, const vector<double>& p_B,
    const vector<double>& p_c, const vector<double>& obs_weights,
    bool apply_prior, double prior_mean, double prior_var,
    double& nuisance_argmax){

    // The binomial log likelihood is concave in the remaining scalar coordinate
    // for fixed values of the other coordinate. Search the complete constrained
    // domain [0,1], refine the interior deterministically, and explicitly compare
    // both exact endpoints. This is required for genuine r=0/r=1 solutions.
    const double lo0 = 0.0;
    const double hi0 = 1.0;
    const double gr = 0.6180339887498948482;
    auto f = [&](double x){
        return fix_c
            ? eval_three_dataset_objective(x, fixed_value, n, k, p_A, p_B, p_c,
                obs_weights, apply_prior, prior_mean, prior_var)
            : eval_three_dataset_objective(fixed_value, x, n, k, p_A, p_B, p_c,
                obs_weights, apply_prior, prior_mean, prior_var);
    };

    double lo=lo0;
    double hi=hi0;
    double x1 = hi - gr * (hi - lo);
    double x2 = lo + gr * (hi - lo);
    double f1 = f(x1);
    double f2 = f(x2);
    for (int iter = 0; iter < 48; ++iter){
        if (f1 < f2){
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + gr * (hi - lo);
            f2 = f(x2);
        } else {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = hi - gr * (hi - lo);
            f1 = f(x1);
        }
    }

    const double f0=f(lo0);
    const double fmid=f(0.5);
    const double fhi=f(hi0);
    nuisance_argmax=x1;
    double best=f1;
    if (f2>best){nuisance_argmax=x2;best=f2;}
    if (f0>best){nuisance_argmax=lo0;best=f0;}
    if (fmid>best){nuisance_argmax=0.5;best=fmid;}
    if (fhi>best){nuisance_argmax=hi0;best=fhi;}

    // On an exactly flat nuisance profile, return the neutral midpoint rather
    // than a golden-section artifact. The downstream continuous r-profile
    // interval still marks the coordinate unidentifiable and suppresses output.
    const double fmin=std::min(std::min(f0,fmid),fhi);
    const double fmax=std::max(std::max(f0,fmid),fhi);
    if (std::isfinite(fmin) && std::isfinite(fmax) && fmax-fmin<=1e-10){
        nuisance_argmax=0.5;
        return fmid;
    }
    return best;
}

static double maximize_c_profile_interval(double lo, double hi,
    const vector<double>& n, const vector<double>& k,
    const vector<double>& p_A, const vector<double>& p_B,
    const vector<double>& p_c, const vector<double>& obs_weights,
    bool apply_prior, double prior_mean, double prior_var,
    double& c_argmax, double& r_argmax){

    lo=clamp_closed_unit_interval(lo);
    hi=clamp_closed_unit_interval(hi);
    if (hi<lo) std::swap(lo,hi);
    const double original_lo=lo;
    const double original_hi=hi;
    auto profile = [&](double c_value, double& nuisance_r){
        return maximize_profile_coordinate(true,c_value,n,k,p_A,p_B,p_c,
            obs_weights,apply_prior,prior_mean,prior_var,nuisance_r);
    };
    if (hi-lo<=1e-12){
        c_argmax=lo;
        return profile(lo,r_argmax);
    }

    const double gr=0.6180339887498948482;
    double x1=hi-gr*(hi-lo);
    double x2=lo+gr*(hi-lo);
    double r1=0.5,r2=0.5;
    double f1=profile(x1,r1);
    double f2=profile(x2,r2);
    for (int iter=0;iter<28;++iter){
        if (f1<f2){
            lo=x1;
            x1=x2; f1=f2; r1=r2;
            x2=lo+gr*(hi-lo);
            f2=profile(x2,r2);
        } else {
            hi=x2;
            x2=x1; f2=f1; r2=r1;
            x1=hi-gr*(hi-lo);
            f1=profile(x1,r1);
        }
    }

    // Compare the refined interior points with both bracket endpoints so a
    // genuine constrained endpoint is not hidden by golden-section interior
    // sampling.
    double rlo=0.5,rhi=0.5;
    const double flo=profile(original_lo,rlo);
    const double fhi=profile(original_hi,rhi);
    c_argmax=x1; r_argmax=r1; double best=f1;
    if (f2>best){c_argmax=x2;r_argmax=r2;best=f2;}
    if (flo>best){c_argmax=original_lo;r_argmax=rlo;best=flo;}
    if (fhi>best){c_argmax=original_hi;r_argmax=rhi;best=fhi;}
    return best;
}

static ContinuousProfileIntervalAssessment continuous_profile_interval(
    bool profile_c, double mle_r, double mle_c, double best_ll,
    const vector<double>& n, const vector<double>& k,
    const vector<double>& p_A, const vector<double>& p_B,
    const vector<double>& p_c, const vector<double>& obs_weights,
    bool apply_prior, double prior_mean, double prior_var){

    ContinuousProfileIntervalAssessment out;
    if (!std::isfinite(mle_r) || !std::isfinite(mle_c) ||
        !std::isfinite(best_ll) || n.empty() || n.size()!=k.size() ||
        n.size()!=p_A.size() || n.size()!=p_B.size() ||
        n.size()!=p_c.size() ||
        (!obs_weights.empty() && obs_weights.size()!=n.size())){
        out.status="invalid_continuous_profile_interval_inputs";
        return out;
    }

    const double mle=clamp_closed_unit_interval(profile_c ? mle_c : mle_r);
    auto profile_ll = [&](double value){
        double nuisance=0.5;
        return maximize_profile_coordinate(profile_c,
            clamp_closed_unit_interval(value),n,k,p_A,p_B,p_c,obs_weights,
            apply_prior,prior_mean,prior_var,nuisance);
    };

    const double ll_at_mle=profile_ll(mle);
    if (!std::isfinite(ll_at_mle)){
        out.status="nonfinite_profile_likelihood_at_optimum";
        return out;
    }
    const double peak_ll=std::max(best_ll,ll_at_mle);
    const double threshold=peak_ll-PROFILE_SUPPORT_LL_DROP;

    auto solve_lower = [&](){
        const double ll0=profile_ll(0.0);
        if (!std::isfinite(ll0)) return std::numeric_limits<double>::quiet_NaN();
        if (ll0>=threshold || mle<=PROFILE_INTERVAL_DOMAIN_TOLERANCE) return 0.0;
        double below=0.0;
        double above=mle;
        for (int iter=0;iter<PROFILE_INTERVAL_BISECTION_STEPS;++iter){
            const double mid=0.5*(below+above);
            const double ll=profile_ll(mid);
            if (!std::isfinite(ll)) return std::numeric_limits<double>::quiet_NaN();
            if (ll>=threshold) above=mid;
            else below=mid;
        }
        return above;
    };
    auto solve_upper = [&](){
        const double ll1=profile_ll(1.0);
        if (!std::isfinite(ll1)) return std::numeric_limits<double>::quiet_NaN();
        if (ll1>=threshold || mle>=1.0-PROFILE_INTERVAL_DOMAIN_TOLERANCE) return 1.0;
        double above=mle;
        double below=1.0;
        for (int iter=0;iter<PROFILE_INTERVAL_BISECTION_STEPS;++iter){
            const double mid=0.5*(above+below);
            const double ll=profile_ll(mid);
            if (!std::isfinite(ll)) return std::numeric_limits<double>::quiet_NaN();
            if (ll>=threshold) above=mid;
            else below=mid;
        }
        return above;
    };

    out.lower=solve_lower();
    out.upper=solve_upper();
    out.width=out.upper-out.lower;
    if (!std::isfinite(out.lower) || !std::isfinite(out.upper) ||
        !std::isfinite(out.width)){
        out.status="nonfinite_continuous_profile_interval";
        return out;
    }
    if (out.lower < -PROFILE_INTERVAL_DOMAIN_TOLERANCE ||
        out.upper > 1.0+PROFILE_INTERVAL_DOMAIN_TOLERANCE ||
        out.lower > mle+PROFILE_INTERVAL_DOMAIN_TOLERANCE ||
        out.upper < mle-PROFILE_INTERVAL_DOMAIN_TOLERANCE ||
        out.width < -PROFILE_INTERVAL_DOMAIN_TOLERANCE ||
        out.width > 1.0+PROFILE_INTERVAL_DOMAIN_TOLERANCE){
        out.status="invalid_continuous_profile_interval_invariants";
        return out;
    }
    out.lower=clamp_closed_unit_interval(out.lower);
    out.upper=clamp_closed_unit_interval(out.upper);
    out.width=out.upper-out.lower;
    out.success=true;
    out.status="continuous_profile_interval_succeeded";
    return out;
}

static void assess_r_profile_identifiability(TwoDimensionalFitAssessment& out,
    const vector<double>& n, const vector<double>& k,
    const vector<double>& p_A, const vector<double>& p_B,
    const vector<double>& p_c, const vector<double>& obs_weights,
    bool apply_prior, double prior_mean, double prior_var){

    out.r_profile_optimization_succeeded=false;
    out.r_profile_identifiable=false;
    out.r_validation_status="r_profile_not_attempted";
    if (!std::isfinite(out.validated_r) || !std::isfinite(out.validated_c) ||
        !std::isfinite(out.validated_objective)){
        out.r_validation_status="invalid_r_profile_inputs";
        return;
    }
    ContinuousProfileIntervalAssessment support=continuous_profile_interval(
        false,out.validated_r,out.validated_c,out.validated_objective,n,k,p_A,
        p_B,p_c,obs_weights,apply_prior,prior_mean,prior_var);
    if (!support.success){
        out.r_validation_status=string("r_")+support.status;
        return;
    }
    out.r_profile_optimization_succeeded=true;
    out.nuisance_r_argmax=out.validated_r;
    out.r_profile_support_low=support.lower;
    out.r_profile_support_high=support.upper;
    out.r_profile_support_width=support.width;
    const bool full_domain = support.lower<=PROFILE_INTERVAL_DOMAIN_TOLERANCE &&
        support.upper>=1.0-PROFILE_INTERVAL_DOMAIN_TOLERANCE;
    if (full_domain){
        out.r_validation_status="unidentifiable_full_domain_r_profile";
        out.validated_r=std::numeric_limits<double>::quiet_NaN();
        return;
    }
    out.r_profile_identifiable=true;
    out.r_validation_status="r_profile_identifiable";
}

static double projected_gradient_violation(double value, double gradient){
    if (!std::isfinite(value) || !std::isfinite(gradient)) return INFINITY;
    if (value <= C_BOUNDARY_TOLERANCE){
        // At a lower bound, a positive derivative points into the feasible
        // interval and therefore violates the maximization KKT condition.
        return std::max(0.0, gradient);
    }
    if (value >= 1.0-C_BOUNDARY_TOLERANCE){
        // At an upper bound, a negative derivative points back into the
        // feasible interval and therefore violates the KKT condition.
        return std::max(0.0, -gradient);
    }
    return fabs(gradient);
}

HeterotypicCBoundaryAssessment assess_unpenalized_heterotypic_c_boundary(
    double candidate_r, double candidate_c,
    const vector<double>& n, const vector<double>& k,
    const vector<double>& p_A, const vector<double>& p_B,
    const vector<double>& p_c, const vector<double>& weights,
    bool regularized_candidate){

    HeterotypicCBoundaryAssessment out;
    if (regularized_candidate){
        out.status = "regularized_two_dimensional_boundary_not_selectable";
        return out;
    }
    if (n.empty() || n.size()!=k.size() || n.size()!=p_A.size() ||
        n.size()!=p_B.size() || n.size()!=p_c.size() ||
        (!weights.empty() && weights.size()!=n.size()) ||
        !std::isfinite(candidate_r) || !std::isfinite(candidate_c)){
        out.status = "invalid_two_dimensional_profile_inputs";
        return out;
    }
    if (candidate_c <= C_BOUNDARY_TOLERANCE) out.boundary = -1;
    else if (candidate_c >= 1.0-C_BOUNDARY_TOLERANCE) out.boundary = 1;
    else {
        out.status = "not_a_two_dimensional_c_boundary_candidate";
        return out;
    }

    auto profile_ll = [&](double c_value, double& nuisance_r){
        return maximize_profile_coordinate(true, c_value, n, k, p_A, p_B,
            p_c, weights, false, -1.0, -1.0, nuisance_r);
    };

    double r_lower=0.5, r_upper=0.5;
    const double ll_lower = profile_ll(0.0, r_lower);
    const double ll_upper = profile_ll(1.0, r_upper);
    out.ll_boundary = out.boundary < 0 ? ll_lower : ll_upper;
    out.ll_opposite = out.boundary < 0 ? ll_upper : ll_lower;
    out.nuisance_r_at_boundary = out.boundary < 0 ? r_lower : r_upper;
    out.ll_candidate = eval_three_dataset_objective(candidate_r,candidate_c,
        n,k,p_A,p_B,p_c,weights,false,-1.0,-1.0);
    out.candidate_profile_gap = out.ll_boundary-out.ll_candidate;

    vector<double> profile_grid;
    for (int i=0;i<=40;++i) profile_grid.push_back((double)i/40.0);
    profile_grid.push_back(PROFILE_BOUNDARY_STEP);
    profile_grid.push_back(1.0-PROFILE_BOUNDARY_STEP);
    sort(profile_grid.begin(),profile_grid.end());
    profile_grid.erase(unique(profile_grid.begin(),profile_grid.end()),
        profile_grid.end());

    vector<pair<double,double> > profiled;
    double ll_min=INFINITY;
    double ll_max=-INFINITY;
    for (size_t i=0;i<profile_grid.size();++i){
        const double cc=profile_grid[i];
        double rr=0.5;
        const double ll=profile_ll(cc,rr);
        if (!std::isfinite(ll)){
            out.status = "nonfinite_two_dimensional_profile_likelihood";
            return out;
        }
        profiled.push_back(make_pair(cc,ll));
        ll_min=std::min(ll_min,ll);
        ll_max=std::max(ll_max,ll);
        if (cc > 0.0 && cc < 1.0){
            out.ll_best_interior=std::max(out.ll_best_interior,ll);
        }
    }
    if (!std::isfinite(out.ll_boundary) || !std::isfinite(out.ll_opposite) ||
        !std::isfinite(out.ll_best_interior) || !std::isfinite(out.ll_candidate)){
        out.status = "nonfinite_two_dimensional_boundary_likelihood";
        return out;
    }
    if (ll_max-ll_min <= 1e-8){
        out.status = "unidentifiable_flat_profile_c_likelihood";
        return out;
    }
    if (out.ll_boundary-out.ll_opposite < MIN_BOUNDARY_VS_OPPOSITE_LL ||
        out.ll_boundary-out.ll_best_interior < MIN_BOUNDARY_VS_INTERIOR_LL){
        out.status = "two_dimensional_boundary_profile_support_too_weak";
        return out;
    }

    double r_near=0.5;
    const double c_near = out.boundary < 0 ? PROFILE_BOUNDARY_STEP :
        1.0-PROFILE_BOUNDARY_STEP;
    const double ll_near = profile_ll(c_near,r_near);
    out.derivative_near_boundary = out.boundary < 0 ?
        (ll_near-out.ll_boundary)/PROFILE_BOUNDARY_STEP :
        (out.ll_boundary-ll_near)/PROFILE_BOUNDARY_STEP;
    const bool kkt_supported = out.boundary < 0 ?
        (out.derivative_near_boundary < -1e-8) :
        (out.derivative_near_boundary > 1e-8);
    if (!std::isfinite(out.derivative_near_boundary) || !kkt_supported){
        out.status = "two_dimensional_boundary_profile_kkt_not_supported";
        return out;
    }

    ContinuousProfileIntervalAssessment support=continuous_profile_interval(
        true,out.nuisance_r_at_boundary,out.boundary<0 ? 0.0 : 1.0,
        out.ll_boundary,n,k,p_A,p_B,p_c,weights,false,-1.0,-1.0);
    if (!support.success){
        out.status=string("two_dimensional_boundary_")+support.status;
        return out;
    }
    out.support_low=support.lower;
    out.support_high=support.upper;
    out.support_width=support.width;
    if (out.support_width > MAX_BOUNDARY_PROFILE_SUPPORT_WIDTH){
        out.status = "two_dimensional_boundary_profile_too_broad";
        return out;
    }

    out.accepted=true;
    out.status=out.boundary < 0 ?
        "profile_validated_two_dimensional_lower_boundary" :
        "profile_validated_two_dimensional_upper_boundary";
    return out;
}

TwoDimensionalFitAssessment fit_two_dimensional_profile_model(
    const vector<pair<double, double> >& raw_candidates,
    const vector<double>& n, const vector<double>& k,
    const vector<double>& p_A, const vector<double>& p_B,
    const vector<double>& p_c, const vector<double>& weights,
    bool apply_prior, double prior_mean, double prior_var){

    TwoDimensionalFitAssessment out;
    out.profile_optimization_attempted=true;
    if (n.empty() || n.size()!=k.size() || n.size()!=p_A.size() ||
        n.size()!=p_B.size() || n.size()!=p_c.size() ||
        (!weights.empty() && weights.size()!=n.size())){
        out.status="invalid_two_dimensional_profile_inputs";
        return out;
    }

    // BFGS is an optional accelerator and diagnostic source only. Recompute the
    // objective for every finite in-domain coordinate and retain the best raw
    // candidate, but continue to the deterministic free-model profile even when
    // this list is empty or every supplied coordinate is malformed.
    double best_raw_obj=-INFINITY;
    double bounded_raw_r=std::numeric_limits<double>::quiet_NaN();
    double bounded_raw_c=std::numeric_limits<double>::quiet_NaN();
    for (size_t i=0;i<raw_candidates.size();++i){
        const double rr=raw_candidates[i].first;
        const double cc=raw_candidates[i].second;
        if (!std::isfinite(rr) || !std::isfinite(cc) ||
            rr < -C_BOUNDARY_TOLERANCE || rr > 1.0+C_BOUNDARY_TOLERANCE ||
            cc < -C_BOUNDARY_TOLERANCE || cc > 1.0+C_BOUNDARY_TOLERANCE){
            continue;
        }
        const double br=clamp_closed_unit_interval(rr);
        const double bc=clamp_closed_unit_interval(cc);
        const double obj=eval_three_dataset_objective(br,bc,n,k,p_A,p_B,p_c,
            weights,apply_prior,prior_mean,prior_var);
        if (!std::isfinite(obj)) continue;
        out.finite_in_domain_raw_candidates++;
        if (!out.raw_candidate_available || obj>best_raw_obj){
            out.raw_candidate_available=true;
            bounded_raw_r=br;
            bounded_raw_c=bc;
            best_raw_obj=obj;
            out.raw_candidate_r=rr;
            out.raw_candidate_c=cc;
            out.raw_candidate_objective=obj;
        }
    }

    double profile_at_raw=std::numeric_limits<double>::quiet_NaN();
    if (out.raw_candidate_available){
        double nuisance_r=0.5;
        profile_at_raw=maximize_profile_coordinate(true,bounded_raw_c,n,k,
            p_A,p_B,p_c,weights,apply_prior,prior_mean,prior_var,nuisance_r);
        out.candidate_profile_gap=profile_at_raw-best_raw_obj;
        if (!std::isfinite(out.candidate_profile_gap) ||
            !std::isfinite(profile_at_raw) || !std::isfinite(nuisance_r)){
            // The raw candidate is not required. Discard its diagnostics and
            // continue with the unconditional deterministic profile.
            out.raw_candidate_available=false;
            out.raw_candidate_r=std::numeric_limits<double>::quiet_NaN();
            out.raw_candidate_c=std::numeric_limits<double>::quiet_NaN();
            out.raw_candidate_objective=-INFINITY;
            out.candidate_profile_gap=std::numeric_limits<double>::quiet_NaN();
            profile_at_raw=std::numeric_limits<double>::quiet_NaN();
        }
    }

    out.c_boundary_candidate=false;
    vector<double> profile_grid;
    profile_grid.reserve(INTERIOR_C_PROFILE_GRID_STEPS+4);
    for (int i=0;i<=INTERIOR_C_PROFILE_GRID_STEPS;++i){
        profile_grid.push_back((double)i/(double)INTERIOR_C_PROFILE_GRID_STEPS);
    }
    if (out.raw_candidate_available) profile_grid.push_back(bounded_raw_c);
    sort(profile_grid.begin(),profile_grid.end());
    profile_grid.erase(unique(profile_grid.begin(),profile_grid.end()),
        profile_grid.end());

    struct ProfilePoint { double c; double ll; double r; };
    vector<ProfilePoint> profiled;
    double profile_min=INFINITY;
    double profile_max=-INFINITY;
    size_t best_index=0;
    for (size_t i=0;i<profile_grid.size();++i){
        double rr=0.5;
        const double ll=maximize_profile_coordinate(true,profile_grid[i],n,k,
            p_A,p_B,p_c,weights,apply_prior,prior_mean,prior_var,rr);
        if (!std::isfinite(ll) || !std::isfinite(rr)){
            out.status="nonfinite_two_dimensional_profile_likelihood";
            return out;
        }
        ProfilePoint pp={profile_grid[i],ll,rr};
        profiled.push_back(pp);
        profile_min=std::min(profile_min,ll);
        if (ll>profile_max){profile_max=ll;best_index=i;}
    }
    out.profile_optimization_succeeded=true;
    out.profile_likelihood_span=profile_max-profile_min;
    if (!std::isfinite(out.profile_likelihood_span) ||
        out.profile_likelihood_span <= MIN_INTERIOR_C_PROFILE_LL_SPAN){
        out.status="unidentifiable_nearly_flat_profile_c_likelihood";
        return out;
    }

    const double bracket_lo = best_index==0 ? profiled[0].c :
        profiled[best_index-1].c;
    const double bracket_hi = best_index+1>=profiled.size() ?
        profiled.back().c : profiled[best_index+1].c;
    double refined_c=profiled[best_index].c;
    double refined_r=profiled[best_index].r;
    double refined_ll=maximize_c_profile_interval(bracket_lo,bracket_hi,n,k,
        p_A,p_B,p_c,weights,apply_prior,prior_mean,prior_var,refined_c,refined_r);
    if (!std::isfinite(refined_ll) || !std::isfinite(refined_c) ||
        !std::isfinite(refined_r)){
        out.profile_optimization_succeeded=false;
        out.status="nonfinite_two_dimensional_profile_refinement";
        return out;
    }
    out.validated_r=clamp_closed_unit_interval(refined_r);
    out.validated_c=clamp_closed_unit_interval(refined_c);
    out.validated_objective=refined_ll;
    if (out.raw_candidate_available && std::isfinite(profile_at_raw)){
        out.profile_global_gain=refined_ll-profile_at_raw;
    }

    if (out.validated_c <= C_BOUNDARY_TOLERANCE ||
        out.validated_c >= 1.0-C_BOUNDARY_TOLERANCE){
        out.c_boundary_candidate=true;
        out.boundary=assess_unpenalized_heterotypic_c_boundary(
            out.validated_r,out.validated_c,n,k,p_A,p_B,p_c,weights,
            apply_prior);
        out.profile_neighbor_gain=out.boundary.ll_best_interior-
            out.boundary.ll_boundary;
        out.profile_support_low=out.boundary.support_low;
        out.profile_support_high=out.boundary.support_high;
        out.profile_support_width=out.boundary.support_width;
        if (std::isfinite(out.boundary.nuisance_r_at_boundary) &&
            std::isfinite(out.boundary.ll_boundary)){
            out.validated_r=clamp_closed_unit_interval(
                out.boundary.nuisance_r_at_boundary);
            out.validated_c=out.boundary.boundary<0 ? 0.0 : 1.0;
            out.validated_objective=out.boundary.ll_boundary;
        }
        out.accepted=out.boundary.accepted;
        out.status=out.boundary.status;
        if (out.accepted){
            assess_r_profile_identifiability(out,n,k,p_A,p_B,p_c,weights,
                apply_prior,prior_mean,prior_var);
        }
        return out;
    }

    ContinuousProfileIntervalAssessment c_support=continuous_profile_interval(
        true,out.validated_r,out.validated_c,out.validated_objective,n,k,p_A,
        p_B,p_c,weights,apply_prior,prior_mean,prior_var);
    if (!c_support.success){
        out.status=c_support.status;
        return out;
    }
    out.profile_support_low=c_support.lower;
    out.profile_support_high=c_support.upper;
    out.profile_support_width=c_support.width;

    const double h=std::min(0.01,
        0.45*std::min(out.validated_c,1.0-out.validated_c));
    double nuisance_left=0.5,nuisance_right=0.5;
    const double ll_left=maximize_profile_coordinate(true,out.validated_c-h,
        n,k,p_A,p_B,p_c,weights,apply_prior,prior_mean,prior_var,nuisance_left);
    const double ll_right=maximize_profile_coordinate(true,out.validated_c+h,
        n,k,p_A,p_B,p_c,weights,apply_prior,prior_mean,prior_var,nuisance_right);
    out.profile_neighbor_gain=std::max(ll_left,ll_right)-refined_ll;
    if (!std::isfinite(out.profile_neighbor_gain) ||
        out.profile_neighbor_gain > MAX_TWO_DIMENSIONAL_PROFILE_NEIGHBOR_GAIN){
        out.status="two_dimensional_profile_refinement_not_locally_maximal";
        return out;
    }

    double grad_r=0.0,grad_c=0.0;
    eval_three_dataset_gradient(out.validated_r,out.validated_c,n,k,p_A,p_B,
        p_c,weights,apply_prior,prior_mean,prior_var,grad_r,grad_c);
    const double vr=projected_gradient_violation(out.validated_r,grad_r);
    const double vc=projected_gradient_violation(out.validated_c,grad_c);
    out.projected_gradient_norm=sqrt(vr*vr+vc*vc);
    if (!std::isfinite(out.projected_gradient_norm)){
        out.status="nonfinite_two_dimensional_projected_gradient";
        return out;
    }

    assess_r_profile_identifiability(out,n,k,p_A,p_B,p_c,weights,
        apply_prior,prior_mean,prior_var);
    out.accepted=true;
    if (!out.raw_candidate_available){
        out.status="two_dimensional_profile_fit_without_bfgs";
    } else {
        const bool materially_refined =
            fabs(out.validated_r-bounded_raw_r)>1e-5 ||
            fabs(out.validated_c-bounded_raw_c)>1e-5;
        out.status=materially_refined ?
            "two_dimensional_same_model_profile_refined" :
            "two_dimensional_converged_profile_validated";
    }
    return out;
}

TwoDimensionalFitAssessment assess_two_dimensional_fit_candidate(
    double candidate_r, double candidate_c,
    const vector<double>& n, const vector<double>& k,
    const vector<double>& p_A, const vector<double>& p_B,
    const vector<double>& p_c, const vector<double>& weights,
    bool apply_prior, double prior_mean, double prior_var){
    vector<pair<double,double> > candidates;
    candidates.push_back(make_pair(candidate_r,candidate_c));
    return fit_two_dimensional_profile_model(candidates,n,k,p_A,p_B,p_c,
        weights,apply_prior,prior_mean,prior_var);
}


static bool numerical_three_hessian(double r, double c,
    const vector<double>& n, const vector<double>& k,
    const vector<double>& p_A, const vector<double>& p_B,
    const vector<double>& p_c, const vector<double>& obs_weights,
    bool apply_prior, double prior_mean, double prior_var,
    double& se_r, double& se_c, double& corr){

    const double h = 1e-4;
    if (r <= 2*h || r >= 1.0-2*h || c <= 2*h || c >= 1.0-2*h){
        return false;
    }
    auto f = [&](double rr, double cc){
        return eval_three_dataset_objective(rr, cc, n, k, p_A, p_B, p_c,
            obs_weights, apply_prior, prior_mean, prior_var);
    };
    const double f00 = f(r,c);
    const double frr = (f(r+h,c) - 2.0*f00 + f(r-h,c))/(h*h);
    const double fcc = (f(r,c+h) - 2.0*f00 + f(r,c-h))/(h*h);
    const double frc = (f(r+h,c+h)-f(r+h,c-h)-f(r-h,c+h)+f(r-h,c-h))/(4.0*h*h);
    const double i11 = -frr;
    const double i22 = -fcc;
    const double i12 = -frc;
    const double det = i11*i22 - i12*i12;
    if (!std::isfinite(det) || det <= 1e-12 || i11 <= 0.0 || i22 <= 0.0){
        return false;
    }
    const double v_r = i22/det;
    const double v_c = i11/det;
    const double cov = -i12/det;
    if (v_r <= 0.0 || v_c <= 0.0){
        return false;
    }
    se_r = sqrt(v_r);
    se_c = sqrt(v_c);
    corr = cov / sqrt(v_r*v_c);
    if (corr < -1.0) corr = -1.0;
    if (corr > 1.0) corr = 1.0;
    return std::isfinite(se_r) && std::isfinite(se_c) && std::isfinite(corr);
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
    set<int>& allowed_ids2)
    : contamFinder3(indv_allelecounts, assn, assn_llr, exp_match_fracs,
        n_samples, allowed_ids, allowed_ids, allowed_ids2){}

contamFinder3::contamFinder3(robin_hood::unordered_map<unsigned long,
        map<pair<int, int>, map<pair<int, int>, pair<float, float> > > >& indv_allelecounts,
    robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    map<pair<int, int>, map<int, float> >& exp_match_fracs,
    int n_samples,
    set<int>& profile_allowed_ids,
    set<int>& reassign_allowed_ids,
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
    this->profile_total_starts = -1;
    this->expfracs = exp_match_fracs;
    this->allowed_ids = profile_allowed_ids;
    this->reassign_allowed_ids = reassign_allowed_ids;
    this->allowed_ids2 = allowed_ids2;
    
    // Reassign cells by default
    this->skip_reassign = false;

    // Don't worry about doublet rate by default
    this->doublet_rate = -1;
    
    this->ef_all_avg = true;
    this->final_loglik = 0.0;
    this->final_loglik_valid = false;

    llrtot = 0.0;

    bool caller_restricted_ids = (profile_allowed_ids.size() > 0);

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
    this->heterotypic_start_mode = START_TOPK;
    this->heterotypic_start_top_k = 5;
    this->adaptive_profile_intervals = true;

    this->contam_prior_mode = PRIOR_FUSION;
    this->fix_r_enabled = false;
    this->fixed_r_by_ident.clear();
    this->contam_prior_min_cells = 20;
    this->contam_prior_max_ci_width = 0.50;
    this->contam_prior_min_informative_weight = 10.0;
    this->heterotypic_contam_prior = make_pair(-1.0, -1.0);
    this->twocomp_contam_prior = make_pair(-1.0, -1.0);
    this->heterotypic_contam_prior_valid = false;
    this->twocomp_contam_prior_valid = false;
    this->dispersion_phi = 1.0;
    this->profile_successful_starts = 0;
    this->profile_near_optimal_count = 0;
    this->profile_best_ll = -DBL_MAX;
    this->profile_second_best_ll = -DBL_MAX;
    this->profile_near_optimal_l1_spread = std::numeric_limits<double>::quiet_NaN();
    this->multistart_attempted = false;
    this->multistart_configured_starts = 0;
    this->multistart_successful_starts = 0;
    this->multistart_near_optimal_count = 0;
    this->multistart_best_ll = -DBL_MAX;
    this->multistart_second_best_ll = -DBL_MAX;
    this->multistart_near_optimal_l1_spread = std::numeric_limits<double>::quiet_NaN();

    // Tetraploid-aware defaults (v1.4)
    this->tetraploid_aware = false;
    this->min_signal_gap = 0.10;
    this->amb_mu_available = false;
    this->ids_restricted = false;

    // O53 step 2 defaults: historical behaviour unless explicitly enabled.
    // The split default of 0.05 is measured, not chosen: swept over seven
    // tetraploid units at 0.0, 0.05, 0.1, 0.2, 0.4 and 0.6, it minimises mean
    // profile L1 against the observed-fraction truth at -0.227, improving six
    // of seven units. Both endpoints are worse, and 0.0 is much worse, which is
    // why the assignment-keyed block is kept rather than replaced.
    this->candidate_keyed_rows = false;
    this->candidate_keyed_split = 0.05;

    // Source-exclusion default: exact historical global-profile behavior.
    this->use_loo = false;
    this->source_exclusion_strength = 0.0;
    this->profile_holdout_barcodes.clear();
    this->profile_holdout_basis_label = "none";
    this->fixed_r_basis_label = "none";
    this->fixed_ambient_basis_label = "none";
    this->profile_variant_label = "fitted_global";
    this->truth_assisted_condition = false;

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

void contamFinder3::set_profile_total_starts(int nt){
    this->profile_total_starts = nt;
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
    if (!enabled && heterotypic_start_mode == START_TOPK){
        heterotypic_start_mode = START_SINGLE;
    }
}

void contamFinder3::set_heterotypic_start_mode(const string& mode, int top_k){
    if (mode == "single") heterotypic_start_mode = START_SINGLE;
    else if (mode == "topk") heterotypic_start_mode = START_TOPK;
    else if (mode == "exhaustive") heterotypic_start_mode = START_EXHAUSTIVE;
    else {
        fprintf(stderr, "ERROR: unknown heterotypic start mode: %s\n", mode.c_str());
        exit(1);
    }
    heterotypic_start_top_k = std::max(1, std::min(25, top_k));
}

void contamFinder3::set_adaptive_profile_intervals(bool enabled){
    adaptive_profile_intervals = enabled;
}

void contamFinder3::set_contam_prior_mode(const string& mode){
    if (mode == "none") contam_prior_mode = PRIOR_NONE;
    else if (mode == "global") contam_prior_mode = PRIOR_GLOBAL;
    else if (mode == "heterotypic") contam_prior_mode = PRIOR_HETEROTYPIC;
    else if (mode == "fusion") contam_prior_mode = PRIOR_FUSION;
    else {
        fprintf(stderr, "ERROR: unknown contamination prior mode: %s\n", mode.c_str());
        exit(1);
    }
}

void contamFinder3::set_contam_prior_support(int min_cells, double max_ci_width,
    double min_informative_weight, double deprecated_max_gradient_norm){
    contam_prior_min_cells = std::max(2, min_cells);
    contam_prior_max_ci_width = max_ci_width;
    contam_prior_min_informative_weight = min_informative_weight;
    // Retained only for source/API compatibility with older callers. Raw
    // objective gradients scale with read depth and are not scientifically
    // valid prior-training eligibility gates.
    (void)deprecated_max_gradient_norm;
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

void contamFinder3::set_candidate_keyed_rows(bool enabled, double split){
    this->candidate_keyed_rows = enabled;
    if (split < 0.0){
        split = 0.0;
    }
    if (split > 1.0){
        split = 1.0;
    }
    this->candidate_keyed_split = split;
}

void contamFinder3::set_use_loo(bool enabled){
    this->use_loo = enabled;
    this->source_exclusion_strength = enabled ? 1.0 : 0.0;
}

void contamFinder3::set_source_exclusion_strength(double strength){
    if (!std::isfinite(strength) || strength < 0.0 || strength > 1.0){
        throw std::invalid_argument("source exclusion strength must be finite and in [0,1]");
    }
    this->source_exclusion_strength = strength;
    this->use_loo = strength > 0.0;
}

void contamFinder3::set_profile_holdout_barcodes(
    const set<unsigned long>& held_out, const string& basis_label){
    this->profile_holdout_barcodes = held_out;
    this->profile_holdout_basis_label = basis_label.empty() ?
        "library_crossfit" : basis_label;
}

void contamFinder3::set_diagnostic_context(const string& fixed_r_basis,
    const string& fixed_ambient_basis, bool truth_assisted){
    this->fixed_r_basis_label = fixed_r_basis.empty() ? "none" : fixed_r_basis;
    this->fixed_ambient_basis_label = fixed_ambient_basis.empty() ? "none" : fixed_ambient_basis;
    this->truth_assisted_condition = truth_assisted;
    const bool oracle = this->fixed_amb_prof || this->fixed_ambient_basis_label == "truth";
    if (source_exclusion_strength >= 1.0 - 1e-12){
        this->profile_variant_label = oracle ? "oracle_hard_source_excluded" :
            "fitted_hard_source_excluded";
    } else if (source_exclusion_strength <= 1e-12){
        this->profile_variant_label = oracle ? "oracle_global" : "fitted_global";
    } else {
        this->profile_variant_label = oracle ? "oracle_partial_source_excluded" :
            "fitted_partial_source_excluded";
    }
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

void contamFinder3::set_fixed_r(const map<int, double>& fr){
    this->fixed_r_by_ident = fr;
    this->fix_r_enabled = !fr.empty();
    if (this->fix_r_enabled){
        fprintf(stderr, "Fixed allele ratio supplied for %lu identity/identities. "
            "r will NOT be fitted for those cells; only c is solved.\n",
            this->fixed_r_by_ident.size());
    }
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


bool contamination_prior_training_supported(
    bool fit_success, double log_likelihood, double informative_weight,
    double min_informative_weight, bool is_heterotypic,
    bool profile_optimization_attempted, bool profile_optimization_succeeded,
    bool profile_validation_passed, double profile_support_width,
    double max_profile_width, bool c_boundary, string& reason){

    if (!fit_success || !std::isfinite(log_likelihood)){
        reason="fit_failed";
        return false;
    }
    if (!std::isfinite(informative_weight) ||
        informative_weight < min_informative_weight){
        reason="low_informative_weight";
        return false;
    }
    if (is_heterotypic){
        if (!profile_optimization_attempted ||
            !profile_optimization_succeeded){
            reason="profile_optimization_unresolved";
            return false;
        }
        if (!profile_validation_passed ||
            !std::isfinite(profile_support_width) ||
            profile_support_width < 0.0 || profile_support_width > 1.0){
            reason="profile_interval_unresolved";
            return false;
        }
        if (!(max_profile_width>0.0) ||
            profile_support_width > max_profile_width){
            reason="wide_contam_profile_interval";
            return false;
        }
        if (c_boundary && profile_support_width >
            std::min(0.20,max_profile_width)){
            reason="boundary_profile_unresolved";
            return false;
        }
    }
    reason="eligible";
    return true;
}

bool contamFinder3::prior_training_supported(const CellContamFitResult& fit,
    string& reason) const{
    return contamination_prior_training_supported(
        fit.fit_success,fit.ll,fit.informative_weight,
        contam_prior_min_informative_weight,fit.is_heterotypic,
        fit.profile_optimization_attempted!=0,
        fit.profile_optimization_succeeded!=0,
        fit.profile_validation_passed!=0,
        fit.c_profile_support_width,contam_prior_max_ci_width,
        fit.hit_c_lower || fit.hit_c_upper,reason);
}

std::pair<double,double> contamination_prior_safe_moments(
    const vector<double>& vals){
    if (vals.size() < 2) return make_pair(-1.0, -1.0);
    vector<double> vals_copy = vals;
    pair<double,double> mv = welford(vals_copy);
    double mu = std::min(1.0-1e-6, std::max(1e-6, mv.first));
    double vmax = mu*(1.0-mu);
    double var = mv.second;
    if (!std::isfinite(var) || var <= 1e-8) var = std::min(1e-4, 0.25*vmax);
    if (var >= vmax) var = 0.95*vmax;
    return make_pair(mu,var);
}


void contamFinder3::learn_contam_priors(const vector<CellContamFitResult>& mle_results){
    fusion_contam_priors.clear();
    heterotypic_contam_prior_valid = false;
    twocomp_contam_prior_valid = false;

    vector<double> hetero_vals;
    vector<double> twocomp_vals;
    map<int, vector<double> > fusion_vals;

    for (vector<CellContamFitResult>::const_iterator it = mle_results.begin();
        it != mle_results.end(); ++it){
        string reason;
        bool ok = prior_training_supported(*it, reason);
        prior_training_eligible[it->barcode] = ok ? 1 : 0;
        prior_training_reason_by_cell[it->barcode] = reason;
        if (!ok) continue;
        if (it->is_heterotypic){
            hetero_vals.push_back(it->c);
            if (assn.count(it->barcode) > 0){
                fusion_vals[assn.at(it->barcode)].push_back(it->c);
            }
        } else {
            twocomp_vals.push_back(it->c);
        }
    }

    if ((int)hetero_vals.size() >= contam_prior_min_cells){
        heterotypic_contam_prior = contamination_prior_safe_moments(hetero_vals);
        heterotypic_contam_prior_valid = heterotypic_contam_prior.second > 0.0;
    }
    if ((int)twocomp_vals.size() >= contam_prior_min_cells){
        twocomp_contam_prior = contamination_prior_safe_moments(twocomp_vals);
        twocomp_contam_prior_valid = twocomp_contam_prior.second > 0.0;
    }
    if (contam_prior_mode == PRIOR_FUSION){
        for (map<int, vector<double> >::const_iterator fv = fusion_vals.begin();
            fv != fusion_vals.end(); ++fv){
            if ((int)fv->second.size() >= contam_prior_min_cells){
                pair<double,double> mv = contamination_prior_safe_moments(fv->second);
                if (mv.second > 0.0) fusion_contam_priors[fv->first] = mv;
            }
        }
    }

    fprintf(stderr, "Empirical contamination prior training: %lu supported heterotypic, "
        "%lu supported two-component, %lu fusion priors\n",
        hetero_vals.size(), twocomp_vals.size(), fusion_contam_priors.size());
}

pair<double,double> contamFinder3::prior_for_cell(unsigned long barcode,
    bool is_heterotypic, string& group) const{
    group = "none";
    if (user_prior_set){
        group = "user_fixed";
        return make_pair(user_prior_mean, user_prior_var);
    }
    if (contam_prior_mode == PRIOR_NONE){
        return make_pair(-1.0,-1.0);
    }
    if (is_heterotypic){
        if (contam_prior_mode == PRIOR_FUSION && assn.count(barcode) > 0){
            map<int,pair<double,double> >::const_iterator fp =
                fusion_contam_priors.find(assn.at(barcode));
            if (fp != fusion_contam_priors.end()){
                group = "fusion";
                return fp->second;
            }
        }
        if (heterotypic_contam_prior_valid){
            group = (contam_prior_mode == PRIOR_GLOBAL) ? "global_heterotypic" : "heterotypic";
            return heterotypic_contam_prior;
        }
        return make_pair(-1.0,-1.0);
    }
    if (twocomp_contam_prior_valid){
        group = "two_component";
        return twocomp_contam_prior;
    }
    return make_pair(-1.0,-1.0);
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
                    double ef_val = 0.5;
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
                                double ef1 = 0.5;
                                double ef2 = 0.5;
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
// Conditional-fraction coverage audit
// =============================================================================

contamFinder3::CondfCoverageReport contamFinder3::write_condf_coverage_report(
    const std::string& filename,
    const std::vector<std::string>& samples) const{

    typedef std::tuple<int, int, int> LookupKey;
    std::map<LookupKey, std::pair<unsigned long long, double> > required;

    for (auto cell = indv_allelecounts.begin(); cell != indv_allelecounts.end(); ++cell){
        for (auto first = cell->second.begin(); first != cell->second.end(); ++first){
            for (auto second = first->second.begin(); second != first->second.end(); ++second){
                double row_weight = (double)second->second.first + (double)second->second.second;
                if (!(row_weight > 0.0)) continue;

                const std::pair<int,int>& key1 = first->first;
                const std::pair<int,int>& key2 = second->first;
                const bool is_doublet = key2.first >= 0;
                for (int source : idx2samp){
                    // In ef_all_avg mode, a receiver constituent uses its observed
                    // dosage directly.  For a doublet, either constituent bypasses
                    // BOTH conditional lookups in the likelihood, so do not create a
                    // false strict-condf requirement for the other parent row.
                    if (ef_all_avg &&
                        (source == key1.first || (is_doublet && source == key2.first))){
                        continue;
                    }
                    LookupKey lk1 = std::make_tuple(key1.first, key1.second, source);
                    required[lk1].first += 1;
                    required[lk1].second += row_weight;
                    if (is_doublet){
                        LookupKey lk2 = std::make_tuple(key2.first, key2.second, source);
                        required[lk2].first += 1;
                        required[lk2].second += row_weight;
                    }
                }
            }
        }
    }

    CondfCoverageReport report;
    report.required_lookups = 0;
    report.missing_lookups = 0;
    report.required_weight = 0.0;
    report.missing_weight = 0.0;

    FILE* out = fopen(filename.c_str(), "w");
    if (out == NULL){
        fprintf(stderr, "ERROR: could not write conditional-fraction coverage report: %s\n",
            filename.c_str());
        report.missing_lookups = 1;
        report.missing_weight = 1.0;
        return report;
    }
    fprintf(out, "condition_sample_index\tcondition_sample\tnalt\tsource_index\tsource_sample\tlookup_rows\tlookup_weight\tstatus\n");

    for (auto it = required.begin(); it != required.end(); ++it){
        int condition_idx = std::get<0>(it->first);
        int nalt = std::get<1>(it->first);
        int source_idx = std::get<2>(it->first);
        unsigned long long rows = it->second.first;
        double weight = it->second.second;
        bool present = false;
        auto row = expfracs.find(std::make_pair(condition_idx, nalt));
        if (row != expfracs.end() && row->second.count(source_idx) > 0){
            present = true;
        }
        report.required_lookups++;
        report.required_weight += weight;
        if (!present){
            report.missing_lookups++;
            report.missing_weight += weight;
        }
        const char* condition_name = (condition_idx >= 0 && condition_idx < (int)samples.size())
            ? samples[condition_idx].c_str() : "NA";
        const char* source_name = (source_idx >= 0 && source_idx < (int)samples.size())
            ? samples[source_idx].c_str() : "NA";
        fprintf(out, "%d\t%s\t%d\t%d\t%s\t%llu\t%.17g\t%s\n",
            condition_idx, condition_name, nalt, source_idx, source_name,
            rows, weight, present ? "present" : "missing");
    }
    fclose(out);

    fprintf(stderr,
        "Conditional-fraction coverage: required=%llu missing=%llu required_weight=%.6g missing_weight=%.6g report=%s\n",
        report.required_lookups, report.missing_lookups,
        report.required_weight, report.missing_weight, filename.c_str());
    return report;
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
    // Empty-matrix guard. A zero-row component matrix makes add_mixcomp return
    // false, leaving the solver with nmixcomp==0 and an empty
    // mixcompfracs_sparse; calling solve() then indexes that empty structure and
    // segfaults. Skip the solve cleanly and keep the current profile. Returning
    // 0.0 matches this function's other failure paths (no-data, no-active-species,
    // all-solves-failed), so downstream consumers already tolerate it.
    if (species_mixfracs.empty()){
        fprintf(stderr, "WARNING: empty species mixture component matrix; "
            "skipping species-level pi solve (keeping current profile)\n");
        return 0.0;
    }
    if (!solver.add_mixcomp(species_mixfracs) || !solver.add_mixcomp_fracs(startprops)){
        fprintf(stderr, "WARNING: species mixture component setup failed "
            "(empty or mismatched matrix); skipping species-level pi solve "
            "(keeping current profile)\n");
        return 0.0;
    }
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

        // Use the same fixed-prior path as an explicit user prior so the
        // per-cell prior dispatcher cannot accidentally ignore the adaptive
        // override.
        set_user_prior(adaptive_prior_mean, current_var);

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

bool contamination_observation_passes_model_filter(
    bool filter_required, bool ambient_profile_available,
    double p_e, double p_c, double min_signal_gap,
    int nalt1, int nalt2){
    if (!filter_required) return true;
    if (ambient_profile_available){
        return category_passes_adaptive_filter(p_e,p_c,min_signal_gap);
    }
    return category_passes_hard_filter(true,nalt1,nalt2);
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

bool weighted_composition_expected_from_row_complete(
    const map<int, double>& comp,
    const pair<int, int>& key1,
    const pair<int, int>& key2,
    const map<pair<int, int>, map<int, float> >& conditional_fractions,
    double& expected){

    double out = 0.0;
    double wsum = 0.0;
    for (auto it = comp.begin(); it != comp.end(); ++it){
        const int idx = it->first;
        const double w = it->second;
        if (!(w > 0.0) || !std::isfinite(w)) continue;

        double value = 0.0;
        bool have = false;
        if (key1.first == idx){
            value = (double)key1.second/2.0;
            have = true;
        } else if (key2.first == idx){
            value = (double)key2.second/2.0;
            have = true;
        } else {
            double acc = 0.0;
            int nhave = 0;
            auto e1 = conditional_fractions.find(key1);
            if (e1 != conditional_fractions.end()){
                auto v = e1->second.find(idx);
                if (v != e1->second.end() && std::isfinite(v->second)){
                    acc += v->second;
                    nhave++;
                }
            }
            if (key2.first != -1){
                auto e2 = conditional_fractions.find(key2);
                if (e2 != conditional_fractions.end()){
                    auto v = e2->second.find(idx);
                    if (v != e2->second.end() && std::isfinite(v->second)){
                        acc += v->second;
                        nhave++;
                    }
                }
            }
            if (nhave > 0){
                value = acc/(double)nhave;
                have = true;
            }
        }

        // A neutral substitute would conceal missing model support and can
        // manufacture a weak or flat contamination likelihood.
        if (!have) return false;
        out += w*value;
        wsum += w;
    }

    if (!(wsum > 0.0) || !std::isfinite(out)) return false;
    expected = out/wsum;
    return std::isfinite(expected);
}

bool contamFinder3::composition_expected_from_row_complete(
    const map<int, double>& comp,
    const pair<int, int>& key1,
    const pair<int, int>& key2, double& expected) const{
    return weighted_composition_expected_from_row_complete(
        comp,key1,key2,expfracs,expected);
}

set<int> weighted_composition_source_exclusions(
    const map<int, double>& composition){
    set<int> excluded;
    for (auto it = composition.begin(); it != composition.end(); ++it){
        if (it->second > 0.0 && std::isfinite(it->second)) excluded.insert(it->first);
    }
    return excluded;
}

SourceExclusionProfileAssessment::SourceExclusionProfileAssessment():
    source_exclusion_strength(0.0), global_profile_mass_sum(0.0),
    scoring_profile_mass_sum(0.0), excluded_source_mass_raw_total(0.0),
    effective_removed_source_mass_total(0.0),
    scoring_profile_renormalization_denominator(0.0),
    status("not_assessed"){}

SourceExclusionProfileAssessment apply_source_exclusion_profile(
    const map<int, double>& global_profile, const set<int>& excluded_sources,
    double strength){

    SourceExclusionProfileAssessment out;
    out.source_exclusion_strength = strength;
    double eligible_nonnegative_mass = 0.0;
    for (auto it = global_profile.begin(); it != global_profile.end(); ++it){
        if (std::isfinite(it->second) && it->second > 0.0){
            out.global_profile_mass_sum += it->second;
            if (it->first >= 0){
                eligible_nonnegative_mass += it->second;
                if (excluded_sources.count(it->first) > 0){
                    out.excluded_source_mass_raw_total += it->second;
                }
            }
        }
    }
    if (!std::isfinite(strength) || strength < 0.0 || strength > 1.0){
        out.scoring_profile = global_profile;
        out.scoring_profile_mass_sum = out.global_profile_mass_sum;
        out.scoring_profile_renormalization_denominator = out.global_profile_mass_sum;
        out.status = "invalid_source_exclusion_strength";
        return out;
    }

    // Exact no-exclusion endpoint: preserve every key and floating-point value.
    if (strength == 0.0 || excluded_sources.empty()){
        out.scoring_profile = global_profile;
        out.scoring_profile_mass_sum = out.global_profile_mass_sum;
        out.scoring_profile_renormalization_denominator = out.global_profile_mass_sum;
        out.effective_removed_source_mass_total = 0.0;
        out.status = excluded_sources.empty() ? "no_eligible_sources" : "global_profile_no_exclusion";
        return out;
    }

    out.effective_removed_source_mass_total =
        strength * out.excluded_source_mass_raw_total;
    if (strength == 1.0){
        // Sum retained sources directly in the same map iteration order as the
        // historical hard-exclusion path. This preserves that endpoint's
        // floating-point denominator rather than subtracting two rounded totals.
        out.scoring_profile_renormalization_denominator = 0.0;
        for (auto it = global_profile.begin(); it != global_profile.end(); ++it){
            if (it->first >= 0 && excluded_sources.count(it->first) == 0 &&
                std::isfinite(it->second) && it->second > 0.0){
                out.scoring_profile_renormalization_denominator += it->second;
            }
        }
    } else {
        out.scoring_profile_renormalization_denominator =
            eligible_nonnegative_mass - out.effective_removed_source_mass_total;
    }

    // Preserve the historical fallback; do not silently substitute a uniform
    // profile. Since the scoring profile remains global, actual removed mass is 0.
    if (!(out.scoring_profile_renormalization_denominator > 1e-12) ||
        !std::isfinite(out.scoring_profile_renormalization_denominator)){
        out.scoring_profile = global_profile;
        out.scoring_profile_mass_sum = out.global_profile_mass_sum;
        out.effective_removed_source_mass_total = 0.0;
        out.scoring_profile_renormalization_denominator = out.global_profile_mass_sum;
        out.status = "degenerate_all_mass_excluded_fallback_global";
        return out;
    }

    const double inv = 1.0/out.scoring_profile_renormalization_denominator;
    for (auto it = global_profile.begin(); it != global_profile.end(); ++it){
        // The legacy hard-exclusion scoring simplex contains concrete source
        // indices only. Use the same source universe for intermediate strengths.
        if (it->first < 0) continue;
        double value = it->second;
        if (!std::isfinite(value) || value < 0.0) value = 0.0;
        if (excluded_sources.count(it->first) > 0){
            value *= (1.0-strength);
        }
        value *= inv;
        out.scoring_profile[it->first] = value;
        out.scoring_profile_mass_sum += value;
    }
    out.status = strength == 1.0 ? "hard_source_exclusion" : "partial_source_exclusion";
    return out;
}

ParentAxisGeometryAssessment::ParentAxisGeometryAssessment():
    alpha(std::numeric_limits<double>::quiet_NaN()),
    orthogonal_norm(std::numeric_limits<double>::quiet_NaN()),
    denominator(std::numeric_limits<double>::quiet_NaN()),
    status("not_assessed"){}

ParentAxisGeometryAssessment assess_parent_axis_geometry(
    const vector<double>& p_A, const vector<double>& p_B,
    const vector<double>& p_ambient, const vector<double>& weights){

    ParentAxisGeometryAssessment out;
    if (p_A.size() != p_B.size() || p_A.size() != p_ambient.size() ||
        (!weights.empty() && weights.size() != p_A.size()) || p_A.empty()){
        out.status = "invalid_geometry_inputs";
        return out;
    }
    double numerator = 0.0;
    double denominator = 0.0;
    double wsum = 0.0;
    for (size_t i = 0; i < p_A.size(); ++i){
        const double w = weights.empty() ? 1.0 : weights[i];
        if (!(w > 0.0) || !std::isfinite(w) || !std::isfinite(p_A[i]) ||
            !std::isfinite(p_B[i]) || !std::isfinite(p_ambient[i])) continue;
        const double axis = p_A[i]-p_B[i];
        numerator += w*(p_ambient[i]-p_B[i])*axis;
        denominator += w*axis*axis;
        wsum += w;
    }
    out.denominator = denominator;
    if (!(wsum > 0.0)){
        out.status = "no_positive_geometry_weight";
        return out;
    }
    if (!(denominator > 1e-14) || !std::isfinite(denominator)){
        out.status = "zero_parent_axis_denominator";
        return out;
    }
    out.alpha = numerator/denominator;
    double orth = 0.0;
    for (size_t i = 0; i < p_A.size(); ++i){
        const double w = weights.empty() ? 1.0 : weights[i];
        if (!(w > 0.0) || !std::isfinite(w) || !std::isfinite(p_A[i]) ||
            !std::isfinite(p_B[i]) || !std::isfinite(p_ambient[i])) continue;
        const double projected = p_B[i] + out.alpha*(p_A[i]-p_B[i]);
        const double residual = p_ambient[i]-projected;
        orth += w*residual*residual;
    }
    out.orthogonal_norm = sqrt(std::max(0.0, orth/wsum));
    out.status = "ok";
    return out;
}

ConditionalInformationAssessment::ConditionalInformationAssessment():
    information_rr(std::numeric_limits<double>::quiet_NaN()),
    information_cc(std::numeric_limits<double>::quiet_NaN()),
    information_rc(std::numeric_limits<double>::quiet_NaN()),
    conditional_information_c_given_r(std::numeric_limits<double>::quiet_NaN()),
    status("not_assessed"){}

ConditionalInformationAssessment assess_conditional_information_rc(
    double r, double c, const vector<double>& n, const vector<double>& k,
    const vector<double>& p_A, const vector<double>& p_B,
    const vector<double>& p_ambient, const vector<double>& weights){

    ConditionalInformationAssessment out;
    if (!std::isfinite(r) || !std::isfinite(c) || n.size() != k.size() ||
        n.size() != p_A.size() || n.size() != p_B.size() ||
        n.size() != p_ambient.size() || (!weights.empty() && weights.size() != n.size()) ||
        n.empty()){
        out.status = "invalid_information_inputs";
        return out;
    }
    double irr = 0.0, icc = 0.0, irc = 0.0;
    int used = 0;
    for (size_t i = 0; i < n.size(); ++i){
        const double w = weights.empty() ? 1.0 : weights[i];
        if (!(w > 0.0) || !(n[i] > 0.0) || !std::isfinite(w) ||
            !std::isfinite(k[i]) || !std::isfinite(p_A[i]) ||
            !std::isfinite(p_B[i]) || !std::isfinite(p_ambient[i])) continue;
        const double axis = p_A[i]-p_B[i];
        const double endogenous = p_B[i] + r*axis;
        double prob = (1.0-c)*endogenous + c*p_ambient[i];
        prob = clamp_unit_interval(prob, 1e-9);
        const double dldr_p = (k[i]-n[i]*prob)/(prob*(1.0-prob));
        const double d2ldp2 = -k[i]/(prob*prob) -
            (n[i]-k[i])/((1.0-prob)*(1.0-prob));
        const double dpdr = (1.0-c)*axis;
        const double dpdc = p_ambient[i]-endogenous;
        irr += -w*d2ldp2*dpdr*dpdr;
        icc += -w*d2ldp2*dpdc*dpdc;
        irc += -w*(d2ldp2*dpdr*dpdc - dldr_p*axis);
        used++;
    }
    out.information_rr = irr;
    out.information_cc = icc;
    out.information_rc = irc;
    if (used == 0){
        out.status = "no_finite_information_rows";
        return out;
    }
    if (!(irr > 1e-12) || !std::isfinite(irr)){
        out.status = "information_rr_nonpositive_or_singular";
        return out;
    }
    const double conditional = icc - irc*irc/irr;
    if (!std::isfinite(conditional) || conditional < -1e-8){
        out.status = "information_schur_complement_invalid";
        return out;
    }
    out.conditional_information_c_given_r = std::max(0.0, conditional);
    out.status = "ok";
    return out;
}

vector<pair<int, int> > rank_weighted_composition_fallback_bases(
    const map<int, double>& comp,
    const map<pair<int, int>, map<pair<int, int>, pair<float, float> > >& allelecounts){

    struct RankedBasis {
        pair<int,int> key;
        int direct_components;
        double support;
    };
    map<pair<int,int>, RankedBasis> support_by_basis;

    for (auto ac = allelecounts.begin(); ac != allelecounts.end(); ++ac){
        for (auto ac2 = ac->second.begin(); ac2 != ac->second.end(); ++ac2){
            const double support = (double)ac2->second.first +
                (double)ac2->second.second;
            if (!(support > 0.0) || !std::isfinite(support)) continue;

            const int i = ac->first.first;
            const int j = ac2->first.first;
            pair<int,int> key = make_pair(i,j);
            int direct = comp.count(i) > 0 ? 1 : 0;
            if (j >= 0 && j != i && comp.count(j) > 0) direct++;
            if (direct <= 0) continue;

            auto found = support_by_basis.find(key);
            if (found == support_by_basis.end()){
                RankedBasis rb = {key,direct,support};
                support_by_basis[key] = rb;
            } else {
                found->second.direct_components =
                    std::max(found->second.direct_components,direct);
                found->second.support += support;
            }
        }
    }

    vector<RankedBasis> ranked_bases;
    for (auto it = support_by_basis.begin(); it != support_by_basis.end(); ++it){
        ranked_bases.push_back(it->second);
    }
    std::sort(ranked_bases.begin(),ranked_bases.end(),
        [](const RankedBasis& a,const RankedBasis& b){
            // Direct biological coverage is primary. Depth breaks ties only
            // among bases representing the same number of native components.
            if (a.direct_components != b.direct_components){
                return a.direct_components > b.direct_components;
            }
            if (a.support != b.support) return a.support > b.support;
            return a.key < b.key;
        });

    vector<pair<int,int> > ranked;
    for (auto it = ranked_bases.begin(); it != ranked_bases.end(); ++it){
        ranked.push_back(it->key);
    }
    return ranked;
}

pair<int, int> select_weighted_composition_fallback_basis(
    const map<int, double>& comp,
    const map<pair<int, int>, map<pair<int, int>, pair<float, float> > >& allelecounts){
    vector<pair<int,int> > ranked =
        rank_weighted_composition_fallback_bases(comp,allelecounts);
    return ranked.empty() ? make_pair(-1,-1) : ranked.front();
}

WeightedCompositionRowSelection::WeightedCompositionRowSelection():
    direct_component_rows(false),fallback_basis(make_pair(-1,-1)),
    status("weighted_composition_no_complete_conditional_rows"){}

WeightedCompositionRowSelection select_weighted_composition_rows(
    const map<int, double>& comp,
    const map<pair<int, int>, map<pair<int, int>, pair<float, float> > >& allelecounts,
    const map<pair<int, int>, map<int, float> >& conditional_fractions){

    WeightedCompositionRowSelection out;
    auto complete_positive_row = [&](const pair<int,int>& key1,
                                     const pair<int,int>& key2,
                                     const pair<float,float>& counts)->bool{
        if (!((double)counts.first+(double)counts.second > 0.0)) return false;
        double expected = 0.0;
        return weighted_composition_expected_from_row_complete(
            comp,key1,key2,conditional_fractions,expected);
    };

    // Prefer all native component-pair rows that have complete conditional
    // support for every remaining positive composition component.
    for (auto ac = allelecounts.begin(); ac != allelecounts.end(); ++ac){
        for (auto ac2 = ac->second.begin(); ac2 != ac->second.end(); ++ac2){
            if (ac2->first.first < 0) continue;
            if (comp.count(ac->first.first) == 0 ||
                comp.count(ac2->first.first) == 0) continue;
            if (!complete_positive_row(ac->first,ac2->first,ac2->second)) continue;
            out.rows.push_back(make_pair(ac->first,ac2->first));
        }
    }
    if (!out.rows.empty()){
        out.direct_component_rows = true;
        out.status = "weighted_composition_component_pair_rows_complete_condf";
        return out;
    }

    // No complete direct component pair exists. Try ranked non-duplicating
    // bases and select the first basis with at least one complete row.
    vector<pair<int,int> > candidates =
        rank_weighted_composition_fallback_bases(comp,allelecounts);
    for (auto candidate = candidates.begin(); candidate != candidates.end(); ++candidate){
        vector<pair<pair<int,int>,pair<int,int> > > candidate_rows;
        for (auto ac = allelecounts.begin(); ac != allelecounts.end(); ++ac){
            if (ac->first.first != candidate->first) continue;
            for (auto ac2 = ac->second.begin(); ac2 != ac->second.end(); ++ac2){
                if (candidate->second < 0){
                    if (ac2->first.first != -1) continue;
                } else if (ac2->first.first != candidate->second){
                    continue;
                }
                if (!complete_positive_row(ac->first,ac2->first,ac2->second)) continue;
                candidate_rows.push_back(make_pair(ac->first,ac2->first));
            }
        }
        if (candidate_rows.empty()) continue;
        out.rows.swap(candidate_rows);
        out.fallback_basis = *candidate;
        if (candidate->second < 0){
            out.status = string(
                "weighted_composition_singlet_fallback_complete_condf:") +
                std::to_string(candidate->first);
        } else {
            out.status = string(
                "weighted_composition_pair_fallback_complete_condf:") +
                std::to_string(candidate->first) + "+" +
                std::to_string(candidate->second);
        }
        return out;
    }
    return out;
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
    compile_status_by_cell.clear();
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
        WeightedCompositionRowSelection selection =
            select_weighted_composition_rows(comp,allelecounts,expfracs);

        for (auto row = selection.rows.begin(); row != selection.rows.end(); ++row){
            auto ac = allelecounts.find(row->first);
            if (ac == allelecounts.end()) continue;
            auto ac2 = ac->second.find(row->second);
            if (ac2 == ac->second.end()) continue;

            double expected_raw = 0.0;
            if (!composition_expected_from_row_complete(
                    comp,row->first,row->second,expected_raw)){
                continue;
            }
            const double ref = ac2->second.first;
            const double alt = ac2->second.second;
            if (!(ref+alt > 0.0)) continue;

            // Weighted species-composition identities do not currently have a
            // generalized multi-r expression model. Use the biologically
            // specified dosage composition as p_e and estimate only c.
            n.push_back(ref+alt);
            k.push_back(alt);
            p_e.push_back(expected_raw);
            p_A.push_back(-1.0);
            p_B.push_back(-1.0);
            type1.push_back(row->first);
            type2.push_back(row->second);
        }

        compile_status_by_cell[barcode] = n.empty() ?
            "weighted_composition_no_complete_conditional_rows" : selection.status;
        return;
    }

    compile_status_by_cell[barcode] = "native_assignment_rows";

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
    amb_mu_loo_cell.clear();

    auto build_excluded_profile = [&](const set<int>& excluded, bool& undefined){
        map<pair<int,int>, map<pair<int,int>, double> > prof;
        SourceExclusionProfileAssessment assessment = apply_source_exclusion_profile(
            contam_prof, excluded, source_exclusion_strength);
        undefined = assessment.status == "degenerate_all_mass_excluded_fallback_global";

        // Exact lambda=0 endpoint: use the existing global expectation table
        // byte-for-byte rather than rebuilding it from the same profile.
        if (source_exclusion_strength == 0.0 || assessment.status == "no_eligible_sources"){
            return amb_mu;
        }
        if (undefined) return amb_mu;

        for (map<pair<int,int>, map<pair<int,int>, double> >::const_iterator t1 = amb_mu.begin();
            t1 != amb_mu.end(); ++t1){
            for (map<pair<int,int>,double>::const_iterator t2 = t1->second.begin();
                t2 != t1->second.end(); ++t2){
                double val = 0.0;
                map<pair<int,int>,map<int,float> >::const_iterator ef = expfracs.find(t1->first);
                if (ef != expfracs.end()){
                    for (map<int,double>::const_iterator cp = assessment.scoring_profile.begin();
                        cp != assessment.scoring_profile.end(); ++cp){
                        if (cp->first < 0 || !(cp->second > 0.0)) continue;
                        map<int,float>::const_iterator ev = ef->second.find(cp->first);
                        if (ev != ef->second.end()) val += cp->second*ev->second;
                    }
                }
                if (std::isfinite(val)){
                    prof[t1->first][t2->first] = std::min(1.0, std::max(0.0,val));
                } else {
                    prof[t1->first][t2->first] = t2->second;
                }
            }
        }
        return prof;
    };

    set<int> receiver_ids;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a) receiver_ids.insert(a->second);

    int n_undefined = 0;
    for (set<int>::iterator rid = receiver_ids.begin(); rid != receiver_ids.end(); ++rid){
        set<int> excluded;
        if (*rid < n_samples){
            excluded.insert(*rid);
        } else {
            pair<int,int> combo = idx_to_hap_comb(*rid, n_samples);
            excluded.insert(combo.first);
            excluded.insert(combo.second);
        }
        bool undefined = false;
        amb_mu_loo[*rid] = build_excluded_profile(excluded, undefined);
        if (undefined) n_undefined++;
    }

    int n_cell_specific = 0;
    int n_cell_undefined = 0;
    for (auto it = cell_composition_overrides.begin();
        it != cell_composition_overrides.end(); ++it){
        set<int> excluded = weighted_composition_source_exclusions(it->second);
        if (excluded.empty()) continue;
        bool undefined = false;
        amb_mu_loo_cell[it->first] = build_excluded_profile(excluded, undefined);
        n_cell_specific++;
        if (undefined) n_cell_undefined++;
    }

    fprintf(stderr, "  Computed source-exclusion profiles at lambda=%.6g for %lu receiver identities",
        source_exclusion_strength, receiver_ids.size());
    if (n_undefined > 0){
        fprintf(stderr, " (%d undefined because exclusion removed all profile mass)", n_undefined);
    }
    if (n_cell_specific > 0){
        fprintf(stderr, "; %d weighted-composition cells use component-complete exclusion",
            n_cell_specific);
        if (n_cell_undefined > 0){
            fprintf(stderr, " (%d cell-specific exclusions undefined)", n_cell_undefined);
        }
    }
    fprintf(stderr, "\n");
}

set<int> contamFinder3::source_exclusions_for_cell(unsigned long barcode) const{
    auto override_it = cell_composition_overrides.find(barcode);
    if (override_it != cell_composition_overrides.end()){
        return weighted_composition_source_exclusions(override_it->second);
    }
    set<int> excluded;
    auto ai = assn.find(barcode);
    if (ai == assn.end()) return excluded;
    if (ai->second < n_samples){
        excluded.insert(ai->second);
    } else {
        pair<int,int> combo = idx_to_hap_comb(ai->second, n_samples);
        excluded.insert(combo.first);
        excluded.insert(combo.second);
    }
    return excluded;
}

bool contamFinder3::compile_cell_model_data(unsigned long barcode,
    const vector<int>& obs_idx, vector<double>& n, vector<double>& k,
    vector<double>& p_e, vector<double>& p_c_scoring, vector<double>& p_c_global,
    vector<double>& cell_p_A, vector<double>& cell_p_B,
    vector<double>& observation_weights, bool& is_heterotypic,
    string& status) const{

    n.clear(); k.clear(); p_e.clear(); p_c_scoring.clear(); p_c_global.clear();
    cell_p_A.clear(); cell_p_B.clear(); observation_weights.clear();
    status = "not_compiled";

    int cell_ident = -1;
    auto ai = assn.find(barcode);
    if (ai != assn.end()) cell_ident = ai->second;
    const bool have_cell_excluded = use_loo && amb_mu_loo_cell.count(barcode) > 0;
    const bool have_identity_excluded = use_loo && cell_ident >= 0 && amb_mu_loo.count(cell_ident) > 0;

    is_heterotypic = false;
    if (cell_ident >= n_samples){
        pair<int,int> combo = idx_to_hap_comb(cell_ident, n_samples);
        is_heterotypic = combo.first != combo.second;
    }

    auto profile_value = [&](const map<pair<int,int>, map<pair<int,int>, double> >& profile,
        const pair<int,int>& t1, const pair<int,int>& t2, double fallback) -> double {
        auto a1 = profile.find(t1);
        if (a1 != profile.end()){
            auto a2 = a1->second.find(t2);
            if (a2 != a1->second.end()) return a2->second;
        }
        return fallback;
    };
    auto get_global_pc = [&](int obs_index) -> double {
        return profile_value(amb_mu, type1_all[obs_index], type2_all[obs_index], 0.5);
    };
    auto get_scoring_pc = [&](int obs_index) -> double {
        const pair<int,int>& t1 = type1_all[obs_index];
        const pair<int,int>& t2 = type2_all[obs_index];
        if (have_cell_excluded){
            auto cell_it = amb_mu_loo_cell.find(barcode);
            if (cell_it != amb_mu_loo_cell.end()){
                const double value = profile_value(cell_it->second, t1, t2,
                    std::numeric_limits<double>::quiet_NaN());
                if (std::isfinite(value)) return value;
            }
        }
        if (have_identity_excluded){
            auto ident_it = amb_mu_loo.find(cell_ident);
            if (ident_it != amb_mu_loo.end()){
                const double value = profile_value(ident_it->second, t1, t2,
                    std::numeric_limits<double>::quiet_NaN());
                if (std::isfinite(value)) return value;
            }
        }
        return get_global_pc(obs_index);
    };

    for (vector<int>::const_iterator i = obs_idx.begin(); i != obs_idx.end(); ++i){
        const double this_p_e = adjust_p_err(p_e_all[*i], e_r, e_a);
        const double this_p_c = get_scoring_pc(*i);
        const double global_p_c = get_global_pc(*i);
        double this_p_A = -1.0, this_p_B = -1.0;
        if (is_heterotypic && p_A_all[*i] >= 0){
            this_p_A = adjust_p_err(p_A_all[*i], e_r, e_a);
            this_p_B = adjust_p_err(p_B_all[*i], e_r, e_a);
        }

        if (use_fi_weight){
            const double w = compute_fi_weight(this_p_e, this_p_c);
            if (w < 1e-12) continue;
            n.push_back(n_all[*i]); k.push_back(k_all[*i]); p_e.push_back(this_p_e);
            p_c_scoring.push_back(this_p_c); p_c_global.push_back(global_p_c);
            cell_p_A.push_back(this_p_A); cell_p_B.push_back(this_p_B);
            observation_weights.push_back(w);
        } else {
            const bool filter_required = !is_heterotypic && tetraploid_aware &&
                type2_all[*i].first != -1;
            if (!contamination_observation_passes_model_filter(filter_required,
                    amb_mu_available, this_p_e, this_p_c, min_signal_gap,
                    type1_all[*i].second, type2_all[*i].second)) continue;
            n.push_back(n_all[*i]); k.push_back(k_all[*i]); p_e.push_back(this_p_e);
            p_c_scoring.push_back(this_p_c); p_c_global.push_back(global_p_c);
            cell_p_A.push_back(this_p_A); cell_p_B.push_back(this_p_B);
        }
    }

    if (n.empty() && (tetraploid_aware || use_fi_weight)){
        observation_weights.clear(); cell_p_A.clear(); cell_p_B.clear();
        p_c_global.clear(); p_c_scoring.clear(); p_e.clear(); k.clear(); n.clear();
        for (vector<int>::const_iterator i = obs_idx.begin(); i != obs_idx.end(); ++i){
            const double this_p_e = adjust_p_err(p_e_all[*i], e_r, e_a);
            const double this_p_c = get_scoring_pc(*i);
            n.push_back(n_all[*i]); k.push_back(k_all[*i]); p_e.push_back(this_p_e);
            p_c_scoring.push_back(this_p_c); p_c_global.push_back(get_global_pc(*i));
            double this_p_A = -1.0, this_p_B = -1.0;
            if (is_heterotypic && p_A_all[*i] >= 0){
                this_p_A = adjust_p_err(p_A_all[*i], e_r, e_a);
                this_p_B = adjust_p_err(p_B_all[*i], e_r, e_a);
            }
            cell_p_A.push_back(this_p_A); cell_p_B.push_back(this_p_B);
            if (use_fi_weight) observation_weights.push_back(compute_fi_weight(this_p_e, this_p_c));
        }
        status = n.empty() ? "no_observations" : "fallback_all_compiled_rows";
    } else {
        status = n.empty() ? "no_observations" : "ok";
    }
    return !n.empty();
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
    const vector<int>& obs_idx,
    bool apply_prior,
    double prior_mean,
    double prior_var,
    const string& prior_group){

    const double nanv = std::numeric_limits<double>::quiet_NaN();
    CellContamFitResult out;
    out.barcode = barcode;
    out.c = nanv;
    out.c_se = nanv;
    out.ll = -INFINITY;
    out.r = nanv;
    out.r_se = nanv;
    out.c_ci_low = nanv;
    out.c_ci_high = nanv;
    out.r_ci_low = nanv;
    out.r_ci_high = nanv;
    out.r_c_correlation = nanv;
    out.ridge_width = nanv;
    out.gradient_norm = nanv;
    out.informative_weight = 0.0;
    out.pearson_ss = 0.0;
    out.pearson_df = 0.0;
    out.diagnostic_fixed_r_fallback_c = nanv;
    out.diagnostic_fixed_r_fallback_ll = nanv;
    out.c_profile_candidate_gap = nanv;
    out.c_profile_neighbor_gain = nanv;
    out.c_profile_global_gain = nanv;
    out.c_profile_likelihood_span = nanv;
    out.c_profile_boundary_slope = nanv;
    out.c_profile_support_low = nanv;
    out.c_profile_support_high = nanv;
    out.c_profile_support_width = nanv;
    out.r_profile_support_low = nanv;
    out.r_profile_support_high = nanv;
    out.r_profile_support_width = nanv;
    out.raw_candidate_r = nanv;
    out.raw_candidate_c = nanv;
    out.raw_candidate_ll = nanv;
    out.validated_candidate_r = nanv;
    out.validated_candidate_c = nanv;
    out.validated_candidate_ll = nanv;
    out.starts_evaluated = 0;
    out.starts_optimized = 0;
    out.solver_successful_starts = 0;
    out.finite_in_domain_bfgs_candidates = 0;
    out.profile_optimization_attempted = 0;
    out.profile_optimization_succeeded = 0;
    out.profile_validation_passed = 0;
    out.r_profile_optimization_succeeded = 0;
    out.r_profile_identifiable = 0;
    out.optimizer_iterations = -1;
    out.is_heterotypic = false;
    out.bfgs_fallback = false;
    out.has_allele_ratio = false;
    out.fit_success = false;
    out.interval_success = false;
    out.hit_c_lower = false;
    out.hit_c_upper = false;
    out.hit_r_lower = false;
    out.hit_r_upper = false;
    out.prior_applied = apply_prior;
    out.interval_method = "not_applicable";
    out.optimizer_status = "not_run";
    out.r_validation_status = "r_profile_not_attempted";
    out.diagnostic_fallback_status = "not_run";
    out.prior_group = prior_group;

    vector<double> n;
    vector<double> k;
    vector<double> p_e;
    vector<double> p_c;
    vector<double> p_c_global;
    vector<double> cell_p_A;
    vector<double> cell_p_B;
    vector<double> fi_weights;
    bool is_heterotypic = false;
    string model_compile_status;
    if (!compile_cell_model_data(barcode, obs_idx, n, k, p_e, p_c, p_c_global,
            cell_p_A, cell_p_B, fi_weights, is_heterotypic, model_compile_status)){
        out.optimizer_status = model_compile_status;
        return out;
    }
    out.is_heterotypic = is_heterotypic;

    vector<double> obs_weights;
    if (use_fi_weight && fi_weights.size() == n.size()) obs_weights = fi_weights;
    for (size_t i = 0; i < n.size(); ++i){
        const double w = obs_weights.empty() ? 1.0 : obs_weights[i];
        out.informative_weight += w * n[i];
    }

    double prior_alpha = 0.0;
    double prior_beta = 0.0;
    bool prior_valid = false;
    if (apply_prior){
        prior_valid = contamination_beta_prior_is_log_concave(
            prior_mean, prior_var, prior_alpha, prior_beta);
        if (!prior_valid){
            out.optimizer_status = "regularized_prior_not_log_concave";
            return out;
        }
    }

    if (is_heterotypic && !cell_p_A.empty() && cell_p_A[0] >= 0){

        // ================= FIXED-r MODE (pooled-r experiment) =================
        // If an allele ratio has been supplied externally for this cell's
        // identity, r is not a free parameter. Build the endogenous expectation
        // from the fixed r and solve c alone, in one dimension.
        //
        // adjust_p_err() is affine in p, so it commutes with the convex
        // combination:  adjust(r*pA + (1-r)*pB) == r*adjust(pA) + (1-r)*adjust(pB).
        // cell_p_A and cell_p_B are already error-adjusted, so mixing them
        // directly is exact and no re-adjustment is required.
        //
        // The point of this mode: contamination absorption requires r to track
        // this cell's own c. A fixed r cannot. Whether that fixes anything
        // depends entirely on the pool the value came from, which is the whole
        // experiment.
        double r_fixed = -1.0;
        if (fix_r_enabled){
            int ident_for_r = -1;
            if (assn.count(barcode) > 0) ident_for_r = assn[barcode];
            map<int, double>::const_iterator fr = fixed_r_by_ident.find(ident_for_r);
            if (fr != fixed_r_by_ident.end()) r_fixed = fr->second;
        }

        if (r_fixed >= 0.0){
            const double rr = clamp_unit_interval(r_fixed);
            vector<double> p_e_fixed;
            p_e_fixed.reserve(n.size());
            for (size_t i = 0; i < n.size(); ++i){
                p_e_fixed.push_back(rr * cell_p_A[i] + (1.0 - rr) * cell_p_B[i]);
            }

            out.is_heterotypic = true;
            out.starts_evaluated = 1;

            optimML::brent_solver c_only(ll_c, dll_dc, d2ll_dc2);
            c_only.add_data("n", n);
            c_only.add_data("k", k);
            c_only.add_data("p_e", p_e_fixed);
            c_only.add_data("p_c", p_c);
            c_only.constrain_01();
            if (!obs_weights.empty()) c_only.add_weights(obs_weights);
            if (prior_valid) c_only.add_beta_prior(prior_alpha, prior_beta);
            c_only.set_maxiter(-1);

            UnivariateCInformationAssessment c_info =
                assess_univariate_c_information(n, p_e_fixed, p_c, obs_weights);
            if (!c_info.identifiable){
                out.optimizer_status = c_info.status;
                return out;
            }

            try {
                double c_hat = c_only.solve(0, 1);
                bool accepted = c_only.root_found;
                UnivariateCBoundaryAssessment boundary;
                if (!accepted){
                    boundary = assess_unpenalized_c_boundary(c_hat, n, k,
                        p_e_fixed, p_c, obs_weights, apply_prior);
                    accepted = boundary.accepted;
                }
                if (accepted){
                    out.r = rr;
                    out.c = clamp_unit_interval(c_hat);
                    out.ll = c_only.log_likelihood;
                    // r was not estimated. Report zero uncertainty on it, which
                    // is honest: it carries no error because it carries no fit.
                    out.r_se = 0.0;
                    out.r_ci_low = rr;
                    out.r_ci_high = rr;
                    out.r_c_correlation = 0.0;
                    out.has_allele_ratio = true;
                    out.r_profile_optimization_succeeded = 1;
                    out.r_profile_identifiable = 1;
                    out.r_profile_support_low = rr;
                    out.r_profile_support_high = rr;
                    out.r_profile_support_width = 0.0;
                    out.r_validation_status = "fixed_r_supplied_not_estimated";
                    out.fit_success = true;
                    out.starts_optimized = 1;
                    out.optimizer_status = c_only.root_found ?
                        "fixed_r_brent" : boundary.status;
                    out.interval_method = c_only.root_found ?
                        "brent_c_only" : "supported_boundary_no_curvature_interval";
                    out.hit_c_lower = out.c <= 1e-4;
                    out.hit_c_upper = out.c >= 1.0 - 1e-4;
                    out.hit_r_lower = false;
                    out.hit_r_upper = false;
                    if (c_only.root_found && c_only.se_found && std::isfinite(c_only.se)){
                        out.c_se = c_only.se;
                        out.c_ci_low  = std::max(0.0, out.c - 1.959963984540054*out.c_se);
                        out.c_ci_high = std::min(1.0, out.c + 1.959963984540054*out.c_se);
                        out.interval_success = true;
                    }
                } else {
                    out.optimizer_status = boundary.status;
                }
            } catch (...) {
                out.optimizer_status = "fixed_r_brent_failed";
            }
            return out;
        }
        // =============== end FIXED-r MODE; free-r path follows ===============

        struct GridStart { double r; double c; double obj; };
        vector<GridStart> grid;
        const double r_grid[] = {0.15, 0.30, 0.50, 0.70, 0.85};
        const double c_grid[] = {0.02, 0.10, 0.25, 0.45, 0.70};

        if (heterotypic_start_mode == START_SINGLE){
            const double c0 = prior_valid ? prior_mean : 0.15;
            grid.push_back({0.5, clamp_unit_interval(c0),
                eval_three_dataset_objective(0.5, clamp_unit_interval(c0), n, k,
                    cell_p_A, cell_p_B, p_c, obs_weights, prior_valid,
                    prior_mean, prior_var)});
        } else {
            for (double rr : r_grid){
                for (double cc : c_grid){
                    grid.push_back({rr, cc,
                        eval_three_dataset_objective(rr, cc, n, k, cell_p_A,
                            cell_p_B, p_c, obs_weights, prior_valid,
                            prior_mean, prior_var)});
                }
            }
        }
        out.starts_evaluated = (int)grid.size();
        std::sort(grid.begin(), grid.end(), [](const GridStart& a, const GridStart& b){
            return a.obj > b.obj;
        });

        int desired = 1;
        if (heterotypic_start_mode == START_EXHAUSTIVE || thorough_multistart){
            desired = (int)grid.size();
        } else if (heterotypic_start_mode == START_TOPK){
            desired = std::min((int)grid.size(), heterotypic_start_top_k);
        }

        vector<GridStart> starts;
        for (size_t gi = 0; gi < grid.size() && (int)starts.size() < desired; ++gi){
            bool duplicate = false;
            for (size_t sj = 0; sj < starts.size(); ++sj){
                const double dr = grid[gi].r - starts[sj].r;
                const double dc = grid[gi].c - starts[sj].c;
                if (sqrt(dr*dr + dc*dc) < 0.08){
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate || heterotypic_start_mode == START_EXHAUSTIVE){
                starts.push_back(grid[gi]);
            }
        }
        if (starts.empty()) starts.push_back(grid.front());

        vector<pair<double,double> > bfgs_candidates;
        int solver_reported_successes = 0;
        for (size_t si = 0; si < starts.size(); ++si){
            vector<double> init = {starts[si].r, starts[si].c};
            optimML::multivar_ml_solver solver(init, ll_three, dll_three);
            solver.add_data("n", n);
            solver.add_data("k", k);
            solver.add_data("p_A", cell_p_A);
            solver.add_data("p_B", cell_p_B);
            solver.add_data("p_c", p_c);
            solver.constrain_01(0);
            solver.constrain_01(1);
            if (!obs_weights.empty()) solver.add_weights(obs_weights);
            if (prior_valid) solver.add_beta_prior(1, prior_alpha, prior_beta);
            solver.set_delta(1e-7);
            solver.set_silent(true);
            try {
                bool success = solver.solve();
                out.starts_optimized++;
                if (success && solver.results.size() >= 2){
                    solver_reported_successes++;
                    // Preserve the raw coordinates even when malformed. The
                    // profile-model helper independently filters them and still
                    // runs the exact free model if every BFGS result is unusable.
                    bfgs_candidates.push_back(make_pair(
                        solver.results[0],solver.results[1]));
                }
            } catch (...) {
                out.starts_optimized++;
            }
        }

        out.solver_successful_starts=solver_reported_successes;
        TwoDimensionalFitAssessment best_assessment=
            fit_two_dimensional_profile_model(bfgs_candidates,n,k,cell_p_A,
                cell_p_B,p_c,obs_weights,prior_valid,prior_mean,prior_var);
        const int profile_validation_passed=best_assessment.accepted ? 1 : 0;
        out.finite_in_domain_bfgs_candidates=
            best_assessment.finite_in_domain_raw_candidates;
        out.profile_optimization_attempted=
            best_assessment.profile_optimization_attempted ? 1 : 0;
        out.profile_optimization_succeeded=
            best_assessment.profile_optimization_succeeded ? 1 : 0;
        out.profile_validation_passed=profile_validation_passed;
        out.r_profile_optimization_succeeded=
            best_assessment.r_profile_optimization_succeeded ? 1 : 0;
        out.r_profile_identifiable=best_assessment.r_profile_identifiable ? 1 : 0;
        out.r_validation_status=best_assessment.r_validation_status;
        if (best_assessment.raw_candidate_available){
            out.raw_candidate_r=best_assessment.raw_candidate_r;
            out.raw_candidate_c=best_assessment.raw_candidate_c;
            out.raw_candidate_ll=best_assessment.raw_candidate_objective;
        }
        if (best_assessment.r_profile_identifiable &&
            std::isfinite(best_assessment.validated_r))
            out.validated_candidate_r=best_assessment.validated_r;
        if (std::isfinite(best_assessment.validated_c))
            out.validated_candidate_c=best_assessment.validated_c;
        if (std::isfinite(best_assessment.validated_objective))
            out.validated_candidate_ll=best_assessment.validated_objective;

        double best_ll=best_assessment.accepted ?
            best_assessment.validated_objective : -INFINITY;
        double best_r=best_assessment.accepted ?
            (std::isfinite(best_assessment.validated_r) ?
                best_assessment.validated_r : best_assessment.nuisance_r_argmax) : nanv;
        double best_c=best_assessment.accepted ?
            best_assessment.validated_c : nanv;
        const string best_rejection_status=best_assessment.status;
        // This flag now means that the exact profile model supplied the result
        // without a usable BFGS candidate. It is not the forbidden fixed-r
        // diagnostic fallback.
        out.bfgs_fallback=best_assessment.accepted &&
            best_assessment.finite_in_domain_raw_candidates==0;
        bool diagnostic_fixed_r_candidate_available = false;
        if (profile_validation_passed == 0){
            out.optimizer_status = best_rejection_status;
            out.c_profile_candidate_gap=best_assessment.candidate_profile_gap;
            out.c_profile_neighbor_gain=best_assessment.profile_neighbor_gain;
            out.c_profile_global_gain=best_assessment.profile_global_gain;
            out.c_profile_likelihood_span=best_assessment.profile_likelihood_span;
            out.c_profile_boundary_slope=
                best_assessment.boundary.derivative_near_boundary;
            out.c_profile_support_low=best_assessment.profile_support_low;
            out.c_profile_support_high=best_assessment.profile_support_high;
            out.c_profile_support_width=best_assessment.profile_support_width;
            out.r_profile_support_low=best_assessment.r_profile_support_low;
            out.r_profile_support_high=best_assessment.r_profile_support_high;
            out.r_profile_support_width=best_assessment.r_profile_support_width;

            // The fixed-r solve is a diagnostic candidate only.  It is a
            // different model from the requested free-(r,c) fit and therefore
            // cannot supply the production estimate after all 2-D starts fail.
            optimML::brent_solver c_fallback(ll_c, dll_dc, d2ll_dc2);
            c_fallback.add_data("n", n);
            c_fallback.add_data("k", k);
            c_fallback.add_data("p_e", p_e);
            c_fallback.add_data("p_c", p_c);
            c_fallback.constrain_01();
            if (!obs_weights.empty()) c_fallback.add_weights(obs_weights);
            c_fallback.set_maxiter(-1);
            try {
                double diagnostic_c = c_fallback.solve(0,1);
                out.diagnostic_fixed_r_fallback_c =
                    clamp_unit_interval(diagnostic_c);
                out.diagnostic_fixed_r_fallback_ll = c_fallback.log_likelihood;
                if (c_fallback.root_found){
                    out.diagnostic_fallback_status =
                        "diagnostic_fixed_r_interior_candidate";
                } else {
                    UnivariateCBoundaryAssessment diagnostic_boundary =
                        assess_unpenalized_c_boundary(diagnostic_c, n, k,
                            p_e, p_c, obs_weights, false);
                    out.diagnostic_fallback_status = string(
                        "diagnostic_fixed_r_") + diagnostic_boundary.status;
                }
                diagnostic_fixed_r_candidate_available = true;
            } catch (...) {
                out.diagnostic_fallback_status =
                    "diagnostic_fixed_r_solver_failed";
            }
        } else {
            out.optimizer_status = best_assessment.status;
            out.c_profile_candidate_gap=best_assessment.candidate_profile_gap;
            out.c_profile_neighbor_gain=best_assessment.profile_neighbor_gain;
            out.c_profile_global_gain=best_assessment.profile_global_gain;
            out.c_profile_likelihood_span=best_assessment.profile_likelihood_span;
            out.c_profile_boundary_slope=
                best_assessment.boundary.derivative_near_boundary;
            out.c_profile_support_low=best_assessment.profile_support_low;
            out.c_profile_support_high=best_assessment.profile_support_high;
            out.c_profile_support_width=best_assessment.profile_support_width;
            out.r_profile_support_low=best_assessment.r_profile_support_low;
            out.r_profile_support_high=best_assessment.r_profile_support_high;
            out.r_profile_support_width=best_assessment.r_profile_support_width;
        }

        HeterotypicFitSource fit_source = choose_heterotypic_fit_source(
            profile_validation_passed,diagnostic_fixed_r_candidate_available);
        out.fit_success = heterotypic_fit_source_is_production_valid(fit_source);

        // A diagnostic fixed-r candidate never changes fit_success.  This
        // explicit policy guard prevents future refactors from silently turning
        // a surrogate model into the selected free-(r,c) production fit.
        if (!heterotypic_fit_source_is_production_valid(fit_source)){
            out.fit_success = false;
        }

        if (out.fit_success){
            const double internal_r=clamp_closed_unit_interval(best_r);
            out.r = best_assessment.r_profile_identifiable ? internal_r : nanv;
            out.c = clamp_closed_unit_interval(best_c);
            out.ll = best_ll;
            out.has_allele_ratio = best_assessment.r_profile_identifiable;
            out.hit_c_lower = out.c <= 1e-4;
            out.hit_c_upper = out.c >= 1.0-1e-4;
            out.hit_r_lower = out.has_allele_ratio && out.r <= 1e-4;
            out.hit_r_upper = out.has_allele_ratio && out.r >= 1.0-1e-4;
            out.c_ci_low=best_assessment.profile_support_low;
            out.c_ci_high=best_assessment.profile_support_high;
            if (std::isfinite(out.c_ci_low) && std::isfinite(out.c_ci_high)){
                out.interval_success=true;
                out.interval_method="continuous_c_profile";
            }
            if (out.has_allele_ratio){
                out.r_ci_low=best_assessment.r_profile_support_low;
                out.r_ci_high=best_assessment.r_profile_support_high;
            }

            // Store the constrained projected-gradient violation as a
            // diagnostic. It is not a production acceptance gate because its
            // absolute scale grows with observation depth.
            out.gradient_norm=best_assessment.projected_gradient_norm;

            double hse_r = nanv;
            double hse_c = nanv;
            double corr = nanv;
            bool hessian_ok = out.has_allele_ratio && numerical_three_hessian(out.r, out.c, n, k,
                cell_p_A, cell_p_B, p_c, obs_weights, prior_valid,
                prior_mean, prior_var, hse_r, hse_c, corr);
            if (hessian_ok){
                out.r_se = hse_r;
                out.c_se = hse_c;
                out.r_c_correlation = corr;
            }
            if (out.interval_success){
                out.ridge_width = out.c_ci_high - out.c_ci_low;
            }

            for (size_t i = 0; i < n.size(); ++i){
                const double pendo = internal_r*cell_p_A[i] + (1.0-internal_r)*cell_p_B[i];
                const double pp = clamp_unit_interval((1.0-out.c)*pendo + out.c*p_c[i], 1e-6);
                const double v = std::max(1e-12, n[i]*pp*(1.0-pp));
                const double resid = (k[i]-n[i]*pp)/sqrt(v);
                const double w = obs_weights.empty() ? 1.0 : obs_weights[i];
                out.pearson_ss += w*resid*resid;
                out.pearson_df += w;
            }
            out.pearson_df = std::max(0.0, out.pearson_df - 2.0);
        }
    } else {
        optimML::brent_solver c_cell(ll_c, dll_dc, d2ll_dc2);
        c_cell.add_data("n", n);
        c_cell.add_data("k", k);
        c_cell.add_data("p_e", p_e);
        c_cell.add_data("p_c", p_c);
        c_cell.constrain_01();
        if (!obs_weights.empty()) c_cell.add_weights(obs_weights);
        if (prior_valid) c_cell.add_beta_prior(prior_alpha, prior_beta);
        c_cell.set_maxiter(-1);
        UnivariateCInformationAssessment c_info =
            assess_univariate_c_information(n, p_e, p_c, obs_weights);
        if (!c_info.identifiable){
            out.optimizer_status = c_info.status;
            return out;
        }

        try {
            out.c = c_cell.solve(0,1);
            bool accepted = c_cell.root_found;
            UnivariateCBoundaryAssessment boundary;
            if (!accepted){
                boundary = assess_unpenalized_c_boundary(out.c, n, k,
                    p_e, p_c, obs_weights, apply_prior);
                accepted = boundary.accepted;
            }
            if (accepted){
                out.c = clamp_unit_interval(out.c);
                out.fit_success = true;
                out.ll = c_cell.log_likelihood;
                out.optimizer_status = c_cell.root_found ?
                    "brent_converged" : boundary.status;
                if (c_cell.root_found && c_cell.se_found && std::isfinite(c_cell.se)){
                    out.c_se = c_cell.se;
                    out.c_ci_low = std::max(0.0, out.c - 1.959963984540054*out.c_se);
                    out.c_ci_high = std::min(1.0, out.c + 1.959963984540054*out.c_se);
                    out.ridge_width = out.c_ci_high - out.c_ci_low;
                    out.interval_success = true;
                    out.interval_method = "brent_curvature";
                } else if (!c_cell.root_found) {
                    out.interval_method = "supported_boundary_no_curvature_interval";
                }
                out.hit_c_lower = out.c <= 1e-4;
                out.hit_c_upper = out.c >= 1.0-1e-4;
                double grad = 0.0;
                for (size_t i = 0; i < n.size(); ++i){
                    map<string,double> dd;
                    map<string,int> di;
                    dd["n"] = n[i]; dd["k"] = k[i]; dd["p_e"] = p_e[i]; dd["p_c"] = p_c[i];
                    const double w = obs_weights.empty() ? 1.0 : obs_weights[i];
                    grad += w*dll_dc(out.c,dd,di);
                    const double pp = clamp_unit_interval((1.0-out.c)*p_e[i] + out.c*p_c[i],1e-6);
                    const double v = std::max(1e-12,n[i]*pp*(1.0-pp));
                    const double resid = (k[i]-n[i]*pp)/sqrt(v);
                    out.pearson_ss += w*resid*resid;
                    out.pearson_df += w;
                }
                if (prior_valid) grad += beta_log_gradient(out.c, prior_mean, prior_var);
                out.gradient_norm = fabs(grad);
                out.pearson_df = std::max(0.0,out.pearson_df-1.0);
            } else {
                out.optimizer_status = boundary.status;
            }
        } catch (...) {
            out.optimizer_status = "brent_failed";
        }
    }

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

    contam_rate_mle.clear();
    contam_rate_regularized.clear();
    allele_ratio_mle.clear();
    allele_ratio_regularized.clear();
    contam_ci_low.clear();
    contam_ci_high.clear();
    allele_ratio_ci_low.clear();
    allele_ratio_ci_high.clear();
    r_c_correlation.clear();
    ridge_width.clear();
    prior_mean_by_cell.clear();
    prior_var_by_cell.clear();
    prior_displacement.clear();
    mle_log_likelihood.clear();
    regularized_log_likelihood.clear();
    informative_weight_by_cell.clear();
    gradient_norm_by_cell.clear();
    diagnostic_fixed_r_fallback_c_by_cell.clear();
    diagnostic_fixed_r_fallback_ll_by_cell.clear();
    c_profile_candidate_gap_by_cell.clear();
    c_profile_neighbor_gain_by_cell.clear();
    c_profile_global_gain_by_cell.clear();
    c_profile_likelihood_span_by_cell.clear();
    c_profile_boundary_slope_by_cell.clear();
    c_profile_support_low_by_cell.clear();
    c_profile_support_high_by_cell.clear();
    c_profile_support_width_by_cell.clear();
    r_profile_support_low_by_cell.clear();
    r_profile_support_high_by_cell.clear();
    r_profile_support_width_by_cell.clear();
    raw_candidate_r_by_cell.clear();
    raw_candidate_c_by_cell.clear();
    raw_candidate_ll_by_cell.clear();
    validated_candidate_r_by_cell.clear();
    validated_candidate_c_by_cell.clear();
    validated_candidate_ll_by_cell.clear();
    starts_evaluated_by_cell.clear();
    starts_optimized_by_cell.clear();
    solver_successful_starts_by_cell.clear();
    finite_in_domain_bfgs_candidates_by_cell.clear();
    profile_optimization_attempted_by_cell.clear();
    profile_optimization_succeeded_by_cell.clear();
    profile_validation_passed_by_cell.clear();
    r_profile_optimization_succeeded_by_cell.clear();
    r_profile_identifiable_by_cell.clear();
    profile_validation_status_by_cell.clear();
    r_validation_status_by_cell.clear();
    diagnostic_fallback_status_by_cell.clear();
    prior_training_eligible.clear();
    fit_boundary_flags.clear();
    prior_group_by_cell.clear();
    interval_method_by_cell.clear();
    optimizer_status_by_cell.clear();
    prior_training_reason_by_cell.clear();

    vector<pair<unsigned long, vector<int> > > cells;
    cells.reserve(cell_to_idx.size());
    for (map<unsigned long, vector<int> >::iterator ci = cell_to_idx.begin();
        ci != cell_to_idx.end(); ++ci){
        cells.push_back(*ci);
    }

    const int nt = (num_threads > 1) ? num_threads : 1;
    vector<CellContamFitResult> mle_results(cells.size());

    #pragma omp parallel for num_threads(nt) schedule(dynamic, 8)
    for (int i = 0; i < (int)cells.size(); ++i){
        mle_results[i] = fit_one_contam_cell(cells[i].first, cells[i].second,
            false, -1.0, -1.0, "none");
    }

    if (!user_prior_set && contam_prior_mode != PRIOR_NONE){
        learn_contam_priors(mle_results);
    } else {
        for (vector<CellContamFitResult>::const_iterator it = mle_results.begin();
            it != mle_results.end(); ++it){
            string reason;
            bool ok = prior_training_supported(*it, reason);
            prior_training_eligible[it->barcode] = ok ? 1 : 0;
            prior_training_reason_by_cell[it->barcode] = reason;
        }
    }

    vector<CellContamFitResult> reg_results = mle_results;
    vector<int> reg_fit_selected(cells.size(), 0);
    vector<string> reg_attempt_status(cells.size(), "not_run");
    vector<double> reg_prior_mean(cells.size(), -1.0);
    vector<double> reg_prior_var(cells.size(), -1.0);
    vector<string> reg_prior_group(cells.size(), "none");

    #pragma omp parallel for num_threads(nt) schedule(dynamic, 8)
    for (int i = 0; i < (int)cells.size(); ++i){
        string group;
        pair<double,double> pv = prior_for_cell(cells[i].first,
            mle_results[i].is_heterotypic, group);
        reg_prior_mean[i] = pv.first;
        reg_prior_var[i] = pv.second;
        reg_prior_group[i] = group;
        if (pv.first > 0.0 && pv.second > 0.0){
            CellContamFitResult rr = fit_one_contam_cell(cells[i].first,
                cells[i].second, true, pv.first, pv.second, group);
            reg_attempt_status[i] = rr.optimizer_status;
            if (rr.fit_success && !regularized_c_candidate_is_selectable(
                    rr.fit_success,rr.c)){
                reg_attempt_status[i] = "regularized_boundary_not_selectable";
            } else if (rr.fit_success && mle_results[i].fit_success){
                reg_results[i] = rr;
                reg_fit_selected[i] = 1;
            }
        }
    }

    double pearson_ss_total = 0.0;
    double pearson_df_total = 0.0;
    for (vector<CellContamFitResult>::const_iterator it = mle_results.begin();
        it != mle_results.end(); ++it){
        if (std::isfinite(it->pearson_ss) && it->pearson_ss >= 0.0 &&
            std::isfinite(it->pearson_df) && it->pearson_df > 0.0){
            pearson_ss_total += it->pearson_ss;
            pearson_df_total += it->pearson_df;
        }
    }
    dispersion_phi = (pearson_df_total > 0.0) ? pearson_ss_total / pearson_df_total : 1.0;
    if (!std::isfinite(dispersion_phi) || dispersion_phi < 1.0) dispersion_phi = 1.0;
    const double dispersion_scale = sqrt(dispersion_phi);

    vector<double> selected_c;
    vector<double> selected_r;
    int n_three_comp = 0;
    int n_two_comp = 0;
    int n_failed = 0;
    int n_fallback = 0;

    for (size_t i = 0; i < cells.size(); ++i){
        const CellContamFitResult& mle = mle_results[i];
        const CellContamFitResult& reg = reg_results[i];
        const CellContamFitResult& selected = reg_fit_selected[i] ? reg : mle;
        const unsigned long bc = cells[i].first;

        starts_evaluated_by_cell[bc]=mle.starts_evaluated;
        starts_optimized_by_cell[bc]=mle.starts_optimized;
        solver_successful_starts_by_cell[bc]=mle.solver_successful_starts;
        finite_in_domain_bfgs_candidates_by_cell[bc]=
            mle.finite_in_domain_bfgs_candidates;
        profile_optimization_attempted_by_cell[bc]=
            mle.profile_optimization_attempted;
        profile_optimization_succeeded_by_cell[bc]=
            mle.profile_optimization_succeeded;
        profile_validation_passed_by_cell[bc]=
            mle.profile_validation_passed;
        r_profile_optimization_succeeded_by_cell[bc]=
            mle.r_profile_optimization_succeeded;
        r_profile_identifiable_by_cell[bc]=mle.r_profile_identifiable;
        profile_validation_status_by_cell[bc]=mle.optimizer_status;
        r_validation_status_by_cell[bc]=mle.r_validation_status;
        if (std::isfinite(mle.raw_candidate_r))
            raw_candidate_r_by_cell[bc]=mle.raw_candidate_r;
        if (std::isfinite(mle.raw_candidate_c))
            raw_candidate_c_by_cell[bc]=mle.raw_candidate_c;
        if (std::isfinite(mle.raw_candidate_ll))
            raw_candidate_ll_by_cell[bc]=mle.raw_candidate_ll;
        // Same authority as the allele-ratio outputs: a validated r is only
        // publishable when the r profile is identifiable.
        if (mle.r_profile_identifiable && std::isfinite(mle.validated_candidate_r))
            validated_candidate_r_by_cell[bc]=mle.validated_candidate_r;
        if (std::isfinite(mle.validated_candidate_c))
            validated_candidate_c_by_cell[bc]=mle.validated_candidate_c;
        if (std::isfinite(mle.validated_candidate_ll))
            validated_candidate_ll_by_cell[bc]=mle.validated_candidate_ll;
        if (std::isfinite(mle.c_profile_support_low))
            c_profile_support_low_by_cell[bc]=mle.c_profile_support_low;
        if (std::isfinite(mle.c_profile_support_high))
            c_profile_support_high_by_cell[bc]=mle.c_profile_support_high;
        if (std::isfinite(mle.r_profile_support_low))
            r_profile_support_low_by_cell[bc]=mle.r_profile_support_low;
        if (std::isfinite(mle.r_profile_support_high))
            r_profile_support_high_by_cell[bc]=mle.r_profile_support_high;
        if (std::isfinite(mle.r_profile_support_width))
            r_profile_support_width_by_cell[bc]=mle.r_profile_support_width;

        if (!mle.fit_success){
            n_failed++;
            // Preserve a status-bearing diagnostic row for every compiled cell,
            // including genuine numerical/data failures that cannot produce a
            // contamination estimate. The launcher copies this sidecar back to
            // durable storage, so failure is explicit rather than a silent drop.
            optimizer_status_by_cell[bc] = string("mle:") + mle.optimizer_status +
                string(";regularized:") + reg_attempt_status[i];
            informative_weight_by_cell[bc] = mle.informative_weight;
            gradient_norm_by_cell[bc] = mle.gradient_norm;
            interval_method_by_cell[bc] = mle.interval_method;
            prior_group_by_cell[bc] = mle.prior_group;
            if (std::isfinite(mle.ll)) mle_log_likelihood[bc] = mle.ll;
            if (std::isfinite(mle.diagnostic_fixed_r_fallback_c)){
                diagnostic_fixed_r_fallback_c_by_cell[bc] =
                    mle.diagnostic_fixed_r_fallback_c;
            }
            if (std::isfinite(mle.diagnostic_fixed_r_fallback_ll)){
                diagnostic_fixed_r_fallback_ll_by_cell[bc] =
                    mle.diagnostic_fixed_r_fallback_ll;
            }
            if (std::isfinite(mle.c_profile_candidate_gap))
                c_profile_candidate_gap_by_cell[bc]=mle.c_profile_candidate_gap;
            if (std::isfinite(mle.c_profile_neighbor_gain))
                c_profile_neighbor_gain_by_cell[bc]=mle.c_profile_neighbor_gain;
            if (std::isfinite(mle.c_profile_global_gain))
                c_profile_global_gain_by_cell[bc]=mle.c_profile_global_gain;
            if (std::isfinite(mle.c_profile_likelihood_span))
                c_profile_likelihood_span_by_cell[bc]=mle.c_profile_likelihood_span;
            if (std::isfinite(mle.c_profile_boundary_slope))
                c_profile_boundary_slope_by_cell[bc]=mle.c_profile_boundary_slope;
            if (std::isfinite(mle.c_profile_support_width))
                c_profile_support_width_by_cell[bc]=mle.c_profile_support_width;
            diagnostic_fallback_status_by_cell[bc] =
                mle.diagnostic_fallback_status;
            int flags = 0;
            if (mle.hit_c_lower) flags |= 1;
            if (mle.hit_c_upper) flags |= 2;
            if (mle.hit_r_lower) flags |= 4;
            if (mle.hit_r_upper) flags |= 8;
            fit_boundary_flags[bc] = flags;
            continue;
        }

        contam_rate_mle[bc] = mle.c;
        contam_rate_regularized[bc] = selected.c;
        contam_rate[bc] = selected.c;
        contam_rate_ll[bc] = selected.ll;
        mle_log_likelihood[bc] = mle.ll;
        regularized_log_likelihood[bc] = selected.ll;
        informative_weight_by_cell[bc] = mle.informative_weight;
        gradient_norm_by_cell[bc] = mle.gradient_norm;

        double selected_c_se = selected.c_se;
        if (std::isfinite(selected_c_se)) selected_c_se *= dispersion_scale;
        contam_rate_se[bc] = selected_c_se;

        contam_ci_low[bc] = mle.c_ci_low;
        contam_ci_high[bc] = mle.c_ci_high;
        ridge_width[bc] = mle.ridge_width;
        r_c_correlation[bc] = mle.r_c_correlation;
        interval_method_by_cell[bc] = mle.interval_method;
        optimizer_status_by_cell[bc] = string("mle:") + mle.optimizer_status +
            string(";regularized:") +
            (reg_fit_selected[i] ? reg.optimizer_status : reg_attempt_status[i]);
        if (std::isfinite(mle.diagnostic_fixed_r_fallback_c)){
            diagnostic_fixed_r_fallback_c_by_cell[bc] =
                mle.diagnostic_fixed_r_fallback_c;
        }
        if (std::isfinite(mle.diagnostic_fixed_r_fallback_ll)){
            diagnostic_fixed_r_fallback_ll_by_cell[bc] =
                mle.diagnostic_fixed_r_fallback_ll;
        }
        if (std::isfinite(mle.c_profile_candidate_gap))
            c_profile_candidate_gap_by_cell[bc]=mle.c_profile_candidate_gap;
        if (std::isfinite(mle.c_profile_neighbor_gain))
            c_profile_neighbor_gain_by_cell[bc]=mle.c_profile_neighbor_gain;
        if (std::isfinite(mle.c_profile_global_gain))
            c_profile_global_gain_by_cell[bc]=mle.c_profile_global_gain;
        if (std::isfinite(mle.c_profile_likelihood_span))
            c_profile_likelihood_span_by_cell[bc]=mle.c_profile_likelihood_span;
        if (std::isfinite(mle.c_profile_boundary_slope))
            c_profile_boundary_slope_by_cell[bc]=mle.c_profile_boundary_slope;
        if (std::isfinite(mle.c_profile_support_width))
            c_profile_support_width_by_cell[bc]=mle.c_profile_support_width;
        diagnostic_fallback_status_by_cell[bc] =
            mle.diagnostic_fallback_status;

        int flags = 0;
        if (mle.hit_c_lower) flags |= 1;
        if (mle.hit_c_upper) flags |= 2;
        if (mle.hit_r_lower) flags |= 4;
        if (mle.hit_r_upper) flags |= 8;
        fit_boundary_flags[bc] = flags;

        prior_mean_by_cell[bc] = reg_prior_mean[i];
        prior_var_by_cell[bc] = reg_prior_var[i];
        prior_group_by_cell[bc] = reg_prior_group[i];
        prior_displacement[bc] = selected.c - mle.c;

        // Allele-ratio reporting is gated on the identifiability of the free
        // (r,c) profile, which is assessed on the unregularized MLE fit and is
        // the value published as r_profile_identifiable.  The regularized fit
        // places a prior on c, not on r, so it cannot add information about r:
        // if the MLE r profile is flat, the regularized r is equally
        // unsupported.  Publishing selected.r on the strength of
        // selected.has_allele_ratio alone let a prior manufacture an allele
        // ratio for a cell whose own diagnostics declare r unidentifiable,
        // which is exactly the surrogate the free-model contract forbids.
        // Model classification below is deliberately left on
        // selected.has_allele_ratio: whether a cell was fit with the three
        // component model is a statement about the model used, not about
        // whether r is resolvable, and the two must not be conflated.
        const bool r_reportable =
            mle.r_profile_identifiable && selected.has_allele_ratio;

        if (mle.has_allele_ratio && mle.r_profile_identifiable){
            allele_ratio_mle[bc] = mle.r;
        }
        if (r_reportable){
            allele_ratio_regularized[bc] = selected.r;
            allele_ratio[bc] = selected.r;
            double selected_r_se = selected.r_se;
            if (std::isfinite(selected_r_se)) selected_r_se *= dispersion_scale;
            allele_ratio_se[bc] = selected_r_se;
            allele_ratio_ci_low[bc] = selected.r_ci_low;
            allele_ratio_ci_high[bc] = selected.r_ci_high;
            selected_r.push_back(selected.r);
        }
        if (selected.has_allele_ratio){
            n_three_comp++;
        } else {
            n_two_comp++;
        }
        if (selected.bfgs_fallback) n_fallback++;
        selected_c.push_back(selected.c);
    }

    if (heterotypic_contam_prior_valid){
        contam_cell_prior = heterotypic_contam_prior.first;
        contam_cell_prior_var = heterotypic_contam_prior.second;
    } else if (twocomp_contam_prior_valid){
        contam_cell_prior = twocomp_contam_prior.first;
        contam_cell_prior_var = twocomp_contam_prior.second;
    }

    if (!selected_c.empty()){
        pair<double,double> cv = contamination_prior_safe_moments(selected_c);
        fprintf(stderr, "Per-cell contamination estimates: mean=%f sd=%f phi=%f\n",
            cv.first, cv.second > 0.0 ? sqrt(cv.second) : 0.0, dispersion_phi);
    }
    fprintf(stderr, "  Three-component cells: %d  Two-component cells: %d  Failed: %d\n",
        n_three_comp, n_two_comp, n_failed);
    if (n_fallback > 0){
        fprintf(stderr, "  Heterotypic profile fits without usable BFGS candidates: %d\n", n_fallback);
    }
    if (!selected_r.empty()){
        pair<double,double> rv = contamination_prior_safe_moments(selected_r);
        fprintf(stderr, "  Allele ratio mean=%f sd=%f\n",
            rv.first, rv.second > 0.0 ? sqrt(rv.second) : 0.0);
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
/**
 * Per-cell data compiler for the individual ambient-profile fit.
 *
 * Two row blocks are emitted.
 *
 * ASSIGNMENT-KEYED BLOCK, historical and always present. Rows keyed on the
 * cell's own assigned identity: the three marginal bins of a singlet's donor,
 * or the nine pairwise bins of a fusion's declared pair. Each physical read
 * falls in exactly one of these rows.
 *
 * CANDIDATE-KEYED BLOCK, O53 step 2, emitted only when candidate_keyed_rows is
 * set. The marginal bin of every ambient candidate in idx2samp, at every
 * dosage present for that candidate in this cell. These are the rows that
 * separate donors from one another; the assignment-keyed block cannot, because
 * for a candidate that is not this cell's receiver its design column collapses
 * to a genome-wide average that depends on site allele frequency and not on the
 * donor. Measured across 40 benchmark units, the assignment-keyed block alone
 * is rank deficient by 13 to 16 of 23 free directions on every unit, while
 * adding this block reaches full rank on every unit.
 *
 * WEIGHT BUDGET. A read enters one marginal bin per candidate, so the
 * candidate-keyed block replicates each read once per candidate. The block is
 * therefore divided by the candidate count, and the two blocks share a single
 * unit of weight per read via candidate_keyed_split, so total evidence per read
 * is unchanged. See the header for why the second correction is not optional.
 */
void contamFinder3::compile_amb_prof_dat(bool solve_for_c, 
    bool use_global_c,
    vector<vector<double> >& mixfracs,
    vector<double>& weights,
    vector<double>& n,
    vector<double>& k,
    vector<double>& p_e,
    vector<double>& c){
    
    // O53 step 2 setup. The candidate set for the marginal block is exactly
    // idx2samp, the ambient candidate roster the mixture columns are built
    // over, so the block cannot introduce a column the solver does not model.
    // A block whose weight budget is zero must not be emitted at all. Pushing
    // zero-weight observations is not the same as omitting them: the optimizer
    // carries them into normalization and they poison the fit rather than
    // vanishing. Skipping them makes split = 1.0 an exact no-op, which is the
    // control that proves the weight plumbing is correct.
    const bool ck_enabled = candidate_keyed_rows && !idx2samp.empty() &&
        candidate_keyed_split < 1.0;
    set<int> ck_candidates;
    double ck_assigned_scale = 1.0;
    double ck_marginal_scale = 0.0;
    if (ck_enabled){
        ck_candidates.insert(idx2samp.begin(), idx2samp.end());
        ck_assigned_scale = candidate_keyed_split;
        ck_marginal_scale = (1.0 - candidate_keyed_split) /
            (double)ck_candidates.size();
    }

    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != 
        assn.end(); ++a){

        // Cross-fitted nuisance profile: held-out receivers are still scored
        // later by est_contam_cells(), but none of their CK rows may train the
        // ambient profile used to score them.
        if (profile_holdout_barcodes.count(a->first) > 0){
            continue;
        }
        
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

        // Historical rows keep the full weight unless the candidate-keyed
        // block is present, in which case the two blocks share it. A zero share
        // here would mean emitting zero-weight assignment rows, so the split is
        // floored to keep that block real; split = 0 is expressed by making the
        // candidate-keyed block dominant rather than by deleting these rows.
        const double weight_assigned = weight * ck_assigned_scale;

        double this_c;
        if (!solve_for_c){
            if (use_global_c){
                this_c = contam_cell_prior;
            }
            else{
                auto rate_it = contam_rate.find(a->first);
                if (rate_it == contam_rate.end()){
                    // A failed per-cell fit has no usable c. Keep the cell in
                    // assignment/diagnostic state, but omit it from this
                    // profile calculation rather than inventing c=0.
                    continue;
                }
                this_c = rate_it->second;
            }
        }

        // ---- O53 step 2: candidate-keyed marginal block ------------------
        //
        // Emitted for every cell, whichever assignment-keyed branch below
        // produced that cell's historical rows. No tetraploid-aware filter is
        // applied, matching the composition-override branch, where the filter
        // is required only for a genuine pair row and a marginal row has
        // second key -1.
        //
        // For a singlet cell the candidate that equals its own identity
        // reproduces rows the assignment-keyed block already emitted. That is
        // intentional and is not double counting: the weight budget above is
        // shared across both blocks, so total evidence per read is one either
        // way, and the self-conditioned bins are the most informative rows for
        // the endogenous term.
        auto emit_candidate_keyed_rows = [&](){
            if (!ck_enabled || !(ck_marginal_scale > 0.0)){
                return;
            }
            auto cell_rows = indv_allelecounts.find(a->first);
            if (cell_rows == indv_allelecounts.end()){
                return;
            }
            for (auto ac1 = cell_rows->second.begin();
                ac1 != cell_rows->second.end(); ++ac1){

                if (ck_candidates.find(ac1->first.first) == ck_candidates.end()){
                    continue;
                }
                auto ac2 = ac1->second.find(make_pair(-1, -1));
                if (ac2 == ac1->second.end()){
                    continue;
                }
                double ref = ac2->second.first;
                double alt = ac2->second.second;
                if (ref + alt <= 0){
                    continue;
                }
                auto cf_it = expfracs.find(ac1->first);
                if (cf_it == expfracs.end()){
                    continue;
                }
                map<int, float>& cf = cf_it->second;

                // Endogenous expectation at sites conditioned on this
                // candidate's dosage, from the same conditional-fraction table
                // the historical rows already read. A fusion mixes its two
                // halves by the receiver's allele ratio, the same quantity and
                // the same guards the pairwise branch applies.
                double endo_raw;
                if (is_comb){
                    double cell_r = 0.5;
                    if (r_feedback_enabled && comb.first != comb.second &&
                        allele_ratio.count(a->first) > 0){
                        cell_r = allele_ratio[a->first];
                        if (cell_r < 0.01) cell_r = 0.01;
                        if (cell_r > 0.99) cell_r = 0.99;
                    }
                    double ef1 = 0.5;
                    double ef2 = 0.5;
                    if (cf.count(comb.first) > 0){
                        ef1 = cf[comb.first];
                    }
                    if (cf.count(comb.second) > 0){
                        ef2 = cf[comb.second];
                    }
                    endo_raw = cell_r * ef1 + (1.0 - cell_r) * ef2;
                }
                else{
                    endo_raw = 0.5;
                    if (cf.count(a->second) > 0){
                        endo_raw = cf[a->second];
                    }
                }
                double expected = adjust_p_err(endo_raw, e_r, e_a);

                double row_weight = weight * ck_marginal_scale;
                if (!(row_weight > 0.0) || !std::isfinite(row_weight)){
                    continue;
                }
                n.push_back(ref + alt);
                k.push_back(alt);
                weights.push_back(row_weight);
                if (!solve_for_c){
                    c.push_back(this_c);
                }
                p_e.push_back(expected);

                // Mixture columns are the candidate-conditioned expectation for
                // every modelled ambient source, matching the singlet branch of
                // the historical block.
                vector<double> mixfrac_row;
                for (int i = 0; i < (int)idx2samp.size(); ++i){
                    int samp = idx2samp[i];
                    double ef_val = 0.5;
                    if (cf.count(samp) > 0){
                        ef_val = cf[samp];
                    }
                    mixfrac_row.push_back(ef_val);
                }
                if (inter_species){
                    mixfrac_row.push_back(adjust_p_err(0.0, e_r, e_a));
                }
                mixfracs.push_back(mixfrac_row);
            }
        };

        // Symmetric to the marginal-block guard above. A zero assignment-keyed
        // share must drop those rows rather than push them at weight zero,
        // because a zero-weight observation is carried into the optimizer's
        // normalization instead of vanishing. split = 0 therefore means the
        // candidate-keyed block alone, which is a legitimate configuration and
        // is the lower endpoint of the split sweep.
        if (!(weight_assigned > 0.0)){
            emit_candidate_keyed_rows();
            continue;
        }

        if (has_composition_override(a->first)){
            const map<int, double>& comp = cell_composition_overrides[a->first];
            WeightedCompositionRowSelection selection =
                select_weighted_composition_rows(
                    comp,indv_allelecounts[a->first],expfracs);

            for (auto row = selection.rows.begin(); row != selection.rows.end(); ++row){
                auto ac1 = indv_allelecounts[a->first].find(row->first);
                if (ac1 == indv_allelecounts[a->first].end()) continue;
                auto ac2 = ac1->second.find(row->second);
                if (ac2 == ac1->second.end()) continue;

                double ref = ac2->second.first;
                double alt = ac2->second.second;
                if (ref + alt <= 0) continue;

                double expected_raw = 0.0;
                if (!composition_expected_from_row_complete(
                        comp,ac1->first,ac2->first,expected_raw)){
                    continue;
                }
                double expected = adjust_p_err(expected_raw,e_r,e_a);

                if (tetraploid_aware && comp.size() > 1){
                    double this_p_c=0.0;
                    if (amb_mu_available && amb_mu.count(ac1->first)>0 &&
                        amb_mu[ac1->first].count(ac2->first)>0){
                        this_p_c=amb_mu[ac1->first][ac2->first];
                    }
                    const bool filter_required=ac2->first.first>=0;
                    if (!contamination_observation_passes_model_filter(
                            filter_required,amb_mu_available,expected,this_p_c,
                            min_signal_gap,ac1->first.second,
                            ac2->first.second)) continue;
                }

                n.push_back(ref + alt);
                k.push_back(alt);
                weights.push_back(weight_assigned);
                if (!solve_for_c) c.push_back(this_c);
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
                                ac1->first.second / 2.0,e_r,e_a));
                        } else if (ef_all_avg && ac2->first.first == samp){
                            mixfrac_row.push_back(adjust_p_err(
                                ac2->first.second / 2.0,e_r,e_a));
                        } else {
                            double ef1 = 0.5, ef2 = 0.5;
                            if (expfracs.count(ac1->first) > 0 &&
                                expfracs[ac1->first].count(samp) > 0){
                                ef1 = expfracs[ac1->first][samp];
                            }
                            if (expfracs.count(ac2->first) > 0 &&
                                expfracs[ac2->first].count(samp) > 0){
                                ef2 = expfracs[ac2->first][samp];
                            }
                            mixfrac_row.push_back(0.5*ef1+0.5*ef2);
                        }
                    }
                }
                if (inter_species){
                    mixfrac_row.push_back(adjust_p_err(0.0,e_r,e_a));
                }
                mixfracs.push_back(mixfrac_row);
            }
            emit_candidate_keyed_rows();
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
                        
                        // R-feedback updates only the receiver's endogenous
                        // parental expectation. Ambient-source genotype columns
                        // remain receiver-independent. When feedback is off or no
                        // allele_ratio exists, cell_r remains the historical 0.5.
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
                            double this_p_c=0.0;
                            if (amb_mu_available && amb_mu.count(ac1->first)>0 &&
                                amb_mu[ac1->first].count(ac2->first)>0){
                                this_p_c=amb_mu[ac1->first][ac2->first];
                            }
                            if (!contamination_observation_passes_model_filter(
                                    true,amb_mu_available,expected,this_p_c,
                                    min_signal_gap,ac1->first.second,
                                    ac2->first.second)) continue;
                        }
                        
                        double ref = ac2->second.first;
                        double alt = ac2->second.second;
                        
                        n.push_back(ref+alt);
                        k.push_back(alt);
                        weights.push_back(weight_assigned);
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
                                    // Ambient-source genotype expectations must not depend
                                    // on the receiver cell's endogenous parental RNA ratio.
                                    // The compressed pair-category representation retains the
                                    // symmetric approximation until the site-level artifact is
                                    // available.
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

        emit_candidate_keyed_rows();
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
    // Empty-matrix guard. A zero-row component matrix (e.g. assn emptied by a
    // reclassification pass on a degenerate unit) makes add_mixcomp return false,
    // leaving the solver with nmixcomp==0 and an empty mixcompfracs_sparse;
    // calling solve() then indexes that empty structure and segfaults. Skip the
    // solve cleanly, leave contam_prof and amb_mu untouched, and return a value
    // the iteration accept/reject logic treats as no improvement. (Accept/reject
    // keys on compute_ll(), not this return; -DBL_MAX is an unambiguous no-op.)
    if (mixfracs.empty()){
        fprintf(stderr, "WARNING: empty ambient mixture component matrix; "
            "skipping ambient profile solve (keeping current profile)\n");
        return -DBL_MAX;
    }
    if (!solver.add_mixcomp(mixfracs) || !solver.add_mixcomp_fracs(startprops)){
        fprintf(stderr, "WARNING: ambient mixture component setup failed "
            "(empty or mismatched matrix); skipping ambient profile solve "
            "(keeping current profile)\n");
        return -DBL_MAX;
    }
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
        // This is the only path that compares starting points; record that it
        // ran so a later refinement solve cannot be mistaken for a multistart.
        multistart_attempted = true;
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
        
        // The caller can request an exact total number of profile starts.  The
        // first start is the supplied/warm-start profile and the second is
        // uniform.  Remaining starts are deterministic Dirichlet(1) draws so
        // repeated benchmark runs are byte-for-byte reproducible with respect
        // to optimizer initialization.  A negative value preserves the legacy
        // candidate-count-scaled behavior.
        const bool use_uniform_start = (profile_total_starts < 0 || profile_total_starts >= 2);
        if (use_uniform_start){
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
        }

        int ntrials = profile_total_starts >= 0
            ? std::max(0, profile_total_starts - 2)
            : (int)contam_prof.size() * n_mixprop_trials;
        // Starts actually configured for this solve: the supplied/warm start,
        // the uniform start when enabled, and the deterministic Dirichlet
        // trials. Recorded so successful/configured is interpretable without
        // reproducing this arithmetic in the driver.
        multistart_configured_starts = 1 + (use_uniform_start ? 1 : 0)
            + std::max(0, ntrials);
        std::mt19937 profile_rng(0x51A7E123u
            + (unsigned int)startprops.size() * 131u
            + (unsigned int)mixfracs.size() * 17u);
        std::gamma_distribution<double> profile_gamma(1.0, 1.0);
        for (int trial_idx = 0; trial_idx < ntrials; ++trial_idx){
            try {
                vector<double> trialprops(startprops.size(), 0.0);
                double trialsum = 0.0;
                for (int j = 0; j < (int)trialprops.size(); ++j){
                    trialprops[j] = profile_gamma(profile_rng);
                    trialsum += trialprops[j];
                }
                if (!(trialsum > 0.0) || !std::isfinite(trialsum)){
                    continue;
                }
                for (int j = 0; j < (int)trialprops.size(); ++j){
                    trialprops[j] /= trialsum;
                }
                solver.add_mixcomp_fracs(trialprops);
                solver.solve();
                if (std::isfinite(solver.log_likelihood) && std::isfinite(solver.results[0])){
                    any_solve_succeeded = true;
                    size_t this_idx = lls.size();
                    lls.push_back(solver.log_likelihood);
                    mcs.push_back(solver.results_mixcomp);
                    cs.push_back(solver.results[0]);
                    if (solver.log_likelihood > maxll){
                        maxll = solver.log_likelihood;
                        maxidx = (int)this_idx;
                    }
                }
            } catch (...) {
                // ignore, try next deterministic starting condition
            }
        }

        if (any_solve_succeeded){
            solver.log_likelihood = maxll;
            solver.results[0] = cs[maxidx];
            solver.results_mixcomp = mcs[maxidx];

            profile_successful_starts = (int)lls.size();
            profile_best_ll = maxll;
            profile_second_best_ll = -DBL_MAX;
            for (int i = 0; i < (int)lls.size(); ++i){
                if (i != maxidx && lls[i] > profile_second_best_ll){
                    profile_second_best_ll = lls[i];
                }
            }
            const double near_tol = std::max(0.1, 1e-6 * std::max(1.0, std::fabs(maxll)));
            profile_near_optimal_count = 0;
            profile_near_optimal_l1_spread = 0.0;
            for (int i = 0; i < (int)lls.size(); ++i){
                if (maxll - lls[i] <= near_tol){
                    profile_near_optimal_count++;
                    double l1 = 0.0;
                    if (mcs[i].size() == mcs[maxidx].size()){
                        for (int j = 0; j < (int)mcs[i].size(); ++j){
                            l1 += std::fabs(mcs[i][j] - mcs[maxidx][j]);
                        }
                    }
                    if (l1 > profile_near_optimal_l1_spread){
                        profile_near_optimal_l1_spread = l1;
                    }
                }
            }

            // Copy the multistart comparison into fields the later refinement
            // solve does not touch. fit() calls update_amb_prof_mixture() again
            // with solve_for_c == false, and that path hardcodes the profile_*
            // fields to a single-solve tuple, so without this copy the
            // multistart result never reaches the serialized diagnostics.
            multistart_successful_starts = profile_successful_starts;
            multistart_best_ll = profile_best_ll;
            multistart_second_best_ll = profile_second_best_ll;
            multistart_near_optimal_count = profile_near_optimal_count;
            multistart_near_optimal_l1_spread = profile_near_optimal_l1_spread;

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
            // Attempted but nothing finite came back. Record zero successes
            // explicitly so this is distinguishable from "multistart not run".
            multistart_successful_starts = 0;
            multistart_near_optimal_count = 0;
            multistart_best_ll = -DBL_MAX;
            multistart_second_best_ll = -DBL_MAX;
            multistart_near_optimal_l1_spread = std::numeric_limits<double>::quiet_NaN();
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
            profile_successful_starts = 1;
            profile_near_optimal_count = 1;
            profile_best_ll = maxll;
            profile_second_best_ll = -DBL_MAX;
            profile_near_optimal_l1_spread = 0.0;
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
                    // Defensive neutral fallback only; production strict-condf validation rejects any required missing lookup
                    double ef_val = 0.5;
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
                                double ef1 = 0.5;
                                double ef2 = 0.5;
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
    // Never leave stale concentration parameters behind after a skipped or
    // failed bootstrap attempt.
    dirichlet_params.clear();
    species_contam_prof_conc.clear();

    if (n_boots <= 0){
        fprintf(stderr, "Bootstrap disabled (n_boots <= 0); skipping concentration fit\n");
        return;
    }

    // Candidate-keyed rows are multiple dependent representations of the same
    // cell/read evidence. Row-wise resampling breaks those blocks apart, so its
    // concentration parameters are not calibrated even when the replicate
    // objective is weighted correctly. Suppress them until cell-clustered
    // resampling is implemented. split == 1.0 emits no candidate-keyed block
    // and remains a valid control.
    const bool candidate_keyed_bootstrap_active =
        !bulk_mode && candidate_keyed_rows && !idx2samp.empty() &&
        candidate_keyed_split < 1.0;
    if (candidate_keyed_bootstrap_active){
        fprintf(stderr,
            "WARNING: candidate-keyed ambient-profile rows are active; "
            "row-wise bootstrap does not preserve cell-level dependence. "
            "Profile concentration parameters will be omitted until clustered "
            "resampling is implemented.\n");
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

    const size_t n_compiled_obs = mixfracs.size();
    if (n_compiled_obs == 0 ||
        weights.size() != n_compiled_obs ||
        n.size() != n_compiled_obs ||
        k.size() != n_compiled_obs ||
        p_e.size() != n_compiled_obs ||
        c.size() != n_compiled_obs){
        fprintf(stderr,
            "Bootstrap samples complete: 0 successful of %d requested\n",
            n_boots);
        fprintf(stderr,
            "WARNING: compiled bootstrap data are empty or misaligned; "
            "profile concentration parameters will be omitted\n");
        return;
    }

    // fit_dirichlet() historically substitutes point proportions when its
    // optimizer fails. Bootstrap uncertainty must never silently become a copy
    // of the point estimate, so use the same Dirichlet objective here with
    // explicit success/failure reporting and no fallback values.
    extern double ll_dirichlet(const vector<double>&,
        const map<string, double>&, const map<string, int>&);
    extern void dll_dirichlet(const vector<double>&,
        const map<string, double>&, const map<string, int>&, vector<double>&);

    auto fit_bootstrap_dirichlet = [&](
        vector<double>& mle_fracs_local,
        vector<vector<double> >& bootstrap_props,
        vector<double>& fitted_params) -> bool {

        fitted_params.clear();
        const size_t n_components = mle_fracs_local.size();
        if (n_components == 0 || bootstrap_props.size() != n_components){
            return false;
        }

        const size_t n_replicates = bootstrap_props[0].size();
        if (n_replicates == 0){
            return false;
        }

        for (size_t j = 0; j < n_components; ++j){
            if (!std::isfinite(mle_fracs_local[j]) ||
                !(mle_fracs_local[j] > 0.0) ||
                bootstrap_props[j].size() != n_replicates){
                return false;
            }
            for (size_t b = 0; b < n_replicates; ++b){
                // The Dirichlet likelihood contains log(f_j), so boundary or
                // non-finite bootstrap proportions are not usable inputs.
                if (!std::isfinite(bootstrap_props[j][b]) ||
                    !(bootstrap_props[j][b] > 0.0)){
                    return false;
                }
            }
        }

        vector<double> dir_init(mle_fracs_local);
        optimML::multivar_ml_solver dirsolver(
            dir_init, ll_dirichlet, dll_dirichlet);
        dirsolver.set_silent(true);

        for (size_t j = 0; j < n_components; ++j){
            dirsolver.constrain_pos((int)j);
            string key = string("f_") + std::to_string(j);
            if (!dirsolver.add_data(key, bootstrap_props[j])){
                return false;
            }
        }

        try {
            if (!dirsolver.solve() ||
                !std::isfinite(dirsolver.log_likelihood) ||
                dirsolver.results.size() != n_components){
                return false;
            }
        } catch (...) {
            return false;
        }

        for (size_t j = 0; j < n_components; ++j){
            double alpha = dirsolver.results[j];
            if (!std::isfinite(alpha) || !(alpha > 0.0)){
                fitted_params.clear();
                return false;
            }
            fitted_params.push_back(alpha);
        }
        return true;
    };

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

        // Species aggregation indexes mixfracs by idx2samp position before the
        // optimizer can validate the matrix, so guard those dimensions here.
        for (size_t obs = 0; obs < mixfracs.size(); ++obs){
            if (mixfracs[obs].size() < idx2samp.size()){
                fprintf(stderr,
                    "Bootstrap samples complete: 0 successful of %d requested\n",
                    n_boots);
                fprintf(stderr,
                    "WARNING: compiled species bootstrap matrix has too few "
                    "columns; profile concentration parameters will be omitted\n");
                return;
            }
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
            fprintf(stderr,
                "Bootstrap samples complete: 0 successful of %d requested\n",
                n_boots);
            fprintf(stderr,
                "WARNING: no active species for bootstrap; profile "
                "concentration parameters will be omitted\n");
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

            const size_t replicate_n = mf_boot.size();
            if (replicate_n == 0 ||
                w_boot.size() != replicate_n ||
                n_boot.size() != replicate_n ||
                k_boot.size() != replicate_n ||
                pe_boot.size() != replicate_n ||
                c_boot.size() != replicate_n){
                continue;
            }

            bool weights_valid = true;
            for (size_t x = 0; x < w_boot.size(); ++x){
                if (!std::isfinite(w_boot[x]) || w_boot[x] < 0.0){
                    weights_valid = false;
                    break;
                }
            }
            if (!weights_valid){
                continue;
            }

            vector<double> params;
            optimML::multivar_ml_solver solver(params, ll_amb_prof_mixture,
                dll_amb_prof_mixture);
            solver.set_silent(true);
            if (!solver.add_mixcomp(mf_boot)){
                continue;
            }
            if (!solver.add_mixcomp_fracs(startprops)){
                continue;
            }
            if (!solver.add_data("n", n_boot)){
                continue;
            }
            if (!solver.add_data("k", k_boot)){
                continue;
            }
            if (!solver.add_data("p_e", pe_boot)){
                continue;
            }
            if (!solver.add_data("c", c_boot)){
                continue;
            }
            if (!solver.add_weights(w_boot)){
                continue;
            }

            try {
                if (!solver.solve() ||
                    !std::isfinite(solver.log_likelihood) ||
                    (int)solver.results_mixcomp.size() != n_active_sp){
                    continue;
                }

                bool result_valid = true;
                double result_sum = 0.0;
                for (int x = 0; x < n_active_sp; ++x){
                    double value = solver.results_mixcomp[x];
                    if (!std::isfinite(value) || value < 0.0 || value > 1.0){
                        result_valid = false;
                        break;
                    }
                    result_sum += value;
                }
                if (!result_valid || !std::isfinite(result_sum) ||
                    fabs(result_sum - 1.0) > 1e-6){
                    continue;
                }
                boot_results[b] = solver.results_mixcomp;
            } catch (...) {
                // Skip failed bootstrap sample.
            }
        }

        int n_boot_success = 0;
        for (int b = 0; b < n_boots; ++b){
            if ((int)boot_results[b].size() != n_active_sp) continue;
            ++n_boot_success;
            for (int x = 0; x < n_active_sp; ++x){
                dirprops[x].push_back(boot_results[b][x]);
            }
        }
        fprintf(stderr,
            "Bootstrap samples complete: %d successful of %d requested\n",
            n_boot_success, n_boots);

        if (n_boot_success == 0){
            fprintf(stderr,
                "WARNING: no complete bootstrap replicates succeeded; "
                "profile concentration parameters will be omitted\n");
            return;
        }

        // Fit Dirichlet at active-species level without substituting point
        // proportions when the concentration fit fails.
        vector<double> dirichlet_soln;
        if (!fit_bootstrap_dirichlet(mle_fracs, dirprops, dirichlet_soln)){
            fprintf(stderr,
                "WARNING: Dirichlet concentration fit failed; profile "
                "concentration parameters will be omitted\n");
            return;
        }

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
    int n_obs = (int)mixfracs.size();
    int n_props = (int)startprops.size();
    if (n_props == 0){
        fprintf(stderr,
            "Bootstrap samples complete: 0 successful of %d requested\n",
            n_boots);
        fprintf(stderr,
            "WARNING: bootstrap mixture has no components; profile "
            "concentration parameters will be omitted\n");
        return;
    }
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

        const size_t replicate_n = mixfracs_boot.size();
        if (replicate_n == 0 ||
            weights_boot.size() != replicate_n ||
            n_boot.size() != replicate_n ||
            k_boot.size() != replicate_n ||
            p_e_boot.size() != replicate_n ||
            c_boot.size() != replicate_n){
            continue;
        }

        bool weights_valid = true;
        for (size_t x = 0; x < weights_boot.size(); ++x){
            if (!std::isfinite(weights_boot[x]) || weights_boot[x] < 0.0){
                weights_valid = false;
                break;
            }
        }
        if (!weights_valid){
            continue;
        }

        vector<double> params;
        optimML::multivar_ml_solver solver(params, ll_amb_prof_mixture,
            dll_amb_prof_mixture);
        solver.set_silent(true);
        if (!solver.add_mixcomp(mixfracs_boot)){
            continue;
        }
        if (!solver.add_mixcomp_fracs(startprops)){
            continue;
        }
        if (!solver.add_data("n", n_boot)){
            continue;
        }
        if (!solver.add_data("k", k_boot)){
            continue;
        }
        if (!solver.add_data("p_e", p_e_boot)){
            continue;
        }
        if (!solver.add_data("c", c_boot)){
            continue;
        }
        if (!solver.add_weights(weights_boot)){
            continue;
        }

        try {
            if (!solver.solve() ||
                !std::isfinite(solver.log_likelihood) ||
                (int)solver.results_mixcomp.size() != n_props){
                continue;
            }

            bool result_valid = true;
            double result_sum = 0.0;
            for (int x = 0; x < n_props; ++x){
                double value = solver.results_mixcomp[x];
                if (!std::isfinite(value) || value < 0.0 || value > 1.0){
                    result_valid = false;
                    break;
                }
                result_sum += value;
            }
            if (!result_valid || !std::isfinite(result_sum) ||
                fabs(result_sum - 1.0) > 1e-6){
                continue;
            }
            boot_results[b] = solver.results_mixcomp;
        } catch (...) {
            // Skip failed bootstrap sample.
        }
    }

    int n_boot_success = 0;
    for (int b = 0; b < n_boots; ++b){
        if ((int)boot_results[b].size() != n_props) continue;
        ++n_boot_success;
        for (int x = 0; x < n_props; ++x){
            dirprops[x].push_back(boot_results[b][x]);
        }
    }
    fprintf(stderr,
        "Bootstrap samples complete: %d successful of %d requested\n",
        n_boot_success, n_boots);

    if (n_boot_success == 0){
        fprintf(stderr,
            "WARNING: no complete bootstrap replicates succeeded; profile "
            "concentration parameters will be omitted\n");
        return;
    }

    vector<double> dirichlet_soln;
    if (!fit_bootstrap_dirichlet(mle_fracs, dirprops, dirichlet_soln)){
        fprintf(stderr,
            "WARNING: Dirichlet concentration fit failed; profile "
            "concentration parameters will be omitted\n");
        return;
    }

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

        // A failed per-cell fit has no accepted c to compare or update.
        // Reclassification is optional, so keep its assignment/diagnostics
        // untouched instead of creating c=0 through map operator[].
        if (contam_rate.find(a->first) == contam_rate.end()){
            continue;
        }
        
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
            reassign_allowed_ids, *cell_allowed_ids2_ptr, doub_rate_table, e_r, e_a, priorweights_ptr,
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
                // Branch (A): no candidate identity beats the null margin.
                // Reclassification is optional, so preserve the previously
                // accepted assignment and fit state.
            }
        }
        else{
            // Branch (B): populate_llr_table could not produce a usable
            // comparison. Preserve the previously accepted state.
        }
    }
    
    fprintf(stderr, "  %d cells reassigned\n", n_reassigned);
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

void contamFinder3::write_cell_source_profile(const string& filename,
    const vector<string>& samples, const string& libname,
    bool cellranger, bool seurat, bool underscore) const{

    ofstream out(filename.c_str());
    if (!out.is_open()){
        fprintf(stderr, "ERROR: could not open cell source-profile output: %s\n",
            filename.c_str());
        return;
    }
    out << std::setprecision(17);
    out << "barcode\tidentity\tsource_index\tsource_label"
        << "\tglobal_profile_mass\tscoring_profile_mass\tsource_excluded"
        << "\tsource_exclusion_strength\tprofile_variant"
        << "\tprofile_training_holdout\tprofile_holdout_basis"
        << "\tall_source_columns_retained\tglobal_profile_mass_sum"
        << "\tscoring_profile_mass_sum\tscoring_profile_status\n";

    vector<string> sample_copy = samples;
    for (robin_hood::unordered_map<unsigned long,int>::const_iterator ai = assn.begin();
        ai != assn.end(); ++ai){
        const unsigned long bc = ai->first;
        string bc_out = bc2str(bc);
        if (seurat && !libname.empty()) bc_out = libname + "_" + bc_out;
        else if (cellranger) bc_out += "-1";
        else if (!libname.empty()) bc_out += (underscore ? "_" : "-") + libname;

        const string ident = idx2name(ai->second, sample_copy);
        const set<int> exclusions = source_exclusions_for_cell(bc);
        const SourceExclusionProfileAssessment assessment =
            apply_source_exclusion_profile(contam_prof, exclusions,
                source_exclusion_strength);
        const bool held_out = profile_holdout_barcodes.count(bc) > 0;
        const bool all_columns_retained = source_exclusion_strength <= 1e-12;

        // contam_prof is the complete active source simplex. Iterate its keys,
        // not just positive values, so an explicit zero source column cannot
        // disappear from the scientific ledger.
        for (map<int,double>::const_iterator cp = contam_prof.begin();
            cp != contam_prof.end(); ++cp){
            const int source = cp->first;
            string label = "other_species";
            if (source >= 0 && source < (int)sample_copy.size()){
                label = sample_copy[source];
            }
            double scoring = 0.0;
            map<int,double>::const_iterator sp = assessment.scoring_profile.find(source);
            if (sp != assessment.scoring_profile.end()) scoring = sp->second;
            out << bc_out << '\t' << ident << '\t' << source << '\t' << label
                << '\t' << cp->second << '\t' << scoring
                << '\t' << ((source_exclusion_strength > 1e-12 &&
                    exclusions.count(source) > 0) ? 1 : 0)
                << '\t' << source_exclusion_strength
                << '\t' << profile_variant_label
                << '\t' << (held_out ? 1 : 0)
                << '\t' << profile_holdout_basis_label
                << '\t' << (all_columns_retained ? 1 : 0)
                << '\t' << assessment.global_profile_mass_sum
                << '\t' << assessment.scoring_profile_mass_sum
                << '\t' << assessment.status << '\n';
        }
    }
    out.close();
}

void contamFinder3::write_cell_fit_diagnostics(const string& filename,
    const vector<string>& samples, const string& libname,
    bool cellranger, bool seurat, bool underscore) const{

    ofstream out(filename.c_str());
    if (!out.is_open()){
        fprintf(stderr,"ERROR: could not open cell fit diagnostics: %s\n",filename.c_str());
        return;
    }
    out << std::setprecision(17);
    out << "barcode\tidentity\tis_heterotypic\tcontam_mle\tcontam_regularized"
        << "\tallele_ratio_mle\tallele_ratio_regularized\tcontam_se\tallele_ratio_se"
        << "\tcontam_ci_low\tcontam_ci_high\tallele_ratio_ci_low\tallele_ratio_ci_high"
        << "\tr_c_correlation\tridge_width\tprior_mean\tprior_variance\tprior_group"
        << "\tprior_displacement\tlog_likelihood_mle\tlog_likelihood_regularized"
        << "\tinformative_weight\tgradient_norm\tprior_training_eligible"
        << "\tprior_training_reason\tboundary_flags\tinterval_method\toptimizer_status"
        << "\tcompile_status\tdiagnostic_fixed_r_fallback_c"
        << "\tdiagnostic_fixed_r_fallback_ll\tdiagnostic_fallback_status"
        << "\tc_profile_candidate_gap\tc_profile_neighbor_gain"
        << "\tc_profile_global_gain\tc_profile_likelihood_span"
        << "\tc_profile_boundary_slope\tc_profile_support_low"
        << "\tc_profile_support_high\tc_profile_support_width"
        << "\tr_profile_optimization_succeeded\tr_profile_identifiable"
        << "\tr_profile_support_low\tr_profile_support_high"
        << "\tr_profile_support_width\tr_validation_status"
        << "\ttwo_dimensional_grid_starts_evaluated"
        << "\ttwo_dimensional_starts_attempted"
        << "\ttwo_dimensional_solver_successful_starts"
        << "\ttwo_dimensional_finite_in_domain_bfgs_candidates"
        << "\tprofile_optimization_attempted\tprofile_optimization_succeeded"
        << "\ttwo_dimensional_profile_validation_passed\tprofile_validation_status"
        << "\traw_candidate_r\traw_candidate_c\traw_candidate_log_likelihood"
        << "\tvalidated_candidate_r\tvalidated_candidate_c"
        << "\tvalidated_candidate_log_likelihood\tdispersion_phi"
        << "\tprofile_variant\tsource_exclusion_strength"
        << "\tfixed_r_enabled\tfixed_r_basis\tfixed_ambient_enabled"
        << "\tfixed_ambient_basis\ttruth_assisted_condition"
        << "\tambient_alt_expectation_global_weighted_mean"
        << "\tambient_alt_expectation_scoring_weighted_mean"
        << "\tambient_alt_expectation_weighted_mean_delta"
        << "\tambient_alt_expectation_weighted_rms_delta"
        << "\tambient_expectation_total_weight\tambient_expectation_status"
        << "\tambient_parent_axis_alpha_global\tambient_parent_axis_alpha_scoring"
        << "\tambient_parent_axis_alpha_delta\tambient_orthogonal_norm_global"
        << "\tambient_orthogonal_norm_scoring\tparent_axis_denominator"
        << "\tparent_axis_geometry_status"
        << "\tinformation_rr\tinformation_cc\tinformation_rc"
        << "\tconditional_information_c_given_r\tconditional_information_status"
        << "\texcluded_source_mass_raw_total\texcluded_parent_A_mass_raw"
        << "\texcluded_parent_B_mass_raw\teffective_removed_source_mass_total"
        << "\teffective_removed_parent_A_mass\teffective_removed_parent_B_mass"
        << "\tscoring_profile_renormalization_denominator"
        << "\tglobal_profile_mass_sum\tscoring_profile_mass_sum"
        << "\tsource_exclusion_status\n";

    const double nanv = std::numeric_limits<double>::quiet_NaN();
    vector<string> sample_copy = samples;
    auto profile_mass = [&](int idx)->double{
        auto it=contam_prof.find(idx);
        return it==contam_prof.end()?0.0:it->second;
    };

    for (robin_hood::unordered_map<unsigned long,int>::const_iterator ai = assn.begin();
        ai != assn.end(); ++ai){
        const unsigned long bc = ai->first;
        string bc_out = bc2str(bc);
        if (seurat && !libname.empty()) bc_out = libname + "_" + bc_out;
        else if (cellranger) bc_out += "-1";
        else if (!libname.empty()) bc_out += (underscore ? "_" : "-") + libname;

        string ident = idx2name(ai->second,sample_copy);
        bool hetero = false;
        int parent_a = -1, parent_b = -1;
        if (ai->second >= n_samples){
            pair<int,int> comb = idx_to_hap_comb(ai->second,n_samples);
            parent_a = comb.first; parent_b = comb.second;
            hetero = comb.first != comb.second;
        } else {
            parent_a = ai->second;
        }
        auto dget = [&](const robin_hood::unordered_map<unsigned long,double>& m)->double{
            auto it=m.find(bc); return it==m.end()?nanv:it->second;
        };
        auto iget = [&](const robin_hood::unordered_map<unsigned long,int>& m)->int{
            auto it=m.find(bc); return it==m.end()?0:it->second;
        };
        auto sget = [&](const map<unsigned long,string>& m,const string& def)->string{
            auto it=m.find(bc); return it==m.end()?def:it->second;
        };

        vector<double> n,k,p_e,p_c_scoring,p_c_global,p_A,p_B,obs_weights;
        bool compiled_hetero=false;
        string model_status;
        auto ci = cell_to_idx.find(bc);
        bool model_ok = ci != cell_to_idx.end() && compile_cell_model_data(
            bc, ci->second, n, k, p_e, p_c_scoring, p_c_global,
            p_A, p_B, obs_weights, compiled_hetero, model_status);

        double mean_global=nanv, mean_scoring=nanv, mean_delta=nanv, rms_delta=nanv;
        double expectation_weight=0.0;
        string expectation_status=model_ok?"ok":model_status;
        vector<double> geometry_weights;
        if (model_ok){
            double sg=0.0, ss=0.0, sd2=0.0;
            geometry_weights.reserve(n.size());
            for(size_t i=0;i<n.size();++i){
                // w_j is the exact observation weight supplied to the model.
                // Read depth n_j remains inside the binomial likelihood and must
                // not be multiplied into these separate weighted summaries.
                const double w=obs_weights.empty()?1.0:obs_weights[i];
                geometry_weights.push_back(w);
                if(!(w>0.0) || !std::isfinite(w)) continue;
                expectation_weight+=w;
                sg+=w*p_c_global[i]; ss+=w*p_c_scoring[i];
                const double delta=p_c_scoring[i]-p_c_global[i];
                sd2+=w*delta*delta;
            }
            if(expectation_weight>0.0){
                mean_global=sg/expectation_weight;
                mean_scoring=ss/expectation_weight;
                mean_delta=mean_scoring-mean_global;
                rms_delta=sqrt(std::max(0.0,sd2/expectation_weight));
            }else expectation_status="no_positive_expectation_weight";
        }

        ParentAxisGeometryAssessment geom_global, geom_scoring;
        string geometry_status="not_applicable_nonheterotypic";
        double alpha_delta=nanv;
        if(model_ok && compiled_hetero){
            geom_global=assess_parent_axis_geometry(p_A,p_B,p_c_global,geometry_weights);
            geom_scoring=assess_parent_axis_geometry(p_A,p_B,p_c_scoring,geometry_weights);
            geometry_status = geom_global.status=="ok" && geom_scoring.status=="ok" ?
                "ok" : ("global="+geom_global.status+";scoring="+geom_scoring.status);
            if(std::isfinite(geom_global.alpha) && std::isfinite(geom_scoring.alpha))
                alpha_delta=geom_scoring.alpha-geom_global.alpha;
        }

        ConditionalInformationAssessment information;
        if(model_ok && compiled_hetero){
            information=assess_conditional_information_rc(dget(allele_ratio_mle),
                dget(contam_rate_mle),n,k,p_A,p_B,p_c_scoring,obs_weights);
        }else information.status="not_applicable_nonheterotypic_or_uncompiled";

        const set<int> exclusions=source_exclusions_for_cell(bc);
        SourceExclusionProfileAssessment exclusion=apply_source_exclusion_profile(
            contam_prof,exclusions,source_exclusion_strength);
        const double raw_a=(parent_a>=0 && exclusions.count(parent_a)>0)
            ? profile_mass(parent_a) : 0.0;
        const double raw_b=(parent_b>=0 && exclusions.count(parent_b)>0)
            ? profile_mass(parent_b) : 0.0;
        const double applied_strength =
            exclusion.status == "degenerate_all_mass_excluded_fallback_global" ?
            0.0 : source_exclusion_strength;
        const double effective_a=exclusions.count(parent_a)>0?applied_strength*raw_a:0.0;
        const double effective_b=exclusions.count(parent_b)>0?applied_strength*raw_b:0.0;

        out << bc_out << '\t' << ident << '\t' << (hetero?1:0)
            << '\t' << dget(contam_rate_mle) << '\t' << dget(contam_rate_regularized)
            << '\t' << dget(allele_ratio_mle) << '\t' << dget(allele_ratio_regularized)
            << '\t' << dget(contam_rate_se) << '\t' << dget(allele_ratio_se)
            << '\t' << dget(contam_ci_low) << '\t' << dget(contam_ci_high)
            << '\t' << dget(allele_ratio_ci_low) << '\t' << dget(allele_ratio_ci_high)
            << '\t' << dget(r_c_correlation) << '\t' << dget(ridge_width)
            << '\t' << dget(prior_mean_by_cell) << '\t' << dget(prior_var_by_cell)
            << '\t' << sget(prior_group_by_cell,"none") << '\t' << dget(prior_displacement)
            << '\t' << dget(mle_log_likelihood) << '\t' << dget(regularized_log_likelihood)
            << '\t' << dget(informative_weight_by_cell) << '\t' << dget(gradient_norm_by_cell)
            << '\t' << iget(prior_training_eligible)
            << '\t' << sget(prior_training_reason_by_cell,"not_evaluated")
            << '\t' << iget(fit_boundary_flags)
            << '\t' << sget(interval_method_by_cell,"not_available")
            << '\t' << sget(optimizer_status_by_cell,"not_available")
            << '\t' << sget(compile_status_by_cell,"not_compiled")
            << '\t' << dget(diagnostic_fixed_r_fallback_c_by_cell)
            << '\t' << dget(diagnostic_fixed_r_fallback_ll_by_cell)
            << '\t' << sget(diagnostic_fallback_status_by_cell,"not_run")
            << '\t' << dget(c_profile_candidate_gap_by_cell)
            << '\t' << dget(c_profile_neighbor_gain_by_cell)
            << '\t' << dget(c_profile_global_gain_by_cell)
            << '\t' << dget(c_profile_likelihood_span_by_cell)
            << '\t' << dget(c_profile_boundary_slope_by_cell)
            << '\t' << dget(c_profile_support_low_by_cell)
            << '\t' << dget(c_profile_support_high_by_cell)
            << '\t' << dget(c_profile_support_width_by_cell)
            << '\t' << iget(r_profile_optimization_succeeded_by_cell)
            << '\t' << iget(r_profile_identifiable_by_cell)
            << '\t' << dget(r_profile_support_low_by_cell)
            << '\t' << dget(r_profile_support_high_by_cell)
            << '\t' << dget(r_profile_support_width_by_cell)
            << '\t' << sget(r_validation_status_by_cell,"not_attempted")
            << '\t' << iget(starts_evaluated_by_cell)
            << '\t' << iget(starts_optimized_by_cell)
            << '\t' << iget(solver_successful_starts_by_cell)
            << '\t' << iget(finite_in_domain_bfgs_candidates_by_cell)
            << '\t' << iget(profile_optimization_attempted_by_cell)
            << '\t' << iget(profile_optimization_succeeded_by_cell)
            << '\t' << iget(profile_validation_passed_by_cell)
            << '\t' << sget(profile_validation_status_by_cell,"not_attempted")
            << '\t' << dget(raw_candidate_r_by_cell)
            << '\t' << dget(raw_candidate_c_by_cell)
            << '\t' << dget(raw_candidate_ll_by_cell)
            << '\t' << dget(validated_candidate_r_by_cell)
            << '\t' << dget(validated_candidate_c_by_cell)
            << '\t' << dget(validated_candidate_ll_by_cell)
            << '\t' << dispersion_phi
            << '\t' << profile_variant_label << '\t' << source_exclusion_strength
            << '\t' << (fix_r_enabled?1:0) << '\t' << fixed_r_basis_label
            << '\t' << (fixed_amb_prof?1:0) << '\t' << fixed_ambient_basis_label
            << '\t' << (truth_assisted_condition?1:0)
            << '\t' << mean_global << '\t' << mean_scoring << '\t' << mean_delta
            << '\t' << rms_delta << '\t' << expectation_weight << '\t' << expectation_status
            << '\t' << geom_global.alpha << '\t' << geom_scoring.alpha << '\t' << alpha_delta
            << '\t' << geom_global.orthogonal_norm << '\t' << geom_scoring.orthogonal_norm
            << '\t' << geom_scoring.denominator << '\t' << geometry_status
            << '\t' << information.information_rr << '\t' << information.information_cc
            << '\t' << information.information_rc
            << '\t' << information.conditional_information_c_given_r
            << '\t' << information.status
            << '\t' << exclusion.excluded_source_mass_raw_total
            << '\t' << raw_a << '\t' << raw_b
            << '\t' << exclusion.effective_removed_source_mass_total
            << '\t' << effective_a << '\t' << effective_b
            << '\t' << exclusion.scoring_profile_renormalization_denominator
            << '\t' << exclusion.global_profile_mass_sum
            << '\t' << exclusion.scoring_profile_mass_sum
            << '\t' << exclusion.status << '\n';
    }
    out.close();
}


void contamFinder3::write_r_c_likelihood_surface(const string& filename,
    const set<string>& selected_barcodes, const vector<string>& samples,
    const string& libname, bool cellranger, bool seurat, bool underscore,
    const string& condition_key, const string& synthetic_id) const{

    ofstream out(filename.c_str());
    if (!out.is_open()){
        fprintf(stderr,"ERROR: could not open r-C likelihood surface: %s\n",filename.c_str());
        return;
    }
    out << std::setprecision(17);
    out << "condition\tsynthetic_id\tbarcode\tprofile_variant\tr_grid\tc_grid"
        << "\tlog_likelihood\tdelta_log_likelihood\tselected_objective"
        << "\tdelta_selected_objective\tis_free_mle_nearest_grid_point"
        << "\tis_selected_fit_nearest_grid_point\tis_fixed_r_nearest_grid_point\n";
    if (selected_barcodes.empty()){
        out.close();
        return;
    }

    vector<string> sample_copy=samples;
    const int grid_steps=50;
    for(auto ai=assn.begin();ai!=assn.end();++ai){
        const unsigned long bc=ai->first;
        string bc_out=bc2str(bc);
        if(seurat && !libname.empty()) bc_out=libname+"_"+bc_out;
        else if(cellranger) bc_out+="-1";
        else if(!libname.empty()) bc_out+=(underscore?"_":"-")+libname;
        if(selected_barcodes.count(bc_out)==0) continue;
        if(ai->second<n_samples) continue;
        pair<int,int> comb=idx_to_hap_comb(ai->second,n_samples);
        if(comb.first==comb.second) continue;
        auto ci=cell_to_idx.find(bc);
        if(ci==cell_to_idx.end()) continue;

        vector<double> n,k,p_e,p_c,p_c_global,p_A,p_B,obs_weights;
        bool hetero=false;
        string status;
        if(!compile_cell_model_data(bc,ci->second,n,k,p_e,p_c,p_c_global,
                p_A,p_B,obs_weights,hetero,status) || !hetero) continue;

        auto dget=[&](const robin_hood::unordered_map<unsigned long,double>& m)->double{
            auto it=m.find(bc); return it==m.end()?std::numeric_limits<double>::quiet_NaN():it->second;
        };
        const double mle_r=dget(allele_ratio_mle);
        const double mle_c=dget(contam_rate_mle);
        const double selected_r=allele_ratio.count(bc)>0?allele_ratio.at(bc):mle_r;
        const double selected_c=contam_rate.count(bc)>0?contam_rate.at(bc):mle_c;
        const double prior_mean=dget(prior_mean_by_cell);
        const double prior_var=dget(prior_var_by_cell);
        double alpha=0.0,beta=0.0;
        const bool selected_uses_prior=std::isfinite(prior_mean) && std::isfinite(prior_var) &&
            contamination_beta_prior_is_log_concave(prior_mean,prior_var,alpha,beta) &&
            std::isfinite(selected_c) && std::isfinite(dget(contam_rate_regularized)) &&
            fabs(selected_c-dget(contam_rate_regularized))<=1e-10;

        double fixed_r=std::numeric_limits<double>::quiet_NaN();
        if(fix_r_enabled){
            auto fr=fixed_r_by_ident.find(ai->second);
            if(fr!=fixed_r_by_ident.end()) fixed_r=fr->second;
        }

        vector<double> ll_values((grid_steps+1)*(grid_steps+1),-INFINITY);
        vector<double> selected_values((grid_steps+1)*(grid_steps+1),-INFINITY);
        double best_ll=-INFINITY,best_selected=-INFINITY;
        int best_ll_ri=-1,best_ll_ci=-1;
        for(int ri=0;ri<=grid_steps;++ri){
            const double rr=(double)ri/(double)grid_steps;
            for(int ci_grid=0;ci_grid<=grid_steps;++ci_grid){
                const double cc=(double)ci_grid/(double)grid_steps;
                const size_t index=(size_t)ri*(grid_steps+1)+(size_t)ci_grid;
                const double ll=eval_three_dataset_objective(rr,cc,n,k,p_A,p_B,p_c,
                    obs_weights,false,-1.0,-1.0);
                const double selected_obj=eval_three_dataset_objective(rr,cc,n,k,p_A,p_B,p_c,
                    obs_weights,selected_uses_prior,prior_mean,prior_var);
                ll_values[index]=ll; selected_values[index]=selected_obj;
                if(std::isfinite(ll) && ll>best_ll){
                    best_ll=ll;
                    best_ll_ri=ri;
                    best_ll_ci=ci_grid;
                }
                if(std::isfinite(selected_obj) && selected_obj>best_selected) best_selected=selected_obj;
            }
        }
        auto nearest_index=[&](double value)->int{
            if(!std::isfinite(value)) return -1;
            int idx=(int)llround(std::max(0.0,std::min(1.0,value))*grid_steps);
            return std::max(0,std::min(grid_steps,idx));
        };
        // The unpenalized grid maximum is the free-MLE reference even in fixed-r
        // mechanism-control conditions, where the selected continuous fit is constrained.
        const int mle_ri=best_ll_ri,mle_ci=best_ll_ci;
        const int sel_ri=nearest_index(selected_r),sel_ci=nearest_index(selected_c);
        const int fix_ri=nearest_index(fixed_r),fix_ci=nearest_index(selected_c);
        for(int ri=0;ri<=grid_steps;++ri){
            const double rr=(double)ri/(double)grid_steps;
            for(int ci_grid=0;ci_grid<=grid_steps;++ci_grid){
                const double cc=(double)ci_grid/(double)grid_steps;
                const size_t index=(size_t)ri*(grid_steps+1)+(size_t)ci_grid;
                const double ll=ll_values[index];
                const double selected_obj=selected_values[index];
                out << condition_key << '\t' << synthetic_id << '\t' << bc_out
                    << '\t' << profile_variant_label << '\t' << rr << '\t' << cc
                    << '\t' << ll << '\t' << (std::isfinite(ll)&&std::isfinite(best_ll)?ll-best_ll:std::numeric_limits<double>::quiet_NaN())
                    << '\t' << selected_obj << '\t' << (std::isfinite(selected_obj)&&std::isfinite(best_selected)?selected_obj-best_selected:std::numeric_limits<double>::quiet_NaN())
                    << '\t' << ((ri==mle_ri && ci_grid==mle_ci)?1:0)
                    << '\t' << ((ri==sel_ri && ci_grid==sel_ci)?1:0)
                    << '\t' << ((ri==fix_ri && ci_grid==fix_ci)?1:0) << '\n';
            }
        }
    }
    out.close();
}


void contamFinder3::write_class_residual_report(const string& filename,
    const vector<string>& samples) const{
    struct Agg {
        double sum;
        double sumsq;
        double abs_sum;
        double obj;
        double nobs;
        double ncells;
        Agg():sum(0),sumsq(0),abs_sum(0),obj(0),nobs(0),ncells(0){}
    };
    map<pair<int,string>,Agg> agg;
    double total_abs = 0.0;

    for (map<unsigned long,vector<int> >::const_iterator ci=cell_to_idx.begin();
        ci!=cell_to_idx.end();++ci){
        unsigned long bc=ci->first;
        if (assn.count(bc)==0 || contam_rate.count(bc)==0) continue;
        double depth=0.0;
        for (size_t z=0;z<ci->second.size();++z) depth+=n_all[ci->second[z]];
        string bin = depth<25?"lt25":(depth<100?"25_99":(depth<500?"100_499":"ge500"));
        pair<int,string> key=make_pair(assn.at(bc),bin);
        agg[key].ncells+=1.0;
        for (size_t z=0;z<ci->second.size();++z){
            int i=ci->second[z];
            double c=contam_rate.at(bc);
            double pc=0.5;
            auto a1=amb_mu.find(type1_all[i]);
            if(a1!=amb_mu.end()){
                auto a2=a1->second.find(type2_all[i]);
                if(a2!=a1->second.end()) pc=a2->second;
            }
            double p;
            if(allele_ratio.count(bc)>0 && p_A_all[i]>=0){
                double r=allele_ratio.at(bc);
                double pa=adjust_p_err(p_A_all[i],e_r,e_a);
                double pb=adjust_p_err(p_B_all[i],e_r,e_a);
                p=(1-c)*(r*pa+(1-r)*pb)+c*pc;
            }else{
                double pe=adjust_p_err(p_e_all[i],e_r,e_a);
                p=(1-c)*pe+c*pc;
            }
            p=clamp_unit_interval(p,1e-6);
            double resid=k_all[i]-n_all[i]*p;
            agg[key].sum+=resid;
            agg[key].sumsq+=resid*resid;
            agg[key].abs_sum+=fabs(resid);
            agg[key].obj+=logbinom(n_all[i],k_all[i],p);
            agg[key].nobs+=1.0;
            total_abs+=fabs(resid);
        }
    }

    FILE* out=fopen(filename.c_str(),"w");
    if(out==NULL){
        fprintf(stderr,"ERROR: could not open class residual report: %s\n",filename.c_str());
        return;
    }
    fprintf(out,"identity_index\tidentity\tdepth_bin\tn_cells\tn_observations\tresidual_mean"
        "\tresidual_variance\tabs_residual_sum\tobjective_contribution\tprofile_leverage\n");
    for(map<pair<int,string>,Agg>::const_iterator it=agg.begin();it!=agg.end();++it){
        const Agg& a=it->second;
        double mean=a.nobs>0?a.sum/a.nobs:0.0;
        double var=a.nobs>1?(a.sumsq-a.sum*a.sum/a.nobs)/(a.nobs-1):0.0;
        double lev=total_abs>0?a.abs_sum/total_abs:0.0;
        vector<string> sample_copy = samples;
        string ident = idx2name(it->first.first, sample_copy);
        fprintf(out,"%d\t%s\t%s\t%.0f\t%.0f\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
            it->first.first,ident.c_str(),it->first.second.c_str(),a.ncells,a.nobs,
            mean,var,a.abs_sum,a.obj,lev);
    }
    fclose(out);
}


/**
 * Compute log likelihood of data set with current parameters
 * Assumes pre-compiled data. This will need to be re-generated after
 * assignments are changed via reclassify_cells().
 */
double contamFinder3::compute_ll(){
    double loglik = 0.0;
    set<unsigned long> prior_counted;
    for (int i = 0; i < (int)n_all.size(); ++i){
        unsigned long cell = idx_to_cell[i];
        if (contam_rate.count(cell) == 0) continue;

        double c = contam_rate[cell];
        double p_c = amb_mu[type1_all[i]][type2_all[i]];
        double binom_p;
        if (allele_ratio.count(cell) > 0 && p_A_all[i] >= 0){
            double r = allele_ratio[cell];
            double p_A = adjust_p_err(p_A_all[i], e_r, e_a);
            double p_B = adjust_p_err(p_B_all[i], e_r, e_a);
            double p_endo = r*p_A + (1.0-r)*p_B;
            binom_p = (1.0-c)*p_endo + c*p_c;
        } else {
            double p_e = adjust_p_err(p_e_all[i], e_r, e_a);
            binom_p = (1.0-c)*p_e + c*p_c;
        }
        binom_p = clamp_unit_interval(binom_p,1e-6);
        loglik += logbinom(n_all[i],k_all[i],binom_p);

        // Priors are per-cell terms. The legacy objective incorrectly added the
        // beta density once per observation, making the acceptance gate depend on
        // a cell's number of count rows. Count it exactly once per cell.
        if (prior_counted.insert(cell).second &&
            prior_mean_by_cell.count(cell) > 0 && prior_var_by_cell.count(cell) > 0 &&
            prior_mean_by_cell.at(cell) > 0.0 && prior_var_by_cell.at(cell) > 0.0){
            loglik += beta_log_density_full(c, prior_mean_by_cell.at(cell),
                prior_var_by_cell.at(cell));
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
        fprintf(stderr, "Leave-source-out diagnostic profiles will be computed\n");
    }
    if (user_prior_set){
        fprintf(stderr, "Fixed prior: mean=%.4f var=%.6f (overrides empirical Bayes)\n",
            contam_cell_prior, contam_cell_prior_var);
    }
    if (r_feedback_enabled){
        fprintf(stderr, "R-feedback: per-cell allele ratios update receiver endogenous "
            "expectations during profile estimation\n");
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

    // -- Phase 2: MLE plus empirical-Bayes regularization --
    // est_contam_cells now preserves the unpenalized MLE, learns priors only
    // from likelihood-supported cells, and emits the regularized fit separately.
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
    // Uses per-cell allele ratios from est_contam_cells to update only the
    // heterotypic receiver endogenous expectation during profile estimation.
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
        map<unsigned long, map<pair<int,int>, map<pair<int,int>, double> > > old_amb_mu_loo_cell = this->amb_mu_loo_cell;
        robin_hood::unordered_map<unsigned long, double> old_contam_rate = this->contam_rate;
        robin_hood::unordered_map<unsigned long, double> old_contam_rate_se = this->contam_rate_se;
        robin_hood::unordered_map<unsigned long, double> old_contam_rate_ll = this->contam_rate_ll;
        robin_hood::unordered_map<unsigned long, double> old_allele_ratio = this->allele_ratio;
        robin_hood::unordered_map<unsigned long, double> old_allele_ratio_se = this->allele_ratio_se;
        robin_hood::unordered_map<unsigned long, double> old_contam_rate_mle = this->contam_rate_mle;
        robin_hood::unordered_map<unsigned long, double> old_contam_rate_regularized = this->contam_rate_regularized;
        robin_hood::unordered_map<unsigned long, double> old_allele_ratio_mle = this->allele_ratio_mle;
        robin_hood::unordered_map<unsigned long, double> old_allele_ratio_regularized = this->allele_ratio_regularized;
        robin_hood::unordered_map<unsigned long, double> old_contam_ci_low = this->contam_ci_low;
        robin_hood::unordered_map<unsigned long, double> old_contam_ci_high = this->contam_ci_high;
        robin_hood::unordered_map<unsigned long, double> old_allele_ratio_ci_low = this->allele_ratio_ci_low;
        robin_hood::unordered_map<unsigned long, double> old_allele_ratio_ci_high = this->allele_ratio_ci_high;
        robin_hood::unordered_map<unsigned long, double> old_r_c_correlation = this->r_c_correlation;
        robin_hood::unordered_map<unsigned long, double> old_ridge_width = this->ridge_width;
        robin_hood::unordered_map<unsigned long, double> old_prior_mean_by_cell = this->prior_mean_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_prior_var_by_cell = this->prior_var_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_prior_displacement = this->prior_displacement;
        robin_hood::unordered_map<unsigned long, double> old_mle_log_likelihood = this->mle_log_likelihood;
        robin_hood::unordered_map<unsigned long, double> old_regularized_log_likelihood = this->regularized_log_likelihood;
        robin_hood::unordered_map<unsigned long, double> old_informative_weight_by_cell = this->informative_weight_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_gradient_norm_by_cell = this->gradient_norm_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_diagnostic_fixed_r_fallback_c_by_cell = this->diagnostic_fixed_r_fallback_c_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_diagnostic_fixed_r_fallback_ll_by_cell = this->diagnostic_fixed_r_fallback_ll_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_c_profile_candidate_gap_by_cell = this->c_profile_candidate_gap_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_c_profile_neighbor_gain_by_cell = this->c_profile_neighbor_gain_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_c_profile_global_gain_by_cell = this->c_profile_global_gain_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_c_profile_likelihood_span_by_cell = this->c_profile_likelihood_span_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_c_profile_boundary_slope_by_cell = this->c_profile_boundary_slope_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_c_profile_support_low_by_cell = this->c_profile_support_low_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_c_profile_support_high_by_cell = this->c_profile_support_high_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_c_profile_support_width_by_cell = this->c_profile_support_width_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_r_profile_support_low_by_cell = this->r_profile_support_low_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_r_profile_support_high_by_cell = this->r_profile_support_high_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_r_profile_support_width_by_cell = this->r_profile_support_width_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_raw_candidate_r_by_cell = this->raw_candidate_r_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_raw_candidate_c_by_cell = this->raw_candidate_c_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_raw_candidate_ll_by_cell = this->raw_candidate_ll_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_validated_candidate_r_by_cell = this->validated_candidate_r_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_validated_candidate_c_by_cell = this->validated_candidate_c_by_cell;
        robin_hood::unordered_map<unsigned long, double> old_validated_candidate_ll_by_cell = this->validated_candidate_ll_by_cell;
        robin_hood::unordered_map<unsigned long, int> old_starts_evaluated_by_cell = this->starts_evaluated_by_cell;
        robin_hood::unordered_map<unsigned long, int> old_starts_optimized_by_cell = this->starts_optimized_by_cell;
        robin_hood::unordered_map<unsigned long, int> old_solver_successful_starts_by_cell = this->solver_successful_starts_by_cell;
        robin_hood::unordered_map<unsigned long, int> old_finite_in_domain_bfgs_candidates_by_cell = this->finite_in_domain_bfgs_candidates_by_cell;
        robin_hood::unordered_map<unsigned long, int> old_profile_optimization_attempted_by_cell = this->profile_optimization_attempted_by_cell;
        robin_hood::unordered_map<unsigned long, int> old_profile_optimization_succeeded_by_cell = this->profile_optimization_succeeded_by_cell;
        robin_hood::unordered_map<unsigned long, int> old_profile_validation_passed_by_cell = this->profile_validation_passed_by_cell;
        robin_hood::unordered_map<unsigned long, int> old_r_profile_optimization_succeeded_by_cell = this->r_profile_optimization_succeeded_by_cell;
        robin_hood::unordered_map<unsigned long, int> old_r_profile_identifiable_by_cell = this->r_profile_identifiable_by_cell;
        robin_hood::unordered_map<unsigned long, int> old_prior_training_eligible = this->prior_training_eligible;
        robin_hood::unordered_map<unsigned long, int> old_fit_boundary_flags = this->fit_boundary_flags;
        map<unsigned long, string> old_prior_group_by_cell = this->prior_group_by_cell;
        map<unsigned long, string> old_interval_method_by_cell = this->interval_method_by_cell;
        map<unsigned long, string> old_optimizer_status_by_cell = this->optimizer_status_by_cell;
        map<unsigned long, string> old_prior_training_reason_by_cell = this->prior_training_reason_by_cell;
        map<unsigned long, string> old_diagnostic_fallback_status_by_cell = this->diagnostic_fallback_status_by_cell;
        map<unsigned long, string> old_profile_validation_status_by_cell = this->profile_validation_status_by_cell;
        map<unsigned long, string> old_r_validation_status_by_cell = this->r_validation_status_by_cell;
        map<int, pair<double,double> > old_fusion_contam_priors = this->fusion_contam_priors;
        pair<double,double> old_heterotypic_contam_prior = this->heterotypic_contam_prior;
        pair<double,double> old_twocomp_contam_prior = this->twocomp_contam_prior;
        bool old_heterotypic_contam_prior_valid = this->heterotypic_contam_prior_valid;
        bool old_twocomp_contam_prior_valid = this->twocomp_contam_prior_valid;
        double old_dispersion_phi = this->dispersion_phi;
        double old_contam_cell_prior = this->contam_cell_prior;
        double old_contam_cell_prior_var = this->contam_cell_prior_var;
        bool old_user_prior_set = this->user_prior_set;
        double old_user_prior_mean = this->user_prior_mean;
        double old_user_prior_var = this->user_prior_var;

        double old_obj = this->compute_ll();
        fprintf(stderr, "R-feedback acceptance gate: current objective before update = %f\n",
            old_obj);

        fprintf(stderr, "Re-estimating ambient profile with r-feedback...\n");
        double dummy2 = -1.0;
        double proposed_profile_ll = this->update_amb_prof_mixture(false, dummy2, false);

        // Source-exclusion expectations are a deterministic function of the
        // proposed global profile. Rebuild them before scoring the proposal;
        // otherwise the proposal is evaluated with stale pre-feedback LOO
        // expectations and an accepted state contains C/r values that were not
        // fit against its published profile.
        if (use_loo){
            fprintf(stderr, "Recomputing source-exclusion profiles for proposed r-informed ambient...\n");
            this->compute_loo_profiles();
        }

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
            // Source-exclusion profiles were rebuilt before the proposal was
            // scored, so the accepted C/r state and published profile are now
            // internally consistent.
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
            this->amb_mu_loo_cell = old_amb_mu_loo_cell;
            this->contam_rate = old_contam_rate;
            this->contam_rate_se = old_contam_rate_se;
            this->contam_rate_ll = old_contam_rate_ll;
            this->allele_ratio = old_allele_ratio;
            this->allele_ratio_se = old_allele_ratio_se;
            this->contam_rate_mle = old_contam_rate_mle;
            this->contam_rate_regularized = old_contam_rate_regularized;
            this->allele_ratio_mle = old_allele_ratio_mle;
            this->allele_ratio_regularized = old_allele_ratio_regularized;
            this->contam_ci_low = old_contam_ci_low;
            this->contam_ci_high = old_contam_ci_high;
            this->allele_ratio_ci_low = old_allele_ratio_ci_low;
            this->allele_ratio_ci_high = old_allele_ratio_ci_high;
            this->r_c_correlation = old_r_c_correlation;
            this->ridge_width = old_ridge_width;
            this->prior_mean_by_cell = old_prior_mean_by_cell;
            this->prior_var_by_cell = old_prior_var_by_cell;
            this->prior_displacement = old_prior_displacement;
            this->mle_log_likelihood = old_mle_log_likelihood;
            this->regularized_log_likelihood = old_regularized_log_likelihood;
            this->informative_weight_by_cell = old_informative_weight_by_cell;
            this->gradient_norm_by_cell = old_gradient_norm_by_cell;
            this->diagnostic_fixed_r_fallback_c_by_cell = old_diagnostic_fixed_r_fallback_c_by_cell;
            this->diagnostic_fixed_r_fallback_ll_by_cell = old_diagnostic_fixed_r_fallback_ll_by_cell;
            this->c_profile_candidate_gap_by_cell = old_c_profile_candidate_gap_by_cell;
            this->c_profile_neighbor_gain_by_cell = old_c_profile_neighbor_gain_by_cell;
            this->c_profile_global_gain_by_cell = old_c_profile_global_gain_by_cell;
            this->c_profile_likelihood_span_by_cell = old_c_profile_likelihood_span_by_cell;
            this->c_profile_boundary_slope_by_cell = old_c_profile_boundary_slope_by_cell;
            this->c_profile_support_low_by_cell = old_c_profile_support_low_by_cell;
            this->c_profile_support_high_by_cell = old_c_profile_support_high_by_cell;
            this->c_profile_support_width_by_cell = old_c_profile_support_width_by_cell;
            this->r_profile_support_low_by_cell = old_r_profile_support_low_by_cell;
            this->r_profile_support_high_by_cell = old_r_profile_support_high_by_cell;
            this->r_profile_support_width_by_cell = old_r_profile_support_width_by_cell;
            this->raw_candidate_r_by_cell = old_raw_candidate_r_by_cell;
            this->raw_candidate_c_by_cell = old_raw_candidate_c_by_cell;
            this->raw_candidate_ll_by_cell = old_raw_candidate_ll_by_cell;
            this->validated_candidate_r_by_cell = old_validated_candidate_r_by_cell;
            this->validated_candidate_c_by_cell = old_validated_candidate_c_by_cell;
            this->validated_candidate_ll_by_cell = old_validated_candidate_ll_by_cell;
            this->starts_evaluated_by_cell = old_starts_evaluated_by_cell;
            this->starts_optimized_by_cell = old_starts_optimized_by_cell;
            this->solver_successful_starts_by_cell = old_solver_successful_starts_by_cell;
            this->finite_in_domain_bfgs_candidates_by_cell = old_finite_in_domain_bfgs_candidates_by_cell;
            this->profile_optimization_attempted_by_cell = old_profile_optimization_attempted_by_cell;
            this->profile_optimization_succeeded_by_cell = old_profile_optimization_succeeded_by_cell;
            this->profile_validation_passed_by_cell = old_profile_validation_passed_by_cell;
            this->r_profile_optimization_succeeded_by_cell = old_r_profile_optimization_succeeded_by_cell;
            this->r_profile_identifiable_by_cell = old_r_profile_identifiable_by_cell;
            this->prior_training_eligible = old_prior_training_eligible;
            this->fit_boundary_flags = old_fit_boundary_flags;
            this->prior_group_by_cell = old_prior_group_by_cell;
            this->interval_method_by_cell = old_interval_method_by_cell;
            this->optimizer_status_by_cell = old_optimizer_status_by_cell;
            this->prior_training_reason_by_cell = old_prior_training_reason_by_cell;
            this->diagnostic_fallback_status_by_cell = old_diagnostic_fallback_status_by_cell;
            this->profile_validation_status_by_cell = old_profile_validation_status_by_cell;
            this->r_validation_status_by_cell = old_r_validation_status_by_cell;
            this->fusion_contam_priors = old_fusion_contam_priors;
            this->heterotypic_contam_prior = old_heterotypic_contam_prior;
            this->twocomp_contam_prior = old_twocomp_contam_prior;
            this->heterotypic_contam_prior_valid = old_heterotypic_contam_prior_valid;
            this->twocomp_contam_prior_valid = old_twocomp_contam_prior_valid;
            this->dispersion_phi = old_dispersion_phi;
            this->contam_cell_prior = old_contam_cell_prior;
            this->contam_cell_prior_var = old_contam_cell_prior_var;
            this->user_prior_set = old_user_prior_set;
            this->user_prior_mean = old_user_prior_mean;
            this->user_prior_var = old_user_prior_var;
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

    // -- Final objective --
    //
    // Retained from V1_R20 so a caller can score this profile against a
    // supplied profile on identical data. The profile is fit over
    // compile_amb_prof_dat rows and evaluated here over compile_data rows, so a
    // fitted profile has no obligation to beat a supplied one in this value.
    // A negative difference is a finding, not an error; do not add a
    // non-negativity gate.
    this->final_loglik = this->compute_ll();
    this->final_loglik_valid = std::isfinite(this->final_loglik);
    fprintf(stderr, "Final model log likelihood: %.6f\n", this->final_loglik);
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
// V1_R12: Degenerate-unit crash fix (no algorithm change on healthy units).
//        Symptom: on low-cell synthetic units, every cell hits one of the two
//        unclassifiable branches in reclassify_cells() ((A) get_max returns
//        a_new==-1 or llr_new<=0; (B) populate_llr_table returns false at
//        tab.n_indvs<2), so all cells are deleted and assn empties. The next
//        outer iteration's mixture solve calls compile_amb_prof_dat (which
//        iterates assn) and builds a zero-row component matrix; add_mixcomp()
//        returns false without setting nmixcomp, add_mixcomp_fracs() returns
//        false, and solve() then indexes an empty mixcompfracs_sparse and
//        segfaults.
//        Three localized guards:
//        (1) reclassify_cells(): MIN_CELLS_FOR_REMOVAL=100 floor evaluated
//            against assn.size() at pass entry. When below the floor,
//            suppress_removal gates both cell_rm.insert sites so unclassifiable
//            cells are kept with their current identity (no reassign, no delete).
//            At or above the floor, removal is byte-for-byte unchanged. Because
//            the gate re-evaluates each iteration against the carried-over count,
//            a unit that drops below the floor cannot erode further.
//        (2) solve_species_level_pi(): guard before the solver. If
//            species_mixfracs is empty, or add_mixcomp / add_mixcomp_fracs
//            return false, warn and return 0.0 (matching the function's other
//            failure paths); the current profile is retained.
//        (3) update_amb_prof_mixture(): same guard on the individual path. If
//            mixfracs is empty, or component setup fails, warn and return
//            -DBL_MAX with contam_prof/amb_mu untouched. Accept/reject keys on
//            compute_ll(), so -DBL_MAX is a non-corrupting no-improvement value.
//        optimML (multivar.cpp) is not modified; it already returns false on the
//        empty/mismatched matrix. The two solver guards make a zero-row matrix
//        non-fatal regardless of cause; the removal floor prevents low-data units
//        from collapsing while preserving result-identity above the floor.
// V1_R13: Individual source-recovery correctness hardening. The profile simplex
//        and receiver-reassignment identity sets are now separate constructor
//        inputs; deterministic --profile_restarts runs an exact total number of
//        supplied/uniform/Dirichlet starts; successful/near-optimal starts and
//        near-optimal profile L1 spread are retained for output diagnostics.
// V1_R14: Historical intermediate revision that accepted any finite endpoint
//        returned when Brent found no derivative root. Superseded by V1_R16
//        because flat, weak, surrogate, or prior-driven endpoints are unsafe.
// V1_R15: Retain weighted hybrid-composition cells that lack a direct pair row
//        by selecting one best-supported singlet or one-component pair basis and
//        deriving the complete endogenous expectation from ordinary condf data.
//        Leave-source-out now excludes all native components of each weighted
//        composition via cell-specific profiles. Cell diagnostics enumerate all
//        assignments, include compile_status, and preserve true fit failures.
// V1_R16: Replace permissive endpoint acceptance with an unpenalized
//        identifiability + likelihood-support + KKT gate. Exactly flat and
//        weak c likelihoods remain explicit failures. Fixed-r fallback after
//        total 2-D optimizer failure is diagnostic only. Empirical Beta priors
//        must be log-concave (alpha,beta >= 1), and regularized boundary
//        candidates can never replace a valid MLE. Weighted-composition fitting
//        and ambient-profile compilation now share one exact row selector; all
//        selected direct/fallback rows require complete conditional support,
//        and candidate ranking prioritizes direct native-component coverage
//        before depth.
// V1_R17: Added deterministic same-model profile refinement and strict
//        profile-based c-boundary validation; superseded as production control
//        flow by V1_R18 because BFGS success still admitted the profile fit.
// V1_R18: Run the exact free-(r,c) profile for every compiled heterotypic cell,
//        independent of BFGS success or coordinate validity. BFGS candidates are
//        optional diagnostics/grid hints only. Remove absolute raw-gradient
//        prior-training gating and use scale-invariant profile support instead.
//        Diagnostics now separate BFGS attempts/successes, finite candidates,
//        profile attempt/success, validation status, and raw/final coordinates.
// V1_R19: Withdrawn marginal self-column experiment; fully reverted in V1_R20.
// V1_R20: Reverted the V1_R19 marginal self-column work after direct
//        measurement refuted the proposed correction. Retained final_loglik and
//        final_loglik_valid as the full-model objective instrument; no point-fit
//        behavior changed relative to the pre-experiment estimator.
// V1_R21: Attach each resampled observation-weight vector to both species and
//        native individual bootstrap solvers, validate replicate dimensions and
//        weights before solving, reject incomplete/non-finite solutions, and
//        report successful/requested counts. Zero-success and failed Dirichlet
//        fits omit concentrations instead of substituting point estimates.
//        Candidate-keyed row-wise uncertainty is suppressed pending clustered
//        cell-level resampling. V1_R20 objective instrumentation is preserved.
// V1_R23: Missing per-cell contamination estimates no longer enter profile
//        compilation or optional reclassification as implicit c=0 values.
//        Reclassification branches without an accepted comparison now retain
//        the prior cell state at every unit size; the arbitrary removal threshold
//        and cell-deletion path are gone.
// ============================================================================
