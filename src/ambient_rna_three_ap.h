// ============================================================================
// ambient_rna_three_ap.h
// Three-component contamination model for tetraploid cells
// with r-feedback and adaptive prior extensions (ambient profile refinements)
//
//
// For heterotypic tetraploid (A+B) cells, decomposes observed allele fractions
// into three sources:
//   p = (1-c) * [r * p_A + (1-r) * p_B] + c * p_c
// where:
//   r  = fraction of endogenous RNA from genome A (vs genome B)
//   c  = ambient RNA contamination fraction
//   p_A = nalt_A / 2  (expected alt fraction from genome A alone)
//   p_B = nalt_B / 2  (expected alt fraction from genome B alone)
//   p_c = ambient pool alt fraction at this SNP category
//
// For diploid singlets and homotypic tetraploids (A+A), collapses to the
// standard two-component model: p = (1-c)*p_e + c*p_c
//
// New features (runtime flags, off by default):
//   --r_feedback: Feed per-cell allele ratio estimates back into the ambient
//     profile mixture solver. Heterotypic cells use their estimated r instead
//     of the hardcoded 0.5 when contributing expected allele fractions.
//   --adaptive_prior: After empirical Bayes estimation, check if the
//     contamination distribution is pathological (boundary pile-up). If so,
//     switch to a fixed prior and iteratively halve variance until the
//     distribution normalizes.
//   --source_exclusion_strength: Scale assigned receiver/parent source mass
//     continuously from the exact global-profile endpoint (0) to the exact
//     historical hard source-exclusion endpoint (1).
//   --r_c_surface_selector/--r_c_surface_out: Emit bounded likelihood geometry
//     only for barcodes selected by the orchestrator-owned signed selector.
// ============================================================================

#ifndef _CELLBOUNCER_AMBIENT_RNA_THREE_AP_H
#define _CELLBOUNCER_AMBIENT_RNA_THREE_AP_H

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
#include <tuple>
#include <math.h>
#include <mixtureDist/functions.h>
#include <optimML/brent.h>
#include <optimML/multivar.h>
#include "common.h"

// One-dimensional contamination identifiability and constrained-boundary
// validation.  Boundary candidates are accepted only for the unpenalized
// model when the data contain meaningful contamination contrast, the
// one-sided derivative satisfies the KKT direction, and the boundary is
// meaningfully better than both the opposite endpoint and an interior grid.
struct UnivariateCInformationAssessment {
    bool identifiable;
    double fisher_information;
    double max_abs_contrast;
    std::string status;
    UnivariateCInformationAssessment();
};

struct UnivariateCBoundaryAssessment {
    bool accepted;
    int boundary;
    double fisher_information;
    double ll_boundary;
    double ll_opposite;
    double ll_best_interior;
    double derivative_near_boundary;
    std::string status;
    UnivariateCBoundaryAssessment();
};

UnivariateCInformationAssessment assess_univariate_c_information(
    const std::vector<double>& n, const std::vector<double>& p_e,
    const std::vector<double>& p_c, const std::vector<double>& weights);

UnivariateCBoundaryAssessment assess_unpenalized_c_boundary(
    double estimate, const std::vector<double>& n,
    const std::vector<double>& k, const std::vector<double>& p_e,
    const std::vector<double>& p_c, const std::vector<double>& weights,
    bool regularized_candidate);

// Empirical Beta priors used for regularization must be log-concave.  Priors
// with alpha < 1 or beta < 1 favor a boundary and are not permitted to replace
// an otherwise valid unpenalized estimate.
bool contamination_beta_prior_is_log_concave(
    double mean, double var, double& alpha, double& beta);

enum HeterotypicFitSource {
    HET_FIT_NONE = 0,
    HET_FIT_TWO_DIMENSIONAL = 1,
    HET_FIT_DIAGNOSTIC_FIXED_R = 2
};

HeterotypicFitSource choose_heterotypic_fit_source(
    int successful_two_dimensional_starts,
    bool diagnostic_fixed_r_candidate_available);

bool heterotypic_fit_source_is_production_valid(HeterotypicFitSource source);

// A regularized result may replace the unpenalized MLE only when it is a
// finite interior c estimate. Boundary regularized candidates remain
// diagnostic and the MLE is retained.
bool regularized_c_candidate_is_selectable(bool fit_success, double c);

// Shared downstream SNP-category filter used after base row selection. The
// caller states whether the model path requires tetraploid category filtering;
// identical inputs then produce identical retained row keys in per-cell and
// ambient-profile compilation.
bool contamination_observation_passes_model_filter(
    bool filter_required, bool ambient_profile_available,
    double p_e, double p_c, double min_signal_gap,
    int nalt1, int nalt2);

// Independent profile-likelihood validation for free heterotypic (r,c) fits.
// A c-boundary candidate is production-valid only when the unpenalized
// profile likelihood L_profile(c)=max_r L(r,c) provides strong endpoint
// support, a KKT-consistent one-sided profile slope, a narrow support region,
// and agreement with the independently profiled endpoint.
struct HeterotypicCBoundaryAssessment {
    bool accepted;
    int boundary;
    double ll_boundary;
    double ll_opposite;
    double ll_best_interior;
    double ll_candidate;
    double candidate_profile_gap;
    double derivative_near_boundary;
    double support_low;
    double support_high;
    double support_width;
    double nuisance_r_at_boundary;
    std::string status;
    HeterotypicCBoundaryAssessment();
};

struct ContinuousProfileIntervalAssessment {
    bool success;
    double lower;
    double upper;
    double width;
    std::string status;
    ContinuousProfileIntervalAssessment();
};

struct TwoDimensionalFitAssessment {
    bool accepted;
    bool c_boundary_candidate;
    bool profile_optimization_attempted;
    bool profile_optimization_succeeded;
    bool raw_candidate_available;
    int finite_in_domain_raw_candidates;
    double raw_candidate_r;
    double raw_candidate_c;
    double raw_candidate_objective;
    double projected_gradient_norm;
    double candidate_profile_gap;
    double profile_neighbor_gain;
    double profile_global_gain;
    double profile_likelihood_span;
    double profile_support_low;
    double profile_support_high;
    double profile_support_width;
    bool r_profile_optimization_succeeded;
    bool r_profile_identifiable;
    double r_profile_support_low;
    double r_profile_support_high;
    double r_profile_support_width;
    std::string r_validation_status;
    double nuisance_r_argmax;
    double validated_r;
    double validated_c;
    double validated_objective;
    HeterotypicCBoundaryAssessment boundary;
    std::string status;
    TwoDimensionalFitAssessment();
};

HeterotypicCBoundaryAssessment assess_unpenalized_heterotypic_c_boundary(
    double candidate_r, double candidate_c,
    const std::vector<double>& n, const std::vector<double>& k,
    const std::vector<double>& p_A, const std::vector<double>& p_B,
    const std::vector<double>& p_c, const std::vector<double>& weights,
    bool regularized_candidate);

TwoDimensionalFitAssessment fit_two_dimensional_profile_model(
    const std::vector<std::pair<double, double> >& raw_candidates,
    const std::vector<double>& n, const std::vector<double>& k,
    const std::vector<double>& p_A, const std::vector<double>& p_B,
    const std::vector<double>& p_c, const std::vector<double>& weights,
    bool apply_prior, double prior_mean, double prior_var);

// Compatibility/test wrapper for one optional raw optimizer candidate. The
// production path calls fit_two_dimensional_profile_model() even when this
// candidate list is empty, so BFGS success is never an admission requirement.
TwoDimensionalFitAssessment assess_two_dimensional_fit_candidate(
    double candidate_r, double candidate_c,
    const std::vector<double>& n, const std::vector<double>& k,
    const std::vector<double>& p_A, const std::vector<double>& p_B,
    const std::vector<double>& p_c, const std::vector<double>& weights,
    bool apply_prior, double prior_mean, double prior_var);

// Scale-invariant empirical-prior eligibility. Raw objective-gradient magnitude
// is intentionally absent because it scales with observation depth.
bool contamination_prior_training_supported(
    bool fit_success, double log_likelihood, double informative_weight,
    double min_informative_weight, bool is_heterotypic,
    bool profile_optimization_attempted, bool profile_optimization_succeeded,
    bool profile_validation_passed, double profile_support_width,
    double max_profile_width, bool c_boundary, std::string& reason);

std::pair<double,double> contamination_prior_safe_moments(
    const std::vector<double>& values);

// Compute the weighted endogenous expectation for one pairwise count row.
// Every positive native component must be represented directly or by a finite
// conditional fraction; no neutral 0.5 substitution is permitted.
bool weighted_composition_expected_from_row_complete(
    const std::map<int, double>& composition,
    const std::pair<int, int>& key1,
    const std::pair<int, int>& key2,
    const std::map<std::pair<int, int>, std::map<int, float> >& conditional_fractions,
    double& expected);

struct WeightedCompositionRowSelection {
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int> > > rows;
    bool direct_component_rows;
    std::pair<int, int> fallback_basis;
    std::string status;
    WeightedCompositionRowSelection();
};

// Select the exact weighted-composition rows used by both per-cell fitting and
// ambient-profile compilation. Direct complete component-pair rows are used
// when available; otherwise exactly one complete conditional basis is chosen.
WeightedCompositionRowSelection select_weighted_composition_rows(
    const std::map<int, double>& composition,
    const std::map<std::pair<int, int>,
        std::map<std::pair<int, int>, std::pair<float, float> > >& allelecounts,
    const std::map<std::pair<int, int>, std::map<int, float> >& conditional_fractions);

// Select one non-duplicating conditional count basis for a weighted native
// composition when no direct complete row is available. Candidates with more
// native components represented directly rank ahead of depth-only alternatives.
std::vector<std::pair<int, int> > rank_weighted_composition_fallback_bases(
    const std::map<int, double>& composition,
    const std::map<std::pair<int, int>,
        std::map<std::pair<int, int>, std::pair<float, float> > >& allelecounts);

std::pair<int, int> select_weighted_composition_fallback_basis(
    const std::map<int, double>& composition,
    const std::map<std::pair<int, int>,
        std::map<std::pair<int, int>, std::pair<float, float> > >& allelecounts);

// Return every positive native component represented by a weighted receiver;
// leave-source-out must remove this complete set, not only the pair label.
std::set<int> weighted_composition_source_exclusions(
    const std::map<int, double>& composition);

// Generalized receiver/source exclusion assessment. The input profile is not
// modified. Endpoints intentionally preserve the legacy behavior: lambda=0
// returns the original global profile and lambda=1 removes all eligible source
// mass before renormalization. Degenerate all-mass exclusion falls back to the
// original global profile with an explicit non-OK status, matching the legacy
// estimator endpoint without concealing the fallback.
struct SourceExclusionProfileAssessment {
    std::map<int, double> scoring_profile;
    double source_exclusion_strength;
    double global_profile_mass_sum;
    double scoring_profile_mass_sum;
    double excluded_source_mass_raw_total;
    double effective_removed_source_mass_total;
    double scoring_profile_renormalization_denominator;
    std::string status;
    SourceExclusionProfileAssessment();
};

SourceExclusionProfileAssessment apply_source_exclusion_profile(
    const std::map<int, double>& global_profile,
    const std::set<int>& excluded_sources,
    double source_exclusion_strength);

struct ParentAxisGeometryAssessment {
    double alpha;
    double orthogonal_norm;
    double denominator;
    std::string status;
    ParentAxisGeometryAssessment();
};

ParentAxisGeometryAssessment assess_parent_axis_geometry(
    const std::vector<double>& p_A, const std::vector<double>& p_B,
    const std::vector<double>& p_ambient, const std::vector<double>& weights);

struct ConditionalInformationAssessment {
    double information_rr;
    double information_cc;
    double information_rc;
    double conditional_information_c_given_r;
    std::string status;
    ConditionalInformationAssessment();
};

ConditionalInformationAssessment assess_conditional_information_rc(
    double r, double c, const std::vector<double>& n,
    const std::vector<double>& k, const std::vector<double>& p_A,
    const std::vector<double>& p_B, const std::vector<double>& p_ambient,
    const std::vector<double>& weights);

// Likelihood callbacks exposed for build-integrated mathematical tests.
double ll_c(double c, const std::map<std::string, double>& data_d,
    const std::map<std::string, int>& data_i);
double dll_dc(double c, const std::map<std::string, double>& data_d,
    const std::map<std::string, int>& data_i);
double d2ll_dc2(double c, const std::map<std::string, double>& data_d,
    const std::map<std::string, int>& data_i);
double ll_three(const std::vector<double>& params,
    const std::map<std::string, double>& data_d,
    const std::map<std::string, int>& data_i);
void dll_three(const std::vector<double>& params,
    const std::map<std::string, double>& data_d,
    const std::map<std::string, int>& data_i, std::vector<double>& results);
double ll_ambmu(const std::vector<double>& params,
    const std::map<std::string, double>& data_d,
    const std::map<std::string, int>& data_i);
void dll_ambmu(const std::vector<double>& params,
    const std::map<std::string, double>& data_d,
    const std::map<std::string, int>& data_i, std::vector<double>& results);
double ll_err_rates(const std::vector<double>& params,
    const std::map<std::string, double>& data_d,
    const std::map<std::string, int>& data_i);
void dll_err_rates(const std::vector<double>& params,
    const std::map<std::string, double>& data_d,
    const std::map<std::string, int>& data_i, std::vector<double>& results);
double ll_amb_prof_mixture(const std::vector<double>& params,
    const std::map<std::string, double>& data_d,
    const std::map<std::string, int>& data_i);
void dll_amb_prof_mixture(const std::vector<double>& params,
    const std::map<std::string, double>& data_d,
    const std::map<std::string, int>& data_i, std::vector<double>& results);

class contamFinder3{
    private:

        // ---- Copies of external data structures ----

        robin_hood::unordered_map<unsigned long,
                std::map<std::pair<int, int>,
                    std::map<std::pair<int, int>,
                        std::pair<float, float> > > > indv_allelecounts;
        // Ambient profile source candidates are separate from reassignment identities.
        std::set<int> allowed_ids;
        std::set<int> reassign_allowed_ids;
        std::set<int> allowed_ids2;

        // Map each ID to the sum of all LLRs of all cells assigned to it
        std::map<int, double> id_llrsum;
        std::map<int, double> id_llrsum2;
        std::map<int, double> id_count;
        double llrtot;

        // Expected doublet rate (for reclassifying cells)
        double doublet_rate;

        // ---- Per-observation data vectors ----

        std::vector<double> n_all;
        std::vector<double> k_all;
        std::vector<double> p_e_all;   // two-component p_e = (nalt_A + nalt_B)/4 or nalt/2
        std::vector<double> p_A_all;   // genome A expected alt frac (nalt_A/2), -1 for singlets
        std::vector<double> p_B_all;   // genome B expected alt frac (nalt_B/2), -1 for singlets
        std::vector<std::pair<int, int> > type1_all;
        std::vector<std::pair<int, int> > type2_all;
        std::vector<double> weight_all;

        // ---- Index mappings ----

        std::map<std::pair<int, int>, std::map<std::pair<int, int>,
            std::vector<int> > > expfrac_to_idx;
        std::map<unsigned long, std::vector<int> > cell_to_idx;
        std::map<int, unsigned long> idx_to_cell;
        std::map<std::pair<int, int>, std::vector<int> > sitecomb_type_to_idx;

        // ---- Sample info ----

        int n_samples;
        std::vector<int> idx2samp;

        // Expected alt allele fractions conditional on genotype category
        std::map<std::pair<int, int>, std::map<int, float> > expfracs;

        // ---- Prior on contamination rate ----

        double contam_cell_prior;
        double contam_cell_prior_var;

        // ---- Error rates ----

        double e_r;
        double e_a;

        // ---- Model settings ----

        bool inter_species;
        int n_mixprop_trials;
        int profile_total_starts;
        double delta_thresh;
        int maxits;
        int nits;
        bool weighted;
        int num_threads;
        bool skip_reassign;

        // Per-cell fit result used by the OpenMP-parallel contamination
        // estimator.  Worker threads fill these POD-like records only; the
        // owning maps are updated serially after the parallel region.
        struct CellContamFitResult {
            unsigned long barcode;
            double c;
            double c_se;
            double ll;
            double r;
            double r_se;
            double c_ci_low;
            double c_ci_high;
            double r_ci_low;
            double r_ci_high;
            double r_c_correlation;
            double ridge_width;
            double gradient_norm;
            double informative_weight;
            double pearson_ss;
            double pearson_df;
            double diagnostic_fixed_r_fallback_c;
            double diagnostic_fixed_r_fallback_ll;
            double c_profile_candidate_gap;
            double c_profile_neighbor_gain;
            double c_profile_global_gain;
            double c_profile_likelihood_span;
            double c_profile_boundary_slope;
            double c_profile_support_low;
            double c_profile_support_high;
            double c_profile_support_width;
            double r_profile_support_low;
            double r_profile_support_high;
            double r_profile_support_width;
            double raw_candidate_r;
            double raw_candidate_c;
            double raw_candidate_ll;
            double validated_candidate_r;
            double validated_candidate_c;
            double validated_candidate_ll;
            int starts_evaluated;
            int starts_optimized;
            int solver_successful_starts;
            int finite_in_domain_bfgs_candidates;
            int profile_optimization_attempted;
            int profile_optimization_succeeded;
            int profile_validation_passed;
            int r_profile_optimization_succeeded;
            int r_profile_identifiable;
            int optimizer_iterations;
            bool is_heterotypic;
            bool bfgs_fallback;
            bool has_allele_ratio;
            bool fit_success;
            bool interval_success;
            bool hit_c_lower;
            bool hit_c_upper;
            bool hit_r_lower;
            bool hit_r_upper;
            bool prior_applied;
            std::string interval_method;
            std::string optimizer_status;
            std::string r_validation_status;
            std::string diagnostic_fallback_status;
            std::string prior_group;
        };

        enum ContamPriorMode {
            PRIOR_NONE = 0,
            PRIOR_GLOBAL = 1,
            PRIOR_HETEROTYPIC = 2,
            PRIOR_FUSION = 3
        };

        enum HeterotypicStartMode {
            START_SINGLE = 0,
            START_TOPK = 1,
            START_EXHAUSTIVE = 2
        };

        bool thorough_multistart;
        bool adaptive_multistart;
        HeterotypicStartMode heterotypic_start_mode;
        int heterotypic_start_top_k;
        bool adaptive_profile_intervals;

        ContamPriorMode contam_prior_mode;
        int contam_prior_min_cells;
        double contam_prior_max_ci_width;
        double contam_prior_min_informative_weight;
        std::map<int, std::pair<double, double> > fusion_contam_priors;
        std::pair<double, double> heterotypic_contam_prior;
        std::pair<double, double> twocomp_contam_prior;
        bool heterotypic_contam_prior_valid;
        bool twocomp_contam_prior_valid;
        double dispersion_phi;

        std::pair<double, double> prior_for_cell(unsigned long barcode,
            bool is_heterotypic, std::string& group) const;
        void learn_contam_priors(const std::vector<CellContamFitResult>& mle_results);
        bool prior_training_supported(const CellContamFitResult& fit,
            std::string& reason) const;

        // ---- Tetraploid-aware mode ----

        // Optional per-cell species composition override.  Native species mode can
        // derive a biologically weighted endogenous composition from the original
        // individual assignment (e.g. Chinobo-mCherry+JOS3C1 ->
        // B:0.25,C:0.25,O:0.5).  Keys are numeric cell barcodes; values are
        // sample/species-index -> dosage weight maps that sum to 1.0.
        std::map<unsigned long, std::map<int, double> > cell_composition_overrides;
        bool has_composition_override(unsigned long barcode) const;
        bool composition_expected_from_row_complete(
            const std::map<int, double>& comp,
            const std::pair<int, int>& key1,
            const std::pair<int, int>& key2, double& expected) const;

        std::set<int> locked_identities;
        std::set<int> safe_singlets;
        bool tetraploid_aware;
        double min_signal_gap;
        bool amb_mu_available;
        bool ids_restricted;

        // ---- O53 step 2: candidate-keyed profile rows --------------------
        //
        // The historical per-cell profile data compiler keys every likelihood
        // row on the cell's own assigned identity: three marginal bins for a
        // singlet, nine pairwise bins for a fusion. Measured 2026-07-26 across
        // 40 benchmark units, that design is rank deficient by 13 to 16 of the
        // 23 free ambient-profile directions on every unit, at every ambient
        // level, in both ploidies. The bins that carry the missing directions
        // are the marginal bins keyed on the other candidate donors. They are
        // already loaded and already have conditional-fraction entries; the
        // compiler simply never visits them.
        //
        // When enabled, the compiler emits a second block of rows per cell:
        // the marginal bin of every ambient candidate in idx2samp, at every
        // genotype dosage present for that candidate in this cell.
        //
        // READ MULTIPLICITY. A physical read sits at one site, and at that
        // site every candidate has a dosage, so the read enters one marginal
        // bin per candidate. Emitting the marginal block therefore replicates
        // each read once per candidate. Two corrections are applied together
        // and neither is optional:
        //
        //   1. the marginal block is divided by the candidate count, so the
        //      block as a whole carries one unit of weight per read;
        //   2. the two blocks share one unit of weight between them, governed
        //      by candidate_keyed_split, so a read contributes exactly the
        //      same total evidence as it did before the flag existed.
        //
        // Without (2), point estimates still improve while every concentration,
        // bootstrap interval, and cross-condition likelihood comparison
        // silently changes meaning, because each read would be counted about
        // twice.
        //
        // candidate_keyed_split is the fraction of the per-read weight budget
        // retained by the historical assignment-keyed block. 1.0 reproduces
        // the historical design exactly and is an exact no-op; 0.0 uses
        // candidate-keyed marginals alone. Neither endpoint emits zero-weight
        // rows: the block whose share is zero is skipped entirely. The default of 0.5 keeps both blocks, because the pairwise
        // block carries the endogenous parental contrast that pins the
        // contamination rate, and the contamination-rate table must not
        // regress.
        //
        // Individual path only. Candidate-keyed expansion of the species path
        // was measured and rejected: worse on 20 of 20 native species units.

        bool candidate_keyed_rows;
        double candidate_keyed_split;

        // ---- LOO ambient profiles ----

        bool use_loo;
        double source_exclusion_strength;
        // Receiver-data cross-fitting: barcodes in this set are scored by the
        // estimator but are excluded from ambient-profile sufficient statistics.
        // This preserves every source column and total-C semantics while
        // preventing a scored receiver from training its own nuisance profile.
        std::set<unsigned long> profile_holdout_barcodes;
        std::string profile_holdout_basis_label;
        std::string fixed_r_basis_label;
        std::string fixed_ambient_basis_label;
        std::string profile_variant_label;
        bool truth_assisted_condition;
        std::map<int, std::map<std::pair<int,int>,
            std::map<std::pair<int,int>, double>>> amb_mu_loo;
        std::map<unsigned long, std::map<std::pair<int,int>,
            std::map<std::pair<int,int>, double>>> amb_mu_loo_cell;
        void compute_loo_profiles();
        std::set<int> source_exclusions_for_cell(unsigned long barcode) const;
        bool compile_cell_model_data(unsigned long barcode,
            const std::vector<int>& obs_idx, std::vector<double>& n,
            std::vector<double>& k, std::vector<double>& p_e,
            std::vector<double>& p_c_scoring, std::vector<double>& p_c_global,
            std::vector<double>& p_A, std::vector<double>& p_B,
            std::vector<double>& observation_weights, bool& is_heterotypic,
            std::string& status) const;

        // ---- Fisher Information weighting ----

        bool use_fi_weight;

        // ---- User-specified prior override ----

        bool user_prior_set;
        double user_prior_mean;
        double user_prior_var;

        // ---- R-feedback: update receiver endogenous expectation only ----

        bool r_feedback_enabled;

        // ---- Fixed allele ratio (pooled-r experiment) --------------------

        bool fix_r_enabled;
        std::map<int, double> fixed_r_by_ident;

        // ---- Adaptive prior: empirical Bayes -> fixed prior fallback ----

        bool adaptive_prior_enabled;
        double adaptive_prior_mean;
        double adaptive_prior_init_var;
        double adaptive_prior_min_var;
        double adaptive_prior_boundary_thresh;
        double adaptive_prior_boundary_low;
        double adaptive_prior_boundary_high;
        int adaptive_prior_max_shrink_steps;

        // ---- Fixed ambient profile (Step 0a) ----

        bool fixed_amb_prof;
        void rebuild_amb_mu_from_contam_prof();

        // ---- Species mode ----

        bool species_mode;
        PanelMetadata panel_meta;
        std::map<std::string, double> species_prior_prof;
        std::map<std::string, double> species_prior_conc;
        void expand_species_prior_to_indiv();
        double solve_species_level_pi();

        // ---- Species-diagnostic counts (--species_counts) ----

        bool has_species_counts;
        bool primary_species_counts_enabled;  // true when native --interspecies uses .species_counts as primary counts
        robin_hood::unordered_map<unsigned long,
                std::map<std::pair<int, int>,
                    std::map<std::pair<int, int>,
                        std::pair<float, float> > > > species_allelecounts;
        std::map<std::pair<int, int>, std::map<int, float> > species_expfracs;

        // ---- Loading weights for within-species individual weighting ----

        std::map<int, double> indiv_loading_weights;
        void compute_loading_weights();

        // ---- Warm-start init profiles ----

        std::map<std::string, double> species_init_prof;  // species-level warm start
        bool species_init_used;  // becomes true after first solve_species_level_pi call

        // ---- Bulk mode (for empty-drops tool) ----

        bool bulk_mode;

        // ---- Internal methods ----

        void compile_data(robin_hood::unordered_map<unsigned long, int>& assn,
            robin_hood::unordered_map<unsigned long,
                std::map<std::pair<int, int>,
                std::map<std::pair<int, int>,
                std::pair<float, float> > > >& indv_allelecounts);

        void clear_data();

        void get_reads_expectations(unsigned long barcode,
            int ident,
            std::map<std::pair<int, int>,
                std::map<std::pair<int, int>,
                std::pair<float, float> > >& allelecounts,
            std::vector<double>& n,
            std::vector<double>& k,
            std::vector<double>& p_e,
            std::vector<double>& p_A,
            std::vector<double>& p_B,
            std::vector<std::pair<int, int> >& type1,
            std::vector<std::pair<int, int> >& type2);

        // Per-cell estimation
        void est_contam_cells();
        void est_contam_cells_global();
        CellContamFitResult fit_one_contam_cell(unsigned long barcode,
            const std::vector<int>& obs_idx, bool apply_prior,
            double prior_mean, double prior_var, const std::string& prior_group);

        // Ambient profile
        double update_ambient_profile(bool global_c = false);
        std::pair<double, double> est_error_rates(bool init);

        double c_init_global;

        void compile_amb_prof_dat(bool solve_for_c,
            bool use_global_c,
            std::vector<std::vector<double> >& mixfracs,
            std::vector<double>& weights,
            std::vector<double>& n,
            std::vector<double>& k,
            std::vector<double>& p_e,
            std::vector<double>& c);

        // Bulk-mode data compiler: includes all count rows regardless of
        // identity assignment. Used by tet_ambient_profile (empty drops).
        void compile_bulk_amb_prof_dat(
            std::vector<std::vector<double> >& mixfracs,
            std::vector<double>& weights,
            std::vector<double>& n,
            std::vector<double>& k,
            std::vector<double>& p_e,
            std::vector<double>& c);

        double update_amb_prof_mixture(bool est_c, double& global_c,
            bool use_global_c = false);
        double est_min_c();
        bool reclassify_cells();
        double init_params(double& c);
        bool contam_prof_initialized;
        bool c_initialized;
        double c_init;

        bool ef_all_avg;

    public:

        void set_init_contam_prof(std::map<int, double>& cp);
        void set_init_c(double c);
        void set_doublet_rate(double d);

        // ---- Public data (accessed by driver for output) ----

        // Full-model log likelihood at the accepted parameters, retained so a
        // caller can score a fitted profile against a supplied one on identical
        // data. Under a correctly specified model the difference between the
        // fitted and the true profile is on the order of half the number of free
        // profile parameters, so tens for a roster of this size. A difference
        // orders of magnitude larger means the likelihood maximum does not sit
        // at the truth, which is misspecification rather than search failure and
        // is not reachable by conditioning or regularization work.
        //
        // ⚠️ The fit and the evaluation do not share a row set: the profile is
        // fit over compile_amb_prof_dat rows and this is evaluated over
        // compile_data rows. A fitted profile therefore has no obligation to
        // beat a fixed one here, and a negative difference is a finding rather
        // than an error. Do not reintroduce a non-negativity check.
        double final_loglik;
        bool final_loglik_valid;

        robin_hood::unordered_map<unsigned long, int> assn;
        robin_hood::unordered_map<unsigned long, double> assn_llr;

        // Ambient profile: type1 -> type2 -> expected alt frac in ambient
        std::map<std::pair<int, int>, std::map<std::pair<int, int>, double> > amb_mu;

        // Per-cell contamination rate, SE, and log likelihood
        robin_hood::unordered_map<unsigned long, double> contam_rate;
        robin_hood::unordered_map<unsigned long, double> contam_rate_se;
        robin_hood::unordered_map<unsigned long, double> contam_rate_ll;

        // Per-cell allele ratio r (genome A fraction) -- only for heterotypic combos
        robin_hood::unordered_map<unsigned long, double> allele_ratio;
        robin_hood::unordered_map<unsigned long, double> allele_ratio_se;

        // Explicit unpenalized and regularized outputs. contam_rate/allele_ratio
        // retain the selected regularized values for backward compatibility.
        robin_hood::unordered_map<unsigned long, double> contam_rate_mle;
        robin_hood::unordered_map<unsigned long, double> contam_rate_regularized;
        robin_hood::unordered_map<unsigned long, double> allele_ratio_mle;
        robin_hood::unordered_map<unsigned long, double> allele_ratio_regularized;
        robin_hood::unordered_map<unsigned long, double> contam_ci_low;
        robin_hood::unordered_map<unsigned long, double> contam_ci_high;
        robin_hood::unordered_map<unsigned long, double> allele_ratio_ci_low;
        robin_hood::unordered_map<unsigned long, double> allele_ratio_ci_high;
        robin_hood::unordered_map<unsigned long, double> r_c_correlation;
        robin_hood::unordered_map<unsigned long, double> ridge_width;
        robin_hood::unordered_map<unsigned long, double> prior_mean_by_cell;
        robin_hood::unordered_map<unsigned long, double> prior_var_by_cell;
        robin_hood::unordered_map<unsigned long, double> prior_displacement;
        robin_hood::unordered_map<unsigned long, double> mle_log_likelihood;
        robin_hood::unordered_map<unsigned long, double> regularized_log_likelihood;
        robin_hood::unordered_map<unsigned long, double> informative_weight_by_cell;
        robin_hood::unordered_map<unsigned long, double> gradient_norm_by_cell;
        robin_hood::unordered_map<unsigned long, double> diagnostic_fixed_r_fallback_c_by_cell;
        robin_hood::unordered_map<unsigned long, double> diagnostic_fixed_r_fallback_ll_by_cell;
        robin_hood::unordered_map<unsigned long, double> c_profile_candidate_gap_by_cell;
        robin_hood::unordered_map<unsigned long, double> c_profile_neighbor_gain_by_cell;
        robin_hood::unordered_map<unsigned long, double> c_profile_global_gain_by_cell;
        robin_hood::unordered_map<unsigned long, double> c_profile_likelihood_span_by_cell;
        robin_hood::unordered_map<unsigned long, double> c_profile_boundary_slope_by_cell;
        robin_hood::unordered_map<unsigned long, double> c_profile_support_low_by_cell;
        robin_hood::unordered_map<unsigned long, double> c_profile_support_high_by_cell;
        robin_hood::unordered_map<unsigned long, double> c_profile_support_width_by_cell;
        robin_hood::unordered_map<unsigned long, double> r_profile_support_low_by_cell;
        robin_hood::unordered_map<unsigned long, double> r_profile_support_high_by_cell;
        robin_hood::unordered_map<unsigned long, double> r_profile_support_width_by_cell;
        robin_hood::unordered_map<unsigned long, double> raw_candidate_r_by_cell;
        robin_hood::unordered_map<unsigned long, double> raw_candidate_c_by_cell;
        robin_hood::unordered_map<unsigned long, double> raw_candidate_ll_by_cell;
        robin_hood::unordered_map<unsigned long, double> validated_candidate_r_by_cell;
        robin_hood::unordered_map<unsigned long, double> validated_candidate_c_by_cell;
        robin_hood::unordered_map<unsigned long, double> validated_candidate_ll_by_cell;
        robin_hood::unordered_map<unsigned long, int> starts_evaluated_by_cell;
        robin_hood::unordered_map<unsigned long, int> starts_optimized_by_cell;
        robin_hood::unordered_map<unsigned long, int> solver_successful_starts_by_cell;
        robin_hood::unordered_map<unsigned long, int> finite_in_domain_bfgs_candidates_by_cell;
        robin_hood::unordered_map<unsigned long, int> profile_optimization_attempted_by_cell;
        robin_hood::unordered_map<unsigned long, int> profile_optimization_succeeded_by_cell;
        robin_hood::unordered_map<unsigned long, int> profile_validation_passed_by_cell;
        robin_hood::unordered_map<unsigned long, int> r_profile_optimization_succeeded_by_cell;
        robin_hood::unordered_map<unsigned long, int> r_profile_identifiable_by_cell;
        robin_hood::unordered_map<unsigned long, int> prior_training_eligible;
        robin_hood::unordered_map<unsigned long, int> fit_boundary_flags;
        std::map<unsigned long, std::string> prior_group_by_cell;
        std::map<unsigned long, std::string> interval_method_by_cell;
        std::map<unsigned long, std::string> optimizer_status_by_cell;
        std::map<unsigned long, std::string> prior_training_reason_by_cell;
        std::map<unsigned long, std::string> compile_status_by_cell;
        std::map<unsigned long, std::string> diagnostic_fallback_status_by_cell;
        std::map<unsigned long, std::string> profile_validation_status_by_cell;
        std::map<unsigned long, std::string> r_validation_status_by_cell;

        // Ambient profile mixture proportions
        std::map<int, double> contam_prof;

        // Ambient-profile optimizer stability diagnostics. These describe the
        // final profile solve, not per-cell contamination fits. Missing values
        // are represented by -DBL_MAX/NaN in the driver output, never by zero.
        //
        // NOTE ON SCOPE: ordinary cell-profile refinement uses a single no-c
        // solve after its initial multistart solve. Bulk empty-droplet mode can
        // explicitly request deterministic multistart on the no-c path through
        // set_profile_total_starts(); in that case these fields describe the
        // selected bulk multistart optimum directly. Consumers should use the
        // multistart_* fields below to distinguish configured, successful, and
        // near-optimal starts.
        int profile_successful_starts;
        int profile_near_optimal_count;
        double profile_best_ll;
        double profile_second_best_ll;
        double profile_near_optimal_l1_spread;

        // Ambient-profile MULTISTART diagnostics. Populated by the
        // solve_for_c path and by an explicitly requested bulk no-c multistart;
        // ordinary single-start refinement does not overwrite them. These
        // retain the comparison across starting points:
        // how many starts were configured, how many returned finite results,
        // how many landed within the near-optimal likelihood tolerance of the
        // best, and the largest L1 distance in mixture space among those
        // near-optimal solutions. A large near_optimal_count with a large
        // l1_spread is the signature of a flat ridge (many equally good and
        // materially different profiles); a count of one indicates a uniquely
        // selected optimum. multistart_attempted distinguishes "not run"
        // (species-mode or fixed-profile bypass) from "run and failed".
        // Missing values are -DBL_MAX/NaN, never zero.
        bool multistart_attempted;
        int multistart_configured_starts;
        int multistart_successful_starts;
        int multistart_near_optimal_count;
        double multistart_best_ll;
        double multistart_second_best_ll;
        double multistart_near_optimal_l1_spread;

        // ---- Constructor ----

        // Backward-compatible constructor used by legacy tools such as
        // quant3_contam_ap. Historically one allowed-ID set controlled both
        // ambient-profile candidates and receiver reassignment. Preserve that
        // behavior by delegating to the explicit three-set constructor below.
        contamFinder3(robin_hood::unordered_map<unsigned long,
                std::map<std::pair<int, int>,
                    std::map<std::pair<int, int>,
                        std::pair<float, float> > > >& indv_allelecounts,
            robin_hood::unordered_map<unsigned long, int>& assn,
            robin_hood::unordered_map<unsigned long, double>& assn_llr,
            std::map<std::pair<int, int>, std::map<int, float> >& exp_match_fracs,
            int n_samples,
            std::set<int>& allowed_ids,
            std::set<int>& allowed_ids2);

        // Explicit constructor used by the corrected benchmark path. The
        // ambient source simplex and receiver reassignment universe are
        // intentionally separate.
        contamFinder3(robin_hood::unordered_map<unsigned long,
                std::map<std::pair<int, int>,
                    std::map<std::pair<int, int>,
                        std::pair<float, float> > > >& indv_allelecounts,
            robin_hood::unordered_map<unsigned long, int>& assn,
            robin_hood::unordered_map<unsigned long, double>& assn_llr,
            std::map<std::pair<int, int>, std::map<int, float> >& exp_match_fracs,
            int n_samples,
            std::set<int>& profile_allowed_ids,
            std::set<int>& reassign_allowed_ids,
            std::set<int>& allowed_ids2);

        // ---- Configuration methods ----

        void set_error_rates(double e_r, double e_a);
        void model_other_species();
        void model_single_species();
        void set_mixprop_trials(int nt);
        void set_profile_total_starts(int nt);
        void set_delta(double d);
        void set_maxiter(int i);
        void use_weights();
        void no_weights();
        void set_num_threads(int nt);
        void set_thorough_multistart(bool enabled);
        void set_adaptive_multistart(bool enabled);
        void set_heterotypic_start_mode(const std::string& mode, int top_k = 5);
        void set_adaptive_profile_intervals(bool enabled);
        void set_contam_prior_mode(const std::string& mode);
        void set_contam_prior_support(int min_cells, double max_ci_width,
            double min_informative_weight, double deprecated_max_gradient_norm);
        void no_reassign();

        // Tetraploid-aware mode
        void set_locked_identities(const std::set<int>& locked);
        void set_safe_singlets(const std::set<int>& safe);
        void set_tetraploid_aware(bool enabled);
        void set_min_signal_gap(double gap);
        void set_ids_restricted(bool restricted);

        // O53 step 2: candidate-keyed profile rows, individual path only.
        // split is the fraction of the per-read weight budget kept by the
        // historical assignment-keyed block; it must lie in [0, 1].
        void set_candidate_keyed_rows(bool enabled, double split);
        void set_cell_composition_overrides(
            const std::map<unsigned long, std::map<int, double> >& overrides);

        // Source exclusion. set_use_loo() is retained as the exact legacy
        // hard-exclusion endpoint; the strength setter supplies the generalized
        // lambda in [0,1].
        void set_use_loo(bool enabled);
        void set_source_exclusion_strength(double strength);
        void set_profile_holdout_barcodes(
            const std::set<unsigned long>& held_out,
            const std::string& basis_label = "library_crossfit");
        void set_diagnostic_context(const std::string& fixed_r_basis,
            const std::string& fixed_ambient_basis, bool truth_assisted);

        // FI weighting
        void set_fi_weight(bool enabled);

        // Prior override
        void set_user_prior(double mean, double var);

        // R-feedback
        void set_r_feedback(bool enabled);

        // ---- Fixed allele ratio (pooled-r experiment) --------------------
        // Supply an externally determined allele ratio, keyed by identity
        // index (the same index that idx2name() renders into the `identity`
        // column of the .allele_ratio output).
        //
        // When an identity has an entry here, r is NOT fitted for cells of
        // that identity. The endogenous expectation is built from the supplied
        // r and only c is solved, in one dimension. r therefore cannot track
        // the cell's own c, which is the mechanism by which a free per-cell r
        // absorbs contamination.
        //
        // The pooling policy lives entirely in whoever builds the table. The
        // estimator does not know or care how the value was derived.
        void set_fixed_r(const std::map<int, double>& fr);
        void report_r_feedback_stats();

        // Adaptive prior
        struct BoundaryDiag {
            double frac_low;
            double frac_high;
            double frac_boundary;
            double median_c;
            int n_cells;
        };
        void set_adaptive_prior(bool enabled,
                                double mean = 0.05,
                                double init_var = 0.04,
                                double min_var = 0.001,
                                double boundary_thresh = 0.20);
        BoundaryDiag diagnose_contam_distribution();
        void run_adaptive_prior();

        // Fixed ambient profile (Step 0a)
        void set_fixed_amb_prof(bool enabled);

        // Conditional-fraction coverage validation.  A .condf row is required
        // for every observed genotype category x active ambient-source donor
        // lookup that is not the donor's own directly observed genotype.
        struct CondfCoverageReport {
            unsigned long long required_lookups;
            unsigned long long missing_lookups;
            double required_weight;
            double missing_weight;
        };
        CondfCoverageReport write_condf_coverage_report(
            const std::string& filename,
            const std::vector<std::string>& samples) const;

        // Species mode
        void set_species_mode(const PanelMetadata& panel);
        void set_species_prior(const std::map<std::string, double>& species_prof,
                               const std::map<std::string, double>& species_prof_conc);

        // Warm-start init profiles
        void set_species_init(const std::map<std::string, double>& init_prof);
        void set_indiv_init(const std::map<int, double>& init_prof);

        // Species-diagnostic counts
        void set_primary_species_counts_enabled(bool enabled);
        void set_species_counts(
            robin_hood::unordered_map<unsigned long,
                std::map<std::pair<int, int>,
                    std::map<std::pair<int, int>,
                        std::pair<float, float> > > >& counts,
            std::map<std::pair<int, int>, std::map<int, float> >& condf);

        // Bulk mode (for quant3_contam_empty_drops)
        void set_bulk_mode(bool enabled);

        // Species output
        std::map<std::string, double> species_contam_prof;
        std::map<std::string, double> species_contam_prof_conc;

        // Bootstrap
        void bootstrap_amb_prof(int n_boots,
            std::map<int, double>& contam_prof_cont);

        // Compute total log likelihood
        double compute_ll();

        // Write extended diagnostics for the accepted estimator state.
        void write_cell_fit_diagnostics(const std::string& filename,
            const std::vector<std::string>& samples, const std::string& libname,
            bool cellranger, bool seurat, bool underscore) const;
        // Long-format, per-cell source profile actually used for scoring.
        // This is required to audit source allocation, gate-selected profiles,
        // and cross-fitted total-C semantics without reconstructing estimator
        // state from aggregate .contam_prof files.
        void write_cell_source_profile(const std::string& filename,
            const std::vector<std::string>& samples, const std::string& libname,
            bool cellranger, bool seurat, bool underscore) const;
        void write_class_residual_report(const std::string& filename,
            const std::vector<std::string>& samples) const;
        void write_r_c_likelihood_surface(const std::string& filename,
            const std::set<std::string>& selected_barcodes,
            const std::vector<std::string>& samples, const std::string& libname,
            bool cellranger, bool seurat, bool underscore,
            const std::string& condition_key, const std::string& synthetic_id) const;

        // Run the full estimation pipeline
        void fit();
};

#endif

// ============================================================================
// Revision History
// V1_R1: Forked from ambient_rna_three.h V1_R3. Added r-feedback flag and
//        setter, adaptive prior members (BoundaryDiag struct, set_adaptive_prior,
//        diagnose_contam_distribution, run_adaptive_prior), report_r_feedback_stats.
// V1_R2: Added fixed_amb_prof support (set_fixed_amb_prof, rebuild_amb_mu_from_contam_prof).
//        Added species mode (set_species_mode, set_species_prior,
//        expand_species_prior_to_indiv, solve_species_level_pi, species output maps).
//        Added bulk mode (set_bulk_mode) for quant3_contam_empty_drops.
// V1_R3: No header changes; see ambient_rna_three_ap.cpp V1_R3 for dimension
//        mismatch fix in expand_species_prior_to_indiv and solve_species_level_pi.
// V1_R4: No header changes; see ambient_rna_three_ap.cpp V1_R4 for est_min_c
//        and solve_species_level_pi failure path fixes.
// V1_R5: No header changes; see ambient_rna_three_ap.cpp V1_R5 for species-level
//        bootstrap in bootstrap_amb_prof.
// V1_R6: No header changes; see ambient_rna_three_ap.cpp V1_R6 for orphan mass
//        redistribution in expand_species_prior_to_indiv.
// V1_R7: No _ap header changes. common.h gains species_sample_weight map and
//        get_weight() on PanelMetadata; load_panel_metadata adds fold_hybrid param.
// V1_R8: Added species_init_prof, species_init_used private members.
//        Added set_species_init(), set_indiv_init() public methods.
// V1_R9: Added has_species_counts, species_allelecounts, species_expfracs
//        private members. Added set_species_counts() public method.
// V1_R14: Historical permissive constrained-endpoint experiment; superseded.
// V1_R15: Added weighted-composition fallback, complete LOO component exclusion,
//        compile-status diagnostics, and explicit failed-cell diagnostic rows.
// V1_R16: Added safe unpenalized boundary validation, diagnostic-only fixed-r
//        fallback after 2-D failure, log-concave-prior enforcement, and a shared
//        complete-conditional weighted row selector for fitting/profile paths.
// V1_R17: Added deterministic same-model profile refinement and profiled
//        c-boundary assessment; superseded as control flow by V1_R18.
// V1_R20: Reverted the V1_R19 marginal self-column work in full, per project
//        rules 3.4a: the mechanism was refuted on 2026-07-27 by measuring 84
//        self rows against nalt/2 at a maximum absolute delta of 0.000000, so
//        the conditional-fraction self entry is already exact and there is
//        nothing to correct. Removed marginal_self_dosage,
//        marginal_ambient_column, marginal_ambient_column_raw,
//        unconditional_dosage_expectation, uncond_dosage_cache,
//        set_marginal_self_dosage and marginal_condf_fallbacks.
//        RETAINED BY EXCEPTION: final_loglik and final_loglik_valid, the
//        likelihood instrument. Justification: it is the only output that
//        distinguishes an estimator that cannot find the answer from an
//        objective whose maximum is not at the answer, and it located both
//        O76 and O77. Cost of retention: one output file per task and one
//        entry in the launcher copy-back list.
// V1_R19: Withdrawn. Marginal self-column construction, refuted before use.
// V1_R18: Added BFGS-independent free-model profile fitting, scale-invariant
//        prior-training eligibility, and separated BFGS/profile diagnostics.
// V1_R22: Added endpoint-compatible continuous source exclusion, source-mass
//        accounting, parent-axis and conditional-information diagnostics,
//        diagnostic provenance, and selector-limited r-C likelihood surfaces.
// ============================================================================
