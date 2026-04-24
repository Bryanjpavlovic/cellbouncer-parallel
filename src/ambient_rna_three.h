// ============================================================================
// ambient_rna_three.h
// Three-component contamination model for tetraploid cells
//
// Version: V1_R3
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
// ============================================================================

#ifndef _CELLBOUNCER_AMBIENT_RNA_THREE_H
#define _CELLBOUNCER_AMBIENT_RNA_THREE_H

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
#include <mixtureDist/functions.h>
#include <optimML/brent.h>
#include <optimML/multivar.h>
#include "common.h"

class contamFinder3{
    private:

        // ---- Copies of external data structures ----

        robin_hood::unordered_map<unsigned long,
                std::map<std::pair<int, int>,
                    std::map<std::pair<int, int>,
                        std::pair<float, float> > > > indv_allelecounts;
        std::set<int> allowed_ids;
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
        double delta_thresh;
        int maxits;
        int nits;
        bool weighted;
        int num_threads;
        bool skip_reassign;

        // ---- Tetraploid-aware mode ----

        std::set<int> locked_identities;
        std::set<int> safe_singlets;
        bool tetraploid_aware;
        double min_signal_gap;
        bool amb_mu_available;
        bool ids_restricted;

        // ---- LOO ambient profiles ----

        bool use_loo;
        std::map<int, std::map<std::pair<int,int>,
            std::map<std::pair<int,int>, double>>> amb_mu_loo;
        void compute_loo_profiles();

        // ---- Fisher Information weighting ----

        bool use_fi_weight;

        // ---- User-specified prior override ----

        bool user_prior_set;
        double user_prior_mean;
        double user_prior_var;

        // ---- Internal methods ----

        void compile_data(robin_hood::unordered_map<unsigned long, int>& assn,
            robin_hood::unordered_map<unsigned long,
                std::map<std::pair<int, int>,
                std::map<std::pair<int, int>,
                std::pair<float, float> > > >& indv_allelecounts);

        void clear_data();

        void get_reads_expectations(int ident,
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

        // Ambient profile mixture proportions
        std::map<int, double> contam_prof;

        // ---- Constructor ----

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

        // ---- Configuration methods ----

        void set_error_rates(double e_r, double e_a);
        void model_other_species();
        void model_single_species();
        void set_mixprop_trials(int nt);
        void set_delta(double d);
        void set_maxiter(int i);
        void use_weights();
        void no_weights();
        void set_num_threads(int nt);
        void no_reassign();

        // Tetraploid-aware mode
        void set_locked_identities(const std::set<int>& locked);
        void set_safe_singlets(const std::set<int>& safe);
        void set_tetraploid_aware(bool enabled);
        void set_min_signal_gap(double gap);
        void set_ids_restricted(bool restricted);

        // LOO
        void set_use_loo(bool enabled);

        // FI weighting
        void set_fi_weight(bool enabled);

        // Prior override
        void set_user_prior(double mean, double var);

        // Bootstrap
        void bootstrap_amb_prof(int n_boots,
            std::map<int, double>& contam_prof_cont);

        // Compute total log likelihood
        double compute_ll();

        // Run the full estimation pipeline
        void fit();
};

#endif

// ============================================================================
// Revision History
// V1_R1: Initial three-component model header
// V1_R2: Fix ll init to -INFINITY, allele_ratio cleanup on reclassify,
//        num_threads default, split Welford prior, BFGS fallback counter
// V1_R3: Fix BFGS line search assertion failure (version bump to match .cpp)
// ============================================================================
