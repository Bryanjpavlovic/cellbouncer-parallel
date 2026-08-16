#ifndef _CELLBOUNCER_GENOTYPE_LLR_H
#define _CELLBOUNCER_GENOTYPE_LLR_H
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
#include <set>
#include <cstdlib>
#include <utility>
#include <limits>
#include "common.h"
#include "vcf_hts.h"  // For CellCounts

// Global verbose flag (set from main)
extern bool g_verbose;

// Global debug flag (set from main) - controls DEBUG spam to stderr
extern bool g_debug;

// ============================================================================
// DIAGNOSTIC STRUCTURES
// ============================================================================

/** Explicit state for a direct pairwise comparison in corrected diagnostics. */
enum class ComparisonState {
    PRESENT_NONZERO = 0,
    PRESENT_ZERO = 1,
    UNAVAILABLE = 2,
    PARTIAL_SUPPORT = 3,
    NOT_APPLICABLE = 4
};

const char* comparison_state_name(ComparisonState state);

struct PairwiseComparison {
    bool present;
    bool partial_support;
    double value;  // LLR(lhs versus rhs) when present

    PairwiseComparison()
        : present(false), partial_support(false),
          value(std::numeric_limits<double>::quiet_NaN()) {}
    PairwiseComparison(bool is_present, double comparison_value,
                       bool has_partial_support = false)
        : present(is_present), partial_support(has_partial_support),
          value(comparison_value) {}

    ComparisonState state() const {
        if (present) {
            return value == 0.0 ? ComparisonState::PRESENT_ZERO
                                : ComparisonState::PRESENT_NONZERO;
        }
        return partial_support ? ComparisonState::PARTIAL_SUPPORT
                               : ComparisonState::UNAVAILABLE;
    }
};

/** Per-cell diagnostic information for downstream refinement. */
struct CellDiagnostics {
    // Margin diagnostics. Missing numeric values are NaN and accompanied by an
    // explicit state; exact numeric zero remains a valid observed comparison.
    double min_margin;                    // Winner's maximin score
    double llr_vs_runnerup;               // Direct LLR(winner versus rank-1 runner)
    ComparisonState runnerup_comparison_state;
    int worst_competitor;                 // Identity that gave the worst present edge
    ComparisonState worst_comparison_state;
    int n_close;
    double total_depth;
    bool selection_resolved;
    int maximin_candidate;
    double maximin_score;
    std::vector<int> missing_comparison_alternatives;

    // Read-only ploidy diagnostics for the frozen accepted identity.
    double het_balance_var;
    int n_het_sites;
    double het_total_depth;
    bool het_diagnostic_available;
    HetBalanceMethod het_method;

    // Softmax summary of maximin margins. These are not calibrated Bayesian
    // posterior probabilities.
    double margin_softmax_score;
    double margin_entropy;

    // Non-mutating selection audit comparator: argmax(maxllr). It is not the
    // historical destructive elimination selector.
    int max_llr_comparator_winner;
    double max_llr_comparator_score;

    CellDiagnostics()
        : min_margin(std::numeric_limits<double>::quiet_NaN()),
          llr_vs_runnerup(std::numeric_limits<double>::quiet_NaN()),
          runnerup_comparison_state(ComparisonState::NOT_APPLICABLE),
          worst_competitor(-1),
          worst_comparison_state(ComparisonState::NOT_APPLICABLE),
          n_close(0),
          total_depth(0.0),
          selection_resolved(false),
          maximin_candidate(-1),
          maximin_score(std::numeric_limits<double>::quiet_NaN()),
          het_balance_var(std::numeric_limits<double>::quiet_NaN()),
          n_het_sites(0),
          het_total_depth(0.0),
          het_diagnostic_available(false),
          het_method(HetBalanceMethod::WELFORD),
          margin_softmax_score(std::numeric_limits<double>::quiet_NaN()),
          margin_entropy(std::numeric_limits<double>::quiet_NaN()),
          max_llr_comparator_winner(-1),
          max_llr_comparator_score(std::numeric_limits<double>::quiet_NaN()) {}
};

/** Runner-up identity information for quad detection. */
struct RunnerUp {
    int identity;
    double llr_vs_winner;  // LLR(runner versus winner) when comparison is present
    ComparisonState comparison_state;
    double min_margin;

    RunnerUp()
        : identity(-1),
          llr_vs_winner(std::numeric_limits<double>::quiet_NaN()),
          comparison_state(ComparisonState::NOT_APPLICABLE),
          min_margin(std::numeric_limits<double>::quiet_NaN()) {}
    RunnerUp(int id, double llr, ComparisonState state, double margin)
        : identity(id), llr_vs_winner(llr), comparison_state(state), min_margin(margin) {}
};

// ============================================================================
// LLR TABLE CLASS
// ============================================================================

/**
 * Pairwise log-likelihood-ratio table for identity assignment.
 *
 * Production selection is maximin: choose the retained identity with the
 * highest minimum present pairwise LLR. The table also retains an authoritative
 * signed pairwise map so diagnostic lookup can distinguish absence from zero.
 */
class llr_table {
    private:
        std::map<double, std::vector<std::pair<int, int> > > lookup_llr;
        std::map<double, std::vector<std::pair<int, int> > >::iterator it;
        std::map<std::pair<int, int>, double> pairwise_llr;
        std::set<std::pair<int, int> > pairwise_partial_support;
        std::vector<double> maxllr;
        std::vector<double> minllr;

    public:
        std::vector<bool> included;
        int n_indvs;

        explicit llr_table(int x);
        ~llr_table();

        void print(std::string& bc_str, std::vector<std::string>& samples);
        void print_ranges(std::string& barcode, std::vector<std::string>& samples);
        void insert(int i1, int i2, double llr);
        void mark_partial_support(int i1, int i2);
        void disallow(int i);
        void recalculate_minmax();
        bool del(int n_keep);
        void get_max(int& best_idx, double& best_llr) const;

        // Audit-only non-mutating argmax(maxllr) comparator.
        void get_max_by_max_llr_comparator(int& best_idx, double& best_maxllr) const;
        void get_max_by_maxllr(int& best_idx, double& best_maxllr) const;

        double get_min_margin(int identity) const;
        PairwiseComparison get_pairwise(int lhs, int rhs) const;
        std::vector<int> retained_identities() const;
        bool winner_has_complete_comparisons(
            int winner,
            std::vector<int>& missing_alternatives) const;

        const std::vector<double>& get_minllr() const { return minllr; }
        const std::vector<double>& get_maxllr() const { return maxllr; }
};

// ============================================================================
// ORIGINAL LLR FUNCTIONS (for legacy nested map structure)
// ============================================================================

bool populate_llr_table(std::map<std::pair<int, int>,
            std::map<std::pair<int, int>, 
                std::pair<float, float> > >& counts,
    std::map<int, std::map<int, double> >& llrs,
    llr_table& tab,
    int n_samples,
    std::set<int>& allowed_assignments,
    std::set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    std::map<int, double>* prior_weights=NULL,
    bool incl_contam=false,
    double contam_rate=0.0,
    double contam_rate_var=0.0,
    std::map<std::pair<int, int>, std::map<std::pair<int, int>, double> >* amb_fracs=NULL,
    int n_target=-1);


// Explicit serial ambient-model variants retained for the shared production
// ambient implementation. The per-identity entry point is the historical
// populate_llr_table model; the pairwise entry point uses the shared-soup term.
bool populate_llr_table_peridentity(
    std::map<std::pair<int, int>, std::map<std::pair<int, int>, std::pair<float, float> > >& counts,
    std::map<int, std::map<int, double> >& llrs,
    llr_table& tab,
    int n_samples,
    std::set<int>& allowed_assignments,
    std::set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    std::map<int, double>* prior_weights=NULL,
    bool incl_contam=false,
    double contam_rate=0.0,
    double contam_rate_var=0.0,
    std::map<std::pair<int, int>, std::map<std::pair<int, int>, double> >* amb_fracs=NULL,
    int n_target=-1);

bool populate_llr_table_pairwise(
    std::map<std::pair<int, int>, std::map<std::pair<int, int>, std::pair<float, float> > >& counts,
    std::map<int, std::map<int, double> >& llrs,
    llr_table& tab,
    int n_samples,
    std::set<int>& allowed_assignments,
    std::set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    std::map<int, double>* prior_weights,
    bool incl_contam,
    double contam_rate,
    double contam_rate_var,
    std::map<std::pair<int, int>, std::map<std::pair<int, int>, double> >* amb_fracs);

// ============================================================================
// NEW LLR FUNCTIONS (for optimized CellCounts structure)
// ============================================================================

/**
 * Populate LLR table from CellCounts structure (optimized format)
 * 
 * @param n_target Singlet pruning control:
 *                 -1 = auto (skip pruning if allowed_assignments provided, else use 10)
 *                  0 = never prune
 *                 >0 = prune to this many singlets
 */
bool populate_llr_table_optimized(
    const CellCounts& counts,
    std::map<int, std::map<int, double> >& llrs,
    llr_table& tab,
    int n_samples,
    std::set<int>& allowed_assignments,
    std::set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    std::map<int, double>* prior_weights=NULL,
    bool incl_contam=false,
    double contam_rate=0.0,
    double contam_rate_var=0.0,
    std::map<std::pair<int, int>, std::map<std::pair<int, int>, double> >* amb_fracs=NULL,
    int n_target=-1);

/**
 * Batch process multiple cells in parallel using CellCounts
 * 
 * @param n_target Singlet pruning control (see populate_llr_table_optimized)
 */
bool assign_ids_parallel(
    robin_hood::unordered_map<unsigned long, CellCounts>& cell_counts,
    std::vector<std::string>& samples,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, double>& assignments_llr,
    std::set<int>& allowed_assignments,
    std::set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    bool use_prior_weights,
    std::map<int, double>& prior_weights,
    int n_threads,
    int n_target = -1);

/**
 * NEW: Batch process with diagnostic collection
 * 
 * @param compute_diagnostics Whether to collect diagnostic info
 * @param n_runner_ups Number of runner-ups to extract per cell
 * @param close_threshold LLR threshold for counting "close" competitors
 * @param het_counts Optional het site counts (NULL if --het_vcf not provided)
 * @param diagnostics Output: per-cell diagnostics
 * @param runner_ups Output: per-cell runner-up lists
 */
bool assign_ids_parallel_with_diagnostics(
    robin_hood::unordered_map<unsigned long, CellCounts>& cell_counts,
    std::vector<std::string>& samples,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, double>& assignments_llr,
    std::set<int>& allowed_assignments,
    std::set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    bool use_prior_weights,
    std::map<int, double>& prior_weights,
    int n_threads,
    int n_target,
    // Diagnostic options
    bool compute_diagnostics,
    int n_runner_ups,
    double close_threshold,
    robin_hood::unordered_map<unsigned long, CellCounts>* het_counts,
    // Diagnostic outputs
    robin_hood::unordered_map<unsigned long, CellDiagnostics>& diagnostics,
    robin_hood::unordered_map<unsigned long, std::vector<RunnerUp> >& runner_ups,
    // Extended evidence (2A/2B/2C) - unused in this variant, accepted for call-site compatibility
    robin_hood::unordered_map<unsigned long, CellCounts>* atac_cell_counts = nullptr,
    const std::map<int, double>* identity_prior = nullptr,
    double z_doublet = 0.0,
    robin_hood::unordered_map<unsigned long, CellCounts>* species_counts_rna = nullptr,
    robin_hood::unordered_map<unsigned long, CellCounts>* species_counts_atac = nullptr);

/**
 * Batch process with Welford or per-site het balance method
 */
bool assign_ids_parallel_with_diagnostics_extended(
    robin_hood::unordered_map<unsigned long, CellCounts>& cell_counts,
    std::vector<std::string>& samples,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, double>& assignments_llr,
    std::set<int>& allowed_assignments,
    std::set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    bool use_prior_weights,
    std::map<int, double>& prior_weights,
    int n_threads,
    int n_target,
    bool compute_diagnostics,
    int n_runner_ups,
    double close_threshold,
    // Het data
    robin_hood::unordered_map<unsigned long, CellHetData>* het_data,
    const robin_hood::unordered_map<int, ChromSNPs>* het_snpdat,
    const std::vector<std::pair<int, int>>* idx_to_site,
    HetBalanceMethod het_method,
    int min_het_sites,
    double min_het_depth,
    // Outputs
    robin_hood::unordered_map<unsigned long, CellDiagnostics>& diagnostics,
    robin_hood::unordered_map<unsigned long, std::vector<RunnerUp> >& runner_ups,
    // Extended evidence (2A/2B/2C)
    robin_hood::unordered_map<unsigned long, CellCounts>* atac_cell_counts = nullptr,
    const std::map<int, double>* identity_prior = nullptr,
    double z_doublet = 0.0,
    robin_hood::unordered_map<unsigned long, CellCounts>* species_counts_rna = nullptr,
    robin_hood::unordered_map<unsigned long, CellCounts>* species_counts_atac = nullptr);

// ============================================================================
// DIAGNOSTIC EXTRACTION FUNCTIONS (NEW)
// ============================================================================

/**
 * Extract diagnostic information from LLR table before destruction
 * 
 * @param llrs The pairwise LLR map
 * @param tab The LLR table (after all processing)
 * @param winner The winning identity
 * @param winner_llr The winner's LLR
 * @param n_samples Number of samples
 * @param n_runner_ups Number of runner-ups to extract
 * @param close_threshold Threshold for counting "close" competitors
 * @param diag Output: diagnostic struct to populate
 * @param runners Output: runner-up list to populate
 */
void get_diagnostics_from_llrs(
    const std::map<int, std::map<int, double> >& llrs,
    const llr_table& tab,
    int winner,
    double winner_llr,
    int n_samples,
    int n_runner_ups,
    double close_threshold,
    CellDiagnostics& diag,
    std::vector<RunnerUp>& runners);

/**
 * Compute het site balance statistics for ploidy detection
 * DEPRECATED: Uses aggregated CellCounts which loses per-site info
 */
void compute_het_balance_stats(
    const CellCounts& het_counts,
    int assigned_id,
    int n_samples,
    CellDiagnostics& diag);

/**
 * Compute het balance using per-site data (PERSITE method)
 */
void compute_het_balance_persite(
    const HetSiteData& persite_data,
    const robin_hood::unordered_map<int, ChromSNPs>& het_snpdat,
    const std::vector<std::pair<int, int>>& idx_to_site,
    int assigned_id,
    int n_samples,
    double min_depth,
    int min_sites,
    CellDiagnostics& diag);

/**
 * Compute het balance using Welford stats (WELFORD method)
 */
void compute_het_balance_welford(
    const CellWelfordStats& welford_stats,
    int assigned_id,
    int n_samples,
    int min_sites,
    CellDiagnostics& diag);

/**
 * Compute total depth from main demux counts
 */
double compute_total_depth(const CellCounts& counts, int n_samples);

/**
 * Compute a softmax score and entropy from retained maximin margins.
 * These summaries are descriptive margin transforms, not calibrated Bayesian
 * posterior probabilities.
 *
 * @param llrs Retained for call-site compatibility; authoritative margins come from tab.
 * @param winner The winning identity index
 * @param n_samples Number of samples
 * @param diag Output: margin_softmax_score and margin_entropy populated
 */
void compute_margin_softmax_scores(
    const std::map<int, std::map<int, double> >& llrs,
    const llr_table& tab,
    int winner,
    int n_samples,
    CellDiagnostics& diag);

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

double adjust_p_err(double p, double e_r, double e_a);

void compute_k_comps(
    const std::map<int, std::map<int, double> >& llrs, 
    llr_table& tab,
    const std::vector<int>& ks,
    int n_samples,
    std::set<int>& allowed_assignments,
    double doublet_rate,
    std::map<int, double>* prior_weights);

#endif
