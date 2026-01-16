#ifndef _CELLBOUNCER_DEMUX_PARALLEL_LLR_H
#define _CELLBOUNCER_DEMUX_PARALLEL_LLR_H
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
#include "common.h"
#include "demux_parallel_hts.h"  // For CellCounts

// Global verbose flag (set from main)
extern bool g_verbose;

// Global debug flag (set from main) - controls DEBUG spam to stderr
extern bool g_debug;

// ============================================================================
// DIAGNOSTIC STRUCTURES (NEW)
// ============================================================================

/**
 * Per-cell diagnostic information for downstream refinement
 */
struct CellDiagnostics {
    // Margin diagnostics (always computed)
    double min_margin;          // Winner's worst pairwise LLR
    int worst_competitor;       // Identity that gave min_margin
    int n_close;               // Count of identities within close_threshold of winner
    double total_depth;        // Total allele counts from main demux VCF
    
    // Ploidy diagnostics (only if --het_vcf provided)
    double het_balance_var;    // Variance of alt_frac at het sites (-1 if not computed)
    int n_het_sites;           // Number of het sites used (0 if not computed)
    double het_total_depth;    // Total depth at het sites (0 if not computed)
    
    // Method used for het_balance computation
    HetBalanceMethod het_method;
    
    CellDiagnostics() 
        : min_margin(0.0), worst_competitor(-1), n_close(0), total_depth(0.0),
          het_balance_var(-1.0), n_het_sites(0), het_total_depth(0.0),
          het_method(HetBalanceMethod::WELFORD) {}
};

/**
 * Runner-up identity information for quad detection
 */
struct RunnerUp {
    int identity;              // Runner-up identity index
    double llr_vs_winner;      // Direct LLR vs winner (negative = winner wins)
    double min_margin;         // Runner-up's own worst comparison
    
    RunnerUp() : identity(-1), llr_vs_winner(0.0), min_margin(0.0) {}
    RunnerUp(int id, double llr, double margin) 
        : identity(id), llr_vs_winner(llr), min_margin(margin) {}
};

// ============================================================================
// LLR TABLE CLASS
// ============================================================================

/**
 * Log likelihood ratio table for cell identity assignment.
 * 
 * Stores pairwise LLRs between all possible identities and provides
 * methods to iteratively eliminate unlikely identities.
 */
class llr_table{
    private:
        std::map<double, std::vector<std::pair<short, short> > > lookup_llr;
        std::map<double, std::vector<std::pair<short, short> > >::iterator it;
        std::vector<double> maxllr;
        std::vector<double> minllr;

    public:
        std::vector<bool> included;
        int n_indvs;
        
        llr_table(int x);
        ~llr_table();
        
        void print(std::string& bc_str, std::vector<std::string>& samples);
        void print_ranges(std::string& barcode, std::vector<std::string>& samples);
        void insert(short i1, short i2, double llr);
        void disallow(short i);
        void recalculate_minmax();
        bool del(int n_keep);
        void get_max(int& best_idx, double& best_llr);
        
        // NEW: Get min_margin for a specific identity
        double get_min_margin(int identity) const;
        
        // NEW: Get the lookup_llr map for diagnostic extraction
        const std::map<double, std::vector<std::pair<short, short> > >& get_lookup_llr() const { 
            return lookup_llr; 
        }
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
void assign_ids_parallel(
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
void assign_ids_parallel_with_diagnostics(
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
    robin_hood::unordered_map<unsigned long, std::vector<RunnerUp> >& runner_ups);

/**
 * Batch process with Welford or per-site het balance method
 */
void assign_ids_parallel_with_diagnostics_extended(
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
    robin_hood::unordered_map<unsigned long, std::vector<RunnerUp> >& runner_ups);

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

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

double adjust_p_err(double p, double e_r, double e_a);

void compute_k_comps(
    std::map<int, std::map<int, double> >& llrs, 
    llr_table& tab,
    std::vector<int>& ks,
    int n_samples,
    std::set<int>& allowed_assignments,
    double doublet_rate,
    std::map<int, double>* prior_weights);

#endif
