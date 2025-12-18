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
