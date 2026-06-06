#ifndef _CELLBOUNCER_COMMON_H
#define _CELLBOUNCER_COMMON_H
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
#include <map>
#include <set>
#include <cstdlib>
#include <utility>
#include <htswrapper/bc.h>
/**
 * Contains functions used by more than one program in this
 * repository.
 */

using std::cout;
using std::endl;
using namespace std;

// Print help screen information about the parameter for
// library names that will be appended to cell barcodes
// (used by multiple programs)
void print_libname_help();

// Parse a file mapping cell barcodes to identities
void parse_barcode_map(std::string& fn, 
    std::map<unsigned long, std::string>& bc2hap,
    std::set<std::string>& barcode_groups,
    double llr_cutoff,
    bool keep_doublets);

// Utility functions for translating back and forth between
// identity indices and names, using a standard way of 
// representing doublet identities as integers
short hap_comb_to_idx(short, short, short);
std::pair<short, short> idx_to_hap_comb(short, short);
std::string idx2name(int, std::vector<std::string>&);

// Represent/print a triangular row-major distance matrix
void init_distmat(vector<vector<float> >& distmat, int dim);
void print_distmat(vector<vector<float> >& distmat);
void print_distmat_square(vector<vector<float> >& distmat);

// Transform a log likelihood ratio table of possible cell identities
// into a single identity and its log likelihood ratio
int collapse_llrs(std::map<int, std::map<int, double> >& llrs, double& llr_final);

// Check for over/underrepresentation of counts of specific 
// doublet types in a data set
double doublet_chisq(std::map<int, int>& idcounts, int n_samples);

// Trim the path off of a file name
std::string filename_nopath(std::string& filename);

// Log PDF of binomial distribution
double logbinom(double n, double k, double p);

bool file_exists(std::string name);
bool is_dir(std::string name);

// Approximate derivative using slope of neighboring points
void derivative(std::map<double, double>& hist, std::map<double, double>& result, int smooth);

// Find inflection point in a histogram
double find_knee(std::map<double, double>& hist, double min_frac_to_allow);

void fit_dirichlet(std::vector<double>& mle_fracs,
    std::vector<std::vector<double> >& dirichlet_bootstraps,
    std::vector<double>& conc_param_results,
    int nthreads=1);

// ============================================================================
// PANEL METADATA (species-to-individual mapping)
// ============================================================================

struct PanelMetadata {
    std::map<std::string, std::string> indiv_to_species;
    std::vector<std::string> species_list;                    // sorted, unique
    std::map<std::string, std::vector<int>> species_to_sample_indices;
    // species_to_sample_indices populated using the VCF sample list passed
    // to the loader; values are indices into that sample list.

    // Per-(species, sample_index) weight. Default 1.0 for normal individuals.
    // Hybrid individuals folded into parent species get 0.5 in each parent.
    // Key: (species_label, sample_index). Missing entries imply weight 1.0.
    std::map<std::pair<std::string, int>, double> species_sample_weight;

    std::vector<std::pair<std::string, std::string>> species_pairs;  // unordered pairs
    std::map<std::pair<std::string, std::string>, int> pair_to_index;
    int n_pairs;

    // Convenience: get weight for a (species, sample_index) pair
    double get_weight(const std::string& sp, int idx) const {
        auto key = std::make_pair(sp, idx);
        auto it = species_sample_weight.find(key);
        return (it != species_sample_weight.end()) ? it->second : 1.0;
    }
};

PanelMetadata load_panel_metadata(const std::string& filename,
                                  const std::vector<std::string>& vcf_samples,
                                  bool fold_hybrid = true);

#endif
