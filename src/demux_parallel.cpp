#include <getopt.h>
#include <argp.h>
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
#include <float.h>
#include <chrono>
#include <omp.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <zlib.h>
#include <htswrapper/bc.h>
#include <htswrapper/bam.h>
#include <htswrapper/gzreader.h>
#include <mixtureDist/mixtureDist.h>
#include <mixtureDist/mixtureModel.h>
#include <mixtureDist/functions.h>
#include <optimML/multivar_ml.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "demux_vcf_io.h"
#include "demux_parallel_hts.h"
#include "demux_parallel_llr.h"

using std::cout;
using std::endl;
using namespace std;

// Version information
const string VERSION = "1.16";
const string VERSION_MESSAGE = "-I flag fix (recalculate_minmax), SNP end boundary fix, SNP iterator order fix, non-ACGT->N, exact V1 read count, fixed-point int64 accumulation, auto htslib_threads, timing output, base extraction debug stats, shared memory VCF support";
const string VERSION_NEW = "v1.16: Fixed-point int64, per-thread storage, chunking by SNP density OR read density (fixes chrM bottleneck)";

// Global verbose flag (defined in demux_vcf_llr.cpp)
extern bool g_verbose;

/**
 * Log likelihood function for computing error rates
 */
double ll_err(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p0 = data_d.at("exp");
    double e_r = params[0];
    double e_a = params[1];
    double p = p0 - p0*e_a + (1.0 - p0)*e_r;
    if (p <= 0){
        p = DBL_MIN*1e6;
    }
    else if (p >= 1){
        p = 1.0-DBL_MIN*1e6;
    }
    double ll = dbinom(n, k, p)/log2(exp(1));
    return ll;
}

void dll_err(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i, vector<double>& results){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double e_r = params[0];
    double e_a = params[1];
    double p0 = data_d.at("exp");
    double p = p0 - p0*e_a + (1.0 - p0)*e_r;
    if (p <= 0){
        p = DBL_MIN*1e6;
    }
    else if (p >= 1.0){
        p = 1.0-DBL_MIN*1e6;
    }
    
    double dy_dp = (k-n*p)/(p-p*p);
    double dp_de_a = -p0;
    double dp_de_r = 1.0 - p0;
    results[0] = dy_dp * dp_de_r;
    results[1] = dy_dp * dp_de_a;
}

/**
 * Re-infer error rates from assignments using optimized structure
 */
pair<double, double> infer_error_rates_optimized(
    robin_hood::unordered_map<unsigned long, CellCounts>& cell_counts,
    int n_samples,
    robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    double error_ref,
    double error_alt,
    double error_sigma){
   
    vector<double> n;
    vector<double> k;
    vector<double> expected;
    vector<double> weights_llr;

    int doublet_points = 0, singlet_points = 0;
    for (auto& a : assn){
        double weight = assn_llr[a.first];
        const CellCounts& counts = cell_counts[a.first];
        
        bool is_combo = a.second >= n_samples;
        pair<int, int> combo;
        if (is_combo){
            combo = idx_to_hap_comb(a.second, n_samples);
        }
        
        // Debug: show assignment details for each cell
        fprintf(stderr, "DEBUG cell: assn=%d n_samples=%d is_combo=%d", 
                a.second, n_samples, is_combo);
        if (is_combo) {
            fprintf(stderr, " combo=(%d,%d)", combo.first, combo.second);
        }
        fprintf(stderr, "\n");
        
        if (is_combo){
            // Doublet - use pairwise counts for the two component individuals
            int i = combo.first;
            int j = combo.second;
            
            for (int nalt_i = 0; nalt_i < 3; ++nalt_i){
                for (int nalt_j = 0; nalt_j < 3; ++nalt_j){
                    auto pair_counts = counts.get(i, nalt_i, j, nalt_j);
                    float ref_count = pair_counts.first;
                    float alt_count = pair_counts.second;
                    
                    if (ref_count + alt_count > 0){
                        double this_expected = (double)(nalt_i + nalt_j) / 4.0;
                        expected.push_back(this_expected);
                        n.push_back(ref_count + alt_count);
                        k.push_back(alt_count);
                        weights_llr.push_back(weight);
                        doublet_points++;
                    }
                }
            }
        }
        else{
            // Singlet - use total counts for this individual
            for (int nalt = 0; nalt < 3; ++nalt){
                auto total = counts.get_total(a.second, nalt);
                float ref_count = total.first;
                float alt_count = total.second;
                
                if (ref_count + alt_count > 0){
                    double this_expected = (double)nalt / 2.0;
                    expected.push_back(this_expected);
                    n.push_back(ref_count + alt_count);
                    k.push_back(alt_count);
                    weights_llr.push_back(weight);
                    singlet_points++;
                }
            }
        }
    }
    fprintf(stderr, "DEBUG: doublet_points=%d singlet_points=%d\n", doublet_points, singlet_points);

    if (n.empty()){
        return make_pair(error_ref, error_alt);
    }

    // Diagnostic: print totals going into optimizer
    double total_n = 0, total_k = 0, total_weight = 0;
    std::map<double, std::pair<double, double>> by_expected;  // expected -> (n, k)
    for (size_t i = 0; i < n.size(); i++) {
        total_n += n[i];
        total_k += k[i];
        total_weight += weights_llr[i];
        by_expected[expected[i]].first += n[i];
        by_expected[expected[i]].second += k[i];
    }
    fprintf(stderr, "DEBUG infer_error_rates_optimized:\n");
    fprintf(stderr, "  data_points=%lu total_n=%.2f total_k=%.2f total_weight=%.2f\n", 
            n.size(), total_n, total_k, total_weight);
    fprintf(stderr, "  overall_alt_frac=%.6f\n", total_k / total_n);
    for (auto& e : by_expected) {
        fprintf(stderr, "  expected=%.2f: n=%.2f k=%.2f alt_frac=%.6f\n",
                e.first, e.second.first, e.second.second, 
                e.second.second / e.second.first);
    }

    optimML::multivar_ml_solver solver({error_ref, error_alt}, ll_err, dll_err);
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("exp", expected);
    solver.add_weights(weights_llr);
    solver.constrain_01(0);
    solver.constrain_01(1);
    solver.add_normal_prior(0, error_ref, error_sigma, 0.0, 1.0);
    solver.add_normal_prior(1, error_alt, error_sigma, 0.0, 1.0);
    
    solver.set_silent(true);

    double sigma_curr = error_sigma;
    bool success = false;
    while (!success){
        success = true;
        try{ 
            solver.solve();
        } 
        catch (const int& err){
            if (err == optimML::OPTIMML_MATH_ERR){
                sigma_curr *= 0.5;
                fprintf(stderr, "Decreasing prior sd to %f...\n", sigma_curr);
                solver.set_prior_param(0, "sigma", sigma_curr);
                solver.set_prior_param(1, "sigma", sigma_curr);
                solver.set_param(0, error_ref);
                solver.set_param(1, error_alt);        
                success = false;
            }
            else{
                fprintf(stderr, "Unknown error encountered while inferring error rates\n");
                return make_pair(error_ref, error_alt);
            }
        }
    }
    return make_pair(solver.results[0], solver.results[1]);
}

/**
 * Mann-Whitney U test for comparing LLR distributions between identities.
 * Returns the p-value for testing whether id1's LLRs are lower than id2's (or all others if id2 == -1).
 */
double mannwhitney_llr(int id1, 
    int id2, 
    robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr){
    
    double n1 = 0.0;
    double n2 = 0.0;

    vector<pair<double, int> > llrs;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){
        if (assn_llr[a->first] > 0){
            if (a->second == id1){
                llrs.push_back(make_pair(assn_llr[a->first], id1));
                n1++;
            }
            else if (id2 == -1 || a->second == id2){
                llrs.push_back(make_pair(assn_llr[a->first], id2));
                n2++;
            }
        }
    }
    
    if (n1 == 0){
        return 0.0;
    }
    else if (n2 == 0){
        return 1.0;
    }

    sort(llrs.begin(), llrs.end());
    
    // Assign ranks
    vector<double> ranks;
    int rank = 1;
    for (int i = 0; i < llrs.size(); ++i){
        ranks.push_back((double)rank);
        rank++;
    }

    // Deal with ties in ranks
    double prevllr = 0.0;
    double ranksum = 0.0;
    int ranktot = 0;
    bool ties = false;

    vector<double> nties;

    for (int i = 0; i < llrs.size(); ++i){
        if (llrs[i].first == prevllr){
            ranksum += ranks[i];
            ranktot++;
        }
        else{
            if (ranktot > 1){
                nties.push_back((double)ranktot);
                ties = true;
                double rankmean = ranksum / (double)ranktot;
                for (int j = i - 1; j >= i - 1 - (ranktot-1); --j){
                    ranks[j] = rankmean;
                }
                ranksum = 0.0;
                ranktot = 0;
            }
            else{
                if (i > 0){
                    nties.push_back(0);
                    ranks[i-1] = ranksum;
                }
                ranksum = 0;
                ranktot = 0;
            }
            ranksum += ranks[i];
            ranktot++;

        }
        prevllr = llrs[i].first;    
    }
    // Handle last one.
    if (ranktot > 1){
        ties = true;
        nties.push_back(ranktot);
        double rankmean = ranksum / (double)ranktot;
        for (int j = llrs.size()-1; j >= (int)llrs.size()-1 - (ranktot-1); --j){
            ranks[j] = rankmean;
        }
    }
    else{
        nties.push_back(0);
        ranks[llrs.size()-1] = ranksum;
    }
    
    double sum_id1 = 0;
    double sum_id2 = 0;
    for (int i = 0; i < ranks.size(); ++i){
        if (llrs[i].second == id1){
            sum_id1 += ranks[i];
        }
        else{
            sum_id2 += ranks[i];
        }
    } 
    
    double U1 = sum_id1 - (n1*(n1+1))/2.0;
    double U2 = sum_id2 - (n2*(n2+1))/2.0;
    double m_u = (n1*n2)/2.0;
    double sigma_u = sqrt((n1*n2*(n1+n2 +1))/(12.0));
    if (ties){
        double tsum = 0;
        for (int i = 0; i < nties.size(); ++i){
            if (nties[i] > 0){
                tsum += (pow(nties[i], 3) - nties[i]);
            }
        }
        double term1 = (n1*n2*(n1+n2+1))/12;
        double term2 = (n1*n2*tsum)/(12*(n1+n2)*(n1+n2-1));
        sigma_u = sqrt(term1-term2);
    }
    
    if (n1 < 3 || n2 < 3){
        // Can't get a reliable sigma here.
        if (U1 < m_u){
            return 0.0;
        }
        else{
            return 1.0;
        }
    }

    // Test whether id1 < id2
    double p = pnorm(U1, m_u, sigma_u);
    return p;
}

/**
 * Does some QC on assignments - for each ID in the assignment file,
 * checks for significantly lower LLR distribution than the rest of
 * the data set. Also checks for significantly lower numbers of cells.
 */
void id_qc(robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    map<int, double>& pois_p,
    map<int, double>& mannwhitney_p){
    
    // Get total num cells for each ID 
    map<int, int> idsizes;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){
        if (assn_llr[a->first] > 0){
            if (idsizes.count(a->second) == 0){
                idsizes.insert(make_pair(a->second, 0));
            }
            idsizes[a->second]++;
        }
    }

    for (map<int, int>::iterator ids = idsizes.begin(); ids != idsizes.end(); ++ids){
        double mean_othersize = 0.0;
        double mean_othersize_tot = 0.0;
        for (map<int, int>::iterator ids2 = idsizes.begin(); ids2 != idsizes.end(); ++ids2){
            if (ids2->first != ids->first){
                mean_othersize += (double)ids2->second;
                mean_othersize_tot++;
            }
            
        }
        double p1 = ppois(ids->second, mean_othersize_tot);
        pois_p.insert(make_pair(ids->first, p1));
        double p2 = mannwhitney_llr(ids->first, -1, assn, assn_llr);
        mannwhitney_p.insert(make_pair(ids->first, p2));
    } 
}

/**
 * Filter identities based on user-specified doublet ID file.
 * Returns true if any identities were removed.
 */
bool filter_identities(robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr, 
    int n_samples,
    set<int>& allowed_ids, 
    set<int>& allowed_ids2){
    
    bool removed_ids = false;
    
    // Get number of cells per ID
    map<int, int> cells_per_id;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){
        if (cells_per_id.count(a->second) == 0){
            cells_per_id.insert(make_pair(a->second, 0));
        }
        cells_per_id[a->second]++;
    }

    // Get a list of identities not in the original filtered list, but added because
    // they're singlet components of allowed doublets
    set<int> candidate_ids;
    for (set<int>::iterator a = allowed_ids.begin(); a != allowed_ids.end(); ++a){
        if (allowed_ids2.find(*a) == allowed_ids2.end()){
            
            // Check whether this individual has significantly fewer cells than 
            // all "parent" individuals (doublets including this individual)
            vector<double> p_ncell;
            
            // Check whether cells assigned to this individual have significantly
            // lower LLRs than cells assigned to "parent" individuals
            vector<double> p_llr;
            
            for (set<int>::iterator a2 = allowed_ids2.begin(); a2 != allowed_ids2.end(); ++a2){
                if (*a2 >= n_samples){
                    pair<int, int> combo = idx_to_hap_comb(*a2, n_samples);
                    if (combo.first == *a || combo.second == *a){
                        double p1 = ppois((double)cells_per_id[*a], (double)cells_per_id[*a2]);
                        double p2 = mannwhitney_llr(*a, *a2, assn, assn_llr);
                        p_ncell.push_back(p1);
                        p_llr.push_back(p2);
                    }
                }
            }
            bool all_p_signif = true;
            for (int i = 0; i < p_ncell.size(); ++i){
                if (p_ncell[i] > 0.05 || p_llr[i] > 0.05){
                    all_p_signif = false;  // At least one comparison not significant
                }
            }
            if (all_p_signif){
                // Remove a.
                removed_ids = true;
            }
            else{
                // Safe to keep.
                allowed_ids2.insert(*a);
            }
        }
    }
    return removed_ids;
}

/**
 * Dump cell counts from optimized CellCounts structure to file.
 * Format matches original: cell_barcode indv1 nalt1 indv2 nalt2 ref_count alt_count
 */
void dump_cellcounts_optimized(gzFile& out_cell,
    robin_hood::unordered_map<unsigned long, CellCounts>& cell_counts,
    int n_samples){
    
    char linebuf[1024];

    for (auto& cell : cell_counts){
        unsigned long bc = cell.first;
        const CellCounts& counts = cell.second;
        
        // Output totals (indv2 = -1, nalt2 = -1)
        for (int indv = 0; indv < n_samples; ++indv){
            for (int nalt = 0; nalt < 3; ++nalt){
                auto total = counts.get_total(indv, nalt);
                if (total.first > 0 || total.second > 0){
                    sprintf(&linebuf[0], "%lu\t%d\t%d\t%d\t%d\t%f\t%f\n", 
                        bc, indv, nalt, -1, -1, total.first, total.second);
                    gzwrite(out_cell, &linebuf[0], strlen(linebuf));
                }
            }
        }
        
        // Output pairwise counts
        for (int indv1 = 0; indv1 < n_samples; ++indv1){
            for (int nalt1 = 0; nalt1 < 3; ++nalt1){
                for (int indv2 = indv1 + 1; indv2 < n_samples; ++indv2){
                    for (int nalt2 = 0; nalt2 < 3; ++nalt2){
                        auto pair_counts = counts.get(indv1, nalt1, indv2, nalt2);
                        if (pair_counts.first > 0 || pair_counts.second > 0){
                            sprintf(&linebuf[0], "%lu\t%d\t%d\t%d\t%d\t%f\t%f\n", 
                                bc, indv1, nalt1, indv2, nalt2, 
                                pair_counts.first, pair_counts.second);
                            gzwrite(out_cell, &linebuf[0], strlen(linebuf));
                        }
                    }
                }
            }
        }
    } 
}

/**
 * Print help message
 */
void help(int code){
    fprintf(stderr, "demux_parallel version %s\n", VERSION.c_str());
    fprintf(stderr, "%s\n\n", VERSION_MESSAGE.c_str());
    fprintf(stderr, "demux_parallel [OPTIONS]\n");
    fprintf(stderr, "Parallel version of demux_vcf - demultiplexes cells based on genotype data.\n");
    fprintf(stderr, "\n[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --bam -b The BAM file of interest\n");
    fprintf(stderr, "    --vcf -v A VCF/BCF file containing genotypes\n");
    fprintf(stderr, "    --output_prefix -o Base name for output files\n");
    fprintf(stderr, "\n===== OPTIONAL =====\n");
    fprintf(stderr, "    --threads -t Number of threads for parallel processing [auto]\n");
    // hts_threads now auto-calculated based on thread count
    fprintf(stderr, "    --barcodes -B A file listing cell barcodes (one per line)\n");
    fprintf(stderr, "    --ids -i A file listing allowed individual IDs (singlets)\n");
    fprintf(stderr, "    --ids_doublet -I A file listing allowed doublet combinations\n");
    fprintf(stderr, "    --qual -q Minimum QUAL score for variants [50]\n");
    fprintf(stderr, "    --doublet_rate -D Prior probability of doublet [0.5]\n");
    fprintf(stderr, "    --error_ref -e Prior error rate for ref allele [0.005]\n");
    fprintf(stderr, "    --error_alt -E Prior error rate for alt allele [0.005]\n");
    fprintf(stderr, "    --error_sigma -s Sigma for error rate prior [0.1]\n");
    fprintf(stderr, "    --no_preload -P Disable VCF preloading (use with --shared_vcf for low memory)\n");
    fprintf(stderr, "    --shared_vcf -S Name of shared memory VCF to attach to\n");
    fprintf(stderr, "    --vcf_chroms -c File listing chromosomes to use from VCF\n");
    fprintf(stderr, "    --libname -n Library name to append to barcodes\n");
    fprintf(stderr, "    --cellranger -C Format barcodes for CellRanger\n");
    fprintf(stderr, "    --seurat -R Format barcodes for Seurat\n");
    fprintf(stderr, "    --underscore -U Use underscore instead of hyphen for libname\n");
    fprintf(stderr, "    --disable_conditional -f Disable computing conditional match fractions\n");
    fprintf(stderr, "    --n_target -N Max singlets to keep before doublet eval [-1=auto, 0=no limit]\n");
    fprintf(stderr, "    --verbose -V Enable verbose output\n");
    fprintf(stderr, "    --help -h Display this message and exit\n");
    exit(code);
}

// Timing helper
void print_elapsed(const std::chrono::steady_clock::time_point& start, const char* step) {
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
    int hours = elapsed / 3600;
    int mins = (elapsed % 3600) / 60;
    int secs = elapsed % 60;
    fprintf(stderr, "[%02d:%02d:%02d] %s\n", hours, mins, secs, step);
}

int main(int argc, char *argv[]) {    
    
    // Start timing
    auto start_time = std::chrono::steady_clock::now();
    
    // Print version info
    fprintf(stderr, "demux_parallel version %s\n", VERSION.c_str());
    fprintf(stderr, "%s\n", VERSION_MESSAGE.c_str());
    fprintf(stderr, "New: %s\n", VERSION_NEW.c_str());
   
    static struct option long_options[] = {
       {"bam", required_argument, 0, 'b'},
       {"vcf", required_argument, 0, 'v'},
       {"output_prefix", required_argument, 0, 'o'},
       {"barcodes", required_argument, 0, 'B'},
       {"doublet_rate", required_argument, 0, 'D'},
       {"ids", required_argument, 0, 'i'},
       {"ids_doublet", required_argument, 0, 'I'},
       {"qual", required_argument, 0, 'q'},
       {"libname", required_argument, 0, 'n'},
       {"cellranger", no_argument, 0, 'C'},
       {"seurat", no_argument, 0, 'R'},
       {"underscore", no_argument, 0, 'U'},
       {"error_ref", required_argument, 0, 'e'},
       {"error_alt", required_argument, 0, 'E'},
       {"error_sigma", required_argument, 0, 's'},
       {"disable_conditional", no_argument, 0, 'f'},
       {"no_preload", no_argument, 0, 'P'},
       {"vcf_chroms", required_argument, 0, 'c'},
       {"threads", required_argument, 0, 't'},
       {"shared_vcf", required_argument, 0, 'S'},
       {"n_target", required_argument, 0, 'N'},
       {"verbose", no_argument, 0, 'V'},
       {"help", no_argument, 0, 'h'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string bamfile = "";
    string vcf_file = "";
    bool cell_barcode = false;
    string cell_barcode_file = "";
    string output_prefix = "";
    int vq = 50;
    string idfile;
    string idfile_doublet;
    bool idfile_given = false;
    bool idfile_doublet_given = false;
    double doublet_rate = 0.5;
    string barcode_group = "";
    double error_ref = 0.005;
    double error_alt = 0.005;
    double error_sigma = 0.1;
    
    bool cellranger = false;
    bool seurat = false;
    bool underscore = false;

    bool disable_conditional = false;

    bool no_preload = false;
    string vcf_chroms_file = "";
    bool vcf_chroms_given = false;
    
    // New parallel processing options
    int n_threads = omp_get_num_procs();
    int htslib_threads = 0;  // 0 = auto-calculate after parsing args
    string shared_vcf_name = "";
    int n_target = -1;    // -1=auto, 0=no prune, >0=use value
    bool verbose = false;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:v:o:B:i:I:q:D:n:e:E:s:c:t:S:N:fPCRUVh", 
        long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                break;
            case 'h':
                help(0);
                break;
            case 'b':
                bamfile = optarg;
                break;
            case 'v':
                vcf_file = optarg;
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 'n':
                barcode_group = optarg;
                break;
            case 'C':
                cellranger = true;
                break;
            case 'R':
                seurat = true;
                break;
            case 'U':
                underscore = true;
                break;
            case 'B':
                cell_barcode = true;
                cell_barcode_file = optarg;
                break;
            case 'D':
                doublet_rate = atof(optarg);
                break;
            case 'i':
                idfile_given = true;
                idfile = optarg;
                break;
            case 'I':
                idfile_doublet_given = true;
                idfile_doublet = optarg;
                break;
            case 'q':
                vq = atoi(optarg);
                break;
            case 'e':
                error_ref = atof(optarg);
                break;
            case 'E':
                error_alt = atof(optarg);
                break;
            case 's':
                error_sigma = atof(optarg);
                break;
            case 'f':
                disable_conditional = true;
                break;
            case 'P':
                no_preload = true;
                break;
            case 'c':
                vcf_chroms_given = true;
                vcf_chroms_file = optarg;
                break;
            case 't':
                n_threads = atoi(optarg);
                break;
            case 'S':
                shared_vcf_name = optarg;
                break;
            case 'N':
                n_target = atoi(optarg);
                break;
            case 'V':
                verbose = true;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Set global verbose flag
    g_verbose = verbose;
    
    // Error check arguments
    if (vq < 0){
        fprintf(stderr, "ERROR: variant quality must be >= 0\n");
        exit(1);
    }
    if (output_prefix.length() == 0){
        fprintf(stderr, "ERROR: output_prefix/-o required\n");
        exit(1);
    }
    if (is_dir(output_prefix)){
        fprintf(stderr, "ERROR: output_prefix %s is a directory\n", output_prefix.c_str());
        exit(1);
    }
    if (doublet_rate < 0 || doublet_rate > 1){
        fprintf(stderr, "ERROR: doublet rate must be between 0 and 1\n");
        exit(1);
    }
    if (idfile_given && idfile_doublet_given){
        fprintf(stderr, "ERROR: only one of -i/-I is allowed\n");
        exit(1);
    }
    if (n_threads < 1){
        n_threads = 1;
    }
    
    // Auto-calculate htslib threads per reader
    // With many workers (>=8), use 1 thread to avoid oversubscription
    // With few workers (<8), use 2 threads since decompression may bottleneck
    if (htslib_threads == 0){
        htslib_threads = (n_threads >= 8) ? 1 : 2;
    }
    
    fprintf(stderr, "demux_parallel: Parallel VCF-based demultiplexing\n");
    fprintf(stderr, "Using %d threads (%d htslib threads per reader)\n", n_threads, htslib_threads);
    
    // Init BAM reader for header info
    bam_reader reader;
    
    // Check if we can load counts from previous run
    bool load_counts = false;
    string countsfilename = output_prefix + ".counts";
    if (file_exists(countsfilename)){
        load_counts = true;
        fprintf(stderr, "Found existing counts file: %s\n", countsfilename.c_str());
    }
    else{
        if (bamfile.length() == 0){
            fprintf(stderr, "ERROR: bam file (--bam) required\n");
            exit(1);
        }    
        reader.set_file(bamfile);
    }

    // Load sample names
    vector<string> samples;
    bool samples_from_vcf = false;

    if (load_counts){
        string samplesfile = output_prefix + ".samples";
        if (file_exists(samplesfile)){
           load_samples(samplesfile, samples); 
        }
        else{
            if (vcf_file == ""){
                fprintf(stderr, "ERROR: vcf file is required\n");
                exit(1);
            }
            read_vcf_samples(vcf_file, samples);         
        }
    }
    else{
        if (vcf_file == "" && shared_vcf_name == ""){
            fprintf(stderr, "ERROR: vcf file is required\n");
            exit(1);
        }
        if (vcf_file != ""){
            read_vcf_samples(vcf_file, samples);
        }
        samples_from_vcf = true;
    }
    
    // Parse allowed IDs - deferred until after VCF loading if using shared memory
    set<int> allowed_ids;
    set<int> allowed_ids2;
    
    // Process samples now if not using shared VCF (otherwise defer until after attach)
    if (shared_vcf_name.length() == 0){
        fprintf(stderr, "Number of individuals in VCF: %lu\n", samples.size());
        
        if (idfile_given){
            parse_idfile(idfile, samples, allowed_ids, allowed_ids2, true);
            if (allowed_ids.size() == 0){
                fprintf(stderr, "No valid individual names found in %s; allowing all\n", idfile.c_str());
            }
        }
        if (idfile_doublet_given){
            parse_idfile(idfile_doublet, samples, allowed_ids, allowed_ids2, false);
            if (allowed_ids.size() == 0){
                fprintf(stderr, "No valid individual names found in %s; allowing all\n", idfile_doublet.c_str());
            }
        }
        
        if (samples_from_vcf){
            string samplesfile = output_prefix + ".samples"; 
            write_samples(samplesfile, samples);
        }
    }
    
    // Load cell barcodes
    set<unsigned long> cell_barcodes;
    if (cell_barcode){
        parse_barcode_file(cell_barcode_file, cell_barcodes);
        if (cell_barcodes.size() == 0){
            fprintf(stderr, "ERROR reading cell barcode list\n");
            exit(1);
        }
        fprintf(stderr, "Loaded %lu cell barcodes\n", cell_barcodes.size());
    }

    // Main data structures
    robin_hood::unordered_map<unsigned long, CellCounts> cell_counts;
    
    map<pair<int, int>, map<int, float> > conditional_match_fracs;
    map<pair<int, int>, map<int, float> > conditional_match_tots;
    
    if (load_counts){
        // Load from file - need to convert format
        fprintf(stderr, "Loading counts from file not yet supported in v3 format\n");
        fprintf(stderr, "Please delete %s and re-run\n", countsfilename.c_str());
        exit(1);
    }
    else{
        // Determine chromosomes to process
        set<string> chroms_bam;
        set<string> chroms_vcf;
        set<string> chroms_to_process;
        
        fprintf(stderr, "Identifying shared chromosomes between BAM and VCF...\n");
        get_bam_chroms(reader, chroms_bam);
        
        if (shared_vcf_name.length() > 0){
            // Using shared memory VCF - chroms will be determined after loading
        }
        else{
            get_vcf_chroms(vcf_file, chroms_vcf);
            
            for (auto& c : chroms_bam){
                if (chroms_vcf.find(c) != chroms_vcf.end()){
                    chroms_to_process.insert(c);
                }
            }
            
            if (vcf_chroms_given){
                set<string> user_chroms;
                ifstream chromfile(vcf_chroms_file.c_str());
                string chrom;
                while (chromfile >> chrom){
                    user_chroms.insert(chrom);
                }
                set<string> filtered_chroms;
                for (auto& c : chroms_to_process){
                    if (user_chroms.find(c) != user_chroms.end()){
                        filtered_chroms.insert(c);
                    }
                }
                chroms_to_process = filtered_chroms;
            }
            
            fprintf(stderr, "Found %lu chromosomes in BAM, %lu in VCF, %lu shared\n",
                chroms_bam.size(), chroms_vcf.size(), chroms_to_process.size());
        }
        
        // Load VCF data
        print_elapsed(start_time, "Starting VCF loading...");
        robin_hood::unordered_map<int, ChromSNPs> snpdat_optimized;
        map<string, int> seq2tid = reader.get_seq2tid();
        
        if (shared_vcf_name.length() > 0){
            // Attach to shared memory VCF
            fprintf(stderr, "Attaching to shared VCF: %s\n", shared_vcf_name.c_str());
            if (!attach_shared_vcf(shared_vcf_name, snpdat_optimized, samples)){
                fprintf(stderr, "ERROR: Could not attach to shared VCF\n");
                exit(1);
            }
            
            // Deferred samples processing (samples now populated from shared memory)
            fprintf(stderr, "Number of individuals in VCF: %lu\n", samples.size());
            
            if (idfile_given){
                parse_idfile(idfile, samples, allowed_ids, allowed_ids2, true);
                if (allowed_ids.size() == 0){
                    fprintf(stderr, "No valid individual names found in %s; allowing all\n", idfile.c_str());
                }
            }
            if (idfile_doublet_given){
                parse_idfile(idfile_doublet, samples, allowed_ids, allowed_ids2, false);
                if (allowed_ids.size() == 0){
                    fprintf(stderr, "No valid individual names found in %s; allowing all\n", idfile_doublet.c_str());
                }
            }
            
            // Write samples file
            string samplesfile = output_prefix + ".samples"; 
            write_samples(samplesfile, samples);
        }
        else{
            // Load VCF into memory (default behavior)
            // Use --no_preload to disable if memory is limited
            if (!no_preload){
                print_elapsed(start_time, "Loading VCF data into memory...");
                fprintf(stderr, "Loading VCF data into memory...\n");
                int nloaded = read_vcf_chroms_optimized(vcf_file, chroms_to_process, 
                    seq2tid, snpdat_optimized, vq);
                print_elapsed(start_time, "VCF loading complete");
                fprintf(stderr, "Loaded %d SNPs from %lu chromosomes\n", 
                    nloaded, chroms_to_process.size());
            }
            else{
                fprintf(stderr, "VCF preloading disabled. Use --shared_vcf for large datasets.\n");
                exit(1);
            }
        }
        
        // Count alleles
        print_elapsed(start_time, "Starting allele counting...");
        fprintf(stderr, "Counting alleles in BAM file...\n");
        
        if (n_threads > 1){
            // Parallel processing
            robin_hood::unordered_map<unsigned long, AlignedCellCounts> parallel_counts;
        
            count_alleles_parallel(bamfile, snpdat_optimized, parallel_counts,
                cell_barcodes, samples.size(), n_threads, htslib_threads);
        
            // Finalize counts
            finalize_parallel_counts(parallel_counts, cell_counts);
        }
        else{
            // Single-threaded fallback
            count_alleles_single_threaded(bamfile, snpdat_optimized, cell_counts,
                cell_barcodes, samples.size(), conditional_match_fracs,
                conditional_match_tots, !disable_conditional);
        }
    
        // Compute conditional match fractions if parallel (didn't do it during counting)
        if (n_threads > 1 && !disable_conditional){
            fprintf(stderr, "Computing conditional match fractions...\n");
            for (auto& kv : snpdat_optimized){
                get_conditional_match_fracs_chrom_optimized(kv.second,
                    conditional_match_fracs, conditional_match_tots, samples.size());
            }
        }
    
        if (!disable_conditional){
            conditional_match_fracs_normalize(conditional_match_fracs,
                conditional_match_tots, samples.size());
        
            string outname = output_prefix + ".condf";
            FILE* outf = fopen(outname.c_str(), "w");
            dump_exp_fracs(outf, conditional_match_fracs);
            fclose(outf);
        }
    } // end else (not load_counts)
    
    // Write counts to disk
    print_elapsed(start_time, "Writing allele counts to disk...");
    {
        string fname = output_prefix + ".counts";
        gzFile outf = gzopen(fname.c_str(), "w");
        fprintf(stderr, "Writing allele counts to disk...\n");
        dump_cellcounts_optimized(outf, cell_counts, samples.size());
        gzclose(outf);
        fprintf(stderr, "Done writing counts\n");
    }
        
    // Assign identities
    print_elapsed(start_time, "Starting identity assignment (round 1)...");
    robin_hood::unordered_map<unsigned long, int> assn;
    robin_hood::unordered_map<unsigned long, double> assn_llr;
    
    fprintf(stderr, "Finding likeliest identities of cells...\n");
    map<int, double> prior_weights;
    
    // First round of assignments
    // NOTE: allowed_ids contains all singlets + doublets (including singlet components)
    //       allowed_ids2 contains ONLY what was in the original ID file
    // When -I is used with doublet combinations, allowed_ids2 has only doublets,
    // which causes singlets to be disallowed in the final assignment step.
    assign_ids_parallel(cell_counts, samples, assn, assn_llr,
        allowed_ids, allowed_ids2, doublet_rate, error_ref, error_alt,
        false, prior_weights, n_threads, n_target);
    
    robin_hood::unordered_map<unsigned long, int> assncpy = assn;
    
    // Re-estimate error rates
    print_elapsed(start_time, "Estimating error rates...");
    fprintf(stderr, "Finding likeliest error rates...\n");
    pair<double, double> err_new = infer_error_rates_optimized(cell_counts, samples.size(),
        assn, assn_llr, error_ref, error_alt, error_sigma);
    
    double error_ref_posterior = err_new.first;
    double error_alt_posterior = err_new.second; 
    
    fprintf(stderr, "Posterior error rates:\n");
    fprintf(stderr, "\tref mismatch: %f\n", error_ref_posterior);
    fprintf(stderr, "\talt mismatch: %f\n", error_alt_posterior);
    
    // Re-assign with posterior error rates
    print_elapsed(start_time, "Re-inferring identities (round 2)...");
    fprintf(stderr, "Re-inferring identities of cells...\n");
    assign_ids_parallel(cell_counts, samples, assn, assn_llr,
        allowed_ids, allowed_ids2, doublet_rate, error_ref_posterior, error_alt_posterior,
        false, prior_weights, n_threads, n_target);

    // Handle doublet-specific ID filtering
    bool do_again = true;
    if (idfile_doublet_given){
        // The user gave an allowed list of specific doublet combinations, and we included
        // all possible singlets from the allowable doublets in the first round. If any of 
        // these singlet identities turned out to be uncommon, we will assume the user was
        // right and that those singlet identities truly do not exist in the pool -- and
        // perform one more round of assignments without them.
        if (allowed_ids.size() > allowed_ids2.size()){
            bool altered = filter_identities(assn, assn_llr, samples.size(), allowed_ids, 
                allowed_ids2);

            if (altered){
                print_elapsed(start_time, "Re-inferring identities (round 3)...");
                fprintf(stderr, "Re-inferring with unlikely singlet identities removed...\n");
                assign_ids_parallel(cell_counts, samples, assn, assn_llr,
                    allowed_ids, allowed_ids2, doublet_rate, error_ref_posterior,
                    error_alt_posterior,
                    false, prior_weights, n_threads, n_target); 
                do_again = false;
            }
        }
    }

    // QC
    print_elapsed(start_time, "Running QC...");
    map<int, double> p_ncell;
    map<int, double> p_llr;
    id_qc(assn, assn_llr, p_ncell, p_llr);

    // Write assignments
    print_elapsed(start_time, "Writing outputs...");
    {
        string fname = output_prefix + ".assignments";
        FILE* outf = fopen(fname.c_str(), "w");
        fprintf(stderr, "Writing cell-individual assignments to disk...\n");
        dump_assignments(outf, assn, assn_llr, samples, barcode_group, cellranger, seurat, underscore);
        fclose(outf);
    }
    
    // Write summary
    {
        string fname = output_prefix + ".summary";
        FILE* outf = fopen(fname.c_str(), "w");
        write_summary(outf, output_prefix, assn, samples, error_ref,
            error_alt, error_sigma, error_ref_posterior,
            error_alt_posterior, vcf_file, vq, doublet_rate,
            p_ncell, p_llr);
        fclose(outf);
    }
    
    print_elapsed(start_time, "Complete!");
    fprintf(stderr, "Done!\n");
    
    return 0;
}
