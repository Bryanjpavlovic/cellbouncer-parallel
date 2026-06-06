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
const string VERSION = "2.01";
const string VERSION_MESSAGE = "ATAC dual-modality, identity prior, species panel";
const string VERSION_NEW = "v2.01: diagnostics `llr` column now reports winner-vs-runner-up (distinct from min_margin)";

// Global verbose flag (defined in demux_parallel_llr.cpp)
extern bool g_verbose;

// Global debug flag (defined in demux_parallel_llr.cpp)
extern bool g_debug;

// Species panel mode (used in main and potentially in helper functions)
enum class SpeciesPanelMode { NONE, COUNT_ONLY, FILTER, AUGMENT, BOTH };

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
        if (g_debug){
            fprintf(stderr, "DEBUG cell: assn=%d n_samples=%d is_combo=%d", 
                    a.second, n_samples, is_combo);
            if (is_combo) {
                fprintf(stderr, " combo=(%d,%d)", combo.first, combo.second);
            }
            fprintf(stderr, "\n");
        }
        
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
    if (g_debug){
        fprintf(stderr, "DEBUG: doublet_points=%d singlet_points=%d\n", doublet_points, singlet_points);
    }

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
    if (g_debug){
        fprintf(stderr, "DEBUG infer_error_rates_optimized:\n");
        fprintf(stderr, "  data_points=%lu total_n=%.2f total_k=%.2f total_weight=%.2f\n", 
                n.size(), total_n, total_k, total_weight);
        fprintf(stderr, "  overall_alt_frac=%.6f\n", total_k / total_n);
        for (auto& e : by_expected) {
            fprintf(stderr, "  expected=%.2f: n=%.2f k=%.2f alt_frac=%.6f\n",
                    e.first, e.second.first, e.second.second, 
                    e.second.second / e.second.first);
        }
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


static bool validate_gzip_readable(const string& filename){
    gzFile test_gz = gzopen(filename.c_str(), "r");
    if (!test_gz){
        return false;
    }
    char buf[1024 * 64];
    bool gz_ok = true;
    while (true){
        int n_read = gzread(test_gz, buf, sizeof(buf));
        if (n_read < 0){
            gz_ok = false;
            break;
        }
        if (n_read == 0) break;
    }
    int close_code = gzclose(test_gz);
    return gz_ok && close_code == Z_OK;
}

/**
 * Load cell counts from gzipped file into CellCounts structures.
 * Reverses dump_cellcounts_optimized: reads 7-column format, scales
 * float values back to int64 fixed-point representation.
 *
 * Format per line: barcode indv1 nalt1 indv2 nalt2 ref_count alt_count
 * When indv2==-1 && nalt2==-1: total row -> add_total(indv1, nalt1, ...)
 * Otherwise: pairwise row -> add(indv1, nalt1, indv2, nalt2, ...)
 */
int load_cellcounts_optimized(const string& filename,
    robin_hood::unordered_map<unsigned long, CellCounts>& cell_counts,
    int n_samples){
    
    gzreader reader(filename);
    int n_lines = 0;
    
    while(reader.next()){
        istringstream splitter(reader.line);
        string field;
        int idx = 0;
        
        unsigned long bc = 0;
        int indv1 = 0, nalt1 = 0, indv2 = 0, nalt2 = 0;
        float ref_val = 0.0f, alt_val = 0.0f;
        
        while(getline(splitter, field, '\t')){
            switch(idx){
                case 0: bc = strtoul(field.c_str(), NULL, 10); break;
                case 1: indv1 = atoi(field.c_str()); break;
                case 2: nalt1 = atoi(field.c_str()); break;
                case 3: indv2 = atoi(field.c_str()); break;
                case 4: nalt2 = atoi(field.c_str()); break;
                case 5: ref_val = atof(field.c_str()); break;
                case 6: alt_val = atof(field.c_str()); break;
            }
            idx++;
        }
        
        if (idx < 7) continue;
        
        // Initialize CellCounts for this barcode if needed
        if (cell_counts.count(bc) == 0){
            cell_counts.emplace(bc, CellCounts(n_samples));
        }
        
        // Convert float values back to fixed-point int64
        int64_t ref_scaled = (int64_t)round((double)ref_val * FIXED_POINT_SCALE);
        int64_t alt_scaled = (int64_t)round((double)alt_val * FIXED_POINT_SCALE);
        
        if (indv2 == -1 && nalt2 == -1){
            // Total row
            cell_counts[bc].add_total(indv1, nalt1, ref_scaled, alt_scaled);
        }
        else{
            // Pairwise row
            cell_counts[bc].add(indv1, nalt1, indv2, nalt2, ref_scaled, alt_scaled);
        }
        
        n_lines++;
    }
    
    return n_lines;
}


// ============================================================================
// Native species artifact helpers (Path Separation Contract V1_R3)
// ============================================================================

static int64_t scaled_from_float(double x){
    return (int64_t)llround(x * (double)FIXED_POINT_SCALE);
}

static vector<vector<pair<int, double>>> build_indiv_to_species_weights(
    const PanelMetadata& pm,
    int n_indiv){

    vector<vector<pair<int, double>>> out(n_indiv);
    for (int sp_idx = 0; sp_idx < (int)pm.species_list.size(); ++sp_idx){
        const string& sp = pm.species_list[sp_idx];
        auto it = pm.species_to_sample_indices.find(sp);
        if (it == pm.species_to_sample_indices.end()) continue;
        for (int indiv_idx : it->second){
            if (indiv_idx < 0 || indiv_idx >= n_indiv) continue;
            double w = pm.get_weight(sp, indiv_idx);
            if (w <= 0.0) continue;
            out[indiv_idx].push_back(make_pair(sp_idx, w));
        }
    }
    return out;
}

static void aggregate_cellcounts_to_species_native(
    const robin_hood::unordered_map<unsigned long, CellCounts>& indiv_counts,
    robin_hood::unordered_map<unsigned long, CellCounts>& species_counts,
    const PanelMetadata& pm,
    int n_indiv){

    const int n_species = (int)pm.species_list.size();
    if (n_species <= 0){
        fprintf(stderr, "ERROR: panel metadata contains zero species; cannot write native .species_counts\n");
        exit(1);
    }

    vector<vector<pair<int, double>>> i2sp = build_indiv_to_species_weights(pm, n_indiv);
    species_counts.clear();
    species_counts.reserve(indiv_counts.size());

    size_t unmapped_count_rows = 0;
    for (const auto& kv : indiv_counts){
        unsigned long bc = kv.first;
        const CellCounts& in = kv.second;
        CellCounts out(n_species);

        // Aggregate singlet/total rows: individual -> species at same genotype state.
        for (int i = 0; i < n_indiv; ++i){
            if (i2sp[i].empty()) continue;
            for (int g = 0; g < GENOTYPE_STATES; ++g){
                auto total = in.get_total(i, g);
                if (total.first == 0.0f && total.second == 0.0f) continue;
                for (const auto& m : i2sp[i]){
                    out.add_total(m.first, g,
                        scaled_from_float((double)total.first * m.second),
                        scaled_from_float((double)total.second * m.second));
                }
            }
        }

        // Aggregate heterotypic rows: individual pairs -> species pairs.
        for (int i = 0; i < n_indiv; ++i){
            if (i2sp[i].empty()) continue;
            for (int j = i + 1; j < n_indiv; ++j){
                if (i2sp[j].empty()) continue;
                for (int gi = 0; gi < GENOTYPE_STATES; ++gi){
                    for (int gj = 0; gj < GENOTYPE_STATES; ++gj){
                        auto pair_counts = in.get(i, gi, j, gj);
                        if (pair_counts.first == 0.0f && pair_counts.second == 0.0f) continue;
                        for (const auto& mi : i2sp[i]){
                            for (const auto& mj : i2sp[j]){
                                int sp_i = mi.first;
                                int sp_j = mj.first;
                                if (sp_i == sp_j) continue; // homospecies doublets are singlet-like at this resolution
                                double w = mi.second * mj.second;
                                if (sp_i < sp_j){
                                    out.add(sp_i, gi, sp_j, gj,
                                        scaled_from_float((double)pair_counts.first * w),
                                        scaled_from_float((double)pair_counts.second * w));
                                } else {
                                    out.add(sp_j, gj, sp_i, gi,
                                        scaled_from_float((double)pair_counts.first * w),
                                        scaled_from_float((double)pair_counts.second * w));
                                }
                            }
                        }
                    }
                }
            }
        }
        species_counts.emplace(bc, std::move(out));
    }

    for (int i = 0; i < n_indiv; ++i){
        if (i2sp[i].empty()) unmapped_count_rows++;
    }
    if (unmapped_count_rows > 0){
        fprintf(stderr, "WARNING: %lu VCF samples had no species mapping and were excluded from native species aggregation\n",
            unmapped_count_rows);
    }
    fprintf(stderr, "Native species aggregation: %lu cells, %d individuals -> %d species\n",
        indiv_counts.size(), n_indiv, n_species);
}

static void compute_species_condf_native(
    robin_hood::unordered_map<int, ChromSNPs>& species_snpdat,
    map<pair<int, int>, map<int, float>>& species_condf,
    const PanelMetadata& pm,
    int n_indiv){

    const int n_species = (int)pm.species_list.size();
    vector<vector<pair<int, double>>> i2sp = build_indiv_to_species_weights(pm, n_indiv);

    for (auto& kv : species_snpdat){
        for (auto& snp : kv.second.snps){
            if ((int)snp.geno.size() != n_indiv){
                snp.precompute_genotypes(n_indiv);
            }
        }
    }

    map<pair<int, int>, map<int, double>> sums;
    map<pair<int, int>, map<int, double>> tots;

    for (auto& chrom_kv : species_snpdat){
        for (auto& snp : chrom_kv.second.snps){
            const int8_t* geno = snp.geno.data();

            vector<double> target_avg(n_species, 0.0);
            vector<double> target_wsum(n_species, 0.0);
            for (int j = 0; j < n_indiv; ++j){
                int8_t gj = geno[j];
                if (gj < 0 || gj >= GENOTYPE_STATES) continue;
                for (const auto& mj : i2sp[j]){
                    target_avg[mj.first] += mj.second * ((double)gj / 2.0);
                    target_wsum[mj.first] += mj.second;
                }
            }
            for (int sp = 0; sp < n_species; ++sp){
                if (target_wsum[sp] > 0.0) target_avg[sp] /= target_wsum[sp];
            }

            for (int i = 0; i < n_indiv; ++i){
                int8_t gi = geno[i];
                if (gi < 0 || gi >= GENOTYPE_STATES) continue;
                if (i2sp[i].empty()) continue;
                for (const auto& mi : i2sp[i]){
                    pair<int, int> row = make_pair(mi.first, (int)gi);
                    for (int sp_t = 0; sp_t < n_species; ++sp_t){
                        if (target_wsum[sp_t] <= 0.0) continue;
                        sums[row][sp_t] += mi.second * target_avg[sp_t];
                        tots[row][sp_t] += mi.second;
                    }
                }
            }
        }
    }

    species_condf.clear();
    for (auto& row_kv : sums){
        for (auto& col_kv : row_kv.second){
            double denom = tots[row_kv.first][col_kv.first];
            if (denom > 0.0){
                species_condf[row_kv.first][col_kv.first] = (float)(col_kv.second / denom);
            }
        }
    }

    fprintf(stderr, "Native species condf: %lu row entries, %d species\n",
        species_condf.size(), n_species);
}

static void load_panel_metadata_if_needed(
    const string& panel_metadata_file,
    const vector<string>& samples,
    PanelMetadata& panel_meta,
    bool& panel_meta_loaded){

    if (!panel_meta_loaded){
        if (panel_metadata_file.empty()){
            fprintf(stderr, "ERROR: panel metadata required for native species artifacts\n");
            exit(1);
        }
        panel_meta = load_panel_metadata(panel_metadata_file, samples);
        panel_meta_loaded = true;
    }
}

static void write_native_species_assignments(
    const string& output_prefix,
    robin_hood::unordered_map<unsigned long, CellCounts>& species_counts_native,
    const vector<string>& species_samples,
    double doublet_rate,
    double error_ref,
    double error_alt,
    int n_threads,
    int n_target,
    const string& barcode_group,
    bool cellranger,
    bool seurat,
    bool underscore){

    if (species_counts_native.empty()){
        fprintf(stderr, "WARNING: native species assignment requested but no species counts are available\n");
        return;
    }

    robin_hood::unordered_map<unsigned long, int> sp_assn;
    robin_hood::unordered_map<unsigned long, double> sp_assn_llr;
    set<int> allowed_ids;
    set<int> allowed_ids2;
    map<int, double> prior_weights;

    assign_ids_parallel(species_counts_native, const_cast<vector<string>&>(species_samples),
        sp_assn, sp_assn_llr, allowed_ids, allowed_ids2, doublet_rate,
        error_ref, error_alt, false, prior_weights, n_threads, n_target);

    string fname = output_prefix + ".species_assignments";
    FILE* outf = fopen(fname.c_str(), "w");
    if (!outf){
        fprintf(stderr, "ERROR: could not open %s for writing\n", fname.c_str());
        exit(1);
    }
    vector<string> tmp_species_samples = species_samples;
    dump_assignments(outf, sp_assn, sp_assn_llr, tmp_species_samples,
        const_cast<string&>(barcode_group), cellranger, seurat, underscore);
    fclose(outf);
    fprintf(stderr, "Wrote native species assignments for %lu cells to %s\n",
        sp_assn.size(), fname.c_str());
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
    fprintf(stderr, "    --shared_het_vcf NAME Shared memory het VCF for ploidy stats (from vcf_loader_daemon)\n");
    fprintf(stderr, "    --vcf_chroms -c File listing chromosomes to use from VCF\n");
    fprintf(stderr, "    --libname -n Library name to append to barcodes\n");
    fprintf(stderr, "    --cellranger -C Format barcodes for CellRanger\n");
    fprintf(stderr, "    --seurat -R Format barcodes for Seurat\n");
    fprintf(stderr, "    --underscore -U Use underscore instead of hyphen for libname\n");
    fprintf(stderr, "    --disable_conditional -f Disable computing conditional match fractions\n");
    fprintf(stderr, "    --dump_conditional -F Load VCF, compute conditional match fractions, write\n");
    fprintf(stderr, "       .condf file, and exit. No BAM required. Use this to generate .condf after\n");
    fprintf(stderr, "       a run that used -f.\n");
    fprintf(stderr, "    --n_target -N Max singlets to keep before doublet eval [-1=auto, 0=no limit]\n");
    fprintf(stderr, "    --verbose -V Enable verbose output\n");
    fprintf(stderr, "\n===== DIAGNOSTIC OUTPUT (NEW) =====\n");
    fprintf(stderr, "    --debug          Enable DEBUG spam to stderr\n");
    fprintf(stderr, "    --diagnostics    Write diagnostic files (default: ON)\n");
    fprintf(stderr, "    --no-diagnostics Disable diagnostic file output\n");
    fprintf(stderr, "    --het_vcf FILE   Het VCF for ploidy stats (from downsample_vcf_parallel)\n");
    fprintf(stderr, "    --het_method M   Het balance method: welford (default) or persite\n");
    fprintf(stderr, "                     welford: Online variance, low memory\n");
    fprintf(stderr, "                     persite: Store per-site counts, more memory\n");
    fprintf(stderr, "    --min_het_sites N  Minimum het sites for variance (default: 100)\n");
    fprintf(stderr, "    --min_het_depth D  Minimum depth per het site for persite method (default: 5.0)\n");
    fprintf(stderr, "    --n_runner_ups N Number of runner-ups to report [8]\n");
    fprintf(stderr, "    --close_threshold F LLR threshold for n_close [20.0]\n");
    fprintf(stderr, "\n===== SKIP ASSIGNMENT =====\n");
    fprintf(stderr, "    --skip_assignment -K Write .counts and exit (no assignment)\n");
    fprintf(stderr, "\n===== PILEUP (variant-consistency benchmark) =====\n");
    fprintf(stderr, "    --dump_pileup PREFIX Emit PREFIX.pileup_sites.tsv.gz and PREFIX.pileup_obs.tsv.gz\n");
    fprintf(stderr, "    --dump_selection_audit  Append maximin-vs-maxllr winner columns to .diagnostics.gz\n");
    fprintf(stderr, "                            (maximin_winner, maximin_score, maxllr_winner, maxllr_score, selection_agree).\n");
    fprintf(stderr, "                            For comparing the current selection criterion against the original maxllr elimination.\n");
    fprintf(stderr, "                         for the interindividual panel. Use with -B (candidate\n");
    fprintf(stderr, "                         barcodes) and -K. Routes through the single-panel path,\n");
    fprintf(stderr, "                         so do not pass a species VCF on the pileup run.\n");
    fprintf(stderr, "\n===== COUNTS FILE SAFETY =====\n");
    fprintf(stderr, "    --reuse_counts       Load existing .counts file after validating integrity\n");
    fprintf(stderr, "    --force_recount      Delete existing .counts and recount from BAM\n");
    fprintf(stderr, "    (default: error if .counts exists, to prevent stale/truncated reuse)\n");
    fprintf(stderr, "\n===== ATAC DUAL-MODALITY (2A) =====\n");
    fprintf(stderr, "    --atac_bam FILE      ATAC BAM (expects RNA-aligned barcodes in CB tag)\n");
    fprintf(stderr, "    --atac_vcf FILE      ATAC-demux VCF panel (required if --atac_bam set)\n");
    fprintf(stderr, "    --atac_het_vcf FILE  ATAC-side het VCF for ATAC het balance diagnostics\n");
    fprintf(stderr, "    --atac_shared_vcf NAME   Shared-memory ATAC VCF (mutually exclusive with --atac_vcf)\n");
    fprintf(stderr, "    --atac_shared_het_vcf NAME   Shared-memory ATAC het VCF\n");
    fprintf(stderr, "\n===== IDENTITY PRIOR (2B) =====\n");
    fprintf(stderr, "    --identity_prior FILE  .contam_prof_empty file for Bayesian prior on identity\n");
    fprintf(stderr, "\n===== SPECIES PANEL (2C) =====\n");
    fprintf(stderr, "    --species_vcf FILE         Species-discrimination VCF panel\n");
    fprintf(stderr, "    --species_shared_vcf NAME  Shared-memory species VCF\n");
    fprintf(stderr, "    --species_panel_mode MODE  count_only (native V1_R3 mode; required if --species_vcf set)\n");
    fprintf(stderr, "       count_only: count/write native species panel artifacts without using species SNPs in individual assignment\n");
    fprintf(stderr, "       filter|augment|both are legacy mixed modes and require --allow_legacy_mixed_species_panel\n");
    fprintf(stderr, "    --allow_legacy_mixed_species_panel  Permit old species-panel filter/augment/both behavior (not V1_R3 native)\n");
    fprintf(stderr, "    --species_assignment_output Produce .species_assignments (species-only scoring pass)\n");
    fprintf(stderr, "    --species_counts_output     Write .species_counts, .species_condf alongside .counts\n");
    fprintf(stderr, "    --panel_metadata FILE  Tab-separated indiv_id/species mapping (required for species output)\n");
    fprintf(stderr, "\n    --help -h Display this message and exit\n");
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

// ============================================================================
// NEW: DIAGNOSTIC OUTPUT FUNCTIONS
// ============================================================================

/**
 * Write per-cell diagnostics to gzipped TSV file
 */
void write_diagnostics_gz(
    const string& filename,
    vector<string>& samples,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, double>& assignments_llr,
    robin_hood::unordered_map<unsigned long, CellDiagnostics>& diagnostics,
    int n_samples,
    const string& libname,
    bool cellranger,
    bool seurat,
    bool underscore,
    // Optional ATAC het data for extra columns
    robin_hood::unordered_map<unsigned long, CellHetData>* atac_het_data = NULL,
    robin_hood::unordered_map<int, ChromSNPs>* atac_het_snpdat_ptr = NULL,
    vector<pair<int, int>>* atac_idx_to_site_ptr = NULL,
    HetBalanceMethod atac_het_method = HetBalanceMethod::WELFORD,
    int atac_min_het_sites = 100,
    double atac_min_het_depth = 5.0,
    // When true, append maximin-vs-maxllr selection-audit columns at the end of
    // each row. Off by default so the standard .diagnostics.gz layout (read
    // positionally by tetra_refine) is unchanged.
    bool dump_selection_audit = false){
    
    gzFile outf = gzopen(filename.c_str(), "w");
    if (!outf){
        fprintf(stderr, "ERROR: Could not open %s for writing\n", filename.c_str());
        return;
    }
    
    bool write_atac_het = (atac_het_data != NULL && !atac_het_data->empty());

    // Header
    gzprintf(outf, "barcode\tassignment\tsinglet_doublet\tllr\tmin_margin\tworst_competitor\t"
                   "n_close\ttotal_depth\thet_balance_var\tn_het_sites\thet_total_depth\t"
                   "posterior\tentropy");
    if (write_atac_het){
        gzprintf(outf, "\tatac_het_balance_var\tatac_n_het_sites\tatac_het_total_depth");
    }
    if (dump_selection_audit){
        gzprintf(outf, "\tmaximin_winner\tmaximin_score\tmaxllr_winner\tmaxllr_score\tselection_agree");
    }
    gzprintf(outf, "\n");
    
    for (const auto& kv : assignments){
        unsigned long bc = kv.first;
        int assn = kv.second;
        
        string bc_str = bc2str(bc);
        mod_bc_libname(bc_str, libname, cellranger, seurat, underscore);
        
        string identity = idx2name(assn, samples);
        char s_d = (assn < n_samples) ? 'S' : 'D';

        // Get diagnostics
        CellDiagnostics diag;
        if (diagnostics.count(bc) > 0){
            diag = diagnostics.at(bc);
        }

        // The `llr` column reports the winner's margin over the rank-1 runner-up
        // (best vs second-best), distinct from min_margin (worst pairwise margin
        // over all competitors). assignments_llr holds the maximin selection
        // score (== min_margin) and is intentionally NOT used here so the two
        // columns differ. assignments_llr still drives the separate .assignments
        // file and winner selection, which are unchanged.
        double llr = diag.llr_vs_runnerup;

        // Format worst competitor
        string worst_comp = ".";
        if (diag.worst_competitor >= 0){
            worst_comp = idx2name(diag.worst_competitor, samples);
        }
        
        gzprintf(outf, "%s\t%s\t%c\t%.4f\t%.4f\t%s\t%d\t%.2f\t%.6f\t%d\t%.2f\t%.6f\t%.4f",
            bc_str.c_str(),
            identity.c_str(),
            s_d,
            llr,
            diag.min_margin,
            worst_comp.c_str(),
            diag.n_close,
            diag.total_depth,
            diag.het_balance_var,
            diag.n_het_sites,
            diag.het_total_depth,
            diag.posterior,
            diag.entropy);

        if (write_atac_het){
            // Compute ATAC het balance for this cell
            CellDiagnostics atac_diag;
            auto atac_it = atac_het_data->find(bc);
            if (atac_it != atac_het_data->end()){
                const CellHetData& cell_atac_het = atac_it->second;
                if (atac_het_method == HetBalanceMethod::PERSITE &&
                    atac_het_snpdat_ptr != NULL && atac_idx_to_site_ptr != NULL){
                    compute_het_balance_persite(cell_atac_het.persite_data,
                        *atac_het_snpdat_ptr, *atac_idx_to_site_ptr,
                        assn, n_samples, atac_min_het_depth, atac_min_het_sites,
                        atac_diag);
                }
                else{
                    compute_het_balance_welford(cell_atac_het.welford_stats,
                        assn, n_samples, atac_min_het_sites, atac_diag);
                }
            }
            gzprintf(outf, "\t%.6f\t%d\t%.2f",
                atac_diag.het_balance_var, atac_diag.n_het_sites, atac_diag.het_total_depth);
        }
        if (dump_selection_audit){
            string maxllr_name = (diag.maxllr_winner >= 0)
                ? idx2name(diag.maxllr_winner, samples) : ".";
            int agree = (diag.maxllr_winner == assn) ? 1 : 0;
            gzprintf(outf, "\t%s\t%.4f\t%s\t%.4f\t%d",
                identity.c_str(),       // maximin_winner == the actual assignment
                diag.min_margin,        // maximin selection score
                maxllr_name.c_str(),    // maxllr criterion winner
                diag.maxllr_winner_score,
                agree);
        }
        gzprintf(outf, "\n");
    }
    
    gzclose(outf);
    fprintf(stderr, "Wrote diagnostics for %lu cells to %s\n", assignments.size(), filename.c_str());
}

/**
 * Write runner-ups to gzipped TSV file
 */
void write_runner_ups_gz(
    const string& filename,
    vector<string>& samples,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, vector<RunnerUp> >& runner_ups,
    int n_samples,
    const string& libname,
    bool cellranger,
    bool seurat,
    bool underscore){
    
    gzFile outf = gzopen(filename.c_str(), "w");
    if (!outf){
        fprintf(stderr, "ERROR: Could not open %s for writing\n", filename.c_str());
        return;
    }
    
    // Header
    gzprintf(outf, "barcode\trank\tidentity\tllr_vs_winner\tmin_margin\n");
    
    for (const auto& kv : runner_ups){
        unsigned long bc = kv.first;
        const vector<RunnerUp>& runners = kv.second;
        
        string bc_str = bc2str(bc);
        mod_bc_libname(bc_str, libname, cellranger, seurat, underscore);
        
        for (size_t i = 0; i < runners.size(); ++i){
            const RunnerUp& runner = runners[i];
            string identity = idx2name(runner.identity, samples);
            
            gzprintf(outf, "%s\t%d\t%s\t%.4f\t%.4f\n",
                bc_str.c_str(),
                (int)(i + 1),  // 1-indexed rank
                identity.c_str(),
                runner.llr_vs_winner,
                runner.min_margin);
        }
    }
    
    gzclose(outf);
    fprintf(stderr, "Wrote runner-ups for %lu cells to %s\n", runner_ups.size(), filename.c_str());
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
       {"dump_conditional", no_argument, 0, 'F'},
       {"no_preload", no_argument, 0, 'P'},
       {"vcf_chroms", required_argument, 0, 'c'},
       {"threads", required_argument, 0, 't'},
       {"shared_vcf", required_argument, 0, 'S'},
       {"shared_het_vcf", required_argument, 0, 1007},
       {"n_target", required_argument, 0, 'N'},
       {"verbose", no_argument, 0, 'V'},
       // NEW diagnostic options
       {"debug", no_argument, 0, 1001},
       {"diagnostics", no_argument, 0, 1002},
       {"no-diagnostics", no_argument, 0, 1003},
       {"het_vcf", required_argument, 0, 1004},
       {"n_runner_ups", required_argument, 0, 1005},
       {"close_threshold", required_argument, 0, 1006},
       {"het_method", required_argument, 0, 1010},
       {"min_het_sites", required_argument, 0, 1011},
       {"min_het_depth", required_argument, 0, 1012},
       // Step 0a: skip assignment
       {"skip_assignment", no_argument, 0, 'K'},
       // Variant-consistency benchmark: per-SNP pileup sidecars
       {"dump_pileup", required_argument, 0, 1099},
       {"dump_selection_audit", no_argument, 0, 1100},
       // 2A: ATAC dual-modality
       {"atac_bam", required_argument, 0, 1020},
       {"atac_vcf", required_argument, 0, 1021},
       {"atac_het_vcf", required_argument, 0, 1022},
       {"atac_shared_vcf", required_argument, 0, 1023},
       {"atac_shared_het_vcf", required_argument, 0, 1024},
       // 2B: identity prior
       {"identity_prior", required_argument, 0, 1030},
       // 2C: species panel
       {"species_vcf", required_argument, 0, 1040},
       {"species_shared_vcf", required_argument, 0, 1041},
       {"species_panel_mode", required_argument, 0, 1042},
       {"allow_legacy_mixed_species_panel", no_argument, 0, 1059},
       {"species_assignment_output", no_argument, 0, 1043},
       {"species_counts_output", no_argument, 0, 1045},
       {"panel_metadata", required_argument, 0, 1044},
       {"reuse_counts", no_argument, 0, 1050},
       {"force_recount", no_argument, 0, 1051},
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
    bool dump_conditional = false;

    bool no_preload = false;
    string vcf_chroms_file = "";
    bool vcf_chroms_given = false;
    
    // New parallel processing options
    int n_threads = omp_get_num_procs();
    int htslib_threads = 0;  // 0 = auto-calculate after parsing args
    string shared_vcf_name = "";
    string shared_het_vcf_name = "";
    int n_target = -1;    // -1=auto, 0=no prune, >0=use value
    bool verbose = false;
    
    // NEW diagnostic options
    bool debug_mode = false;
    bool write_diagnostics = true;  // Default ON
    string het_vcf_file = "";
    HetBalanceMethod het_method = HetBalanceMethod::WELFORD;  // Default to Welford
    int min_het_sites = 100;   // Minimum sites for reliable variance
    double min_het_depth = 5.0;  // Minimum depth per site (persite method only)
    int n_runner_ups = 8;
    double close_threshold = 20.0;

    // Step 0a: skip assignment
    bool skip_assignment = false;
    bool dump_selection_audit = false;

    // Variant-consistency benchmark: emit per-SNP pileup sidecars
    bool dump_pileup = false;
    string dump_pileup_prefix = "";

    // 2A: ATAC dual-modality
    string atac_bamfile = "";
    string atac_vcf_file = "";
    string atac_het_vcf_file = "";
    string atac_shared_vcf_name = "";
    string atac_shared_het_vcf_name = "";

    // 2B: identity prior
    string identity_prior_file = "";

    // 2C: species panel
    string species_vcf_file = "";
    string species_shared_vcf_name = "";
    string species_panel_mode_str = "";
    bool allow_legacy_mixed_species_panel = false;
    bool species_assignment_output = false;
    bool species_counts_output = false;
    string panel_metadata_file = "";

    // Counts file reuse safety
    bool reuse_counts = false;
    bool force_recount = false;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:v:o:B:i:I:q:D:n:e:E:s:c:t:S:N:fFPCRUVhK", 
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
            case 'F':
                dump_conditional = true;
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
            // NEW diagnostic options
            case 1001:  // --debug
                debug_mode = true;
                break;
            case 1002:  // --diagnostics
                write_diagnostics = true;
                break;
            case 1003:  // --no-diagnostics
                write_diagnostics = false;
                break;
            case 1004:  // --het_vcf
                het_vcf_file = optarg;
                break;
            case 1005:  // --n_runner_ups
                n_runner_ups = atoi(optarg);
                break;
            case 1006:  // --close_threshold
                close_threshold = atof(optarg);
                break;
            case 1007:  // --shared_het_vcf
                shared_het_vcf_name = optarg;
                break;
            case 1010:  // --het_method
                if (strcmp(optarg, "persite") == 0) {
                    het_method = HetBalanceMethod::PERSITE;
                } else if (strcmp(optarg, "welford") == 0) {
                    het_method = HetBalanceMethod::WELFORD;
                } else {
                    fprintf(stderr, "ERROR: Unknown het_method: %s\n", optarg);
                    fprintf(stderr, "Valid options: welford, persite\n");
                    exit(1);
                }
                break;
            case 1011:  // --min_het_sites
                min_het_sites = atoi(optarg);
                break;
            case 1012:  // --min_het_depth
                min_het_depth = atof(optarg);
                break;
            case 'K':  // --skip_assignment
                skip_assignment = true;
                break;
            case 1099:  // --dump_pileup
                dump_pileup = true;
                dump_pileup_prefix = optarg;
                break;
            case 1100:  // --dump_selection_audit
                dump_selection_audit = true;
                break;
            case 1020:  // --atac_bam
                atac_bamfile = optarg;
                break;
            case 1021:  // --atac_vcf
                atac_vcf_file = optarg;
                break;
            case 1022:  // --atac_het_vcf
                atac_het_vcf_file = optarg;
                break;
            case 1023:  // --atac_shared_vcf
                atac_shared_vcf_name = optarg;
                break;
            case 1024:  // --atac_shared_het_vcf
                atac_shared_het_vcf_name = optarg;
                break;
            case 1030:  // --identity_prior
                identity_prior_file = optarg;
                break;
            case 1040:  // --species_vcf
                species_vcf_file = optarg;
                break;
            case 1041:  // --species_shared_vcf
                species_shared_vcf_name = optarg;
                break;
            case 1042:  // --species_panel_mode
                species_panel_mode_str = optarg;
                break;
            case 1059:  // --allow_legacy_mixed_species_panel
                allow_legacy_mixed_species_panel = true;
                break;
            case 1043:  // --species_assignment_output
                species_assignment_output = true;
                break;
            case 1045:  // --species_counts_output
                species_counts_output = true;
                break;
            case 1044:  // --panel_metadata
                panel_metadata_file = optarg;
                break;
            case 1050:  // --reuse_counts
                reuse_counts = true;
                break;
            case 1051:  // --force_recount
                force_recount = true;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Set global verbose flag
    g_verbose = verbose;
    
    // Set global debug flag
    g_debug = debug_mode;
    
    // Error check arguments
    if (reuse_counts && force_recount){
        fprintf(stderr, "ERROR: --reuse_counts and --force_recount are mutually exclusive\n");
        exit(1);
    }
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

    // Validate new flags
    bool atac_mode = (atac_bamfile.length() > 0);
    if (atac_mode){
        if (atac_vcf_file.length() == 0 && atac_shared_vcf_name.length() == 0){
            fprintf(stderr, "ERROR: --atac_vcf or --atac_shared_vcf required when --atac_bam is set\n");
            exit(1);
        }
        if (atac_vcf_file.length() > 0 && atac_shared_vcf_name.length() > 0){
            fprintf(stderr, "ERROR: --atac_vcf and --atac_shared_vcf are mutually exclusive\n");
            exit(1);
        }
    }
    if (atac_het_vcf_file.length() > 0 && !atac_mode){
        fprintf(stderr, "WARNING: --atac_het_vcf ignored without --atac_bam\n");
    }

    // Species panel mode validation
    SpeciesPanelMode species_mode = SpeciesPanelMode::NONE;
    bool has_species_vcf = (species_vcf_file.length() > 0 || species_shared_vcf_name.length() > 0);
    if (has_species_vcf){
        if (species_panel_mode_str.length() == 0){
            fprintf(stderr, "ERROR: --species_panel_mode required when --species_vcf is set\n");
            exit(1);
        }
        if (species_panel_mode_str == "count_only" || species_panel_mode_str == "native" || species_panel_mode_str == "none") species_mode = SpeciesPanelMode::COUNT_ONLY;
        else if (species_panel_mode_str == "filter" || species_panel_mode_str == "augment" || species_panel_mode_str == "both"){
            if (!allow_legacy_mixed_species_panel){
                fprintf(stderr, "ERROR: --species_panel_mode %s is a legacy mixed species-panel mode.\n",
                    species_panel_mode_str.c_str());
                fprintf(stderr, "       Under V1_R3 path separation, use --species_panel_mode count_only.\n");
                fprintf(stderr, "       Use --allow_legacy_mixed_species_panel only for intentional legacy debugging.\n");
                exit(1);
            }
            if (species_panel_mode_str == "filter") species_mode = SpeciesPanelMode::FILTER;
            else if (species_panel_mode_str == "augment") species_mode = SpeciesPanelMode::AUGMENT;
            else species_mode = SpeciesPanelMode::BOTH;
        }
        else{
            fprintf(stderr, "ERROR: --species_panel_mode must be count_only (or legacy filter/augment/both with --allow_legacy_mixed_species_panel)\n");
            exit(1);
        }
        if (species_vcf_file.length() > 0 && species_shared_vcf_name.length() > 0){
            fprintf(stderr, "ERROR: --species_vcf and --species_shared_vcf are mutually exclusive\n");
            exit(1);
        }
    }
    if (species_assignment_output && !has_species_vcf){
        fprintf(stderr, "ERROR: --species_assignment_output requires --species_vcf\n");
        exit(1);
    }
    if ((species_counts_output || species_assignment_output) && panel_metadata_file.length() == 0){
        fprintf(stderr, "ERROR: --species_counts_output/--species_assignment_output requires --panel_metadata so .species_* artifacts are species-native\n");
        exit(1);
    }
    
    // BAM header/TID mapping is read directly through HTSlib below when needed.
    
    // Check if we can load counts from previous run
    bool load_counts = false;
    string countsfilename = output_prefix + ".counts";
    if (file_exists(countsfilename)){
        if (force_recount){
            fprintf(stderr, "Found existing counts file: %s\n", countsfilename.c_str());
            fprintf(stderr, "  --force_recount set: deleting stale counts and recounting from BAM\n");
            remove(countsfilename.c_str());
            // Also remove companion files that may be stale
            string stale_samples = output_prefix + ".samples";
            string stale_condf = output_prefix + ".condf";
            string stale_sp_counts = output_prefix + ".species_counts";
            string stale_sp_condf = output_prefix + ".species_condf";
            if (file_exists(stale_samples)) remove(stale_samples.c_str());
            if (file_exists(stale_condf)) remove(stale_condf.c_str());
            if (file_exists(stale_sp_counts)) remove(stale_sp_counts.c_str());
            if (file_exists(stale_sp_condf)) remove(stale_sp_condf.c_str());
        } else if (reuse_counts){
            // Validate the existing file before trusting it
            fprintf(stderr, "Found existing counts file: %s\n", countsfilename.c_str());
            fprintf(stderr, "  --reuse_counts set: validating before loading\n");

            // Check .samples companion
            string samplesfile = output_prefix + ".samples";
            if (!file_exists(samplesfile)){
                fprintf(stderr, "ERROR: --reuse_counts: .samples file missing alongside .counts\n");
                fprintf(stderr, "  Delete %s and rerun, or use --force_recount\n", countsfilename.c_str());
                exit(1);
            }

            // Check .counts is gzip-readable (not truncated)
            if (!validate_gzip_readable(countsfilename)){
                fprintf(stderr, "ERROR: --reuse_counts: .counts file appears truncated or corrupt: %s\n",
                    countsfilename.c_str());
                fprintf(stderr, "  Delete it and rerun, or use --force_recount\n");
                exit(1);
            }

            // If species output is requested, the species companion outputs must
            // be present and readable too.  Otherwise a reused main .counts can
            // silently coexist with stale or missing .species_counts/.species_condf.
            if (species_counts_output){
                string sp_counts = output_prefix + ".species_counts";
                string sp_condf = output_prefix + ".species_condf";
                if (!file_exists(sp_counts)){
                    fprintf(stderr, "ERROR: --reuse_counts with --species_counts_output: missing %s\n",
                        sp_counts.c_str());
                    fprintf(stderr, "  Use --force_recount to regenerate the full dual-panel output set.\n");
                    exit(1);
                }
                if (!file_exists(sp_condf)){
                    fprintf(stderr, "ERROR: --reuse_counts with --species_counts_output: missing %s\n",
                        sp_condf.c_str());
                    fprintf(stderr, "  Use --force_recount to regenerate the full dual-panel output set.\n");
                    exit(1);
                }
                if (!validate_gzip_readable(sp_counts)){
                    fprintf(stderr, "ERROR: --reuse_counts: .species_counts appears truncated or corrupt: %s\n",
                        sp_counts.c_str());
                    fprintf(stderr, "  Use --force_recount to regenerate the full dual-panel output set.\n");
                    exit(1);
                }
            }

            load_counts = true;
            fprintf(stderr, "  Validation passed; loading existing counts\n");
        } else {
            // Default: refuse to silently reuse stale counts
            fprintf(stderr, "ERROR: existing .counts file found: %s\n", countsfilename.c_str());
            fprintf(stderr, "  A previous run may have left this file truncated or stale.\n");
            fprintf(stderr, "  Use --reuse_counts to load it after validation,\n");
            fprintf(stderr, "  or --force_recount to delete it and recount from BAM.\n");
            exit(1);
        }
    }
    else{
        if (bamfile.length() == 0 && !dump_conditional){
            fprintf(stderr, "ERROR: bam file (--bam) required\n");
            exit(1);
        }
        // BAM readability/header validation is performed by direct HTSlib helpers
        // when chromosome/TID mapping is needed.
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
    
    // Process samples now if not using shared VCF (otherwise defer until after attach).
    // EXCEPTION: when load_counts is true, sample names are already populated from
    // the .samples companion file above (or read_vcf_samples), so we can and must
    // parse the idfile here even when a shared VCF is in use. The deferred parse
    // further below only runs inside the (!load_counts) branch, so without this the
    // combination (--reuse_counts + --shared_vcf + -i/-I) would leave allowed_ids
    // empty and silently disable the identity restriction, letting the assignment
    // step consider every possible singlet/doublet combination.
    if (shared_vcf_name.length() == 0 || load_counts){
        fprintf(stderr, "Number of individuals in VCF: %lu\n", samples.size());
        
        if (idfile_given){
            parse_idfile(idfile, samples, allowed_ids, allowed_ids2, true);
            if (allowed_ids.size() == 0){
                fprintf(stderr, "ERROR: No valid individual names found in %s; refusing to ignore -i\n", idfile.c_str());
                exit(1);
            }
        }
        if (idfile_doublet_given){
            parse_idfile(idfile_doublet, samples, allowed_ids, allowed_ids2, false);
            if (allowed_ids.size() == 0){
                fprintf(stderr, "ERROR: No valid individual names found in %s; refusing to ignore -I\n", idfile_doublet.c_str());
                exit(1);
            }
        }

        if (idfile_given || idfile_doublet_given){
            fprintf(stderr, "Identity restriction active: allowed_ids=%lu allowed_ids2=%lu\n",
                allowed_ids.size(), allowed_ids2.size());
        }
        
        if (samples_from_vcf){
            string samplesfile = output_prefix + ".samples"; 
            write_samples(samplesfile, samples);
        }
    }
    
    // Handle --dump_conditional / -F mode: VCF-only .condf generation
    if (dump_conditional){
        if (vcf_file.length() == 0 && shared_vcf_name.length() == 0){
            fprintf(stderr, "ERROR: --vcf/-v or --shared_vcf/-S required for --dump_conditional/-F\n");
            exit(1);
        }
        
        robin_hood::unordered_map<int, ChromSNPs> snpdat_optimized;
        
        if (shared_vcf_name.length() > 0){
            // Fast path: attach to shared memory daemon
            print_elapsed(start_time, "dump_conditional mode: attaching to shared VCF...");
            fprintf(stderr, "Attaching to shared VCF: %s\n", shared_vcf_name.c_str());
            
            if (!attach_shared_vcf(shared_vcf_name, snpdat_optimized, samples)){
                fprintf(stderr, "ERROR: Could not attach to shared VCF: %s\n", shared_vcf_name.c_str());
                exit(1);
            }
            fprintf(stderr, "Attached. %lu chromosomes, %lu samples\n",
                snpdat_optimized.size(), samples.size());
        }
        else{
            // Disk path: load VCF from file
            print_elapsed(start_time, "dump_conditional mode: loading VCF from disk...");
            
            set<string> chroms_vcf;
            get_vcf_chroms(vcf_file, chroms_vcf);
            
            // Build synthetic seq2tid (no BAM header available)
            map<string, int> seq2tid_synthetic;
            int tid_counter = 0;
            for (const auto& c : chroms_vcf){
                seq2tid_synthetic[c] = tid_counter++;
            }
            
            int nloaded = read_vcf_chroms_optimized(vcf_file, chroms_vcf,
                seq2tid_synthetic, snpdat_optimized, vq);
            fprintf(stderr, "Loaded %d SNPs from %lu chromosomes\n", nloaded, chroms_vcf.size());
        }
        
        // Compute conditional match fractions. The parallel CONDF function now
        // performs genotype-only precompute internally, avoiding unnecessary
        // per-SNP pair_targets allocation in VCF-only mode.
        print_elapsed(start_time, "Computing conditional match fractions...");
        map<pair<int, int>, map<int, float>> conditional_match_fracs;
        
        compute_conditional_match_fracs_parallel(snpdat_optimized,
            conditional_match_fracs, samples.size(), n_threads);
        
        // Write .condf
        string outname = output_prefix + ".condf";
        FILE* outf = fopen(outname.c_str(), "w");
        dump_exp_fracs(outf, conditional_match_fracs);
        fclose(outf);
        
        print_elapsed(start_time, "Done. Wrote .condf file.");
        fprintf(stderr, "Wrote conditional match fractions to %s\n", outname.c_str());
        return 0;
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

    // Position set for species panel dedup (§5.5)
    // Populated from demux VCF snpdat before it goes out of scope;
    // used later to filter overlapping sites from the species VCF.
    // Key: (tid, position)
    set<pair<int, int>> demux_positions;
    
    // Species panel data - declared early so dual counting path can populate them
    robin_hood::unordered_map<int, ChromSNPs> species_snpdat;
    robin_hood::unordered_map<unsigned long, CellCounts> species_cell_counts_rna;
    robin_hood::unordered_map<unsigned long, CellCounts> species_cell_counts_atac;
    PanelMetadata panel_meta;
    bool panel_meta_loaded = false;
    bool species_counted_dual = false;  // Set by WP3 dual counting path

    // Main BAM seqname->TID mapping, populated once for BAM-backed runs.
    // Keep it outside the counting block because later disk-based species/het
    // VCF loaders need the same TID mapping.
    map<string, int> bam_seq2tid;
    
    if (load_counts){
        // Load counts from previous run
        print_elapsed(start_time, "Loading counts from existing file...");
        fprintf(stderr, "Loading allele counts from %s\n", countsfilename.c_str());
        int n_lines = load_cellcounts_optimized(countsfilename, cell_counts, samples.size());
        fprintf(stderr, "Loaded %d count records for %lu cells\n", n_lines, cell_counts.size());
        
        // Generate .condf if missing and not disabled
        string condf_file = output_prefix + ".condf";
        if (!disable_conditional && !file_exists(condf_file)){
            if (vcf_file.length() == 0 && shared_vcf_name.length() == 0){
                fprintf(stderr, "WARNING: No VCF or shared VCF provided and no .condf file found. "
                                "Skipping conditional match fraction computation.\n");
            }
            else{
                robin_hood::unordered_map<int, ChromSNPs> snpdat_optimized;
                
                if (shared_vcf_name.length() > 0){
                    // Attach to shared memory daemon
                    print_elapsed(start_time, "Generating .condf: attaching to shared VCF...");
                    fprintf(stderr, "Attaching to shared VCF: %s\n", shared_vcf_name.c_str());
                    
                    vector<string> shm_samples;
                    if (!attach_shared_vcf(shared_vcf_name, snpdat_optimized, shm_samples)){
                        fprintf(stderr, "WARNING: Could not attach to shared VCF. Skipping .condf.\n");
                    }
                }
                else{
                    // Load from disk
                    print_elapsed(start_time, "Generating .condf from VCF on disk...");
                    
                    set<string> chroms_vcf;
                    get_vcf_chroms(vcf_file, chroms_vcf);
                    
                    map<string, int> seq2tid_synthetic;
                    int tid_counter = 0;
                    for (const auto& c : chroms_vcf){
                        seq2tid_synthetic[c] = tid_counter++;
                    }
                    
                    int nloaded = read_vcf_chroms_optimized(vcf_file, chroms_vcf,
                        seq2tid_synthetic, snpdat_optimized, vq);
                    fprintf(stderr, "Loaded %d SNPs for .condf computation\n", nloaded);
                }
                
                if (!snpdat_optimized.empty()){
                    // The parallel CONDF function performs genotype-only precompute
                    // internally, avoiding unnecessary pair_targets allocation here.
                    compute_conditional_match_fracs_parallel(snpdat_optimized,
                        conditional_match_fracs, samples.size(), n_threads);
                    
                    FILE* outf = fopen(condf_file.c_str(), "w");
                    dump_exp_fracs(outf, conditional_match_fracs);
                    fclose(outf);
                    fprintf(stderr, "Wrote .condf to %s\n", condf_file.c_str());
                }
            }
        }
    }
    else{
        // Determine chromosomes to process
        set<string> chroms_bam;
        set<string> chroms_vcf;
        set<string> chroms_to_process;
        
        fprintf(stderr, "Identifying shared chromosomes between BAM and VCF...\n");
        get_bam_header_chroms_and_seq2tid(bamfile, chroms_bam, bam_seq2tid);
        fprintf(stderr, "BAM header loaded: %lu chromosomes\n", chroms_bam.size());
        
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
                    fprintf(stderr, "ERROR: No valid individual names found in %s; refusing to ignore -i\n", idfile.c_str());
                    exit(1);
                }
            }
            if (idfile_doublet_given){
                parse_idfile(idfile_doublet, samples, allowed_ids, allowed_ids2, false);
                if (allowed_ids.size() == 0){
                    fprintf(stderr, "ERROR: No valid individual names found in %s; refusing to ignore -I\n", idfile_doublet.c_str());
                    exit(1);
                }
            }

            if (idfile_given || idfile_doublet_given){
                fprintf(stderr, "Identity restriction active: allowed_ids=%lu allowed_ids2=%lu\n",
                    allowed_ids.size(), allowed_ids2.size());
            }
            
            // Write samples file
            string samplesfile = output_prefix + ".samples"; 
            write_samples(samplesfile, samples);
            precompute_all_genotypes(snpdat_optimized, samples.size());
        }
        else{
            // Load VCF into memory (default behavior)
            // Use --no_preload to disable if memory is limited
            if (!no_preload){
                print_elapsed(start_time, "Loading VCF data into memory...");
                fprintf(stderr, "Loading VCF data into memory...\n");
                int nloaded = read_vcf_chroms_optimized(vcf_file, chroms_to_process, 
                    bam_seq2tid, snpdat_optimized, vq);
                print_elapsed(start_time, "VCF loading complete");
                fprintf(stderr, "Loaded %d SNPs from %lu chromosomes\n", 
                    nloaded, chroms_to_process.size());
                precompute_all_genotypes(snpdat_optimized, samples.size());
            }
            else{
                fprintf(stderr, "VCF preloading disabled. Use --shared_vcf for large datasets.\n");
                exit(1);
            }
        }
        
        // Count alleles
        print_elapsed(start_time, "Starting allele counting...");
        fprintf(stderr, "Counting alleles in BAM file...\n");
        
        // WP3 optimization: detect when both interindividual and species panels
        // need counting from the same RNA BAM. If so, load species VCF early,
        // merge SNP sets, and do a single BAM pass via count_alleles_parallel_dual.
        bool use_dual_counting = false;
        if (n_threads > 1 && has_species_vcf 
            && (species_mode == SpeciesPanelMode::AUGMENT 
                || species_mode == SpeciesPanelMode::BOTH
                || species_counts_output)){
            use_dual_counting = true;
        }
        
        if (use_dual_counting){
            // === DUAL-PANEL PATH (WP3) ===
            // Step 1: Build demux_positions for dedup BEFORE species VCF load
            for (auto& kv : snpdat_optimized){
                for (auto& snp : kv.second.snps){
                    demux_positions.insert(make_pair(kv.first, snp.pos));
                }
            }
            fprintf(stderr, "Built demux position set: %lu sites for species dedup\n",
                demux_positions.size());
            
            // Step 2: Load species VCF early
            print_elapsed(start_time, "Loading species panel VCF for dual counting...");
            vector<string> species_samples_early;
            if (species_shared_vcf_name.length() > 0){
                if (!attach_shared_vcf(species_shared_vcf_name, species_snpdat, species_samples_early)){
                    fprintf(stderr, "ERROR: Could not attach to shared species VCF: %s\n",
                        species_shared_vcf_name.c_str());
                    exit(1);
                }
            }
            else{
                read_vcf_samples(species_vcf_file, species_samples_early);
                set<string> sp_chroms;
                get_vcf_chroms(species_vcf_file, sp_chroms);
                int nloaded = read_vcf_chroms_optimized(species_vcf_file, sp_chroms,
                    bam_seq2tid, species_snpdat, vq);
                fprintf(stderr, "Loaded %d species panel SNPs\n", nloaded);
            }
            
            // Step 3: Sample set alignment check
            if (species_samples_early.size() != samples.size()){
                fprintf(stderr, "ERROR: Species VCF has %lu samples but demux VCF has %lu samples. "
                    "Sample sets must match.\n", species_samples_early.size(), samples.size());
                exit(1);
            }
            for (size_t si = 0; si < samples.size(); ++si){
                if (species_samples_early[si] != samples[si]){
                    fprintf(stderr, "ERROR: Species VCF sample[%lu]='%s' differs from demux VCF sample[%lu]='%s'. "
                        "Sample sets must match in same order.\n",
                        si, species_samples_early[si].c_str(), si, samples[si].c_str());
                    exit(1);
                }
            }
            fprintf(stderr, "Species VCF sample set matches demux VCF (%lu samples)\n", samples.size());
            
            // Step 4: Dedup species SNPs against demux positions
            if ((species_mode == SpeciesPanelMode::AUGMENT || species_mode == SpeciesPanelMode::BOTH)
                && !demux_positions.empty()){
                int overlap_count = 0;
                int total_species_sites = 0;
                for (auto& sp_chrom : species_snpdat){
                    int tid = sp_chrom.first;
                    vector<SNPData> filtered_snps;
                    for (auto& sp_snp : sp_chrom.second.snps){
                        total_species_sites++;
                        if (demux_positions.count(make_pair(tid, sp_snp.pos)) > 0){
                            overlap_count++;
                        }
                        else{
                            filtered_snps.push_back(sp_snp);
                        }
                    }
                    sp_chrom.second.snps = filtered_snps;
                }
                if (overlap_count > 0){
                    fprintf(stderr, "Species/demux overlap: %d of %d species sites removed "
                        "(demux VCF retains priority)\n", overlap_count, total_species_sites);
                }
                else{
                    fprintf(stderr, "Species/demux overlap: 0 of %d species sites overlap\n",
                        total_species_sites);
                }
            }
            
            // Step 5: Tag species SNPs with panel_id = 1
            for (auto& kv : species_snpdat){
                for (auto& snp : kv.second.snps){
                    snp.panel_id = 1;
                }
            }
            
            // Step 6: Merge into combined SNP set
            robin_hood::unordered_map<int, ChromSNPs> combined_snpdat = snpdat_optimized;
            for (auto& kv : species_snpdat){
                int tid = kv.first;
                auto& combined_chrom = combined_snpdat[tid];
                for (auto& snp : kv.second.snps){
                    combined_chrom.snps.push_back(snp);
                }
                combined_chrom.sort_snps();
            }
            
            long n_interindiv = 0, n_species = 0;
            for (auto& kv : combined_snpdat){
                for (auto& snp : kv.second.snps){
                    if (snp.panel_id == 0) n_interindiv++;
                    else n_species++;
                }
            }
            fprintf(stderr, "Combined SNP set: %ld interindividual + %ld species = %ld total\n",
                n_interindiv, n_species, n_interindiv + n_species);
            
            // Precompute genotypes and targets for the combined set
            precompute_all_genotypes(combined_snpdat, samples.size());
            
            // Step 7: Single BAM pass with dual output
            robin_hood::unordered_map<unsigned long, AlignedCellCounts> dual_panel0;
            robin_hood::unordered_map<unsigned long, AlignedCellCounts> dual_panel1;
            count_alleles_parallel_dual(bamfile, combined_snpdat,
                dual_panel0, dual_panel1,
                cell_barcodes, samples.size(), n_threads, htslib_threads);
            
            // Step 8: Finalize each panel separately
            finalize_parallel_counts(dual_panel0, cell_counts);
            finalize_parallel_counts(dual_panel1, species_cell_counts_rna);
            
            fprintf(stderr, "Dual counting complete: %lu interindiv cells, %lu species cells\n",
                cell_counts.size(), species_cell_counts_rna.size());
            
            // Step 9: Load panel metadata if needed (for later use)
            if (panel_metadata_file.length() > 0){
                panel_meta = load_panel_metadata(panel_metadata_file, samples);
                panel_meta_loaded = true;
            }
            
            // Mark species RNA counting as done so the later block skips it
            species_counted_dual = true;
        }
        else if (n_threads > 1){
            // Parallel processing
            robin_hood::unordered_map<unsigned long, AlignedCellCounts> parallel_counts;
        
            count_alleles_parallel(bamfile, snpdat_optimized, parallel_counts,
                cell_barcodes, samples.size(), n_threads, htslib_threads,
                dump_pileup, dump_pileup_prefix);
        
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
            compute_conditional_match_fracs_parallel(snpdat_optimized,
                conditional_match_fracs, samples.size(), n_threads);
        }
    
        if (!disable_conditional){
            // For single-threaded path, condf was accumulated during counting
            // and still needs normalization. For parallel path, the new function
            // already normalized internally, so this is a no-op (tots map is empty).
            conditional_match_fracs_normalize(conditional_match_fracs,
                conditional_match_tots, samples.size());
        
            string outname = output_prefix + ".condf";
            FILE* outf = fopen(outname.c_str(), "w");
            dump_exp_fracs(outf, conditional_match_fracs);
            fclose(outf);
        }

        // Build position set from demux VCF for species panel dedup (§5.5)
        if (has_species_vcf){
            for (auto& kv : snpdat_optimized){
                for (auto& snp : kv.second.snps){
                    demux_positions.insert(make_pair(kv.first, snp.pos));
                }
            }
            fprintf(stderr, "Built demux position set: %lu sites for species dedup\n",
                demux_positions.size());
        }
    } // end else (not load_counts)
    
    // Write counts to disk (skip if we loaded from existing file)
    if (!load_counts){
        print_elapsed(start_time, "Writing allele counts to disk...");
        {
            string fname = output_prefix + ".counts";
            gzFile outf = gzopen(fname.c_str(), "w");
            fprintf(stderr, "Writing allele counts to disk...\n");
            dump_cellcounts_optimized(outf, cell_counts, samples.size());
            gzclose(outf);
            fprintf(stderr, "Done writing counts\n");
        }
    }

    // ================================================================
    // 2A: ATAC dual-modality counting
    // ================================================================
    robin_hood::unordered_map<unsigned long, CellCounts> atac_cell_counts;

    if (atac_mode && !load_counts){
        print_elapsed(start_time, "Starting ATAC allele counting...");
        robin_hood::unordered_map<int, ChromSNPs> atac_snpdat;

        if (atac_shared_vcf_name.length() > 0){
            vector<string> atac_samples;
            if (!attach_shared_vcf(atac_shared_vcf_name, atac_snpdat, atac_samples)){
                fprintf(stderr, "ERROR: Could not attach to shared ATAC VCF: %s\n",
                    atac_shared_vcf_name.c_str());
                exit(1);
            }
            fprintf(stderr, "Attached to shared ATAC VCF with %lu chromosomes\n", atac_snpdat.size());
        }
        else{
            set<string> atac_chroms;
            get_vcf_chroms(atac_vcf_file, atac_chroms);
            map<string, int> seq2tid_atac;
            set<string> chroms_atac_header;
            get_bam_header_chroms_and_seq2tid(atac_bamfile, chroms_atac_header, seq2tid_atac);
            int nloaded = read_vcf_chroms_optimized(atac_vcf_file, atac_chroms,
                seq2tid_atac, atac_snpdat, vq);
            fprintf(stderr, "Loaded %d ATAC SNPs\n", nloaded);
        }

        precompute_all_genotypes(atac_snpdat, samples.size());

        if (n_threads > 1){
            robin_hood::unordered_map<unsigned long, AlignedCellCounts> atac_parallel_counts;
            count_alleles_parallel(atac_bamfile, atac_snpdat, atac_parallel_counts,
                cell_barcodes, samples.size(), n_threads, htslib_threads);
            finalize_parallel_counts(atac_parallel_counts, atac_cell_counts);
        }
        else{
            map<pair<int, int>, map<int, float>> atac_cond_fracs, atac_cond_tots;
            count_alleles_single_threaded(atac_bamfile, atac_snpdat, atac_cell_counts,
                cell_barcodes, samples.size(), atac_cond_fracs, atac_cond_tots, false);
        }

        // CB tag intersection check
        int n_intersect = 0;
        int n_rna = (int)cell_counts.size();
        int n_atac = (int)atac_cell_counts.size();
        for (auto& kv : atac_cell_counts){
            if (cell_counts.count(kv.first) > 0) n_intersect++;
        }
        int min_set = (n_rna < n_atac) ? n_rna : n_atac;
        if (min_set > 0){
            double overlap_frac = (double)n_intersect / (double)min_set;
            fprintf(stderr, "ATAC/RNA barcode overlap: %d / %d = %.3f\n",
                n_intersect, min_set, overlap_frac);
            if (overlap_frac < 0.10){
                fprintf(stderr, "ERROR: ATAC/RNA barcode overlap (%.3f) below threshold (0.10).\n"
                    "Check that the ATAC BAM CB tag holds RNA-aligned barcodes.\n", overlap_frac);
                exit(1);
            }
        }

        // Write ATAC counts
        {
            string fname = output_prefix + ".atac.counts";
            gzFile outf = gzopen(fname.c_str(), "w");
            dump_cellcounts_optimized(outf, atac_cell_counts, samples.size());
            gzclose(outf);
            fprintf(stderr, "Wrote ATAC counts for %lu cells\n", atac_cell_counts.size());
        }
    }

    // ================================================================
    // 2C: Species panel loading (needed before assignment)
    // ================================================================
    // NOTE: species_snpdat, species_cell_counts_rna, species_cell_counts_atac,
    // panel_meta, and panel_meta_loaded are declared earlier (before counting)
    // to support the WP3 dual counting optimization path.

    if (has_species_vcf){
        // WP3: if dual counting already loaded species VCF, validated samples,
        // performed dedup, and counted RNA, skip those steps here.
        if (!species_counted_dual){
            // Species panel counts are never cached in .counts, so this block must run
            // even when load_counts is true (main demux counts were loaded from disk).
            print_elapsed(start_time, "Loading species panel VCF...");
            vector<string> species_samples;  // Sample names from species VCF
            if (species_shared_vcf_name.length() > 0){
                if (!attach_shared_vcf(species_shared_vcf_name, species_snpdat, species_samples)){
                    fprintf(stderr, "ERROR: Could not attach to shared species VCF: %s\n",
                        species_shared_vcf_name.c_str());
                    exit(1);
                }
            }
            else if (!load_counts){
                // Disk-based species VCF requires the BAM reader for tid mapping,
                // which is only initialized when load_counts is false.
                read_vcf_samples(species_vcf_file, species_samples);
                set<string> sp_chroms;
                get_vcf_chroms(species_vcf_file, sp_chroms);
                int nloaded = read_vcf_chroms_optimized(species_vcf_file, sp_chroms,
                    bam_seq2tid, species_snpdat, vq);
                fprintf(stderr, "Loaded %d species panel SNPs\n", nloaded);
            }
            else{
                // load_counts + disk-based species VCF: not currently supported
                fprintf(stderr, "ERROR: Disk-based --species_vcf with cached .counts is not supported. "
                    "Use --species_shared_vcf or delete the .counts file to re-run from BAM.\n");
                exit(1);
            }

            // §5.3: Sample set alignment check
            if (species_samples.size() != samples.size()){
                fprintf(stderr, "ERROR: Species VCF has %lu samples but demux VCF has %lu samples. "
                    "Sample sets must match.\n", species_samples.size(), samples.size());
                exit(1);
            }
            for (size_t si = 0; si < samples.size(); ++si){
                if (species_samples[si] != samples[si]){
                    fprintf(stderr, "ERROR: Species VCF sample[%lu]='%s' differs from demux VCF sample[%lu]='%s'. "
                        "Sample sets must match in same order.\n",
                        si, species_samples[si].c_str(), si, samples[si].c_str());
                    exit(1);
                }
            }
            fprintf(stderr, "Species VCF sample set matches demux VCF (%lu samples)\n", samples.size());

            // §5.5: SNP overlap dedup - remove species sites that overlap with demux VCF
            // Demux VCF retains priority at overlapping positions.
            // When load_counts is true, demux_positions is empty (demux SNP data was not
            // loaded), so dedup is skipped. This is safe: species and demux panels are
            // designed to be non-overlapping, and any rare overlap would only cause minor
            // double-counting in the species-only scoring pass.
            if ((species_mode == SpeciesPanelMode::AUGMENT || species_mode == SpeciesPanelMode::BOTH)
                && !demux_positions.empty()){
                int overlap_count = 0;
                int total_species_sites = 0;
                for (auto& sp_chrom : species_snpdat){
                    int tid = sp_chrom.first;
                    vector<SNPData> filtered_snps;
                    for (auto& sp_snp : sp_chrom.second.snps){
                        total_species_sites++;
                        if (demux_positions.count(make_pair(tid, sp_snp.pos)) > 0){
                            overlap_count++;
                        }
                        else{
                            filtered_snps.push_back(sp_snp);
                        }
                    }
                    sp_chrom.second.snps = filtered_snps;
                }
                if (overlap_count > 0){
                    fprintf(stderr, "Species/demux overlap: %d of %d species sites removed "
                        "(demux VCF retains priority)\n", overlap_count, total_species_sites);
                }
                else{
                    fprintf(stderr, "Species/demux overlap: 0 of %d species sites overlap\n",
                        total_species_sites);
                }
            }

            // Load panel metadata if needed
            if (panel_metadata_file.length() > 0){
                panel_meta = load_panel_metadata(panel_metadata_file, samples);
                panel_meta_loaded = true;
            }

            // Count reads at species-panel sites.
            // Always needed for augment/both (LLR augmentation).
            // Also needed for filter mode when --species_counts_output requests persisting counts.
            if (species_mode == SpeciesPanelMode::AUGMENT || species_mode == SpeciesPanelMode::BOTH
                || species_counts_output){
                precompute_all_genotypes(species_snpdat, samples.size());
                print_elapsed(start_time, "Counting alleles at species panel sites (RNA)...");
                if (n_threads > 1){
                    robin_hood::unordered_map<unsigned long, AlignedCellCounts> sp_parallel;
                    count_alleles_parallel(bamfile, species_snpdat, sp_parallel,
                        cell_barcodes, samples.size(), n_threads, htslib_threads);
                    finalize_parallel_counts(sp_parallel, species_cell_counts_rna);
                }
                else{
                    map<pair<int, int>, map<int, float>> sp_cond, sp_tots;
                    count_alleles_single_threaded(bamfile, species_snpdat, species_cell_counts_rna,
                        cell_barcodes, samples.size(), sp_cond, sp_tots, false);
                }
                fprintf(stderr, "Species RNA counts for %lu cells\n", species_cell_counts_rna.size());

                if (atac_mode){
                    print_elapsed(start_time, "Counting alleles at species panel sites (ATAC)...");
                    if (n_threads > 1){
                        robin_hood::unordered_map<unsigned long, AlignedCellCounts> sp_atac_parallel;
                        count_alleles_parallel(atac_bamfile, species_snpdat, sp_atac_parallel,
                            cell_barcodes, samples.size(), n_threads, htslib_threads);
                        finalize_parallel_counts(sp_atac_parallel, species_cell_counts_atac);
                    }
                    else{
                        map<pair<int, int>, map<int, float>> sp_cond2, sp_tots2;
                        count_alleles_single_threaded(atac_bamfile, species_snpdat, species_cell_counts_atac,
                            cell_barcodes, samples.size(), sp_cond2, sp_tots2, false);
                    }
                }
            }
        }
        else{
            // Dual counting path already handled VCF loading, validation, dedup, 
            // and RNA counting. Only ATAC counting remains (if needed).
            fprintf(stderr, "Species RNA counting already done via dual-panel path\n");
            
            if (atac_mode && (species_mode == SpeciesPanelMode::AUGMENT
                || species_mode == SpeciesPanelMode::BOTH || species_counts_output)){
                print_elapsed(start_time, "Counting alleles at species panel sites (ATAC)...");
                if (n_threads > 1){
                    robin_hood::unordered_map<unsigned long, AlignedCellCounts> sp_atac_parallel;
                    count_alleles_parallel(atac_bamfile, species_snpdat, sp_atac_parallel,
                        cell_barcodes, samples.size(), n_threads, htslib_threads);
                    finalize_parallel_counts(sp_atac_parallel, species_cell_counts_atac);
                }
                else{
                    map<pair<int, int>, map<int, float>> sp_cond2, sp_tots2;
                    count_alleles_single_threaded(atac_bamfile, species_snpdat, species_cell_counts_atac,
                        cell_barcodes, samples.size(), sp_cond2, sp_tots2, false);
                }
            }
        }
    }

    // ================================================================
    // Write native species counts / condf / samples to disk (if requested)
    // ================================================================
    robin_hood::unordered_map<unsigned long, CellCounts> species_cell_counts_native;
    vector<string> species_samples_native;

    if ((species_counts_output || species_assignment_output) && !species_cell_counts_rna.empty()){
        load_panel_metadata_if_needed(panel_metadata_file, samples, panel_meta, panel_meta_loaded);
        species_samples_native = panel_meta.species_list;
        aggregate_cellcounts_to_species_native(species_cell_counts_rna,
            species_cell_counts_native, panel_meta, samples.size());
    }

    if (species_counts_output && !species_cell_counts_native.empty()){
        {
            string fname = output_prefix + ".species_counts";
            gzFile outf = gzopen(fname.c_str(), "w");
            dump_cellcounts_optimized(outf, species_cell_counts_native,
                species_samples_native.size());
            gzclose(outf);
            fprintf(stderr, "Wrote native species counts for %lu cells to %s (%lu species columns)\n",
                species_cell_counts_native.size(), fname.c_str(), species_samples_native.size());
        }

        {
            map<pair<int, int>, map<int, float>> sp_condf;
            compute_species_condf_native(species_snpdat, sp_condf, panel_meta, samples.size());
            string fname = output_prefix + ".species_condf";
            FILE* outf = fopen(fname.c_str(), "w");
            dump_exp_fracs(outf, sp_condf);
            fclose(outf);
            fprintf(stderr, "Wrote native species conditional match fracs to %s\n", fname.c_str());
        }

        {
            string fname = output_prefix + ".species_samples";
            write_samples(fname, species_samples_native);
            fprintf(stderr, "Wrote native species sample list to %s\n", fname.c_str());
        }
    }

    // ================================================================
    // Step 0a: --skip_assignment bypass
    // ================================================================
    if (skip_assignment){
        print_elapsed(start_time, "Skipping identity assignment (--skip_assignment set)");
        fprintf(stderr, "Wrote %s.counts; exiting (no assignment performed).\n",
                output_prefix.c_str());
        if (identity_prior_file.length() > 0){
            fprintf(stderr, "WARNING: --identity_prior was loaded but unused (--skip_assignment)\n");
        }
        return 0;
    }

    // ================================================================
    // 2B: Load identity prior if provided
    // ================================================================
    map<int, double> identity_prior_map;
    map<int, double> identity_prior_conc;
    double z_doublet_prior = 0.0;
    bool has_identity_prior = false;

    if (identity_prior_file.length() > 0){
        print_elapsed(start_time, "Loading identity prior...");
        load_contam_prof(identity_prior_file, identity_prior_map, identity_prior_conc,
            samples, false);  // false: missing samples get warning, not error
        has_identity_prior = true;

        // Compute Z_doublet for doublet prior normalization
        // Sum over allowed doublets of 2 * pi[i] * pi[j]
        set<int>& aa = allowed_ids;
        for (int i = 0; i < (int)samples.size() - 1; ++i){
            if (!aa.empty() && aa.find(i) == aa.end()) continue;
            for (int j = i + 1; j < (int)samples.size(); ++j){
                if (!aa.empty() && aa.find(j) == aa.end()) continue;
                int k = hap_comb_to_idx(i, j, samples.size());
                if (!aa.empty() && aa.find(k) == aa.end()) continue;

                double pi_i = 0.0, pi_j = 0.0;
                if (identity_prior_map.count(i) > 0) pi_i = identity_prior_map[i];
                if (identity_prior_map.count(j) > 0) pi_j = identity_prior_map[j];
                z_doublet_prior += 2.0 * pi_i * pi_j;
            }
        }
        fprintf(stderr, "Identity prior: %lu entries, Z_doublet=%.6f\n",
            identity_prior_map.size(), z_doublet_prior);
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

    // NEW: Final assignment with diagnostic collection
    robin_hood::unordered_map<unsigned long, CellDiagnostics> cell_diagnostics;
    robin_hood::unordered_map<unsigned long, vector<RunnerUp> > cell_runner_ups;
    
    // Het data structures
    robin_hood::unordered_map<unsigned long, CellHetData> het_data;
    robin_hood::unordered_map<int, ChromSNPs> het_snpdat;
    vector<pair<int, int>> idx_to_site;  // Site index to (tid, pos) for persite method
    
    bool het_vcf_available = (het_vcf_file.length() > 0 || shared_het_vcf_name.length() > 0);
    int n_het_loaded = 0;
    
    if (het_vcf_available && write_diagnostics){
        if (load_counts && bamfile.length() == 0){
            fprintf(stderr, "WARNING: Het VCF diagnostics require a BAM file. Skipping het balance "
                            "computation in counts-only mode.\n");
            fprintf(stderr, "Re-run with --bam to include ploidy diagnostics.\n");
        }
        else{
            vector<string> het_samples;  // Not used but needed for attach_shared_vcf
            
            if (shared_het_vcf_name.length() > 0){
                // Load from shared memory
                print_elapsed(start_time, "Attaching to shared het VCF for ploidy diagnostics...");
                fprintf(stderr, "Attaching to shared het VCF: %s\n", shared_het_vcf_name.c_str());
                
                if (!attach_shared_vcf(shared_het_vcf_name, het_snpdat, het_samples)){
                    fprintf(stderr, "WARNING: Could not attach to shared het VCF: %s\n", shared_het_vcf_name.c_str());
                    fprintf(stderr, "Ploidy-related diagnostics will not be computed.\n");
                }
                else{
                    // Count SNPs loaded
                    for (auto& kv : het_snpdat){
                        n_het_loaded += kv.second.snps.size();
                    }
                    fprintf(stderr, "Attached to shared het VCF with %d sites\n", n_het_loaded);
                }
            }
            else{
                // Load from file
                print_elapsed(start_time, "Loading het VCF for ploidy diagnostics...");
                fprintf(stderr, "Loading het VCF: %s\n", het_vcf_file.c_str());
                
                map<string, int> seq2tid_het;
                
                // Get chromosome names from BAM header
                set<string> chroms_for_het;
                get_bam_header_chroms_and_seq2tid(bamfile, chroms_for_het, seq2tid_het);
                
                n_het_loaded = load_het_vcf(het_vcf_file, chroms_for_het, seq2tid_het, het_snpdat, vq);
                fprintf(stderr, "Loaded %d het sites\n", n_het_loaded);
            }
            
            if (n_het_loaded > 0){
                const char* method_name = (het_method == HetBalanceMethod::PERSITE) ? "per-site" : "Welford";
                fprintf(stderr, "Using %s method (min_het_sites=%d", method_name, min_het_sites);
                if (het_method == HetBalanceMethod::PERSITE) {
                    fprintf(stderr, ", min_het_depth=%.1f", min_het_depth);
                }
                fprintf(stderr, ")\n");
                print_elapsed(start_time, "Counting alleles at het sites...");
                count_het_alleles_extended(bamfile, het_snpdat, het_data, idx_to_site,
                    cell_barcodes, samples.size(), n_threads, htslib_threads, het_method);
            }
            else{
                fprintf(stderr, "WARNING: No het sites available\n");
            }
        }
    }
    else if (write_diagnostics && !het_vcf_available){
        fprintf(stderr, "WARNING: --het_vcf/--shared_het_vcf not provided. Ploidy-related diagnostics "
                        "(het_balance_var, n_het_sites, het_total_depth) will not be computed.\n"
                        "For full tetraploid analysis, provide het VCF from downsample_vcf_parallel.\n");
    }

    // ================================================================
    // 3.6: ATAC het balance diagnostics
    // ================================================================
    robin_hood::unordered_map<unsigned long, CellHetData> atac_het_data;
    robin_hood::unordered_map<int, ChromSNPs> atac_het_snpdat;
    vector<pair<int, int>> atac_idx_to_site;
    bool atac_het_available = atac_mode &&
        (atac_het_vcf_file.length() > 0 || atac_shared_het_vcf_name.length() > 0);
    int n_atac_het_loaded = 0;

    if (atac_het_available && write_diagnostics && !load_counts){
        if (atac_shared_het_vcf_name.length() > 0){
            print_elapsed(start_time, "Attaching to shared ATAC het VCF...");
            vector<string> atac_het_samples;
            if (!attach_shared_vcf(atac_shared_het_vcf_name, atac_het_snpdat, atac_het_samples)){
                fprintf(stderr, "WARNING: Could not attach to shared ATAC het VCF: %s\n",
                    atac_shared_het_vcf_name.c_str());
            }
            else{
                for (auto& kv : atac_het_snpdat){
                    n_atac_het_loaded += kv.second.snps.size();
                }
                fprintf(stderr, "Attached to shared ATAC het VCF with %d sites\n", n_atac_het_loaded);
            }
        }
        else{
            print_elapsed(start_time, "Loading ATAC het VCF...");
            map<string, int> seq2tid_atac_het;
            set<string> chroms_atac_het;
            get_bam_header_chroms_and_seq2tid(atac_bamfile, chroms_atac_het, seq2tid_atac_het);
            n_atac_het_loaded = load_het_vcf(atac_het_vcf_file, chroms_atac_het,
                seq2tid_atac_het, atac_het_snpdat, vq);
            fprintf(stderr, "Loaded %d ATAC het sites\n", n_atac_het_loaded);
        }

        if (n_atac_het_loaded > 0){
            print_elapsed(start_time, "Counting alleles at ATAC het sites...");
            count_het_alleles_extended(atac_bamfile, atac_het_snpdat, atac_het_data,
                atac_idx_to_site, cell_barcodes, samples.size(), n_threads,
                htslib_threads, het_method);
            fprintf(stderr, "ATAC het data for %lu cells\n", atac_het_data.size());
        }
    }
    
    // Prepare extended evidence pointers for diagnostic calls
    robin_hood::unordered_map<unsigned long, CellCounts>* atac_ptr =
        atac_mode ? &atac_cell_counts : nullptr;
    const map<int, double>* prior_ptr =
        has_identity_prior ? &identity_prior_map : nullptr;
    robin_hood::unordered_map<unsigned long, CellCounts>* sp_rna_ptr =
        (species_mode == SpeciesPanelMode::AUGMENT || species_mode == SpeciesPanelMode::BOTH)
        ? &species_cell_counts_rna : nullptr;
    robin_hood::unordered_map<unsigned long, CellCounts>* sp_atac_ptr =
        (atac_mode && sp_rna_ptr != nullptr) ? &species_cell_counts_atac : nullptr;

    // Do final assignment with diagnostics if requested
    if (write_diagnostics){
        print_elapsed(start_time, "Final assignment with diagnostic collection...");
        
        if (het_data.empty()){
            // No het data - use original function
            assign_ids_parallel_with_diagnostics(
                cell_counts, samples, assn, assn_llr,
                allowed_ids, allowed_ids2, doublet_rate, error_ref_posterior, error_alt_posterior,
                false, prior_weights, n_threads, n_target,
                true,  // compute_diagnostics
                n_runner_ups,
                close_threshold,
                NULL,  // no het counts
                cell_diagnostics,
                cell_runner_ups,
                atac_ptr, prior_ptr, z_doublet_prior,
                sp_rna_ptr, sp_atac_ptr);
        }
        else{
            // Use extended method (Welford or per-site)
            assign_ids_parallel_with_diagnostics_extended(
                cell_counts, samples, assn, assn_llr,
                allowed_ids, allowed_ids2, doublet_rate, error_ref_posterior, error_alt_posterior,
                false, prior_weights, n_threads, n_target,
                true,  // compute_diagnostics
                n_runner_ups,
                close_threshold,
                &het_data, &het_snpdat, &idx_to_site,
                het_method, min_het_sites, min_het_depth,
                cell_diagnostics,
                cell_runner_ups,
                atac_ptr, prior_ptr, z_doublet_prior,
                sp_rna_ptr, sp_atac_ptr);
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
    
    // NEW: Write diagnostic files
    if (write_diagnostics){
        // Write .diagnostics.gz
        {
            string fname = output_prefix + ".diagnostics.gz";
            write_diagnostics_gz(fname, samples, assn, assn_llr, cell_diagnostics,
                samples.size(), barcode_group, cellranger, seurat, underscore,
                atac_het_available ? &atac_het_data : NULL,
                atac_het_available ? &atac_het_snpdat : NULL,
                atac_het_available ? &atac_idx_to_site : NULL,
                het_method, min_het_sites, min_het_depth,
                dump_selection_audit);
        }
        
        // Write .runner_ups.gz
        {
            string fname = output_prefix + ".runner_ups.gz";
            write_runner_ups_gz(fname, samples, assn, cell_runner_ups,
                samples.size(), barcode_group, cellranger, seurat, underscore);
        }
        
        fprintf(stderr, "Diagnostic output complete.\n");
        if (!het_vcf_available){
            fprintf(stderr, "Note: het_balance_var columns contain -1 (no --het_vcf provided)\n");
        }
    }
    
    // ================================================================
    // 2C: Native species assignment output
    // ================================================================
    // This is deliberately independent from .assignments.  .assignments remains
    // the individual-native demux call from the individual SNP panel.  The species
    // panel writes its own species-native assignment file with species labels only.
    if (species_assignment_output){
        if (species_cell_counts_native.empty() && !species_cell_counts_rna.empty()){
            load_panel_metadata_if_needed(panel_metadata_file, samples, panel_meta, panel_meta_loaded);
            species_samples_native = panel_meta.species_list;
            aggregate_cellcounts_to_species_native(species_cell_counts_rna,
                species_cell_counts_native, panel_meta, samples.size());
        }
        write_native_species_assignments(output_prefix, species_cell_counts_native,
            species_samples_native, doublet_rate, error_ref_posterior, error_alt_posterior,
            n_threads, n_target, barcode_group, cellranger, seurat, underscore);
    }

    print_elapsed(start_time, "Complete!");
    fprintf(stderr, "Done!\n");
    
    return 0;
}
