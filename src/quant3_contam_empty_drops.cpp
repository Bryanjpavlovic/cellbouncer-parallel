// ============================================================================
// quant3_contam_empty_drops.cpp
// Estimate ambient RNA mixture profile (pi) from empty-droplet allele counts.
//
// Produces .contam_prof_empty (per-individual) or .species_prof_empty
// (per-species, with --aggregate_to_species) files consumed by
// quant3_contam_ap --fixed_ambient_prof or --species_mode fixed --species_prior.
//
// Algorithm:
//   1. Load .counts file (same format as demux_parallel output)
//   2. Load .condf (conditional matching probabilities)
//   3. Aggregate all cell allele counts into a single pseudo-bulk
//   4. Set c=1 for the pseudo-bulk (all RNA is "ambient" in empty drops)
//   5. Solve for pi using contamFinder3 in bulk mode
//   6. Bootstrap pi for Dirichlet concentration parameters
//   7. Write output
// ============================================================================

#include <getopt.h>
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
#include <zlib.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "ambient_rna_three_ap.h"
#include "demux_vcf_io.h"

using std::cout;
using std::endl;
using namespace std;

const string TOOL_VERSION = "1.0";

void help(int code){
    fprintf(stderr, "quant3_contam_empty_drops [OPTIONS]\n");
    fprintf(stderr, "Estimates the ambient RNA mixture profile (pi) from allele counts.\n");
    fprintf(stderr, "Uses all cells as a pseudo-bulk to estimate the mixture of individuals\n");
    fprintf(stderr, "contributing to the ambient RNA pool.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --output_prefix -o The output prefix used in a prior run of demux_parallel\n");
    fprintf(stderr, "===== OPTIONAL =====\n");
    fprintf(stderr, "    --output -O Output file path (default: <output_prefix>.contam_prof_empty\n");
    fprintf(stderr, "       or <output_prefix>.species_prof_empty with --aggregate_to_species)\n");
    fprintf(stderr, "    --num_threads -T Number of parallel threads (default = all available)\n");
    fprintf(stderr, "    --n_bootstrap -b Number of bootstrap replicates (default = 100)\n");
    fprintf(stderr, "    --error_ref -e Reference error rate (default 0.001)\n");
    fprintf(stderr, "    --error_alt -E Alt error rate (default 0.001)\n");
    fprintf(stderr, "    --ids -i Restrict to listed individuals (one name per line)\n");
    fprintf(stderr, "    --libname -n Library name for diagnostics output\n");
    fprintf(stderr, "    --condf -F Pre-computed .condf file (if absent, loaded from output_prefix)\n");
    fprintf(stderr, "    --aggregate_to_species  Compute one pi per species instead of per individual.\n");
    fprintf(stderr, "       Output file extension becomes .species_prof_empty.\n");
    fprintf(stderr, "       Requires --panel_metadata.\n");
    fprintf(stderr, "    --panel_metadata -P TSV mapping individual to species (header required)\n");
    fprintf(stderr, "    --species_counts FILE  Species-diagnostic allele counts file.\n");
    fprintf(stderr, "       Same format as .counts, at species-diagnostic SNP sites.\n");
    fprintf(stderr, "       Used with --aggregate_to_species for cleaner species separation.\n");
    fprintf(stderr, "       Requires a companion .species_condf file alongside.\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

int main(int argc, char *argv[]){

    static struct option long_options[] = {
        {"output_prefix", required_argument, 0, 'o'},
        {"output", required_argument, 0, 'O'},
        {"num_threads", required_argument, 0, 'T'},
        {"n_bootstrap", required_argument, 0, 'b'},
        {"error_ref", required_argument, 0, 'e'},
        {"error_alt", required_argument, 0, 'E'},
        {"ids", required_argument, 0, 'i'},
        {"libname", required_argument, 0, 'n'},
        {"condf", required_argument, 0, 'F'},
        {"aggregate_to_species", no_argument, 0, 1001},
        {"panel_metadata", required_argument, 0, 'P'},
        {"species_counts", required_argument, 0, 1002},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    string output_prefix = "";
    string output_file = "";
    int num_threads = 0;
    int n_bootstrap = 100;
    double error_ref = 0.001;
    double error_alt = 0.001;
    string idfile = "";
    bool idfile_given = false;
    string libname = "";
    string condf_file = "";
    bool aggregate_to_species = false;
    string panel_metadata_file = "";
    string species_counts_file = "";

    int option_index = 0;
    int ch;

    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "o:O:T:b:e:E:i:n:F:P:h",
            long_options, &option_index)) != -1){
        switch(ch){
            case 'o':
                output_prefix = optarg;
                break;
            case 'O':
                output_file = optarg;
                break;
            case 'T':
                num_threads = atoi(optarg);
                break;
            case 'b':
                n_bootstrap = atoi(optarg);
                break;
            case 'e':
                error_ref = atof(optarg);
                break;
            case 'E':
                error_alt = atof(optarg);
                break;
            case 'i':
                idfile = optarg;
                idfile_given = true;
                break;
            case 'n':
                libname = optarg;
                break;
            case 'F':
                condf_file = optarg;
                break;
            case 1001:
                aggregate_to_species = true;
                break;
            case 'P':
                panel_metadata_file = optarg;
                break;
            case 1002:
                species_counts_file = optarg;
                break;
            case 'h':
                help(0);
                break;
            default:
                help(1);
                break;
        }
    }

    // Validate arguments
    if (output_prefix.empty()){
        fprintf(stderr, "ERROR: --output_prefix/-o is required\n");
        exit(1);
    }
    if (aggregate_to_species && panel_metadata_file.empty()){
        fprintf(stderr, "ERROR: --aggregate_to_species requires --panel_metadata\n");
        exit(1);
    }
    if (!species_counts_file.empty() && !aggregate_to_species){
        fprintf(stderr, "ERROR: --species_counts requires --aggregate_to_species\n");
        exit(1);
    }
    if (num_threads <= 1){
        num_threads = 0;
    }

    // Default output file
    if (output_file.empty()){
        if (aggregate_to_species){
            output_file = output_prefix + ".species_prof_empty";
        } else {
            output_file = output_prefix + ".contam_prof_empty";
        }
    }

    fprintf(stderr, "quant3_contam_empty_drops v%s\n", TOOL_VERSION.c_str());
    fprintf(stderr, "  Output prefix: %s\n", output_prefix.c_str());
    fprintf(stderr, "  Output file: %s\n", output_file.c_str());
    if (aggregate_to_species){
        fprintf(stderr, "  Mode: aggregate to species\n");
    }

    // ---- Load samples ----
    string sample_name = output_prefix + ".samples";
    vector<string> samples;
    if (file_exists(sample_name)){
        load_samples(sample_name, samples);
    } else {
        fprintf(stderr, "ERROR: no .samples file found for %s\n", output_prefix.c_str());
        exit(1);
    }
    fprintf(stderr, "  Loaded %lu samples\n", samples.size());

    // ---- Load conditional matching probabilities ----
    map<pair<int, int>, map<int, float> > exp_match_fracs;
    string expfrac_name = condf_file.empty() ? output_prefix + ".condf" : condf_file;
    if (file_exists(expfrac_name)){
        load_exp_fracs(expfrac_name, exp_match_fracs);
    } else {
        fprintf(stderr, "ERROR: no .condf file found: %s\n", expfrac_name.c_str());
        exit(1);
    }

    // ---- Load filtered ID list, if given ----
    set<int> allowed_ids;
    set<int> allowed_ids2;
    if (idfile_given){
        parse_idfile(idfile, samples, allowed_ids, allowed_ids2, true);
    }

    // ---- Load allele counts ----
    robin_hood::unordered_map<unsigned long, map<pair<int, int>,
        map<pair<int, int>, pair<float, float> > > > indv_allelecounts;
    string counts_name = output_prefix + ".counts";
    if (file_exists(counts_name)){
        fprintf(stderr, "Loading counts from %s...\n", counts_name.c_str());
        load_counts_from_file(indv_allelecounts, samples, counts_name, allowed_ids);
    } else {
        fprintf(stderr, "ERROR: no .counts file found: %s\n", counts_name.c_str());
        exit(1);
    }
    fprintf(stderr, "  Loaded counts for %lu cells/barcodes\n", indv_allelecounts.size());

    // ---- Load assignments ----
    // For the pseudo-bulk, we need cell->identity assignments to set up the
    // contamFinder3 (it requires assn). Load from .assignments file.
    robin_hood::unordered_map<unsigned long, int> assn;
    robin_hood::unordered_map<unsigned long, double> assn_llr;
    string assn_name = output_prefix + ".assignments";
    if (file_exists(assn_name)){
        load_assignments_from_file(assn_name, assn, assn_llr, samples);
    } else {
        fprintf(stderr, "ERROR: no .assignments file found: %s\n", assn_name.c_str());
        exit(1);
    }

    // ---- Create contamFinder3 in bulk mode ----
    fprintf(stderr, "Setting up bulk-mode contamFinder3...\n");
    contamFinder3 cf(indv_allelecounts, assn, assn_llr, exp_match_fracs,
        samples.size(), allowed_ids, allowed_ids, allowed_ids2);
    cf.set_error_rates(error_ref, error_alt);
    cf.set_num_threads(num_threads);
    cf.set_bulk_mode(true);
    cf.no_reassign();

    // Species mode for species-aggregated estimation
    PanelMetadata pm;
    if (aggregate_to_species){
        pm = load_panel_metadata(panel_metadata_file, samples);
        cf.set_species_mode(pm);
        fprintf(stderr, "  Species-aggregated mode: %lu species\n", pm.species_list.size());
    }

    // Species-diagnostic counts for species-level solver
    if (!species_counts_file.empty()){
        robin_hood::unordered_map<unsigned long, map<pair<int, int>,
            map<pair<int, int>, pair<float, float> > > > sp_counts;
        if (file_exists(species_counts_file)){
            fprintf(stderr, "Loading species-diagnostic counts from %s...\n",
                species_counts_file.c_str());
            load_counts_from_file(sp_counts, samples, species_counts_file, allowed_ids);
            fprintf(stderr, "  Loaded species counts for %lu cells\n", sp_counts.size());
        } else {
            fprintf(stderr, "ERROR: species counts file not found: %s\n",
                species_counts_file.c_str());
            exit(1);
        }

        // Load companion .species_condf
        string sp_condf_name = species_counts_file;
        size_t pos = sp_condf_name.rfind(".species_counts");
        if (pos != string::npos){
            sp_condf_name = sp_condf_name.substr(0, pos) + ".species_condf";
        } else {
            sp_condf_name = sp_condf_name + ".species_condf";
        }
        map<pair<int, int>, map<int, float> > sp_expfracs;
        if (file_exists(sp_condf_name)){
            load_exp_fracs(sp_condf_name, sp_expfracs);
            fprintf(stderr, "  Loaded species condf from %s\n", sp_condf_name.c_str());
        } else {
            fprintf(stderr, "ERROR: species condf file not found: %s\n",
                sp_condf_name.c_str());
            exit(1);
        }

        cf.set_species_counts(sp_counts, sp_expfracs);
    }

    // ---- Solve for pi ----
    // In bulk mode, we fix c=1 and solve only for the mixture proportions.
    // The contamFinder3::fit() method handles bulk mode internally:
    // it skips per-cell estimation and calls update_amb_prof_mixture once.
    fprintf(stderr, "Estimating ambient profile...\n");

    // Set c_init to 1.0 for bulk mode (all RNA is ambient)
    cf.set_init_c(1.0);
    cf.fit();

    fprintf(stderr, "Ambient profile estimation complete.\n");

    // Print estimated profile
    for (int i = 0; i < (int)samples.size(); i++){
        if (cf.contam_prof.count(i) > 0){
            fprintf(stderr, "  %s: %.6f\n", samples[i].c_str(), cf.contam_prof[i]);
        }
    }

    // ---- Bootstrap ----
    map<int, double> contam_prof_conc;
    if (n_bootstrap > 0){
        fprintf(stderr, "Bootstrapping (%d replicates)...\n", n_bootstrap);
        cf.bootstrap_amb_prof(n_bootstrap, contam_prof_conc);
        fprintf(stderr, "Bootstrap complete.\n");
    }

    // ---- Write output ----
    if (aggregate_to_species){
        // Species-level output
        map<string, double> species_prof;
        map<string, double> species_prof_conc;

        // If species mode was active, species_contam_prof has the results
        if (!cf.species_contam_prof.empty()){
            species_prof = cf.species_contam_prof;
        } else {
            // Aggregate individual-level to species-level
            for (const auto& sp : pm.species_list){
                double sum = 0.0;
                if (pm.species_to_sample_indices.count(sp) > 0){
                    for (int idx : pm.species_to_sample_indices.at(sp)){
                        if (cf.contam_prof.count(idx) > 0){
                            sum += cf.contam_prof[idx];
                        }
                    }
                }
                species_prof[sp] = sum;
            }
        }

        // Aggregate bootstrap concentrations to species level if available
        // (simplified: use individual-level Dirichlet params summed within species)
        for (const auto& sp : pm.species_list){
            double alpha_sum = 0.0;
            if (pm.species_to_sample_indices.count(sp) > 0){
                for (int idx : pm.species_to_sample_indices.at(sp)){
                    if (contam_prof_conc.count(idx) > 0){
                        alpha_sum += contam_prof_conc[idx];
                    }
                }
            }
            if (alpha_sum > 0){
                species_prof_conc[sp] = alpha_sum;
            }
        }

        FILE* outf = fopen(output_file.c_str(), "w");
        if (!outf){
            fprintf(stderr, "ERROR: cannot open output file: %s\n", output_file.c_str());
            exit(1);
        }
        dump_species_prof(outf, species_prof, species_prof_conc);
        fclose(outf);
        fprintf(stderr, "Species profile written to: %s\n", output_file.c_str());

        // Print species profile
        for (const auto& sp : species_prof){
            fprintf(stderr, "  %s: %.6f", sp.first.c_str(), sp.second);
            if (species_prof_conc.count(sp.first) > 0){
                fprintf(stderr, " (alpha=%.2f)", species_prof_conc[sp.first]);
            }
            fprintf(stderr, "\n");
        }
    } else {
        // Individual-level output
        FILE* outf = fopen(output_file.c_str(), "w");
        if (!outf){
            fprintf(stderr, "ERROR: cannot open output file: %s\n", output_file.c_str());
            exit(1);
        }
        dump_contam_prof(outf, cf.contam_prof, contam_prof_conc, samples);
        fclose(outf);
        fprintf(stderr, "Individual profile written to: %s\n", output_file.c_str());
    }

    // ---- Write diagnostics ----
    {
        string diag_name = output_prefix + ".empty_diagnostics.tsv";
        FILE* diagf = fopen(diag_name.c_str(), "w");
        if (diagf){
            fprintf(diagf, "metric\tvalue\n");
            fprintf(diagf, "n_barcodes\t%lu\n", indv_allelecounts.size());
            fprintf(diagf, "n_samples\t%d\n", (int)samples.size());
            fprintf(diagf, "n_bootstrap\t%d\n", n_bootstrap);
            fprintf(diagf, "error_ref\t%.6f\n", error_ref);
            fprintf(diagf, "error_alt\t%.6f\n", error_alt);
            if (!libname.empty()){
                fprintf(diagf, "libname\t%s\n", libname.c_str());
            }
            if (aggregate_to_species){
                fprintf(diagf, "mode\taggregate_to_species\n");
                fprintf(diagf, "n_species\t%lu\n", pm.species_list.size());
            } else {
                fprintf(diagf, "mode\tindividual\n");
            }

            // Per-individual proportions
            for (int i = 0; i < (int)samples.size(); i++){
                if (cf.contam_prof.count(i) > 0){
                    fprintf(diagf, "pi_%s\t%.10f\n", samples[i].c_str(), cf.contam_prof[i]);
                }
            }
            fclose(diagf);
            fprintf(stderr, "Diagnostics written to: %s\n", diag_name.c_str());
        }
    }

    fprintf(stderr, "Done.\n");
    return 0;
}
