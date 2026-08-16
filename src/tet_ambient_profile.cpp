// ============================================================================
// tet_ambient_profile.cpp
// Estimate ambient RNA mixture profile (pi) from truly empty droplets.
//
// Replaces quant3_contam_empty_drops with three improvements:
//   1. Barcode filtering: exclude cell-containing droplets via --filtered_barcodes
//   2. Pre-aggregation: collapse barcodes into per-identity pseudo-bulk entries
//      before solving (runtime drops from ~16h to <1 min)
//   3. Explicit --interindividual / --interspecies panel selection (mutually
//      exclusive) replaces --aggregate_to_species + --species_counts
//
// Algorithm:
//   1. Load .counts (or .species_counts) and .condf (or .species_condf)
//   2. Filter out cell-containing barcodes; optionally subsample empties
//   3. Aggregate kept barcodes into per-identity pseudo-bulk entries
//   4. Construct contamFinder3 in bulk mode (c fixed at 1.0)
//   5. Solve for pi (mixture proportions)
//   6. Bootstrap pi for Dirichlet concentration parameters
//   7. Write output
//
// Outputs:
//   --interindividual: {prefix}.contam_prof_empty + {prefix}.empty_diagnostics.tsv
//   --interspecies:    {prefix}.species_prof_empty + {prefix}.empty_diagnostics.tsv
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
#include <random>
#include <htswrapper/robin_hood/robin_hood.h>
#include <htswrapper/bc.h>
#include <htswrapper/gzreader.h>
#include "common.h"
#include "ambient_rna_three_ap.h"
#include "io.h"

using std::cout;
using std::endl;
using namespace std;

const string TOOL_VERSION = "2.1-exact-include";
const string TOOL_NAME = "tet_ambient_profile";


struct CountRow {
    unsigned long cell;
    int indv1;
    int type1;
    int indv2;
    int type2;
    float ref;
    float alt;
};

struct PseudoBulkLoadStats {
    size_t rows_read;
    size_t rows_kept;
    size_t empty_barcodes_available;
    size_t cell_barcodes_available;
    size_t empty_barcodes_kept;
    size_t cell_barcodes_kept;
    size_t categories_kept;
    size_t informative_barcodes;
    PseudoBulkLoadStats() : rows_read(0), rows_kept(0),
        empty_barcodes_available(0), cell_barcodes_available(0),
        empty_barcodes_kept(0), cell_barcodes_kept(0), categories_kept(0),
        informative_barcodes(0) {}
};

static bool parse_count_row(const string& line, CountRow& row){
    // Count files are written as:
    // barcode_encoded  indv1  type1  indv2  type2  ref_count  alt_count
    int n = sscanf(line.c_str(), "%lu\t%d\t%d\t%d\t%d\t%f\t%f",
        &row.cell, &row.indv1, &row.type1, &row.indv2, &row.type2,
        &row.ref, &row.alt);
    if (n == 7) return true;
    // Be permissive if a file was space-delimited by an older writer.
    n = sscanf(line.c_str(), "%lu %d %d %d %d %f %f",
        &row.cell, &row.indv1, &row.type1, &row.indv2, &row.type2,
        &row.ref, &row.alt);
    return n == 7;
}

static bool parse_count_barcode_only(const string& line, unsigned long& cell){
    int n = sscanf(line.c_str(), "%lu", &cell);
    return n == 1;
}

enum class NativePanelMode { INTERINDIVIDUAL_NATIVE, INTERSPECIES_NATIVE };

static const char* panel_mode_name(NativePanelMode mode){
    return mode == NativePanelMode::INTERSPECIES_NATIVE
        ? "INTERSPECIES_NATIVE" : "INTERINDIVIDUAL_NATIVE";
}

static int max_index_in_counts_file(const string& filename){
    gzFile gz = gzopen(filename.c_str(), "r");
    if (!gz){
        fprintf(stderr, "ERROR: could not open counts file for dimension guard: %s\n", filename.c_str());
        exit(1);
    }
    char buf[1 << 16];
    int max_idx = -1;
    while (gzgets(gz, buf, sizeof(buf)) != Z_NULL){
        unsigned long cell = 0;
        int i = -1, g1 = -1, j = -1, g2 = -1;
        float ref = 0.0f, alt = 0.0f;
        int n = sscanf(buf, "%lu\t%d\t%d\t%d\t%d\t%f\t%f", &cell, &i, &g1, &j, &g2, &ref, &alt);
        if (n != 7){
            n = sscanf(buf, "%lu %d %d %d %d %f %f", &cell, &i, &g1, &j, &g2, &ref, &alt);
        }
        if (n != 7) continue;
        if (i > max_idx) max_idx = i;
        if (j > max_idx) max_idx = j;
    }
    gzclose(gz);
    return max_idx;
}

static int max_index_in_condf_file(const string& filename){
    ifstream infile(filename.c_str());
    if (!infile.good()){
        fprintf(stderr, "ERROR: could not open condf file for dimension guard: %s\n", filename.c_str());
        exit(1);
    }
    int idx1 = -1, type = -1, idx2 = -1;
    float frac = 0.0f;
    int max_idx = -1;
    while (infile >> idx1 >> type >> idx2 >> frac){
        if (idx1 > max_idx) max_idx = idx1;
        if (idx2 > max_idx) max_idx = idx2;
    }
    return max_idx;
}

static void enforce_native_dimensions(
    NativePanelMode mode,
    const string& counts_name,
    const string& condf_name,
    const vector<string>& samples){

    int n = (int)samples.size();
    int max_counts = max_index_in_counts_file(counts_name);
    int max_condf = max_index_in_condf_file(condf_name);
    fprintf(stderr, "Dimension guard [%s]: n_samples=%d max_counts_index=%d max_condf_index=%d\n",
        panel_mode_name(mode), n, max_counts, max_condf);
    if (max_counts >= n){
        fprintf(stderr, "ERROR: %s has index %d but %s loaded only %d samples.\n",
            counts_name.c_str(), max_counts, panel_mode_name(mode), n);
        fprintf(stderr, "       This indicates cross-resolution input, e.g. an individual-shaped species file.\n");
        exit(1);
    }
    if (max_condf >= n){
        fprintf(stderr, "ERROR: %s has index %d but %s loaded only %d samples.\n",
            condf_name.c_str(), max_condf, panel_mode_name(mode), n);
        fprintf(stderr, "       This indicates cross-resolution input, e.g. an individual-shaped species condf.\n");
        exit(1);
    }
}

static void load_filtered_barcodes_streaming(const string& filename,
    robin_hood::unordered_set<unsigned long>& cell_barcodes){

    gzreader reader(filename);
    while (reader.next()){
        string bc_line = reader.line;
        if (bc_line.empty()) continue;
        while (!bc_line.empty() && (bc_line.back() == '\r' || bc_line.back() == '\n' ||
            bc_line.back() == ' ' || bc_line.back() == '\t')){
            bc_line.pop_back();
        }
        if (bc_line.empty()) continue;
        unsigned long ul = bc_ul(bc_line);
        cell_barcodes.insert(ul);
    }
}

static void load_positive_barcodes_streaming(const string& filename,
    vector<unsigned long>& requested,
    robin_hood::unordered_set<unsigned long>& requested_set){

    gzreader reader(filename);
    size_t line_number = 0;
    while (reader.next()){
        line_number++;
        string bc_line = reader.line;
        if (bc_line.empty()) continue;
        while (!bc_line.empty() && (bc_line.back() == '\r' || bc_line.back() == '\n' ||
            bc_line.back() == ' ' || bc_line.back() == '\t')){
            bc_line.pop_back();
        }
        if (bc_line.empty()) continue;
        unsigned long ul = bc_ul(bc_line);
        if (!requested_set.insert(ul).second){
            fprintf(stderr, "ERROR: duplicate barcode in --include_barcodes at line %lu: %s\n",
                line_number, bc_line.c_str());
            exit(1);
        }
        requested.push_back(ul);
    }
    if (requested.empty()){
        fprintf(stderr, "ERROR: --include_barcodes contains no barcodes: %s\n", filename.c_str());
        exit(1);
    }
}

static string file_signature_fnv1a64(const string& filename){
    std::ifstream in(filename.c_str(), std::ios::binary);
    if (!in){
        fprintf(stderr, "ERROR: cannot open include-list for signature: %s\n", filename.c_str());
        exit(1);
    }
    unsigned long long hash = 1469598103934665603ULL;
    char buffer[8192];
    while (in.good()){
        in.read(buffer, sizeof(buffer));
        std::streamsize n = in.gcount();
        for (std::streamsize i = 0; i < n; ++i){
            hash ^= static_cast<unsigned char>(buffer[i]);
            hash *= 1099511628211ULL;
        }
    }
    char out[64];
    snprintf(out, sizeof(out), "fnv1a64:%016llx", hash);
    return string(out);
}

static void discover_count_barcodes(const string& counts_name,
    const robin_hood::unordered_set<unsigned long>& cell_barcodes,
    vector<unsigned long>& empty_bc_list,
    vector<unsigned long>& cell_bc_list,
    size_t& rows_read){

    robin_hood::unordered_set<unsigned long> seen;
    gzreader reader(counts_name);
    while (reader.next()){
        rows_read++;
        unsigned long cell;
        if (!parse_count_barcode_only(reader.line, cell)) continue;
        if (seen.insert(cell).second){
            if (cell_barcodes.count(cell) > 0){
                cell_bc_list.push_back(cell);
            } else {
                empty_bc_list.push_back(cell);
            }
        }
    }
}

static void aggregate_selected_counts_streaming(const string& counts_name,
    const robin_hood::unordered_set<unsigned long>& kept_barcodes,
    const set<int>& allowed_ids,
    map<pair<int, int>, map<pair<int, int>, pair<float, float> > >& bulk_agg,
    PseudoBulkLoadStats& stats){

    gzreader reader(counts_name);
    CountRow row;
    robin_hood::unordered_set<unsigned long> informative_barcodes;
    while (reader.next()){
        stats.rows_read++;
        if (!parse_count_row(reader.line, row)) continue;
        if (kept_barcodes.count(row.cell) == 0) continue;

        pair<int, int> indv1key = make_pair(row.indv1, row.type1);
        pair<int, int> indv2key = make_pair(row.indv2, row.type2);
        if (!allowed_ids.empty()){
            if (allowed_ids.find(indv1key.first) == allowed_ids.end()) continue;
            if (indv2key.first != -1 && allowed_ids.find(indv2key.first) == allowed_ids.end()) continue;
        }
        bulk_agg[indv1key][indv2key].first += row.ref;
        bulk_agg[indv1key][indv2key].second += row.alt;
        informative_barcodes.insert(row.cell);
        stats.rows_kept++;
    }
    stats.informative_barcodes = informative_barcodes.size();
}

void help(int code){
    fprintf(stderr, "%s v%s\n", TOOL_NAME.c_str(), TOOL_VERSION.c_str());
    fprintf(stderr, "Estimates the ambient RNA mixture profile (pi) from empty droplets.\n");
    fprintf(stderr, "Filters out cell-containing barcodes, optionally subsamples empties,\n");
    fprintf(stderr, "aggregates to pseudo-bulk, and solves for mixture proportions.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "REQUIRED:\n");
    fprintf(stderr, "    --output_prefix, -o   Output prefix. Interindividual mode finds .counts/.condf/.samples;\n");
    fprintf(stderr, "                          interspecies mode finds .species_counts/.species_condf/.species_samples.\n");
    fprintf(stderr, "    --filtered_barcodes   Path to filtered barcode list (one per line,\n");
    fprintf(stderr, "                          cell-containing barcodes to EXCLUDE)\n");
    fprintf(stderr, "    --include_barcodes    Experimental exact positive barcode list; requested\n");
    fprintf(stderr, "                          barcodes must be present in counts and non-cell.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "SNP PANEL (exactly one required, mutually exclusive):\n");
    fprintf(stderr, "    --interindividual     Use .counts and .condf at the output prefix\n");
    fprintf(stderr, "    --interspecies        Native species path: use .species_counts/.species_condf/.species_samples.\n");
    fprintf(stderr, "                          Species inputs must already be species-shaped.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL:\n");
    fprintf(stderr, "    --panel_metadata, -P  Deprecated for native --interspecies solve time;\n");
    fprintf(stderr, "                          species aggregation is producer-side in demux_parallel.\n");
    fprintf(stderr, "    --max_empty, -m       Subsample empty barcodes to at most N\n");
    fprintf(stderr, "                          (default: 50000). Set to 0 for no subsampling.\n");
    fprintf(stderr, "    --cell_fraction       Fraction of cell-containing barcodes to include\n");
    fprintf(stderr, "                          (default: 0.0 = only truly empty droplets).\n");
    fprintf(stderr, "                          1.0 reproduces old behavior of using all barcodes.\n");
    fprintf(stderr, "    --seed                RNG seed for subsampling (default: 42)\n");
    fprintf(stderr, "    --bootstrap, -b       Bootstrap replicates for Dirichlet concentrations\n");
    fprintf(stderr, "                          (default: 100)\n");
    fprintf(stderr, "    --num_threads, -T     Parallel threads (default: all available)\n");
    fprintf(stderr, "    --error_ref, -e       Ref->alt error rate (default: 0.001)\n");
    fprintf(stderr, "    --error_alt, -E       Alt->ref error rate (default: 0.001)\n");
    fprintf(stderr, "    --libname, -n         Library name for diagnostics output\n");
    fprintf(stderr, "    --condf, -F           Pre-computed .condf file\n");
    fprintf(stderr, "                          (if absent, derived from output_prefix)\n");
    fprintf(stderr, "    --ids, -i             Filtered individual list (one name per line)\n");
    fprintf(stderr, "    --output, -O          Override output profile file path\n");
    fprintf(stderr, "    --diagnostics         Override diagnostics TSV path (recommended for exact lists)\n");
    fprintf(stderr, "    --version             Display the tool version and exit\n");
    fprintf(stderr, "    --help, -h            Display this message and exit\n");
    exit(code);
}

int main(int argc, char *argv[]){

    static struct option long_options[] = {
        {"output_prefix",     required_argument, 0, 'o'},
        {"output",            required_argument, 0, 'O'},
        {"num_threads",       required_argument, 0, 'T'},
        {"bootstrap",         required_argument, 0, 'b'},
        {"n_bootstrap",       required_argument, 0, 'b'},
        {"error_ref",         required_argument, 0, 'e'},
        {"error_alt",         required_argument, 0, 'E'},
        {"ids",               required_argument, 0, 'i'},
        {"libname",           required_argument, 0, 'n'},
        {"condf",             required_argument, 0, 'F'},
        {"panel_metadata",    required_argument, 0, 'P'},
        {"filtered_barcodes", required_argument, 0, 1003},
        {"include_barcodes",  required_argument, 0, 1006},
        {"diagnostics",       required_argument, 0, 1007},
        {"max_empty",         required_argument, 0, 'm'},
        {"cell_fraction",     required_argument, 0, 1004},
        {"seed",              required_argument, 0, 1005},
        {"interindividual",   no_argument,       0, 1001},
        {"interspecies",      no_argument,       0, 1002},
        {"version",           no_argument,       0, 1008},
        {"help",              no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    string output_prefix = "";
    string output_file = "";
    string diagnostics_file = "";
    int num_threads = 0;
    int n_bootstrap = 100;
    double error_ref = 0.001;
    double error_alt = 0.001;
    string idfile = "";
    bool idfile_given = false;
    string libname = "";
    string condf_file = "";
    string panel_metadata_file = "";
    string filtered_barcodes_file = "";
    string include_barcodes_file = "";
    int max_empty = 50000;
    double cell_fraction = 0.0;
    int seed = 42;
    bool use_interindividual = false;
    bool use_interspecies = false;

    int option_index = 0;
    int ch;

    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "o:O:T:b:e:E:i:n:F:P:m:h",
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
            case 'P':
                panel_metadata_file = optarg;
                break;
            case 'm':
                max_empty = atoi(optarg);
                break;
            case 1001:
                use_interindividual = true;
                break;
            case 1002:
                use_interspecies = true;
                break;
            case 1003:
                filtered_barcodes_file = optarg;
                break;
            case 1004:
                cell_fraction = atof(optarg);
                break;
            case 1005:
                seed = atoi(optarg);
                break;
            case 1006:
                include_barcodes_file = optarg;
                break;
            case 1007:
                diagnostics_file = optarg;
                break;
            case 1008:
                fprintf(stdout, "%s v%s\n", TOOL_NAME.c_str(), TOOL_VERSION.c_str());
                return 0;
            case 'h':
                help(0);
                break;
            default:
                help(1);
                break;
        }
    }

    // ---- Validate arguments ----

    if (output_prefix.empty()){
        fprintf(stderr, "ERROR: --output_prefix/-o is required\n");
        exit(1);
    }
    if (filtered_barcodes_file.empty()){
        fprintf(stderr, "ERROR: --filtered_barcodes is required\n");
        exit(1);
    }
    if (!use_interindividual && !use_interspecies){
        fprintf(stderr, "ERROR: exactly one of --interindividual or --interspecies is required\n");
        exit(1);
    }
    if (use_interindividual && use_interspecies){
        fprintf(stderr, "ERROR: --interindividual and --interspecies are mutually exclusive\n");
        exit(1);
    }
    if (num_threads <= 1){
        num_threads = 0;
    }

    // ---- Default output file ----

    if (output_file.empty()){
        if (use_interspecies){
            output_file = output_prefix + ".species_prof_empty";
        } else {
            output_file = output_prefix + ".contam_prof_empty";
        }
    }

    fprintf(stderr, "%s v%s\n", TOOL_NAME.c_str(), TOOL_VERSION.c_str());
    fprintf(stderr, "  Output prefix: %s\n", output_prefix.c_str());
    fprintf(stderr, "  Output file: %s\n", output_file.c_str());
    fprintf(stderr, "  SNP panel: %s\n", use_interspecies ? "interspecies" : "interindividual");
    fprintf(stderr, "  Selection mode: %s\n",
        include_barcodes_file.empty() ? "legacy_complement" : "exact_include_list");
    fprintf(stderr, "  Max empty barcodes: %d%s%s\n", max_empty,
        max_empty == 0 ? " (no subsampling)" : "",
        include_barcodes_file.empty() ? "" : " (disabled by exact include-list)");
    fprintf(stderr, "  Cell include fraction: %.3f\n", cell_fraction);
    fprintf(stderr, "  Seed: %d\n", seed);

    // ---- Determine counts/condf file names based on panel selection ----

    string counts_name;
    string condf_name;
    if (use_interspecies){
        counts_name = output_prefix + ".species_counts";
        condf_name = condf_file.empty() ? output_prefix + ".species_condf" : condf_file;
    } else {
        counts_name = output_prefix + ".counts";
        condf_name = condf_file.empty() ? output_prefix + ".condf" : condf_file;
    }

    // ---- Load native sample set ----

    NativePanelMode panel_mode = use_interspecies
        ? NativePanelMode::INTERSPECIES_NATIVE
        : NativePanelMode::INTERINDIVIDUAL_NATIVE;
    string sample_name = output_prefix + (use_interspecies ? ".species_samples" : ".samples");
    vector<string> samples;
    if (file_exists(sample_name)){
        load_samples(sample_name, samples);
    } else {
        fprintf(stderr, "ERROR: no %s samples file found for %s. Expected %s\n",
            use_interspecies ? "species-native" : "individual-native",
            output_prefix.c_str(), sample_name.c_str());
        exit(1);
    }
    fprintf(stderr, "  Panel mode: %s\n", panel_mode_name(panel_mode));
    fprintf(stderr, "  Loaded %lu samples from %s\n", samples.size(), sample_name.c_str());

    if (use_interspecies){
        fprintf(stderr, "  Interspecies mode: using species-native SNP panel\n");
    }

    // ---- Load conditional matching probabilities ----

    map<pair<int, int>, map<int, float> > exp_match_fracs;
    if (file_exists(condf_name)){
        load_exp_fracs(condf_name, exp_match_fracs);
    } else {
        fprintf(stderr, "ERROR: no .condf file found: %s\n", condf_name.c_str());
        exit(1);
    }

    enforce_native_dimensions(panel_mode, counts_name, condf_name, samples);

    // ---- Load filtered ID list, if given ----

    set<int> allowed_ids;
    set<int> allowed_ids2;
    if (idfile_given){
        parse_idfile(idfile, samples, allowed_ids, allowed_ids2, true);
    }

    // ---- Load filtered (cell-containing) barcodes ----

    robin_hood::unordered_set<unsigned long> cell_barcodes;
    load_filtered_barcodes_streaming(filtered_barcodes_file, cell_barcodes);
    fprintf(stderr, "  Loaded %lu filtered (cell) barcodes\n", cell_barcodes.size());

    // ---- Streaming two-pass pseudo-bulk aggregation ----
    // Pass 1 discovers unique barcode IDs only.  Pass 2 parses full count rows
    // only for the selected barcode set and aggregates directly into one bulk
    // count map.  This avoids materializing hundreds of thousands of raw
    // per-barcode count maps before filtering.

    if (!file_exists(counts_name)){
        fprintf(stderr, "ERROR: no counts file found: %s\n", counts_name.c_str());
        exit(1);
    }

    vector<unsigned long> empty_bc_list;
    vector<unsigned long> cell_bc_list;
    size_t pass1_rows = 0;
    fprintf(stderr, "Streaming pass 1 over %s: discovering count barcodes...\n",
        counts_name.c_str());
    discover_count_barcodes(counts_name, cell_barcodes, empty_bc_list,
        cell_bc_list, pass1_rows);

    size_t n_empty_available = empty_bc_list.size();
    size_t n_cell_available = cell_bc_list.size();
    fprintf(stderr, "  Pass 1 rows scanned: %lu\n", pass1_rows);
    fprintf(stderr, "  Barcodes in counts: %lu empty, %lu cell\n",
        n_empty_available, n_cell_available);

    size_t include_requested = 0;
    size_t include_found = 0;
    string include_signature = "not_applicable";
    if (!include_barcodes_file.empty()){
        if (cell_fraction != 0.0){
            fprintf(stderr, "ERROR: --cell_fraction must be 0 when --include_barcodes is supplied\n");
            exit(1);
        }
        vector<unsigned long> requested;
        robin_hood::unordered_set<unsigned long> requested_set;
        load_positive_barcodes_streaming(include_barcodes_file, requested, requested_set);
        include_requested = requested.size();
        include_signature = file_signature_fnv1a64(include_barcodes_file);
        robin_hood::unordered_set<unsigned long> empty_available_set;
        robin_hood::unordered_set<unsigned long> cell_available_set;
        for (auto bc : empty_bc_list) empty_available_set.insert(bc);
        for (auto bc : cell_bc_list) cell_available_set.insert(bc);
        vector<unsigned long> exact;
        vector<unsigned long> missing;
        vector<unsigned long> filtered_overlap;
        for (auto bc : requested){
            if (cell_barcodes.count(bc) > 0 || cell_available_set.count(bc) > 0){
                filtered_overlap.push_back(bc);
            } else if (empty_available_set.count(bc) == 0){
                missing.push_back(bc);
            } else {
                exact.push_back(bc);
            }
        }
        if (!filtered_overlap.empty()){
            fprintf(stderr, "ERROR: --include_barcodes contains %lu filtered-cell barcode(s); first encoded barcode=%lu\n",
                filtered_overlap.size(), filtered_overlap.front());
            exit(1);
        }
        if (!missing.empty()){
            fprintf(stderr, "ERROR: %lu requested include-list barcode(s) are absent from counts; first encoded barcode=%lu\n",
                missing.size(), missing.front());
            exit(1);
        }
        empty_bc_list.swap(exact);
        cell_bc_list.clear();
        include_found = empty_bc_list.size();
        fprintf(stderr, "  Exact include-list: requested=%lu found=%lu signature=%s\n",
            include_requested, include_found, include_signature.c_str());
    }

    // ---- Subsample empties (legacy complement mode only) ----

    if (include_barcodes_file.empty() && max_empty > 0 && (int)empty_bc_list.size() > max_empty){
        std::mt19937 rng(seed);
        std::shuffle(empty_bc_list.begin(), empty_bc_list.end(), rng);
        empty_bc_list.resize(max_empty);
        fprintf(stderr, "  Subsampled empties: %d (from %lu)\n", max_empty, n_empty_available);
    }

    // ---- Subsample cells ----

    size_t n_cell_kept = 0;
    if (cell_fraction > 0.0 && !cell_bc_list.empty()){
        size_t n_keep = (size_t)(cell_fraction * cell_bc_list.size());
        if (n_keep < 1) n_keep = 1;
        if (n_keep < cell_bc_list.size()){
            std::mt19937 rng2(seed + 1);
            std::shuffle(cell_bc_list.begin(), cell_bc_list.end(), rng2);
            cell_bc_list.resize(n_keep);
        }
        n_cell_kept = cell_bc_list.size();
        fprintf(stderr, "  Including %lu cell barcodes (fraction %.3f)\n",
            n_cell_kept, cell_fraction);
    } else {
        cell_bc_list.clear();
        n_cell_kept = 0;
    }

    robin_hood::unordered_set<unsigned long> kept_barcodes;
    for (auto bc : empty_bc_list) kept_barcodes.insert(bc);
    for (auto bc : cell_bc_list) kept_barcodes.insert(bc);

    map<pair<int, int>, map<pair<int, int>, pair<float, float> > > bulk_agg;
    PseudoBulkLoadStats pb_stats;
    pb_stats.empty_barcodes_available = n_empty_available;
    pb_stats.cell_barcodes_available = n_cell_available;
    pb_stats.empty_barcodes_kept = empty_bc_list.size();
    pb_stats.cell_barcodes_kept = n_cell_kept;

    fprintf(stderr, "Streaming pass 2: aggregating %lu selected barcodes...\n",
        kept_barcodes.size());
    aggregate_selected_counts_streaming(counts_name, kept_barcodes, allowed_ids,
        bulk_agg, pb_stats);

    size_t n_categories = 0;
    for (auto& t1 : bulk_agg) n_categories += t1.second.size();
    pb_stats.categories_kept = n_categories;

    if (bulk_agg.empty()){
        fprintf(stderr, "ERROR: pseudo-bulk aggregation produced no count categories.\n");
        fprintf(stderr, "  Counts: %s\n", counts_name.c_str());
        fprintf(stderr, "  Selected barcodes: %lu\n", kept_barcodes.size());
        exit(1);
    }

    robin_hood::unordered_map<unsigned long, map<pair<int, int>,
        map<pair<int, int>, pair<float, float> > > > indv_allelecounts;
    robin_hood::unordered_map<unsigned long, int> assn;
    robin_hood::unordered_map<unsigned long, double> assn_llr;
    unsigned long synth_bc = 1;
    indv_allelecounts[synth_bc] = std::move(bulk_agg);
    assn[synth_bc] = 0;         // placeholder identity (unused in bulk compiler)
    assn_llr[synth_bc] = 100.0; // placeholder LLR

    fprintf(stderr, "Pseudo-bulk aggregation complete:\n");
    fprintf(stderr, "  Rows scanned in pass 2: %lu\n", pb_stats.rows_read);
    fprintf(stderr, "  Rows retained: %lu\n", pb_stats.rows_kept);
    fprintf(stderr, "  Barcodes retained: %lu empty + %lu cell = %lu\n",
        pb_stats.empty_barcodes_kept, pb_stats.cell_barcodes_kept,
        pb_stats.empty_barcodes_kept + pb_stats.cell_barcodes_kept);
    fprintf(stderr, "  Categories retained: %lu\n", pb_stats.categories_kept);

    // Ensure the allowed ID set covers all VCF samples so idx2samp is fully
    // populated. Without this, an empty allowed_ids set from the constructor
    // would build idx2samp from assignments only (which is just the placeholder).
    if (allowed_ids.empty()){
        for (int i = 0; i < (int)samples.size(); ++i){
            allowed_ids.insert(i);
        }
    }

    // ---- Create contamFinder3 in bulk mode ----

    fprintf(stderr, "Setting up bulk-mode contamFinder3...\n");
    contamFinder3 cf(indv_allelecounts, assn, assn_llr, exp_match_fracs,
        samples.size(), allowed_ids, allowed_ids, allowed_ids2);
    cf.set_error_rates(error_ref, error_alt);
    cf.set_num_threads(num_threads);
    cf.set_bulk_mode(true);
    cf.no_reassign();

    // ---- Solve for pi ----

    fprintf(stderr, "Estimating ambient profile...\n");
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

    if (use_interspecies){
        // Native species output: samples are already species labels, so the
        // ordinary contamination profile is the species profile.
        FILE* outf = fopen(output_file.c_str(), "w");
        if (!outf){
            fprintf(stderr, "ERROR: cannot open output file: %s\n", output_file.c_str());
            exit(1);
        }
        dump_contam_prof(outf, cf.contam_prof, contam_prof_conc, samples);
        fclose(outf);
        fprintf(stderr, "Native species profile written to: %s\n", output_file.c_str());

        for (int i = 0; i < (int)samples.size(); i++){
            if (cf.contam_prof.count(i) > 0){
                fprintf(stderr, "  %s: %.6f", samples[i].c_str(), cf.contam_prof[i]);
                if (contam_prof_conc.count(i) > 0){
                    fprintf(stderr, " (alpha=%.2f)", contam_prof_conc[i]);
                }
                fprintf(stderr, "\n");
            }
        }
    } else {
        // Individual-level output: .contam_prof_empty with individual names
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
        string diag_name = diagnostics_file.empty()
            ? output_prefix + ".empty_diagnostics.tsv"
            : diagnostics_file;
        FILE* diagf = fopen(diag_name.c_str(), "w");
        if (diagf){
            fprintf(diagf, "metric\tvalue\n");
            fprintf(diagf, "tool_name\t%s\n", TOOL_NAME.c_str());
            fprintf(diagf, "tool_version\t%s\n", TOOL_VERSION.c_str());
            fprintf(diagf, "selection_mode\t%s\n",
                include_barcodes_file.empty() ? "legacy_complement" : "exact_include_list");
            fprintf(diagf, "include_list_path\t%s\n",
                include_barcodes_file.empty() ? "not_applicable" : include_barcodes_file.c_str());
            fprintf(diagf, "include_list_signature\t%s\n", include_signature.c_str());
            fprintf(diagf, "barcodes_requested\t%lu\n",
                include_barcodes_file.empty() ? empty_bc_list.size() : include_requested);
            fprintf(diagf, "barcodes_found_in_counts\t%lu\n",
                include_barcodes_file.empty() ? empty_bc_list.size() : include_found);
            fprintf(diagf, "barcodes_with_informative_individual_observations\t%lu\n",
                use_interindividual ? pb_stats.informative_barcodes : 0UL);
            fprintf(diagf, "barcodes_with_informative_species_observations\t%lu\n",
                use_interspecies ? pb_stats.informative_barcodes : 0UL);
            fprintf(diagf, "count_rows_used\t%lu\n", pb_stats.rows_kept);
            fprintf(diagf, "labels_supported\t%d\n", (int)samples.size());
            fprintf(diagf, "profile_identifiability_status\t%s\n",
                (pb_stats.rows_kept > 0 && pb_stats.categories_kept > 0) ? "identifiable_input" : "non_identifiable");
            fprintf(diagf, "n_empty_available\t%lu\n", n_empty_available);
            fprintf(diagf, "n_empty_kept\t%lu\n", empty_bc_list.size());
            fprintf(diagf, "n_cell_available\t%lu\n", n_cell_available);
            fprintf(diagf, "n_cell_kept\t%lu\n", n_cell_kept);
            fprintf(diagf, "cell_include_fraction\t%.2f\n", cell_fraction);
            fprintf(diagf, "max_empty_barcodes\t%d\n", max_empty);
            fprintf(diagf, "seed\t%d\n", seed);
            fprintf(diagf, "snp_panel\t%s\n",
                use_interspecies ? "interspecies" : "interindividual");
            fprintf(diagf, "n_samples\t%d\n", (int)samples.size());
            fprintf(diagf, "n_bootstrap\t%d\n", n_bootstrap);
            fprintf(diagf, "error_ref\t%.6f\n", error_ref);
            fprintf(diagf, "error_alt\t%.6f\n", error_alt);
            if (!libname.empty()){
                fprintf(diagf, "libname\t%s\n", libname.c_str());
            }

            // Per-component proportions
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
