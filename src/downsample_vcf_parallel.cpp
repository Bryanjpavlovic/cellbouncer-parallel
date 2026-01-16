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
#include <unistd.h>
#include <map>
#include <unordered_map>
#include <set>
#include <bitset>
#include <cstdlib>
#include <random>
#include <functional>
#include <utility>
#include <cmath>
#include <float.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include <htslib/kseq.h>
#include <zlib.h>
#include <htswrapper/bc.h>
#include <htswrapper/bam.h>
#include <htswrapper/gzreader.h>
#include <mixtureDist/functions.h>
#include <optimML/multivar_ml.h>
#include <optimML/brent.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include <iomanip>
#include "common.h"
#include "downsample_vcf.h"
#include <chrono>
#include <omp.h>
#include <mutex>
#include <atomic>

using std::cout;
using std::endl;
using namespace std;

// ============================================================================
// BJPs ADDITIONS: GTF and Coverage support structures
// ============================================================================

struct Interval {
    int start;
    int end;
    bool operator<(const Interval& other) const {
        return start < other.start;
    }
};

void load_gtf(const string& gtf_file, map<string, vector<Interval>>& annotations) {
    ifstream file(gtf_file);
    string line;
    while (getline(file, line)) {
        if (line[0] == '#') continue;
        stringstream ss(line);
        string seqname, source, feature, start_s, end_s;
        ss >> seqname >> source >> feature >> start_s >> end_s;
        
        if (feature == "gene" || feature == "exon") {
            try {
                Interval iv = {stoi(start_s), stoi(end_s)};
                annotations[seqname].push_back(iv);
            } catch (...) { continue; }
        }
    }
    for (auto& pair : annotations) {
        sort(pair.second.begin(), pair.second.end());
    }
    fprintf(stderr, "Loaded annotations for %lu chromosomes.\n", annotations.size());
}

float get_annot_score(const string& chrom, int pos, const map<string, vector<Interval>>& annotations) {
    if (annotations.count(chrom) == 0) return 0.1;
    
    const vector<Interval>& ivs = annotations.at(chrom);
    auto it = lower_bound(ivs.begin(), ivs.end(), Interval{pos, pos});
    
    if (it != ivs.end() && it->start <= pos && it->end >= pos) return 1.0;
    if (it != ivs.begin()) {
        auto prev = std::prev(it);
        if (prev->start <= pos && prev->end >= pos) return 1.0;
    }
    
    return 0.1;
}

// Coverage interval with score
struct CoverageInterval {
    int start;
    int end;
    float score;
    bool operator<(const CoverageInterval& other) const {
        return start < other.start;
    }
};

// Load coverage bedgraph into memory
map<string, vector<CoverageInterval>> load_coverage(const string& filename) {
    map<string, vector<CoverageInterval>> cov_map;
    
    htsFile* fp = hts_open(filename.c_str(), "r");
    if (!fp) return cov_map;
    
    kstring_t str = {0, 0, 0};
    long count = 0;
    
    while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
        if (str.l == 0 || str.s[0] == '#') continue;
        
        char* chrom = str.s;
        char* p = str.s;
        int col = 0;
        int start = 0, end = 0;
        float score = 0;
        
        while (*p) {
            if (*p == '\t') {
                *p = '\0';
                col++;
                if (col == 1) start = atoi(p + 1);
                else if (col == 2) end = atoi(p + 1);
                else if (col == 3) {
                    score = atof(p + 1);
                    break;
                }
            }
            p++;
        }
        
        if (col >= 3) {
            cov_map[string(chrom)].push_back({start, end, score});
            count++;
            if (count % 10000000 == 0) {
                fprintf(stderr, "  Loaded %ld coverage intervals...\r", count);
            }
        }
    }
    
    free(str.s);
    hts_close(fp);
    
    // Sort each chromosome's intervals
    for (auto& kv : cov_map) {
        sort(kv.second.begin(), kv.second.end());
    }
    
    fprintf(stderr, "  Loaded %ld coverage intervals from %zu chromosomes\n", count, cov_map.size());
    return cov_map;
}

// Fast in-memory coverage lookup
float get_coverage_score_fast(const string& chrom, int pos, const map<string, vector<CoverageInterval>>& coverage) {
    auto it = coverage.find(chrom);
    if (it == coverage.end()) return 0.0;
    
    const vector<CoverageInterval>& ivs = it->second;
    auto cit = lower_bound(ivs.begin(), ivs.end(), CoverageInterval{pos, 0, 0});
    
    // Check current interval
    if (cit != ivs.end() && cit->start <= pos && cit->end > pos) {
        return cit->score;
    }
    // Check previous interval
    if (cit != ivs.begin()) {
        auto prev = std::prev(cit);
        if (prev->start <= pos && prev->end > pos) {
            return prev->score;
        }
    }
    
    return 0.0;
}

float get_coverage_score(const string& chrom, int pos, tbx_t* tbx, htsFile* fp) {
    if (!tbx || !fp) return 0.0;
    
    int tid = tbx_name2id(tbx, chrom.c_str());
    if (tid < 0) return 0.0;

    hts_itr_t* itr = tbx_itr_queryi(tbx, tid, pos, pos + 1);
    if (!itr) return 0.0;

    kstring_t str = {0, 0, 0};
    float cov = 0.0;
    
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
        char* p = str.s;
        int col = 0;
        while (*p) {
            if (col == 3) {
                cov = atof(p);
                break;
            }
            if (*p == '\t') col++;
            p++;
        }
    }
    
    tbx_itr_destroy(itr);
    free(str.s);
    return cov;
}

// ============================================================================
// HET SITE CANDIDATE STRUCTURE
// ============================================================================

struct HetCandidate {
    string chrom;
    int64_t pos;
    int chrom_idx;          // Index into chroms vector for ordering
    int bin_idx;            // Bin index within chromosome
    float het_score;        // Combined het score
    float het_freq;         // Fraction of called individuals that are het
    float missing_freq;     // Fraction of individuals with missing genotype
    float annot_score;      // Annotation score (if GTF provided)
    float cov_score;        // Coverage score (if coverage provided)
    float demux_score;      // Clade rarity score (computed after Pass 1)
    string ref_allele;
    string alt_allele;
    vector<int32_t> genotypes;  // Store genotypes for output
    int n_samples;
    
    // For sorting by het_score descending
    bool operator<(const HetCandidate& other) const {
        return het_score > other.het_score;  // Higher score = better
    }
};

// ============================================================================
// STORED VARIANT STRUCTURE - for in-memory VCF processing
// ============================================================================

struct StoredVariant {
    int chrom_idx;              // bcf_record->rid
    int64_t pos;                // bcf_record->pos
    int n_alleles;              // bcf_record->n_allele
    string alleles;             // All alleles comma-separated: "A,T" or "A,T,C"
    vector<int32_t> genotypes;  // Raw genotypes (num_samples * 2)
    
    // Precomputed by proc_bcf_record logic (only valid for biallelic)
    bool is_biallelic;          // n_alleles == 2
    bool passes_qc;             // What proc_bcf_record would return
    bitset<NBITS> alt;          // From proc_bcf_record
    bitset<NBITS> alt_flip;     // From proc_bcf_record  
    bitset<NBITS> present;      // From proc_bcf_record
    
    // For stratified selection (Option C)
    size_t variant_idx;         // Index into all_variants for back-reference
};

// ============================================================================
// PARALLEL CODE
// ============================================================================

bool operator<(const std::bitset<NBITS>& a, const std::bitset<NBITS>& b) {
    for (int i = NBITS-1; i >= 0; --i) {
        if (a[i] != b[i]) return b[i];
    }
    return false;
}

void help(int code){
    fprintf(stderr, "downsample_vcf_parallel [OPTIONS]\n");
    fprintf(stderr, "Parallel version of downsample_vcf - uses OpenMP for chromosome-level parallelization.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Given a VCF file containing variants segregating within a panel of\n");
    fprintf(stderr, "deeply divergent genomes, intelligently downsamples the variant panel.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --vcf -v     A VCF/BCF file listing variants (must be indexed).\n");
    fprintf(stderr, "    --output -o  Output VCF file name. Will be gzipped.\n");
    fprintf(stderr, "    --num -n     Desired final number of SNPs (not required with --annotate_only).\n");
    fprintf(stderr, "===== OPTIONAL =====\n");
    fprintf(stderr, "    --threads -t Number of threads (default: all available).\n");
    fprintf(stderr, "    --gtf -g     GTF annotation file. Adds ANNOT_SCORE INFO field.\n");
    fprintf(stderr, "    --cov -c     Tabix-indexed bedgraph coverage file. Adds COV_SCORE INFO field.\n");
    fprintf(stderr, "    --min_cov -m Minimum coverage threshold (default: 1.0). Sites below this are excluded.\n");
    fprintf(stderr, "    --seed -s    Random seed for reproducibility.\n");
    fprintf(stderr, "===== HET SITE SELECTION =====\n");
    fprintf(stderr, "    --het_output -H  Output file for het-enriched SNP set (optional).\n");
    fprintf(stderr, "    --het_num -N     Target number of het sites (default: half of --num).\n");
    fprintf(stderr, "    --bin_size -b    Bin size (bp) for evenness scoring (default: 100000).\n");
    fprintf(stderr, "===== ANNOTATION MODES =====\n");
    fprintf(stderr, "    --annotate_only  Annotate ALL sites with scores and exit (no downsampling).\n");
    fprintf(stderr, "    --annotate_all FILE  Annotate ALL sites first, save to FILE, then downsample.\n");
    fprintf(stderr, "                         Produces both full annotated VCF and downsampled outputs.\n");
    fprintf(stderr, "    INFO fields added: HET_FREQ, MISSING_FREQ, HET_SCORE, N_HET, N_CALLED,\n");
    fprintf(stderr, "                       ANNOT_SCORE (if -g), COV_SCORE (if -c), BIN_ID\n");
    fprintf(stderr, "    --help -h    Display this message and exit.\n");
    exit(code);
}

bool proc_bcf_record(bcf1_t* bcf_record,
    bcf_hdr_t* bcf_header,
    int num_samples,
    bitset<NBITS>& alt,
    bitset<NBITS>& alt_flip,
    bitset<NBITS>& present,
    bool& pass){

    alt.reset();
    alt_flip.reset();
    present.reset();

    bcf_unpack(bcf_record, BCF_UN_STR);
    pass = true;
    
    for (int i = 0; i < bcf_record->n_allele; ++i){
        if (strcmp(bcf_record->d.allele[i], "A") != 0 &&
            strcmp(bcf_record->d.allele[i], "C") != 0 &&
            strcmp(bcf_record->d.allele[i], "G") != 0 && 
            strcmp(bcf_record->d.allele[i], "T") != 0){
            pass = false;
            break;
        }
    }
    if (bcf_record->d.allele[0][0] == bcf_record->d.allele[1][0]){
        pass = false;
    }
    if (pass){
        int32_t* gts = NULL;
        int n_gts = 0;
        int num_loaded = bcf_get_genotypes(bcf_header, bcf_record, &gts, &n_gts);
        if (num_loaded <= 0){
            free(gts);
            return false;
        }
        
        int ploidy = 2;
        int nref = 0;
        int nalt = 0;
         
        for (int i = 0; i < num_samples; ++i){
            int32_t* gtptr = gts + i*ploidy;
            
            if (!bcf_gt_is_missing(gtptr[0])){
                present.set(i);
                alt_flip.set(i);
                
                if (bcf_gt_allele(gtptr[0]) != 0 || bcf_gt_allele(gtptr[1]) != 0){
                    alt.set(i);
                    alt_flip.reset(i);
                    nalt++;
                } 
                else{
                    nref++;
                }
            }    
        }
        
        free(gts);
     
        if (nref == 0 || nalt == 0){
            pass = false;
            return false;
        } 
        return true;
    }
    return false;
}

// Process BCF record and compute het statistics
// Returns true if site has het individuals, false otherwise
bool proc_bcf_record_het(bcf1_t* bcf_record,
    bcf_hdr_t* bcf_header,
    int num_samples,
    int& n_het,          // Output: count of het individuals
    int& n_called,       // Output: count of called (non-missing) individuals
    int& n_missing,      // Output: count of missing individuals
    vector<int32_t>& gt_copy, // Output: copy of genotypes for later use
    bool& pass){

    bcf_unpack(bcf_record, BCF_UN_STR);
    pass = true;
    n_het = 0;
    n_called = 0;
    n_missing = 0;
    gt_copy.clear();
    
    // Check alleles are valid ACGT
    for (int i = 0; i < bcf_record->n_allele; ++i){
        if (strcmp(bcf_record->d.allele[i], "A") != 0 &&
            strcmp(bcf_record->d.allele[i], "C") != 0 &&
            strcmp(bcf_record->d.allele[i], "G") != 0 && 
            strcmp(bcf_record->d.allele[i], "T") != 0){
            pass = false;
            break;
        }
    }
    if (bcf_record->d.allele[0][0] == bcf_record->d.allele[1][0]){
        pass = false;
    }
    
    if (!pass) return false;
    
    int32_t* gts = NULL;
    int n_gts = 0;
    int num_loaded = bcf_get_genotypes(bcf_header, bcf_record, &gts, &n_gts);
    if (num_loaded <= 0){
        free(gts);
        return false;
    }
    
    int ploidy = 2;
    
    // Copy genotypes and count het
    gt_copy.resize(num_samples * ploidy);
    for (int i = 0; i < num_samples * ploidy; ++i) {
        gt_copy[i] = gts[i];
    }
    
    for (int i = 0; i < num_samples; ++i){
        int32_t* gtptr = gts + i*ploidy;
        
        if (bcf_gt_is_missing(gtptr[0]) || bcf_gt_is_missing(gtptr[1])){
            n_missing++;
        } else {
            n_called++;
            int a1 = bcf_gt_allele(gtptr[0]);
            int a2 = bcf_gt_allele(gtptr[1]);
            // Het = one ref and one alt (0/1 or 1/0)
            if ((a1 == 0 && a2 != 0) || (a1 != 0 && a2 == 0)){
                n_het++;
            }
        }
    }
    
    free(gts);
    return (n_het > 0);  // Only return true if at least one het individual
}

void count_branch(unordered_map<bitset<NBITS>, double>& branchcounts,
    bitset<NBITS>& clade,
    bitset<NBITS>& clade_flip,
    int cladecount,
    int cladecount_flip,
    double count){
   
    if (cladecount < cladecount_flip || (cladecount == cladecount_flip && clade < clade_flip)){
        if (branchcounts.count(clade) == 0){
            branchcounts.insert(make_pair(clade, count));
        }
        else{
            branchcounts[clade] += count;
        }
    }
    else{
        if (branchcounts.count(clade_flip) == 0){
            branchcounts.insert(make_pair(clade_flip, count));
        }
        else{
            branchcounts[clade_flip] += count;
        }
    }
}

void count_branch_missing(unordered_map<pair<bitset<NBITS>, bitset<NBITS> >, double>& branchcounts,
    unordered_map<bitset<NBITS>, bitset<NBITS> >& miss2flip,
    bitset<NBITS>& present,
    bitset<NBITS>& clade,
    bitset<NBITS>& clade_flip,
    int cladecount,
    int cladecount_flip,
    double count){
    
    pair<bitset<NBITS>, bitset<NBITS> > key;
    
    if (cladecount < cladecount_flip || (cladecount == cladecount_flip && clade < clade_flip)){
        key = make_pair(present, clade);
        if (miss2flip.count(clade) == 0){
            miss2flip.insert(make_pair(clade, clade_flip));
        }
    }
    else{
        key = make_pair(present, clade_flip);
        if (miss2flip.count(clade_flip) == 0){
            miss2flip.insert(make_pair(clade_flip, clade));
        }
    }
    if (branchcounts.count(key) == 0){
        branchcounts.insert(make_pair(key, count));
    }
    else{
        branchcounts[key] += count;
    }
}

// Merge thread-local branchcounts into global
void merge_branchcounts(
    unordered_map<bitset<NBITS>, double>& global,
    unordered_map<bitset<NBITS>, double>& local,
    mutex& mtx){
    
    lock_guard<mutex> lock(mtx);
    for (auto& kv : local){
        if (global.count(kv.first) == 0){
            global.insert(kv);
        }
        else{
            global[kv.first] += kv.second;
        }
    }
}

void merge_branchcounts_missing(
    unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, double>& global,
    unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, double>& local,
    mutex& mtx){
    
    lock_guard<mutex> lock(mtx);
    for (auto& kv : local){
        if (global.count(kv.first) == 0){
            global.insert(kv);
        }
        else{
            global[kv.first] += kv.second;
        }
    }
}

void merge_miss2flip(
    unordered_map<bitset<NBITS>, bitset<NBITS>>& global,
    unordered_map<bitset<NBITS>, bitset<NBITS>>& local,
    mutex& mtx){
    
    lock_guard<mutex> lock(mtx);
    for (auto& kv : local){
        if (global.count(kv.first) == 0){
            global.insert(kv);
        }
    }
}

// Merge thread-local het candidates into global
void merge_het_candidates(
    vector<HetCandidate>& global,
    vector<HetCandidate>& local,
    mutex& mtx){
    
    lock_guard<mutex> lock(mtx);
    global.insert(global.end(), 
        make_move_iterator(local.begin()), 
        make_move_iterator(local.end()));
}

double brentfun(double param, const map<string, double>& data_d, const map<string, int>& data_i){
    int num = data_i.at("num");
    char buf[50];
    string bufstr;
    double sum = 0.0;
    for (int i = 0; i < num; ++i){
        sprintf(&buf[0], "x%d", i);
        bufstr = buf;
        double count = data_d.at(bufstr);
        sum += pow(count, param);
    }
    return sum - (double)data_i.at("target");
}

double dbrentfun(double param, const map<string, double>& data_d, const map<string, int>& data_i){
    int num = data_i.at("num");
    char buf[50];
    string bufstr;
    double d_df = 0.0;
    for (int i = 0; i < num; ++i){
        sprintf(&buf[0], "x%d", i);
        bufstr = buf;
        double count = data_d.at(bufstr);
        d_df += pow(count, param) * log(count);
    }
    return d_df;
}

// Helper to get bin index from position
inline int get_bin_idx(int64_t pos, int bin_size) {
    return (int)(pos / bin_size);
}

// Helper to create BIN_ID string
inline string make_bin_id(const string& chrom, int bin_idx) {
    return chrom + ":" + to_string(bin_idx);
}

// Compute within-clade selection score for Option C stratified selection
// Coverage dominates, annotation is a small tiebreaker bonus
inline float compute_selection_score(float cov_score, float annot_score) {
    // annot_score is 1.0 for genic, 0.1 for intergenic
    // Convert to annot_boost: 1.0 if genic, 0.0 if intergenic
    float annot_boost = (annot_score >= 1.0f) ? 1.0f : 0.0f;
    return log2f(cov_score + 1.0f) + (0.1f * annot_boost);
}

// Structure for tracking SNPs within a clade for stratified selection
struct CladeCandidate {
    size_t variant_idx;      // Index into all_variants
    float selection_score;   // log2(cov+1) + 0.1*annot_boost
    
    bool operator<(const CladeCandidate& other) const {
        return selection_score > other.selection_score;  // Descending order
    }
};


int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"vcf", required_argument, 0, 'v'},
       {"output", required_argument, 0, 'o'},
       {"num", required_argument, 0, 'n'},
       {"threads", required_argument, 0, 't'},
       {"gtf", required_argument, 0, 'g'},
       {"cov", required_argument, 0, 'c'},
       {"seed", required_argument, 0, 's'},
       {"het_output", required_argument, 0, 'H'},
       {"het_num", required_argument, 0, 'N'},
       {"bin_size", required_argument, 0, 'b'},
       {"min_cov", required_argument, 0, 'm'},
       {"annotate_only", no_argument, 0, 'A'},
       {"annotate_all", required_argument, 0, 'a'},
       {0, 0, 0, 0} 
    };
    
    string vcf_file = "";
    string outfile = "";
    string gtf_file = "";
    string cov_file = "";
    string het_outfile = "";
    string annotate_all_file = "";  // Output file for full annotated VCF
    int num = -1;
    int het_num = -1;  // Will default to num/2 if not specified
    int bin_size = 100000;  // 100kb default
    float min_cov = 1.0f;  // Minimum coverage threshold (hard filter)
    int n_threads = omp_get_max_threads();
    int seed = -1;
    bool annotate_only = false;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:n:o:t:g:c:s:H:N:b:m:Aa:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                break;
            case 'h':
                help(0);
                break;
            case 'v':
                vcf_file = optarg;
                break;
            case 'n':
                num = atoi(optarg);
                break;
            case 'o':
                outfile = optarg;
                break;
            case 't':
                n_threads = atoi(optarg);
                break;
            case 'g':
                gtf_file = optarg;
                break;
            case 'c':
                cov_file = optarg;
                break;
            case 's':
                seed = atoi(optarg);
                break;
            case 'H':
                het_outfile = optarg;
                break;
            case 'N':
                het_num = atoi(optarg);
                break;
            case 'b':
                bin_size = atoi(optarg);
                break;
            case 'm':
                min_cov = atof(optarg);
                break;
            case 'A':
                annotate_only = true;
                break;
            case 'a':
                annotate_all_file = optarg;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    if (vcf_file == ""){
        fprintf(stderr, "ERROR: VCF file required\n");
        exit(1);
    }
    if (num <= 0 && !annotate_only){
        fprintf(stderr, "ERROR: --num/-n must be positive (Recommend > 1M)\n");
        exit(1);
    }
    if (outfile == ""){
        fprintf(stderr, "ERROR: --output/-o is required.\n");
        exit(1);
    }
    if (outfile.rfind(".bcf") == string::npos && outfile.rfind(".vcf") == string::npos){
        outfile += ".bcf";
    }
    
    // Handle annotate_all output file extension
    if (!annotate_all_file.empty()) {
        if (annotate_all_file.rfind(".bcf") == string::npos && annotate_all_file.rfind(".vcf") == string::npos){
            annotate_all_file += ".bcf";
        }
    }
    
    // Set het_num default to half of num if not specified
    if (het_num < 0) {
        het_num = num / 2;
    }
    
    // Handle het output file extension
    if (!het_outfile.empty()) {
        if (het_outfile.rfind(".bcf") == string::npos && het_outfile.rfind(".vcf") == string::npos){
            het_outfile += ".het.bcf";
        }
    }

    fprintf(stderr, "Using %d threads\n", n_threads);
    fprintf(stderr, "Bin size for evenness scoring: %d bp\n", bin_size);
    fprintf(stderr, "Minimum coverage threshold: %.1f\n", min_cov);
    if (annotate_only) {
        fprintf(stderr, "Mode: ANNOTATE ONLY (all sites will be annotated, no downsampling)\n");
    } else if (!annotate_all_file.empty()) {
        fprintf(stderr, "Mode: ANNOTATE ALL first -> %s, then downsample\n", annotate_all_file.c_str());
    }
    if (!het_outfile.empty() && !annotate_only) {
        fprintf(stderr, "Het output enabled: %s (target: %d sites)\n", het_outfile.c_str(), het_num);
    }
    omp_set_num_threads(n_threads);

    // ========================================================================
    // Load optional annotation and coverage files
    // ========================================================================
    map<string, vector<Interval>> annotations;
    if (!gtf_file.empty()) {
        load_gtf(gtf_file, annotations);
    }

    // Load coverage into memory for fast lookup
    map<string, vector<CoverageInterval>> coverage_map;
    if (!cov_file.empty()) {
        fprintf(stderr, "Loading coverage from %s...\n", cov_file.c_str());
        coverage_map = load_coverage(cov_file);
    }

    // ========================================================================
    // ANNOTATE_ONLY MODE: Stream through VCF, add scores to all sites, exit
    // ========================================================================
    if (annotate_only) {
        fprintf(stderr, "\n========================================\n");
        fprintf(stderr, "ANNOTATE ONLY MODE\n");
        fprintf(stderr, "========================================\n");
        
        auto start_time = chrono::high_resolution_clock::now();
        
        // Open input VCF
        htsFile* in_fp = bcf_open(vcf_file.c_str(), "r");
        if (!in_fp) {
            fprintf(stderr, "ERROR: Could not open %s\n", vcf_file.c_str());
            exit(1);
        }
        hts_set_threads(in_fp, n_threads);
        
        bcf_hdr_t* hdr = bcf_hdr_read(in_fp);
        if (!hdr) {
            fprintf(stderr, "ERROR: Could not read VCF header\n");
            exit(1);
        }
        
        int n_samples = bcf_hdr_nsamples(hdr);
        fprintf(stderr, "Samples: %d\n", n_samples);
        
        // Add INFO headers
        bcf_hdr_append(hdr, "##INFO=<ID=HET_FREQ,Number=1,Type=Float,Description=\"Fraction of called individuals that are heterozygous\">");
        bcf_hdr_append(hdr, "##INFO=<ID=MISSING_FREQ,Number=1,Type=Float,Description=\"Fraction of individuals with missing genotype\">");
        bcf_hdr_append(hdr, "##INFO=<ID=HET_SCORE,Number=1,Type=Float,Description=\"Combined het score: het_freq * (1 - missing_freq * 0.5)\">");
        bcf_hdr_append(hdr, "##INFO=<ID=N_HET,Number=1,Type=Integer,Description=\"Number of heterozygous individuals\">");
        bcf_hdr_append(hdr, "##INFO=<ID=N_CALLED,Number=1,Type=Integer,Description=\"Number of individuals with called genotype\">");
        if (!annotations.empty()) {
            bcf_hdr_append(hdr, "##INFO=<ID=ANNOT_SCORE,Number=1,Type=Float,Description=\"Annotation weight (1.0=genic, 0.1=intergenic)\">");
        }
        if (!coverage_map.empty()) {
            bcf_hdr_append(hdr, "##INFO=<ID=COV_SCORE,Number=1,Type=Float,Description=\"Coverage score from bedgraph\">");
        }
        bcf_hdr_append(hdr, "##INFO=<ID=BIN_ID,Number=1,Type=String,Description=\"Chromosome:bin_index for evenness tracking\">");
        
        if (bcf_hdr_sync(hdr) < 0) {
            fprintf(stderr, "ERROR: Failed to sync header\n");
            exit(1);
        }
        
        // Open output BCF
        htsFile* out_fp = hts_open(outfile.c_str(), "wb");
        if (!out_fp) {
            fprintf(stderr, "ERROR: Could not open output file %s\n", outfile.c_str());
            exit(1);
        }
        hts_set_threads(out_fp, n_threads);
        
        if (bcf_hdr_write(out_fp, hdr) < 0) {
            fprintf(stderr, "ERROR: Failed to write header\n");
            exit(1);
        }
        
        // Process records
        bcf1_t* rec = bcf_init();
        int32_t* gts = NULL;
        int n_gts = 0;
        
        long n_records = 0;
        long n_with_het = 0;
        
        fprintf(stderr, "Processing variants...\n");
        
        while (bcf_read(in_fp, hdr, rec) >= 0) {
            
            // Get genotypes
            int ret = bcf_get_genotypes(hdr, rec, &gts, &n_gts);
            
            int n_het = 0;
            int n_called = 0;
            int n_missing = 0;
            
            if (ret > 0) {
                int ploidy = ret / n_samples;
                
                for (int i = 0; i < n_samples; i++) {
                    int32_t* gt_ptr = gts + i * ploidy;
                    
                    // Check for missing
                    if (bcf_gt_is_missing(gt_ptr[0]) || 
                        (ploidy > 1 && bcf_gt_is_missing(gt_ptr[1]))) {
                        n_missing++;
                        continue;
                    }
                    
                    n_called++;
                    
                    // Check for het (different alleles)
                    if (ploidy >= 2) {
                        int a1 = bcf_gt_allele(gt_ptr[0]);
                        int a2 = bcf_gt_allele(gt_ptr[1]);
                        if (a1 != a2) {
                            n_het++;
                        }
                    }
                }
            } else {
                // No genotype data - all missing
                n_missing = n_samples;
            }
            
            // Compute metrics
            float het_freq = (n_called > 0) ? (float)n_het / n_called : 0.0f;
            float missing_freq = (float)n_missing / n_samples;
            float missing_penalty = 1.0f - (missing_freq * 0.5f);
            float het_score = het_freq * missing_penalty;
            
            // Add INFO fields
            bcf_update_info_float(hdr, rec, "HET_FREQ", &het_freq, 1);
            bcf_update_info_float(hdr, rec, "MISSING_FREQ", &missing_freq, 1);
            bcf_update_info_float(hdr, rec, "HET_SCORE", &het_score, 1);
            bcf_update_info_int32(hdr, rec, "N_HET", &n_het, 1);
            bcf_update_info_int32(hdr, rec, "N_CALLED", &n_called, 1);
            
            // Add ANNOT_SCORE if GTF provided
            if (!annotations.empty()) {
                const char* chrom = bcf_hdr_id2name(hdr, rec->rid);
                float annot_score = get_annot_score(chrom, rec->pos, annotations);
                bcf_update_info_float(hdr, rec, "ANNOT_SCORE", &annot_score, 1);
            }
            
            // Add COV_SCORE if coverage provided
            if (!coverage_map.empty()) {
                const char* chrom = bcf_hdr_id2name(hdr, rec->rid);
                float cov_score = get_coverage_score_fast(chrom, rec->pos, coverage_map);
                bcf_update_info_float(hdr, rec, "COV_SCORE", &cov_score, 1);
            }
            
            // Add BIN_ID
            const char* chrom = bcf_hdr_id2name(hdr, rec->rid);
            int bin_idx_val = get_bin_idx(rec->pos, bin_size);
            string bin_id = make_bin_id(chrom, bin_idx_val);
            bcf_update_info_string(hdr, rec, "BIN_ID", bin_id.c_str());
            
            // Write record
            if (bcf_write(out_fp, hdr, rec) < 0) {
                fprintf(stderr, "ERROR: Failed to write record\n");
                exit(1);
            }
            
            n_records++;
            if (n_het > 0) n_with_het++;
            
            if (n_records % 1000000 == 0) {
                fprintf(stderr, "  Processed %ld million variants...\r", n_records / 1000000);
                fflush(stderr);
            }
        }
        
        fprintf(stderr, "\n");
        
        // Cleanup
        free(gts);
        bcf_destroy(rec);
        bcf_hdr_destroy(hdr);
        hts_close(in_fp);
        hts_close(out_fp);
        
        // Index output BCF
        fprintf(stderr, "Indexing output BCF...\n");
        if (bcf_index_build(outfile.c_str(), 14) < 0) {
            fprintf(stderr, "WARNING: Failed to create index for %s\n", outfile.c_str());
        }
        
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
        
        fprintf(stderr, "\n========================================\n");
        fprintf(stderr, "ANNOTATION COMPLETE\n");
        fprintf(stderr, "========================================\n");
        fprintf(stderr, "  Total variants: %ld\n", n_records);
        fprintf(stderr, "  With het individuals: %ld (%.1f%%)\n", 
                n_with_het, 100.0 * n_with_het / n_records);
        fprintf(stderr, "  Time: %ld seconds\n", duration.count());
        fprintf(stderr, "  Output: %s\n", outfile.c_str());
        
        return 0;
    }

    // Note: --annotate_all is now handled in Pass 2 (parallel, with all scores including DEMUX_SCORE)

    // ========================================================================
    // Get chromosome list from VCF
    // ========================================================================
    htsFile* bcf_reader = bcf_open(vcf_file.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR interpreting %s as BCF format.\n", vcf_file.c_str());
        exit(1);
    }
    bcf_hdr_t* bcf_header = bcf_hdr_read(bcf_reader);
    int num_samples = bcf_hdr_nsamples(bcf_header);
    
    vector<string> samples;
    for (int i = 0; i < num_samples; ++i){
        samples.push_back(bcf_header->samples[i]);
    }
    
    if (samples.size() * 2 > NBITS){
        fprintf(stderr, "ERROR: too many samples in VCF. Please recompile with NBITS = %ld or higher.\n",
            samples.size() * 2);
        exit(1); 
    }

    // Get chromosome names
    vector<string> chroms;
    map<string, int> chrom_to_idx;  // For ordering chromosomes
    int n_seqs = 0;
    const char** seqnames = bcf_hdr_seqnames(bcf_header, &n_seqs);
    for (int i = 0; i < n_seqs; ++i){
        chrom_to_idx[seqnames[i]] = chroms.size();
        chroms.push_back(seqnames[i]);
    }
    free(seqnames);
    
    bcf_hdr_destroy(bcf_header);
    hts_close(bcf_reader);

    fprintf(stderr, "Found %d chromosomes, %d samples\n", (int)chroms.size(), num_samples);

    // ========================================================================
    // Load VCF index (CSI required for parallel region queries)
    // ========================================================================
    
    fprintf(stderr, "Loading VCF index...\n");
    hts_idx_t* vcf_idx = bcf_index_load(vcf_file.c_str());
    if (!vcf_idx) {
        fprintf(stderr, "ERROR: Could not load VCF index (.csi required for parallel queries)\n");
        fprintf(stderr, "  Create with: bcftools index -c %s\n", vcf_file.c_str());
        exit(1);
    }
    fprintf(stderr, "  Index loaded successfully\n");

    // ========================================================================
    // LOAD ENTIRE VCF INTO MEMORY
    // ========================================================================
    
    fprintf(stderr, "\n========================================\n");
    fprintf(stderr, "Loading VCF into memory...\n");
    fprintf(stderr, "========================================\n");
    
    auto load_start = chrono::high_resolution_clock::now();
    
    htsFile* load_fp = bcf_open(vcf_file.c_str(), "r");
    if (!load_fp) {
        fprintf(stderr, "ERROR: Could not open %s for loading\n", vcf_file.c_str());
        exit(1);
    }
    hts_set_threads(load_fp, n_threads);
    
    bcf_hdr_t* load_hdr = bcf_hdr_read(load_fp);
    bcf1_t* load_rec = bcf_init();
    int32_t* load_gts = NULL;
    int load_n_gts = 0;
    
    vector<StoredVariant> all_variants;
    all_variants.reserve(150000000);  // Pre-allocate for ~150M variants
    
    long n_loaded = 0;
    long n_biallelic = 0;
    long n_pass_qc = 0;
    
    while (bcf_read(load_fp, load_hdr, load_rec) >= 0) {
        StoredVariant sv;
        sv.chrom_idx = load_rec->rid;
        sv.pos = load_rec->pos;
        sv.n_alleles = load_rec->n_allele;
        sv.is_biallelic = (load_rec->n_allele == 2);
        sv.passes_qc = false;
        
        // Get alleles - store ALL alleles for ALL variants
        bcf_unpack(load_rec, BCF_UN_STR);
        if (load_rec->n_allele > 0) {
            sv.alleles = load_rec->d.allele[0];
            for (int i = 1; i < load_rec->n_allele; i++) {
                sv.alleles += ",";
                sv.alleles += load_rec->d.allele[i];
            }
        }
        if (sv.is_biallelic) {
            n_biallelic++;
        }
        
        // Get genotypes
        int num_loaded_gt = bcf_get_genotypes(load_hdr, load_rec, &load_gts, &load_n_gts);
        if (num_loaded_gt > 0) {
            // Store raw genotypes
            sv.genotypes.assign(load_gts, load_gts + num_loaded_gt);
            
            // For biallelic variants, precompute what proc_bcf_record would compute
            // EXACTLY matching proc_bcf_record logic
            if (sv.is_biallelic) {
                bool pass = true;
                
                // Check alleles are valid A/C/G/T (exactly as proc_bcf_record does)
                for (int i = 0; i < load_rec->n_allele; ++i) {
                    if (strcmp(load_rec->d.allele[i], "A") != 0 &&
                        strcmp(load_rec->d.allele[i], "C") != 0 &&
                        strcmp(load_rec->d.allele[i], "G") != 0 && 
                        strcmp(load_rec->d.allele[i], "T") != 0) {
                        pass = false;
                        break;
                    }
                }
                // Check ref != alt (exactly as proc_bcf_record does)
                if (load_rec->d.allele[0][0] == load_rec->d.allele[1][0]) {
                    pass = false;
                }
                
                if (pass) {
                    // Compute bitsets exactly as proc_bcf_record does
                    sv.alt.reset();
                    sv.alt_flip.reset();
                    sv.present.reset();
                    
                    int ploidy = 2;  // Hardcoded as in proc_bcf_record
                    int nref = 0;
                    int nalt = 0;
                    
                    for (int i = 0; i < num_samples; ++i) {
                        int32_t* gtptr = load_gts + i * ploidy;
                        
                        // proc_bcf_record only checks gtptr[0] for missing
                        if (!bcf_gt_is_missing(gtptr[0])) {
                            sv.present.set(i);
                            sv.alt_flip.set(i);
                            
                            // proc_bcf_record: alt if EITHER allele is non-zero
                            if (bcf_gt_allele(gtptr[0]) != 0 || bcf_gt_allele(gtptr[1]) != 0) {
                                sv.alt.set(i);
                                sv.alt_flip.reset(i);
                                nalt++;
                            } else {
                                nref++;
                            }
                        }
                    }
                    
                    // proc_bcf_record returns true only if nref > 0 AND nalt > 0
                    if (nref > 0 && nalt > 0) {
                        sv.passes_qc = true;
                        n_pass_qc++;
                    }
                }
            }
        }
        
        sv.variant_idx = n_loaded;  // Store index for back-reference
        
        all_variants.push_back(std::move(sv));
        n_loaded++;
        
        if (n_loaded % 5000000 == 0) {
            fprintf(stderr, "  Loaded %ld variants (%ld biallelic, %ld pass QC)...\r", 
                n_loaded, n_biallelic, n_pass_qc);
            fflush(stderr);
        }
    }
    
    free(load_gts);
    bcf_destroy(load_rec);
    bcf_hdr_destroy(load_hdr);
    hts_close(load_fp);
    
    all_variants.shrink_to_fit();
    
    auto load_end = chrono::high_resolution_clock::now();
    auto load_duration = chrono::duration_cast<chrono::seconds>(load_end - load_start);
    
    fprintf(stderr, "\n  Loaded %ld variants in %ld seconds\n", n_loaded, load_duration.count());
    fprintf(stderr, "  Biallelic: %ld (%.1f%%)\n", n_biallelic, 100.0 * n_biallelic / n_loaded);
    fprintf(stderr, "  Pass QC: %ld (%.1f%%)\n", n_pass_qc, 100.0 * n_pass_qc / n_loaded);
    
    // Build chromosome -> variant range index
    vector<pair<size_t, size_t>> chrom_ranges(chroms.size(), {0, 0});
    {
        size_t start_idx = 0;
        int current_chrom = -1;
        for (size_t i = 0; i < all_variants.size(); i++) {
            if (all_variants[i].chrom_idx != current_chrom) {
                if (current_chrom >= 0 && current_chrom < (int)chroms.size()) {
                    chrom_ranges[current_chrom] = {start_idx, i};
                }
                current_chrom = all_variants[i].chrom_idx;
                start_idx = i;
            }
        }
        if (current_chrom >= 0 && current_chrom < (int)chroms.size()) {
            chrom_ranges[current_chrom] = {start_idx, all_variants.size()};
        }
    }

    // ========================================================================
    // PASS 1: Count clades AND collect het candidates (PARALLEL BY CHROMOSOME)
    // Now using in-memory variants instead of disk I/O
    // ========================================================================
    
    unordered_map<bitset<NBITS>, double> branchcounts;
    unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, double> branchcounts_missing;
    unordered_map<bitset<NBITS>, bitset<NBITS>> miss2flip;
    vector<HetCandidate> all_het_candidates;  // Global het candidates
    
    mutex bc_mutex, bcm_mutex, m2f_mutex, het_mutex;
    atomic<long> total_snps(0);
    atomic<long> total_het_candidates(0);
    atomic<int> chroms_done(0);

    auto start_time = chrono::high_resolution_clock::now();
    
    bool collect_het = !het_outfile.empty();
    fprintf(stderr, "Pass 1: Counting mutations on branches%s (parallel, %d threads)...\n", 
        collect_het ? " and collecting het candidates" : "", n_threads);

    #pragma omp parallel num_threads(n_threads)
    {
        // Thread-local data structures (SAME as before)
        unordered_map<bitset<NBITS>, double> local_branchcounts;
        unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, double> local_branchcounts_missing;
        unordered_map<bitset<NBITS>, bitset<NBITS>> local_miss2flip;
        vector<HetCandidate> local_het_candidates;
        
        #pragma omp for schedule(dynamic, 1)
        for (size_t c = 0; c < chroms.size(); c++){
            const string& chrom = chroms[c];
            
            // Get range of variants for this chromosome from memory
            size_t range_start = chrom_ranges[c].first;
            size_t range_end = chrom_ranges[c].second;
            
            if (range_start >= range_end) {
                chroms_done++;
                continue;
            }
            
            long local_snps = 0;
            long local_het = 0;
            
            // Iterate over stored variants for this chromosome
            for (size_t vi = range_start; vi < range_end; vi++) {
                const StoredVariant& sv = all_variants[vi];
                
                // Skip non-biallelic (same as original: if n_allele != 2 continue)
                if (!sv.is_biallelic) continue;
                
                // Use precomputed proc_bcf_record result
                if (sv.passes_qc) {
                    // Make copies since count_branch/count_branch_missing take non-const refs
                    bitset<NBITS> alt = sv.alt;
                    bitset<NBITS> alt_flip = sv.alt_flip;
                    bitset<NBITS> present = sv.present;
                    
                    int ac = alt.count();
                    int afc = alt_flip.count();
                
                    if (present.count() == (size_t)num_samples){
                        count_branch(local_branchcounts, alt, alt_flip, ac, afc, 1.0);
                    }
                    else{
                        count_branch_missing(local_branchcounts_missing, local_miss2flip, 
                            present, alt, alt_flip, ac, afc, 1.0);
                    }
                    local_snps++;
                }
                
                // Collect het candidates if enabled
                // Using EXACT same logic as proc_bcf_record_het
                if (collect_het && !sv.genotypes.empty() && sv.is_biallelic) {
                    // Parse ref and alt from alleles string
                    size_t comma_pos = sv.alleles.find(',');
                    string ref_allele = sv.alleles.substr(0, comma_pos);
                    string alt_allele = sv.alleles.substr(comma_pos + 1);
                    
                    // Check alleles are valid (same check as proc_bcf_record_het)
                    bool het_pass = true;
                    if (ref_allele.length() != 1 || alt_allele.length() != 1) {
                        het_pass = false;
                    } else {
                        char r = ref_allele[0];
                        char a = alt_allele[0];
                        if ((r != 'A' && r != 'C' && r != 'G' && r != 'T') ||
                            (a != 'A' && a != 'C' && a != 'G' && a != 'T') ||
                            r == a) {
                            het_pass = false;
                        }
                    }
                    
                    if (het_pass) {
                        int ploidy = 2;  // Hardcoded as in proc_bcf_record_het
                        int n_het = 0;
                        int n_called = 0;
                        int n_missing = 0;
                        
                        // EXACT same loop as proc_bcf_record_het
                        for (int i = 0; i < num_samples; ++i) {
                            const int32_t* gtptr = sv.genotypes.data() + i * ploidy;
                            
                            // proc_bcf_record_het checks BOTH alleles for missing
                            if (bcf_gt_is_missing(gtptr[0]) || bcf_gt_is_missing(gtptr[1])) {
                                n_missing++;
                            } else {
                                n_called++;
                                int a1 = bcf_gt_allele(gtptr[0]);
                                int a2 = bcf_gt_allele(gtptr[1]);
                                // Het = one ref and one alt (0/1 or 1/0)
                                if ((a1 == 0 && a2 != 0) || (a1 != 0 && a2 == 0)) {
                                    n_het++;
                                }
                            }
                        }
                        
                        // Only proceed if at least one het (same as proc_bcf_record_het return)
                        if (n_het > 0) {
                            float het_freq = (n_called > 0) ? (float)n_het / n_called : 0.0f;
                            float missing_freq = (float)n_missing / num_samples;
                            
                            // Get annotation and coverage scores FRESH (same as original)
                            float annot_score = 0.1f;
                            if (!annotations.empty()) {
                                annot_score = get_annot_score(chrom, sv.pos, annotations);
                            }
                            float cov_score = 0.0f;
                            if (!coverage_map.empty()) {
                                cov_score = get_coverage_score_fast(chrom, sv.pos, coverage_map);
                            }
                            
                            // Apply min_cov hard filter for het sites
                            if (cov_score < min_cov) {
                                continue;  // Skip sites below coverage threshold
                            }
                            
                            // New het scoring formula: het_freq dominates, then coverage, annotation is bonus
                            // het_score = het_freq × (log2(cov+1) + 0.1 × annot_boost)
                            float annot_boost = (annot_score >= 1.0f) ? 1.0f : 0.0f;
                            float selection_component = log2f(cov_score + 1.0f) + (0.1f * annot_boost);
                            float het_score = het_freq * selection_component;
                            
                            HetCandidate hc;
                            hc.chrom = chrom;
                            hc.pos = sv.pos;
                            hc.chrom_idx = c;
                            hc.bin_idx = get_bin_idx(sv.pos, bin_size);
                            hc.het_score = het_score;
                            hc.het_freq = het_freq;
                            hc.missing_freq = missing_freq;
                            hc.annot_score = annot_score;
                            hc.cov_score = cov_score;  // Store RAW coverage, same as demux
                            hc.ref_allele = ref_allele;
                            hc.alt_allele = alt_allele;
                            hc.genotypes = sv.genotypes;  // Copy genotypes
                            hc.n_samples = num_samples;
                            
                            local_het_candidates.push_back(std::move(hc));
                            local_het++;
                        }
                    }
                }
            }
            
            total_snps += local_snps;
            total_het_candidates += local_het;
            int done = ++chroms_done;
            
            if (local_snps > 0 || done % 100 == 0) {
                fprintf(stderr, "  Chromosome %s: %ld SNPs, %ld het [%d/%zu]\n", 
                    chrom.c_str(), local_snps, local_het, done, chroms.size());
            }
        }
        
        // Merge thread-local counts into global (SAME as original)
        merge_branchcounts(branchcounts, local_branchcounts, bc_mutex);
        merge_branchcounts_missing(branchcounts_missing, local_branchcounts_missing, bcm_mutex);
        merge_miss2flip(miss2flip, local_miss2flip, m2f_mutex);
        if (collect_het) {
            merge_het_candidates(all_het_candidates, local_het_candidates, het_mutex);
        }
    }
    
    // vcf_idx still available for any code that needs it
    
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    
    fprintf(stderr, "\nPass 1 complete: %ld SNPs in %ld seconds\n", total_snps.load(), duration.count());
    fprintf(stderr, "Found %lu unique clade patterns\n", branchcounts.size());
    if (collect_het) {
        fprintf(stderr, "Found %ld het candidates\n", total_het_candidates.load());
    }
    
    if (branchcounts.size() >= (size_t)num){
        fprintf(stderr, "ERROR: %zu distinct allele sharing patterns found, but %d SNPs requested.\n", 
            branchcounts.size(), num);
        fprintf(stderr, "  Please set sample size to at least %zu.\n", branchcounts.size());
        exit(1);
    }

    // ========================================================================
    // HET SITE SELECTION WITH BIN BALANCING
    // (Output writing happens after DEMUX_SCORE computation)
    // ========================================================================
    
    set<pair<string, int64_t>> selected_het_sites;  // (chrom, pos) of selected het sites
    vector<size_t> selected_het_indices;  // Indices into all_het_candidates
    
    if (collect_het && !all_het_candidates.empty()) {
        fprintf(stderr, "\n========================================\n");
        fprintf(stderr, "Selecting het sites with bin balancing...\n");
        start_time = chrono::high_resolution_clock::now();
        
        // Sort candidates by het_score descending
        fprintf(stderr, "  Sorting %zu het candidates by score...\n", all_het_candidates.size());
        sort(all_het_candidates.begin(), all_het_candidates.end());
        
        // Bin tracking: map from (chrom, bin_idx) -> count of selected SNPs
        map<pair<string, int>, int> bins_selected;
        
        // Greedy selection with bin balancing
        selected_het_indices.reserve(het_num);
        
        fprintf(stderr, "  Selecting up to %d het sites...\n", het_num);
        
        for (size_t i = 0; i < all_het_candidates.size() && (int)selected_het_indices.size() < het_num; ++i) {
            const HetCandidate& hc = all_het_candidates[i];
            pair<string, int> bin_key = make_pair(hc.chrom, hc.bin_idx);
            
            // Compute evenness score
            int snps_in_bin = bins_selected[bin_key];
            float evenness_score = 1.0f / (1.0f + snps_in_bin);
            
            // Adjusted score (for potential replacement logic, not used currently)
            float adjusted_score = hc.het_score * evenness_score;
            (void)adjusted_score;  // Suppress unused warning
            
            // Select this candidate
            selected_het_indices.push_back(i);
            bins_selected[bin_key]++;
            selected_het_sites.insert(make_pair(hc.chrom, hc.pos));
            
            if ((selected_het_indices.size() % 100000) == 0) {
                fprintf(stderr, "    Selected %zu het sites...\r", selected_het_indices.size());
                fflush(stderr);
            }
        }
        
        end_time = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
        fprintf(stderr, "  Selected %zu het sites from %zu bins in %ld seconds\n", 
            selected_het_indices.size(), bins_selected.size(), duration.count());
        fprintf(stderr, "  (Het output will be written after DEMUX_SCORE computation)\n");
    }

    // ========================================================================
    // Process missing data and compute downsampling probabilities
    // (Same as original - this is fast)
    // ========================================================================
    
    map<int, vector<bitset<NBITS>>> size2clades;
    for (auto& bc : branchcounts){
        for (int i = 1; i <= (int)bc.first.count(); ++i){
            if (size2clades.count(i) == 0){
                vector<bitset<NBITS>> v;
                size2clades.insert(make_pair(i, v));
            }
            size2clades[i].push_back(bc.first);
        }
    }
    
    fprintf(stderr, "Adjusting counts using data from sites with missing genotypes (%zu entries)...\n", branchcounts_missing.size());
    fflush(stderr);
    
    unordered_map<bitset<NBITS>, vector<bitset<NBITS>>> keyuniqparent;
    unordered_map<bitset<NBITS>, float> to_add;
    unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, unordered_map<bitset<NBITS>, float>> miss_clade_probs;
    
    long adj_count = 0;
    for (auto& bcm : branchcounts_missing){
        adj_count++;
        if (adj_count % 100000 == 0) {
            fprintf(stderr, "  Adjusting: %ld/%zu\r", adj_count, branchcounts_missing.size());
            fflush(stderr);
        }
        unordered_map<bitset<NBITS>, float> mx;
        miss_clade_probs.insert(make_pair(bcm.first, mx));
        double parent_tot = 0.0;

        int k = bcm.first.second.count();
        if (keyuniqparent.count(bcm.first.second) > 0){
            for (auto& v : keyuniqparent[bcm.first.second]){
                if ((v & bcm.first.first).count() == k){
                    double cladecount = branchcounts[v];
                    parent_tot += cladecount;
                    miss_clade_probs[bcm.first].insert(make_pair(v, (float)cladecount));
                }
            }
        }
        else{
            for (auto& v : size2clades[k]){
                if ((v & bcm.first.second).count() == k){
                    if (keyuniqparent.count(bcm.first.second) == 0){
                        vector<bitset<NBITS>> vec;
                        keyuniqparent.insert(make_pair(bcm.first.second, vec));
                    }
                    keyuniqparent[bcm.first.second].push_back(v);
                    if ((v & bcm.first.first).count() == k){
                        double cladecount = branchcounts[v];
                        parent_tot += cladecount;
                        miss_clade_probs[bcm.first].insert(make_pair(v, (float)cladecount));
                    }
                }
            }
        }

        if (miss2flip.count(bcm.first.second) > 0){
            bitset<NBITS> flipped = (miss2flip[bcm.first.second] & bcm.first.first);
            k = flipped.count();
            if (keyuniqparent.count(flipped) > 0){
                for (auto& v : keyuniqparent[flipped]){
                    if ((v & bcm.first.first).count() == k){
                        double cladecount = branchcounts[v];
                        parent_tot += cladecount;
                        miss_clade_probs[bcm.first].insert(make_pair(v, (float)cladecount));
                    }
                }
            }
            else{
                for (auto& v : size2clades[k]){
                    if ((v & flipped).count() == k){
                        if (keyuniqparent.count(flipped) == 0){
                            vector<bitset<NBITS>> vec;
                            keyuniqparent.insert(make_pair(flipped, vec));
                        }
                        keyuniqparent[flipped].push_back(v);
                        if ((v & bcm.first.first).count() == k){
                            double cladecount = branchcounts[v];
                            parent_tot += cladecount;
                            miss_clade_probs[bcm.first].insert(make_pair(v, (float)cladecount));
                        }
                    }
                }
            }
        }
        
        if (parent_tot > 0){
            for (auto& mcp2 : miss_clade_probs[bcm.first]){
                float frac = mcp2.second / parent_tot;
                if (to_add.count(mcp2.first) == 0){
                    to_add.insert(make_pair(mcp2.first, frac * bcm.second));
                }
                else{
                    to_add[mcp2.first] += frac * bcm.second;
                }
                miss_clade_probs[bcm.first][mcp2.first] = frac;
            }
        }
        else{
            miss_clade_probs.erase(bcm.first);
        }
    }
    
    for (auto& ta : to_add){
        branchcounts[ta.first] += ta.second;
    }
    to_add.clear();
    
    // ========================================================================
    // Compute slot allocations using Brent's method (Option C: Stratified Deterministic)
    // ========================================================================
    
    // Map from clade bitset -> number of slots allocated
    unordered_map<bitset<NBITS>, long> clade_slot_allocation;
    // Keep downsample_prob for DEMUX_SCORE computation (still useful for reporting)
    unordered_map<bitset<NBITS>, double> downsample_prob;
    
    // Store the Brent exponent for later use
    double brent_exponent = 1.0;
    
    fprintf(stderr, "Determine number of sites to downsample...\n");
    fflush(stderr);

    double allbc2 = 0.0;
    vector<bitset<NBITS>> bs_idx;
    vector<pair<double, int>> clcountsort;  // (negative_count, index)
    
    // Build sorted list
    fprintf(stderr, "Building clade count list from %zu entries...\n", branchcounts.size());
    fflush(stderr);
    for (auto bc = branchcounts.begin(); bc != branchcounts.end(); ){
        allbc2 += bc->second;
        clcountsort.push_back(make_pair(-bc->second, bs_idx.size()));
        bs_idx.push_back(bc->first);
        branchcounts.erase(bc++);
    }
    fprintf(stderr, "  Built list with %zu entries\n", clcountsort.size());
    fflush(stderr);
    
    fprintf(stderr, "%f total clade counts\n", allbc2);
    fflush(stderr);
    
    if (allbc2 > (double)num){
        fprintf(stderr, "Sorting %zu clade counts...\n", clcountsort.size());
        fflush(stderr);
        sort(clcountsort.begin(), clcountsort.end());
        fprintf(stderr, "Sort complete.\n");
        fflush(stderr);
        
        // Extract counts into a simple vector for fast access
        vector<double> counts(clcountsort.size());
        for (size_t i = 0; i < clcountsort.size(); i++) {
            counts[i] = -clcountsort[i].first;
        }
        
        // Pre-compute cumulative sums
        fprintf(stderr, "Computing cumulative sums...\n");
        fflush(stderr);
        vector<double> cumsum(counts.size() + 1, 0.0);
        for (size_t i = 0; i < counts.size(); i++) {
            cumsum[i + 1] = cumsum[i] + counts[i];
        }
        
        // Find starting point where sum_lo < num
        double threshold = allbc2 - (double)num;
        int start_n_hi = 1;
        for (size_t i = 1; i <= counts.size(); i++) {
            if (cumsum[i] > threshold) {
                start_n_hi = i;
                break;
            }
        }
        fprintf(stderr, "Starting optimization at n_hi=%d (sum_hi=%.0f, need > %.0f)\n", 
            start_n_hi, cumsum[start_n_hi], threshold);
        fprintf(stderr, "Top 10 clade counts: ");
        for (int i = 0; i < min(10, (int)counts.size()); i++) {
            fprintf(stderr, "%.0f ", counts[i]);
        }
        fprintf(stderr, "\n");
        fflush(stderr);
        
        auto brent_start = chrono::high_resolution_clock::now();
        
        // Lambda to compute sum of counts[0..n-1]^x (parallelized)
        auto power_sum = [&counts](int n, double x) -> double {
            double sum = 0.0;
            #pragma omp parallel for reduction(+:sum)
            for (int i = 0; i < n; i++) {
                sum += pow(counts[i], x);
            }
            return sum;
        };
        
        // For each n_hi, find x such that power_sum(n_hi, x) = target
        // Using bisection (simpler and robust)
        int solutions_tried = 0;
        int final_n_hi = -1;
        for (int n_hi = start_n_hi; n_hi <= (int)counts.size(); ++n_hi){
            double sum_hi = cumsum[n_hi];
            double sum_lo = allbc2 - sum_hi;
            
            if (n_hi == (int)counts.size() || sum_lo < (double)num){
                solutions_tried++;
                
                if (solutions_tried % 1000 == 0) {
                    auto now = chrono::high_resolution_clock::now();
                    auto elapsed = chrono::duration_cast<chrono::seconds>(now - brent_start);
                    fprintf(stderr, "  Tried %d solutions, n_hi=%d, elapsed=%lds\r", 
                        solutions_tried, n_hi, elapsed.count());
                    fflush(stderr);
                }
                
                double target = (double)num - sum_lo;
                double lastcount = counts[n_hi - 1];
                double nextcount = (n_hi < (int)counts.size()) ? counts[n_hi] : -1;
                
                // Bisection to find x where power_sum(n_hi, x) = target
                // At x=0: power_sum = n_hi (each count^0 = 1)
                // At x=1: power_sum = sum_hi
                // We need target, which is between n_hi and sum_hi
                
                double x_lo = 0.0, x_hi = 1.0;
                double res = 0.5;
                
                // Quick check: if target is out of range, skip
                if (target < n_hi || target > sum_hi) {
                    continue;
                }
                
                // Bisection with 20 iterations (gives ~6 decimal places)
                for (int iter = 0; iter < 20; iter++) {
                    res = (x_lo + x_hi) / 2.0;
                    double ps = power_sum(n_hi, res);
                    if (ps < target) {
                        x_lo = res;
                    } else {
                        x_hi = res;
                    }
                }
                
                if (res < 1.0){
                    double transformed_last = pow(lastcount, res);
                    
                    if (transformed_last > nextcount){
                        auto brent_end = chrono::high_resolution_clock::now();
                        auto brent_dur = chrono::duration_cast<chrono::seconds>(brent_end - brent_start);
                        fprintf(stderr, "\nFound solution: x = %f at n_hi=%d (after %ld seconds, %d solutions tried)\n", 
                            res, n_hi, brent_dur.count(), solutions_tried);
                        fflush(stderr);
                        
                        brent_exponent = res;
                        final_n_hi = n_hi;
                        
                        // Build slot allocation map (Option C: integer slot counts)
                        // Also keep probability for DEMUX_SCORE reporting
                        fprintf(stderr, "Building slot allocation map for %d clades...\n", n_hi);
                        fflush(stderr);
                        
                        long total_slots_allocated = 0;
                        for (int i = 0; i < n_hi; ++i){
                            double slots_float = pow(counts[i], res);
                            long slots = (long)floor(slots_float);  // Floor to get integer slots
                            if (slots < 1) slots = 1;  // Minimum 1 slot per clade
                            
                            clade_slot_allocation.insert(make_pair(bs_idx[clcountsort[i].second], slots));
                            total_slots_allocated += slots;
                            
                            // Also compute probability for DEMUX_SCORE
                            double dsp = pow(counts[i], res) / counts[i];
                            downsample_prob.insert(make_pair(bs_idx[clcountsort[i].second], dsp));
                        }
                        fprintf(stderr, "  Allocated %ld initial slots across %zu clades\n", 
                            total_slots_allocated, clade_slot_allocation.size());
                        fflush(stderr);
                        break;
                    }
                    else if (solutions_tried <= 10) {
                        fprintf(stderr, "  n_hi=%d: x=%.4f, transformed_last=%.2f <= nextcount=%.2f (failed)\n",
                            n_hi, res, transformed_last, nextcount);
                        fflush(stderr);
                    }
                }
                else if (solutions_tried <= 10) {
                    fprintf(stderr, "  n_hi=%d: x=%.4f >= 1.0 (failed)\n", n_hi, res);
                    fflush(stderr);
                }
            }
        }
    }
    else{
        fprintf(stderr, "%f total clade counts; no downsampling needed\n", allbc2);
        return 0;
    }
    
    // Compute downsample probabilities for missing genotype sites
    fprintf(stderr, "Computing downsample probabilities for missing genotype sites (%zu entries)...\n", miss_clade_probs.size());
    fflush(stderr);
    
    unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, double> downsample_miss_prob;

    long miss_count = 0;
    for (auto mcp = miss_clade_probs.begin(); mcp != miss_clade_probs.end(); ){
        double p = 0.0;
        for (auto mcp2 = mcp->second.begin(); mcp2 != mcp->second.end(); ){
            if (downsample_prob.count(mcp2->first) > 0){
                p += mcp2->second * downsample_prob[mcp2->first];
            }
            else{
                p += mcp2->second;
            }
            mcp->second.erase(mcp2++);
        }
        if (p < 1.0){
            downsample_miss_prob.insert(make_pair(mcp->first, p));
        }
        miss_clade_probs.erase(mcp++);
        miss_count++;
        if (miss_count % 100000 == 0) {
            fprintf(stderr, "  Processed %ld missing entries...\r", miss_count);
            fflush(stderr);
        }
    }
    fprintf(stderr, "  Computed %zu missing genotype probabilities\n", downsample_miss_prob.size());
    fflush(stderr);
    
    if (clade_slot_allocation.size() == 0){
        fprintf(stderr, "No clades were found with more than %d occurrences.\n", num);
        fprintf(stderr, "There is nothing to downsample.\n");
        return 0;
    }

    // ========================================================================
    // STRATIFIED SELECTION (Option C): Select top SNPs within each clade
    // ========================================================================
    
    fprintf(stderr, "\n========================================\n");
    fprintf(stderr, "Stratified Selection: Selecting top SNPs within each clade...\n");
    start_time = chrono::high_resolution_clock::now();
    
    // Step 1: Group all QC-passing biallelic SNPs by their clade
    // Map from clade bitset -> vector of (variant_idx, selection_score)
    unordered_map<bitset<NBITS>, vector<CladeCandidate>> clade_to_variants;
    
    // Also need to handle missing genotype sites - map to their probable parent clades
    // For simplicity, we'll use the same parent-mapping logic as the probability computation
    
    fprintf(stderr, "  Grouping %zu variants by clade...\n", all_variants.size());
    fflush(stderr);
    
    long grouped_count = 0;
    long skipped_low_cov = 0;
    
    for (size_t vi = 0; vi < all_variants.size(); vi++) {
        const StoredVariant& sv = all_variants[vi];
        
        // Only consider biallelic SNPs that pass QC
        if (!sv.is_biallelic || !sv.passes_qc) continue;
        
        // Get chromosome name for coverage lookup
        string chrom = chroms[sv.chrom_idx];
        
        // Get coverage score FRESH (same as het candidate collection)
        float cov_score = 0.0f;
        if (!coverage_map.empty()) {
            cov_score = get_coverage_score_fast(chrom, sv.pos, coverage_map);
        }
        
        // Get annotation score FRESH
        float annot_score = 0.1f;
        if (!annotations.empty()) {
            annot_score = get_annot_score(chrom, sv.pos, annotations);
        }
        
        // Apply min_cov hard filter
        if (cov_score < min_cov) {
            skipped_low_cov++;
            continue;
        }
        
        // Determine the canonical clade bitset
        int ac = sv.alt.count();
        int afc = sv.alt_flip.count();
        bitset<NBITS> clade_key;
        
        if (sv.present.count() == (size_t)num_samples) {
            // No missing data - use direct clade
            if (ac < afc || (ac == afc && sv.alt < sv.alt_flip)) {
                clade_key = sv.alt;
            } else {
                clade_key = sv.alt_flip;
            }
        } else {
            // Has missing data - need to find parent clade
            // For now, skip sites with missing data in stratified selection
            // (they contribute fractionally to multiple parent clades, complex to handle)
            continue;
        }
        
        // Only add if this clade has a slot allocation
        if (clade_slot_allocation.count(clade_key) > 0) {
            CladeCandidate cc;
            cc.variant_idx = vi;
            cc.selection_score = compute_selection_score(cov_score, annot_score);
            clade_to_variants[clade_key].push_back(cc);
            grouped_count++;
        }
        
        if (grouped_count % 5000000 == 0 && grouped_count > 0) {
            fprintf(stderr, "    Grouped %ld variants...\r", grouped_count);
            fflush(stderr);
        }
    }
    
    fprintf(stderr, "  Grouped %ld variants into %zu clades (skipped %ld below min_cov=%.1f)\n", 
        grouped_count, clade_to_variants.size(), skipped_low_cov, min_cov);
    fflush(stderr);
    
    // Step 2: Sort each clade's variants by selection score and select top N
    fprintf(stderr, "  Selecting top variants within each clade...\n");
    fflush(stderr);
    
    // Set of selected variant indices
    set<size_t> selected_variant_indices;
    
    long total_shortfall = 0;
    long clades_with_shortfall = 0;
    vector<pair<bitset<NBITS>, long>> clades_needing_more;  // (clade, extra_capacity)
    
    for (auto& kv : clade_to_variants) {
        const bitset<NBITS>& clade = kv.first;
        vector<CladeCandidate>& candidates = kv.second;
        
        long slots = clade_slot_allocation[clade];
        long available = candidates.size();
        
        // Sort by selection score (descending)
        sort(candidates.begin(), candidates.end());
        
        // Select top min(slots, available)
        long to_select = min(slots, available);
        for (long i = 0; i < to_select; i++) {
            selected_variant_indices.insert(candidates[i].variant_idx);
        }
        
        if (available < slots) {
            // Shortfall: couldn't fill all allocated slots
            long shortfall = slots - available;
            total_shortfall += shortfall;
            clades_with_shortfall++;
        } else if (available > slots) {
            // This clade has extra capacity for redistribution
            clades_needing_more.push_back(make_pair(clade, available - slots));
        }
    }
    
    fprintf(stderr, "  Initial selection: %zu variants\n", selected_variant_indices.size());
    fprintf(stderr, "  Shortfall: %ld slots from %ld clades\n", total_shortfall, clades_with_shortfall);
    fflush(stderr);
    
    // Step 3: Redistribute shortfall proportionally to clades with extra capacity
    if (total_shortfall > 0 && !clades_needing_more.empty()) {
        fprintf(stderr, "  Redistributing %ld shortfall slots...\n", total_shortfall);
        fflush(stderr);
        
        // Calculate total extra capacity
        long total_extra_capacity = 0;
        for (auto& cn : clades_needing_more) {
            total_extra_capacity += cn.second;
        }
        
        // Redistribute proportionally
        long redistributed = 0;
        for (auto& cn : clades_needing_more) {
            const bitset<NBITS>& clade = cn.first;
            long extra_capacity = cn.second;
            
            // Proportional share of shortfall
            long extra_slots = (long)round((double)total_shortfall * extra_capacity / total_extra_capacity);
            if (extra_slots > extra_capacity) extra_slots = extra_capacity;
            
            // Get this clade's candidates (already sorted)
            vector<CladeCandidate>& candidates = clade_to_variants[clade];
            long already_selected = clade_slot_allocation[clade];
            
            // Add more from this clade
            for (long i = already_selected; i < already_selected + extra_slots && i < (long)candidates.size(); i++) {
                if (selected_variant_indices.insert(candidates[i].variant_idx).second) {
                    redistributed++;
                }
            }
        }
        
        fprintf(stderr, "  Redistributed %ld additional slots\n", redistributed);
        fprintf(stderr, "  Final selection: %zu variants\n", selected_variant_indices.size());
    }
    
    end_time = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    fprintf(stderr, "Stratified selection complete in %ld seconds\n", duration.count());
    fprintf(stderr, "  Target: %d, Selected: %zu, Shortfall redistributed: %ld\n", 
        num, selected_variant_indices.size(), total_shortfall);
    fflush(stderr);

    // ========================================================================
    // COMPUTE DEMUX_SCORE FOR HET CANDIDATES AND WRITE HET OUTPUT
    // ========================================================================
    
    if (collect_het && !all_het_candidates.empty() && !selected_het_indices.empty()) {
        fprintf(stderr, "\n========================================\n");
        fprintf(stderr, "Computing DEMUX_SCORE for het candidates and writing output...\n");
        start_time = chrono::high_resolution_clock::now();
        
        // Compute DEMUX_SCORE for each het candidate from stored genotypes
        fprintf(stderr, "  Computing DEMUX_SCORE for %zu selected het candidates...\n", selected_het_indices.size());
        
        #pragma omp parallel for schedule(dynamic, 1000)
        for (size_t i = 0; i < selected_het_indices.size(); ++i) {
            HetCandidate& hc = all_het_candidates[selected_het_indices[i]];
            
            // Reconstruct alt/alt_flip/present bitsets from stored genotypes
            bitset<NBITS> alt, alt_flip, present;
            int ploidy = hc.genotypes.size() / hc.n_samples;
            
            for (int s = 0; s < hc.n_samples; ++s) {
                int32_t* gt_ptr = hc.genotypes.data() + s * ploidy;
                
                if (bcf_gt_is_missing(gt_ptr[0]) || (ploidy > 1 && bcf_gt_is_missing(gt_ptr[1]))) {
                    continue;  // Missing - not present
                }
                
                present.set(s);
                
                int a1 = bcf_gt_allele(gt_ptr[0]);
                int a2 = (ploidy > 1) ? bcf_gt_allele(gt_ptr[1]) : a1;
                
                if (a1 > 0 || a2 > 0) {
                    alt.set(s);
                }
                if (a1 == 0 || a2 == 0) {
                    alt_flip.set(s);
                }
            }
            
            // Look up downsample probability
            double samp_prob = 1.0;
            int ac = alt.count();
            int afc = alt_flip.count();
            
            if (present.count() == (size_t)hc.n_samples) {
                // No missing data
                if (ac < afc || (ac == afc && alt < alt_flip)) {
                    if (downsample_prob.count(alt) > 0) {
                        samp_prob = downsample_prob.at(alt);
                    }
                } else {
                    if (downsample_prob.count(alt_flip) > 0) {
                        samp_prob = downsample_prob.at(alt_flip);
                    }
                }
            } else {
                // Has missing data
                pair<bitset<NBITS>, bitset<NBITS>> key;
                if (ac < afc || (ac == afc && alt < alt_flip)) {
                    key = make_pair(present, alt);
                } else {
                    key = make_pair(present, alt_flip);
                }
                if (downsample_miss_prob.count(key) > 0) {
                    samp_prob = downsample_miss_prob.at(key);
                }
            }
            
            hc.demux_score = (samp_prob < 1.0) ? (float)samp_prob : 1.0f;
        }
        
        fprintf(stderr, "  Writing het output to %s...\n", het_outfile.c_str());
        
        // Create output header
        htsFile* hdr_reader = hts_open(vcf_file.c_str(), "r");
        bcf_hdr_t* het_header = bcf_hdr_read(hdr_reader);
        hts_close(hdr_reader);
        
        // Add INFO headers - ALL fields
        bcf_hdr_append(het_header, "##INFO=<ID=DEMUX_SCORE,Number=1,Type=Float,Description=\"Clade rarity score\">");
        bcf_hdr_append(het_header, "##INFO=<ID=HET_SCORE,Number=1,Type=Float,Description=\"Combined het site quality score\">");
        bcf_hdr_append(het_header, "##INFO=<ID=HET_FREQ,Number=1,Type=Float,Description=\"Fraction of called individuals that are het\">");
        bcf_hdr_append(het_header, "##INFO=<ID=MISSING_FREQ,Number=1,Type=Float,Description=\"Fraction of individuals with missing genotype\">");
        bcf_hdr_append(het_header, "##INFO=<ID=ANNOT_SCORE,Number=1,Type=Float,Description=\"Annotation weight (1.0=genic, 0.1=intergenic)\">");
        bcf_hdr_append(het_header, "##INFO=<ID=COV_SCORE,Number=1,Type=Float,Description=\"Normalized coverage score\">");
        bcf_hdr_append(het_header, "##INFO=<ID=BIN_ID,Number=1,Type=String,Description=\"Chromosome:bin_index for evenness tracking\">");
        if (bcf_hdr_sync(het_header) < 0) {
            fprintf(stderr, "ERROR: Failed to sync het header\n");
            exit(1);
        }
        
        // Open output BCF file
        htsFile* het_out = hts_open(het_outfile.c_str(), "wb");
        hts_set_threads(het_out, n_threads);
        if (bcf_hdr_write(het_out, het_header) < 0) {
            fprintf(stderr, "ERROR: Failed to write het header\n");
            exit(1);
        }
        
        // Sort selected candidates by chromosome then position for ordered output
        vector<HetCandidate*> selected_sorted;
        for (size_t idx : selected_het_indices) {
            selected_sorted.push_back(&all_het_candidates[idx]);
        }
        sort(selected_sorted.begin(), selected_sorted.end(), 
            [](const HetCandidate* a, const HetCandidate* b) {
                if (a->chrom_idx != b->chrom_idx) return a->chrom_idx < b->chrom_idx;
                return a->pos < b->pos;
            });
        
        // Write records
        bcf1_t* het_rec = bcf_init();
        for (const HetCandidate* hc : selected_sorted) {
            bcf_clear(het_rec);
            
            // Set chromosome and position
            het_rec->rid = bcf_hdr_name2id(het_header, hc->chrom.c_str());
            het_rec->pos = hc->pos;
            het_rec->qual = 0;
            
            // Set alleles
            string alleles = hc->ref_allele + "," + hc->alt_allele;
            bcf_update_alleles_str(het_header, het_rec, alleles.c_str());
            
            // Set genotypes
            bcf_update_genotypes(het_header, het_rec, hc->genotypes.data(), hc->genotypes.size());
            
            // Set INFO fields - ALL scores
            float demux_score = hc->demux_score;
            float het_score = hc->het_score;
            float het_freq = hc->het_freq;
            float missing_freq = hc->missing_freq;
            float annot = hc->annot_score;
            float cov = hc->cov_score;
            
            bcf_update_info_float(het_header, het_rec, "DEMUX_SCORE", &demux_score, 1);
            bcf_update_info_float(het_header, het_rec, "HET_SCORE", &het_score, 1);
            bcf_update_info_float(het_header, het_rec, "HET_FREQ", &het_freq, 1);
            bcf_update_info_float(het_header, het_rec, "MISSING_FREQ", &missing_freq, 1);
            bcf_update_info_float(het_header, het_rec, "ANNOT_SCORE", &annot, 1);
            bcf_update_info_float(het_header, het_rec, "COV_SCORE", &cov, 1);
            
            string bin_id = make_bin_id(hc->chrom, hc->bin_idx);
            bcf_update_info_string(het_header, het_rec, "BIN_ID", bin_id.c_str());
            
            if (bcf_write(het_out, het_header, het_rec) < 0) {
                fprintf(stderr, "ERROR: Failed to write het record\n");
                exit(1);
            }
        }
        
        bcf_destroy(het_rec);
        hts_close(het_out);
        bcf_hdr_destroy(het_header);
        
        // Index het output BCF
        fprintf(stderr, "  Indexing het output BCF...\n");
        if (bcf_index_build(het_outfile.c_str(), 14) < 0) {
            fprintf(stderr, "WARNING: Failed to create index for %s\n", het_outfile.c_str());
        }
        
        // Free het candidates (large memory)
        all_het_candidates.clear();
        all_het_candidates.shrink_to_fit();
        selected_het_indices.clear();
        selected_het_indices.shrink_to_fit();
        
        end_time = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
        fprintf(stderr, "Het output complete: %zu sites in %ld seconds\n", 
            selected_het_sites.size(), duration.count());
        fprintf(stderr, "Het output: %s\n", het_outfile.c_str());
    }

    // ========================================================================
    // PASS 2: Write selected variants (Option C: Deterministic from stratified selection)
    // Also writes annotate_all output if requested
    // ========================================================================
    
    fprintf(stderr, "\n========================================\n");
    fprintf(stderr, "Pass 2: Writing selected variants (deterministic, parallel, %d threads)...\n", n_threads);
    if (!annotate_all_file.empty()) {
        fprintf(stderr, "  Also writing full annotated output to: %s\n", annotate_all_file.c_str());
    }
    fflush(stderr);
    start_time = chrono::high_resolution_clock::now();
    
    // vcf_idx is reused from Pass 1
    
    // Get header for INFO field additions
    htsFile* hdr_reader = hts_open(vcf_file.c_str(), "r");
    bcf_hdr_t* out_header = bcf_hdr_read(hdr_reader);
    hts_close(hdr_reader);
    
    // Add ALL INFO headers (used for both demux and annotate_all outputs)
    bcf_hdr_append(out_header, "##INFO=<ID=DEMUX_SCORE,Number=1,Type=Float,Description=\"Clade rarity score\">");
    bcf_hdr_append(out_header, "##INFO=<ID=HET_FREQ,Number=1,Type=Float,Description=\"Fraction of called individuals that are heterozygous\">");
    bcf_hdr_append(out_header, "##INFO=<ID=MISSING_FREQ,Number=1,Type=Float,Description=\"Fraction of individuals with missing genotype\">");
    bcf_hdr_append(out_header, "##INFO=<ID=HET_SCORE,Number=1,Type=Float,Description=\"Combined het score: het_freq * (1 - missing_freq * 0.5)\">");
    bcf_hdr_append(out_header, "##INFO=<ID=ANNOT_SCORE,Number=1,Type=Float,Description=\"Annotation weight (1.0=genic, 0.1=intergenic)\">");
    bcf_hdr_append(out_header, "##INFO=<ID=COV_SCORE,Number=1,Type=Float,Description=\"Coverage score from bedgraph\">");
    bcf_hdr_append(out_header, "##INFO=<ID=BIN_ID,Number=1,Type=String,Description=\"Chromosome:bin_index for evenness tracking\">");
    bcf_hdr_append(out_header, "##INFO=<ID=SELECTION_SCORE,Number=1,Type=Float,Description=\"Within-clade selection score: log2(cov+1) + 0.1*annot_boost\">");
    if (bcf_hdr_sync(out_header) < 0) {
        fprintf(stderr, "ERROR: Failed to sync header\n");
        exit(1);
    }
    
    // Create temp directory for per-chromosome demux outputs
    string temp_dir = outfile + ".tmp";
    mkdir(temp_dir.c_str(), 0755);
    fprintf(stderr, "  Demux temp directory: %s\n", temp_dir.c_str());
    
    // Create temp directory for per-chromosome annotate_all outputs (if needed)
    string annot_temp_dir;
    if (!annotate_all_file.empty()) {
        annot_temp_dir = annotate_all_file + ".tmp";
        mkdir(annot_temp_dir.c_str(), 0755);
        fprintf(stderr, "  Annotate temp directory: %s\n", annot_temp_dir.c_str());
    }
    fflush(stderr);
    
    atomic<long> total_kept(0);
    atomic<long> total_elim(0);
    atomic<long> total_annotated(0);  // For annotate_all tracking
    atomic<int> chroms_written(0);
    
    // Vector to track which chromosomes have output (for concatenation order)
    vector<string> chrom_files(chroms.size());
    vector<string> annot_chrom_files(chroms.size());  // For annotate_all
    vector<bool> chrom_has_data(chroms.size(), false);
    vector<bool> annot_chrom_has_data(chroms.size(), false);
    
    #pragma omp parallel num_threads(n_threads)
    {
        // Thread-local copy of output header
        bcf_hdr_t* thread_out_header = bcf_hdr_dup(out_header);
        
        // Thread-local output record
        bcf1_t* out_rec = bcf_init();
        
        #pragma omp for schedule(dynamic, 1)
        for (size_t c = 0; c < chroms.size(); c++) {
            const string& chrom = chroms[c];
            
            // Get range of variants for this chromosome from memory
            size_t range_start = chrom_ranges[c].first;
            size_t range_end = chrom_ranges[c].second;
            
            if (range_start >= range_end) {
                chroms_written++;
                continue;
            }
            
            // Create temp output file for demux (filtered) output
            string temp_file = temp_dir + "/chr_" + to_string(c) + ".bcf";
            chrom_files[c] = temp_file;
            htsFile* temp_out = hts_open(temp_file.c_str(), "wb");
            (void)bcf_hdr_write(temp_out, thread_out_header);
            
            // Create temp output file for annotate_all (all sites) if requested
            htsFile* annot_out = NULL;
            if (!annotate_all_file.empty()) {
                string annot_file = annot_temp_dir + "/chr_" + to_string(c) + ".bcf";
                annot_chrom_files[c] = annot_file;
                annot_out = hts_open(annot_file.c_str(), "wb");
                (void)bcf_hdr_write(annot_out, thread_out_header);
            }
            
            long local_kept = 0;
            long local_elim = 0;
            long local_annotated = 0;
            
            // Iterate over stored variants for this chromosome
            for (size_t vi = range_start; vi < range_end; vi++) {
                const StoredVariant& sv = all_variants[vi];
                
                // For annotate_all, we process ALL records (including multi-allelic)
                // For demux, we only process biallelic
                bool is_biallelic = sv.is_biallelic;
                
                // Compute DEMUX_SCORE (for reporting) and filtering decision
                // Option C: Use pre-computed selected_variant_indices set (deterministic)
                double samp_prob = 1.0;
                bool passes_filter = false;
                
                if (sv.passes_qc) {
                    // Use precomputed bitsets (same as what proc_bcf_record would return)
                    int ac = sv.alt.count();
                    int afc = sv.alt_flip.count();
                    
                    if (sv.present.count() == (size_t)num_samples) {
                        if (ac < afc || (ac == afc && sv.alt < sv.alt_flip)) {
                            if (downsample_prob.count(sv.alt) > 0) {
                                samp_prob = downsample_prob[sv.alt];
                            }
                        } else {
                            if (downsample_prob.count(sv.alt_flip) > 0) {
                                samp_prob = downsample_prob[sv.alt_flip];
                            }
                        }
                    } else {
                        pair<bitset<NBITS>, bitset<NBITS>> key;
                        if (ac < afc || (ac == afc && sv.alt < sv.alt_flip)) {
                            key = make_pair(sv.present, sv.alt);
                        } else {
                            key = make_pair(sv.present, sv.alt_flip);
                        }
                        if (downsample_miss_prob.count(key) > 0) {
                            samp_prob = downsample_miss_prob[key];
                        }
                    }
                    
                    // Option C: Deterministic selection using pre-computed set
                    // (replaces random sampling)
                    passes_filter = (selected_variant_indices.count(vi) > 0);
                }
                
                // Compute ALL scores for this record
                float demux_score = (samp_prob < 1.0) ? (float)samp_prob : 1.0f;
                // Get annotation and coverage scores FRESH (consistent with het and stratified selection)
                float annot_score = 0.1f;
                if (!annotations.empty()) {
                    annot_score = get_annot_score(chrom, sv.pos, annotations);
                }
                float cov_score = 0.0f;
                if (!coverage_map.empty()) {
                    cov_score = get_coverage_score_fast(chrom, sv.pos, coverage_map);
                }
                int bin_idx = get_bin_idx(sv.pos, bin_size);
                string bin_id = make_bin_id(chrom, bin_idx);
                
                // Compute HET scores from genotypes (EXACT same logic as original)
                float het_freq = 0.0f;
                float missing_freq = 0.0f;
                float het_score = 0.0f;
                
                if (!sv.genotypes.empty()) {
                    int ploidy = sv.genotypes.size() / num_samples;  // Same as: ret / num_samples
                    int n_het = 0;
                    int n_called = 0;
                    int n_missing = 0;
                    
                    for (int s = 0; s < num_samples; s++) {
                        const int32_t* gt_ptr = sv.genotypes.data() + s * ploidy;
                        
                        if (bcf_gt_is_missing(gt_ptr[0]) || 
                            (ploidy > 1 && bcf_gt_is_missing(gt_ptr[1]))) {
                            n_missing++;
                            continue;
                        }
                        
                        n_called++;
                        
                        if (ploidy >= 2) {
                            int a1 = bcf_gt_allele(gt_ptr[0]);
                            int a2 = bcf_gt_allele(gt_ptr[1]);
                            if (a1 != a2) {
                                n_het++;
                            }
                        }
                    }
                    
                    het_freq = (n_called > 0) ? (float)n_het / n_called : 0.0f;
                    missing_freq = (float)n_missing / num_samples;
                    float missing_penalty = 1.0f - (missing_freq * 0.5f);
                    het_score = het_freq * missing_penalty;
                }
                
                // Build output record from stored data
                bcf_clear(out_rec);
                out_rec->rid = sv.chrom_idx;
                out_rec->pos = sv.pos;
                out_rec->qual = 0;
                
                // Set alleles
                if (!sv.alleles.empty()) {
                    bcf_update_alleles_str(thread_out_header, out_rec, sv.alleles.c_str());
                }
                
                // Set genotypes
                if (!sv.genotypes.empty()) {
                    bcf_update_genotypes(thread_out_header, out_rec, sv.genotypes.data(), sv.genotypes.size());
                }
                
                // Add ALL INFO fields
                bcf_update_info_float(thread_out_header, out_rec, "DEMUX_SCORE", &demux_score, 1);
                bcf_update_info_float(thread_out_header, out_rec, "HET_FREQ", &het_freq, 1);
                bcf_update_info_float(thread_out_header, out_rec, "MISSING_FREQ", &missing_freq, 1);
                bcf_update_info_float(thread_out_header, out_rec, "HET_SCORE", &het_score, 1);
                bcf_update_info_float(thread_out_header, out_rec, "ANNOT_SCORE", &annot_score, 1);
                bcf_update_info_float(thread_out_header, out_rec, "COV_SCORE", &cov_score, 1);
                bcf_update_info_string(thread_out_header, out_rec, "BIN_ID", bin_id.c_str());
                
                // Write to annotate_all output (ALL records)
                if (annot_out) {
                    bcf_write1(annot_out, thread_out_header, out_rec);
                    local_annotated++;
                }
                
                // Write to demux output (filtered biallelic records only)
                if (is_biallelic && passes_filter) {
                    bcf_write1(temp_out, thread_out_header, out_rec);
                    local_kept++;
                } else if (is_biallelic) {
                    local_elim++;
                }
            }
            
            hts_close(temp_out);
            if (annot_out) {
                hts_close(annot_out);
            }
            
            if (local_kept > 0) {
                chrom_has_data[c] = true;
            } else {
                // Remove empty file
                remove(temp_file.c_str());
            }
            
            if (local_annotated > 0) {
                annot_chrom_has_data[c] = true;
            } else if (annot_out) {
                remove(annot_chrom_files[c].c_str());
            }
            
            total_kept += local_kept;
            total_elim += local_elim;
            total_annotated += local_annotated;
            int done = ++chroms_written;
            
            if (local_kept > 0 || done % 100 == 0) {
                fprintf(stderr, "  Chromosome %s: kept %ld [%d/%zu]\n", 
                    chrom.c_str(), local_kept, done, chroms.size());
                fflush(stderr);
            }
        }
        
        bcf_destroy(out_rec);
        bcf_hdr_destroy(thread_out_header);
    }
    
    hts_idx_destroy(vcf_idx);
    
    end_time = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    
    fprintf(stderr, "\nPass 2 complete: %ld demux kept, %ld eliminated", 
        total_kept.load(), total_elim.load());
    if (!annotate_all_file.empty()) {
        fprintf(stderr, ", %ld annotated", total_annotated.load());
    }
    fprintf(stderr, " in %ld seconds\n", duration.count());
    
    // ========================================================================
    // Concatenate chromosome files in order
    // ========================================================================
    
    fprintf(stderr, "Concatenating chromosome files...\n");
    start_time = chrono::high_resolution_clock::now();
    
    htsFile* outf = hts_open(outfile.c_str(), "wb");
    hts_set_threads(outf, n_threads);
    if (bcf_hdr_write(outf, out_header) < 0) {
        fprintf(stderr, "ERROR: Failed to write output header\n");
        exit(1);
    }
    
    bcf1_t* rec = bcf_init();
    long final_count = 0;
    
    for (size_t c = 0; c < chroms.size(); c++) {
        if (!chrom_has_data[c]) continue;
        
        htsFile* temp_in = hts_open(chrom_files[c].c_str(), "r");
        bcf_hdr_t* temp_hdr = bcf_hdr_read(temp_in);
        
        while (bcf_read(temp_in, temp_hdr, rec) >= 0) {
            bcf_write1(outf, out_header, rec);
            final_count++;
        }
        
        bcf_hdr_destroy(temp_hdr);
        hts_close(temp_in);
        
        // Remove temp file
        remove(chrom_files[c].c_str());
    }
    
    bcf_destroy(rec);
    hts_close(outf);
    
    // Index demux output BCF
    fprintf(stderr, "Indexing demux output BCF...\n");
    if (bcf_index_build(outfile.c_str(), 14) < 0) {
        fprintf(stderr, "WARNING: Failed to create index for %s\n", outfile.c_str());
    }
    
    // Remove demux temp directory
    rmdir(temp_dir.c_str());
    
    fprintf(stderr, "Demux output: %ld SNPs written to %s\n", final_count, outfile.c_str());
    
    // ========================================================================
    // Concatenate annotate_all chromosome files (if requested)
    // ========================================================================
    
    if (!annotate_all_file.empty()) {
        fprintf(stderr, "\nConcatenating annotate_all chromosome files...\n");
        
        htsFile* annot_outf = hts_open(annotate_all_file.c_str(), "wb");
        hts_set_threads(annot_outf, n_threads);
        if (bcf_hdr_write(annot_outf, out_header) < 0) {
            fprintf(stderr, "ERROR: Failed to write annotate_all header\n");
            exit(1);
        }
        
        bcf1_t* annot_rec = bcf_init();
        long annot_final_count = 0;
        
        for (size_t c = 0; c < chroms.size(); c++) {
            if (!annot_chrom_has_data[c]) continue;
            
            htsFile* temp_in = hts_open(annot_chrom_files[c].c_str(), "r");
            bcf_hdr_t* temp_hdr = bcf_hdr_read(temp_in);
            
            while (bcf_read(temp_in, temp_hdr, annot_rec) >= 0) {
                bcf_write1(annot_outf, out_header, annot_rec);
                annot_final_count++;
            }
            
            bcf_hdr_destroy(temp_hdr);
            hts_close(temp_in);
            
            // Remove temp file
            remove(annot_chrom_files[c].c_str());
        }
        
        bcf_destroy(annot_rec);
        hts_close(annot_outf);
        
        // Index annotate_all output BCF
        fprintf(stderr, "Indexing annotate_all output BCF...\n");
        if (bcf_index_build(annotate_all_file.c_str(), 14) < 0) {
            fprintf(stderr, "WARNING: Failed to create index for %s\n", annotate_all_file.c_str());
        }
        
        // Remove annotate_all temp directory
        rmdir(annot_temp_dir.c_str());
        
        fprintf(stderr, "Annotate_all output: %ld variants written to %s\n", annot_final_count, annotate_all_file.c_str());
    }
    
    bcf_hdr_destroy(out_header);
    
    end_time = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    
    fprintf(stderr, "\n========================================\n");
    fprintf(stderr, "Pipeline complete in %ld seconds\n", duration.count());
    fprintf(stderr, "========================================\n");
    fprintf(stderr, "Demux output: %s (%ld SNPs)\n", outfile.c_str(), final_count);
    if (!annotate_all_file.empty()) {
        fprintf(stderr, "Annotated output: %s (all variants)\n", annotate_all_file.c_str());
    }
    if (!het_outfile.empty()) {
        fprintf(stderr, "Het output: %s (%zu sites)\n", het_outfile.c_str(), selected_het_sites.size());
    }

    return 0;
}
