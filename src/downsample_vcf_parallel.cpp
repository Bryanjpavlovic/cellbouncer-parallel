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
#include <cstdlib>
#include <random>
#include <functional>
#include <utility>
#include <math.h>
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
    fprintf(stderr, "    --num -n     Desired final number of SNPs.\n");
    fprintf(stderr, "    --output -o  Output VCF file name. Will be gzipped.\n");
    fprintf(stderr, "===== OPTIONAL =====\n");
    fprintf(stderr, "    --threads -t Number of threads (default: all available).\n");
    fprintf(stderr, "    --gtf -g     GTF annotation file. Adds ANNOT_SCORE INFO field.\n");
    fprintf(stderr, "    --cov -c     Tabix-indexed bedgraph coverage file. Adds COV_SCORE INFO field.\n");
    fprintf(stderr, "    --seed -s    Random seed for reproducibility.\n");
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


int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"vcf", required_argument, 0, 'v'},
       {"output", required_argument, 0, 'o'},
       {"num", required_argument, 0, 'n'},
       {"threads", required_argument, 0, 't'},
       {"gtf", required_argument, 0, 'g'},
       {"cov", required_argument, 0, 'c'},
       {"seed", required_argument, 0, 's'},
       {0, 0, 0, 0} 
    };
    
    string vcf_file = "";
    string outfile = "";
    string gtf_file = "";
    string cov_file = "";
    int num = -1;
    int n_threads = omp_get_max_threads();
    int seed = -1;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:n:o:t:g:c:s:h", long_options, &option_index )) != -1){
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
            default:
                help(0);
                break;
        }    
    }
    
    if (vcf_file == ""){
        fprintf(stderr, "ERROR: VCF file required\n");
        exit(1);
    }
    if (num <= 0){
        fprintf(stderr, "ERROR: --num/-n must be positive (Recommend > 1M)\n");
        exit(1);
    }
    if (outfile == ""){
        fprintf(stderr, "ERROR: --output/-o is required.\n");
        exit(1);
    }
    if (outfile.rfind(".vcf") == string::npos){
        outfile += ".vcf.gz";
    }
    else if (outfile.rfind(".gz") == string::npos){
        outfile += ".gz";
    }

    fprintf(stderr, "Using %d threads\n", n_threads);
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
    int n_seqs = 0;
    const char** seqnames = bcf_hdr_seqnames(bcf_header, &n_seqs);
    for (int i = 0; i < n_seqs; ++i){
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
    // PASS 1: Count clades (PARALLEL BY CHROMOSOME)
    // ========================================================================
    
    unordered_map<bitset<NBITS>, double> branchcounts;
    unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, double> branchcounts_missing;
    unordered_map<bitset<NBITS>, bitset<NBITS>> miss2flip;
    
    mutex bc_mutex, bcm_mutex, m2f_mutex;
    atomic<long> total_snps(0);
    atomic<int> chroms_done(0);

    auto start_time = chrono::high_resolution_clock::now();
    
    fprintf(stderr, "Pass 1: Counting mutations on branches (parallel, %d threads)...\n", n_threads);

    #pragma omp parallel num_threads(n_threads)
    {
        // Thread-local data structures
        unordered_map<bitset<NBITS>, double> local_branchcounts;
        unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, double> local_branchcounts_missing;
        unordered_map<bitset<NBITS>, bitset<NBITS>> local_miss2flip;
        
        // Thread-local VCF reader (each thread needs its own)
        htsFile* thread_reader = bcf_open(vcf_file.c_str(), "r");
        bcf_hdr_t* thread_header = bcf_hdr_read(thread_reader);
        
        bcf1_t* bcf_record = bcf_init();
        bitset<NBITS> alt, alt_flip, present;
        
        #pragma omp for schedule(dynamic, 1)
        for (size_t c = 0; c < chroms.size(); c++){
            const string& chrom = chroms[c];
            
            // Query this chromosome using shared index
            hts_itr_t* iter = bcf_itr_querys(vcf_idx, thread_header, chrom.c_str());
            if (!iter) {
                chroms_done++;
                continue;
            }
            
            long local_snps = 0;
            
            while (bcf_itr_next(thread_reader, iter, bcf_record) >= 0){
                if (bcf_record->n_allele != 2) continue;
                
                bool pass;
                if (proc_bcf_record(bcf_record, thread_header, num_samples,
                    alt, alt_flip, present, pass)){
                    
                    int ac = alt.count();
                    int afc = alt_flip.count();
                
                    if (present.count() == num_samples){
                        count_branch(local_branchcounts, alt, alt_flip, ac, afc, 1.0);
                    }
                    else{
                        count_branch_missing(local_branchcounts_missing, local_miss2flip, 
                            present, alt, alt_flip, ac, afc, 1.0);
                    }
                    local_snps++;
                }
            }
            
            hts_itr_destroy(iter);
            total_snps += local_snps;
            int done = ++chroms_done;
            
            if (local_snps > 0 || done % 100 == 0) {
                fprintf(stderr, "  Chromosome %s: %ld SNPs [%d/%zu]\n", 
                    chrom.c_str(), local_snps, done, chroms.size());
            }
        }
        
        // Merge thread-local counts into global
        merge_branchcounts(branchcounts, local_branchcounts, bc_mutex);
        merge_branchcounts_missing(branchcounts_missing, local_branchcounts_missing, bcm_mutex);
        merge_miss2flip(miss2flip, local_miss2flip, m2f_mutex);
        
        bcf_destroy(bcf_record);
        bcf_hdr_destroy(thread_header);
        hts_close(thread_reader);
    }
    
    // Clean up shared index
    hts_idx_destroy(vcf_idx);
    
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    
    fprintf(stderr, "\nPass 1 complete: %ld SNPs in %ld seconds\n", total_snps.load(), duration.count());
    fprintf(stderr, "Found %lu unique clade patterns\n", branchcounts.size());
    
    if (branchcounts.size() >= (size_t)num){
        fprintf(stderr, "ERROR: %zu distinct allele sharing patterns found, but %d SNPs requested.\n", 
            branchcounts.size(), num);
        fprintf(stderr, "  Please set sample size to at least %zu.\n", branchcounts.size());
        exit(1);
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
    // Compute downsampling probabilities using Brent's method
    // ========================================================================
    
    unordered_map<bitset<NBITS>, double> downsample_prob;
    
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
                        
                        fprintf(stderr, "Building downsample probability map for %d clades...\n", n_hi);
                        fflush(stderr);
                        for (int i = 0; i < n_hi; ++i){
                            double dsp = pow(counts[i], res) / counts[i];
                            downsample_prob.insert(make_pair(bs_idx[clcountsort[i].second], dsp));
                        }
                        fprintf(stderr, "  Added %zu probabilities\n", downsample_prob.size());
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
    
    if (downsample_prob.size() == 0){
        fprintf(stderr, "No clades were found with more than %d occurrences.\n", num);
        fprintf(stderr, "There is nothing to downsample.\n");
        return 0;
    }

    // ========================================================================
    // PASS 2: Filter VCF (PARALLEL BY CHROMOSOME)
    // ========================================================================
    
    fprintf(stderr, "Pass 2: Filtering VCF (parallel, %d threads)...\n", n_threads);
    fflush(stderr);
    start_time = chrono::high_resolution_clock::now();
    
    // Reload index for Pass 2
    vcf_idx = bcf_index_load(vcf_file.c_str());
    if (!vcf_idx) {
        fprintf(stderr, "ERROR: Could not reload VCF index for Pass 2\n");
        exit(1);
    }
    
    // Get header for INFO field additions
    htsFile* hdr_reader = hts_open(vcf_file.c_str(), "r");
    bcf_hdr_t* out_header = bcf_hdr_read(hdr_reader);
    hts_close(hdr_reader);
    
    // Add INFO headers
    bcf_hdr_append(out_header, "##INFO=<ID=DEMUX_SCORE,Number=1,Type=Float,Description=\"Clade rarity score\">");
    if (!gtf_file.empty()) {
        bcf_hdr_append(out_header, "##INFO=<ID=ANNOT_SCORE,Number=1,Type=Float,Description=\"Annotation weight (1.0=genic, 0.1=intergenic)\">");
    }
    if (!coverage_map.empty()) {
        bcf_hdr_append(out_header, "##INFO=<ID=COV_SCORE,Number=1,Type=Float,Description=\"Coverage score from bedgraph\">");
    }
    if (bcf_hdr_sync(out_header) < 0) {
        fprintf(stderr, "ERROR: Failed to sync header\n");
        exit(1);
    }
    
    // Create temp directory for per-chromosome outputs
    string temp_dir = outfile + ".tmp";
    mkdir(temp_dir.c_str(), 0755);
    fprintf(stderr, "  Temp directory: %s\n", temp_dir.c_str());
    fflush(stderr);
    
    atomic<long> total_kept(0);
    atomic<long> total_elim(0);
    atomic<int> chroms_written(0);
    
    // Vector to track which chromosomes have output (for concatenation order)
    vector<string> chrom_files(chroms.size());
    vector<bool> chrom_has_data(chroms.size(), false);
    
    #pragma omp parallel num_threads(n_threads)
    {
        // Thread-local VCF reader
        htsFile* thread_reader = bcf_open(vcf_file.c_str(), "r");
        bcf_hdr_t* thread_header = bcf_hdr_read(thread_reader);
        
        // Thread-local copy of output header (bcf_update_info_float is NOT thread-safe with shared header)
        bcf_hdr_t* thread_out_header = bcf_hdr_dup(out_header);
        
        // Thread-local random generator with unique seed per thread
        int thread_id = omp_get_thread_num();
        mt19937 gen(seed >= 0 ? seed + thread_id : random_device{}());
        uniform_real_distribution<> rand_dist(0.0, 1.0);
        
        bcf1_t* bcf_record = bcf_init();
        bitset<NBITS> alt, alt_flip, present;
        
        #pragma omp for schedule(dynamic, 1)
        for (size_t c = 0; c < chroms.size(); c++) {
            const string& chrom = chroms[c];
            
            hts_itr_t* iter = bcf_itr_querys(vcf_idx, thread_header, chrom.c_str());
            if (!iter) {
                chroms_written++;
                continue;
            }
            
            // Create temp output file for this chromosome
            string temp_file = temp_dir + "/chr_" + to_string(c) + ".bcf";
            chrom_files[c] = temp_file;
            
            htsFile* temp_out = hts_open(temp_file.c_str(), "wb");
            
            // Write header (needed for valid BCF)
            (void)bcf_hdr_write(temp_out, thread_out_header);
            
            long local_kept = 0;
            long local_elim = 0;
            
            while (bcf_itr_next(thread_reader, iter, bcf_record) >= 0) {
                if (bcf_record->n_allele != 2) continue;
                
                bool pass;
                if (proc_bcf_record(bcf_record, thread_header, num_samples,
                    alt, alt_flip, present, pass)) {
                    
                    double samp_prob = 1.0;
                    int ac = alt.count();
                    int afc = alt_flip.count();
                    
                    if (present.count() == num_samples) {
                        if (ac < afc || (ac == afc && alt < alt_flip)) {
                            if (downsample_prob.count(alt) > 0) {
                                samp_prob = downsample_prob[alt];
                            }
                        } else {
                            if (downsample_prob.count(alt_flip) > 0) {
                                samp_prob = downsample_prob[alt_flip];
                            }
                        }
                    } else {
                        pair<bitset<NBITS>, bitset<NBITS>> key;
                        if (ac < afc || (ac == afc && alt < alt_flip)) {
                            key = make_pair(present, alt);
                        } else {
                            key = make_pair(present, alt_flip);
                        }
                        if (downsample_miss_prob.count(key) > 0) {
                            samp_prob = downsample_miss_prob[key];
                        }
                    }
                    
                    bool print_rec = true;
                    if (samp_prob < 1.0) {
                        double r = rand_dist(gen);
                        if (r > samp_prob) {
                            print_rec = false;
                        }
                    }
                    
                    if (print_rec) {
                        float demux_score = (samp_prob < 1.0) ? (float)samp_prob : 1.0f;
                        bcf_update_info_float(thread_out_header, bcf_record, "DEMUX_SCORE", &demux_score, 1);
                        
                        if (!gtf_file.empty()) {
                            float annot_score = get_annot_score(
                                bcf_hdr_id2name(thread_header, bcf_record->rid),
                                bcf_record->pos,
                                annotations);
                            bcf_update_info_float(thread_out_header, bcf_record, "ANNOT_SCORE", &annot_score, 1);
                        }
                        
                        if (!coverage_map.empty()) {
                            float cov_score = get_coverage_score_fast(
                                bcf_hdr_id2name(thread_header, bcf_record->rid),
                                bcf_record->pos,
                                coverage_map);
                            bcf_update_info_float(thread_out_header, bcf_record, "COV_SCORE", &cov_score, 1);
                        }
                        
                        bcf_write1(temp_out, thread_out_header, bcf_record);
                        local_kept++;
                    } else {
                        local_elim++;
                    }
                } else {
                    local_elim++;
                }
            }
            
            hts_itr_destroy(iter);
            hts_close(temp_out);
            
            if (local_kept > 0) {
                chrom_has_data[c] = true;
            } else {
                // Remove empty file
                remove(temp_file.c_str());
            }
            
            total_kept += local_kept;
            total_elim += local_elim;
            int done = ++chroms_written;
            
            if (local_kept > 0 || done % 100 == 0) {
                fprintf(stderr, "  Chromosome %s: kept %ld [%d/%zu]\n", 
                    chrom.c_str(), local_kept, done, chroms.size());
                fflush(stderr);
            }
        }
        
        bcf_destroy(bcf_record);
        bcf_hdr_destroy(thread_header);
        bcf_hdr_destroy(thread_out_header);
        hts_close(thread_reader);
    }
    
    hts_idx_destroy(vcf_idx);
    
    end_time = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    
    fprintf(stderr, "\nPass 2 filtering complete: %ld kept, %ld eliminated in %ld seconds\n", 
        total_kept.load(), total_elim.load(), duration.count());
    
    // ========================================================================
    // Concatenate chromosome files in order
    // ========================================================================
    
    fprintf(stderr, "Concatenating chromosome files...\n");
    start_time = chrono::high_resolution_clock::now();
    
    htsFile* outf = hts_open(outfile.c_str(), "wz");
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
    bcf_hdr_destroy(out_header);
    
    // Remove temp directory
    rmdir(temp_dir.c_str());
    
    end_time = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    
    fprintf(stderr, "Concatenation complete: %ld SNPs in %ld seconds\n", final_count, duration.count());
    fprintf(stderr, "\nOutput: %s\n", outfile.c_str());

    return 0;
}
