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

    htsFile* cov_fp = NULL;
    tbx_t* cov_tbx = NULL;
    if (!cov_file.empty()) {
        cov_fp = hts_open(cov_file.c_str(), "r");
        cov_tbx = tbx_index_load(cov_file.c_str());
        if (!cov_tbx) {
            fprintf(stderr, "ERROR: Coverage file must be bgzipped and tabix indexed.\n");
            exit(1);
        }
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
    // PASS 1: Count clades (PARALLEL BY CHROMOSOME)
    // ========================================================================
    
    unordered_map<bitset<NBITS>, double> branchcounts;
    unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, double> branchcounts_missing;
    unordered_map<bitset<NBITS>, bitset<NBITS>> miss2flip;
    
    mutex bc_mutex, bcm_mutex, m2f_mutex;
    atomic<long> total_snps(0);
    atomic<int> chroms_done(0);

    auto start_time = chrono::high_resolution_clock::now();
    
    fprintf(stderr, "Pass 1: Counting mutations on branches (parallel)...\n");

    #pragma omp parallel
    {
        // Thread-local data structures
        unordered_map<bitset<NBITS>, double> local_branchcounts;
        unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, double> local_branchcounts_missing;
        unordered_map<bitset<NBITS>, bitset<NBITS>> local_miss2flip;
        
        // Thread-local VCF reader
        htsFile* thread_reader = bcf_open(vcf_file.c_str(), "r");
        bcf_hdr_t* thread_header = bcf_hdr_read(thread_reader);
        
        // Load index for region queries
        hts_idx_t* idx = bcf_index_load(vcf_file.c_str());
        if (!idx){
            // Try tabix index
            tbx_t* tbx = tbx_index_load(vcf_file.c_str());
            if (!tbx){
                fprintf(stderr, "ERROR: VCF must be indexed (.tbi or .csi)\n");
                exit(1);
            }
            tbx_destroy(tbx);
            idx = bcf_index_load(vcf_file.c_str());
        }
        
        bcf1_t* bcf_record = bcf_init();
        bitset<NBITS> alt, alt_flip, present;
        
        #pragma omp for schedule(dynamic, 1)
        for (size_t c = 0; c < chroms.size(); c++){
            const string& chrom = chroms[c];
            
            // Query this chromosome
            hts_itr_t* iter = bcf_itr_querys(idx, thread_header, chrom.c_str());
            if (!iter) continue;
            
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
                        count_branch_missing(local_branchcounts_missing,
                            local_miss2flip, present, alt, alt_flip, ac, afc, 1.0);
                    }
                }
                local_snps++;
            }
            
            hts_itr_destroy(iter);
            total_snps += local_snps;
            chroms_done++;
            
            fprintf(stderr, "  Chromosome %s done (%ld SNPs) [%d/%lu]\r", 
                chrom.c_str(), local_snps, chroms_done.load(), chroms.size());
        }
        
        // Merge thread-local counts into global
        merge_branchcounts(branchcounts, local_branchcounts, bc_mutex);
        merge_branchcounts_missing(branchcounts_missing, local_branchcounts_missing, bcm_mutex);
        merge_miss2flip(miss2flip, local_miss2flip, m2f_mutex);
        
        bcf_destroy(bcf_record);
        hts_idx_destroy(idx);
        bcf_hdr_destroy(thread_header);
        hts_close(thread_reader);
    }
    
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    
    fprintf(stderr, "\nPass 1 complete: %ld SNPs in %ld seconds\n", total_snps.load(), duration.count());
    fprintf(stderr, "Found %lu unique clade patterns\n", branchcounts.size());
    
    if (branchcounts.size() >= (size_t)num){
        fprintf(stderr, "ERROR: %ld distinct allele sharing patterns found, but %d SNPs requested.\n", 
            branchcounts.size(), num);
        fprintf(stderr, "  Please set sample size to at least %ld.\n", branchcounts.size());
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
    
    fprintf(stderr, "Adjusting counts using data from sites with missing genotypes...\n");
    unordered_map<bitset<NBITS>, vector<bitset<NBITS>>> keyuniqparent;
    unordered_map<bitset<NBITS>, float> to_add;
    unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, unordered_map<bitset<NBITS>, float>> miss_clade_probs;
    
    for (auto& bcm : branchcounts_missing){
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

    double allbc2 = 0.0;
    vector<bitset<NBITS>> bs_idx;
    vector<pair<double, int>> clcountsort;  // (negative_count, index)
    
    // Build sorted list
    for (auto bc = branchcounts.begin(); bc != branchcounts.end(); ){
        allbc2 += bc->second;
        clcountsort.push_back(make_pair(-bc->second, bs_idx.size()));
        bs_idx.push_back(bc->first);
        branchcounts.erase(bc++);
    }
    
    fprintf(stderr, "%f total clade counts\n", allbc2);
    
    if (allbc2 > (double)num){
        sort(clcountsort.begin(), clcountsort.end());
        
        for (int n_hi = 1; n_hi <= (int)clcountsort.size(); ++n_hi){
            double sum_hi = 0;
            
            optimML::brent_solver solver(brentfun, dbrentfun);
            solver.set_root();
            char buf[50];
            vector<vector<double>> vecs;

            for (int i = 0; i < n_hi; ++i){
                vecs.push_back(vector<double>{ -clcountsort[i].first });
                sprintf(&buf[0], "x%d", i);
                string bufstr = buf;
                solver.add_data(bufstr, vecs[vecs.size()-1]);
                sum_hi += -clcountsort[i].first;
            }

            double sum_lo = allbc2 - sum_hi;
            if (n_hi == (int)clcountsort.size() || sum_lo < (double)num){
                double lastcount = -clcountsort[n_hi-1].first;
                double nextcount = -1;
                if (n_hi < (int)clcountsort.size()){
                    nextcount = -clcountsort[n_hi].first;
                }
                
                solver.add_data_fixed("target", (int)round((double)num - sum_lo));
                solver.add_data_fixed("num", n_hi);
                double res = solver.solve(0, 1);

                if (res < 1.0){
                    lastcount = pow(lastcount, res);
                    
                    if (lastcount > nextcount){
                        fprintf(stderr, "Brent solver found x = %f\n", res);
                        for (int i = 0; i < n_hi; ++i){
                            double dsp = pow(-clcountsort[i].first, res) / -clcountsort[i].first;
                            downsample_prob.insert(make_pair(bs_idx[clcountsort[i].second], dsp));
                        }
                        break;
                    }
                }
            }
        }
    }
    else{
        fprintf(stderr, "%f total clade counts; no downsampling needed\n", allbc2);
        return 0;
    }
    
    // Compute downsample probabilities for missing genotype sites
    unordered_map<pair<bitset<NBITS>, bitset<NBITS>>, double> downsample_miss_prob;

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
    }
    
    if (downsample_prob.size() == 0){
        fprintf(stderr, "No clades were found with more than %d occurrences.\n", num);
        fprintf(stderr, "There is nothing to downsample.\n");
        return 0;
    }

    // ========================================================================
    // PASS 2: Filter VCF (sequential - output must be ordered)
    // ========================================================================
    
    fprintf(stderr, "Pass 2: Filtering VCF...\n");
    start_time = chrono::high_resolution_clock::now();
    
    bcf_reader = hts_open(vcf_file.c_str(), "r");
    bcf_header = bcf_hdr_read(bcf_reader);
    bcf1_t* bcf_record = bcf_init();
    
    // Add INFO headers
    bcf_hdr_append(bcf_header, "##INFO=<ID=DEMUX_SCORE,Number=1,Type=Float,Description=\"Clade rarity score\">");
    if (!gtf_file.empty()) {
        bcf_hdr_append(bcf_header, "##INFO=<ID=ANNOT_SCORE,Number=1,Type=Float,Description=\"Annotation weight (1.0=genic, 0.1=intergenic)\">");
    }
    if (cov_tbx) {
        bcf_hdr_append(bcf_header, "##INFO=<ID=COV_SCORE,Number=1,Type=Float,Description=\"Coverage score from bedgraph\">");
    }
    
    htsFile* outf = hts_open(outfile.c_str(), "wz");
    if (bcf_hdr_write(outf, bcf_header) < 0) {
        fprintf(stderr, "ERROR: Failed to write VCF header\n");
        exit(1);
    }

    // Random number generator
    random_device rd;
    mt19937 gen(seed >= 0 ? seed : rd());
    uniform_real_distribution<> rand_dist(0.0, 1.0);
    
    bitset<NBITS> alt, alt_flip, present;
    long nsnp = 0;
    int elim = 0;
    int kept = 0;
    
    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        
        if (bcf_record->n_allele == 2){ 
            bool pass;
            if (proc_bcf_record(bcf_record, bcf_header, num_samples,
                alt, alt_flip, present, pass)){
                
                double samp_prob = 1.0;
                int ac = alt.count();
                int afc = alt_flip.count();
                
                if (present.count() == num_samples){
                    if (ac < afc || (ac == afc && alt < alt_flip)){
                        if (downsample_prob.count(alt) > 0){
                            samp_prob = downsample_prob[alt];
                        }
                    }
                    else{
                        if (downsample_prob.count(alt_flip) > 0){
                            samp_prob = downsample_prob[alt_flip];
                        }
                    }
                }
                else{
                    pair<bitset<NBITS>, bitset<NBITS>> key;
                    if (ac < afc || (ac == afc && alt < alt_flip)){
                        key = make_pair(present, alt);
                    }
                    else{
                        key = make_pair(present, alt_flip);
                    }
                    if (downsample_miss_prob.count(key) > 0){
                        samp_prob = downsample_miss_prob[key];
                    }
                }

                bool print_rec = true;
                if (samp_prob < 1.0){
                    double r = rand_dist(gen);
                    if (r > samp_prob){
                        print_rec = false;
                    }
                }
                
                if (print_rec){
                    float demux_score = (samp_prob < 1.0) ? (float)samp_prob : 1.0f;
                    bcf_update_info_float(bcf_header, bcf_record, "DEMUX_SCORE", &demux_score, 1);
                    
                    if (!gtf_file.empty()) {
                        float annot_score = get_annot_score(
                            bcf_hdr_id2name(bcf_header, bcf_record->rid),
                            bcf_record->pos,
                            annotations);
                        bcf_update_info_float(bcf_header, bcf_record, "ANNOT_SCORE", &annot_score, 1);
                    }
                    
                    if (cov_tbx) {
                        float cov_score = get_coverage_score(
                            bcf_hdr_id2name(bcf_header, bcf_record->rid),
                            bcf_record->pos,
                            cov_tbx,
                            cov_fp);
                        bcf_update_info_float(bcf_header, bcf_record, "COV_SCORE", &cov_score, 1);
                    }
                    
                    bcf_write1(outf, bcf_header, bcf_record);
                    kept++;
                }
                else{
                    elim++;
                }
            }
            else{
                ++elim;
            }
        }
        ++nsnp;
        if (nsnp % 100000 == 0){
            fprintf(stderr, "  Processed %ld SNPs, kept %d, eliminated %d\r", nsnp, kept, elim);
        }
    }
    
    end_time = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    
    fprintf(stderr, "\nPass 2 complete: %ld SNPs in %ld seconds\n", nsnp, duration.count());
    fprintf(stderr, "Kept %d SNPs, eliminated %d SNPs\n", kept, elim);

    bcf_destroy(bcf_record);
    bcf_hdr_destroy(bcf_header);
    hts_close(bcf_reader);
    hts_close(outf);
    
    if (cov_tbx) {
        tbx_destroy(cov_tbx);
        hts_close(cov_fp);
    }

    return 0;
}
