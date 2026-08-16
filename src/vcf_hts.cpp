// =============================================================================
// vcf_hts.cpp
// Unified VCF/BCF, shared-memory panel, BAM counting, and HTS helpers.
// =============================================================================

#include <string>
#include <climits>
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
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <map>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <zlib.h>
#include <atomic>
#include <mutex>
#include <chrono>
#include <cmath>
#include <exception>
#include <limits>
#include <omp.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htswrapper/bc.h>
#include <htswrapper/bam.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "vcf_hts.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * ===== Contains functions relating to processing HTSlib-format files =====
 */

namespace {

class ParallelOperationStatus {
  public:
    ParallelOperationStatus() : ok_(true) {}

    bool ok() const { return ok_.load(std::memory_order_acquire); }

    void fail(const std::string& message) {
        bool expected = true;
        if (ok_.compare_exchange_strong(expected, false, std::memory_order_acq_rel)) {
            std::lock_guard<std::mutex> lock(message_mutex_);
            message_ = message;
        }
    }

    std::string message() const {
        std::lock_guard<std::mutex> lock(message_mutex_);
        return message_;
    }

  private:
    std::atomic<bool> ok_;
    mutable std::mutex message_mutex_;
    std::string message_;
};

static std::string format_worker_error(const char* operation, int thread_id, int tid = -1) {
    std::ostringstream out;
    out << operation << " failed in worker " << thread_id;
    if (tid >= 0) out << " for BAM target id " << tid;
    return out.str();
}

}  // namespace

size_t estimate_cellcounts_bytes(int n_samples) {
    if (n_samples <= 0) return 0;
    const size_t state_count = (size_t)n_samples * (size_t)GENOTYPE_STATES;
    if (state_count > std::numeric_limits<size_t>::max() / state_count) {
        return std::numeric_limits<size_t>::max();
    }
    const size_t pair_slots = state_count * state_count;
    if (pair_slots > (std::numeric_limits<size_t>::max() / sizeof(int64_t) - state_count) / 2) {
        return std::numeric_limits<size_t>::max();
    }
    return (2 * pair_slots + 2 * state_count) * sizeof(int64_t);
}

bool validate_identity_and_allocation_request(
    int n_samples,
    size_t* n_identity_states,
    size_t* bytes_per_cell,
    std::string* error_message) {

    auto fail = [&](const std::string& message) {
        if (error_message) *error_message = message;
        return false;
    };

    if (n_samples <= 0) {
        return fail("identity universe must contain at least one sample");
    }
    if (n_samples > MAX_INDIVIDUALS) {
        std::ostringstream out;
        out << "identity universe contains " << n_samples
            << " samples, exceeding the bitset-backed limit of " << MAX_INDIVIDUALS;
        return fail(out.str());
    }

    const size_t n = (size_t)n_samples;
    if (n > 1 && (n - 1) > std::numeric_limits<size_t>::max() / n) {
        return fail("identity-pair count overflows size_t");
    }
    const size_t pair_count = n * (n - 1) / 2;
    if (pair_count > std::numeric_limits<size_t>::max() - n) {
        return fail("identity-state count overflows size_t");
    }
    const size_t identity_states = n + pair_count;
    if (identity_states > (size_t)std::numeric_limits<int>::max()) {
        return fail("identity-state count exceeds the supported integer index range");
    }
    if (n_samples > MAX_COMBINATION_SAFE_INDIVIDUALS ||
        identity_states - 1 > (size_t)std::numeric_limits<short>::max()) {
        std::ostringstream out;
        out << "identity universe contains " << n_samples
            << " samples and " << identity_states
            << " singlet/doublet states, exceeding the signed-short range used by "
               "the shared haplotype-combination mapping; maximum supported sample count is "
            << MAX_COMBINATION_SAFE_INDIVIDUALS;
        return fail(out.str());
    }

    const size_t bytes = estimate_cellcounts_bytes(n_samples);
    if (bytes == std::numeric_limits<size_t>::max()) {
        return fail("dense CellCounts size overflows size_t");
    }
    if (bytes > MAX_CELLCOUNTS_BYTES_PER_CELL) {
        std::ostringstream out;
        out << "dense CellCounts request requires " << bytes
            << " bytes per cell, exceeding the supported "
            << MAX_CELLCOUNTS_BYTES_PER_CELL << "-byte safety limit";
        return fail(out.str());
    }

    if (n_identity_states) *n_identity_states = identity_states;
    if (bytes_per_cell) *bytes_per_cell = bytes;
    if (error_message) error_message->clear();
    return true;
}

const ReadFilterPolicy& default_production_read_filter() {
    static const ReadFilterPolicy policy(
        BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP,
        SupplementaryReadHandling::INCLUDE);
    return policy;
}

bool read_passes_filter(const bam1_t* record, const ReadFilterPolicy& policy) {
    if (record == nullptr) return false;
    uint16_t excluded = policy.excluded_flags;
    if (policy.supplementary == SupplementaryReadHandling::EXCLUDE) {
        excluded |= BAM_FSUPPLEMENTARY;
    }
    return (record->core.flag & excluded) == 0;
}

// ============================================================================
// VCF READING FUNCTIONS
// ============================================================================

void read_vcf_samples(string& filename, 
    vector<string>& samples){
    bcf_srs_t* sr = bcf_sr_init();
    if (!sr){
        fprintf(stderr, "Could not init VCF/BCF reader.\n");
        exit(1);
    }
    if (bcf_sr_add_reader(sr, filename.c_str()) < 0){
        fprintf(stderr, "ERROR: could not open VCF/BCF file %s\n", filename.c_str());
        bcf_sr_destroy(sr);
        exit(1);
    }
    bcf_hdr_t* bcf_header = bcf_sr_get_header(sr, 0);
    for (int i = 0; i < bcf_hdr_nsamples(bcf_header); ++i){
        samples.push_back(bcf_header->samples[i]);           
    }
    bcf_sr_destroy(sr);
}

int read_vcf_chrom(string& vcf_file,
    string& chrom,
    map<int, var>& snps,
    int min_vq,
    bool allow_missing){
    
    // First, check whether an index exists.
    htsFile* test = hts_open(vcf_file.c_str(), "r");
    if (test->format.format == vcf){
        tbx_t* idxptr = tbx_index_load(vcf_file.c_str());
        if (idxptr == NULL){
            fprintf(stderr, "Index not found for %s. Creating...\n", vcf_file.c_str());
            if (tbx_index_build(vcf_file.c_str(), 14, &tbx_conf_vcf) != 0){
                fprintf(stderr, "ERROR writing index for %s\n", vcf_file.c_str());
                exit(1);
            }
        }
        else{
            tbx_destroy(idxptr);
        }
    }
    else if (test->format.format == bcf){
        hts_idx_t* idxptr = bcf_index_load(vcf_file.c_str());
        if (idxptr == NULL){
            fprintf(stderr, "Index not found for %s. Creating...\n", vcf_file.c_str());
            if (bcf_index_build(vcf_file.c_str(), 14) != 0){
                fprintf(stderr, "ERROR: writing index for %s\n", vcf_file.c_str());
                exit(1);
            }
        }
        else{
            hts_idx_destroy(idxptr);
        }
    }
    hts_close(test);

    bcf_srs_t* sr = bcf_sr_init();
    if (!sr){
        fprintf(stderr, "Could not init VCF/BCF reader.\n");
        exit(1);
    }
    if (bcf_sr_set_regions(sr, chrom.c_str(), 0) < 0){
        fprintf(stderr, "ERROR: unable to set region %s\n", chrom.c_str());
        exit(1);
    }
    if (bcf_sr_add_reader(sr, vcf_file.c_str()) < 0){
        fprintf(stderr, "ERROR: could not open VCF/BCF file %s\n", vcf_file.c_str());
        bcf_sr_destroy(sr);
        exit(1);
    }
    int num_samples = bcf_hdr_nsamples(bcf_sr_get_header(sr, 0));

    long int nvar = 0;
    long int skipped_missing_gt = 0;
    long int skipped_malformed_gt = 0;
    set<int> bl; 
    
    bcf_hdr_t* bcf_header = bcf_sr_get_header(sr, 0);
    
    while (bcf_sr_next_line(sr)){
        if (bcf_sr_has_line(sr, 0)){
            bcf1_t* bcf_record = bcf_sr_get_line(sr, 0);

            int pos = bcf_record->pos;
            if (bl.find(pos) != bl.end()){
                continue;
            }
            if (snps.count(pos) > 0){
                fprintf(stderr, "WARNING: duplicate variants at site %s:%d\n", chrom.c_str(), pos+1);
                snps.erase(pos);
                bl.insert(pos);
            }
            if (bcf_record->n_allele == 2){ 
                bcf_unpack(bcf_record, BCF_UN_STR);

                bool pass = true;
                for (int i = 0; i < 2; ++i){
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
                else if (bcf_record->qual < min_vq){
                    pass = false;
                }
                if (pass){
                    var v;
                    v.ref = bcf_record->d.allele[0][0];
                    v.alt = bcf_record->d.allele[1][0];
                    v.vq = bcf_record->qual;

                    int32_t* gts = NULL;
                    int n_gts = 0;
                    int nmiss = 0;
                    int num_loaded = bcf_get_genotypes(bcf_header, bcf_record, &gts, &n_gts);
                    if (num_loaded <= 0){
                        ++skipped_missing_gt;
                        free(gts);
                        continue;
                    }
                    if (num_samples <= 0 || num_loaded % num_samples != 0 ||
                        num_loaded / num_samples < 2){
                        ++skipped_malformed_gt;
                        free(gts);
                        continue;
                    }
                    const int ploidy = num_loaded / num_samples;

                    // GQ is intentionally ignored. Historical and current panels
                    // use different GQ header types, and genotype quality is not
                    // required by the demultiplexing model. GT alone determines
                    // whether a donor genotype is available at this site.
                    for (int i = 0; i < num_samples; ++i){
                        int32_t* gtptr = gts + i*ploidy;
                        const bool missing_gt =
                            bcf_gt_is_missing(gtptr[0]) ||
                            gtptr[0] == bcf_int32_vector_end ||
                            bcf_gt_is_missing(gtptr[1]) ||
                            gtptr[1] == bcf_int32_vector_end;
                        if (missing_gt){
                            ++nmiss;
                            continue;
                        }
                        const int allele0 = bcf_gt_allele(gtptr[0]);
                        const int allele1 = bcf_gt_allele(gtptr[1]);
                        if (allele0 < 0 || allele0 > 1 || allele1 < 0 || allele1 > 1){
                            ++nmiss;
                            continue;
                        }
                        v.haps_covered.set(i);
                        if (allele0 == 1) v.haps1.set(i);
                        if (allele1 == 1) v.haps2.set(i);
                    }
                    free(gts);

                    if (allow_missing || nmiss == 0){
                        snps.insert(make_pair(pos, v));
                        ++nvar;
                    }
                }
            }
        }
    }
    
    bcf_sr_destroy(sr);
    if (skipped_missing_gt > 0 || skipped_malformed_gt > 0){
        fprintf(stderr,
            "WARNING: skipped %ld site(s) without usable GT and %ld site(s) with unsupported GT width on %s\n",
            skipped_missing_gt, skipped_malformed_gt, chrom.c_str());
    }
    return nvar;
}

void get_vcf_chroms(string& vcf_file, set<string>& chroms){
    bcf_srs_t* sr = bcf_sr_init();
    if (!sr){
        fprintf(stderr, "ERROR: could not init VCF/BCF reader\n");
        exit(1);
    }
    if (bcf_sr_add_reader(sr, vcf_file.c_str()) < 0){
        fprintf(stderr, "ERROR: could not open VCF/BCF file %s\n", vcf_file.c_str());
        bcf_sr_destroy(sr);
        exit(1);
    }
    bcf_hdr_t* bcf_header = bcf_sr_get_header(sr, 0);
    for (int i = 0; i < bcf_header->n[BCF_DT_CTG]; ++i){
        string chrom = bcf_hdr_id2name(bcf_header, i);
        chroms.insert(chrom);
    }
    bcf_sr_destroy(sr);
}

void get_bam_chroms(bam_reader& reader, set<string>& chroms){
    map<string, int> seq2tid = reader.get_seq2tid();
    for (map<string, int>::iterator it = seq2tid.begin(); it != seq2tid.end(); ++it){
        chroms.insert(it->first);
    }
}

bool get_bam_header_chroms_and_seq2tid(const string& bamfile,
    set<string>& chroms,
    map<string, int>& seq2tid,
    string* error_message){

    chroms.clear();
    seq2tid.clear();

    htsFile* fp = hts_open(bamfile.c_str(), "r");
    if (!fp){
        if (error_message) *error_message = "could not open BAM file for header read: " + bamfile;
        return false;
    }

    bam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr){
        if (error_message) *error_message = "could not read BAM header from: " + bamfile;
        hts_close(fp);
        return false;
    }

    bool ok = true;
    for (int tid = 0; tid < hdr->n_targets; ++tid){
        const char* name = hdr->target_name[tid];
        if (!name || name[0] == '\0'){
            ok = false;
            if (error_message) *error_message = "BAM header contains an invalid contig name";
            break;
        }
        string chrom(name);
        chroms.insert(chrom);
        seq2tid[chrom] = tid;
    }

    bam_hdr_destroy(hdr);
    hts_close(fp);
    if (ok && error_message) error_message->clear();
    return ok;
}

long int count_vcf_snps(string& vcf_file, set<string>& chroms_to_include, int min_vq){
    htsFile* bcf_reader = bcf_open(vcf_file.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR: could not open VCF/BCF file %s\n", vcf_file.c_str());
        exit(1);
    }
    bcf_hdr_t* bcf_header = bcf_hdr_read(bcf_reader);
    bcf1_t* bcf_record = bcf_init();
    
    long int nvar = 0;
    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        string chrom = bcf_hdr_id2name(bcf_header, bcf_record->rid);
        if (chroms_to_include.find(chrom) == chroms_to_include.end()){
            continue;
        }
        if (bcf_record->n_allele == 2){
            bcf_unpack(bcf_record, BCF_UN_STR);
            bool pass = true;
            for (int i = 0; i < 2; ++i){
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
            else if (bcf_record->qual < min_vq){
                pass = false;
            }
            if (pass){
                nvar++;
            }
        }
    }
    
    bcf_destroy(bcf_record);
    bcf_hdr_destroy(bcf_header);
    hts_close(bcf_reader);
    
    return nvar;
}

int read_vcf_chroms(string& vcf_file,
    set<string>& chroms_to_include,
    map<string, int>& seq2tid,
    map<int, map<int, var> >& snps,
    int min_vq,
    bool allow_missing){
    
    htsFile* bcf_reader = bcf_open(vcf_file.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR: could not open VCF/BCF file %s\n", vcf_file.c_str());
        exit(1);
    }
    bcf_hdr_t* bcf_header = bcf_hdr_read(bcf_reader);
    bcf1_t* bcf_record = bcf_init();
    int num_samples = bcf_hdr_nsamples(bcf_header);
    
    long int nvar = 0;
    long int skipped_missing_gt = 0;
    long int skipped_malformed_gt = 0;
    set<pair<int, int> > bl;
    
    int progress = 1000000;
    long int last_print = 0;
    
    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        string chrom = bcf_hdr_id2name(bcf_header, bcf_record->rid);
        
        // Progress indicator
        if (nvar - last_print >= progress){
            fprintf(stderr, "Loaded %ld SNPs\r", nvar);
            last_print = nvar;
        }
        
        if (chroms_to_include.find(chrom) == chroms_to_include.end()){
            continue;
        }
        if (seq2tid.find(chrom) == seq2tid.end()){
            continue;
        }
        
        int tid = seq2tid[chrom];
        int pos = bcf_record->pos;
        
        pair<int, int> key = make_pair(tid, pos);
        if (bl.find(key) != bl.end()){
            continue;
        }
        
        if (snps.count(tid) > 0 && snps[tid].count(pos) > 0){
            fprintf(stderr, "WARNING: duplicate variants at site %s:%d\n", chrom.c_str(), pos+1);
            snps[tid].erase(pos);
            bl.insert(key);
            continue;
        }
        
        if (bcf_record->n_allele == 2){
            bcf_unpack(bcf_record, BCF_UN_STR);
            
            bool pass = true;
            for (int i = 0; i < 2; ++i){
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
            else if (bcf_record->qual < min_vq){
                pass = false;
            }
            
            if (pass){
                if (snps.count(tid) == 0){
                    map<int, var> m;
                    snps.insert(make_pair(tid, m));
                }
                var v;
                v.ref = bcf_record->d.allele[0][0];
                v.alt = bcf_record->d.allele[1][0];
                v.vq = bcf_record->qual;

                int32_t* gts = NULL;
                int n_gts = 0;
                int nmiss = 0;
                int num_loaded = bcf_get_genotypes(bcf_header, bcf_record, &gts, &n_gts);
                if (num_loaded <= 0){
                    ++skipped_missing_gt;
                    free(gts);
                    continue;
                }
                if (num_samples <= 0 || num_loaded % num_samples != 0 ||
                    num_loaded / num_samples < 2){
                    ++skipped_malformed_gt;
                    free(gts);
                    continue;
                }
                const int ploidy = num_loaded / num_samples;

                // GQ is intentionally ignored. Historical and current panels
                // use different GQ header types, and genotype quality is not
                // required by the demultiplexing model. GT alone determines
                // whether a donor genotype is available at this site.
                for (int i = 0; i < num_samples; ++i){
                    int32_t* gtptr = gts + i*ploidy;
                    const bool missing_gt =
                        bcf_gt_is_missing(gtptr[0]) ||
                        gtptr[0] == bcf_int32_vector_end ||
                        bcf_gt_is_missing(gtptr[1]) ||
                        gtptr[1] == bcf_int32_vector_end;
                    if (missing_gt){
                        ++nmiss;
                        continue;
                    }
                    const int allele0 = bcf_gt_allele(gtptr[0]);
                    const int allele1 = bcf_gt_allele(gtptr[1]);
                    if (allele0 < 0 || allele0 > 1 || allele1 < 0 || allele1 > 1){
                        ++nmiss;
                        continue;
                    }
                    v.haps_covered.set(i);
                    if (allele0 == 1) v.haps1.set(i);
                    if (allele1 == 1) v.haps2.set(i);
                }
                free(gts);

                if (allow_missing || nmiss == 0){
                    snps[tid].insert(make_pair(pos, v));
                    ++nvar;
                }
            }
        }
    }
    
    bcf_destroy(bcf_record);
    bcf_hdr_destroy(bcf_header);
    hts_close(bcf_reader);
    if (skipped_missing_gt > 0 || skipped_malformed_gt > 0){
        fprintf(stderr,
            "WARNING: skipped %ld site(s) without usable GT and %ld site(s) with unsupported GT width while loading %s\n",
            skipped_missing_gt, skipped_malformed_gt, vcf_file.c_str());
    }
    
    return nvar;
}

// ============================================================================
// OPTIMIZED VCF READING FUNCTIONS
// ============================================================================

int read_vcf_chroms_optimized(string& vcf_file,
    set<string>& chroms_to_include,
    map<string, int>& seq2tid,
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_optimized,
    int min_vq,
    bool allow_missing){
    
    auto t1 = std::chrono::steady_clock::now();
    
    // First read into old format
    map<int, map<int, var> > snps_old;
    int nvar = read_vcf_chroms(vcf_file, chroms_to_include, seq2tid, snps_old, min_vq, allow_missing);
    if (nvar < 0){
        snpdat_optimized.clear();
        return -1;
    }
    
    auto t2 = std::chrono::steady_clock::now();
    auto read_secs = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    fprintf(stderr, "  VCF parsing took %ld seconds\n", read_secs);
    
    // Convert to optimized format
    convert_snpdat_to_optimized(snps_old, snpdat_optimized);
    
    auto t3 = std::chrono::steady_clock::now();
    auto convert_secs = std::chrono::duration_cast<std::chrono::seconds>(t3 - t2).count();
    fprintf(stderr, "  Conversion took %ld seconds\n", convert_secs);
    
    return nvar;
}

void convert_snpdat_to_optimized(
    map<int, map<int, var> >& snpdat_old,
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_new){
    
    snpdat_new.clear();
    
    for (auto& kv : snpdat_old){
        int tid = kv.first;
        ChromSNPs& cs = snpdat_new[tid];
        cs.snps.reserve(kv.second.size());
        
        for (auto& snp_kv : kv.second){
            cs.snps.push_back(SNPData(snp_kv.first, snp_kv.second));
        }
        
        // Should already be sorted since map iterates in order, but ensure it
        cs.sort_snps();
    }
}

void precompute_all_genotypes(
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_all,
    int n_samples){
    
    long count = 0;
    for (auto& kv : snpdat_all){
        for (auto& snp : kv.second.snps){
            snp.precompute_genotypes(n_samples);
            snp.precompute_targets(n_samples);
            count++;
        }
    }
    fprintf(stderr, "Precomputed genotypes and targets for %ld SNPs (%d samples)\n",
        count, n_samples);
}

void build_bin_indices(
    const robin_hood::unordered_map<int, ChromSNPs>& snpdat_all,
    const vector<int64_t>& chrom_lengths,
    robin_hood::unordered_map<int, ChromBinIndex>& bin_indices){
    
    bin_indices.clear();
    long hot_bins = 0, total_bins = 0;
    for (const auto& kv : snpdat_all){
        int tid = kv.first;
        int64_t len = (tid < (int)chrom_lengths.size()) ? chrom_lengths[tid] : 0;
        if (len == 0) continue;
        
        ChromBinIndex& idx = bin_indices[tid];
        idx.build(kv.second, len);
        for (int b = 0; b < idx.n_bins; ++b){
            if (idx.hot[b]) hot_bins++;
            total_bins++;
        }
    }
    fprintf(stderr, "Bin index: %ld/%ld bins hot (%.1f%% coverage)\n",
        hot_bins, total_bins, 100.0 * hot_bins / (total_bins > 0 ? total_bins : 1));
}

// ============================================================================
// CONDITIONAL MATCH FRACTION FUNCTIONS
// ============================================================================

void get_conditional_match_fracs_chrom(map<int, var>& snpdat,
    map<pair<int, int>, map<int, float> >& conditional_match_fracs,
    map<pair<int, int>, map<int, float> >& conditional_match_tots,
    int n_samples){
    
    for (map<int, var>::iterator s = snpdat.begin(); s != snpdat.end(); ++s){
        for (int i = 0; i < n_samples; ++i){
            if (s->second.haps_covered.test(i)){
                int nalt_i = 0;
                if (s->second.haps1.test(i)) nalt_i++;
                if (s->second.haps2.test(i)) nalt_i++;
                pair<int, int> key_i = make_pair(i, nalt_i);
                
                for (int j = 0; j < n_samples; ++j){
                    if (s->second.haps_covered.test(j)){
                        int nalt_j = 0;
                        if (s->second.haps1.test(j)) nalt_j++;
                        if (s->second.haps2.test(j)) nalt_j++;
                        
                        if (conditional_match_fracs.count(key_i) == 0){
                            map<int, float> m;
                            conditional_match_fracs.insert(make_pair(key_i, m));
                            conditional_match_tots.insert(make_pair(key_i, m));
                        }
                        if (conditional_match_fracs[key_i].count(j) == 0){
                            conditional_match_fracs[key_i].insert(make_pair(j, 0.0));
                            conditional_match_tots[key_i].insert(make_pair(j, 0.0));
                        }
                        conditional_match_fracs[key_i][j] += (float)nalt_j / 2.0;
                        conditional_match_tots[key_i][j] += 1.0;
                    }
                }
            }
        }
    }
}

void get_conditional_match_fracs_chrom_optimized(ChromSNPs& snpdat,
    map<pair<int, int>, map<int, float> >& conditional_match_fracs,
    map<pair<int, int>, map<int, float> >& conditional_match_tots,
    int n_samples){
    
    for (const auto& s : snpdat.snps){
        for (int i = 0; i < n_samples; ++i){
            if (s.data.haps_covered.test(i)){
                int nalt_i = 0;
                if (s.data.haps1.test(i)) nalt_i++;
                if (s.data.haps2.test(i)) nalt_i++;
                pair<int, int> key_i = make_pair(i, nalt_i);
                
                for (int j = 0; j < n_samples; ++j){
                    if (s.data.haps_covered.test(j)){
                        int nalt_j = 0;
                        if (s.data.haps1.test(j)) nalt_j++;
                        if (s.data.haps2.test(j)) nalt_j++;
                        
                        if (conditional_match_fracs.count(key_i) == 0){
                            map<int, float> m;
                            conditional_match_fracs.insert(make_pair(key_i, m));
                            conditional_match_tots.insert(make_pair(key_i, m));
                        }
                        if (conditional_match_fracs[key_i].count(j) == 0){
                            conditional_match_fracs[key_i].insert(make_pair(j, 0.0));
                            conditional_match_tots[key_i].insert(make_pair(j, 0.0));
                        }
                        conditional_match_fracs[key_i][j] += (float)nalt_j / 2.0;
                        conditional_match_tots[key_i][j] += 1.0;
                    }
                }
            }
        }
    }
}

void conditional_match_fracs_normalize(map<pair<int, int>, map<int, float> >& conditional_match_fracs,
    map<pair<int, int>, map<int, float> >& conditional_match_tots,
    int n_samples){
    
    for (auto& kv : conditional_match_fracs){
        for (auto& kv2 : kv.second){
            if (conditional_match_tots[kv.first][kv2.first] > 0){
                kv2.second /= conditional_match_tots[kv.first][kv2.first];
            }
        }
    }
}

// ============================================================================
// PARALLEL CONDITIONAL MATCH FRACTION COMPUTATION
// ============================================================================
//
// The key space for conditional match fractions is:
//   row = (individual i, genotype nalt_i)  where i in [0, n_samples), nalt in {0,1,2}
//   col = individual j                     where j in [0, n_samples)
//
// So the flat index is: row = i * 3 + nalt_i, col = j
// Total array size: (n_samples * 3) * n_samples
//
// Each thread accumulates into its own flat array, then we sum across threads.

void compute_conditional_match_fracs_parallel(
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_all,
    map<pair<int, int>, map<int, float> >& conditional_match_fracs,
    int n_samples,
    int n_threads){
    
    int n_rows = n_samples * 3;
    int n_cols = n_samples;
    size_t arr_size = (size_t)n_rows * n_cols;
    
    // Flatten the chromosome map into a vector for OpenMP indexing
    vector<ChromSNPs*> chrom_ptrs;
    chrom_ptrs.reserve(snpdat_all.size());
    for (auto& kv : snpdat_all){
        chrom_ptrs.push_back(&kv.second);
    }
    int n_chroms = (int)chrom_ptrs.size();
    
    fprintf(stderr, "Parallel condf: %d chromosomes, %d samples, %d threads\n",
        n_chroms, n_samples, n_threads);

    // CONDF only needs SNPData::geno. Some callers (notably the dual-panel
    // species path) have already precomputed targets on a copied/combined SNP
    // set, but not on the original species_snpdat that is later used here for
    // .species_condf. Indexing an empty geno vector segfaults. Do a cheap,
    // genotype-only preflight here so this function is safe for every caller
    // without allocating the very large per-SNP pair_targets arrays.
    long geno_precomputed = 0;
    for (auto& kv : snpdat_all){
        for (auto& snp : kv.second.snps){
            if ((int)snp.geno.size() != n_samples){
                snp.precompute_genotypes(n_samples);
                ++geno_precomputed;
            }
        }
    }
    if (geno_precomputed > 0){
        fprintf(stderr, "Parallel condf: genotype-only precompute for %ld SNPs\n",
            geno_precomputed);
    }
    
    // Global accumulators (sum of all threads)
    vector<double> global_fracs(arr_size, 0.0);
    vector<double> global_tots(arr_size, 0.0);
    
    omp_set_num_threads(n_threads);
    
    #pragma omp parallel
    {
        // Thread-local accumulators
        vector<double> local_fracs(arr_size, 0.0);
        vector<double> local_tots(arr_size, 0.0);
        
        #pragma omp for schedule(dynamic, 1)
        for (int c = 0; c < n_chroms; c++){
            ChromSNPs& chrom = *chrom_ptrs[c];
            
            for (const auto& s : chrom.snps){
                const int8_t* geno = s.geno.data();
                for (int i = 0; i < n_samples; ++i){
                    int8_t nalt_i = geno[i];
                    if (nalt_i < 0 || nalt_i >= 3) continue;
                    int row = i * 3 + nalt_i;
                    
                    for (int j = 0; j < n_samples; ++j){
                        int8_t nalt_j = geno[j];
                        if (nalt_j < 0 || nalt_j >= 3) continue;
                        
                        size_t idx = (size_t)row * n_cols + j;
                        local_fracs[idx] += (double)nalt_j / 2.0;
                        local_tots[idx] += 1.0;
                    }
                }
            }
        }
        
        // Merge thread-local into global (critical section)
        #pragma omp critical
        {
            for (size_t k = 0; k < arr_size; k++){
                global_fracs[k] += local_fracs[k];
                global_tots[k] += local_tots[k];
            }
        }
    }
    
    // Normalize and convert to map format
    conditional_match_fracs.clear();
    for (int i = 0; i < n_samples; ++i){
        for (int nalt = 0; nalt < 3; ++nalt){
            int row = i * 3 + nalt;
            pair<int, int> key_i = make_pair(i, nalt);
            
            bool has_any = false;
            for (int j = 0; j < n_samples; ++j){
                size_t idx = (size_t)row * n_cols + j;
                if (global_tots[idx] > 0){
                    has_any = true;
                    break;
                }
            }
            if (!has_any) continue;
            
            map<int, float>& frac_map = conditional_match_fracs[key_i];
            for (int j = 0; j < n_samples; ++j){
                size_t idx = (size_t)row * n_cols + j;
                if (global_tots[idx] > 0){
                    frac_map[j] = (float)(global_fracs[idx] / global_tots[idx]);
                }
            }
        }
    }
    
    fprintf(stderr, "Parallel condf complete: %lu entries\n", conditional_match_fracs.size());
}

ConditionalWeightStats compute_conditional_match_fracs_weighted(
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_all,
    const AcceptedSiteWeightMap& accepted_site_weights,
    map<pair<int, int>, map<int, float> >& conditional_match_fracs,
    int n_samples,
    int n_threads){

    ConditionalWeightStats stats;
    const int n_rows = n_samples * 3;
    const int n_cols = n_samples;
    const size_t arr_size = (size_t)n_rows * n_cols;

    vector<pair<int, ChromSNPs*> > chrom_ptrs;
    chrom_ptrs.reserve(snpdat_all.size());
    for (auto& kv : snpdat_all){
        chrom_ptrs.push_back(make_pair(kv.first, &kv.second));
        for (auto& snp : kv.second.snps){
            if ((int)snp.geno.size() != n_samples) snp.precompute_genotypes(n_samples);
        }
    }

    vector<double> global_fracs(arr_size, 0.0);
    vector<double> global_tots(arr_size, 0.0);
    uint64_t global_sites = 0;
    long double global_weight = 0.0L;

    omp_set_num_threads(n_threads);
    #pragma omp parallel
    {
        vector<double> local_fracs(arr_size, 0.0);
        vector<double> local_tots(arr_size, 0.0);
        uint64_t local_sites = 0;
        long double local_weight = 0.0L;

        #pragma omp for schedule(dynamic, 1)
        for (int c = 0; c < (int)chrom_ptrs.size(); ++c){
            const int tid = chrom_ptrs[c].first;
            ChromSNPs& chrom = *chrom_ptrs[c].second;
            for (const auto& snp : chrom.snps){
                auto wit = accepted_site_weights.find(accepted_site_weight_key(tid, snp.pos));
                if (wit == accepted_site_weights.end() || wit->second <= 0) continue;
                const double weight = (double)wit->second;
                ++local_sites;
                local_weight += (long double)weight;
                const int8_t* geno = snp.geno.data();
                for (int i = 0; i < n_samples; ++i){
                    const int8_t nalt_i = geno[i];
                    if (nalt_i < 0 || nalt_i >= 3) continue;
                    const int row = i * 3 + nalt_i;
                    for (int j = 0; j < n_samples; ++j){
                        const int8_t nalt_j = geno[j];
                        if (nalt_j < 0 || nalt_j >= 3) continue;
                        const size_t idx = (size_t)row * n_cols + j;
                        local_fracs[idx] += weight * ((double)nalt_j / 2.0);
                        local_tots[idx] += weight;
                    }
                }
            }
        }
        #pragma omp critical
        {
            for (size_t k = 0; k < arr_size; ++k){
                global_fracs[k] += local_fracs[k];
                global_tots[k] += local_tots[k];
            }
            global_sites += local_sites;
            global_weight += local_weight;
        }
    }

    conditional_match_fracs.clear();
    for (int i = 0; i < n_samples; ++i){
        for (int nalt = 0; nalt < 3; ++nalt){
            const int row = i * 3 + nalt;
            const pair<int, int> key_i = make_pair(i, nalt);
            for (int j = 0; j < n_samples; ++j){
                const size_t idx = (size_t)row * n_cols + j;
                if (global_tots[idx] > 0.0){
                    conditional_match_fracs[key_i][j] =
                        (float)(global_fracs[idx] / global_tots[idx]);
                }
            }
        }
    }

    stats.observed_sites = global_sites;
    stats.accepted_weight = (double)(global_weight / (long double)FIXED_POINT_SCALE);
    fprintf(stderr,
        "Accepted-observation-weighted condf complete: %llu sites, %.6f accepted weight, %lu row entries\n",
        (unsigned long long)stats.observed_sites, stats.accepted_weight,
        conditional_match_fracs.size());
    return stats;
}

// ============================================================================
// ORIGINAL BAM PROCESSING FUNCTIONS
// ============================================================================

void process_bam_record(bam_reader& reader,
    int snppos,
    var& vardat,
    map<int, robin_hood::unordered_map<unsigned long, 
        pair<float, float> > >& varcounts_site,
    bool has_bc_list,
    set<unsigned long>& bcs_valid){

    if (!reader.unmapped() && !reader.secondary() && 
        !reader.dup() && reader.has_cb_z){
        
        bc bc_bits;
        str2bc(reader.cb_z, bc_bits);
        unsigned long bc_key = bc_bits.to_ulong();
        
        if (!has_bc_list || bcs_valid.find(bc_key) != bcs_valid.end()){
            int tid = reader.tid();
            
            float prob_corr = 1.0 - pow(10, -(float)reader.mapq/10.0);
            
            if (varcounts_site.count(snppos) == 0){
                robin_hood::unordered_map<unsigned long, pair<float, float> > m;
                varcounts_site.insert(make_pair(snppos, m));
            }
            if (varcounts_site[snppos].count(bc_key) == 0){
                varcounts_site[snppos].emplace(bc_key, make_pair(0.0f, 0.0f));
            }
            
            // Note: get_base_at expects 1-based position
            char allele = reader.get_base_at(snppos + 1);
            
            if (allele != 'N' && allele != '-'){
                if (allele == vardat.ref){
                    varcounts_site[snppos][bc_key].first += prob_corr;
                }
                else if (allele == vardat.alt){
                    varcounts_site[snppos][bc_key].second += prob_corr;
                }
            }
        }
    }
}

void dump_vcs_counts(robin_hood::unordered_map<unsigned long, pair<float, float> >& varcounts_site,
    robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > >& indv_allelecounts,
    var& snpdat,
    int n_samples){
    
    for (auto& vcs : varcounts_site){
        if (indv_allelecounts.count(vcs.first) == 0){
            map<pair<int, int>, map<pair<int, int>, pair<float, float> > > m;
            indv_allelecounts.emplace(vcs.first, m);
            
            for (int i = 0; i < n_samples; ++i){
                map<pair<int, int>, pair<float, float> > m2;
                for (int j = 0; j < 3; ++j){
                    pair<int, int> key = make_pair(i, j);
                    indv_allelecounts[vcs.first].insert(make_pair(key, m2));
                }
            }
        } 
        
        if (vcs.second.first + vcs.second.second > 0){
            for (int i = 0; i < n_samples; ++i){
                int n_alt_chroms = 0;
                if (snpdat.haps_covered.test(i)){
                    if (snpdat.haps1.test(i)){
                        n_alt_chroms++;
                    }
                    if (snpdat.haps2.test(i)){
                        n_alt_chroms++;
                    }
                    
                    pair<int, int> key = make_pair(i, n_alt_chroms);
                    
                    pair<int, int> nullkey = make_pair(-1, -1);
                    if (indv_allelecounts[vcs.first][key].count(nullkey) == 0){
                        indv_allelecounts[vcs.first][key].insert(make_pair(nullkey, 
                            make_pair(0.0,0.0)));    
                    }
                    indv_allelecounts[vcs.first][key][nullkey].first += vcs.second.first;
                    indv_allelecounts[vcs.first][key][nullkey].second += vcs.second.second; 
                    
                    for (int j = i + 1; j < n_samples; ++j){
                        if (snpdat.haps_covered.test(j)){
                            int n_alt_chroms_j = 0;
                            if (snpdat.haps1.test(j)){
                                n_alt_chroms_j++;
                            }
                            if (snpdat.haps2.test(j)){
                                n_alt_chroms_j++;
                            }
                            pair<int, int> key_j = make_pair(j, n_alt_chroms_j);
                            
                            if (key.first > key_j.first){
                                pair<int, int> tmp = key;
                                key = key_j;
                                key_j = tmp;
                            }

                            if (indv_allelecounts[vcs.first][key].count(key_j) == 0){
                                indv_allelecounts[vcs.first][key].insert(make_pair(key_j, 
                                    make_pair(0.0,0.0)));
                            }                       
                            indv_allelecounts[vcs.first][key][key_j].first += 
                                vcs.second.first;
                            indv_allelecounts[vcs.first][key][key_j].second += 
                                vcs.second.second;
                        }
                    }       
                }
            }
        }
    }
}

// ============================================================================
// PARALLEL BAM PROCESSING FUNCTIONS
// ============================================================================

/**
 * Resolve one 0-based reference position against a BAM record using strict
 * half-open intervals. Insertions and soft clips consume query only and can
 * never satisfy a reference coordinate; hard clips and pads consume neither.
 */
ReferenceCoordinateResult query_reference_coordinate(const bam1_t* record, int pos){
    if (record == nullptr){
        return ReferenceCoordinateResult(ReferenceCoordinateState::MALFORMED_CIGAR);
    }

    const int64_t alignment_start = record->core.pos;
    const int64_t alignment_end = bam_endpos(const_cast<bam1_t*>(record));
    if (pos < alignment_start || pos >= alignment_end){
        return ReferenceCoordinateResult(ReferenceCoordinateState::OUTSIDE_ALIGNMENT);
    }

    uint32_t* cigar = bam_get_cigar(const_cast<bam1_t*>(record));
    uint8_t* seq = bam_get_seq(const_cast<bam1_t*>(record));
    int64_t ref_pos = alignment_start;
    int64_t query_pos = 0;

    for (uint32_t i = 0; i < record->core.n_cigar; ++i){
        const int op = bam_cigar_op(cigar[i]);
        const int64_t len = bam_cigar_oplen(cigar[i]);
        if (len < 0){
            return ReferenceCoordinateResult(ReferenceCoordinateState::MALFORMED_CIGAR);
        }

        switch (op){
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF: {
                const int64_t block_end = ref_pos + len;
                if ((int64_t)pos >= ref_pos && (int64_t)pos < block_end){
                    const int64_t query_index = query_pos + ((int64_t)pos - ref_pos);
                    if (query_index < 0 || query_index >= record->core.l_qseq){
                        return ReferenceCoordinateResult(
                            ReferenceCoordinateState::NO_QUERY_BASE, 'N', -1);
                    }
                    const int base_code = bam_seqi(seq, (int)query_index);
                    const char base = seq_nt16_str[base_code];
                    if (base != 'A' && base != 'C' && base != 'G' && base != 'T'){
                        return ReferenceCoordinateResult(
                            ReferenceCoordinateState::NO_QUERY_BASE, 'N', (int)query_index);
                    }
                    return ReferenceCoordinateResult(
                        ReferenceCoordinateState::BASE, base, (int)query_index);
                }
                ref_pos = block_end;
                query_pos += len;
                break;
            }
            case BAM_CDEL: {
                const int64_t block_end = ref_pos + len;
                if ((int64_t)pos >= ref_pos && (int64_t)pos < block_end){
                    return ReferenceCoordinateResult(ReferenceCoordinateState::DELETION, '-', -1);
                }
                ref_pos = block_end;
                break;
            }
            case BAM_CREF_SKIP: {
                const int64_t block_end = ref_pos + len;
                if ((int64_t)pos >= ref_pos && (int64_t)pos < block_end){
                    return ReferenceCoordinateResult(
                        ReferenceCoordinateState::REFERENCE_SKIP, '-', -1);
                }
                ref_pos = block_end;
                break;
            }
            case BAM_CINS:
            case BAM_CSOFT_CLIP:
                query_pos += len;
                break;
            case BAM_CHARD_CLIP:
            case BAM_CPAD:
                break;
            default:
                return ReferenceCoordinateResult(ReferenceCoordinateState::MALFORMED_CIGAR);
        }

        if (query_pos < 0 || query_pos > record->core.l_qseq){
            return ReferenceCoordinateResult(ReferenceCoordinateState::MALFORMED_CIGAR);
        }
    }

    return ReferenceCoordinateResult(ReferenceCoordinateState::NO_QUERY_BASE);
}

char get_base_at_pos(const bam1_t* record, int pos){
    const ReferenceCoordinateResult result = query_reference_coordinate(record, pos);
    if (result.state == ReferenceCoordinateState::BASE) return result.base;
    if (result.state == ReferenceCoordinateState::DELETION ||
        result.state == ReferenceCoordinateState::REFERENCE_SKIP) return '-';
    return 'N';
}


static inline int64_t apply_species_target_weight(
    int64_t value,
    int32_t weight_scaled){

    if (value == 0 || weight_scaled == 0) return 0;
    const int64_t half = FIXED_POINT_SCALE / 2;
    return (value * (int64_t)weight_scaled + half) / FIXED_POINT_SCALE;
}

static inline bool accumulate_species_native_targets(
    CellCounts& counts,
    const NativeSpeciesChromTargets& chrom_targets,
    size_t snp_index,
    int64_t ref_add,
    int64_t alt_add){

    if (snp_index >= chrom_targets.site_offsets.size()) return false;
    const uint64_t offset = chrom_targets.site_offsets[snp_index];
    if (offset == UINT64_MAX) return false;
    if ((uint64_t)offset + chrom_targets.weights_per_site >
        (uint64_t)chrom_targets.weights.size()) return false;

    const int n_species = chrom_targets.n_species;
    const int state_count = n_species * GENOTYPE_STATES;
    const int32_t* weights = chrom_targets.weights.data() + offset;
    size_t cursor = 0;

    for (int sp = 0; sp < n_species; ++sp){
        for (int g = 0; g < GENOTYPE_STATES; ++g, ++cursor){
            const int32_t weight = weights[cursor];
            if (weight == 0) continue;
            const int total_idx = sp * GENOTYPE_STATES + g;
            counts.total_ref[total_idx] +=
                apply_species_target_weight(ref_add, weight);
            counts.total_alt[total_idx] +=
                apply_species_target_weight(alt_add, weight);
        }
    }

    for (int a = 0; a < n_species; ++a){
        for (int b = a + 1; b < n_species; ++b){
            for (int ga = 0; ga < GENOTYPE_STATES; ++ga){
                const int idx_a = a * GENOTYPE_STATES + ga;
                for (int gb = 0; gb < GENOTYPE_STATES; ++gb, ++cursor){
                    const int32_t weight = weights[cursor];
                    if (weight == 0) continue;
                    const int idx_b = b * GENOTYPE_STATES + gb;
                    const size_t pair_idx =
                        (size_t)idx_a * (size_t)state_count + (size_t)idx_b;
                    counts.ref_counts[pair_idx] +=
                        apply_species_target_weight(ref_add, weight);
                    counts.alt_counts[pair_idx] +=
                        apply_species_target_weight(alt_add, weight);
                }
            }
        }
    }
    return true;
}

bool count_alleles_parallel(
    const string& bamfile,
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_all,
    robin_hood::unordered_map<unsigned long, AlignedCellCounts>& cell_counts,
    const set<unsigned long>& valid_barcodes,
    int n_samples,
    int n_threads,
    int htslib_threads,
    bool dump_pileup,
    const string& pileup_prefix,
    AcceptedSiteWeightMap* accepted_site_weights,
    const NativeSpeciesTargetTable* species_native_targets,
    robin_hood::unordered_map<unsigned long, AlignedCellCounts>* species_native_counts,
    int species_native_n_samples){
    size_t n_identity_states = 0;
    size_t bytes_per_cell = 0;
    string request_error;
    if (!validate_identity_and_allocation_request(
            n_samples, &n_identity_states, &bytes_per_cell, &request_error)){
        fprintf(stderr, "ERROR: invalid individual identity universe: %s\n", request_error.c_str());
        return false;
    }
    if (species_native_n_samples > 0 && !validate_identity_and_allocation_request(
            species_native_n_samples, nullptr, nullptr, &request_error)){
        fprintf(stderr, "ERROR: invalid native-species identity universe: %s\n", request_error.c_str());
        return false;
    }
    if (n_threads < 1 || htslib_threads < 1){
        fprintf(stderr, "ERROR: thread counts must be positive\n");
        return false;
    }
    
    bool has_bc_list = !valid_barcodes.empty();
    const bool collect_species_native =
        species_native_targets != nullptr && species_native_counts != nullptr &&
        species_native_n_samples > 0;
    
    // Pre-allocate count structure for known barcodes after validating the
    // per-cell dense allocation and the total multiplication.
    if (has_bc_list){
        if (bytes_per_cell > 0 && valid_barcodes.size() >
            std::numeric_limits<size_t>::max() / bytes_per_cell){
            fprintf(stderr, "ERROR: projected CellCounts allocation overflows size_t\n");
            return false;
        }
        fprintf(stderr, "Pre-allocating counts for %lu cells (%lu bytes/cell)...\n",
            valid_barcodes.size(), (unsigned long)bytes_per_cell);
        try {
            for (unsigned long bc : valid_barcodes){
                cell_counts.emplace(std::piecewise_construct,
                    std::forward_as_tuple(bc),
                    std::forward_as_tuple(n_samples));
                if (collect_species_native){
                    species_native_counts->emplace(std::piecewise_construct,
                        std::forward_as_tuple(bc),
                        std::forward_as_tuple(species_native_n_samples));
                }
            }
        }
        catch (const std::exception& e){
            fprintf(stderr, "ERROR: CellCounts pre-allocation failed: %s\n", e.what());
            return false;
        }
    }
    
    // Get total number of chromosomes and their lengths from BAM header
    // Also get read counts per chromosome from BAM index
    htsFile* bam_tmp = hts_open(bamfile.c_str(), "r");
    if (!bam_tmp){
        fprintf(stderr, "ERROR: Could not open BAM file to get header: %s\n", bamfile.c_str());
        return false;
    }
    bam_hdr_t* hdr_tmp = sam_hdr_read(bam_tmp);
    if (!hdr_tmp){
        fprintf(stderr, "ERROR: Could not read BAM header: %s\n", bamfile.c_str());
        hts_close(bam_tmp);
        return false;
    }
    hts_idx_t* idx_tmp = sam_index_load(bam_tmp, bamfile.c_str());
    if (!idx_tmp){
        fprintf(stderr, "ERROR: Could not load required BAM index: %s\n", bamfile.c_str());
        bam_hdr_destroy(hdr_tmp);
        hts_close(bam_tmp);
        return false;
    }
    int n_chroms = hdr_tmp->n_targets;
    for (const auto& kv : snpdat_all){
        if (kv.first < 0 || kv.first >= n_chroms){
            fprintf(stderr, "ERROR: SNP panel references invalid BAM target id %d\n", kv.first);
            hts_idx_destroy(idx_tmp);
            bam_hdr_destroy(hdr_tmp);
            hts_close(bam_tmp);
            return false;
        }
    }

    // --dump_pileup: emit the per-SNP genotype sidecar for the variant-consistency
    // metric (interindividual panel only).  geno[] is already populated by
    // precompute_all_genotypes(), and hdr_tmp is still valid here (destroyed below).
    // Columns: tid  chrom  pos  ref  alt  geno_0 .. geno_{n_samples-1}  (0/1/2/-1).
    // The (tid,pos) pair is the producer's SNP join key; ref/alt are the allele
    // bases (informational; the metric is allele-orientation based).
    if (dump_pileup){
        string sites_path = pileup_prefix + ".pileup_sites.tsv.gz";
        gzFile sf = gzopen(sites_path.c_str(), "w");
        if (!sf){
            fprintf(stderr, "ERROR: could not open %s for writing\n", sites_path.c_str());
            hts_idx_destroy(idx_tmp);
            bam_hdr_destroy(hdr_tmp);
            hts_close(bam_tmp);
            return false;
        } else {
            long n_sites_written = 0;
            for (auto& kv : snpdat_all){
                int tid_s = kv.first;
                const char* cname = (tid_s >= 0 && tid_s < n_chroms) ?
                    hdr_tmp->target_name[tid_s] : ".";
                for (auto& snp : kv.second.snps){
                    if (snp.panel_id != 0) continue;
                    gzprintf(sf, "%d\t%s\t%d\t%c\t%c", tid_s, cname, snp.pos,
                        snp.data.ref, snp.data.alt);
                    for (int s = 0; s < n_samples; s++){
                        gzprintf(sf, "\t%d", (int)snp.geno[s]);
                    }
                    gzprintf(sf, "\n");
                    n_sites_written++;
                }
            }
            if (gzclose(sf) != Z_OK){
                fprintf(stderr, "ERROR: failed while closing %s\n", sites_path.c_str());
                hts_idx_destroy(idx_tmp);
                bam_hdr_destroy(hdr_tmp);
                hts_close(bam_tmp);
                return false;
            }
            fprintf(stderr, "Wrote %ld pileup sites to %s\n", n_sites_written, sites_path.c_str());
        }
    }

    // Get read counts per chromosome from index
    vector<uint64_t> chrom_read_counts(n_chroms, 0);
    vector<int64_t> chrom_lengths(n_chroms);
    vector<string> chrom_names(n_chroms);
    for (int i = 0; i < n_chroms; ++i){
        chrom_names[i] = hdr_tmp->target_name[i] ? hdr_tmp->target_name[i] : std::to_string(i);
    }
    
    const int n_index_targets = hts_idx_nseq(idx_tmp);
    if (n_index_targets < 0){
        fprintf(stderr, "ERROR: Could not determine BAM-index target count\n");
        hts_idx_destroy(idx_tmp);
        bam_hdr_destroy(hdr_tmp);
        hts_close(bam_tmp);
        return false;
    }
    int n_missing_index_stats = 0;
    int first_missing_index_stat = -1;
    for (int i = 0; i < n_chroms; i++){
        uint64_t mapped = 0, unmapped = 0;
        if (i >= n_index_targets ||
            hts_idx_get_stat(idx_tmp, i, &mapped, &unmapped) < 0){
            // A valid BAI/CSI may omit the metadata bin for a target with no
            // alignment records. These counts are used only to estimate work
            // unit size and ordering; the iterator still queries the complete
            // target below. Use an unchunked zero-read scheduling estimate.
            mapped = 0;
            unmapped = 0;
            if (first_missing_index_stat < 0) first_missing_index_stat = i;
            n_missing_index_stats++;
        }
        chrom_read_counts[i] = mapped;
        chrom_lengths[i] = hdr_tmp->target_len[i];
    }
    if (n_missing_index_stats > 0){
        const char* first_name =
            (first_missing_index_stat >= 0 && first_missing_index_stat < n_chroms &&
             hdr_tmp->target_name[first_missing_index_stat])
                ? hdr_tmp->target_name[first_missing_index_stat] : ".";
        fprintf(stderr,
            "WARNING: BAM-index mapped/unmapped statistics unavailable for %d of %d header targets "
            "(first target %d: %s); using zero only for work-unit scheduling estimates\n",
            n_missing_index_stats, n_chroms, first_missing_index_stat, first_name);
    }
    hts_idx_destroy(idx_tmp);
    bam_hdr_destroy(hdr_tmp);
    hts_close(bam_tmp);
    
    // Build bin index for read skipping (Change 3)
    robin_hood::unordered_map<int, ChromBinIndex> bin_indices;
    build_bin_indices(snpdat_all, chrom_lengths, bin_indices);
    
    // Work unit: represents a region to process (whole chromosome or chunk)
    struct WorkUnit {
        int tid;
        int start_pos;      // BAM region start (0 for whole chrom)
        int end_pos;        // BAM region end (INT_MAX for whole chrom)
        size_t snp_start;   // Index into ChromSNPs vector
        size_t snp_end;     // Index into ChromSNPs vector (exclusive)
        bool has_snps;
        uint64_t est_reads; // Estimated reads for sorting
    };
    
    // Build work units - chunk based on SNP density OR read density
    vector<WorkUnit> work_units;
    long total_snps = 0;
    int chroms_with_snps = 0;
    int chroms_chunked_by_snp = 0;
    int chroms_chunked_by_reads = 0;
    
    // Thresholds for chunking
    const size_t CHUNK_SNP_THRESHOLD = 100000;       // Chunk if >100k SNPs
    const uint64_t CHUNK_READ_THRESHOLD = 10000000;  // Chunk if >10M iterator records
    
    for (int tid = 0; tid < n_chroms; tid++){
        auto it = snpdat_all.find(tid);
        int64_t chrom_len = chrom_lengths[tid];
        uint64_t chrom_reads = chrom_read_counts[tid];
        
        if (it == snpdat_all.end() || it->second.empty()){
            // Targets with no active panel SNPs cannot contribute scientific
            // evidence. Do not create iterators merely to count records for
            // progress reporting; that adds unnecessary I/O and failure modes.
            continue;
        }
        else{
            ChromSNPs& chrom_snps = it->second;
            size_t n_snps = chrom_snps.snps.size();
            total_snps += n_snps;
            chroms_with_snps++;
            
            // Decide chunking strategy: SNP density or read density, whichever is higher
            size_t chunks_by_snp = (n_snps > CHUNK_SNP_THRESHOLD) ? 
                (n_snps + CHUNK_SNP_THRESHOLD - 1) / CHUNK_SNP_THRESHOLD : 1;
            size_t chunks_by_reads = (chrom_reads > CHUNK_READ_THRESHOLD) ?
                (chrom_reads + CHUNK_READ_THRESHOLD - 1) / CHUNK_READ_THRESHOLD : 1;
            
            // Cap read-based chunks
            chunks_by_reads = std::min(chunks_by_reads, (size_t)20);
            
            if (chunks_by_snp >= chunks_by_reads && chunks_by_snp > 1){
                // Chunk by SNP count
                size_t snps_per_chunk = (n_snps + chunks_by_snp - 1) / chunks_by_snp;
                uint64_t reads_per_chunk = chrom_reads / chunks_by_snp;
                
                for (size_t c = 0; c < chunks_by_snp; c++){
                    size_t snp_start = c * snps_per_chunk;
                    size_t snp_end = std::min(snp_start + snps_per_chunk, n_snps);
                    
                    int start_pos = (snp_start == 0) ? 0 : chrom_snps.snps[snp_start].pos;
                    int end_pos = (snp_end >= n_snps) ? INT_MAX : chrom_snps.snps[snp_end - 1].pos + 1000;
                    
                    work_units.push_back({tid, start_pos, end_pos, snp_start, snp_end, true, reads_per_chunk});
                }
                chroms_chunked_by_snp++;
            }
            else if (chunks_by_reads > 1){
                // Chunk by read density - split position space
                int64_t chunk_size = (chrom_len + chunks_by_reads - 1) / chunks_by_reads;
                uint64_t reads_per_chunk = chrom_reads / chunks_by_reads;
                
                size_t snp_idx = 0;
                for (size_t c = 0; c < chunks_by_reads; c++){
                    int start_pos = c * chunk_size;
                    int end_pos = (c == chunks_by_reads - 1) ? INT_MAX : (int)((c + 1) * chunk_size);
                    
                    // Find SNPs in this position range
                    size_t snp_start = snp_idx;
                    while (snp_idx < n_snps && chrom_snps.snps[snp_idx].pos < (c + 1) * chunk_size){
                        snp_idx++;
                    }
                    size_t snp_end = snp_idx;
                    if (snp_start == snp_end) continue;
                    
                    work_units.push_back({tid, start_pos, end_pos, snp_start, snp_end, true, reads_per_chunk});
                }
                chroms_chunked_by_reads++;
            }
            else{
                // No chunking needed - single unit
                work_units.push_back({tid, 0, INT_MAX, 0, n_snps, true, chrom_reads});
            }
        }
    }
    
    // Sort work units by estimated reads (highest first for best load balancing)
    std::sort(work_units.begin(), work_units.end(),
              [](const WorkUnit& a, const WorkUnit& b){
                  return a.est_reads > b.est_reads;
              });
    
    fprintf(stderr, "BAM header has %d targets; processing %d SNP-bearing targets (%ld total SNPs) using %d threads...\n",
        n_chroms, chroms_with_snps, total_snps, n_threads);
    fprintf(stderr, "  Split into %lu work units (%d by SNP density, %d by read density)\n", 
        work_units.size(), chroms_chunked_by_snp, chroms_chunked_by_reads);
    if (work_units.size() > 0){
        const WorkUnit& largest = work_units[0];
        fprintf(stderr, "  Largest work unit: %lu SNPs, ~%luM iterator records\n", 
            largest.snp_end - largest.snp_start, largest.est_reads / 1000000);
    }
    
    // Progress tracking
    atomic<long> snps_processed(0);
    atomic<int> units_done(0);
    atomic<long> reads_processed(0);
    
    // Per-thread cell counts for memory-efficient accumulation
    // With fixed-point int64 arithmetic, merge order doesn't matter (integer addition is associative)
    omp_set_num_threads(n_threads);
    vector<robin_hood::unordered_map<unsigned long, CellCounts>> thread_counts(n_threads);
    vector<robin_hood::unordered_map<unsigned long, CellCounts>> thread_species_native_counts(
        collect_species_native ? n_threads : 0);

    // --dump_pileup: per-thread per-(cell,SNP) allele evidence (interindividual
    // only).  Inner key packs (tid<<32 | pos); value is (ref_scaled, alt_scaled).
    // Empty and untouched unless dump_pileup is set.
    vector<robin_hood::unordered_map<unsigned long,
        robin_hood::unordered_map<int64_t, std::pair<int64_t, int64_t> > > > thread_pileup(n_threads);
    vector<AcceptedSiteWeightMap> thread_site_weights(n_threads);
    ParallelOperationStatus operation_status;
    std::atomic<bool> hts_thread_warning_emitted(false);
    
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        auto& local_counts = thread_counts[thread_id];
        robin_hood::unordered_map<unsigned long, CellCounts>* local_species_native =
            collect_species_native ? &thread_species_native_counts[thread_id] : nullptr;
        auto& local_site_weights = thread_site_weights[thread_id];
        
        // Each thread gets its own BAM reader. All workers still encounter the
        // same OpenMP work-sharing construct; a failed worker sets shared status.
        htsFile* bam_fp = hts_open(bamfile.c_str(), "r");
        bam_hdr_t* header = nullptr;
        hts_idx_t* idx = nullptr;
        bam1_t* record = nullptr;
        if (!bam_fp){
            operation_status.fail(format_worker_error("BAM open", thread_id));
        }
        else{
            if (htslib_threads > 1 && hts_set_threads(bam_fp, htslib_threads) < 0){
                bool expected = false;
                if (hts_thread_warning_emitted.compare_exchange_strong(expected, true)){
                    fprintf(stderr,
                        "WARNING: HTSlib helper-thread setup failed; continuing with synchronous BAM I/O\n");
                }
            }
            header = sam_hdr_read(bam_fp);
            if (!header){
                operation_status.fail(format_worker_error("BAM header read", thread_id));
            }
            idx = sam_index_load(bam_fp, bamfile.c_str());
            if (!idx){
                operation_status.fail(format_worker_error("BAM index load", thread_id));
            }
            record = bam_init1();
            if (!record){
                operation_status.fail(format_worker_error("BAM record allocation", thread_id));
            }
        }

        // Process work units with dynamic scheduling.
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < work_units.size(); i++){
            if (!operation_status.ok() || !bam_fp || !header || !idx || !record) continue;
            WorkUnit& wu = work_units[i];
            int tid = wu.tid;
            if (tid < 0 || tid >= header->n_targets){
                operation_status.fail(format_worker_error("invalid contig", thread_id, tid));
                continue;
            }
                    
                    // Query BAM for this region
            hts_itr_t* iter = sam_itr_queryi(idx, tid, wu.start_pos, wu.end_pos);
            if (!iter){
                operation_status.fail(format_worker_error("iterator creation", thread_id, tid));
                continue;
            }
                    
                    long local_snps = 0;
                    long local_reads = 0;
                    long local_all_reads = 0;
                    
                    if (!wu.has_snps){
                        // No SNPs on this chromosome - just count reads
                        int iterator_result = 0;
                        while ((iterator_result = sam_itr_next(bam_fp, iter, record)) >= 0){
                            if (!read_passes_filter(record, default_production_read_filter())){
                                continue;
                            }
                            local_all_reads++;
                        }
                        if (iterator_result < -1){
                            operation_status.fail(format_worker_error("iterator read", thread_id, tid));
                        }
                        hts_itr_destroy(iter);
                        reads_processed += local_all_reads;
                        
                        int done = ++units_done;
                        if (done % 100 == 0 || done == (int)work_units.size()){
                            fprintf(stderr, "\rProgress: %d/%lu units, %ld/%ld SNPs, %ld iterator records",
                                done, work_units.size(), snps_processed.load(), total_snps,
                                reads_processed.load());
                        }
                        continue;
                    }
                    
                    // Has SNPs - do full processing for this chunk
                    auto snp_it = snpdat_all.find(tid);
                    ChromSNPs& chrom_snps = snp_it->second;
                    const NativeSpeciesChromTargets* native_chrom_targets = nullptr;
                    if (collect_species_native){
                        auto native_it = species_native_targets->find(tid);
                        if (native_it != species_native_targets->end()){
                            native_chrom_targets = &native_it->second;
                        }
                    }
                    
                    // Get iterators for just our chunk of SNPs
                    auto snp_iter = chrom_snps.snps.begin() + wu.snp_start;
                    auto snp_chunk_end = chrom_snps.snps.begin() + wu.snp_end;
                    
                    // Get chromosome name for progress reporting
                    const char* chrom_name = header->target_name[tid];
                    long chunk_snp_count = wu.snp_end - wu.snp_start;
                    
                    int iterator_result = 0;
                    while ((iterator_result = sam_itr_next(bam_fp, iter, record)) >= 0){
                        // Apply the named production read-filter policy.
                        // V1 filters: unmapped, secondary, qcfail, dup
                        if (!read_passes_filter(record, default_production_read_filter())){
                            continue;
                        }
                        
                        // Count ALL reads that pass flag filter (before CB check) - matches V1
                        local_all_reads++;
                        
                        // Progress within large chunks (every 5M iterator records)
                        if (local_all_reads % 5000000 == 0){
                            fprintf(stderr, "\r  [%s:%d-%d] %ldM iterator records, %ld/%ld SNPs...          ",
                                chrom_name, wu.start_pos, wu.end_pos, 
                                local_all_reads / 1000000, local_snps, chunk_snp_count);
                        }
                        
                        int read_start = record->core.pos;
                        int read_end = bam_endpos(record);
                        
                        // Advance SNP iterator past SNPs before this read (within our chunk)
                        while (snp_iter != snp_chunk_end && snp_iter->pos < read_start){
                            ++snp_iter;
                            ++local_snps;
                        }
                        
                        // Bin-skip: if no SNP bin overlaps this read, skip CB extraction
                        {
                            auto bin_it = bin_indices.find(tid);
                            if (bin_it != bin_indices.end() &&
                                !bin_it->second.might_overlap(read_start, read_end)){
                                continue;
                            }
                        }
                        
                        // Extract cell barcode
                        uint8_t* cb_tag = bam_aux_get(record, "CB");
                        if (!cb_tag) continue;
                        
                        const char* cb_str = bam_aux2Z(cb_tag);
                        bc cb_bits;
                        str2bc(cb_str, cb_bits);
                        unsigned long bc_key = cb_bits.to_ulong();
                        
                        // Skip if not in whitelist
                        if (has_bc_list && valid_barcodes.find(bc_key) == valid_barcodes.end()){
                            continue;
                        }
                        
                        local_reads++;
                        
                        // Get mapping quality probability and scale to fixed-point
                        float prob_correct = 1.0f - powf(10.0f, -(float)record->core.qual / 10.0f);
                        int64_t prob_scaled = (int64_t)(prob_correct * FIXED_POINT_SCALE);
                        
                        // Process all SNPs overlapping this read (within our chunk)
                        for (auto snp_check = snp_iter; 
                             snp_check != snp_chunk_end && snp_check->pos < read_end; 
                             ++snp_check){
                            
                            char allele = get_base_at_pos(record, snp_check->pos);
                            if (allele == 'N' || allele == '-') continue;
                            
                            int64_t ref_add = 0, alt_add = 0;
                            
                            if (allele == snp_check->data.ref){
                                ref_add = prob_scaled;
                            }
                            else if (allele == snp_check->data.alt){
                                alt_add = prob_scaled;
                            }
                            else {
                            }
                            
                            if (ref_add > 0 || alt_add > 0){
                                if (accepted_site_weights != nullptr){
                                    local_site_weights[accepted_site_weight_key(tid, snp_check->pos)] +=
                                        ref_add + alt_add;
                                }
                                // Accumulate to per-thread storage (no locks needed)
                                auto it = local_counts.find(bc_key);
                                if (it == local_counts.end()){
                                    local_counts.emplace(bc_key, CellCounts(n_samples));
                                    it = local_counts.find(bc_key);
                                }
                                CellCounts& cc = it->second;
                                
                                // Precomputed targets: linear traversal, no branches
                                const auto& ttargets = snp_check->total_targets;
                                const auto& ptargets = snp_check->pair_targets;
                                
                                for (const auto& t : ttargets){
                                    cc.total_ref[t.total_idx] += ref_add;
                                    cc.total_alt[t.total_idx] += alt_add;
                                }
                                for (const auto& p : ptargets){
                                    cc.ref_counts[p.pair_idx] += ref_add;
                                    cc.alt_counts[p.pair_idx] += alt_add;
                                }

                                if (local_species_native != nullptr &&
                                    native_chrom_targets != nullptr){
                                    const size_t snp_index = (size_t)(
                                        snp_check - chrom_snps.snps.begin());
                                    if (snp_index < native_chrom_targets->site_offsets.size() &&
                                        native_chrom_targets->site_offsets[snp_index] != UINT64_MAX){
                                        auto native_it = local_species_native->find(bc_key);
                                        if (native_it == local_species_native->end()){
                                            local_species_native->emplace(
                                                bc_key, CellCounts(species_native_n_samples));
                                            native_it = local_species_native->find(bc_key);
                                        }
                                        accumulate_species_native_targets(
                                            native_it->second, *native_chrom_targets, snp_index,
                                            ref_add, alt_add);
                                    }
                                }

                                // --dump_pileup: record per-(cell,SNP) evidence for
                                // interindividual SNPs.  Summed within a thread; the
                                // producer sums any cross-thread duplicates.
                                if (dump_pileup && snp_check->panel_id == 0){
                                    int64_t pkey = ((int64_t)tid << 32) |
                                        (int64_t)(uint32_t)snp_check->pos;
                                    auto& slot = thread_pileup[thread_id][bc_key][pkey];
                                    slot.first += ref_add;
                                    slot.second += alt_add;
                                }
                            }
                        }
                    }
                    
                    if (iterator_result < -1){
                        operation_status.fail(format_worker_error("iterator read", thread_id, tid));
                    }

                    // Count remaining SNPs in this chunk
                    while (snp_iter != snp_chunk_end){
                        ++snp_iter;
                        ++local_snps;
                    }
                    
                    snps_processed += local_snps;
                    reads_processed += local_all_reads;
                    int done = ++units_done;
                    
                    if (done % 10 == 0 || done == (int)work_units.size()){
                        fprintf(stderr, "\rProgress: %d/%lu units, %ld/%ld SNPs, %ld iterator records          ",
                            done, work_units.size(), snps_processed.load(), total_snps,
                            reads_processed.load());
                    }
                    
                    hts_itr_destroy(iter);
        }

        if (record) bam_destroy1(record);
        if (idx) hts_idx_destroy(idx);
        if (header) bam_hdr_destroy(header);
        if (bam_fp) hts_close(bam_fp);
    }

    if (!operation_status.ok()){
        fprintf(stderr, "ERROR: parallel allele counting failed: %s\n",
            operation_status.message().c_str());
        return false;
    }
    
    // Merge per-thread counts (order doesn't matter with fixed-point int64 arithmetic)
    fprintf(stderr, "\nMerging per-thread counts (fixed-point)...\n");
    
    // Efficient merge: iterate each thread's map once and merge directly
    // This is O(total_entries) instead of O(cells * threads) lookups
    for (int t = 0; t < n_threads; t++){
        for (auto& kv : thread_counts[t]){
            unsigned long bc_key = kv.first;
            auto it = cell_counts.find(bc_key);
            if (it == cell_counts.end()){
                cell_counts.emplace(std::piecewise_construct,
                    std::forward_as_tuple(bc_key),
                    std::forward_as_tuple(n_samples));
                it = cell_counts.find(bc_key);
            }
            it->second.counts.merge(kv.second);
        }
        // Free memory incrementally
        thread_counts[t].clear();
        
        if ((t + 1) % 10 == 0 || t == n_threads - 1){
            fprintf(stderr, "  Merged thread %d/%d\n", t + 1, n_threads);
        }
    }
    
    if (collect_species_native){
        fprintf(stderr, "Merging per-thread native species counts...\n");
        for (int t = 0; t < n_threads; ++t){
            for (auto& kv : thread_species_native_counts[t]){
                const unsigned long bc_key = kv.first;
                auto it = species_native_counts->find(bc_key);
                if (it == species_native_counts->end()){
                    species_native_counts->emplace(std::piecewise_construct,
                        std::forward_as_tuple(bc_key),
                        std::forward_as_tuple(species_native_n_samples));
                    it = species_native_counts->find(bc_key);
                }
                it->second.counts.merge(kv.second);
            }
            thread_species_native_counts[t].clear();
        }
        thread_species_native_counts.clear();
        thread_species_native_counts.shrink_to_fit();
        fprintf(stderr, "Native species counts: %lu cells, %d species\n",
            species_native_counts->size(), species_native_n_samples);
    }

    if (accepted_site_weights != nullptr){
        accepted_site_weights->clear();
        for (int t = 0; t < n_threads; ++t){
            for (const auto& kv : thread_site_weights[t]){
                (*accepted_site_weights)[kv.first] += kv.second;
            }
            thread_site_weights[t].clear();
        }
        fprintf(stderr, "Accepted-site weight map: %lu observed sites\n",
            accepted_site_weights->size());
    }
    thread_site_weights.clear();
    thread_site_weights.shrink_to_fit();

    // Final cleanup
    thread_counts.clear();
    thread_counts.shrink_to_fit();

    // --dump_pileup: flush per-thread per-(cell,SNP) observations.  Per-thread
    // rows are pre-summed; the same (bc,tid,pos) may still appear across threads
    // and is summed downstream by the producer.  Columns: bc_hash tid pos ref alt
    // (ref/alt are prob-scaled float evidence, matching the .counts units).
    if (dump_pileup){
        string obs_path = pileup_prefix + ".pileup_obs.tsv.gz";
        gzFile pf = gzopen(obs_path.c_str(), "w");
        if (!pf){
            fprintf(stderr, "ERROR: could not open %s for writing\n", obs_path.c_str());
            return false;
        } else {
            long n_obs_written = 0;
            for (int t = 0; t < n_threads; t++){
                for (auto& cell : thread_pileup[t]){
                    unsigned long bc = cell.first;
                    for (auto& kv : cell.second){
                        int tid_o = (int)(kv.first >> 32);
                        int pos_o = (int)(kv.first & 0xFFFFFFFF);
                        gzprintf(pf, "%lu\t%d\t%d\t%f\t%f\n", bc, tid_o, pos_o,
                            (double)kv.second.first / FIXED_POINT_SCALE,
                            (double)kv.second.second / FIXED_POINT_SCALE);
                        n_obs_written++;
                    }
                }
                thread_pileup[t].clear();
            }
            if (gzclose(pf) != Z_OK){
                fprintf(stderr, "ERROR: failed while closing %s\n", obs_path.c_str());
                return false;
            }
            fprintf(stderr, "Wrote %ld pileup observations to %s\n", n_obs_written, obs_path.c_str());
        }
    }
    thread_pileup.clear();
    thread_pileup.shrink_to_fit();
    
    fprintf(stderr, "Completed: %d chromosomes (%lu work units), %ld SNPs, %ld iterator records, %lu cells\n",
        n_chroms, work_units.size(), snps_processed.load(), reads_processed.load(), 
        cell_counts.size());
    return true;
}

// ============================================================================
// SYNTHETIC SOURCE-PROVENANCE OBSERVATION SUMMARY
// ============================================================================

struct SourceObservationStats {
    uint64_t n_reads = 0;
    uint64_t n_observations = 0;
    double ref_weight = 0.0;
    double alt_weight = 0.0;

    void merge(const SourceObservationStats& other){
        n_reads += other.n_reads;
        n_observations += other.n_observations;
        ref_weight += other.ref_weight;
        alt_weight += other.alt_weight;
    }
};

struct SourceReceiverInfo {
    std::string identity;
    std::string name_a;
    std::string name_b;
    int idx_a = -1;
    int idx_b = -1;
};

using SourceReceiverMap = std::unordered_map<unsigned long, SourceReceiverInfo>;

static std::string trim_copy(const std::string& value){
    const size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

static std::string barcode_core_copy(const std::string& value){
    std::string out = trim_copy(value);
    const size_t dash = out.find('-');
    if (dash != std::string::npos) out.resize(dash);
    return out;
}

static SourceReceiverMap load_source_receiver_map(
    const std::string& path,
    const std::map<std::string, int>& sample_to_idx){
    SourceReceiverMap result;
    if (path.empty()) return result;
    std::ifstream in(path.c_str());
    if (!in){
        fprintf(stderr, "ERROR: could not open source receiver map: %s\n", path.c_str());
        exit(1);
    }
    std::string line;
    size_t line_no = 0;
    while (std::getline(in, line)){
        line_no++;
        line = trim_copy(line);
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string barcode, identity;
        if (!(ss >> barcode >> identity)){
            fprintf(stderr, "ERROR: %s:%lu: expected barcode and receiver_identity\n",
                path.c_str(), (unsigned long)line_no);
            exit(1);
        }
        if (barcode == "barcode" && (identity == "receiver_identity" || identity == "identity")){
            continue;
        }
        const size_t plus = identity.find('+');
        if (plus == std::string::npos || identity.find('+', plus + 1) != std::string::npos){
            fprintf(stderr, "ERROR: %s:%lu: receiver_identity must be exactly A+B, observed '%s'\n",
                path.c_str(), (unsigned long)line_no, identity.c_str());
            exit(1);
        }
        const std::string name_a = identity.substr(0, plus);
        const std::string name_b = identity.substr(plus + 1);
        auto ita = sample_to_idx.find(name_a);
        auto itb = sample_to_idx.find(name_b);
        if (ita == sample_to_idx.end() || itb == sample_to_idx.end()){
            fprintf(stderr, "ERROR: %s:%lu: receiver '%s' contains a name absent from the individual VCF sample set\n",
                path.c_str(), (unsigned long)line_no, identity.c_str());
            exit(1);
        }
        const std::string core = barcode_core_copy(barcode);
        if (core.empty()){
            fprintf(stderr, "ERROR: %s:%lu: empty barcode\n", path.c_str(), (unsigned long)line_no);
            exit(1);
        }
        bc bits;
        str2bc(core.c_str(), bits);
        const unsigned long key = bits.to_ulong();
        SourceReceiverInfo info;
        info.identity = identity;
        info.name_a = name_a;
        info.name_b = name_b;
        info.idx_a = ita->second;
        info.idx_b = itb->second;
        auto existing = result.find(key);
        if (existing != result.end() && existing->second.identity != identity){
            fprintf(stderr, "ERROR: %s:%lu: barcode '%s' maps to conflicting identities '%s' and '%s'\n",
                path.c_str(), (unsigned long)line_no, core.c_str(),
                existing->second.identity.c_str(), identity.c_str());
            exit(1);
        }
        result[key] = info;
    }
    if (result.empty()){
        fprintf(stderr, "ERROR: source receiver map is empty after parsing: %s\n", path.c_str());
        exit(1);
    }
    fprintf(stderr, "Loaded %lu barcode-to-receiver mappings from %s\n",
        (unsigned long)result.size(), path.c_str());
    return result;
}

// barcode, receiver identity/A/B/category, explicit origin bucket,
// resolved source, raw YI source, typing status, panel. The receiver category
// is populated for the individual panel when --source_receiver_map is supplied;
// species rows use -1/-1.
using SourceObservationKey = std::tuple<
    unsigned long,
    std::string, std::string, std::string, int, int,
    std::string, std::string, std::string, std::string, int>;
using SourceObservationMap = std::map<SourceObservationKey, SourceObservationStats>;

using SourceCategoryKey = std::tuple<
    unsigned long,
    std::string, std::string, std::string, int, int, int>;
using SourceCategoryMap = std::map<SourceCategoryKey, SourceObservationStats>;

// Complete audit aggregate.  Every accepted injected individual-panel
// observation is accounted, but rows are compactly aggregated by receiver,
// source label, exact source genotypes, and typing state rather than emitted
// read-by-read.
struct DonorAuditStats {
    uint64_t n_observations = 0;
    double ref_weight = 0.0;
    double alt_weight = 0.0;
    double direct_supported_weight = 0.0;
    double direct_unsupported_weight = 0.0;
    double flipped_supported_weight = 0.0;
    double unique_component_a_weight = 0.0;
    double unique_component_b_weight = 0.0;
    double ambiguous_component_weight = 0.0;
    double no_component_weight = 0.0;
    double raw_equalmix_expected_alt_num = 0.0;
    double raw_equalmix_expected_weight = 0.0;
    double resolved_expected_alt_num = 0.0;
    double resolved_expected_weight = 0.0;

    void merge(const DonorAuditStats& other){
        n_observations += other.n_observations;
        ref_weight += other.ref_weight;
        alt_weight += other.alt_weight;
        direct_supported_weight += other.direct_supported_weight;
        direct_unsupported_weight += other.direct_unsupported_weight;
        flipped_supported_weight += other.flipped_supported_weight;
        unique_component_a_weight += other.unique_component_a_weight;
        unique_component_b_weight += other.unique_component_b_weight;
        ambiguous_component_weight += other.ambiguous_component_weight;
        no_component_weight += other.no_component_weight;
        raw_equalmix_expected_alt_num += other.raw_equalmix_expected_alt_num;
        raw_equalmix_expected_weight += other.raw_equalmix_expected_weight;
        resolved_expected_alt_num += other.resolved_expected_alt_num;
        resolved_expected_weight += other.resolved_expected_weight;
    }
};

// receiver identity/A/B/category, raw YI, resolved source/typing state,
// raw-source component names and exact site genotypes, resolved genotype.
using DonorAuditKey = std::tuple<
    std::string, std::string, std::string, int, int,
    std::string, std::string, std::string,
    std::string, std::string, int, int, int>;
using DonorAuditMap = std::map<DonorAuditKey, DonorAuditStats>;

// Deterministic exact-site evidence sample.  The complete aggregate above uses
// every accepted observation; this sampled table provides inspectable loci
// without generating a multi-gigabyte read/site dump.
using DonorSiteAuditKey = std::tuple<
    int, int, char, char,
    std::string, std::string, std::string, int, int,
    std::string, std::string, std::string,
    std::string, std::string, int, int, int>;
using DonorSiteAuditMap = std::map<DonorSiteAuditKey, DonorAuditStats>;

static bool genotype_supports_observed_allele(int genotype, bool is_ref){
    if (genotype < 0 || genotype > 2) return false;
    if (genotype == 1) return true;
    return is_ref ? genotype == 0 : genotype == 2;
}

static bool donor_site_is_sampled(int tid, int pos, int modulus){
    if (modulus <= 1) return true;
    uint64_t x = ((uint64_t)(uint32_t)tid << 32) ^ (uint64_t)(uint32_t)pos;
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return (x % (uint64_t)modulus) == 0;
}

static void describe_source_genotypes(
    const std::string& raw_source,
    const std::string& resolved_source,
    const SNPData& snp,
    const std::map<std::string, int>& sample_to_idx,
    std::string& source_a,
    std::string& source_b,
    int& source_nalt_a,
    int& source_nalt_b,
    int& resolved_nalt){
    source_a.clear();
    source_b.clear();
    source_nalt_a = -1;
    source_nalt_b = -1;
    resolved_nalt = -1;

    auto genotype_for = [&](const std::string& name) -> int {
        auto hit = sample_to_idx.find(name);
        if (hit == sample_to_idx.end()) return -1;
        const int idx = hit->second;
        if (idx < 0 || idx >= (int)snp.geno.size()) return -1;
        return (int)snp.geno[idx];
    };

    const size_t plus = raw_source.find('+');
    if (plus == std::string::npos){
        source_a = raw_source;
        source_nalt_a = genotype_for(source_a);
    }
    else if (raw_source.find('+', plus + 1) == std::string::npos){
        source_a = raw_source.substr(0, plus);
        source_b = raw_source.substr(plus + 1);
        source_nalt_a = genotype_for(source_a);
        source_nalt_b = genotype_for(source_b);
        if (source_a == source_b){
            source_b.clear();
            source_nalt_b = -1;
        }
    }
    resolved_nalt = genotype_for(resolved_source);
}

static void update_donor_audit_stats(
    DonorAuditStats& stats,
    bool is_ref,
    double weight,
    int source_nalt_a,
    int source_nalt_b,
    int resolved_nalt){
    stats.n_observations += 1;
    if (is_ref) stats.ref_weight += weight;
    else stats.alt_weight += weight;

    const bool a_supports = genotype_supports_observed_allele(source_nalt_a, is_ref);
    const bool b_present = source_nalt_b >= 0 && source_nalt_b <= 2;
    const bool b_supports = b_present && genotype_supports_observed_allele(source_nalt_b, is_ref);
    const bool any_supports = a_supports || b_supports;
    if (any_supports) stats.direct_supported_weight += weight;
    else stats.direct_unsupported_weight += weight;

    const bool flipped_a = genotype_supports_observed_allele(source_nalt_a, !is_ref);
    const bool flipped_b = b_present && genotype_supports_observed_allele(source_nalt_b, !is_ref);
    if (flipped_a || flipped_b) stats.flipped_supported_weight += weight;

    if (a_supports && !b_supports) stats.unique_component_a_weight += weight;
    else if (b_supports && !a_supports) stats.unique_component_b_weight += weight;
    else if (a_supports && b_supports) stats.ambiguous_component_weight += weight;
    else stats.no_component_weight += weight;

    if (source_nalt_a >= 0 && source_nalt_a <= 2){
        double p = (double)source_nalt_a / 2.0;
        if (b_present) p = ((double)source_nalt_a + (double)source_nalt_b) / 4.0;
        stats.raw_equalmix_expected_alt_num += weight * p;
        stats.raw_equalmix_expected_weight += weight;
    }
    if (resolved_nalt >= 0 && resolved_nalt <= 2){
        stats.resolved_expected_alt_num += weight * ((double)resolved_nalt / 2.0);
        stats.resolved_expected_weight += weight;
    }
}

static bool split_pair_source_label(
    const std::string& raw,
    const std::map<std::string, int>& sample_to_idx,
    int& idx_a,
    int& idx_b,
    std::string& name_a,
    std::string& name_b){
    const size_t plus = raw.find('+');
    if (plus == std::string::npos || raw.find('+', plus + 1) != std::string::npos) return false;
    name_a = raw.substr(0, plus);
    name_b = raw.substr(plus + 1);
    auto a = sample_to_idx.find(name_a);
    auto b = sample_to_idx.find(name_b);
    if (name_a.empty() || name_b.empty() || a == sample_to_idx.end() || b == sample_to_idx.end()) return false;
    idx_a = a->second;
    idx_b = b->second;
    return true;
}

static bool source_label_is_known(
    const std::string& raw,
    const std::map<std::string, int>& sample_to_idx){
    if (raw.empty()) return false;
    if (raw.find('+') == std::string::npos){
        return sample_to_idx.find(raw) != sample_to_idx.end();
    }
    int idx_a = -1, idx_b = -1;
    std::string name_a, name_b;
    return split_pair_source_label(raw, sample_to_idx, idx_a, idx_b, name_a, name_b);
}

static void resolve_source_observation(
    const std::string& raw_source,
    const SNPData& snp,
    bool is_ref,
    const std::map<std::string, int>& sample_to_idx,
    std::string& resolved_source,
    std::string& typing_status){
    int idx_a = -1, idx_b = -1;
    std::string name_a, name_b;
    if (raw_source.find('+') == std::string::npos){
        resolved_source = raw_source;
        typing_status = "singleton";
        return;
    }
    if (!split_pair_source_label(raw_source, sample_to_idx, idx_a, idx_b, name_a, name_b)){
        resolved_source = raw_source;
        typing_status = "composite_unmapped";
        return;
    }
    if (idx_a < 0 || idx_b < 0 || idx_a >= (int)snp.geno.size() || idx_b >= (int)snp.geno.size()){
        resolved_source = raw_source;
        typing_status = "composite_unmapped";
        return;
    }
    // A+A is a homotypic cell label, not an ambiguous source mixture: every
    // emitted read originates from the same biological individual A.
    if (name_a == name_b){
        resolved_source = name_a;
        typing_status = "homotypic_composite";
        return;
    }
    const int ga = (int)snp.geno[idx_a];
    const int gb = (int)snp.geno[idx_b];
    // B35 exact primitive: only homozygous-opposite constituent sites type the
    // source observation. Heterozygous/agreement sites remain explicitly
    // ambiguous and are never forced to a 50:50 split.
    if (ga == 0 && gb == 2){
        resolved_source = is_ref ? name_a : name_b;
        typing_status = "typed_composite";
    }
    else if (ga == 2 && gb == 0){
        resolved_source = is_ref ? name_b : name_a;
        typing_status = "typed_composite";
    }
    else{
        resolved_source = raw_source;
        typing_status = "composite_ambiguous";
    }
}

static void write_source_observation_summary(
    const std::string& prefix,
    const SourceObservationMap& observations){
    if (prefix.empty()) return;
    const std::string out_path = prefix + ".source_observations.tsv.gz";
    gzFile out = gzopen(out_path.c_str(), "wb");
    if (!out){
        fprintf(stderr, "ERROR: could not open %s for source-observation output\n", out_path.c_str());
        exit(1);
    }
    gzprintf(out, "barcode\treceiver_identity\treceiver_A\treceiver_B\tnalt_A\tnalt_B\torigin\tsource\traw_source\ttyped_status\tpanel\tn_reads\tn_observations\tref_weight\talt_weight\tweighted_observations\n");
    for (const auto& kv : observations){
        const unsigned long bc_key = std::get<0>(kv.first);
        const std::string& receiver_identity = std::get<1>(kv.first);
        const std::string& receiver_a = std::get<2>(kv.first);
        const std::string& receiver_b = std::get<3>(kv.first);
        const int nalt_a = std::get<4>(kv.first);
        const int nalt_b = std::get<5>(kv.first);
        const std::string& origin = std::get<6>(kv.first);
        const std::string& source = std::get<7>(kv.first);
        const std::string& raw_source = std::get<8>(kv.first);
        const std::string& typed_status = std::get<9>(kv.first);
        const int panel_id = std::get<10>(kv.first);
        const auto& st = kv.second;
        const std::string barcode = bc2str(bc_key);
        const char* panel = panel_id == 0 ? "individual" : "species";
        gzprintf(out, "%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%llu\t%llu\t%.10g\t%.10g\t%.10g\n",
            barcode.c_str(), receiver_identity.c_str(), receiver_a.c_str(), receiver_b.c_str(),
            nalt_a, nalt_b, origin.c_str(), source.c_str(), raw_source.c_str(), typed_status.c_str(), panel,
            (unsigned long long)st.n_reads,
            (unsigned long long)st.n_observations,
            st.ref_weight, st.alt_weight, st.ref_weight + st.alt_weight);
    }
    gzclose(out);
    fprintf(stderr, "Wrote %lu source-provenance rows to %s\n", observations.size(), out_path.c_str());
}

static void write_donor_genotype_audit(
    const std::string& prefix,
    const DonorAuditMap& audit){
    if (prefix.empty()) return;
    const std::string out_path = prefix + ".donor_genotype_audit.tsv.gz";
    gzFile out = gzopen(out_path.c_str(), "wb");
    if (!out){
        fprintf(stderr, "ERROR: could not open %s for donor-genotype audit output\n", out_path.c_str());
        exit(1);
    }
    gzprintf(out,
        "receiver_identity\treceiver_A\treceiver_B\tnalt_A\tnalt_B"
        "\traw_source\tresolved_source\ttyped_status\tsource_A\tsource_B"
        "\tsource_nalt_A\tsource_nalt_B\tresolved_nalt\tn_observations"
        "\tref_weight\talt_weight\tweighted_observations\tempirical_alt_fraction"
        "\tdirect_supported_weight\tdirect_unsupported_weight\tdirect_unsupported_fraction"
        "\tflipped_supported_weight\tunique_component_A_weight\tunique_component_B_weight"
        "\tambiguous_component_weight\tno_component_weight"
        "\traw_equalmix_expected_alt_fraction\tresolved_expected_alt_fraction\n");
    for (const auto& kv : audit){
        const auto& key = kv.first;
        const auto& st = kv.second;
        const double total = st.ref_weight + st.alt_weight;
        const double empirical = total > 0 ? st.alt_weight / total : NAN;
        const double unsupported = total > 0 ? st.direct_unsupported_weight / total : NAN;
        const double raw_expected = st.raw_equalmix_expected_weight > 0
            ? st.raw_equalmix_expected_alt_num / st.raw_equalmix_expected_weight : NAN;
        const double resolved_expected = st.resolved_expected_weight > 0
            ? st.resolved_expected_alt_num / st.resolved_expected_weight : NAN;
        gzprintf(out,
            "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d"
            "\t%llu\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g"
            "\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
            std::get<0>(key).c_str(), std::get<1>(key).c_str(), std::get<2>(key).c_str(),
            std::get<3>(key), std::get<4>(key), std::get<5>(key).c_str(),
            std::get<6>(key).c_str(), std::get<7>(key).c_str(), std::get<8>(key).c_str(),
            std::get<9>(key).c_str(), std::get<10>(key), std::get<11>(key), std::get<12>(key),
            (unsigned long long)st.n_observations,
            st.ref_weight, st.alt_weight, total, empirical,
            st.direct_supported_weight, st.direct_unsupported_weight, unsupported,
            st.flipped_supported_weight, st.unique_component_a_weight,
            st.unique_component_b_weight, st.ambiguous_component_weight,
            st.no_component_weight, raw_expected, resolved_expected);
    }
    gzclose(out);
    fprintf(stderr, "Wrote %lu donor-genotype audit rows to %s\n",
        (unsigned long)audit.size(), out_path.c_str());
}

static void write_donor_site_sample(
    const std::string& prefix,
    const DonorSiteAuditMap& audit,
    const std::vector<std::string>& chrom_names,
    int sample_mod){
    if (prefix.empty()) return;
    const std::string out_path = prefix + ".donor_site_sample.tsv.gz";
    gzFile out = gzopen(out_path.c_str(), "wb");
    if (!out){
        fprintf(stderr, "ERROR: could not open %s for donor-site sample output\n", out_path.c_str());
        exit(1);
    }
    gzprintf(out,
        "chrom\ttid\tpos0\tpos1\tref\talt\treceiver_identity\treceiver_A\treceiver_B"
        "\tnalt_A\tnalt_B\traw_source\tresolved_source\ttyped_status\tsource_A\tsource_B"
        "\tsource_nalt_A\tsource_nalt_B\tresolved_nalt\tn_observations\tref_weight\talt_weight"
        "\tweighted_observations\tempirical_alt_fraction\tdirect_unsupported_fraction"
        "\traw_equalmix_expected_alt_fraction\tresolved_expected_alt_fraction\tsample_mod\n");
    for (const auto& kv : audit){
        const auto& key = kv.first;
        const auto& st = kv.second;
        const int tid = std::get<0>(key);
        const int pos = std::get<1>(key);
        const std::string chrom = tid >= 0 && tid < (int)chrom_names.size()
            ? chrom_names[tid] : std::to_string(tid);
        const double total = st.ref_weight + st.alt_weight;
        const double empirical = total > 0 ? st.alt_weight / total : NAN;
        const double unsupported = total > 0 ? st.direct_unsupported_weight / total : NAN;
        const double raw_expected = st.raw_equalmix_expected_weight > 0
            ? st.raw_equalmix_expected_alt_num / st.raw_equalmix_expected_weight : NAN;
        const double resolved_expected = st.resolved_expected_weight > 0
            ? st.resolved_expected_alt_num / st.resolved_expected_weight : NAN;
        gzprintf(out,
            "%s\t%d\t%d\t%d\t%c\t%c\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s"
            "\t%d\t%d\t%d\t%llu\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%d\n",
            chrom.c_str(), tid, pos, pos + 1, std::get<2>(key), std::get<3>(key),
            std::get<4>(key).c_str(), std::get<5>(key).c_str(), std::get<6>(key).c_str(),
            std::get<7>(key), std::get<8>(key), std::get<9>(key).c_str(),
            std::get<10>(key).c_str(), std::get<11>(key).c_str(), std::get<12>(key).c_str(),
            std::get<13>(key).c_str(), std::get<14>(key), std::get<15>(key), std::get<16>(key),
            (unsigned long long)st.n_observations, st.ref_weight, st.alt_weight, total,
            empirical, unsupported, raw_expected, resolved_expected, sample_mod);
    }
    gzclose(out);
    fprintf(stderr, "Wrote %lu deterministic donor-site sample rows (mod=%d) to %s\n",
        (unsigned long)audit.size(), sample_mod, out_path.c_str());
}

static SourceCategoryKey observation_category_key(const SourceObservationKey& key){
    return SourceCategoryKey(
        std::get<0>(key), std::get<1>(key), std::get<2>(key), std::get<3>(key),
        std::get<4>(key), std::get<5>(key), std::get<10>(key));
}

static void write_source_reconciliation_summary(
    const std::string& prefix,
    const SourceObservationMap& observations,
    const SourceCategoryMap& raw_categories){
    if (prefix.empty()) return;
    const std::string out_path = prefix + ".source_reconciliation.tsv.gz";
    gzFile out = gzopen(out_path.c_str(), "wb");
    if (!out){
        fprintf(stderr, "ERROR: could not open %s for source-reconciliation output\n", out_path.c_str());
        exit(1);
    }

    using OriginMap = std::map<std::string, SourceObservationStats>;
    std::map<SourceCategoryKey, OriginMap> by_origin;
    for (const auto& kv : observations){
        if (std::get<10>(kv.first) != 0) continue;
        by_origin[observation_category_key(kv.first)][std::get<6>(kv.first)].merge(kv.second);
    }

    std::set<SourceCategoryKey> keys;
    for (const auto& kv : raw_categories){
        if (std::get<6>(kv.first) == 0) keys.insert(kv.first);
    }
    for (const auto& kv : by_origin) keys.insert(kv.first);

    const std::vector<std::string> origins = {
        "native", "native_untagged", "injected_expected_ys",
        "injected_missing_yi", "injected_missing_ys", "injected_wrong_ys",
        "missing_yi", "invalid_yi", "unknown_yi"
    };
    gzprintf(out,
        "barcode\treceiver_identity\treceiver_A\treceiver_B\tnalt_A\tnalt_B\tpanel"
        "\traw_n_observations\traw_ref\traw_alt");
    for (const auto& origin : origins){
        gzprintf(out, "\t%s_n_observations\t%s_ref\t%s_alt",
            origin.c_str(), origin.c_str(), origin.c_str());
    }
    gzprintf(out,
        "\taccounted_ref\taccounted_alt\tref_difference_raw_minus_accounted"
        "\talt_difference_raw_minus_accounted\treconciliation_pass\n");

    uint64_t failed_rows = 0;
    double abs_ref_error = 0.0;
    double abs_alt_error = 0.0;
    for (const auto& key : keys){
        const auto raw_it = raw_categories.find(key);
        const SourceObservationStats empty;
        const SourceObservationStats& raw = raw_it == raw_categories.end() ? empty : raw_it->second;
        const auto origin_it = by_origin.find(key);
        double accounted_ref = 0.0;
        double accounted_alt = 0.0;
        const std::string barcode = bc2str(std::get<0>(key));
        gzprintf(out, "%s\t%s\t%s\t%s\t%d\t%d\tindividual\t%llu\t%.17g\t%.17g",
            barcode.c_str(), std::get<1>(key).c_str(), std::get<2>(key).c_str(),
            std::get<3>(key).c_str(), std::get<4>(key), std::get<5>(key),
            (unsigned long long)raw.n_observations, raw.ref_weight, raw.alt_weight);
        for (const auto& origin : origins){
            SourceObservationStats st;
            if (origin_it != by_origin.end()){
                auto hit = origin_it->second.find(origin);
                if (hit != origin_it->second.end()) st = hit->second;
            }
            accounted_ref += st.ref_weight;
            accounted_alt += st.alt_weight;
            gzprintf(out, "\t%llu\t%.17g\t%.17g",
                (unsigned long long)st.n_observations, st.ref_weight, st.alt_weight);
        }
        const double ref_diff = raw.ref_weight - accounted_ref;
        const double alt_diff = raw.alt_weight - accounted_alt;
        const double ref_tol = 1e-7 + 1e-12 * std::max(1.0, std::fabs(raw.ref_weight));
        const double alt_tol = 1e-7 + 1e-12 * std::max(1.0, std::fabs(raw.alt_weight));
        const bool pass = std::fabs(ref_diff) <= ref_tol && std::fabs(alt_diff) <= alt_tol;
        if (!pass) failed_rows++;
        abs_ref_error += std::fabs(ref_diff);
        abs_alt_error += std::fabs(alt_diff);
        gzprintf(out, "\t%.17g\t%.17g\t%.17g\t%.17g\t%s\n",
            accounted_ref, accounted_alt, ref_diff, alt_diff, pass ? "PASS" : "FAIL");
    }
    gzclose(out);
    fprintf(stderr,
        "Source reconciliation: rows=%lu failed_rows=%llu absolute_ref_error=%.10g absolute_alt_error=%.10g -> %s\n",
        (unsigned long)keys.size(), (unsigned long long)failed_rows,
        abs_ref_error, abs_alt_error, out_path.c_str());
}

// ============================================================================
// DUAL-OUTPUT PARALLEL COUNTING (WP3: single BAM pass for two panels)
// ============================================================================

bool count_alleles_parallel_dual(
    const string& bamfile,
    robin_hood::unordered_map<int, ChromSNPs>& combined_snpdat,
    robin_hood::unordered_map<unsigned long, AlignedCellCounts>& counts_panel0,
    robin_hood::unordered_map<unsigned long, AlignedCellCounts>& counts_panel1,
    const set<unsigned long>& valid_barcodes,
    int n_samples,
    int n_threads,
    int htslib_threads){
    const vector<string> no_sample_names;
    return count_alleles_parallel_dual(
        bamfile, combined_snpdat, counts_panel0, counts_panel1, valid_barcodes,
        n_samples, no_sample_names, n_threads, htslib_threads,
        false, "", "YI", "YS", "", "", false, false, 256,
        nullptr, nullptr, nullptr, nullptr, 0,
        false, "");
}

bool count_alleles_parallel_dual(
    const string& bamfile,
    robin_hood::unordered_map<int, ChromSNPs>& combined_snpdat,
    robin_hood::unordered_map<unsigned long, AlignedCellCounts>& counts_panel0,
    robin_hood::unordered_map<unsigned long, AlignedCellCounts>& counts_panel1,
    const set<unsigned long>& valid_barcodes,
    int n_samples,
    const vector<string>& sample_names,
    int n_threads,
    int htslib_threads,
    bool dump_source_observations,
    const string& source_observations_prefix,
    const string& source_provenance_tag,
    const string& synthetic_id_tag,
    const string& expected_synthetic_id,
    const string& source_receiver_map_path,
    bool source_reconciliation_mode,
    bool source_donor_site_audit,
    int source_donor_site_sample_mod,
    AcceptedSiteWeightMap* accepted_site_weights_panel0,
    AcceptedSiteWeightMap* accepted_site_weights_panel1,
    const NativeSpeciesTargetTable* species_native_targets,
    robin_hood::unordered_map<unsigned long, AlignedCellCounts>* species_native_counts,
    int species_native_n_samples,
    bool dump_pileup,
    const string& pileup_prefix){

    size_t bytes_per_cell = 0;
    string request_error;
    if (!validate_identity_and_allocation_request(
            n_samples, nullptr, &bytes_per_cell, &request_error)){
        fprintf(stderr, "ERROR: invalid dual-panel identity universe: %s\n",
            request_error.c_str());
        return false;
    }
    if (species_native_n_samples > 0 && !validate_identity_and_allocation_request(
            species_native_n_samples, nullptr, nullptr, &request_error)){
        fprintf(stderr, "ERROR: invalid native-species identity universe: %s\n",
            request_error.c_str());
        return false;
    }
    if (n_threads < 1 || htslib_threads < 1){
        fprintf(stderr, "ERROR: thread counts must be positive\n");
        return false;
    }

    const bool has_bc_list = !valid_barcodes.empty();
    const bool collect_species_native =
        species_native_targets != nullptr && species_native_counts != nullptr &&
        species_native_n_samples > 0;

    if (has_bc_list){
        if (bytes_per_cell > 0 && valid_barcodes.size() >
            std::numeric_limits<size_t>::max() / (2 * bytes_per_cell)){
            fprintf(stderr, "ERROR: projected dual-panel CellCounts allocation overflows size_t\n");
            return false;
        }
        fprintf(stderr, "Pre-allocating dual counts for %lu cells (%lu bytes/cell/panel)...\n",
            valid_barcodes.size(), (unsigned long)bytes_per_cell);
        try {
            for (unsigned long bc : valid_barcodes){
                counts_panel0.emplace(std::piecewise_construct,
                    std::forward_as_tuple(bc),
                    std::forward_as_tuple(n_samples));
                counts_panel1.emplace(std::piecewise_construct,
                    std::forward_as_tuple(bc),
                    std::forward_as_tuple(n_samples));
                if (collect_species_native){
                    species_native_counts->emplace(std::piecewise_construct,
                        std::forward_as_tuple(bc),
                        std::forward_as_tuple(species_native_n_samples));
                }
            }
        }
        catch (const std::exception& e){
            fprintf(stderr, "ERROR: dual-panel pre-allocation failed: %s\n", e.what());
            return false;
        }
    }

    htsFile* bam_tmp = hts_open(bamfile.c_str(), "r");
    if (!bam_tmp){
        fprintf(stderr, "ERROR: Could not open BAM file to get header: %s\n", bamfile.c_str());
        return false;
    }
    bam_hdr_t* hdr_tmp = sam_hdr_read(bam_tmp);
    if (!hdr_tmp){
        fprintf(stderr, "ERROR: Could not read BAM header: %s\n", bamfile.c_str());
        hts_close(bam_tmp);
        return false;
    }
    hts_idx_t* idx_tmp = sam_index_load(bam_tmp, bamfile.c_str());
    if (!idx_tmp){
        fprintf(stderr, "ERROR: Could not load required BAM index: %s\n", bamfile.c_str());
        bam_hdr_destroy(hdr_tmp);
        hts_close(bam_tmp);
        return false;
    }
    const int n_chroms = hdr_tmp->n_targets;
    for (const auto& kv : combined_snpdat){
        if (kv.first < 0 || kv.first >= n_chroms){
            fprintf(stderr, "ERROR: combined SNP panel references invalid BAM target id %d\n", kv.first);
            hts_idx_destroy(idx_tmp);
            bam_hdr_destroy(hdr_tmp);
            hts_close(bam_tmp);
            return false;
        }
    }

    // --dump_pileup: emit the per-SNP genotype sidecar for the variant-consistency
    // metric (interindividual panel only; panel_id != 0 sites are skipped so the
    // dual-panel sidecar is byte-compatible with the single-panel sidecar).
    // geno[] is populated by precompute_all_genotypes() before this call, and
    // hdr_tmp is still valid here (destroyed below).
    // Columns: tid  chrom  pos  ref  alt  geno_0 .. geno_{n_samples-1}  (0/1/2/-1).
    // The (tid,pos) pair is the producer's SNP join key; ref/alt are the allele
    // bases (informational; the metric is allele-orientation based).
    if (dump_pileup){
        string sites_path = pileup_prefix + ".pileup_sites.tsv.gz";
        gzFile sf = gzopen(sites_path.c_str(), "w");
        if (!sf){
            fprintf(stderr, "ERROR: could not open %s for writing\n", sites_path.c_str());
            hts_idx_destroy(idx_tmp);
            bam_hdr_destroy(hdr_tmp);
            hts_close(bam_tmp);
            return false;
        } else {
            long n_sites_written = 0;
            for (auto& kv : combined_snpdat){
                int tid_s = kv.first;
                const char* cname = (tid_s >= 0 && tid_s < n_chroms) ?
                    hdr_tmp->target_name[tid_s] : ".";
                for (auto& snp : kv.second.snps){
                    if (snp.panel_id != 0) continue;
                    gzprintf(sf, "%d\t%s\t%d\t%c\t%c", tid_s, cname, snp.pos,
                        snp.data.ref, snp.data.alt);
                    for (int s = 0; s < n_samples; s++){
                        gzprintf(sf, "\t%d", (int)snp.geno[s]);
                    }
                    gzprintf(sf, "\n");
                    n_sites_written++;
                }
            }
            if (gzclose(sf) != Z_OK){
                fprintf(stderr, "ERROR: failed while closing %s\n", sites_path.c_str());
                hts_idx_destroy(idx_tmp);
                bam_hdr_destroy(hdr_tmp);
                hts_close(bam_tmp);
                return false;
            }
            fprintf(stderr, "Wrote %ld pileup sites to %s\n", n_sites_written, sites_path.c_str());
        }
    }

    vector<uint64_t> chrom_read_counts(n_chroms, 0);
    vector<int64_t> chrom_lengths(n_chroms);
    vector<string> chrom_names(n_chroms);
    const int n_index_targets = hts_idx_nseq(idx_tmp);
    if (n_index_targets < 0){
        fprintf(stderr, "ERROR: Could not determine BAM-index target count\n");
        hts_idx_destroy(idx_tmp);
        bam_hdr_destroy(hdr_tmp);
        hts_close(bam_tmp);
        return false;
    }
    int n_missing_index_stats = 0;
    int first_missing_index_stat = -1;
    for (int i = 0; i < n_chroms; ++i){
        chrom_names[i] = hdr_tmp->target_name[i]
            ? hdr_tmp->target_name[i] : std::to_string(i);
        uint64_t mapped = 0, unmapped = 0;
        if (i >= n_index_targets ||
            hts_idx_get_stat(idx_tmp, i, &mapped, &unmapped) < 0){
            // A valid BAI/CSI may omit the metadata bin for a target with no
            // alignment records. These counts are used only to estimate work
            // unit size and ordering; the iterator still queries the complete
            // target below. Use an unchunked zero-read scheduling estimate.
            mapped = 0;
            unmapped = 0;
            if (first_missing_index_stat < 0) first_missing_index_stat = i;
            n_missing_index_stats++;
        }
        chrom_read_counts[i] = mapped;
        chrom_lengths[i] = hdr_tmp->target_len[i];
    }
    if (n_missing_index_stats > 0){
        const char* first_name =
            (first_missing_index_stat >= 0 && first_missing_index_stat < n_chroms &&
             hdr_tmp->target_name[first_missing_index_stat])
                ? hdr_tmp->target_name[first_missing_index_stat] : ".";
        fprintf(stderr,
            "WARNING: BAM-index mapped/unmapped statistics unavailable for %d of %d header targets "
            "(first target %d: %s); using zero only for work-unit scheduling estimates\n",
            n_missing_index_stats, n_chroms, first_missing_index_stat, first_name);
    }
    hts_idx_destroy(idx_tmp);
    bam_hdr_destroy(hdr_tmp);
    hts_close(bam_tmp);

    // Build bin index for read skipping (Change 3)
    robin_hood::unordered_map<int, ChromBinIndex> bin_indices;
    build_bin_indices(combined_snpdat, chrom_lengths, bin_indices);
    
    // Work unit structure (same as count_alleles_parallel)
    struct WorkUnit {
        int tid;
        int start_pos;
        int end_pos;
        size_t snp_start;
        size_t snp_end;
        bool has_snps;
        uint64_t est_reads;
    };
    
    // Build work units (identical logic to count_alleles_parallel)
    vector<WorkUnit> work_units;
    long total_snps = 0;
    int chroms_with_snps = 0;
    int chroms_chunked_by_snp = 0;
    int chroms_chunked_by_reads = 0;
    
    const size_t CHUNK_SNP_THRESHOLD = 100000;
    const uint64_t CHUNK_READ_THRESHOLD = 10000000;
    
    for (int tid = 0; tid < n_chroms; tid++){
        auto it = combined_snpdat.find(tid);
        int64_t chrom_len = chrom_lengths[tid];
        uint64_t chrom_reads = chrom_read_counts[tid];
        
        if (it == combined_snpdat.end() || it->second.empty()){
            // Neither active panel has a SNP on this target. Skip it entirely;
            // record-count-only iterators do not affect either output panel.
            continue;
        }
        else{
            ChromSNPs& chrom_snps = it->second;
            size_t n_snps = chrom_snps.snps.size();
            total_snps += n_snps;
            chroms_with_snps++;
            
            size_t chunks_by_snp = (n_snps > CHUNK_SNP_THRESHOLD) ? 
                (n_snps + CHUNK_SNP_THRESHOLD - 1) / CHUNK_SNP_THRESHOLD : 1;
            size_t chunks_by_reads = (chrom_reads > CHUNK_READ_THRESHOLD) ?
                (chrom_reads + CHUNK_READ_THRESHOLD - 1) / CHUNK_READ_THRESHOLD : 1;
            
            chunks_by_reads = std::min(chunks_by_reads, (size_t)20);
            
            if (chunks_by_snp >= chunks_by_reads && chunks_by_snp > 1){
                size_t snps_per_chunk = (n_snps + chunks_by_snp - 1) / chunks_by_snp;
                uint64_t reads_per_chunk = chrom_reads / chunks_by_snp;
                
                for (size_t c = 0; c < chunks_by_snp; c++){
                    size_t snp_start = c * snps_per_chunk;
                    size_t snp_end = std::min(snp_start + snps_per_chunk, n_snps);
                    
                    int start_pos = (snp_start == 0) ? 0 : chrom_snps.snps[snp_start].pos;
                    int end_pos = (snp_end >= n_snps) ? INT_MAX : chrom_snps.snps[snp_end - 1].pos + 1000;
                    
                    work_units.push_back({tid, start_pos, end_pos, snp_start, snp_end, true, reads_per_chunk});
                }
                chroms_chunked_by_snp++;
            }
            else if (chunks_by_reads > 1){
                int64_t chunk_size = (chrom_len + chunks_by_reads - 1) / chunks_by_reads;
                uint64_t reads_per_chunk = chrom_reads / chunks_by_reads;
                
                size_t snp_idx = 0;
                for (size_t c = 0; c < chunks_by_reads; c++){
                    int start_pos = c * chunk_size;
                    int end_pos = (c == chunks_by_reads - 1) ? INT_MAX : (int)((c + 1) * chunk_size);
                    
                    size_t snp_start = snp_idx;
                    while (snp_idx < n_snps && chrom_snps.snps[snp_idx].pos < (c + 1) * chunk_size){
                        snp_idx++;
                    }
                    size_t snp_end = snp_idx;
                    if (snp_start == snp_end) continue;
                    
                    work_units.push_back({tid, start_pos, end_pos, snp_start, snp_end, true, reads_per_chunk});
                }
                chroms_chunked_by_reads++;
            }
            else{
                work_units.push_back({tid, 0, INT_MAX, 0, n_snps, true, chrom_reads});
            }
        }
    }
    
    // Sort work units by estimated reads (highest first)
    std::sort(work_units.begin(), work_units.end(),
              [](const WorkUnit& a, const WorkUnit& b){
                  return a.est_reads > b.est_reads;
              });
    
    fprintf(stderr, "DUAL-PANEL: BAM header has %d targets; processing %d SNP-bearing targets (%ld total combined SNPs) using %d threads...\n",
        n_chroms, chroms_with_snps, total_snps, n_threads);
    fprintf(stderr, "  Split into %lu work units (%d by SNP density, %d by read density)\n", 
        work_units.size(), chroms_chunked_by_snp, chroms_chunked_by_reads);
    
    // Progress tracking
    atomic<long> snps_processed(0);
    atomic<int> units_done(0);
    atomic<long> reads_processed(0);
    
    // Per-thread cell counts: separate maps for panel 0 and panel 1
    omp_set_num_threads(n_threads);
    vector<robin_hood::unordered_map<unsigned long, CellCounts>> thread_counts_p0(n_threads);
    vector<robin_hood::unordered_map<unsigned long, CellCounts>> thread_counts_p1(n_threads);
    vector<robin_hood::unordered_map<unsigned long, CellCounts>> thread_species_native_counts(
        collect_species_native ? n_threads : 0);

    // --dump_pileup: per-thread per-(cell,SNP) allele evidence (interindividual
    // only).  Inner key packs (tid<<32 | pos); value is (ref_scaled, alt_scaled).
    // Empty and untouched unless dump_pileup is set.  Identical structure and
    // downstream contract to the single-panel path.
    vector<robin_hood::unordered_map<unsigned long,
        robin_hood::unordered_map<int64_t, std::pair<int64_t, int64_t> > > > thread_pileup(n_threads);
    vector<AcceptedSiteWeightMap> thread_site_weights_p0(n_threads);
    vector<AcceptedSiteWeightMap> thread_site_weights_p1(n_threads);
    vector<SourceObservationMap> thread_source_observations(n_threads);
    vector<SourceCategoryMap> thread_source_raw_categories(n_threads);
    vector<DonorAuditMap> thread_donor_audit(n_threads);
    vector<DonorSiteAuditMap> thread_donor_site_audit(n_threads);
    std::atomic<uint64_t> source_reads_native_sentinel{0};
    std::atomic<uint64_t> source_reads_missing_yi{0};
    std::atomic<uint64_t> source_reads_invalid_yi{0};
    std::atomic<uint64_t> source_reads_unknown_yi{0};
    std::atomic<uint64_t> source_reads_missing_synthetic_id{0};
    std::atomic<uint64_t> source_reads_mismatched_synthetic_id{0};
    map<string, int> source_sample_to_idx;
    for (int i = 0; i < (int)sample_names.size(); ++i){
        source_sample_to_idx[sample_names[i]] = i;
    }
    const SourceReceiverMap source_receiver_map =
        load_source_receiver_map(source_receiver_map_path, source_sample_to_idx);
    std::atomic<uint64_t> source_observations_missing_receiver_map{0};
    std::atomic<uint64_t> source_observations_invalid_receiver_genotype{0};
    ParallelOperationStatus operation_status;
    std::atomic<bool> hts_thread_warning_emitted(false);
    
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        auto& local_p0 = thread_counts_p0[thread_id];
        auto& local_p1 = thread_counts_p1[thread_id];
        robin_hood::unordered_map<unsigned long, CellCounts>* local_species_native =
            collect_species_native ? &thread_species_native_counts[thread_id] : nullptr;
        auto& local_site_weights_p0 = thread_site_weights_p0[thread_id];
        auto& local_site_weights_p1 = thread_site_weights_p1[thread_id];
        auto& local_source_observations = thread_source_observations[thread_id];
        auto& local_source_raw_categories = thread_source_raw_categories[thread_id];
        auto& local_donor_audit = thread_donor_audit[thread_id];
        auto& local_donor_site_audit = thread_donor_site_audit[thread_id];
        
        htsFile* bam_fp = hts_open(bamfile.c_str(), "r");
        bam_hdr_t* header = nullptr;
        hts_idx_t* idx = nullptr;
        bam1_t* record = nullptr;
        if (!bam_fp){
            operation_status.fail(format_worker_error("BAM open", thread_id));
        }
        else{
            if (htslib_threads > 1 && hts_set_threads(bam_fp, htslib_threads) < 0){
                bool expected = false;
                if (hts_thread_warning_emitted.compare_exchange_strong(expected, true)){
                    fprintf(stderr,
                        "WARNING: HTSlib helper-thread setup failed; continuing with synchronous BAM I/O\n");
                }
            }
            header = sam_hdr_read(bam_fp);
            if (!header){
                operation_status.fail(format_worker_error("BAM header read", thread_id));
            }
            idx = sam_index_load(bam_fp, bamfile.c_str());
            if (!idx){
                operation_status.fail(format_worker_error("BAM index load", thread_id));
            }
            record = bam_init1();
            if (!record){
                operation_status.fail(format_worker_error("BAM record allocation", thread_id));
            }
        }

        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < work_units.size(); i++){
            if (!operation_status.ok() || !bam_fp || !header || !idx || !record) continue;
            WorkUnit& wu = work_units[i];
            int tid = wu.tid;
            if (tid < 0 || tid >= header->n_targets){
                operation_status.fail(format_worker_error("invalid contig", thread_id, tid));
                continue;
            }

            hts_itr_t* iter = sam_itr_queryi(idx, tid, wu.start_pos, wu.end_pos);
            if (!iter){
                operation_status.fail(format_worker_error("iterator creation", thread_id, tid));
                continue;
            }
                    
                    long local_snps = 0;
                    long local_reads = 0;
                    long local_all_reads = 0;
                    
                    if (!wu.has_snps){
                        int iterator_result = 0;
                        while ((iterator_result = sam_itr_next(bam_fp, iter, record)) >= 0){
                            if (!read_passes_filter(record, default_production_read_filter())){
                                continue;
                            }
                            local_all_reads++;
                        }
                        if (iterator_result < -1){
                            operation_status.fail(format_worker_error("iterator read", thread_id, tid));
                        }
                        hts_itr_destroy(iter);
                        reads_processed += local_all_reads;
                        
                        int done = ++units_done;
                        if (done % 100 == 0 || done == (int)work_units.size()){
                            fprintf(stderr, "\rDUAL progress: %d/%lu units, %ld/%ld SNPs, %ld iterator records",
                                done, work_units.size(), snps_processed.load(), total_snps,
                                reads_processed.load());
                        }
                        continue;
                    }
                    
                    auto snp_it = combined_snpdat.find(tid);
                    ChromSNPs& chrom_snps = snp_it->second;
                    const NativeSpeciesChromTargets* native_chrom_targets = nullptr;
                    if (collect_species_native){
                        auto native_it = species_native_targets->find(tid);
                        if (native_it != species_native_targets->end()){
                            native_chrom_targets = &native_it->second;
                        }
                    }
                    
                    auto snp_iter = chrom_snps.snps.begin() + wu.snp_start;
                    auto snp_chunk_end = chrom_snps.snps.begin() + wu.snp_end;
                    
                    int iterator_result = 0;
                    while ((iterator_result = sam_itr_next(bam_fp, iter, record)) >= 0){
                        if (!read_passes_filter(record, default_production_read_filter())){
                            continue;
                        }
                        
                        local_all_reads++;
                        
                        int read_start = record->core.pos;
                        int read_end = bam_endpos(record);
                        
                        while (snp_iter != snp_chunk_end && snp_iter->pos < read_start){
                            ++snp_iter;
                            ++local_snps;
                        }
                        
                        // Bin-skip: if no SNP bin overlaps this read, skip CB extraction
                        {
                            auto bin_it = bin_indices.find(tid);
                            if (bin_it != bin_indices.end() &&
                                !bin_it->second.might_overlap(read_start, read_end)){
                                continue;
                            }
                        }
                        
                        uint8_t* cb_tag = bam_aux_get(record, "CB");
                        if (!cb_tag) continue;
                        
                        const char* cb_str = bam_aux2Z(cb_tag);
                        bc cb_bits;
                        str2bc(cb_str, cb_bits);
                        unsigned long bc_key = cb_bits.to_ulong();
                        
                        if (has_bc_list && valid_barcodes.find(bc_key) == valid_barcodes.end()){
                            continue;
                        }
                        
                        local_reads++;

                        std::string source_label;
                        std::string source_origin;
                        bool source_label_resolvable = false;
                        if (dump_source_observations){
                            uint8_t* source_tag = bam_aux_get(record, source_provenance_tag.c_str());
                            if (!source_tag){
                                source_label = "__MISSING_YI__";
                                uint8_t* sid_tag = bam_aux_get(record, synthetic_id_tag.c_str());
                                const char* sid_z = sid_tag ? bam_aux2Z(sid_tag) : NULL;
                                if (!sid_tag){
                                    // Native reads do not carry YS. This salvages the existing
                                    // even BAM without modifying it; the sidecar is never a
                                    // production estimator input.
                                    source_origin = "native_untagged";
                                }
                                else if (sid_z && sid_z[0] != '\0' && expected_synthetic_id == sid_z){
                                    // Injected reads always carry YS in the benchmark contract.
                                    // Their donor identity is unavailable, but their injected
                                    // origin is still known for evaluator-only rate truth.
                                    source_origin = "injected_missing_yi";
                                }
                                else{
                                    // Present-but-empty, non-string, or wrong-unit YS is malformed
                                    // provenance and remains unclassified instead of becoming native.
                                    source_origin = "missing_yi";
                                }
                                source_reads_missing_yi.fetch_add(1, std::memory_order_relaxed);
                            }
                            else{
                                const char* source_z = bam_aux2Z(source_tag);
                                if (!source_z){
                                    source_origin = "invalid_yi";
                                    source_label = "__INVALID_YI__";
                                    source_reads_invalid_yi.fetch_add(1, std::memory_order_relaxed);
                                }
                                else if (source_z[0] == '\0'){
                                    source_label = "__EMPTY_YI__";
                                    uint8_t* sid_tag = bam_aux_get(record, synthetic_id_tag.c_str());
                                    const char* sid_z = sid_tag ? bam_aux2Z(sid_tag) : NULL;
                                    if (!sid_tag){
                                        source_origin = "native_untagged";
                                    }
                                    else if (sid_z && sid_z[0] != '\0' && expected_synthetic_id == sid_z){
                                        source_origin = "injected_missing_yi";
                                    }
                                    else{
                                        source_origin = "missing_yi";
                                    }
                                    source_reads_missing_yi.fetch_add(1, std::memory_order_relaxed);
                                }
                                else if (strcmp(source_z, "__NATIVE__") == 0){
                                    uint8_t* sid_tag = bam_aux_get(record, synthetic_id_tag.c_str());
                                    if (sid_tag){
                                        // Native provenance and a synthetic-unit injection tag are
                                        // mutually inconsistent. Keep the read explicit but
                                        // unclassified for evaluator-side truth.
                                        source_origin = "invalid_yi";
                                        source_label = "__INVALID_YI__";
                                        source_reads_invalid_yi.fetch_add(1, std::memory_order_relaxed);
                                    }
                                    else{
                                        source_origin = "native";
                                        source_label.assign(source_z);
                                        source_label_resolvable = true;
                                        source_reads_native_sentinel.fetch_add(
                                            1, std::memory_order_relaxed);
                                    }
                                }
                                else{
                                    source_label.assign(source_z);
                                    source_label_resolvable = source_label_is_known(
                                        source_label, source_sample_to_idx);
                                    uint8_t* sid_tag = bam_aux_get(record, synthetic_id_tag.c_str());
                                    const char* sid_z = sid_tag ? bam_aux2Z(sid_tag) : NULL;
                                    if (!sid_z || sid_z[0] == '\0'){
                                        source_origin = "injected_missing_ys";
                                        source_reads_missing_synthetic_id.fetch_add(1, std::memory_order_relaxed);
                                    }
                                    else if (expected_synthetic_id != sid_z){
                                        source_origin = "injected_wrong_ys";
                                        source_reads_mismatched_synthetic_id.fetch_add(1, std::memory_order_relaxed);
                                    }
                                    else if (!source_label_resolvable){
                                        source_origin = "unknown_yi";
                                        source_reads_unknown_yi.fetch_add(1, std::memory_order_relaxed);
                                    }
                                    else{
                                        source_origin = "injected_expected_ys";
                                    }
                                }
                            }
                        }
                        std::set<SourceObservationKey> source_read_keys;
                        
                        float prob_correct = 1.0f - powf(10.0f, -(float)record->core.qual / 10.0f);
                        int64_t prob_scaled = (int64_t)(prob_correct * FIXED_POINT_SCALE);
                        
                        for (auto snp_check = snp_iter; 
                             snp_check != snp_chunk_end && snp_check->pos < read_end; 
                             ++snp_check){
                            
                            char allele = get_base_at_pos(record, snp_check->pos);
                            if (allele == 'N' || allele == '-') continue;
                            
                            int64_t ref_add = 0, alt_add = 0;
                            
                            if (allele == snp_check->data.ref){
                                ref_add = prob_scaled;
                            }
                            else if (allele == snp_check->data.alt){
                                alt_add = prob_scaled;
                            }
                            else {
                            }
                            
                            if (ref_add > 0 || alt_add > 0){
                                const int panel_id = snp_check->panel_id == 0 ? 0 : 1;
                                if (panel_id == 0 && accepted_site_weights_panel0 != nullptr){
                                    local_site_weights_p0[accepted_site_weight_key(tid, snp_check->pos)] +=
                                        ref_add + alt_add;
                                }
                                else if (panel_id == 1 && accepted_site_weights_panel1 != nullptr){
                                    local_site_weights_p1[accepted_site_weight_key(tid, snp_check->pos)] +=
                                        ref_add + alt_add;
                                }
                                std::string receiver_identity;
                                std::string receiver_a;
                                std::string receiver_b;
                                int receiver_nalt_a = -1;
                                int receiver_nalt_b = -1;
                                bool receiver_category_valid = false;
                                if (panel_id == 0 && !source_receiver_map.empty()){
                                    auto receiver_it = source_receiver_map.find(bc_key);
                                    if (receiver_it == source_receiver_map.end()){
                                        source_observations_missing_receiver_map.fetch_add(
                                            1, std::memory_order_relaxed);
                                    }
                                    else{
                                        const SourceReceiverInfo& receiver = receiver_it->second;
                                        receiver_identity = receiver.identity;
                                        receiver_a = receiver.name_a;
                                        receiver_b = receiver.name_b;
                                        if (receiver.idx_a >= 0 && receiver.idx_b >= 0 &&
                                            receiver.idx_a < (int)snp_check->geno.size() &&
                                            receiver.idx_b < (int)snp_check->geno.size()){
                                            receiver_nalt_a = (int)snp_check->geno[receiver.idx_a];
                                            receiver_nalt_b = (int)snp_check->geno[receiver.idx_b];
                                            receiver_category_valid = receiver_nalt_a >= 0 && receiver_nalt_b >= 0;
                                        }
                                        else{
                                            source_observations_invalid_receiver_genotype.fetch_add(
                                                1, std::memory_order_relaxed);
                                        }
                                    }
                                }

                                if (dump_source_observations && source_reconciliation_mode &&
                                    panel_id == 0 && receiver_category_valid){
                                    SourceCategoryKey raw_key(
                                        bc_key, receiver_identity, receiver_a, receiver_b,
                                        receiver_nalt_a, receiver_nalt_b, panel_id);
                                    auto& raw_stats = local_source_raw_categories[raw_key];
                                    raw_stats.n_observations += 1;
                                    raw_stats.ref_weight += (double)ref_add / (double)FIXED_POINT_SCALE;
                                    raw_stats.alt_weight += (double)alt_add / (double)FIXED_POINT_SCALE;
                                }

                                // Record every accepted observation when a source sidecar is
                                // requested, including missing/invalid provenance.  The evaluator
                                // needs those explicit buckets to determine whether rate truth is
                                // complete; silently dropping them recreates the even native-YI
                                // failure.  This sidecar is never consumed by production fitting.
                                const bool record_source = dump_source_observations;
                                if (record_source){
                                    std::string resolved_source = source_label;
                                    std::string typing_status = "unavailable";
                                    if (source_label_resolvable){
                                        resolve_source_observation(
                                            source_label, *snp_check, ref_add > 0, source_sample_to_idx,
                                            resolved_source, typing_status);
                                    }
                                    SourceObservationKey source_key(
                                        bc_key, receiver_identity, receiver_a, receiver_b,
                                        receiver_nalt_a, receiver_nalt_b,
                                        source_origin, resolved_source, source_label, typing_status, panel_id);
                                    auto& source_stats = local_source_observations[source_key];
                                    source_stats.n_observations += 1;
                                    source_stats.ref_weight += (double)ref_add / (double)FIXED_POINT_SCALE;
                                    source_stats.alt_weight += (double)alt_add / (double)FIXED_POINT_SCALE;
                                    source_read_keys.insert(source_key);

                                    if (source_donor_site_audit && panel_id == 0 &&
                                        receiver_category_valid && source_origin == "injected_expected_ys"){
                                        std::string source_a;
                                        std::string source_b;
                                        int source_nalt_a = -1;
                                        int source_nalt_b = -1;
                                        int resolved_nalt = -1;
                                        describe_source_genotypes(
                                            source_label, resolved_source, *snp_check, source_sample_to_idx,
                                            source_a, source_b, source_nalt_a, source_nalt_b, resolved_nalt);
                                        DonorAuditKey donor_key(
                                            receiver_identity, receiver_a, receiver_b,
                                            receiver_nalt_a, receiver_nalt_b,
                                            source_label, resolved_source, typing_status,
                                            source_a, source_b, source_nalt_a, source_nalt_b, resolved_nalt);
                                        const bool is_ref = ref_add > 0;
                                        const double audit_weight = (double)(ref_add + alt_add) /
                                            (double)FIXED_POINT_SCALE;
                                        update_donor_audit_stats(
                                            local_donor_audit[donor_key], is_ref, audit_weight,
                                            source_nalt_a, source_nalt_b, resolved_nalt);
                                        if (donor_site_is_sampled(tid, snp_check->pos,
                                                source_donor_site_sample_mod)){
                                            DonorSiteAuditKey site_key(
                                                tid, snp_check->pos, snp_check->data.ref, snp_check->data.alt,
                                                receiver_identity, receiver_a, receiver_b,
                                                receiver_nalt_a, receiver_nalt_b,
                                                source_label, resolved_source, typing_status,
                                                source_a, source_b, source_nalt_a, source_nalt_b, resolved_nalt);
                                            update_donor_audit_stats(
                                                local_donor_site_audit[site_key], is_ref, audit_weight,
                                                source_nalt_a, source_nalt_b, resolved_nalt);
                                        }
                                    }
                                }
                                // Route to panel-specific thread-local map based on panel_id
                                auto& target_counts = (snp_check->panel_id == 0) ? local_p0 : local_p1;
                                
                                auto it = target_counts.find(bc_key);
                                if (it == target_counts.end()){
                                    target_counts.emplace(bc_key, CellCounts(n_samples));
                                    it = target_counts.find(bc_key);
                                }
                                CellCounts& cc = it->second;
                                
                                // Precomputed targets: linear traversal, no branches
                                const auto& ttargets = snp_check->total_targets;
                                const auto& ptargets = snp_check->pair_targets;
                                
                                for (const auto& t : ttargets){
                                    cc.total_ref[t.total_idx] += ref_add;
                                    cc.total_alt[t.total_idx] += alt_add;
                                }
                                for (const auto& p : ptargets){
                                    cc.ref_counts[p.pair_idx] += ref_add;
                                    cc.alt_counts[p.pair_idx] += alt_add;
                                }

                                // --dump_pileup: record per-(cell,SNP) evidence for
                                // interindividual SNPs.  Summed within a thread; the
                                // producer sums any cross-thread duplicates.
                                if (dump_pileup && snp_check->panel_id == 0){
                                    int64_t pkey = ((int64_t)tid << 32) |
                                        (int64_t)(uint32_t)snp_check->pos;
                                    auto& slot = thread_pileup[thread_id][bc_key][pkey];
                                    slot.first += ref_add;
                                    slot.second += alt_add;
                                }

                                if (local_species_native != nullptr && snp_check->panel_id == 1 &&
                                    native_chrom_targets != nullptr){
                                    const size_t snp_index = (size_t)(
                                        snp_check - chrom_snps.snps.begin());
                                    if (snp_index < native_chrom_targets->site_offsets.size() &&
                                        native_chrom_targets->site_offsets[snp_index] != UINT64_MAX){
                                        auto native_it = local_species_native->find(bc_key);
                                        if (native_it == local_species_native->end()){
                                            local_species_native->emplace(
                                                bc_key, CellCounts(species_native_n_samples));
                                            native_it = local_species_native->find(bc_key);
                                        }
                                        accumulate_species_native_targets(
                                            native_it->second, *native_chrom_targets, snp_index,
                                            ref_add, alt_add);
                                    }
                                }
                            }
                        }
                        if (dump_source_observations){
                            for (const auto& source_key : source_read_keys){
                                local_source_observations[source_key].n_reads += 1;
                            }
                        }
                    }
                    
                    if (iterator_result < -1){
                        operation_status.fail(format_worker_error("iterator read", thread_id, tid));
                    }

                    while (snp_iter != snp_chunk_end){
                        ++snp_iter;
                        ++local_snps;
                    }
                    
                    snps_processed += local_snps;
                    reads_processed += local_all_reads;
                    int done = ++units_done;
                    
                    if (done % 10 == 0 || done == (int)work_units.size()){
                        fprintf(stderr, "\rDUAL progress: %d/%lu units, %ld/%ld SNPs, %ld iterator records          ",
                            done, work_units.size(), snps_processed.load(), total_snps,
                            reads_processed.load());
                    }
                    
                    hts_itr_destroy(iter);
        }

        if (record) bam_destroy1(record);
        if (idx) hts_idx_destroy(idx);
        if (header) bam_hdr_destroy(header);
        if (bam_fp) hts_close(bam_fp);
    }

    if (!operation_status.ok()){
        fprintf(stderr, "ERROR: dual-panel allele counting failed: %s\n",
            operation_status.message().c_str());
        return false;
    }
    
    // Merge per-thread counts for panel 0
    fprintf(stderr, "\nMerging per-thread counts (panel 0)...\n");
    for (int t = 0; t < n_threads; t++){
        for (auto& kv : thread_counts_p0[t]){
            unsigned long bc_key = kv.first;
            auto it = counts_panel0.find(bc_key);
            if (it == counts_panel0.end()){
                counts_panel0.emplace(std::piecewise_construct,
                    std::forward_as_tuple(bc_key),
                    std::forward_as_tuple(n_samples));
                it = counts_panel0.find(bc_key);
            }
            it->second.counts.merge(kv.second);
        }
        thread_counts_p0[t].clear();
    }
    
    // Merge per-thread counts for panel 1
    fprintf(stderr, "Merging per-thread counts (panel 1)...\n");
    for (int t = 0; t < n_threads; t++){
        for (auto& kv : thread_counts_p1[t]){
            unsigned long bc_key = kv.first;
            auto it = counts_panel1.find(bc_key);
            if (it == counts_panel1.end()){
                counts_panel1.emplace(std::piecewise_construct,
                    std::forward_as_tuple(bc_key),
                    std::forward_as_tuple(n_samples));
                it = counts_panel1.find(bc_key);
            }
            it->second.counts.merge(kv.second);
        }
        thread_counts_p1[t].clear();
    }
    
    if (collect_species_native){
        fprintf(stderr, "Merging per-thread native species counts (panel 1)...\n");
        for (int t = 0; t < n_threads; ++t){
            for (auto& kv : thread_species_native_counts[t]){
                const unsigned long bc_key = kv.first;
                auto it = species_native_counts->find(bc_key);
                if (it == species_native_counts->end()){
                    species_native_counts->emplace(std::piecewise_construct,
                        std::forward_as_tuple(bc_key),
                        std::forward_as_tuple(species_native_n_samples));
                    it = species_native_counts->find(bc_key);
                }
                it->second.counts.merge(kv.second);
            }
            thread_species_native_counts[t].clear();
        }
        thread_species_native_counts.clear();
        thread_species_native_counts.shrink_to_fit();
        fprintf(stderr, "Native species counts (panel 1): %lu cells, %d species\n",
            species_native_counts->size(), species_native_n_samples);
    }

    thread_counts_p0.clear();
    thread_counts_p0.shrink_to_fit();
    thread_counts_p1.clear();
    thread_counts_p1.shrink_to_fit();

    // --dump_pileup: flush per-thread per-(cell,SNP) observations.  Per-thread
    // rows are pre-summed; the same (bc,tid,pos) may still appear across threads
    // and is summed downstream by the producer.  Columns: bc_hash tid pos ref alt
    // (ref/alt are prob-scaled float evidence, matching the .counts units).
    if (dump_pileup){
        string obs_path = pileup_prefix + ".pileup_obs.tsv.gz";
        gzFile pf = gzopen(obs_path.c_str(), "w");
        if (!pf){
            fprintf(stderr, "ERROR: could not open %s for writing\n", obs_path.c_str());
            return false;
        } else {
            long n_obs_written = 0;
            for (int t = 0; t < n_threads; t++){
                for (auto& cell : thread_pileup[t]){
                    unsigned long bc = cell.first;
                    for (auto& kv : cell.second){
                        int tid_o = (int)(kv.first >> 32);
                        int pos_o = (int)(kv.first & 0xFFFFFFFF);
                        gzprintf(pf, "%lu\t%d\t%d\t%f\t%f\n", bc, tid_o, pos_o,
                            (double)kv.second.first / FIXED_POINT_SCALE,
                            (double)kv.second.second / FIXED_POINT_SCALE);
                        n_obs_written++;
                    }
                }
                thread_pileup[t].clear();
            }
            if (gzclose(pf) != Z_OK){
                fprintf(stderr, "ERROR: failed while closing %s\n", obs_path.c_str());
                return false;
            }
            fprintf(stderr, "Wrote %ld pileup observations to %s\n", n_obs_written, obs_path.c_str());
        }
    }
    thread_pileup.clear();
    thread_pileup.shrink_to_fit();

    auto merge_site_weights = [n_threads](
        vector<AcceptedSiteWeightMap>& per_thread,
        AcceptedSiteWeightMap* destination,
        const char* panel_name){
        if (destination != nullptr){
            destination->clear();
            for (int t = 0; t < n_threads; ++t){
                for (const auto& kv : per_thread[t]){
                    (*destination)[kv.first] += kv.second;
                }
                per_thread[t].clear();
            }
            fprintf(stderr, "Accepted-site weight map (%s): %lu observed sites\n",
                panel_name, destination->size());
        }
        per_thread.clear();
        per_thread.shrink_to_fit();
    };
    merge_site_weights(thread_site_weights_p0, accepted_site_weights_panel0, "individual");
    merge_site_weights(thread_site_weights_p1, accepted_site_weights_panel1, "species");

    if (dump_source_observations){
        const uint64_t native_sentinel = source_reads_native_sentinel.load();
        if (native_sentinel > 0){
            fprintf(stderr,
                "Source-provenance validation: retained %llu native sentinel reads "
                "(YI='__NATIVE__') without requiring %s; all non-native YI reads "
                "must match %s='%s'.\n",
                (unsigned long long)native_sentinel, synthetic_id_tag.c_str(),
                synthetic_id_tag.c_str(), expected_synthetic_id.c_str());
        }
        const uint64_t missing_yi = source_reads_missing_yi.load();
        const uint64_t invalid_yi = source_reads_invalid_yi.load();
        const uint64_t unknown_yi = source_reads_unknown_yi.load();
        const uint64_t missing_sid = source_reads_missing_synthetic_id.load();
        const uint64_t mismatched_sid = source_reads_mismatched_synthetic_id.load();
        const uint64_t missing_receiver = source_observations_missing_receiver_map.load();
        const uint64_t invalid_receiver_genotype = source_observations_invalid_receiver_genotype.load();
        if (source_reconciliation_mode){
            fprintf(stderr,
                "Source provenance read buckets: native=%llu missing_yi=%llu invalid_yi=%llu "
                "unknown_yi=%llu injected_missing_ys=%llu injected_wrong_ys=%llu.\n",
                (unsigned long long)native_sentinel,
                (unsigned long long)missing_yi,
                (unsigned long long)invalid_yi,
                (unsigned long long)unknown_yi,
                (unsigned long long)missing_sid,
                (unsigned long long)mismatched_sid);
        }
        else if (missing_sid > 0 || mismatched_sid > 0){
            fprintf(stderr, "ERROR: source-provenance validation failed: %llu YI-tagged reads missing %s and %llu carrying a non-matching %s (expected '%s').\n",
                (unsigned long long)missing_sid, synthetic_id_tag.c_str(),
                (unsigned long long)mismatched_sid, synthetic_id_tag.c_str(),
                expected_synthetic_id.c_str());
            return false;
        }
        if (!source_receiver_map_path.empty() &&
            (missing_receiver > 0 || invalid_receiver_genotype > 0)){
            fprintf(stderr,
                "ERROR: category-resolved source provenance failed: %llu accepted individual-panel observations lacked a receiver-map entry and %llu had an invalid receiver genotype index.\n",
                (unsigned long long)missing_receiver,
                (unsigned long long)invalid_receiver_genotype);
            return false;
        }
        if (!source_receiver_map_path.empty()){
            fprintf(stderr,
                "Category-resolved source provenance enabled with %lu receiver mappings; individual-panel rows carry authored receiver (nalt_A,nalt_B).\n",
                (unsigned long)source_receiver_map.size());
        }
        SourceObservationMap merged_source_observations;
        SourceCategoryMap merged_source_raw_categories;
        for (int t = 0; t < n_threads; t++){
            for (const auto& kv : thread_source_observations[t]){
                merged_source_observations[kv.first].merge(kv.second);
            }
            for (const auto& kv : thread_source_raw_categories[t]){
                merged_source_raw_categories[kv.first].merge(kv.second);
            }
            thread_source_observations[t].clear();
            thread_source_raw_categories[t].clear();
        }
        write_source_observation_summary(source_observations_prefix, merged_source_observations);
        if (source_reconciliation_mode){
            write_source_reconciliation_summary(
                source_observations_prefix,
                merged_source_observations,
                merged_source_raw_categories);
        }
        if (source_donor_site_audit){
            DonorAuditMap merged_donor_audit;
            DonorSiteAuditMap merged_donor_site_audit;
            for (int t = 0; t < n_threads; ++t){
                for (const auto& kv : thread_donor_audit[t]){
                    merged_donor_audit[kv.first].merge(kv.second);
                }
                for (const auto& kv : thread_donor_site_audit[t]){
                    merged_donor_site_audit[kv.first].merge(kv.second);
                }
                thread_donor_audit[t].clear();
                thread_donor_site_audit[t].clear();
            }
            write_donor_genotype_audit(source_observations_prefix, merged_donor_audit);
            write_donor_site_sample(
                source_observations_prefix, merged_donor_site_audit,
                chrom_names, source_donor_site_sample_mod);
        }
    }
    thread_source_observations.clear();
    thread_source_observations.shrink_to_fit();
    thread_source_raw_categories.clear();
    thread_source_raw_categories.shrink_to_fit();
    thread_donor_audit.clear();
    thread_donor_audit.shrink_to_fit();
    thread_donor_site_audit.clear();
    thread_donor_site_audit.shrink_to_fit();
    
    fprintf(stderr, "DUAL completed: %d chromosomes (%lu work units), %ld SNPs, %ld iterator records\n",
        n_chroms, work_units.size(), snps_processed.load(), reads_processed.load());
    fprintf(stderr, "  Panel 0 (interindiv): %lu cells\n", counts_panel0.size());
    fprintf(stderr, "  Panel 1 (individual-shaped species evidence): %lu cells\n", counts_panel1.size());
    if (collect_species_native){
        fprintf(stderr, "  Panel 1 (native species evidence):            %lu cells\n",
            species_native_counts->size());
    }
    return true;
}

bool count_alleles_single_threaded(
    const string& bamfile,
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_all,
    robin_hood::unordered_map<unsigned long, CellCounts>& cell_counts,
    const set<unsigned long>& valid_barcodes,
    int n_samples,
    map<pair<int, int>, map<int, float> >& conditional_match_fracs,
    map<pair<int, int>, map<int, float> >& conditional_match_tots,
    bool compute_conditional,
    const NativeSpeciesTargetTable* species_native_targets,
    robin_hood::unordered_map<unsigned long, CellCounts>* species_native_counts,
    int species_native_n_samples){
    robin_hood::unordered_map<unsigned long, AlignedCellCounts> parallel_counts;
    robin_hood::unordered_map<unsigned long, AlignedCellCounts> parallel_species_counts;
    const bool collect_species_native =
        species_native_targets != nullptr && species_native_counts != nullptr &&
        species_native_n_samples > 0;

    if (!count_alleles_parallel(
            bamfile, snpdat_all, parallel_counts, valid_barcodes, n_samples,
            1, 1, false, "", nullptr, species_native_targets,
            collect_species_native ? &parallel_species_counts : nullptr,
            species_native_n_samples)){
        return false;
    }

    finalize_parallel_counts(parallel_counts, cell_counts);
    if (collect_species_native){
        finalize_parallel_counts(parallel_species_counts, *species_native_counts);
    }
    if (compute_conditional){
        conditional_match_tots.clear();
        compute_conditional_match_fracs_parallel(
            snpdat_all, conditional_match_fracs, n_samples, 1);
    }
    return true;
}

void finalize_parallel_counts(
    robin_hood::unordered_map<unsigned long, AlignedCellCounts>& parallel_counts,
    robin_hood::unordered_map<unsigned long, CellCounts>& final_counts){
    
    final_counts.clear();
    
    for (auto& kv : parallel_counts){
        final_counts.emplace(kv.first, std::move(kv.second.counts));
    }
}

// ============================================================================
// SHARED MEMORY VCF FUNCTIONS
// ============================================================================

bool create_shared_vcf(
    const string& vcf_file,
    const string& shm_name,
    set<string>& chroms_to_include,
    map<string, int>& seq2tid,
    int min_vq){
    
    // Load VCF data
    robin_hood::unordered_map<int, ChromSNPs> snpdat;
    string vcf_file_copy = vcf_file;  // read_vcf_chroms_optimized needs non-const
    int nvar = read_vcf_chroms_optimized(vcf_file_copy, chroms_to_include, seq2tid, snpdat, min_vq);
    
    if (nvar <= 0){
        fprintf(stderr, "ERROR: %s variants loaded from VCF\n",
            nvar < 0 ? "failed to load" : "no usable");
        return false;
    }
    
    // Read sample names from VCF
    vector<string> samples;
    read_vcf_samples(vcf_file_copy, samples);
    fprintf(stderr, "Found %lu samples in VCF\n", samples.size());
    
    // Calculate total size needed
    size_t total_size = sizeof(SharedVCFHeader);
    int n_chroms = 0;
    
    for (auto& kv : snpdat){
        total_size += kv.second.snps.size() * sizeof(SNPData);
        n_chroms++;
    }
    
    fprintf(stderr, "Creating shared memory segment: %s (%.2f GB, %d SNPs, %d chroms)\n",
        shm_name.c_str(), (double)total_size / (1024.0 * 1024.0 * 1024.0), nvar, n_chroms);
    
    // Create shared memory
    int shm_fd = shm_open(shm_name.c_str(), O_CREAT | O_RDWR, 0666);
    if (shm_fd < 0){
        perror("shm_open");
        return false;
    }
    
    if (ftruncate(shm_fd, total_size) < 0){
        perror("ftruncate");
        close(shm_fd);
        shm_unlink(shm_name.c_str());
        return false;
    }
    
    void* ptr = mmap(NULL, total_size, PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0);
    if (ptr == MAP_FAILED){
        perror("mmap");
        close(shm_fd);
        shm_unlink(shm_name.c_str());
        return false;
    }
    
    // Serialize VCF data
    SharedVCFHeader* header = (SharedVCFHeader*)ptr;
    header->total_size = total_size;
    header->n_chromosomes = n_chroms;
    header->n_snps_total = nvar;
    
    // Store sample names
    header->n_samples = samples.size();
    for (size_t i = 0; i < samples.size() && i < 512; i++){
        strncpy(header->sample_names[i], samples[i].c_str(), 63);
        header->sample_names[i][63] = '\0';
    }
    fprintf(stderr, "Stored %d sample names in shared memory\n", header->n_samples);
    
    char* data_ptr = (char*)ptr + sizeof(SharedVCFHeader);
    size_t offset = sizeof(SharedVCFHeader);
    int chrom_idx = 0;
    
    for (auto& kv : snpdat){
        header->chrom_offsets[chrom_idx] = offset;
        header->chrom_snp_counts[chrom_idx] = kv.second.snps.size();
        header->chrom_tids[chrom_idx] = kv.first;
        
        // Serialize SNP data field-by-field into shared memory.
        // SNPData contains var::gqs (a std::vector<float>) which is
        // non-trivially-copyable. A raw memcpy would write the vector's
        // internal heap pointer into shared memory, causing UB when the
        // reader process later interprets those bytes. Instead, we zero
        // each slot first (clearing the gqs region to a null/zero state),
        // then copy only the POD fields.
        for (size_t s = 0; s < kv.second.snps.size(); s++){
            SNPData* dest = (SNPData*)(data_ptr + s * sizeof(SNPData));
            // Zero the entire slot so the gqs vector bytes are null
            memset(dest, 0, sizeof(SNPData));
            // Copy POD fields only
            const SNPData& src = kv.second.snps[s];
            dest->pos = src.pos;
            dest->panel_id = src.panel_id;
            dest->data.ref = src.data.ref;
            dest->data.alt = src.data.alt;
            dest->data.haps1 = src.data.haps1;
            dest->data.haps2 = src.data.haps2;
            dest->data.haps_covered = src.data.haps_covered;
            dest->data.vq = src.data.vq;
            // data.gqs left as zeroed bytes (empty vector representation)
        }
        size_t copy_size = kv.second.snps.size() * sizeof(SNPData);
        data_ptr += copy_size;
        offset += copy_size;
        chrom_idx++;
    }
    
    // Sync to ensure data is written
    msync(ptr, total_size, MS_SYNC);
    
    fprintf(stderr, "Shared memory created successfully\n");
    
    // Keep mapping but close fd (other processes will open separately)
    close(shm_fd);
    
    return true;
}

bool attach_shared_vcf(
    const string& shm_name,
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_all,
    vector<string>& samples){
    
    int shm_fd = shm_open(shm_name.c_str(), O_RDONLY, 0);
    if (shm_fd < 0){
        perror("shm_open");
        return false;
    }
    
    struct stat sb;
    if (fstat(shm_fd, &sb) < 0){
        perror("fstat");
        close(shm_fd);
        return false;
    }
    
    void* ptr = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, shm_fd, 0);
    if (ptr == MAP_FAILED){
        perror("mmap");
        close(shm_fd);
        return false;
    }
    
    SharedVCFHeader* header = (SharedVCFHeader*)ptr;
    
    fprintf(stderr, "Attached to shared VCF: %s (%d SNPs, %d chromosomes)\n",
        shm_name.c_str(), header->n_snps_total, header->n_chromosomes);
    
    // Retrieve sample names
    samples.clear();
    for (int i = 0; i < header->n_samples && i < 512; i++){
        samples.push_back(string(header->sample_names[i]));
    }
    fprintf(stderr, "Loaded %lu sample names from shared VCF\n", samples.size());
    
    // Deserialize
    snpdat_all.clear();
    
    for (int i = 0; i < header->n_chromosomes; i++){
        int tid = header->chrom_tids[i];
        size_t n_snps = header->chrom_snp_counts[i];
        
        SNPData* snp_ptr = (SNPData*)((char*)ptr + header->chrom_offsets[i]);
        
        ChromSNPs& cs = snpdat_all[tid];
        cs.snps.reserve(n_snps);
        
        for (size_t j = 0; j < n_snps; j++){
            // Cannot use push_back(snp_ptr[j]) because the SNPData in shared
            // memory was written via memcpy, and the var struct contains a
            // std::vector<float> gqs whose internal heap pointer is from the
            // daemon's address space. Copying it would segfault.
            // Instead, construct each SNPData from the safe POD fields only.
            SNPData sd;
            sd.pos = snp_ptr[j].pos;
            sd.panel_id = snp_ptr[j].panel_id;
            sd.data.ref = snp_ptr[j].data.ref;
            sd.data.alt = snp_ptr[j].data.alt;
            sd.data.haps1 = snp_ptr[j].data.haps1;
            sd.data.haps2 = snp_ptr[j].data.haps2;
            sd.data.haps_covered = snp_ptr[j].data.haps_covered;
            sd.data.vq = snp_ptr[j].data.vq;
            // gqs is left empty (default-constructed); it's not needed
            // for counting or conditional match fraction computation.
            cs.snps.push_back(std::move(sd));
        }
    }
    
    // Note: We keep the mapping active - caller should call detach when done
    close(shm_fd);
    
    return true;
}

void detach_shared_vcf(const string& shm_name){
    // The mmap is automatically unmapped when the process exits
    // This function is a placeholder for explicit cleanup if needed
}

void destroy_shared_vcf(const string& shm_name){
    if (shm_unlink(shm_name.c_str()) < 0){
        perror("shm_unlink");
    }
    else{
        fprintf(stderr, "Shared memory destroyed: %s\n", shm_name.c_str());
    }
}

// ============================================================================
// HET VCF FUNCTIONS (NEW - for ploidy detection)
// ============================================================================

int load_het_vcf(
    const string& het_vcf_file,
    const set<string>& chroms_to_include,
    map<string, int>& seq2tid,
    robin_hood::unordered_map<int, ChromSNPs>& het_snpdat,
    int min_vq){
    
    // This is essentially read_vcf_chroms_optimized but for het VCF
    // We can reuse that function since the format is the same
    string vcf_file_copy = het_vcf_file;
    set<string> chroms_copy = chroms_to_include;
    
    return read_vcf_chroms_optimized(vcf_file_copy, chroms_copy, seq2tid, het_snpdat, min_vq);
}

bool count_het_alleles_parallel(
    const string& bamfile,
    robin_hood::unordered_map<int, ChromSNPs>& het_snpdat,
    robin_hood::unordered_map<unsigned long, CellCounts>& het_counts,
    const set<unsigned long>& valid_barcodes,
    int n_samples,
    int n_threads,
    int htslib_threads){
    
    // This is essentially count_alleles_parallel but with a different output structure
    // We can use the same parallel counting infrastructure
    
    robin_hood::unordered_map<unsigned long, AlignedCellCounts> parallel_counts;
    
    // Use the same counting function - it works with any SNP set
    if (!count_alleles_parallel(bamfile, het_snpdat, parallel_counts,
            valid_barcodes, n_samples, n_threads, htslib_threads)){
        return false;
    }
    
    // Finalize into the output structure
    finalize_parallel_counts(parallel_counts, het_counts);
    
    fprintf(stderr, "Het allele counting complete: %lu cells\n", het_counts.size());
    return true;
}

/**
 * Extended het allele counting that collects per-site data and/or Welford stats.
 * 
 * Unlike count_het_alleles_parallel which aggregates by genotype pair (losing per-site info),
 * this preserves per-site information needed for accurate het balance variance.
 */
bool count_het_alleles_extended(
    const string& bamfile,
    robin_hood::unordered_map<int, ChromSNPs>& het_snpdat,
    robin_hood::unordered_map<unsigned long, CellHetData>& het_data,
    vector<pair<int, int>>& idx_to_site,
    const set<unsigned long>& valid_barcodes,
    int n_samples,
    int n_threads,
    int htslib_threads,
    HetBalanceMethod method) {
    size_t ignored_bytes = 0;
    string request_error;
    if (!validate_identity_and_allocation_request(
            n_samples, nullptr, &ignored_bytes, &request_error)){
        fprintf(stderr, "ERROR: invalid het/ploidy identity universe: %s\n",
            request_error.c_str());
        return false;
    }
    if (n_threads < 1 || htslib_threads < 1){
        fprintf(stderr, "ERROR: thread counts must be positive\n");
        return false;
    }

    const bool collect_persite = (method == HetBalanceMethod::PERSITE);
    const bool collect_welford = (method == HetBalanceMethod::WELFORD);
    const bool has_bc_list = !valid_barcodes.empty();
    if (!collect_persite && !collect_welford){
        fprintf(stderr, "ERROR: unsupported het-balance method\n");
        return false;
    }

    map<pair<int, int>, int32_t> site_to_idx;
    idx_to_site.clear();
    size_t global_idx = 0;
    for (const auto& kv : het_snpdat){
        const int tid = kv.first;
        for (const auto& snp : kv.second.snps){
            if (global_idx > (size_t)std::numeric_limits<int32_t>::max()){
                fprintf(stderr, "ERROR: het site index exceeds int32_t representation\n");
                return false;
            }
            site_to_idx[make_pair(tid, snp.pos)] = (int32_t)global_idx;
            idx_to_site.push_back(make_pair(tid, snp.pos));
            ++global_idx;
        }
    }

    const char* method_name = collect_persite ? "per-site" : "Welford";
    fprintf(stderr, "Processing %lu het sites with %s method\n",
        (unsigned long)global_idx, method_name);

    if (has_bc_list){
        fprintf(stderr, "Pre-allocating het data for %lu cells...\n", valid_barcodes.size());
        try {
            for (unsigned long bc : valid_barcodes){
                het_data.emplace(piecewise_construct,
                    forward_as_tuple(bc), forward_as_tuple(n_samples));
            }
        }
        catch (const std::exception& e){
            fprintf(stderr, "ERROR: het/ploidy pre-allocation failed: %s\n", e.what());
            return false;
        }
    }

    htsFile* bam_tmp = hts_open(bamfile.c_str(), "r");
    if (!bam_tmp){
        fprintf(stderr, "ERROR: Could not open BAM file: %s\n", bamfile.c_str());
        return false;
    }
    bam_hdr_t* hdr_tmp = sam_hdr_read(bam_tmp);
    if (!hdr_tmp){
        fprintf(stderr, "ERROR: Could not read BAM header: %s\n", bamfile.c_str());
        hts_close(bam_tmp);
        return false;
    }
    hts_idx_t* idx_tmp = sam_index_load(bam_tmp, bamfile.c_str());
    if (!idx_tmp){
        fprintf(stderr, "ERROR: Could not load required BAM index: %s\n", bamfile.c_str());
        bam_hdr_destroy(hdr_tmp);
        hts_close(bam_tmp);
        return false;
    }
    const int n_chroms = hdr_tmp->n_targets;
    for (const auto& kv : het_snpdat){
        if (kv.first < 0 || kv.first >= n_chroms){
            fprintf(stderr, "ERROR: het panel references invalid BAM target id %d\n", kv.first);
            hts_idx_destroy(idx_tmp);
            bam_hdr_destroy(hdr_tmp);
            hts_close(bam_tmp);
            return false;
        }
    }
    hts_idx_destroy(idx_tmp);
    bam_hdr_destroy(hdr_tmp);
    hts_close(bam_tmp);

    vector<int> chrom_work;
    for (int tid = 0; tid < n_chroms; ++tid){
        auto it = het_snpdat.find(tid);
        if (it != het_snpdat.end() && !it->second.empty()) chrom_work.push_back(tid);
    }
    fprintf(stderr, "Processing %lu chromosomes with %d threads...\n",
        chrom_work.size(), n_threads);

    atomic<int> chroms_done(0);
    atomic<long> records_processed(0);
    atomic<long> sites_hit(0);
    ParallelOperationStatus operation_status;
    std::atomic<bool> hts_thread_warning_emitted(false);
    vector<robin_hood::unordered_map<unsigned long, CellHetData> > thread_het_data(n_threads);

    #pragma omp parallel num_threads(n_threads)
    {
        const int thread_id = omp_get_thread_num();
        auto& local_het_data = thread_het_data[thread_id];
        if (has_bc_list){
            try {
                for (unsigned long bc : valid_barcodes){
                    local_het_data.emplace(piecewise_construct,
                        forward_as_tuple(bc), forward_as_tuple(n_samples));
                }
            }
            catch (const std::exception& e){
                operation_status.fail(std::string("het worker allocation failed: ") + e.what());
            }
        }

        htsFile* bam_fp = hts_open(bamfile.c_str(), "r");
        bam_hdr_t* header = nullptr;
        hts_idx_t* idx = nullptr;
        bam1_t* record = nullptr;
        if (!bam_fp){
            operation_status.fail(format_worker_error("BAM open", thread_id));
        }
        else{
            if (htslib_threads > 1 && hts_set_threads(bam_fp, htslib_threads) < 0){
                bool expected = false;
                if (hts_thread_warning_emitted.compare_exchange_strong(expected, true)){
                    fprintf(stderr,
                        "WARNING: HTSlib helper-thread setup failed; continuing with synchronous BAM I/O\n");
                }
            }
            header = sam_hdr_read(bam_fp);
            if (!header) operation_status.fail(format_worker_error("BAM header read", thread_id));
            idx = sam_index_load(bam_fp, bamfile.c_str());
            if (!idx) operation_status.fail(format_worker_error("BAM index load", thread_id));
            record = bam_init1();
            if (!record) operation_status.fail(format_worker_error("BAM record allocation", thread_id));
        }

        #pragma omp for schedule(dynamic, 1)
        for (size_t wi = 0; wi < chrom_work.size(); ++wi){
            if (!operation_status.ok() || !bam_fp || !header || !idx || !record) continue;
            const int tid = chrom_work[wi];
            if (tid < 0 || tid >= header->n_targets){
                operation_status.fail(format_worker_error("invalid contig", thread_id, tid));
                continue;
            }
            auto snp_map_it = het_snpdat.find(tid);
            if (snp_map_it == het_snpdat.end()){
                operation_status.fail(format_worker_error("missing het contig data", thread_id, tid));
                continue;
            }
            ChromSNPs& chrom_snps = snp_map_it->second;
            hts_itr_t* iter = sam_itr_queryi(idx, tid, 0, INT_MAX);
            if (!iter){
                operation_status.fail(format_worker_error("iterator creation", thread_id, tid));
                continue;
            }

            long local_records = 0;
            long local_sites = 0;
            auto snp_iter = chrom_snps.snps.begin();
            auto snp_end = chrom_snps.snps.end();
            int iterator_result = 0;
            while ((iterator_result = sam_itr_next(bam_fp, iter, record)) >= 0){
                ++local_records;
                if (!read_passes_filter(record, default_production_read_filter())) continue;

                uint8_t* cb_ptr = bam_aux_get(record, "CB");
                if (!cb_ptr) continue;
                const char* cb_str = bam_aux2Z(cb_ptr);
                if (!cb_str) continue;
                const unsigned long bc_key = bc_ul((char*)cb_str);
                if (has_bc_list && valid_barcodes.find(bc_key) == valid_barcodes.end()) continue;

                const int read_start = record->core.pos;
                const int read_end = bam_endpos(record);
                while (snp_iter != snp_end && snp_iter->pos < read_start) ++snp_iter;

                const float prob_correct =
                    1.0f - powf(10.0f, -(float)record->core.qual / 10.0f);
                for (auto snp_check = snp_iter;
                     snp_check != snp_end && snp_check->pos < read_end;
                     ++snp_check){
                    const char allele = get_base_at_pos(record, snp_check->pos);
                    if (allele == 'N' || allele == '-') continue;
                    float ref_add = 0.0f;
                    float alt_add = 0.0f;
                    if (allele == snp_check->data.ref) ref_add = prob_correct;
                    else if (allele == snp_check->data.alt) alt_add = prob_correct;
                    else continue;
                    ++local_sites;

                    auto cell_it = local_het_data.find(bc_key);
                    if (cell_it == local_het_data.end()){
                        if (has_bc_list) continue;
                        try {
                            local_het_data.emplace(bc_key, CellHetData(n_samples));
                        }
                        catch (const std::exception& e){
                            operation_status.fail(std::string("het worker allocation failed: ") + e.what());
                            break;
                        }
                        cell_it = local_het_data.find(bc_key);
                    }
                    CellHetData& cell_data = cell_it->second;

                    if (collect_persite){
                        const auto site_it = site_to_idx.find(make_pair(tid, snp_check->pos));
                        if (site_it == site_to_idx.end()){
                            operation_status.fail(format_worker_error("het site lookup", thread_id, tid));
                            break;
                        }
                        cell_data.persite_data.add_site(site_it->second, ref_add, alt_add);
                    }
                    if (collect_welford){
                        const float depth = ref_add + alt_add;
                        if (depth > 0.0f){
                            const float alt_frac = alt_add / depth;
                            for (int indiv = 0; indiv < n_samples; ++indiv){
                                if (snp_check->data.is_het(indiv)){
                                    cell_data.welford_stats.add(indiv, alt_frac, depth);
                                }
                            }
                        }
                    }
                }
                if (!operation_status.ok()) break;
            }
            if (iterator_result < -1){
                operation_status.fail(format_worker_error("iterator read", thread_id, tid));
            }
            hts_itr_destroy(iter);
            records_processed += local_records;
            sites_hit += local_sites;
            const int done = ++chroms_done;
            if (done % 5 == 0 || done == (int)chrom_work.size()){
                fprintf(stderr,
                    "\rHet counting: %d/%lu chroms, %ld iterator records, %ld site hits    ",
                    done, chrom_work.size(), records_processed.load(), sites_hit.load());
            }
        }

        if (record) bam_destroy1(record);
        if (idx) hts_idx_destroy(idx);
        if (header) bam_hdr_destroy(header);
        if (bam_fp) hts_close(bam_fp);
    }

    if (!operation_status.ok()){
        fprintf(stderr, "ERROR: het/ploidy counting failed: %s\n",
            operation_status.message().c_str());
        return false;
    }

    for (int t = 0; t < n_threads; ++t){
        for (auto& kv : thread_het_data[t]){
            auto it = het_data.find(kv.first);
            if (it != het_data.end()) it->second.merge(kv.second);
            else het_data.emplace(kv.first, std::move(kv.second));
        }
        thread_het_data[t].clear();
    }

    fprintf(stderr,
        "\nHet counting complete: %lu cells, %ld iterator records, %ld site hits\n",
        het_data.size(), records_processed.load(), sites_hit.load());
    if (collect_persite){
        long total_sites_stored = 0;
        long max_sites = 0;
        for (const auto& kv : het_data){
            const long n = (long)kv.second.persite_data.size();
            total_sites_stored += n;
            if (n > max_sites) max_sites = n;
        }
        fprintf(stderr, "Per-site: %ld total entries, max %ld/cell, avg %.1f/cell\n",
            total_sites_stored, max_sites,
            het_data.empty() ? 0.0 : (double)total_sites_stored / het_data.size());
    }
    if (collect_welford){
        long max_sites = 0;
        for (const auto& kv : het_data){
            for (int i = 0; i < n_samples; ++i){
                const long n = (long)kv.second.welford_stats.get(i).n;
                if (n > max_sites) max_sites = n;
            }
        }
        fprintf(stderr, "Max het sites per individual per cell: %ld\n", max_sites);
    }
    return true;
}


// Serial compatibility helper retained in the unified HTS module.
void process_bam_record_bulk(bam_reader& reader,
    int snppos,
    var& vardat,
    map<int, pair<float, float> >& snp_ref_alt,
    map<int, float>& snp_err,
    bool genes,
    map<pair<int, int>, set<string> >& gene_ids,
    map<string, string>& gene_id2name){

    if (!reader.unmapped() && !reader.secondary() && 
        !reader.dup() && reader.has_cb_z){
         
        int tid = reader.tid();

        // Instead of storing actual read counts, store the probability
        // that the mapping was correct.
        float prob_corr = 1.0 - pow(10, -(float)reader.mapq/10.0);
        
        if (snp_ref_alt.count(snppos) == 0){
            snp_ref_alt.insert(make_pair(snppos, make_pair(0.0, 0.0)));
            snp_err.insert(make_pair(snppos, 0.0));
        }

        // Note: this function expects positions to be 1-based, but 
        // BCF/BAM functions store as 0-based
        char allele = reader.get_base_at(snppos + 1);
        
        if (allele != 'N' && allele != '-'){
            if (allele == vardat.ref){
                snp_ref_alt[snppos].first += prob_corr;
            }
            else if (allele == vardat.alt){
                snp_ref_alt[snppos].second += prob_corr;
            }
            else{
                snp_err[snppos] += prob_corr;
            }
        }

        if (genes){
            pair<int, int> key = make_pair(tid, snppos);
            if (gene_ids.count(key) == 0){
                set<string> s;
                gene_ids.insert(make_pair(key, s));
            }
            if (reader.gene_names.size() > 0 && reader.gene_ids.size() < reader.gene_names.size()){
                // Use names as IDs
                for (int i = 0; i < reader.gene_names.size(); ++i){
                    gene_ids[key].insert(reader.gene_names[i]);
                }
            }
            else if (reader.gene_ids.size() > 0 && reader.gene_names.size() < reader.gene_ids.size()){
                // Ignore gene names (can't map from IDs -> names)
                for (int i = 0; i < reader.gene_ids.size(); ++i){
                    gene_ids[key].insert(reader.gene_ids[i]);
                }
            }
            else if (reader.gene_ids.size() == reader.gene_names.size()){
                // Store gene IDs and map to names
                for (int i = 0; i < reader.gene_ids.size(); ++i){
                    gene_ids[key].insert(reader.gene_ids[i]);
                    if (gene_id2name.count(reader.gene_ids[i]) == 0){
                        gene_id2name.insert(make_pair(reader.gene_ids[i], reader.gene_names[i]));
                    }
                }
            }
        }
    }
}


// Serial compatibility helper retained in the unified HTS module.
void process_bam_record_bysnp(bam_reader& reader,
    int snppos,
    var& vardat,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    map<int, pair<float, float> >& snp_var_counts){

    if (!reader.unmapped() && !reader.secondary() && 
        !reader.dup() && reader.has_cb_z){
                        
        // Get BC key
        bc bc_bits;
        str2bc(reader.cb_z, bc_bits);
        unsigned long bc_key = bc_bits.to_ulong();
        
        if (assignments.count(bc_key) > 0){
            
            // Instead of storing actual read counts, store the probability
            // that the mapping was correct.
            float prob_corr = 1.0 - pow(10, -(float)reader.mapq/10.0);
            
            if (prob_corr > 0.001){
                int a = assignments[bc_key];
                if (snp_var_counts.count(a) == 0){
                    snp_var_counts.insert(make_pair(a, make_pair(0,0)));
                }
                
                // Note: this function expects positions to be 1-based, but 
                // BCF/BAM functions store as 0-based
                char allele = reader.get_base_at(snppos + 1);
                
                if (allele != 'N' && allele != '-'){
                    if (allele == vardat.ref){
                        snp_var_counts[a].first += prob_corr;
                    }
                    else if (allele == vardat.alt){
                        snp_var_counts[a].second += prob_corr;
                    }
                }
            }
        }
    }
}
