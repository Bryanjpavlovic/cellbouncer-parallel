#ifndef _CELLBOUNCER_DEMUX_PARALLEL_HTS_H
#define _CELLBOUNCER_DEMUX_PARALLEL_HTS_H
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
#include <zlib.h>
#include <mutex>
#include <atomic>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htswrapper/bam.h>
#include <htswrapper/bc.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"

using std::cout;
using std::endl;
using namespace std;

// ============================================================================
// CONSTANTS FOR OPTIMIZED DATA STRUCTURES
// ============================================================================

// Maximum number of individuals supported (increase if needed)
constexpr int MAX_INDIVIDUALS = 64;
// Number of genotype states (0, 1, 2 alt alleles)
constexpr int GENOTYPE_STATES = 3;
// Total number of (individual, genotype) state combinations
constexpr int STATE_COUNT = MAX_INDIVIDUALS * GENOTYPE_STATES;

// ============================================================================
// ORIGINAL VAR STRUCTURE (unchanged for compatibility)
// ============================================================================

/**
 * Structure to represent population-level genotype information at a 
 * single SNP.
 */
struct var{
    char ref;
    char alt;
    
    // This limits us to 500 individuals in the input VCF
    bitset<500> haps1;
    bitset<500> haps2;
    bitset<500> haps_covered;
    vector<float> gqs;
    float vq;
    var(){
        this->ref = 'N';
        this->alt = 'N';
        this->haps1.reset();
        this->haps2.reset();
        this->haps_covered.reset();
        this->vq = 0.0;
    }
    var(const var& v){
        this->ref = v.ref;
        this->alt = v.alt;
        this->haps1 = v.haps1;
        this->haps2 = v.haps2;
        this->haps_covered = v.haps_covered;
        this->gqs = v.gqs;
        this->vq = v.vq;
    }
};

// ============================================================================
// NEW OPTIMIZED DATA STRUCTURES
// ============================================================================

/**
 * SNPData - combines position and variant data for cache-efficient iteration
 */
struct SNPData {
    int pos;
    var data;
    
    SNPData() : pos(0) {}
    SNPData(int p, const var& v) : pos(p), data(v) {}
    
    bool operator<(const SNPData& other) const {
        return pos < other.pos;
    }
};

/**
 * ChromSNPs - per-chromosome SNP storage with sorted vector for fast lookup
 */
struct ChromSNPs {
    std::vector<SNPData> snps;  // Sorted by position
    
    // Binary search for first SNP >= pos
    std::vector<SNPData>::iterator lower_bound(int pos) {
        SNPData target;
        target.pos = pos;
        return std::lower_bound(snps.begin(), snps.end(), target);
    }
    
    // Binary search for first SNP > pos  
    std::vector<SNPData>::iterator upper_bound(int pos) {
        SNPData target;
        target.pos = pos;
        return std::upper_bound(snps.begin(), snps.end(), target);
    }
    
    // Const versions
    std::vector<SNPData>::const_iterator lower_bound(int pos) const {
        SNPData target;
        target.pos = pos;
        return std::lower_bound(snps.begin(), snps.end(), target);
    }
    
    std::vector<SNPData>::const_iterator upper_bound(int pos) const {
        SNPData target;
        target.pos = pos;
        return std::upper_bound(snps.begin(), snps.end(), target);
    }
    
    void sort_snps() {
        std::sort(snps.begin(), snps.end());
    }
    
    bool empty() const { return snps.empty(); }
    size_t size() const { return snps.size(); }
};

/**
 * CellCounts - Flat, cache-friendly structure for allele counts per cell
 * 
 * For n individuals with 3 genotype states each, we store counts in a 
 * flattened 2D array. Index [indv1 * 3 + nalt1][indv2 * 3 + nalt2] gives
 * counts for SNPs where individual 1 has nalt1 alt alleles and individual 2
 * has nalt2 alt alleles.
 * 
 * We use a dynamic size based on actual number of samples to avoid
 * allocating huge arrays.
 */
struct CellCounts {
    int n_samples;
    int state_count;
    
    // Flattened storage for pairwise counts: ref_counts[idx1 * state_count + idx2]
    std::vector<float> ref_counts;
    std::vector<float> alt_counts;
    
    // SEPARATE storage for totals per individual/genotype (fixes index collision bug)
    std::vector<float> total_ref;
    std::vector<float> total_alt;
    
    CellCounts() : n_samples(0), state_count(0) {}
    
    CellCounts(int n_samp) : n_samples(n_samp), state_count(n_samp * GENOTYPE_STATES) {
        size_t total_size = (size_t)state_count * state_count;
        ref_counts.resize(total_size, 0.0f);
        alt_counts.resize(total_size, 0.0f);
        total_ref.resize(state_count, 0.0f);
        total_alt.resize(state_count, 0.0f);
    }
    
    CellCounts(const CellCounts& other) = default;
    CellCounts& operator=(const CellCounts& other) = default;
    CellCounts(CellCounts&& other) = default;
    CellCounts& operator=(CellCounts&& other) = default;
    
    inline void add(int indv1, int nalt1, int indv2, int nalt2, float ref, float alt) {
        int idx1 = indv1 * GENOTYPE_STATES + nalt1;
        int idx2 = indv2 * GENOTYPE_STATES + nalt2;
        size_t flat_idx = (size_t)idx1 * state_count + idx2;
        ref_counts[flat_idx] += ref;
        alt_counts[flat_idx] += alt;
    }
    
    inline std::pair<float, float> get(int indv1, int nalt1, int indv2, int nalt2) const {
        int idx1 = indv1 * GENOTYPE_STATES + nalt1;
        int idx2 = indv2 * GENOTYPE_STATES + nalt2;
        size_t flat_idx = (size_t)idx1 * state_count + idx2;
        return {ref_counts[flat_idx], alt_counts[flat_idx]};
    }
    
    // Totals stored in separate vectors - no collision with pairwise data
    inline void add_total(int indv, int nalt, float ref, float alt) {
        int idx = indv * GENOTYPE_STATES + nalt;
        total_ref[idx] += ref;
        total_alt[idx] += alt;
    }
    
    inline std::pair<float, float> get_total(int indv, int nalt) const {
        int idx = indv * GENOTYPE_STATES + nalt;
        return {total_ref[idx], total_alt[idx]};
    }
    
    // Merge another CellCounts into this one (for deterministic accumulation)
    void merge(const CellCounts& other) {
        if (other.n_samples == 0) return;
        
        // Initialize if we're empty
        if (n_samples == 0) {
            n_samples = other.n_samples;
            state_count = other.state_count;
            ref_counts = other.ref_counts;
            alt_counts = other.alt_counts;
            total_ref = other.total_ref;
            total_alt = other.total_alt;
            return;
        }
        
        // Add other's counts to ours
        for (size_t i = 0; i < ref_counts.size(); ++i) {
            ref_counts[i] += other.ref_counts[i];
            alt_counts[i] += other.alt_counts[i];
        }
        for (size_t i = 0; i < total_ref.size(); ++i) {
            total_ref[i] += other.total_ref[i];
            total_alt[i] += other.total_alt[i];
        }
    }
    
    void clear() {
        std::fill(ref_counts.begin(), ref_counts.end(), 0.0f);
        std::fill(alt_counts.begin(), alt_counts.end(), 0.0f);
        std::fill(total_ref.begin(), total_ref.end(), 0.0f);
        std::fill(total_alt.begin(), total_alt.end(), 0.0f);
    }
    
    bool is_empty() const {
        for (size_t i = 0; i < total_ref.size(); ++i) {
            if (total_ref[i] > 0 || total_alt[i] > 0) return false;
        }
        return true;
    }
};

/**
 * AlignedCellCounts - Cache-line aligned wrapper with per-cell mutex
 * The alignas(64) ensures each cell's data starts on a cache line boundary,
 * preventing false sharing between threads.
 */
struct alignas(64) AlignedCellCounts {
    CellCounts counts;
    std::mutex lock;
    
    AlignedCellCounts() = default;
    AlignedCellCounts(int n_samples) : counts(n_samples) {}
    
    // Disable copy (mutex is not copyable)
    AlignedCellCounts(const AlignedCellCounts&) = delete;
    AlignedCellCounts& operator=(const AlignedCellCounts&) = delete;
    
    // Allow move
    AlignedCellCounts(AlignedCellCounts&& other) noexcept 
        : counts(std::move(other.counts)) {}
    AlignedCellCounts& operator=(AlignedCellCounts&& other) noexcept {
        counts = std::move(other.counts);
        return *this;
    }
};

// ============================================================================
// SHARED MEMORY STRUCTURES FOR VCF SHARING
// ============================================================================

/**
 * Header structure for shared memory VCF data
 */
struct SharedVCFHeader {
    size_t total_size;
    int n_chromosomes;
    int n_samples;
    int n_snps_total;
    // Offsets to chromosome data (max 8192 chromosomes)
    size_t chrom_offsets[8192];
    size_t chrom_snp_counts[8192];
    int chrom_tids[8192];
};

// ============================================================================
// VCF READING FUNCTIONS (original interface preserved)
// ============================================================================

void read_vcf_samples(std::string& filename,
    std::vector<std::string>& samples);

int read_vcf_chrom(std::string& vcf_file,
    std::string& chrom,
    std::map<int, var>& snps,
    int min_vq,
    bool allow_missing=true);

void get_vcf_chroms(std::string& vcf_file,
    std::set<std::string>& chroms);

void get_bam_chroms(bam_reader& reader,
    std::set<std::string>& chroms);

long int count_vcf_snps(std::string& vcf_file,
    std::set<std::string>& chroms_to_include,
    int min_vq);

int read_vcf_chroms(std::string& vcf_file,
    std::set<std::string>& chroms_to_include,
    std::map<std::string, int>& seq2tid,
    std::map<int, std::map<int, var> >& snps,
    int min_vq,
    bool allow_missing=true);

int read_vcf(std::string& filename, 
    bam_reader& reader,
    std::vector<std::string>& samples,
    std::map<int, std::map<int, var> >& snps,
    int min_vq,
    bool hdr_only,
    bool skip_seq2tid,
    bool allow_missing=true);

// ============================================================================
// NEW VCF READING FUNCTIONS (optimized structures)
// ============================================================================

/**
 * Read VCF data into optimized ChromSNPs structure
 */
int read_vcf_chroms_optimized(std::string& vcf_file,
    std::set<std::string>& chroms_to_include,
    std::map<std::string, int>& seq2tid,
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_optimized,
    int min_vq,
    bool allow_missing=true);

/**
 * Convert from old map structure to new optimized structure
 */
void convert_snpdat_to_optimized(
    std::map<int, std::map<int, var> >& snpdat_old,
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_new);

// ============================================================================
// CONDITIONAL MATCH FRACTION FUNCTIONS
// ============================================================================

void get_conditional_match_fracs(std::map<int, std::map<int, var> >& snpdat,
    std::map<std::pair<int, int>, std::map<int, float> >& conditional_match_fracs, 
    int n_samples);

void get_conditional_match_fracs_chrom(std::map<int, var>& snpdat,
    std::map<std::pair<int, int>, std::map<int, float> >& conditional_match_fracs,
    std::map<std::pair<int, int>, std::map<int, float> >& conditional_match_tots,
    int n_samples);

void get_conditional_match_fracs_chrom_optimized(ChromSNPs& snpdat,
    std::map<std::pair<int, int>, std::map<int, float> >& conditional_match_fracs,
    std::map<std::pair<int, int>, std::map<int, float> >& conditional_match_tots,
    int n_samples);

void conditional_match_fracs_normalize(std::map<std::pair<int, int>, std::map<int, float> >& conditional_match_fracs,
    std::map<std::pair<int, int>, std::map<int, float> >& conditional_match_tots,
    int n_samples);

// ============================================================================
// BAM PROCESSING FUNCTIONS (original interface)
// ============================================================================

void process_bam_record(bam_reader& reader,
    int snppos,
    var& vardat,
    std::map<int, robin_hood::unordered_map<unsigned long, 
        std::pair<float, float> > >& var_counts,
    bool has_bc_list,
    std::set<unsigned long>& bcs_valid);

void process_bam_record_bulk(bam_reader& reader,
    int snppos,
    var& vardat,
    std::map<int, std::pair<float, float> >& snp_ref_alt,
    std::map<int, float>& snp_err,
    bool genes,
    std::map<std::pair<int, int>, std::set<std::string> >& snp_gene_ids,
    std::map<std::string, std::string>& gene_id2name);

void process_bam_record_bysnp(bam_reader& reader,
    int snppos,
    var& vardat,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    std::map<int, std::pair<float, float> >& snp_var_counts);

void dump_vcs_counts(robin_hood::unordered_map<unsigned long, 
        std::pair<float, float> >& varcounts_site,
    robin_hood::unordered_map<unsigned long, std::map<std::pair<int, int>, 
        std::map<std::pair<int, int>, std::pair<float, float> > > >& indv_allelecounts,
    var& snpdat,
    int n_samples);

// ============================================================================
// NEW PARALLEL BAM PROCESSING FUNCTIONS
// ============================================================================

/**
 * Get base at a specific position from a BAM record, accounting for CIGAR
 */
char get_base_at_pos(bam1_t* record, int pos);

/**
 * Main parallel counting function - processes BAM file using multiple threads
 * 
 * @param bamfile Path to BAM file
 * @param snpdat_all Optimized SNP data structure (tid -> ChromSNPs)
 * @param cell_counts Output: per-cell allele counts
 * @param valid_barcodes Set of barcodes to process (empty = all)
 * @param n_samples Number of individuals in VCF
 * @param n_threads Number of OpenMP threads to use
 * @param htslib_threads Number of decompression threads per BAM reader
 */
void count_alleles_parallel(
    const std::string& bamfile,
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_all,
    robin_hood::unordered_map<unsigned long, AlignedCellCounts>& cell_counts,
    const std::set<unsigned long>& valid_barcodes,
    int n_samples,
    int n_threads,
    int htslib_threads);

/**
 * Single-threaded counting using optimized data structures
 * Fallback for --threads 1 or when parallel processing is disabled
 */
void count_alleles_single_threaded(
    const std::string& bamfile,
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_all,
    robin_hood::unordered_map<unsigned long, CellCounts>& cell_counts,
    const std::set<unsigned long>& valid_barcodes,
    int n_samples,
    std::map<std::pair<int, int>, std::map<int, float> >& conditional_match_fracs,
    std::map<std::pair<int, int>, std::map<int, float> >& conditional_match_tots,
    bool compute_conditional);

/**
 * Convert parallel cell counts to non-mutex version for LLR computation
 */
void finalize_parallel_counts(
    robin_hood::unordered_map<unsigned long, AlignedCellCounts>& parallel_counts,
    robin_hood::unordered_map<unsigned long, CellCounts>& final_counts);

// ============================================================================
// SHARED MEMORY VCF FUNCTIONS
// ============================================================================

/**
 * Create shared memory segment with VCF data
 */
bool create_shared_vcf(
    const std::string& vcf_file,
    const std::string& shm_name,
    std::set<std::string>& chroms_to_include,
    std::map<std::string, int>& seq2tid,
    int min_vq);

/**
 * Attach to existing shared memory VCF
 */
bool attach_shared_vcf(
    const std::string& shm_name,
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_all,
    std::vector<std::string>& samples);

/**
 * Detach from shared memory VCF (does not destroy)
 */
void detach_shared_vcf(const std::string& shm_name);

/**
 * Destroy shared memory segment
 */
void destroy_shared_vcf(const std::string& shm_name);

#endif 
