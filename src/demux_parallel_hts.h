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
#include <cstdint>
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
// HET BALANCE COMPUTATION OPTIONS
// ============================================================================

/**
 * Het balance computation method:
 * - WELFORD: Online Welford variance per individual during counting (default, low memory)
 * - PERSITE: Store per-site (ref, alt) counts, compute variance after assignment (more memory)
 */
enum class HetBalanceMethod {
    WELFORD,      // Online Welford variance (default)
    PERSITE       // Store per-site counts
};

// ============================================================================
// PER-SITE HET DATA STORAGE
// ============================================================================

/**
 * Per-site het count data - stores (ref, alt) for each site with coverage.
 */
struct HetSiteCount {
    int32_t site_idx;      // Global index into het VCF sites
    float ref;             // Ref allele count (quality-weighted)
    float alt;             // Alt allele count (quality-weighted)
    
    HetSiteCount() : site_idx(-1), ref(0), alt(0) {}
    HetSiteCount(int32_t idx, float r, float a) : site_idx(idx), ref(r), alt(a) {}
};

/**
 * Per-cell het site storage for PERSITE method.
 */
struct HetSiteData {
    std::vector<HetSiteCount> sites;
    
    void add_site(int32_t site_idx, float ref, float alt) {
        for (auto& s : sites) {
            if (s.site_idx == site_idx) {
                s.ref += ref;
                s.alt += alt;
                return;
            }
        }
        sites.emplace_back(site_idx, ref, alt);
    }
    
    void clear() { sites.clear(); }
    size_t size() const { return sites.size(); }
    bool empty() const { return sites.empty(); }
    
    void merge(const HetSiteData& other) {
        for (const auto& s : other.sites) {
            add_site(s.site_idx, s.ref, s.alt);
        }
    }
};

// ============================================================================
// WELFORD ONLINE VARIANCE STATISTICS
// ============================================================================

/**
 * Welford online variance for a single individual.
 * Tracks running mean and variance without storing per-site data.
 */
struct WelfordStats {
    int64_t n;           // Number of sites
    double mean;         // Running mean of alt_frac
    double M2;           // Sum of squared deviations
    double total_depth;  // Total depth at het sites
    
    WelfordStats() : n(0), mean(0.0), M2(0.0), total_depth(0.0) {}
    
    void add(double alt_frac, double depth) {
        n++;
        total_depth += depth;
        double delta = alt_frac - mean;
        mean += delta / n;
        double delta2 = alt_frac - mean;
        M2 += delta * delta2;
    }
    
    double variance(int min_sites = 10) const {
        if (n < min_sites) return -1.0;
        return M2 / (n - 1);
    }
    
    // Chan's parallel merge algorithm
    void merge(const WelfordStats& other) {
        if (other.n == 0) return;
        if (n == 0) { *this = other; return; }
        
        int64_t n_combined = n + other.n;
        double delta = other.mean - mean;
        double mean_combined = mean + delta * other.n / n_combined;
        double M2_combined = M2 + other.M2 + delta * delta * n * other.n / n_combined;
        
        n = n_combined;
        mean = mean_combined;
        M2 = M2_combined;
        total_depth += other.total_depth;
    }
    
    void clear() { n = 0; mean = 0.0; M2 = 0.0; total_depth = 0.0; }
};

/**
 * Per-cell Welford stats for all individuals.
 */
struct CellWelfordStats {
    int n_samples;
    std::vector<WelfordStats> indiv_stats;
    
    CellWelfordStats() : n_samples(0) {}
    CellWelfordStats(int n_samp) : n_samples(n_samp), indiv_stats(n_samp) {}
    
    void add(int indiv, double alt_frac, double depth) {
        if (indiv >= 0 && indiv < n_samples) {
            indiv_stats[indiv].add(alt_frac, depth);
        }
    }
    
    const WelfordStats& get(int indiv) const {
        static WelfordStats empty;
        if (indiv >= 0 && indiv < n_samples) return indiv_stats[indiv];
        return empty;
    }
    
    void merge(const CellWelfordStats& other) {
        if (other.n_samples == 0) return;
        if (n_samples == 0) { *this = other; return; }
        for (int i = 0; i < n_samples; ++i) {
            indiv_stats[i].merge(other.indiv_stats[i]);
        }
    }
    
    void clear() { for (auto& s : indiv_stats) s.clear(); }
};

/**
 * Per-cell het data (supports both methods).
 */
struct CellHetData {
    HetSiteData persite_data;       // For PERSITE method
    CellWelfordStats welford_stats; // For WELFORD method
    
    CellHetData() {}
    CellHetData(int n_samples) : welford_stats(n_samples) {}
    
    void merge(const CellHetData& other) {
        persite_data.merge(other.persite_data);
        welford_stats.merge(other.welford_stats);
    }
    
    void clear() {
        persite_data.clear();
        welford_stats.clear();
    }
};

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
    
    // Check if individual is het at this site
    bool is_het(int indiv) const {
        if (!haps_covered.test(indiv)) return false;
        return haps1.test(indiv) != haps2.test(indiv);
    }
    
    // Get genotype (0, 1, 2) for individual, -1 if not covered
    int get_genotype(int indiv) const {
        if (!haps_covered.test(indiv)) return -1;
        int nalt = 0;
        if (haps1.test(indiv)) nalt++;
        if (haps2.test(indiv)) nalt++;
        return nalt;
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
 * 
 * FIXED-POINT ARITHMETIC: Values are stored as int64_t scaled by FIXED_POINT_SCALE.
 * This makes accumulation order-independent (integer addition is associative),
 * enabling deterministic results regardless of thread scheduling.
 */

// Scale factor for fixed-point arithmetic
// 1,000,000 gives 6 decimal places of precision, well within int64 range
constexpr int64_t FIXED_POINT_SCALE = 1000000;

struct CellCounts {
    int n_samples;
    int state_count;
    
    // Flattened storage for pairwise counts (fixed-point int64)
    std::vector<int64_t> ref_counts;
    std::vector<int64_t> alt_counts;
    
    // SEPARATE storage for totals per individual/genotype (fixed-point int64)
    std::vector<int64_t> total_ref;
    std::vector<int64_t> total_alt;
    
    CellCounts() : n_samples(0), state_count(0) {}
    
    CellCounts(int n_samp) : n_samples(n_samp), state_count(n_samp * GENOTYPE_STATES) {
        size_t total_size = (size_t)state_count * state_count;
        ref_counts.resize(total_size, 0);
        alt_counts.resize(total_size, 0);
        total_ref.resize(state_count, 0);
        total_alt.resize(state_count, 0);
    }
    
    CellCounts(const CellCounts& other) = default;
    CellCounts& operator=(const CellCounts& other) = default;
    CellCounts(CellCounts&& other) = default;
    CellCounts& operator=(CellCounts&& other) = default;
    
    // Add using pre-scaled int64 values (caller scales by FIXED_POINT_SCALE)
    inline void add(int indv1, int nalt1, int indv2, int nalt2, int64_t ref, int64_t alt) {
        int idx1 = indv1 * GENOTYPE_STATES + nalt1;
        int idx2 = indv2 * GENOTYPE_STATES + nalt2;
        size_t flat_idx = (size_t)idx1 * state_count + idx2;
        ref_counts[flat_idx] += ref;
        alt_counts[flat_idx] += alt;
    }
    
    // Get returns float (converts from fixed-point)
    inline std::pair<float, float> get(int indv1, int nalt1, int indv2, int nalt2) const {
        int idx1 = indv1 * GENOTYPE_STATES + nalt1;
        int idx2 = indv2 * GENOTYPE_STATES + nalt2;
        size_t flat_idx = (size_t)idx1 * state_count + idx2;
        return {(float)ref_counts[flat_idx] / FIXED_POINT_SCALE, 
                (float)alt_counts[flat_idx] / FIXED_POINT_SCALE};
    }
    
    // Add totals using pre-scaled int64 values
    inline void add_total(int indv, int nalt, int64_t ref, int64_t alt) {
        int idx = indv * GENOTYPE_STATES + nalt;
        total_ref[idx] += ref;
        total_alt[idx] += alt;
    }
    
    // Get totals returns float (converts from fixed-point)
    inline std::pair<float, float> get_total(int indv, int nalt) const {
        int idx = indv * GENOTYPE_STATES + nalt;
        return {(float)total_ref[idx] / FIXED_POINT_SCALE, 
                (float)total_alt[idx] / FIXED_POINT_SCALE};
    }
    
    // Merge another CellCounts into this one
    // Integer addition is associative - order doesn't matter!
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
        
        // Add other's counts to ours (exact integer arithmetic)
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
        std::fill(ref_counts.begin(), ref_counts.end(), 0);
        std::fill(alt_counts.begin(), alt_counts.end(), 0);
        std::fill(total_ref.begin(), total_ref.end(), 0);
        std::fill(total_alt.begin(), total_alt.end(), 0);
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
    // Sample names (max 512 samples, 64 chars each)
    char sample_names[512][64];
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

// ============================================================================
// HET VCF FUNCTIONS (NEW - for ploidy detection)
// ============================================================================

/**
 * Load het VCF into ChromSNPs structure for parallel processing
 * 
 * @param het_vcf_file Path to het-enriched VCF from downsample_vcf_parallel
 * @param chroms_to_include Set of chromosomes to load
 * @param seq2tid Chromosome name to tid mapping
 * @param het_snpdat Output: SNP data structure
 * @param min_vq Minimum variant quality
 * @return Number of variants loaded
 */
int load_het_vcf(
    const std::string& het_vcf_file,
    const std::set<std::string>& chroms_to_include,
    std::map<std::string, int>& seq2tid,
    robin_hood::unordered_map<int, ChromSNPs>& het_snpdat,
    int min_vq);

/**
 * Count alleles at het sites in parallel
 * Similar to count_alleles_parallel but for het VCF
 * 
 * @param bamfile Path to BAM file
 * @param het_snpdat Het SNP data (tid -> ChromSNPs)
 * @param het_counts Output: per-cell allele counts at het sites
 * @param valid_barcodes Set of barcodes to process
 * @param n_samples Number of individuals in VCF
 * @param n_threads Number of OpenMP threads
 * @param htslib_threads Number of decompression threads per reader
 */
void count_het_alleles_parallel(
    const std::string& bamfile,
    robin_hood::unordered_map<int, ChromSNPs>& het_snpdat,
    robin_hood::unordered_map<unsigned long, CellCounts>& het_counts,
    const std::set<unsigned long>& valid_barcodes,
    int n_samples,
    int n_threads,
    int htslib_threads);

/**
 * Extended het allele counting with Welford online variance or per-site storage.
 * 
 * @param het_data Output: per-cell het data
 * @param idx_to_site Output: site index to (tid, pos) mapping (only for PERSITE)
 * @param method Which method to use (WELFORD or PERSITE)
 */
void count_het_alleles_extended(
    const std::string& bamfile,
    robin_hood::unordered_map<int, ChromSNPs>& het_snpdat,
    robin_hood::unordered_map<unsigned long, CellHetData>& het_data,
    std::vector<std::pair<int, int>>& idx_to_site,
    const std::set<unsigned long>& valid_barcodes,
    int n_samples,
    int n_threads,
    int htslib_threads,
    HetBalanceMethod method);

#endif 
