// =============================================================================
// vcf_hts.cpp
// Unified VCF/BCF, shared-memory panel, BAM counting, and HTS helpers.
// =============================================================================

#include <string>
#include <climits>
#include <algorithm>
#include <array>
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
#include <memory>
#include <chrono>
#include <cmath>
#include <exception>
#include <limits>
#include <queue>
#include <functional>
#include <omp.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htswrapper/bc.h>
#include <htswrapper/bam.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "vcf_hts.h"

#ifndef CELLBOUNCER_VCF_HTS_INTERFACE_REVISION
#error "vcf_hts.cpp requires the matching src/vcf_hts.h (interface revision 21901)"
#elif CELLBOUNCER_VCF_HTS_INTERFACE_REVISION != 21901
#error "vcf_hts.cpp and src/vcf_hts.h are from different CellBouncer source revisions"
#else

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

static const std::array<int64_t, 256>& mapq_probability_scaled_table() {
    static const std::array<int64_t, 256> table = []() {
        std::array<int64_t, 256> values{};
        for (size_t mapq = 0; mapq < values.size(); ++mapq) {
            const float probability_correct = 1.0f - powf(
                10.0f, -(float)mapq / 10.0f);
            values[mapq] =
                (int64_t)(probability_correct * FIXED_POINT_SCALE);
        }
        for (size_t mapq = 0; mapq < values.size(); ++mapq) {
            const float old_probability_correct = 1.0f - powf(
                10.0f, -(float)mapq / 10.0f);
            const int64_t old_scaled =
                (int64_t)(old_probability_correct * FIXED_POINT_SCALE);
            if (values[mapq] != old_scaled) {
                fprintf(stderr,
                    "ERROR: MAPQ lookup validation failed at value %lu\n",
                    (unsigned long)mapq);
                abort();
            }
        }
        return values;
    }();
    return table;
}

static inline int64_t mapq_probability_scaled(uint8_t mapq) {
    return mapq_probability_scaled_table()[(size_t)mapq];
}

// Molecule-aware identity scoring sidecar.  Corrected UMI plus gene is the
// preferred RNA molecule key.  A stable query-name key is retained as an
// explicit fallback for BAMs that do not carry 10x UB/GX (or UB/GN) tags.
// Cell barcode is emitted separately and therefore is not folded into the
// molecule hash.  FNV-1a is used instead of std::hash so files are reproducible
// across processes and compiler/library upgrades.
enum PileupMoleculeBasis : uint8_t {
    PILEUP_MOLECULE_UB_GX = 1,
    PILEUP_MOLECULE_UB_GN = 2,
    PILEUP_MOLECULE_QNAME = 3,
};

struct PileupMoleculeObservation {
    unsigned long barcode = 0;
    uint64_t molecule_hash = 0;
    int64_t site = 0;
    int64_t ref_scaled = 0;
    int64_t alt_scaled = 0;
    uint8_t basis = PILEUP_MOLECULE_QNAME;
};

using PileupObservationMap = robin_hood::unordered_map<unsigned long,
    robin_hood::unordered_map<int64_t, std::pair<int64_t, int64_t> > >;

// Pileup evidence used to remain resident until the complete BAM scan ended.
// Large multiome libraries can produce hundreds of millions of observations,
// so keep only bounded per-worker chunks and stream those chunks to independent
// gzip members. The members are concatenated after a successful scan; gzip
// readers transparently consume concatenated members and downstream scoring
// already merges duplicate cell/site and molecule/site rows.
constexpr size_t PILEUP_SITE_CHUNK_ENTRIES = 100000;
constexpr size_t PILEUP_MOLECULE_CHUNK_ENTRIES = 200000;

// Raw-barcode passes cannot preallocate cell keys. A fixed set of independently
// locked shards permits concurrent insertion while retaining only one dense
// matrix per discovered barcode, rather than one copy per barcode per worker.
// The shard count is deliberately much larger than the worker count so updates
// to unrelated barcodes rarely contend.
constexpr size_t RAW_COUNT_SHARDS = 4096;
static_assert((RAW_COUNT_SHARDS & (RAW_COUNT_SHARDS - 1)) == 0,
    "RAW_COUNT_SHARDS must remain a power of two");

struct RawCountShard {
    std::mutex lock;
    std::unordered_map<unsigned long, CellCounts> panel0;
    std::unordered_map<unsigned long, CellCounts> panel1;
    std::unordered_map<unsigned long, CellCounts> native;
};

static vector<std::unique_ptr<RawCountShard>> make_raw_count_shards(bool enabled) {
    vector<std::unique_ptr<RawCountShard>> shards;
    if (!enabled) return shards;
    shards.reserve(RAW_COUNT_SHARDS);
    for (size_t i = 0; i < RAW_COUNT_SHARDS; ++i) {
        shards.emplace_back(new RawCountShard());
    }
    return shards;
}

static size_t raw_count_shard_index(unsigned long barcode) {
    return std::hash<unsigned long>()(barcode) & (RAW_COUNT_SHARDS - 1);
}

static void move_sharded_counts(
        vector<std::unique_ptr<RawCountShard>>& shards,
        int panel,
        robin_hood::unordered_map<unsigned long, AlignedCellCounts>& destination) {
    size_t source_entries = 0;
    for (const auto& shard_ptr : shards) {
        if (panel == 0) source_entries += shard_ptr->panel0.size();
        else if (panel == 1) source_entries += shard_ptr->panel1.size();
        else source_entries += shard_ptr->native.size();
    }
    destination.reserve(destination.size() + source_entries);

    size_t transferred = 0;
    for (auto& shard_ptr : shards) {
        std::unordered_map<unsigned long, CellCounts>* source = nullptr;
        if (panel == 0) source = &shard_ptr->panel0;
        else if (panel == 1) source = &shard_ptr->panel1;
        else source = &shard_ptr->native;

        // Reproduce the long-standing, validated merge path: construct the
        // destination matrix with its final dimensions and merge into it.  Do
        // not move-assign a CellCounts through an over-aligned map node.  Erase
        // each source node immediately so the one-time transfer does not hold
        // duplicate dense matrices and peak memory remains bounded.
        for (auto source_it = source->begin(); source_it != source->end();) {
            auto destination_it = destination.find(source_it->first);
            if (destination_it == destination.end()) {
                destination.emplace(std::piecewise_construct,
                    std::forward_as_tuple(source_it->first),
                    std::forward_as_tuple(source_it->second.n_samples));
                destination_it = destination.find(source_it->first);
            }
            destination_it->second.counts.merge(source_it->second);
            source_it = source->erase(source_it);
            ++transferred;
        }
    }
    fprintf(stderr, "Transferred %lu sharded raw count matrices for panel %d\n",
        (unsigned long)transferred, panel);
}

static uint64_t pileup_fnv1a_update(uint64_t hash, const char* value) {
    static const uint64_t prime = 1099511628211ULL;
    if (value != NULL) {
        for (const unsigned char* p =
                 reinterpret_cast<const unsigned char*>(value); *p; ++p) {
            hash ^= static_cast<uint64_t>(*p);
            hash *= prime;
        }
    }
    // An explicit separator prevents ("AB","C") and ("A","BC") from
    // sharing a key when two tags are combined.
    hash ^= 0xffULL;
    hash *= prime;
    return hash;
}

static std::pair<uint64_t, uint8_t> pileup_molecule_key(const bam1_t* record) {
    static const uint64_t offset = 1469598103934665603ULL;
    const uint8_t* ub_tag = bam_aux_get(record, "UB");
    const char* ub = ub_tag ? bam_aux2Z(ub_tag) : NULL;
    const uint8_t* gx_tag = bam_aux_get(record, "GX");
    const char* gx = gx_tag ? bam_aux2Z(gx_tag) : NULL;
    if (ub != NULL && ub[0] != '\0' && gx != NULL && gx[0] != '\0') {
        uint64_t hash = pileup_fnv1a_update(offset, "UB_GX");
        hash = pileup_fnv1a_update(hash, ub);
        hash = pileup_fnv1a_update(hash, gx);
        return std::make_pair(hash, PILEUP_MOLECULE_UB_GX);
    }
    const uint8_t* gn_tag = bam_aux_get(record, "GN");
    const char* gn = gn_tag ? bam_aux2Z(gn_tag) : NULL;
    if (ub != NULL && ub[0] != '\0' && gn != NULL && gn[0] != '\0') {
        uint64_t hash = pileup_fnv1a_update(offset, "UB_GN");
        hash = pileup_fnv1a_update(hash, ub);
        hash = pileup_fnv1a_update(hash, gn);
        return std::make_pair(hash, PILEUP_MOLECULE_UB_GN);
    }
    uint64_t hash = pileup_fnv1a_update(offset, "QNAME");
    hash = pileup_fnv1a_update(hash, bam_get_qname(record));
    return std::make_pair(hash, PILEUP_MOLECULE_QNAME);
}

static const char* pileup_molecule_basis_name(uint8_t basis) {
    if (basis == PILEUP_MOLECULE_UB_GX) return "UB_GX";
    if (basis == PILEUP_MOLECULE_UB_GN) return "UB_GN";
    return "QNAME_FALLBACK";
}

static long write_collapsed_pileup_molecules(
        gzFile output, vector<PileupMoleculeObservation>& observations) {
    sort(observations.begin(), observations.end(),
        [](const PileupMoleculeObservation& a,
           const PileupMoleculeObservation& b) {
            if (a.barcode != b.barcode) return a.barcode < b.barcode;
            if (a.molecule_hash != b.molecule_hash)
                return a.molecule_hash < b.molecule_hash;
            if (a.site != b.site) return a.site < b.site;
            return a.basis < b.basis;
        });
    long rows = 0;
    size_t begin = 0;
    while (begin < observations.size()) {
        size_t end = begin + 1;
        int64_t ref_scaled = observations[begin].ref_scaled;
        int64_t alt_scaled = observations[begin].alt_scaled;
        uint8_t basis = observations[begin].basis;
        while (end < observations.size() &&
                observations[end].barcode == observations[begin].barcode &&
                observations[end].molecule_hash ==
                    observations[begin].molecule_hash &&
                observations[end].site == observations[begin].site) {
            ref_scaled += observations[end].ref_scaled;
            alt_scaled += observations[end].alt_scaled;
            basis = std::min(basis, observations[end].basis);
            ++end;
        }
        const int tid = (int)(observations[begin].site >> 32);
        const int pos = (int)(observations[begin].site & 0xFFFFFFFF);
        if (gzprintf(output, "%lu\t%llu\t%s\t%d\t%d\t%f\t%f\n",
            observations[begin].barcode,
            (unsigned long long)observations[begin].molecule_hash,
            pileup_molecule_basis_name(basis), tid, pos,
            (double)ref_scaled / FIXED_POINT_SCALE,
            (double)alt_scaled / FIXED_POINT_SCALE) <= 0) {
            observations.clear();
            return -1;
        }
        ++rows;
        begin = end;
    }
    observations.clear();
    return rows;
}

static long write_collapsed_pileup_observations(
        gzFile output, PileupObservationMap& observations) {
    long rows = 0;
    for (auto& cell : observations) {
        const unsigned long barcode = cell.first;
        for (auto& site : cell.second) {
            const int tid = (int)(site.first >> 32);
            const int pos = (int)(site.first & 0xFFFFFFFF);
            if (gzprintf(output, "%lu\t%d\t%d\t%f\t%f\n",
                    barcode, tid, pos,
                    (double)site.second.first / FIXED_POINT_SCALE,
                    (double)site.second.second / FIXED_POINT_SCALE) <= 0) {
                observations.clear();
                return -1;
            }
            ++rows;
        }
    }
    observations.clear();
    return rows;
}

static std::string pileup_part_path(
        const std::string& prefix, const char* kind, int worker) {
    std::ostringstream path;
    path << prefix << "." << kind << ".part." << getpid() << "." << worker
         << ".tsv.gz";
    return path.str();
}

static void remove_files(const vector<string>& paths) {
    for (const string& path : paths) unlink(path.c_str());
}

static bool open_parallel_pileup_parts(
        const string& prefix,
        int n_threads,
        vector<string>& observation_paths,
        vector<gzFile>& observation_files,
        vector<string>& molecule_paths,
        vector<gzFile>& molecule_files,
        string& error_message) {
    observation_paths.resize(n_threads);
    molecule_paths.resize(n_threads);
    observation_files.assign(n_threads, nullptr);
    molecule_files.assign(n_threads, nullptr);
    for (int worker = 0; worker < n_threads; ++worker) {
        observation_paths[worker] = pileup_part_path(prefix, "pileup_obs", worker);
        molecule_paths[worker] = pileup_part_path(prefix, "pileup_molecules", worker);
        observation_files[worker] = gzopen(observation_paths[worker].c_str(), "wb");
        molecule_files[worker] = gzopen(molecule_paths[worker].c_str(), "wb");
        if (!observation_files[worker] || !molecule_files[worker]) {
            error_message = "could not open bounded pileup part files for worker " +
                std::to_string(worker);
            for (int i = 0; i <= worker; ++i) {
                if (observation_files[i]) gzclose(observation_files[i]);
                if (molecule_files[i]) gzclose(molecule_files[i]);
            }
            remove_files(observation_paths);
            remove_files(molecule_paths);
            return false;
        }
    }
    return true;
}

static bool close_parallel_pileup_parts(
        vector<gzFile>& files, string& error_message) {
    bool ok = true;
    for (size_t i = 0; i < files.size(); ++i) {
        if (files[i] && gzclose(files[i]) != Z_OK) {
            if (ok) {
                error_message = "failed closing bounded pileup part for worker " +
                    std::to_string(i);
            }
            ok = false;
        }
        files[i] = nullptr;
    }
    return ok;
}

static bool publish_concatenated_gzip_members(
        const vector<string>& part_paths,
        const string& final_path,
        string& error_message) {
    const string staged_path = final_path + ".tmp." +
        std::to_string((long long)getpid());
    std::ofstream output(staged_path.c_str(), std::ios::binary | std::ios::trunc);
    if (!output) {
        error_message = "could not open staged pileup output: " + staged_path;
        return false;
    }
    vector<char> buffer(1 << 20);
    for (const string& part_path : part_paths) {
        std::ifstream input(part_path.c_str(), std::ios::binary);
        if (!input) {
            error_message = "could not reopen bounded pileup part: " + part_path;
            output.close();
            unlink(staged_path.c_str());
            return false;
        }
        while (input) {
            input.read(buffer.data(), (std::streamsize)buffer.size());
            const std::streamsize bytes_read = input.gcount();
            if (bytes_read > 0) output.write(buffer.data(), bytes_read);
        }
        if (!input.eof() || !output) {
            error_message = "failed concatenating bounded pileup part: " + part_path;
            output.close();
            unlink(staged_path.c_str());
            return false;
        }
    }
    output.close();
    if (!output) {
        error_message = "failed closing staged pileup output: " + staged_path;
        unlink(staged_path.c_str());
        return false;
    }
    if (rename(staged_path.c_str(), final_path.c_str()) != 0) {
        error_message = "failed publishing pileup output: " + final_path;
        unlink(staged_path.c_str());
        return false;
    }
    return true;
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

namespace {

template <typename Policy>
auto set_policy_min_mapq(Policy& policy, uint8_t min_mapq, int)
    -> decltype(policy.min_mapq = min_mapq, void()) {
    policy.min_mapq = min_mapq;
}

template <typename Policy>
void set_policy_min_mapq(Policy&, uint8_t, long) {}

template <typename Policy>
auto get_policy_min_mapq(const Policy& policy, int)
    -> decltype(static_cast<uint8_t>(policy.min_mapq)) {
    return static_cast<uint8_t>(policy.min_mapq);
}

template <typename Policy>
uint8_t get_policy_min_mapq(const Policy&, long) {
    return 0;
}

ReadFilterPolicy& mutable_production_read_filter() {
    // Construct with defaults and assign the fields shared by both the older
    // two-argument policy and the current policy that also accepts min_mapq.
    // The small SFINAE helpers preserve source compatibility with both header
    // revisions while using min_mapq whenever the current field is available.
    static ReadFilterPolicy policy = []() {
        ReadFilterPolicy configured;
        configured.excluded_flags =
            BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
        configured.supplementary = SupplementaryReadHandling::INCLUDE;
        set_policy_min_mapq(configured, 0, 0);
        return configured;
    }();
    return policy;
}

}  // namespace

const ReadFilterPolicy& default_production_read_filter() {
    return mutable_production_read_filter();
}

void configure_read_filter(uint8_t min_mapq, uint16_t excluded_flags) {
    // demux_parallel v2.17 configures this once before any worker threads are
    // launched. Keep the effective policy in the object returned by
    // default_production_read_filter(), which is consumed by every BAM path.
    ReadFilterPolicy& policy = mutable_production_read_filter();
    policy.excluded_flags = excluded_flags;
    set_policy_min_mapq(policy, min_mapq, 0);
}

bool read_passes_filter(const bam1_t* record, const ReadFilterPolicy& policy) {
    if (record == nullptr) return false;
    uint16_t excluded = policy.excluded_flags;
    if (policy.supplementary == SupplementaryReadHandling::EXCLUDE) {
        excluded |= BAM_FSUPPLEMENTARY;
    }
    return (record->core.flag & excluded) == 0 &&
           record->core.qual >= get_policy_min_mapq(policy, 0);
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

    vector<std::pair<int, ChromSNPs*> > chromosomes;
    chromosomes.reserve(snpdat_all.size());
    long count = 0;
    for (auto& kv : snpdat_all) {
        chromosomes.push_back(std::make_pair(kv.first, &kv.second));
        count += (long)kv.second.snps.size();
    }
    ParallelOperationStatus status;
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t chromosome_index = 0;
            chromosome_index < chromosomes.size(); ++chromosome_index) {
        if (!status.ok()) continue;
        try {
            ChromSNPs& chromosome = *chromosomes[chromosome_index].second;
            for (SNPData& snp : chromosome.snps) {
                snp.precompute_genotypes(n_samples);
                snp.precompute_targets(n_samples);
            }
        } catch (const std::exception& error) {
            std::ostringstream message;
            message << "genotype/target preparation failed for TID "
                    << chromosomes[chromosome_index].first << ": "
                    << error.what();
            status.fail(message.str());
        }
    }
    if (!status.ok()) {
        fprintf(stderr, "ERROR: %s\n", status.message().c_str());
        throw std::runtime_error(status.message());
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

static inline bool checked_add_i64(int64_t& destination, int64_t value) {
    int64_t result = 0;
    if (__builtin_add_overflow(destination, value, &result)) return false;
    destination = result;
    return true;
}

static inline bool checked_species_multiplicity(
        int64_t rounded_value, uint64_t multiplicity, int64_t& result) {
    const __int128 product =
        (__int128)rounded_value * (__int128)multiplicity;
    if (product > std::numeric_limits<int64_t>::max() ||
        product < std::numeric_limits<int64_t>::min()) return false;
    result = (int64_t)product;
    return true;
}

static inline bool accumulate_species_native_targets(
    CellCounts& counts,
    const NativeSpeciesChromTargets& chrom_targets,
    size_t snp_index,
    int64_t ref_add,
    int64_t alt_add,
    uint64_t* total_target_updates = nullptr,
    uint64_t* pair_target_updates = nullptr){

    if (snp_index >= chrom_targets.site_offsets.size()) return false;
    const uint64_t offset = chrom_targets.site_offsets[snp_index];
    if (offset == UINT64_MAX) return false;
    if (snp_index >= chrom_targets.site_target_counts.size()) return false;
    const uint64_t count = chrom_targets.site_target_counts[snp_index];
    if (offset + count > chrom_targets.targets.size()) return false;
    if ((ref_add == 0) == (alt_add == 0)) return ref_add == 0;

    const bool is_ref = ref_add != 0;
    const int64_t value = is_ref ? ref_add : alt_add;
    for (uint64_t cursor = offset; cursor < offset + count; ++cursor) {
        const NativeSpeciesTargetEntry& target =
            chrom_targets.targets[(size_t)cursor];
        const int64_t add = apply_species_target_weight(value, target.weight);
        if ((target.encoded_index & NATIVE_SPECIES_PAIR_TARGET) != 0) {
            const uint32_t index =
                target.encoded_index & ~NATIVE_SPECIES_PAIR_TARGET;
            if (index >= counts.ref_counts.size()) return false;
            if (is_ref) counts.ref_counts[index] += add;
            else counts.alt_counts[index] += add;
            if (pair_target_updates) ++(*pair_target_updates);
        } else {
            const uint32_t index = target.encoded_index;
            if (index >= counts.total_ref.size()) return false;
            if (is_ref) counts.total_ref[index] += add;
            else counts.total_alt[index] += add;
            if (total_target_updates) ++(*total_target_updates);
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
    int species_native_n_samples,
    const BarcodeRemap* barcode_remap,
    BarcodeRemapStats* barcode_remap_stats){
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
    
    // A whitelist lets workers update the one preallocated matrix per cell under
    // its existing mutex. Raw-barcode discovery uses independently locked shards
    // so it can insert unknown keys without replicating matrices per worker.
    omp_set_num_threads(n_threads);
    vector<std::unique_ptr<RawCountShard>> raw_count_shards =
        make_raw_count_shards(!has_bc_list);

    std::unordered_map<unsigned long, AlignedCellCounts*> shared_count_lookup;
    std::unordered_map<unsigned long, AlignedCellCounts*> shared_species_native_lookup;
    if (has_bc_list) {
        shared_count_lookup.reserve(valid_barcodes.size());
        if (collect_species_native) {
            shared_species_native_lookup.reserve(valid_barcodes.size());
        }
        for (unsigned long barcode : valid_barcodes) {
            auto count_it = cell_counts.find(barcode);
            if (count_it == cell_counts.end()) {
                fprintf(stderr, "ERROR: internal filtered-cell pre-allocation failure\n");
                return false;
            }
            shared_count_lookup.emplace(barcode, &count_it->second);
            if (collect_species_native) {
                auto native_it = species_native_counts->find(barcode);
                if (native_it == species_native_counts->end()) {
                    fprintf(stderr, "ERROR: internal native-species pre-allocation failure\n");
                    return false;
                }
                shared_species_native_lookup.emplace(barcode, &native_it->second);
            }
        }
        fprintf(stderr,
            "Filtered-cell accumulation uses one shared count matrix per cell; "
            "per-thread dense matrices are disabled.\n");
    }
    else {
        fprintf(stderr,
            "Raw-barcode accumulation uses %lu shared count shards; per-thread "
            "dense matrices are disabled.\n",
            (unsigned long)RAW_COUNT_SHARDS);
    }

    // --dump_pileup: per-thread per-(cell,SNP) allele evidence (interindividual
    // only).  Inner key packs (tid<<32 | pos); value is (ref_scaled, alt_scaled).
    // Empty and untouched unless dump_pileup is set.
    vector<PileupObservationMap> thread_pileup(n_threads);
    vector<vector<PileupMoleculeObservation> > thread_pileup_molecules(n_threads);
    vector<size_t> thread_pileup_entries(n_threads, 0);
    vector<long> thread_pileup_rows(n_threads, 0);
    vector<long> thread_molecule_rows(n_threads, 0);
    vector<string> pileup_observation_part_paths;
    vector<string> pileup_molecule_part_paths;
    vector<gzFile> pileup_observation_part_files;
    vector<gzFile> pileup_molecule_part_files;
    string pileup_stream_error;
    if (dump_pileup && !open_parallel_pileup_parts(
            pileup_prefix, n_threads,
            pileup_observation_part_paths, pileup_observation_part_files,
            pileup_molecule_part_paths, pileup_molecule_part_files,
            pileup_stream_error)) {
        fprintf(stderr, "ERROR: %s\n", pileup_stream_error.c_str());
        return false;
    }
    if (dump_pileup) {
        fprintf(stderr,
            "Pileup memory mode: bounded worker chunks (%lu cell/site, %lu molecule/site).\n",
            (unsigned long)PILEUP_SITE_CHUNK_ENTRIES,
            (unsigned long)PILEUP_MOLECULE_CHUNK_ENTRIES);
    }
    vector<AcceptedSiteWeightMap> thread_site_weights(n_threads);
    // Optional ATAC->RNA barcode namespace remap.  Tracking is per-thread so
    // the hot counting loop remains lock-free; unique barcode sets are merged
    // after counting.
    vector<set<unsigned long>> thread_raw_barcodes(n_threads);
    vector<set<unsigned long>> thread_direct_target_barcodes(n_threads);
    vector<set<unsigned long>> thread_mapped_target_barcodes(n_threads);
    vector<set<unsigned long>> thread_map_entries_used(n_threads);
    ParallelOperationStatus operation_status;
    std::atomic<bool> hts_thread_warning_emitted(false);
    
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        auto& local_site_weights = thread_site_weights[thread_id];
        auto& local_raw_barcodes = thread_raw_barcodes[thread_id];
        auto& local_direct_barcodes = thread_direct_target_barcodes[thread_id];
        auto& local_mapped_barcodes = thread_mapped_target_barcodes[thread_id];
        auto& local_map_entries = thread_map_entries_used[thread_id];
        
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
                        const unsigned long raw_bc_key = cb_bits.to_ulong();
                        unsigned long bc_key = raw_bc_key;
                        local_raw_barcodes.insert(raw_bc_key);
                        if (!has_bc_list || valid_barcodes.find(raw_bc_key) != valid_barcodes.end()){
                            local_direct_barcodes.insert(raw_bc_key);
                        }
                        if (barcode_remap != nullptr){
                            auto remap_it = barcode_remap->find(raw_bc_key);
                            if (remap_it != barcode_remap->end()){
                                bc_key = remap_it->second;
                                local_map_entries.insert(raw_bc_key);
                            }
                        }
                        if (!has_bc_list || valid_barcodes.find(bc_key) != valid_barcodes.end()){
                            local_mapped_barcodes.insert(bc_key);
                        }
                        
                        // Skip if the direct/remapped barcode is not in the RNA-space whitelist.
                        if (has_bc_list && valid_barcodes.find(bc_key) == valid_barcodes.end()){
                            continue;
                        }
                        
                        local_reads++;

                        std::pair<uint64_t, uint8_t> molecule_key;
                        if (dump_pileup) {
                            molecule_key = pileup_molecule_key(record);
                        }
                        
                        // Get mapping quality probability and scale to fixed-point
                        int64_t prob_scaled =
                            mapq_probability_scaled(record->core.qual);
                        
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
                                // Precomputed targets: linear traversal, no branches
                                const auto& ttargets = snp_check->total_targets;
                                const auto& ptargets = snp_check->pair_targets;
                                auto add_precomputed = [&](CellCounts& counts) {
                                    const bool is_ref = ref_add != 0;
                                    for (const auto& t : ttargets){
                                        if (is_ref) counts.total_ref[t.total_idx] += ref_add;
                                        else counts.total_alt[t.total_idx] += alt_add;
                                    }
                                    for (const auto& p : ptargets){
                                        if (is_ref) counts.ref_counts[p.pair_idx] += ref_add;
                                        else counts.alt_counts[p.pair_idx] += alt_add;
                                    }
                                };

                                if (has_bc_list) {
                                    auto shared_it = shared_count_lookup.find(bc_key);
                                    if (shared_it == shared_count_lookup.end()) {
                                        operation_status.fail(format_worker_error(
                                            "filtered-cell count lookup", thread_id, tid));
                                        continue;
                                    }
                                    AlignedCellCounts& shared = *shared_it->second;
                                    std::lock_guard<std::mutex> guard(shared.lock);
                                    add_precomputed(shared.counts);
                                }
                                else {
                                    RawCountShard& shard = *raw_count_shards[
                                        raw_count_shard_index(bc_key)];
                                    std::lock_guard<std::mutex> guard(shard.lock);
                                    auto it = shard.panel0.find(bc_key);
                                    if (it == shard.panel0.end()){
                                        shard.panel0.emplace(bc_key, CellCounts(n_samples));
                                        it = shard.panel0.find(bc_key);
                                    }
                                    add_precomputed(it->second);
                                }

                                if (collect_species_native && native_chrom_targets != nullptr){
                                    const size_t snp_index = (size_t)(
                                        snp_check - chrom_snps.snps.begin());
                                    if (snp_index < native_chrom_targets->site_offsets.size() &&
                                        native_chrom_targets->site_offsets[snp_index] != UINT64_MAX){
                                        if (has_bc_list) {
                                            auto native_it = shared_species_native_lookup.find(bc_key);
                                            if (native_it == shared_species_native_lookup.end()) {
                                                operation_status.fail(format_worker_error(
                                                    "filtered native-species count lookup",
                                                    thread_id, tid));
                                                continue;
                                            }
                                            AlignedCellCounts& shared_native = *native_it->second;
                                            std::lock_guard<std::mutex> guard(shared_native.lock);
                                            accumulate_species_native_targets(
                                                shared_native.counts, *native_chrom_targets,
                                                snp_index, ref_add, alt_add);
                                        }
                                        else {
                                            RawCountShard& shard = *raw_count_shards[
                                                raw_count_shard_index(bc_key)];
                                            std::lock_guard<std::mutex> guard(shard.lock);
                                            auto native_it = shard.native.find(bc_key);
                                            if (native_it == shard.native.end()){
                                                shard.native.emplace(
                                                    bc_key, CellCounts(species_native_n_samples));
                                                native_it = shard.native.find(bc_key);
                                            }
                                            accumulate_species_native_targets(
                                                native_it->second, *native_chrom_targets, snp_index,
                                                ref_add, alt_add);
                                        }
                                    }
                                }

                                // --dump_pileup: record per-(cell,SNP) evidence for
                                // interindividual SNPs.  Summed within a thread; the
                                // producer sums any cross-thread duplicates.
                                if (dump_pileup && snp_check->panel_id == 0){
                                    int64_t pkey = ((int64_t)tid << 32) |
                                        (int64_t)(uint32_t)snp_check->pos;
                                    auto& cell_sites = thread_pileup[thread_id][bc_key];
                                    auto site_it = cell_sites.find(pkey);
                                    if (site_it == cell_sites.end()) {
                                        cell_sites.emplace(
                                            pkey, std::make_pair(ref_add, alt_add));
                                        ++thread_pileup_entries[thread_id];
                                    }
                                    else {
                                        site_it->second.first += ref_add;
                                        site_it->second.second += alt_add;
                                    }
                                    PileupMoleculeObservation observation;
                                    observation.barcode = bc_key;
                                    observation.molecule_hash = molecule_key.first;
                                    observation.site = pkey;
                                    observation.ref_scaled = ref_add;
                                    observation.alt_scaled = alt_add;
                                    observation.basis = molecule_key.second;
                                    thread_pileup_molecules[thread_id].push_back(
                                        observation);
                                    if (thread_pileup_entries[thread_id] >=
                                            PILEUP_SITE_CHUNK_ENTRIES) {
                                        const long written = write_collapsed_pileup_observations(
                                            pileup_observation_part_files[thread_id],
                                            thread_pileup[thread_id]);
                                        thread_pileup_entries[thread_id] = 0;
                                        if (written < 0) {
                                            operation_status.fail(format_worker_error(
                                                "pileup observation write", thread_id, tid));
                                        }
                                        else thread_pileup_rows[thread_id] += written;
                                    }
                                    if (thread_pileup_molecules[thread_id].size() >=
                                            PILEUP_MOLECULE_CHUNK_ENTRIES) {
                                        const long written = write_collapsed_pileup_molecules(
                                            pileup_molecule_part_files[thread_id],
                                            thread_pileup_molecules[thread_id]);
                                        if (written < 0) {
                                            operation_status.fail(format_worker_error(
                                                "pileup molecule write", thread_id, tid));
                                        }
                                        else thread_molecule_rows[thread_id] += written;
                                    }
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

        if (dump_pileup && operation_status.ok()) {
            const long obs_written = write_collapsed_pileup_observations(
                pileup_observation_part_files[thread_id], thread_pileup[thread_id]);
            thread_pileup_entries[thread_id] = 0;
            const long molecule_written = write_collapsed_pileup_molecules(
                pileup_molecule_part_files[thread_id],
                thread_pileup_molecules[thread_id]);
            if (obs_written < 0 || molecule_written < 0) {
                operation_status.fail(format_worker_error(
                    "final bounded pileup write", thread_id));
            }
            else {
                thread_pileup_rows[thread_id] += obs_written;
                thread_molecule_rows[thread_id] += molecule_written;
            }
        }
    }

    if (dump_pileup) {
        string close_error;
        const bool obs_closed = close_parallel_pileup_parts(
            pileup_observation_part_files, close_error);
        const bool molecule_closed = close_parallel_pileup_parts(
            pileup_molecule_part_files, close_error);
        if (!obs_closed || !molecule_closed) operation_status.fail(close_error);
    }

    if (!operation_status.ok()){
        fprintf(stderr, "ERROR: parallel allele counting failed: %s\n",
            operation_status.message().c_str());
        remove_files(pileup_observation_part_paths);
        remove_files(pileup_molecule_part_paths);
        return false;
    }
    
    if (!has_bc_list) {
        fprintf(stderr, "\nMoving sharded raw-barcode counts...\n");
        move_sharded_counts(raw_count_shards, 0, cell_counts);
    }
    
    if (collect_species_native && !has_bc_list){
        fprintf(stderr, "Moving sharded native species counts...\n");
        move_sharded_counts(raw_count_shards, 2, *species_native_counts);
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

    if (barcode_remap_stats != nullptr){
        set<unsigned long> raw_all, direct_all, mapped_all, used_all;
        for (int t = 0; t < n_threads; ++t){
            raw_all.insert(thread_raw_barcodes[t].begin(), thread_raw_barcodes[t].end());
            direct_all.insert(thread_direct_target_barcodes[t].begin(), thread_direct_target_barcodes[t].end());
            mapped_all.insert(thread_mapped_target_barcodes[t].begin(), thread_mapped_target_barcodes[t].end());
            used_all.insert(thread_map_entries_used[t].begin(), thread_map_entries_used[t].end());
        }
        barcode_remap_stats->observed_raw_barcodes = raw_all.size();
        barcode_remap_stats->direct_target_barcodes = direct_all.size();
        barcode_remap_stats->mapped_target_barcodes = mapped_all.size();
        barcode_remap_stats->map_entries_used = used_all.size();
    }
    thread_raw_barcodes.clear();
    thread_direct_target_barcodes.clear();
    thread_mapped_target_barcodes.clear();
    thread_map_entries_used.clear();

    raw_count_shards.clear();
    raw_count_shards.shrink_to_fit();

    // Publish the bounded per-worker gzip members. Rows can be duplicated across
    // chunks or workers; every downstream consumer merges those exact keys.
    if (dump_pileup){
        string obs_path = pileup_prefix + ".pileup_obs.tsv.gz";
        if (!publish_concatenated_gzip_members(
                pileup_observation_part_paths, obs_path, pileup_stream_error)) {
            fprintf(stderr, "ERROR: %s\n", pileup_stream_error.c_str());
            remove_files(pileup_observation_part_paths);
            remove_files(pileup_molecule_part_paths);
            return false;
        }
        string molecule_path = pileup_prefix + ".pileup_molecules.tsv.gz";
        if (!publish_concatenated_gzip_members(
                pileup_molecule_part_paths, molecule_path, pileup_stream_error)) {
            fprintf(stderr, "ERROR: %s\n", pileup_stream_error.c_str());
            remove_files(pileup_observation_part_paths);
            remove_files(pileup_molecule_part_paths);
            return false;
        }
        long n_obs_written = 0;
        long n_molecule_rows = 0;
        for (int t = 0; t < n_threads; ++t) {
            n_obs_written += thread_pileup_rows[t];
            n_molecule_rows += thread_molecule_rows[t];
        }
        remove_files(pileup_observation_part_paths);
        remove_files(pileup_molecule_part_paths);
        fprintf(stderr, "Wrote %ld bounded pileup observation rows to %s\n",
            n_obs_written, obs_path.c_str());
        fprintf(stderr, "Wrote %ld bounded molecule/SNP rows to %s\n",
            n_molecule_rows, molecule_path.c_str());
    }
    thread_pileup.clear();
    thread_pileup.shrink_to_fit();
    thread_pileup_molecules.clear();
    thread_pileup_molecules.shrink_to_fit();
    
    fprintf(stderr, "Completed: %d chromosomes (%lu work units), %ld SNPs, %ld iterator records, %lu cells\n",
        n_chroms, work_units.size(), snps_processed.load(), reads_processed.load(), 
        cell_counts.size());
    return true;
}

// ============================================================================
// FUSED FILTERED/RAW RNA COUNTING
// ============================================================================

namespace {

struct FusedAlignedObservation {
    const SNPData* snp = nullptr;
    size_t snp_index = 0;
    int64_t ref_scaled = 0;
    int64_t alt_scaled = 0;
};

enum FusedObservationPanel : uint8_t {
    FUSED_MAIN_PANEL = 0,
    FUSED_SPECIES_PANEL = 1,
};

struct FusedCountObservation {
    uint64_t barcode;
    int64_t probability_scaled;
    uint32_t snp_index;
    uint32_t tid_and_flags;

    int tid() const { return (int)(tid_and_flags & UINT32_C(0x3fffffff)); }
    uint8_t panel() const { return (uint8_t)((tid_and_flags >> 31) & 1U); }
    uint8_t allele() const { return (uint8_t)((tid_and_flags >> 30) & 1U); }
};

static_assert(sizeof(FusedCountObservation) == 24,
    "fused count observations must remain compact");

constexpr size_t FUSED_RAW_PARTITIONS = 256;
constexpr size_t FUSED_OBSERVATION_BUFFER_RECORDS = 65536;
constexpr size_t FUSED_SORT_CHUNK_RECORDS = 500000;
constexpr size_t FUSED_MERGE_FAN_IN = 24;
static_assert((FUSED_RAW_PARTITIONS & (FUSED_RAW_PARTITIONS - 1)) == 0,
    "raw partition count must remain a power of two");

static inline uint64_t fused_barcode_hash(uint64_t value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

static inline size_t fused_raw_partition(uint64_t barcode) {
    return (size_t)(fused_barcode_hash(barcode) &
        (FUSED_RAW_PARTITIONS - 1));
}

static bool fused_observation_less(
        const FusedCountObservation& left,
        const FusedCountObservation& right) {
    if (left.barcode != right.barcode) return left.barcode < right.barcode;
    if (left.panel() != right.panel()) return left.panel() < right.panel();
    if (left.tid() != right.tid()) return left.tid() < right.tid();
    if (left.snp_index != right.snp_index)
        return left.snp_index < right.snp_index;
    if (left.allele() != right.allele()) return left.allele() < right.allele();
    return left.probability_scaled < right.probability_scaled;
}

struct FusedRawPartitionFile {
    std::mutex lock;
    FILE* file = nullptr;
    string path;
    uint64_t bytes = 0;

    ~FusedRawPartitionFile() {
        if (file) fclose(file);
    }
};

class FusedTemporaryFiles {
  public:
    ~FusedTemporaryFiles() {
        std::lock_guard<std::mutex> guard(lock_);
        for (const string& path : paths_) unlink(path.c_str());
    }

    void add(const string& path) {
        std::lock_guard<std::mutex> guard(lock_);
        paths_.push_back(path);
    }

  private:
    std::mutex lock_;
    vector<string> paths_;
};

struct FusedWorkUnit {
    int tid = -1;
    int owner_start = 0;
    int owner_end = INT_MAX;
    uint64_t estimated_records = 0;
};

using FusedPerfClock = std::chrono::steady_clock;

static double fused_perf_seconds(
        const FusedPerfClock::time_point& start,
        const FusedPerfClock::time_point& end) {
    return std::chrono::duration<double>(end - start).count();
}

static uint64_t fused_perf_nanoseconds(
        const FusedPerfClock::time_point& start,
        const FusedPerfClock::time_point& end) {
    return (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(
        end - start).count();
}

static void print_fused_perf_phase(
        const char* name,
        const FusedPerfClock::time_point& start,
        const FusedPerfClock::time_point& end) {
    fprintf(stderr, "PERF_PHASE name=%s seconds=%.6f\n",
        name, fused_perf_seconds(start, end));
}

struct FusedLockSampleCounters {
    uint64_t applicable_records = 0;
    uint64_t samples = 0;
    uint64_t wait_nanoseconds = 0;
    uint64_t held_nanoseconds = 0;
};

struct FusedThreadPerfCounters {
    uint64_t iterator_records = 0;
    uint64_t records_passing_read_policy = 0;
    uint64_t cb_tagged_records = 0;
    uint64_t filtered_cell_records = 0;
    uint64_t raw_only_records = 0;
    uint64_t main_allele_observations = 0;
    uint64_t species_allele_observations = 0;
    uint64_t main_total_target_updates = 0;
    uint64_t main_pair_target_updates = 0;
    uint64_t species_total_target_updates = 0;
    uint64_t species_pair_target_updates = 0;
    uint64_t raw_shard_acquisitions = 0;
    uint64_t raw_partition_bulk_writes = 0;
    FusedLockSampleCounters filtered_cell_lock;
    FusedLockSampleCounters raw_shard_lock;
    FusedLockSampleCounters raw_partition_lock;
};

struct FusedWorkUnitProfile {
    uint64_t elapsed_nanoseconds = 0;
    uint64_t iterator_records = 0;
    uint64_t records_passing_read_policy = 0;
    uint64_t cb_tagged_records = 0;
    uint64_t filtered_cell_records = 0;
    uint64_t raw_only_records = 0;
    uint64_t main_allele_observations = 0;
    uint64_t species_allele_observations = 0;
};

constexpr uint64_t FUSED_LOCK_SAMPLE_INTERVAL = 4096;

template <typename Operation>
static inline bool with_fused_sampled_lock(
        std::mutex& lock,
        bool sample,
        FusedLockSampleCounters& counters,
        bool& sample_recorded,
        Operation operation) {
    if (!sample) {
        std::lock_guard<std::mutex> guard(lock);
        return operation();
    }

    const FusedPerfClock::time_point wait_start = FusedPerfClock::now();
    std::unique_lock<std::mutex> guard(lock);
    const FusedPerfClock::time_point held_start = FusedPerfClock::now();
    const bool ok = operation();
    const FusedPerfClock::time_point held_end = FusedPerfClock::now();
    counters.wait_nanoseconds += fused_perf_nanoseconds(wait_start, held_start);
    counters.held_nanoseconds += fused_perf_nanoseconds(held_start, held_end);
    if (!sample_recorded) {
        ++counters.samples;
        sample_recorded = true;
    }
    return ok;
}

#pragma pack(push, 1)
struct FusedPileupSpoolRecord {
    uint64_t barcode;
    uint64_t molecule_hash;
    uint64_t site;
    int64_t ref_scaled;
    int64_t alt_scaled;
    uint8_t basis;
};
#pragma pack(pop)

static_assert(sizeof(FusedPileupSpoolRecord) == 41,
    "fused pileup spool records must remain compact");

constexpr size_t FUSED_GZIP_BUFFER_BYTES = 1U << 20;
constexpr size_t FUSED_PILEUP_SORT_RECORDS = 500000;

class FusedGzipWriter {
  public:
    FusedGzipWriter() : file_(nullptr), ok_(true) {
        buffer_.reserve(FUSED_GZIP_BUFFER_BYTES + 4096);
    }

    ~FusedGzipWriter() {
        if (file_) gzclose(file_);
    }

    FusedGzipWriter(const FusedGzipWriter&) = delete;
    FusedGzipWriter& operator=(const FusedGzipWriter&) = delete;

    bool open(const string& path) {
        path_ = path;
        file_ = gzopen(path.c_str(), "wb1");
        if (!file_) {
            ok_ = false;
            return false;
        }
        if (gzbuffer(file_, (unsigned int)FUSED_GZIP_BUFFER_BYTES) != 0) {
            gzclose(file_);
            file_ = nullptr;
            ok_ = false;
            return false;
        }
        return true;
    }

    bool append(const string& text) {
        if (!ok_ || !file_) return false;
        if (buffer_.size() + text.size() > FUSED_GZIP_BUFFER_BYTES && !flush()) {
            return false;
        }
        if (text.size() > FUSED_GZIP_BUFFER_BYTES) {
            const int written = gzwrite(file_, text.data(), (unsigned int)text.size());
            if (written != (int)text.size()) ok_ = false;
            return ok_;
        }
        buffer_.append(text);
        return true;
    }

    bool close() {
        if (!file_) return ok_;
        const bool flushed = flush();
        const int close_status = gzclose(file_);
        file_ = nullptr;
        if (!flushed || close_status != Z_OK) ok_ = false;
        return ok_;
    }

    const string& path() const { return path_; }

  private:
    bool flush() {
        if (!ok_ || !file_) return false;
        if (buffer_.empty()) return true;
        const int written = gzwrite(
            file_, buffer_.data(), (unsigned int)buffer_.size());
        if (written != (int)buffer_.size()) {
            ok_ = false;
            return false;
        }
        buffer_.clear();
        return true;
    }

    gzFile file_;
    bool ok_;
    string path_;
    string buffer_;
};

static void append_scaled_decimal(string& line, int64_t value) {
    const bool negative = value < 0;
    const uint64_t magnitude = negative
        ? (uint64_t)(-(value + 1)) + 1ULL
        : (uint64_t)value;
    const uint64_t whole = magnitude / (uint64_t)FIXED_POINT_SCALE;
    const uint64_t fraction = magnitude % (uint64_t)FIXED_POINT_SCALE;
    char number[96];
    snprintf(number, sizeof(number), "%s%llu.%06llu",
        negative ? "-" : "",
        (unsigned long long)whole,
        (unsigned long long)fraction);
    line.append(number);
}

static bool append_count_row(
        FusedGzipWriter& writer,
        string& line,
        unsigned long barcode,
        int indv1,
        int nalt1,
        int indv2,
        int nalt2,
        int64_t ref_scaled,
        int64_t alt_scaled) {
    char prefix[160];
    const int length = snprintf(prefix, sizeof(prefix), "%lu\t%d\t%d\t%d\t%d\t",
        barcode, indv1, nalt1, indv2, nalt2);
    if (length < 0 || (size_t)length >= sizeof(prefix)) return false;
    line.clear();
    line.append(prefix, (size_t)length);
    append_scaled_decimal(line, ref_scaled);
    line.push_back('\t');
    append_scaled_decimal(line, alt_scaled);
    line.push_back('\n');
    return writer.append(line);
}

static bool append_dense_count_rows(
        FusedGzipWriter& writer,
        unsigned long barcode,
        const CellCounts& counts,
        int n_samples) {
    const int state_count = n_samples * GENOTYPE_STATES;
    string line;
    line.reserve(256);
    if (counts.state_count != state_count ||
        counts.total_ref.size() != (size_t)state_count ||
        counts.total_alt.size() != (size_t)state_count ||
        counts.ref_counts.size() != (size_t)state_count * (size_t)state_count ||
        counts.alt_counts.size() != (size_t)state_count * (size_t)state_count) {
        return false;
    }
    for (int indv = 0; indv < n_samples; ++indv) {
        for (int nalt = 0; nalt < GENOTYPE_STATES; ++nalt) {
            const size_t index = (size_t)indv * GENOTYPE_STATES + (size_t)nalt;
            const int64_t ref = counts.total_ref[index];
            const int64_t alt = counts.total_alt[index];
            if (ref == 0 && alt == 0) continue;
            if (!append_count_row(
                    writer, line, barcode, indv, nalt, -1, -1, ref, alt)) return false;
        }
    }
    for (int indv1 = 0; indv1 < n_samples; ++indv1) {
        for (int nalt1 = 0; nalt1 < GENOTYPE_STATES; ++nalt1) {
            const size_t idx1 =
                (size_t)indv1 * GENOTYPE_STATES + (size_t)nalt1;
            for (int indv2 = indv1 + 1; indv2 < n_samples; ++indv2) {
                for (int nalt2 = 0; nalt2 < GENOTYPE_STATES; ++nalt2) {
                    const size_t idx2 =
                        (size_t)indv2 * GENOTYPE_STATES + (size_t)nalt2;
                    const size_t index = idx1 * (size_t)state_count + idx2;
                    const int64_t ref = counts.ref_counts[index];
                    const int64_t alt = counts.alt_counts[index];
                    if (ref == 0 && alt == 0) continue;
                    if (!append_count_row(
                            writer, line, barcode, indv1, nalt1, indv2, nalt2,
                            ref, alt)) return false;
                }
            }
        }
    }
    return true;
}

static bool add_main_value_dense(
        CellCounts& counts,
        const SNPData& snp,
        bool is_alt,
        int64_t value,
        uint64_t multiplicity,
        vector<uint32_t>* touched_totals,
        vector<uint32_t>* touched_pairs,
        uint64_t& total_target_updates,
        uint64_t& pair_target_updates) {
    for (const SNPTotalTarget& target : snp.total_targets) {
        if (target.total_idx >= counts.total_ref.size()) return false;
        if (touched_totals &&
            counts.total_ref[target.total_idx] == 0 &&
            counts.total_alt[target.total_idx] == 0) {
            touched_totals->push_back(target.total_idx);
        }
        int64_t& destination = is_alt
            ? counts.total_alt[target.total_idx]
            : counts.total_ref[target.total_idx];
        if (!checked_add_i64(destination, value)) return false;
        if (__builtin_add_overflow(
                total_target_updates, multiplicity,
                &total_target_updates)) return false;
    }
    for (const SNPPairTarget& target : snp.pair_targets) {
        if (target.pair_idx >= counts.ref_counts.size()) return false;
        if (touched_pairs &&
            counts.ref_counts[target.pair_idx] == 0 &&
            counts.alt_counts[target.pair_idx] == 0) {
            touched_pairs->push_back(target.pair_idx);
        }
        int64_t& destination = is_alt
            ? counts.alt_counts[target.pair_idx]
            : counts.ref_counts[target.pair_idx];
        if (!checked_add_i64(destination, value)) return false;
        if (__builtin_add_overflow(
                pair_target_updates, multiplicity,
                &pair_target_updates)) return false;
    }
    return true;
}

static bool add_species_value_dense(
        CellCounts& counts,
        const NativeSpeciesChromTargets& targets,
        uint32_t snp_index,
        bool is_alt,
        int64_t probability_scaled,
        uint64_t multiplicity,
        vector<uint32_t>* touched_totals,
        vector<uint32_t>* touched_pairs,
        uint64_t& total_target_updates,
        uint64_t& pair_target_updates) {
    if (snp_index >= targets.site_offsets.size() ||
        snp_index >= targets.site_target_counts.size()) return false;
    const uint64_t offset = targets.site_offsets[snp_index];
    if (offset == UINT64_MAX) return true;
    const uint64_t target_count = targets.site_target_counts[snp_index];
    if (offset + target_count > targets.targets.size()) return false;
    for (uint64_t cursor = offset; cursor < offset + target_count; ++cursor) {
        const NativeSpeciesTargetEntry& target = targets.targets[(size_t)cursor];
        int64_t add = 0;
        if (!checked_species_multiplicity(
                apply_species_target_weight(probability_scaled, target.weight),
                multiplicity, add)) return false;
        if ((target.encoded_index & NATIVE_SPECIES_PAIR_TARGET) != 0) {
            const uint32_t index =
                target.encoded_index & ~NATIVE_SPECIES_PAIR_TARGET;
            if (index >= counts.ref_counts.size()) return false;
            if (add != 0 && touched_pairs && counts.ref_counts[index] == 0 &&
                counts.alt_counts[index] == 0) touched_pairs->push_back(index);
            int64_t& destination =
                is_alt ? counts.alt_counts[index] : counts.ref_counts[index];
            if (!checked_add_i64(destination, add)) return false;
            if (__builtin_add_overflow(
                    pair_target_updates, multiplicity,
                    &pair_target_updates)) return false;
        } else {
            const uint32_t index = target.encoded_index;
            if (index >= counts.total_ref.size()) return false;
            if (add != 0 && touched_totals && counts.total_ref[index] == 0 &&
                counts.total_alt[index] == 0) touched_totals->push_back(index);
            int64_t& destination =
                is_alt ? counts.total_alt[index] : counts.total_ref[index];
            if (!checked_add_i64(destination, add)) return false;
            if (__builtin_add_overflow(
                    total_target_updates, multiplicity,
                    &total_target_updates)) return false;
        }
    }
    return true;
}

static void collect_aligned_block_observations(
        const bam1_t* record,
        const ChromSNPs* panel,
        size_t& cursor,
        int64_t ref_start,
        int64_t query_start,
        int64_t length,
        int64_t probability_scaled,
        vector<FusedAlignedObservation>& observations) {
    if (!panel || length <= 0) return;
    const int64_t ref_end = ref_start + length;
    const vector<SNPData>& snps = panel->snps;
    while (cursor < snps.size() && (int64_t)snps[cursor].pos < ref_start) ++cursor;
    size_t index = cursor;
    uint8_t* sequence = bam_get_seq(const_cast<bam1_t*>(record));
    while (index < snps.size() && (int64_t)snps[index].pos < ref_end) {
        const int64_t query_index =
            query_start + ((int64_t)snps[index].pos - ref_start);
        if (query_index >= 0 && query_index < record->core.l_qseq) {
            const char allele = seq_nt16_str[bam_seqi(sequence, (int)query_index)];
            int64_t ref_scaled = 0;
            int64_t alt_scaled = 0;
            if (allele == snps[index].data.ref) ref_scaled = probability_scaled;
            else if (allele == snps[index].data.alt) alt_scaled = probability_scaled;
            if (ref_scaled != 0 || alt_scaled != 0) {
                FusedAlignedObservation observation;
                observation.snp = &snps[index];
                observation.snp_index = index;
                observation.ref_scaled = ref_scaled;
                observation.alt_scaled = alt_scaled;
                observations.push_back(observation);
            }
        }
        ++index;
    }
    cursor = index;
}

static void collect_fused_alignment_observations(
        const bam1_t* record,
        const ChromSNPs* main_panel,
        const ChromSNPs* species_panel,
        size_t main_start_cursor,
        size_t species_start_cursor,
        int64_t probability_scaled,
        vector<FusedAlignedObservation>& main_observations,
        vector<FusedAlignedObservation>& species_observations) {
    main_observations.clear();
    species_observations.clear();
    if (!record || record->core.pos < 0) return;

    uint32_t* cigar = bam_get_cigar(const_cast<bam1_t*>(record));
    int64_t ref_position = record->core.pos;
    int64_t query_position = 0;
    size_t main_cursor = main_start_cursor;
    size_t species_cursor = species_start_cursor;

    for (uint32_t cigar_index = 0;
            cigar_index < record->core.n_cigar; ++cigar_index) {
        const int operation = bam_cigar_op(cigar[cigar_index]);
        const int64_t length = bam_cigar_oplen(cigar[cigar_index]);
        if (length < 0) return;
        switch (operation) {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                if (ref_position > INT64_MAX - length ||
                    query_position > INT64_MAX - length) return;
                collect_aligned_block_observations(
                    record, main_panel, main_cursor, ref_position,
                    query_position, length, probability_scaled,
                    main_observations);
                collect_aligned_block_observations(
                    record, species_panel, species_cursor, ref_position,
                    query_position, length, probability_scaled,
                    species_observations);
                ref_position += length;
                query_position += length;
                break;
            case BAM_CINS:
            case BAM_CSOFT_CLIP:
                if (query_position > INT64_MAX - length) return;
                query_position += length;
                break;
            case BAM_CDEL:
            case BAM_CREF_SKIP:
                if (ref_position > INT64_MAX - length) return;
                ref_position += length;
                break;
            case BAM_CHARD_CLIP:
            case BAM_CPAD:
                break;
            default:
                return;
        }
        if (query_position < 0 || query_position > record->core.l_qseq) return;
    }
}

static const ChromSNPs* find_chrom_panel(
        const robin_hood::unordered_map<int, ChromSNPs>& panel,
        int tid) {
    auto found = panel.find(tid);
    return found == panel.end() ? nullptr : &found->second;
}

static bool append_compact_observation(
        vector<FusedCountObservation>& destination,
        uint64_t barcode,
        int tid,
        uint8_t panel,
        const FusedAlignedObservation& observation) {
    if (tid < 0 || (uint64_t)tid > UINT32_C(0x3fffffff) || panel > 1 ||
        observation.snp_index > std::numeric_limits<uint32_t>::max() ||
        ((observation.ref_scaled == 0) == (observation.alt_scaled == 0))) {
        return false;
    }
    FusedCountObservation compact;
    compact.barcode = barcode;
    compact.snp_index = (uint32_t)observation.snp_index;
    compact.probability_scaled = observation.ref_scaled != 0
        ? observation.ref_scaled : observation.alt_scaled;
    compact.tid_and_flags = (uint32_t)tid |
        ((uint32_t)panel << 31) |
        ((uint32_t)(observation.alt_scaled != 0 ? 1 : 0) << 30);
    destination.push_back(compact);
    return true;
}

static bool same_observation_site(
        const FusedCountObservation& left,
        const FusedCountObservation& right) {
    return left.barcode == right.barcode && left.panel() == right.panel() &&
        left.tid() == right.tid() && left.snp_index == right.snp_index &&
        left.allele() == right.allele();
}

static bool same_species_observation(
        const FusedCountObservation& left,
        const FusedCountObservation& right) {
    return same_observation_site(left, right) &&
        left.probability_scaled == right.probability_scaled;
}

static bool flush_filtered_observations(
        vector<FusedCountObservation>& observations,
        const robin_hood::unordered_map<int, ChromSNPs>& main_snpdat,
        const NativeSpeciesTargetTable& species_targets,
        const std::unordered_map<unsigned long, AlignedCellCounts*>& main_cells,
        const std::unordered_map<unsigned long, AlignedCellCounts*>& species_cells,
        FusedThreadPerfCounters& perf,
        string& error_message) {
    if (observations.empty()) return true;
    std::sort(observations.begin(), observations.end(), fused_observation_less);
    size_t begin = 0;
    while (begin < observations.size()) {
        const FusedCountObservation& first = observations[begin];
        size_t end = begin + 1;
        if (first.panel() == FUSED_MAIN_PANEL) {
            int64_t sum = first.probability_scaled;
            uint64_t multiplicity = 1;
            while (end < observations.size() &&
                    same_observation_site(first, observations[end])) {
                if (!checked_add_i64(sum, observations[end].probability_scaled)) {
                    error_message = "filtered main observation sum overflow";
                    return false;
                }
                ++multiplicity;
                ++end;
            }
            auto chromosome = main_snpdat.find(first.tid());
            auto cell = main_cells.find((unsigned long)first.barcode);
            if (chromosome == main_snpdat.end() ||
                first.snp_index >= chromosome->second.snps.size() ||
                cell == main_cells.end()) {
                error_message = "filtered main observation lookup failed";
                return false;
            }
            ++perf.filtered_cell_lock.applicable_records;
            const bool sample = perf.filtered_cell_lock.applicable_records %
                FUSED_LOCK_SAMPLE_INTERVAL == 0;
            bool sampled = false;
            const bool ok = with_fused_sampled_lock(
                cell->second->lock, sample, perf.filtered_cell_lock, sampled,
                [&]() {
                    return add_main_value_dense(
                        cell->second->counts,
                        chromosome->second.snps[first.snp_index],
                        first.allele() != 0, sum, multiplicity,
                        nullptr, nullptr,
                        perf.main_total_target_updates,
                        perf.main_pair_target_updates);
                });
            if (!ok) {
                error_message = "filtered main target accumulation overflow";
                return false;
            }
        } else {
            while (end < observations.size() &&
                    same_species_observation(first, observations[end])) ++end;
            auto chromosome = species_targets.find(first.tid());
            auto cell = species_cells.find((unsigned long)first.barcode);
            if (chromosome == species_targets.end() ||
                cell == species_cells.end()) {
                error_message = "filtered species observation lookup failed";
                return false;
            }
            ++perf.filtered_cell_lock.applicable_records;
            const bool sample = perf.filtered_cell_lock.applicable_records %
                FUSED_LOCK_SAMPLE_INTERVAL == 0;
            bool sampled = false;
            const bool ok = with_fused_sampled_lock(
                cell->second->lock, sample, perf.filtered_cell_lock, sampled,
                [&]() {
                    return add_species_value_dense(
                        cell->second->counts, chromosome->second,
                        first.snp_index, first.allele() != 0,
                        first.probability_scaled, (uint64_t)(end - begin),
                        nullptr, nullptr,
                        perf.species_total_target_updates,
                        perf.species_pair_target_updates);
                });
            if (!ok) {
                error_message = "filtered species target accumulation overflow";
                return false;
            }
        }
        begin = end;
    }
    observations.clear();
    return true;
}

static string fused_temporary_base(const string& output_path) {
    const size_t slash = output_path.find_last_of('/');
    const string directory = slash == string::npos ? "." :
        (slash == 0 ? "/" : output_path.substr(0, slash));
    const long long stamp = std::chrono::duration_cast<std::chrono::nanoseconds>(
        FusedPerfClock::now().time_since_epoch()).count();
    const string name = ".demux_raw_partition." +
        std::to_string((long long)getpid()) + "." + std::to_string(stamp);
    return directory == "/" ? directory + name : directory + "/" + name;
}

static bool open_raw_partitions(
        const string& base,
        vector<std::unique_ptr<FusedRawPartitionFile>>& partitions,
        FusedTemporaryFiles& temporary_files,
        string& error_message) {
    partitions.reserve(FUSED_RAW_PARTITIONS);
    for (size_t partition_index = 0;
            partition_index < FUSED_RAW_PARTITIONS; ++partition_index) {
        std::unique_ptr<FusedRawPartitionFile> partition(
            new FusedRawPartitionFile());
        partition->path = base + "." + std::to_string(partition_index) + ".bin";
        temporary_files.add(partition->path);
        partition->file = fopen(partition->path.c_str(), "wb");
        if (!partition->file) {
            error_message = "could not open raw observation partition " +
                partition->path + ": " + strerror(errno);
            for (auto& opened : partitions) {
                if (opened->file) fclose(opened->file);
                opened->file = nullptr;
            }
            return false;
        }
        setvbuf(partition->file, nullptr, _IOFBF, 1U << 18);
        partitions.push_back(std::move(partition));
    }
    return true;
}

static bool close_raw_partitions(
        vector<std::unique_ptr<FusedRawPartitionFile>>& partitions,
        string& error_message) {
    bool ok = true;
    for (auto& partition : partitions) {
        if (partition->file && fclose(partition->file) != 0) {
            if (ok) error_message = "failed closing raw observation partition " +
                partition->path;
            ok = false;
        }
        partition->file = nullptr;
    }
    return ok;
}

static bool spill_raw_observations(
        vector<FusedCountObservation>& observations,
        vector<FusedCountObservation>& partition_order,
        vector<std::unique_ptr<FusedRawPartitionFile>>& partitions,
        FusedThreadPerfCounters& perf,
        string& error_message) {
    if (observations.empty()) return true;
    std::array<size_t, FUSED_RAW_PARTITIONS> counts{};
    std::array<size_t, FUSED_RAW_PARTITIONS> offsets{};
    std::array<size_t, FUSED_RAW_PARTITIONS> cursors{};
    for (const FusedCountObservation& observation : observations)
        ++counts[fused_raw_partition(observation.barcode)];
    size_t prefix = 0;
    for (size_t i = 0; i < FUSED_RAW_PARTITIONS; ++i) {
        offsets[i] = prefix;
        cursors[i] = prefix;
        prefix += counts[i];
    }
    partition_order.resize(observations.size());
    for (const FusedCountObservation& observation : observations) {
        const size_t owner = fused_raw_partition(observation.barcode);
        partition_order[cursors[owner]++] = observation;
    }
    for (size_t partition_index = 0;
            partition_index < FUSED_RAW_PARTITIONS; ++partition_index) {
        if (counts[partition_index] == 0) continue;
        const size_t begin = offsets[partition_index];
        const size_t end = begin + counts[partition_index];
        FusedRawPartitionFile& partition = *partitions[partition_index];
        ++perf.raw_partition_lock.applicable_records;
        const bool sample = perf.raw_partition_lock.applicable_records %
            FUSED_LOCK_SAMPLE_INTERVAL == 0;
        bool sampled = false;
        const bool ok = with_fused_sampled_lock(
            partition.lock, sample, perf.raw_partition_lock, sampled, [&]() {
                const size_t count = end - begin;
                if (fwrite(partition_order.data() + begin,
                        sizeof(FusedCountObservation), count,
                        partition.file) != count) return false;
                partition.bytes +=
                    (uint64_t)count * sizeof(FusedCountObservation);
                return true;
            });
        if (!ok) {
            error_message = "failed writing raw observation partition";
            return false;
        }
        ++perf.raw_partition_bulk_writes;
    }
    observations.clear();
    partition_order.clear();
    return true;
}

static vector<int> fused_rank_boundaries(
        const ChromSNPs* main_panel,
        const ChromSNPs* species_panel,
        size_t chunks,
        int64_t chrom_length) {
    vector<int> boundaries;
    if (chunks <= 1) return boundaries;
    const size_t main_size = main_panel ? main_panel->snps.size() : 0;
    const size_t species_size = species_panel ? species_panel->snps.size() : 0;
    const size_t total = main_size + species_size;
    if (total == 0) return boundaries;

    size_t main_index = 0;
    size_t species_index = 0;
    size_t rank = 0;
    size_t boundary_number = 1;
    size_t target_rank = (total * boundary_number) / chunks;
    int previous = 0;
    while ((main_index < main_size || species_index < species_size) &&
            boundary_number < chunks) {
        int position = 0;
        if (species_index >= species_size ||
            (main_index < main_size &&
             main_panel->snps[main_index].pos <=
                species_panel->snps[species_index].pos)) {
            position = main_panel->snps[main_index++].pos;
        } else {
            position = species_panel->snps[species_index++].pos;
        }
        if (rank >= target_rank) {
            if (position > previous && position > 0 && position < chrom_length) {
                boundaries.push_back(position);
                previous = position;
            }
            ++boundary_number;
            target_rank = (total * boundary_number) / chunks;
        }
        ++rank;
    }
    return boundaries;
}

static vector<int> fused_position_boundaries(size_t chunks, int64_t chrom_length) {
    vector<int> boundaries;
    if (chunks <= 1 || chrom_length <= 1) return boundaries;
    int previous = 0;
    for (size_t chunk = 1; chunk < chunks; ++chunk) {
        const int64_t boundary =
            (chrom_length * (int64_t)chunk) / (int64_t)chunks;
        if (boundary > previous && boundary > 0 && boundary < chrom_length &&
            boundary < INT_MAX) {
            boundaries.push_back((int)boundary);
            previous = (int)boundary;
        }
    }
    return boundaries;
}

static bool write_fused_pileup_sites(
        const string& path,
        const robin_hood::unordered_map<int, ChromSNPs>& main_snpdat,
        const vector<string>& chromosome_names,
        int n_samples) {
    FusedGzipWriter writer;
    if (!writer.open(path)) return false;
    for (const auto& chromosome : main_snpdat) {
        const int tid = chromosome.first;
        const string name =
            tid >= 0 && tid < (int)chromosome_names.size()
                ? chromosome_names[(size_t)tid] : ".";
        for (const SNPData& snp : chromosome.second.snps) {
            if ((int)snp.geno.size() < n_samples) {
                writer.close();
                unlink(path.c_str());
                return false;
            }
            char prefix[256];
            const int length = snprintf(
                prefix, sizeof(prefix), "%d\t%s\t%d\t%c\t%c",
                tid, name.c_str(), snp.pos, snp.data.ref, snp.data.alt);
            if (length < 0 || (size_t)length >= sizeof(prefix)) {
                writer.close();
                unlink(path.c_str());
                return false;
            }
            string line(prefix, (size_t)length);
            for (int sample = 0; sample < n_samples; ++sample) {
                char genotype[32];
                const int genotype_length = snprintf(
                    genotype, sizeof(genotype), "\t%d", (int)snp.geno[sample]);
                if (genotype_length < 0 ||
                    (size_t)genotype_length >= sizeof(genotype)) {
                    writer.close();
                    unlink(path.c_str());
                    return false;
                }
                line.append(genotype, (size_t)genotype_length);
            }
            line.push_back('\n');
            if (!writer.append(line)) {
                writer.close();
                unlink(path.c_str());
                return false;
            }
        }
    }
    if (!writer.close()) {
        unlink(path.c_str());
        return false;
    }
    return true;
}

static string fused_spool_path(const string& prefix, int worker) {
    std::ostringstream path;
    path << prefix << ".pileup_spool." << getpid() << "." << worker << ".bin";
    return path.str();
}

static bool open_fused_spools(
        const string& prefix,
        int n_threads,
        vector<string>& paths,
        vector<FILE*>& files) {
    paths.resize((size_t)n_threads);
    files.assign((size_t)n_threads, nullptr);
    for (int worker = 0; worker < n_threads; ++worker) {
        paths[(size_t)worker] = fused_spool_path(prefix, worker);
        files[(size_t)worker] = fopen(paths[(size_t)worker].c_str(), "wb");
        if (!files[(size_t)worker] ||
            setvbuf(files[(size_t)worker], nullptr, _IOFBF,
                FUSED_GZIP_BUFFER_BYTES) != 0) {
            for (int opened = 0; opened <= worker; ++opened) {
                if (files[(size_t)opened]) fclose(files[(size_t)opened]);
                files[(size_t)opened] = nullptr;
            }
            remove_files(paths);
            return false;
        }
    }
    return true;
}

static bool close_fused_spools(vector<FILE*>& files) {
    bool ok = true;
    for (FILE*& file : files) {
        if (file && fclose(file) != 0) ok = false;
        file = nullptr;
    }
    return ok;
}

static bool append_pileup_observation_row(
        FusedGzipWriter& writer,
        const FusedPileupSpoolRecord& key,
        int64_t ref_scaled,
        int64_t alt_scaled) {
    const int tid = (int)(key.site >> 32);
    const int pos = (int)(key.site & 0xFFFFFFFFULL);
    char prefix[160];
    const int length = snprintf(prefix, sizeof(prefix), "%llu\t%d\t%d\t",
        (unsigned long long)key.barcode, tid, pos);
    if (length < 0 || (size_t)length >= sizeof(prefix)) return false;
    string line(prefix, (size_t)length);
    append_scaled_decimal(line, ref_scaled);
    line.push_back('\t');
    append_scaled_decimal(line, alt_scaled);
    line.push_back('\n');
    return writer.append(line);
}

static bool append_pileup_molecule_row(
        FusedGzipWriter& writer,
        const FusedPileupSpoolRecord& key,
        uint8_t basis,
        int64_t ref_scaled,
        int64_t alt_scaled) {
    const int tid = (int)(key.site >> 32);
    const int pos = (int)(key.site & 0xFFFFFFFFULL);
    char prefix[256];
    const int length = snprintf(
        prefix, sizeof(prefix), "%llu\t%llu\t%s\t%d\t%d\t",
        (unsigned long long)key.barcode,
        (unsigned long long)key.molecule_hash,
        pileup_molecule_basis_name(basis), tid, pos);
    if (length < 0 || (size_t)length >= sizeof(prefix)) return false;
    string line(prefix, (size_t)length);
    append_scaled_decimal(line, ref_scaled);
    line.push_back('\t');
    append_scaled_decimal(line, alt_scaled);
    line.push_back('\n');
    return writer.append(line);
}

static bool derive_fused_pileup_outputs(
        const vector<string>& spool_paths,
        const string& observation_path,
        const string& molecule_path) {
    FusedGzipWriter observation_writer;
    FusedGzipWriter molecule_writer;
    if (!observation_writer.open(observation_path) ||
        !molecule_writer.open(molecule_path)) {
        observation_writer.close();
        molecule_writer.close();
        unlink(observation_path.c_str());
        unlink(molecule_path.c_str());
        remove_files(spool_paths);
        return false;
    }

    vector<FusedPileupSpoolRecord> records(FUSED_PILEUP_SORT_RECORDS);
    bool ok = true;
    for (const string& spool_path : spool_paths) {
        FILE* input = fopen(spool_path.c_str(), "rb");
        if (!input) {
            ok = false;
            break;
        }
        while (ok) {
            const size_t n_records = fread(
                records.data(), sizeof(FusedPileupSpoolRecord),
                records.size(), input);
            if (n_records == 0) break;
            std::sort(records.begin(), records.begin() + n_records,
                [](const FusedPileupSpoolRecord& left,
                   const FusedPileupSpoolRecord& right) {
                    if (left.barcode != right.barcode)
                        return left.barcode < right.barcode;
                    if (left.site != right.site) return left.site < right.site;
                    if (left.molecule_hash != right.molecule_hash)
                        return left.molecule_hash < right.molecule_hash;
                    return left.basis < right.basis;
                });

            size_t observation_begin = 0;
            while (observation_begin < n_records && ok) {
                size_t observation_end = observation_begin + 1;
                int64_t observation_ref = records[observation_begin].ref_scaled;
                int64_t observation_alt = records[observation_begin].alt_scaled;
                while (observation_end < n_records &&
                        records[observation_end].barcode ==
                            records[observation_begin].barcode &&
                        records[observation_end].site ==
                            records[observation_begin].site) {
                    observation_ref += records[observation_end].ref_scaled;
                    observation_alt += records[observation_end].alt_scaled;
                    ++observation_end;
                }
                ok = append_pileup_observation_row(
                    observation_writer, records[observation_begin],
                    observation_ref, observation_alt);

                size_t molecule_begin = observation_begin;
                while (molecule_begin < observation_end && ok) {
                    size_t molecule_end = molecule_begin + 1;
                    int64_t molecule_ref = records[molecule_begin].ref_scaled;
                    int64_t molecule_alt = records[molecule_begin].alt_scaled;
                    uint8_t basis = records[molecule_begin].basis;
                    while (molecule_end < observation_end &&
                            records[molecule_end].molecule_hash ==
                                records[molecule_begin].molecule_hash) {
                        molecule_ref += records[molecule_end].ref_scaled;
                        molecule_alt += records[molecule_end].alt_scaled;
                        basis = std::min(basis, records[molecule_end].basis);
                        ++molecule_end;
                    }
                    ok = append_pileup_molecule_row(
                        molecule_writer, records[molecule_begin], basis,
                        molecule_ref, molecule_alt);
                    molecule_begin = molecule_end;
                }
                observation_begin = observation_end;
            }
        }
        if (ferror(input)) ok = false;
        if (fclose(input) != 0) ok = false;
        if (!ok) break;
    }

    if (!observation_writer.close()) ok = false;
    if (!molecule_writer.close()) ok = false;
    remove_files(spool_paths);
    if (!ok) {
        unlink(observation_path.c_str());
        unlink(molecule_path.c_str());
    }
    return ok;
}

class FusedRawBarcodeAccumulator {
  public:
    FusedRawBarcodeAccumulator(
            int n_samples,
            int n_species,
            const robin_hood::unordered_map<int, ChromSNPs>& main_panel,
            const NativeSpeciesTargetTable& species_targets,
            FusedGzipWriter& main_writer,
            FusedGzipWriter& species_writer,
            FusedThreadPerfCounters& perf,
            size_t owner_partition)
        : main_counts_(n_samples), species_counts_(n_species),
          main_panel_(main_panel), species_targets_(species_targets),
          main_writer_(main_writer), species_writer_(species_writer), perf_(perf),
          owner_partition_(owner_partition) {
        line_.reserve(256);
    }

    bool consume(const FusedCountObservation& observation) {
        if (fused_raw_partition(observation.barcode) != owner_partition_ ||
            observation.panel() > FUSED_SPECIES_PANEL || observation.allele() > 1 ||
            observation.probability_scaled <= 0) return false;
        if (!has_pending_) {
            pending_ = observation;
            multiplicity_ = 1;
            has_pending_ = true;
            return true;
        }
        const bool same_group = pending_.panel() == FUSED_MAIN_PANEL
            ? same_observation_site(pending_, observation)
            : same_species_observation(pending_, observation);
        if (same_group) {
            if (pending_.panel() == FUSED_MAIN_PANEL) {
                if (!checked_add_i64(
                        pending_.probability_scaled,
                        observation.probability_scaled)) return false;
            }
            if (multiplicity_ == std::numeric_limits<uint64_t>::max())
                return false;
            ++multiplicity_;
            return true;
        }

        const uint64_t completed_barcode = pending_.barcode;
        if (!flush_group()) return false;
        if (observation.barcode != completed_barcode &&
            !write_and_reset(completed_barcode)) return false;
        pending_ = observation;
        multiplicity_ = 1;
        has_pending_ = true;
        return true;
    }

    bool finish() {
        if (!has_pending_) return true;
        const uint64_t barcode = pending_.barcode;
        if (!flush_group() || !write_and_reset(barcode)) return false;
        has_pending_ = false;
        return true;
    }

    uint64_t barcode_count() const { return barcode_count_; }

  private:
    bool flush_group() {
        if (!has_pending_) return true;
        bool ok = false;
        if (pending_.panel() == FUSED_MAIN_PANEL) {
            auto chromosome = main_panel_.find(pending_.tid());
            if (chromosome == main_panel_.end() ||
                pending_.snp_index >= chromosome->second.snps.size()) return false;
            ok = add_main_value_dense(
                main_counts_, chromosome->second.snps[pending_.snp_index],
                pending_.allele() != 0, pending_.probability_scaled,
                multiplicity_, &main_touched_totals_, &main_touched_pairs_,
                perf_.main_total_target_updates,
                perf_.main_pair_target_updates);
        } else {
            auto chromosome = species_targets_.find(pending_.tid());
            if (chromosome == species_targets_.end()) return false;
            ok = add_species_value_dense(
                species_counts_, chromosome->second, pending_.snp_index,
                pending_.allele() != 0, pending_.probability_scaled,
                multiplicity_, &species_touched_totals_,
                &species_touched_pairs_, perf_.species_total_target_updates,
                perf_.species_pair_target_updates);
        }
        has_pending_ = false;
        return ok;
    }

    bool append_touched(
            FusedGzipWriter& writer,
            uint64_t barcode,
            CellCounts& counts,
            vector<uint32_t>& totals,
            vector<uint32_t>& pairs) {
        for (uint32_t index : totals) {
            if (index >= counts.total_ref.size()) return false;
            const int64_t ref = counts.total_ref[index];
            const int64_t alt = counts.total_alt[index];
            if ((ref != 0 || alt != 0) && !append_count_row(
                    writer, line_, (unsigned long)barcode,
                    (int)(index / GENOTYPE_STATES),
                    (int)(index % GENOTYPE_STATES), -1, -1, ref, alt)) {
                return false;
            }
        }
        const size_t state_count = (size_t)counts.state_count;
        for (uint32_t index : pairs) {
            if (index >= counts.ref_counts.size()) return false;
            const int64_t ref = counts.ref_counts[index];
            const int64_t alt = counts.alt_counts[index];
            if (ref == 0 && alt == 0) continue;
            const size_t idx1 = index / state_count;
            const size_t idx2 = index % state_count;
            if (!append_count_row(
                    writer, line_, (unsigned long)barcode,
                    (int)(idx1 / GENOTYPE_STATES),
                    (int)(idx1 % GENOTYPE_STATES),
                    (int)(idx2 / GENOTYPE_STATES),
                    (int)(idx2 % GENOTYPE_STATES), ref, alt)) return false;
        }
        return true;
    }

    static void reset_touched(
            CellCounts& counts,
            vector<uint32_t>& totals,
            vector<uint32_t>& pairs) {
        for (uint32_t index : totals) {
            counts.total_ref[index] = 0;
            counts.total_alt[index] = 0;
        }
        for (uint32_t index : pairs) {
            counts.ref_counts[index] = 0;
            counts.alt_counts[index] = 0;
        }
        totals.clear();
        pairs.clear();
    }

    bool write_and_reset(uint64_t barcode) {
        const bool ok = append_touched(
                main_writer_, barcode, main_counts_, main_touched_totals_,
                main_touched_pairs_) &&
            append_touched(
                species_writer_, barcode, species_counts_,
                species_touched_totals_, species_touched_pairs_);
        if (!ok) return false;
        reset_touched(
            main_counts_, main_touched_totals_, main_touched_pairs_);
        reset_touched(
            species_counts_, species_touched_totals_, species_touched_pairs_);
        ++barcode_count_;
        return true;
    }

    CellCounts main_counts_;
    CellCounts species_counts_;
    vector<uint32_t> main_touched_totals_;
    vector<uint32_t> main_touched_pairs_;
    vector<uint32_t> species_touched_totals_;
    vector<uint32_t> species_touched_pairs_;
    const robin_hood::unordered_map<int, ChromSNPs>& main_panel_;
    const NativeSpeciesTargetTable& species_targets_;
    FusedGzipWriter& main_writer_;
    FusedGzipWriter& species_writer_;
    FusedThreadPerfCounters& perf_;
    size_t owner_partition_;
    string line_;
    FusedCountObservation pending_{};
    uint64_t multiplicity_ = 0;
    uint64_t barcode_count_ = 0;
    bool has_pending_ = false;
};

static bool write_sorted_observation_run(
        const string& path, vector<FusedCountObservation>& records) {
    std::sort(records.begin(), records.end(), fused_observation_less);
    FILE* output = fopen(path.c_str(), "wb");
    if (!output) return false;
    const bool ok = records.empty() || fwrite(
        records.data(), sizeof(FusedCountObservation), records.size(), output) ==
            records.size();
    const bool closed = fclose(output) == 0;
    return ok && closed;
}

static bool merge_observation_runs(
        const vector<string>& paths,
        const string* output_path,
        FusedRawBarcodeAccumulator* reducer) {
    struct Cursor { FILE* file = nullptr; FusedCountObservation record{}; };
    struct HeapItem { FusedCountObservation record; size_t cursor; };
    struct HeapGreater {
        bool operator()(const HeapItem& left, const HeapItem& right) const {
            if (fused_observation_less(right.record, left.record)) return true;
            if (fused_observation_less(left.record, right.record)) return false;
            return left.cursor > right.cursor;
        }
    };
    vector<Cursor> cursors(paths.size());
    std::priority_queue<HeapItem, vector<HeapItem>, HeapGreater> heap;
    FILE* output = nullptr;
    bool ok = output_path != nullptr || reducer != nullptr;
    if (output_path) output = fopen(output_path->c_str(), "wb");
    if (output_path && !output) ok = false;
    for (size_t i = 0; i < paths.size() && ok; ++i) {
        cursors[i].file = fopen(paths[i].c_str(), "rb");
        if (!cursors[i].file) { ok = false; break; }
        if (fread(&cursors[i].record, sizeof(FusedCountObservation), 1,
                cursors[i].file) == 1) {
            heap.push({cursors[i].record, i});
        } else if (ferror(cursors[i].file)) { ok = false; break; }
    }
    while (ok && !heap.empty()) {
        const HeapItem item = heap.top();
        heap.pop();
        ok = output
            ? fwrite(&item.record, sizeof(item.record), 1, output) == 1
            : reducer->consume(item.record);
        Cursor& cursor = cursors[item.cursor];
        if (ok && fread(&cursor.record, sizeof(cursor.record), 1,
                cursor.file) == 1) {
            heap.push({cursor.record, item.cursor});
        } else if (ok && ferror(cursor.file)) ok = false;
    }
    for (Cursor& cursor : cursors) {
        if (cursor.file && fclose(cursor.file) != 0) ok = false;
    }
    if (output && fclose(output) != 0) ok = false;
    if (!output && ok) ok = reducer->finish();
    return ok;
}

static bool reduce_raw_partition(
        const FusedRawPartitionFile& partition,
        size_t partition_index,
        const string& temporary_base,
        FusedTemporaryFiles& temporary_files,
        const robin_hood::unordered_map<int, ChromSNPs>& main_panel,
        const NativeSpeciesTargetTable& species_targets,
        int n_samples,
        int n_species,
        const string& main_member,
        const string& species_member,
        FusedThreadPerfCounters& perf,
        uint64_t& barcode_count,
        string& error_message) {
    FusedGzipWriter main_writer;
    FusedGzipWriter species_writer;
    if (!main_writer.open(main_member) || !species_writer.open(species_member)) {
        main_writer.close(); species_writer.close();
        error_message = "could not open raw gzip member";
        return false;
    }
    FusedRawBarcodeAccumulator reducer(
        n_samples, n_species, main_panel, species_targets,
        main_writer, species_writer, perf, partition_index);
    FILE* input = fopen(partition.path.c_str(), "rb");
    if (!input) {
        error_message = "could not reopen raw observation partition";
        main_writer.close(); species_writer.close();
        return false;
    }
    if (partition.bytes % sizeof(FusedCountObservation) != 0) {
        fclose(input);
        main_writer.close(); species_writer.close();
        error_message = "raw observation partition has a partial record";
        return false;
    }
    if (partition.bytes <=
            FUSED_SORT_CHUNK_RECORDS * sizeof(FusedCountObservation)) {
        vector<FusedCountObservation> in_memory(
            (size_t)(partition.bytes / sizeof(FusedCountObservation)));
        bool ok = in_memory.empty() || fread(
            in_memory.data(), sizeof(FusedCountObservation), in_memory.size(),
            input) == in_memory.size();
        if (fclose(input) != 0) ok = false;
        if (ok) {
            std::sort(in_memory.begin(), in_memory.end(), fused_observation_less);
            for (const FusedCountObservation& observation : in_memory) {
                if (!reducer.consume(observation)) { ok = false; break; }
            }
            if (ok) ok = reducer.finish();
        }
        if (!main_writer.close()) ok = false;
        if (!species_writer.close()) ok = false;
        if (!ok) {
            error_message = "in-memory raw observation reduction failed";
            unlink(main_member.c_str()); unlink(species_member.c_str());
            return false;
        }
        barcode_count = reducer.barcode_count();
        unlink(partition.path.c_str());
        return true;
    }
    vector<string> runs;
    vector<FusedCountObservation> records(FUSED_SORT_CHUNK_RECORDS);
    size_t run_index = 0;
    bool ok = true;
    while (ok) {
        const size_t count = fread(
            records.data(), sizeof(FusedCountObservation), records.size(), input);
        if (count == 0) { if (ferror(input)) ok = false; break; }
        records.resize(count);
        const string run_path = temporary_base + ".p" +
            std::to_string(partition_index) + ".r" +
            std::to_string(run_index++) + ".bin";
        temporary_files.add(run_path);
        if (!write_sorted_observation_run(run_path, records)) { ok = false; break; }
        runs.push_back(run_path);
        records.resize(FUSED_SORT_CHUNK_RECORDS);
    }
    if (fclose(input) != 0) ok = false;
    if (ok) unlink(partition.path.c_str());
    size_t generation = 0;
    while (ok && runs.size() > FUSED_MERGE_FAN_IN) {
        vector<string> merged_runs;
        for (size_t begin = 0; begin < runs.size();
                begin += FUSED_MERGE_FAN_IN) {
            const size_t end = std::min(
                runs.size(), begin + FUSED_MERGE_FAN_IN);
            vector<string> batch(runs.begin() + begin, runs.begin() + end);
            const string merged_path = temporary_base + ".p" +
                std::to_string(partition_index) + ".g" +
                std::to_string(generation) + ".m" +
                std::to_string(merged_runs.size()) + ".bin";
            temporary_files.add(merged_path);
            if (!merge_observation_runs(batch, &merged_path, nullptr)) {
                ok = false; break;
            }
            for (const string& path : batch) unlink(path.c_str());
            merged_runs.push_back(merged_path);
        }
        runs.swap(merged_runs);
        ++generation;
    }
    if (ok) ok = runs.empty()
        ? reducer.finish() : merge_observation_runs(runs, nullptr, &reducer);
    for (const string& path : runs) unlink(path.c_str());
    if (!main_writer.close()) ok = false;
    if (!species_writer.close()) ok = false;
    if (!ok) {
        error_message = "raw observation reduction or serialization failed";
        unlink(main_member.c_str()); unlink(species_member.c_str());
        return false;
    }
    barcode_count = reducer.barcode_count();
    unlink(partition.path.c_str());
    return true;
}

static bool write_fused_raw_counts(
        const string& main_path,
        const string& species_path,
        const robin_hood::unordered_map<unsigned long, AlignedCellCounts>& filtered_main,
        const robin_hood::unordered_map<unsigned long, AlignedCellCounts>& filtered_species,
        vector<std::unique_ptr<FusedRawPartitionFile>>& partitions,
        FusedTemporaryFiles& temporary_files,
        const string& temporary_base,
        const robin_hood::unordered_map<int, ChromSNPs>& main_panel,
        const NativeSpeciesTargetTable& species_targets,
        int n_samples,
        int n_species,
        int n_threads,
        FusedThreadPerfCounters& reduction_perf,
        uint64_t& raw_barcode_count) {
    vector<string> main_members(FUSED_RAW_PARTITIONS + 1);
    vector<string> species_members(FUSED_RAW_PARTITIONS + 1);
    main_members[0] = temporary_base + ".filtered.main.gz";
    species_members[0] = temporary_base + ".filtered.species.gz";
    temporary_files.add(main_members[0]);
    temporary_files.add(species_members[0]);
    FusedGzipWriter filtered_main_writer;
    FusedGzipWriter filtered_species_writer;
    bool ok = filtered_main_writer.open(main_members[0]) &&
        filtered_species_writer.open(species_members[0]);
    const FusedPerfClock::time_point filtered_main_start = FusedPerfClock::now();
    if (ok) {
        for (const auto& cell : filtered_main) {
            if (!append_dense_count_rows(
                    filtered_main_writer, cell.first, cell.second.counts,
                    n_samples)) { ok = false; break; }
        }
    }
    print_fused_perf_phase(
        "filtered_main_raw_serialization", filtered_main_start,
        FusedPerfClock::now());

    const FusedPerfClock::time_point filtered_species_start = FusedPerfClock::now();
    if (ok) {
        for (const auto& cell : filtered_species) {
            if (!append_dense_count_rows(
                    filtered_species_writer, cell.first, cell.second.counts,
                    n_species)) { ok = false; break; }
        }
    }
    print_fused_perf_phase(
        "filtered_species_raw_serialization", filtered_species_start,
        FusedPerfClock::now());

    if (!filtered_main_writer.close()) ok = false;
    if (!filtered_species_writer.close()) ok = false;
    if (!ok) {
        unlink(main_path.c_str()); unlink(species_path.c_str());
        return false;
    }

    vector<FusedThreadPerfCounters> partition_perf(FUSED_RAW_PARTITIONS);
    vector<uint64_t> partition_barcodes(FUSED_RAW_PARTITIONS, 0);
    vector<string> partition_errors(FUSED_RAW_PARTITIONS);
    std::atomic<bool> reduction_ok(true);
    for (size_t partition_index = 0;
            partition_index < FUSED_RAW_PARTITIONS; ++partition_index) {
        main_members[partition_index + 1] = temporary_base + ".p" +
            std::to_string(partition_index) + ".main.gz";
        species_members[partition_index + 1] = temporary_base + ".p" +
            std::to_string(partition_index) + ".species.gz";
        temporary_files.add(main_members[partition_index + 1]);
        temporary_files.add(species_members[partition_index + 1]);
    }

    const FusedPerfClock::time_point sparse_raw_start = FusedPerfClock::now();
    const int reduction_threads = std::max(1, std::min(n_threads, 16));
    #pragma omp parallel for schedule(dynamic, 1) num_threads(reduction_threads)
    for (size_t partition_index = 0;
            partition_index < FUSED_RAW_PARTITIONS; ++partition_index) {
        if (!reduction_ok.load(std::memory_order_acquire)) continue;
        bool partition_ok = false;
        try {
            partition_ok = reduce_raw_partition(
                *partitions[partition_index], partition_index, temporary_base,
                temporary_files, main_panel, species_targets,
                n_samples, n_species, main_members[partition_index + 1],
                species_members[partition_index + 1],
                partition_perf[partition_index],
                partition_barcodes[partition_index],
                partition_errors[partition_index]);
        } catch (const std::exception& error) {
            partition_errors[partition_index] =
                string("raw partition worker exception: ") + error.what();
        }
        if (!partition_ok) {
            reduction_ok.store(false, std::memory_order_release);
        }
    }
    print_fused_perf_phase(
        "sparse_raw_only_main_species_serialization", sparse_raw_start,
        FusedPerfClock::now());
    if (!reduction_ok.load(std::memory_order_acquire)) {
        for (const string& error : partition_errors) {
            if (!error.empty()) { fprintf(stderr, "ERROR: %s\n", error.c_str()); break; }
        }
        unlink(main_path.c_str()); unlink(species_path.c_str());
        return false;
    }

    for (size_t i = 0; i < FUSED_RAW_PARTITIONS; ++i) {
        reduction_perf.main_total_target_updates +=
            partition_perf[i].main_total_target_updates;
        reduction_perf.main_pair_target_updates +=
            partition_perf[i].main_pair_target_updates;
        reduction_perf.species_total_target_updates +=
            partition_perf[i].species_total_target_updates;
        reduction_perf.species_pair_target_updates +=
            partition_perf[i].species_pair_target_updates;
        raw_barcode_count += partition_barcodes[i];
    }

    const FusedPerfClock::time_point gzip_close_start = FusedPerfClock::now();
    string error;
    ok = publish_concatenated_gzip_members(main_members, main_path, error) &&
        publish_concatenated_gzip_members(species_members, species_path, error);
    print_fused_perf_phase(
        "gzip_flush_close", gzip_close_start, FusedPerfClock::now());
    if (!ok) {
        fprintf(stderr, "ERROR: %s\n", error.c_str());
        unlink(main_path.c_str()); unlink(species_path.c_str());
    }
    return ok;
}

}  // namespace

bool count_alleles_parallel_fused(
    const string& bamfile,
    const robin_hood::unordered_map<int, ChromSNPs>& main_snpdat,
    const robin_hood::unordered_map<int, ChromSNPs>& species_snpdat,
    robin_hood::unordered_map<unsigned long, AlignedCellCounts>& filtered_main_counts,
    robin_hood::unordered_map<unsigned long, AlignedCellCounts>& filtered_species_counts,
    const set<unsigned long>& filtered_barcodes,
    int n_samples,
    const NativeSpeciesTargetTable& species_native_targets,
    int species_native_n_samples,
    int n_threads,
    int htslib_threads,
    const string& raw_counts_path,
    const string& raw_species_counts_path,
    bool dump_pileup,
    const string& pileup_prefix,
    AcceptedSiteWeightMap* accepted_site_weights_main,
    AcceptedSiteWeightMap* accepted_site_weights_species) {
    string validation_error;
    size_t main_bytes_per_cell = 0;
    size_t species_bytes_per_cell = 0;
    const bool write_raw_outputs =
        !raw_counts_path.empty() || !raw_species_counts_path.empty();
    if (!validate_identity_and_allocation_request(
            n_samples, nullptr, &main_bytes_per_cell, &validation_error)) {
        fprintf(stderr, "ERROR: invalid fused main identity universe: %s\n",
            validation_error.c_str());
        return false;
    }
    if (!validate_identity_and_allocation_request(
            species_native_n_samples, nullptr, &species_bytes_per_cell,
            &validation_error)) {
        fprintf(stderr, "ERROR: invalid fused native-species universe: %s\n",
            validation_error.c_str());
        return false;
    }
    if (n_threads < 1 || htslib_threads < 1 ||
        (raw_counts_path.empty() != raw_species_counts_path.empty()) ||
        (!write_raw_outputs && filtered_barcodes.empty())) {
        fprintf(stderr,
            "ERROR: fused counting requires positive worker/HTSlib thread counts and either both raw output paths or a filtered barcode set\n");
        return false;
    }
    if (dump_pileup && (pileup_prefix.empty() || main_snpdat.empty())) {
        fprintf(stderr, "ERROR: fused pileup output requires a non-empty prefix\n");
        return false;
    }
    if (!filtered_main_counts.empty() || !filtered_species_counts.empty()) {
        fprintf(stderr, "ERROR: fused count destinations must be empty\n");
        return false;
    }
    (void)mapq_probability_scaled_table();

    const FusedPerfClock::time_point dense_allocation_start = FusedPerfClock::now();
    if (!filtered_barcodes.empty()) {
        const size_t active_main_bytes = main_snpdat.empty() ? 0 : main_bytes_per_cell;
        const size_t combined_bytes = active_main_bytes + species_bytes_per_cell;
        if (combined_bytes < active_main_bytes ||
            (combined_bytes > 0 && filtered_barcodes.size() >
                std::numeric_limits<size_t>::max() / combined_bytes)) {
            fprintf(stderr, "ERROR: fused filtered-cell allocation overflows size_t\n");
            return false;
        }
        fprintf(stderr,
            "Pre-allocating fused dense main/native-species stores for %lu STARsolo-filtered cells...\n",
            (unsigned long)filtered_barcodes.size());
        try {
            for (unsigned long barcode : filtered_barcodes) {
                if (!main_snpdat.empty()) {
                    filtered_main_counts.emplace(
                        std::piecewise_construct,
                        std::forward_as_tuple(barcode),
                        std::forward_as_tuple(n_samples));
                }
                filtered_species_counts.emplace(
                    std::piecewise_construct,
                    std::forward_as_tuple(barcode),
                    std::forward_as_tuple(species_native_n_samples));
            }
        } catch (const std::exception& error) {
            fprintf(stderr, "ERROR: fused dense pre-allocation failed: %s\n",
                error.what());
            return false;
        }
    }
    print_fused_perf_phase(
        "filtered_dense_store_allocation", dense_allocation_start,
        FusedPerfClock::now());

    const FusedPerfClock::time_point bam_inspection_start = FusedPerfClock::now();
    htsFile* bam_header_file = hts_open(bamfile.c_str(), "r");
    if (!bam_header_file) {
        fprintf(stderr, "ERROR: could not open BAM for fused counting: %s\n",
            bamfile.c_str());
        return false;
    }
    bam_hdr_t* bam_header = sam_hdr_read(bam_header_file);
    hts_idx_t* bam_index =
        bam_header ? sam_index_load(bam_header_file, bamfile.c_str()) : nullptr;
    if (!bam_header || !bam_index) {
        fprintf(stderr, "ERROR: could not read BAM header/index for fused counting\n");
        if (bam_index) hts_idx_destroy(bam_index);
        if (bam_header) bam_hdr_destroy(bam_header);
        hts_close(bam_header_file);
        return false;
    }

    const int n_chromosomes = bam_header->n_targets;
    vector<int64_t> chromosome_lengths((size_t)n_chromosomes, 0);
    vector<uint64_t> chromosome_records((size_t)n_chromosomes, 0);
    vector<string> chromosome_names((size_t)n_chromosomes);
    const int n_index_targets = hts_idx_nseq(bam_index);
    set<int> panel_tids;
    for (const auto& chromosome : main_snpdat) panel_tids.insert(chromosome.first);
    for (const auto& chromosome : species_snpdat) panel_tids.insert(chromosome.first);
    for (int tid : panel_tids) {
        if (tid < 0 || tid >= n_chromosomes) {
            fprintf(stderr, "ERROR: fused panel references invalid BAM target id %d\n", tid);
            hts_idx_destroy(bam_index);
            bam_hdr_destroy(bam_header);
            hts_close(bam_header_file);
            return false;
        }
    }
    for (int tid = 0; tid < n_chromosomes; ++tid) {
        chromosome_lengths[(size_t)tid] = bam_header->target_len[tid];
        chromosome_names[(size_t)tid] = bam_header->target_name[tid]
            ? bam_header->target_name[tid] : std::to_string(tid);
        uint64_t mapped = 0;
        uint64_t unmapped = 0;
        if (tid < n_index_targets &&
            hts_idx_get_stat(bam_index, tid, &mapped, &unmapped) >= 0) {
            chromosome_records[(size_t)tid] = mapped;
        }
    }
    print_fused_perf_phase(
        "bam_header_index_inspection", bam_inspection_start,
        FusedPerfClock::now());

    if (dump_pileup) {
        const string site_path = pileup_prefix + ".pileup_sites.tsv.gz";
        if (!write_fused_pileup_sites(
                site_path, main_snpdat, chromosome_names, n_samples)) {
            fprintf(stderr, "ERROR: failed writing fused pileup site table: %s\n",
                site_path.c_str());
            hts_idx_destroy(bam_index);
            bam_hdr_destroy(bam_header);
            hts_close(bam_header_file);
            return false;
        }
    }

    const FusedPerfClock::time_point work_unit_start = FusedPerfClock::now();
    const size_t SNP_CHUNK_THRESHOLD = 100000;
    const uint64_t READ_CHUNK_THRESHOLD = 10000000;
    vector<FusedWorkUnit> work_units;
    for (int tid : panel_tids) {
        const ChromSNPs* main_panel = find_chrom_panel(main_snpdat, tid);
        const ChromSNPs* species_panel = find_chrom_panel(species_snpdat, tid);
        const size_t main_sites = main_panel ? main_panel->snps.size() : 0;
        const size_t species_sites = species_panel ? species_panel->snps.size() : 0;
        const size_t site_views = main_sites + species_sites;
        if (site_views == 0) continue;

        const size_t chunks_by_snp = site_views > SNP_CHUNK_THRESHOLD
            ? (site_views + SNP_CHUNK_THRESHOLD - 1) / SNP_CHUNK_THRESHOLD : 1;
        size_t chunks_by_reads =
            chromosome_records[(size_t)tid] > READ_CHUNK_THRESHOLD
                ? (size_t)((chromosome_records[(size_t)tid] +
                    READ_CHUNK_THRESHOLD - 1) / READ_CHUNK_THRESHOLD) : 1;
        chunks_by_reads = std::min(chunks_by_reads, (size_t)20);
        const size_t requested_chunks = std::max(chunks_by_snp, chunks_by_reads);
        vector<int> boundaries = chunks_by_snp >= chunks_by_reads
            ? fused_rank_boundaries(
                main_panel, species_panel, requested_chunks,
                chromosome_lengths[(size_t)tid])
            : fused_position_boundaries(
                requested_chunks, chromosome_lengths[(size_t)tid]);

        int start = 0;
        const uint64_t estimate = chromosome_records[(size_t)tid] /
            (uint64_t)std::max((size_t)1, boundaries.size() + 1);
        for (int boundary : boundaries) {
            if (boundary <= start) continue;
            FusedWorkUnit unit;
            unit.tid = tid;
            unit.owner_start = start;
            unit.owner_end = boundary;
            unit.estimated_records = estimate;
            work_units.push_back(unit);
            start = boundary;
        }
        FusedWorkUnit final_unit;
        final_unit.tid = tid;
        final_unit.owner_start = start;
        final_unit.owner_end = INT_MAX;
        final_unit.estimated_records = estimate;
        work_units.push_back(final_unit);
    }
    std::sort(work_units.begin(), work_units.end(),
        [](const FusedWorkUnit& left, const FusedWorkUnit& right) {
            return left.estimated_records > right.estimated_records;
        });
    print_fused_perf_phase(
        "work_unit_construction", work_unit_start, FusedPerfClock::now());

    hts_idx_destroy(bam_index);
    bam_hdr_destroy(bam_header);
    hts_close(bam_header_file);
    if (work_units.empty()) {
        fprintf(stderr, "ERROR: fused main/species panels contain no schedulable sites\n");
        return false;
    }
    fprintf(stderr,
        "Starting fused BAM traversal with non-overlapping start ownership and %d readers/workers\n",
        n_threads);

    const FusedPerfClock::time_point sparse_initialization_start =
        FusedPerfClock::now();
    FusedTemporaryFiles raw_temporary_files;
    const string raw_temporary_base = write_raw_outputs
        ? fused_temporary_base(raw_counts_path) : string();
    vector<std::unique_ptr<FusedRawPartitionFile>> raw_partitions;
    if (write_raw_outputs) {
        string error;
        if (!open_raw_partitions(
                raw_temporary_base, raw_partitions,
                raw_temporary_files, error)) {
            fprintf(stderr, "ERROR: %s\n", error.c_str());
            return false;
        }
    }

    std::unordered_map<unsigned long, AlignedCellCounts*> filtered_main_lookup;
    std::unordered_map<unsigned long, AlignedCellCounts*> filtered_species_lookup;
    filtered_main_lookup.reserve(filtered_barcodes.size());
    filtered_species_lookup.reserve(filtered_barcodes.size());
    for (unsigned long barcode : filtered_barcodes) {
        auto species_found = filtered_species_counts.find(barcode);
        if (species_found == filtered_species_counts.end()) {
            fprintf(stderr, "ERROR: fused filtered-cell lookup initialization failed\n");
            return false;
        }
        if (!main_snpdat.empty()) {
            auto main_found = filtered_main_counts.find(barcode);
            if (main_found == filtered_main_counts.end()) {
                fprintf(stderr, "ERROR: fused filtered-main lookup initialization failed\n");
                return false;
            }
            filtered_main_lookup.emplace(barcode, &main_found->second);
        }
        filtered_species_lookup.emplace(barcode, &species_found->second);
    }
    print_fused_perf_phase(
        "sparse_shard_lookup_initialization", sparse_initialization_start,
        FusedPerfClock::now());

    vector<string> spool_paths;
    vector<FILE*> spool_files;
    if (dump_pileup &&
        !open_fused_spools(pileup_prefix, n_threads, spool_paths, spool_files)) {
        fprintf(stderr, "ERROR: could not create fused pileup spool files\n");
        return false;
    }

    vector<AcceptedSiteWeightMap> thread_main_site_weights((size_t)n_threads);
    vector<AcceptedSiteWeightMap> thread_species_site_weights((size_t)n_threads);
    vector<FusedThreadPerfCounters> thread_perf_counters((size_t)n_threads);
    vector<FusedWorkUnitProfile> work_unit_profiles(work_units.size());
    ParallelOperationStatus operation_status;
    std::atomic<bool> hts_thread_warning_emitted(false);
    const FusedPerfClock::time_point traversal_start = FusedPerfClock::now();
    omp_set_num_threads(n_threads);

    #pragma omp parallel
    {
        const int thread_id = omp_get_thread_num();
        vector<FusedAlignedObservation> main_observations;
        vector<FusedAlignedObservation> species_observations;
        main_observations.reserve(32);
        species_observations.reserve(32);
        vector<FusedCountObservation> filtered_observation_buffer;
        vector<FusedCountObservation> raw_observation_buffer;
        vector<FusedCountObservation> raw_partition_order_buffer;
        filtered_observation_buffer.reserve(FUSED_OBSERVATION_BUFFER_RECORDS);
        raw_observation_buffer.reserve(FUSED_OBSERVATION_BUFFER_RECORDS);
        raw_partition_order_buffer.reserve(FUSED_OBSERVATION_BUFFER_RECORDS);
        AcceptedSiteWeightMap& local_main_weights =
            thread_main_site_weights[(size_t)thread_id];
        AcceptedSiteWeightMap& local_species_weights =
            thread_species_site_weights[(size_t)thread_id];
        FusedThreadPerfCounters local_perf;

        htsFile* bam_file = hts_open(bamfile.c_str(), "r");
        bam_hdr_t* header = nullptr;
        hts_idx_t* index = nullptr;
        bam1_t* record = nullptr;
        if (!bam_file) {
            operation_status.fail(format_worker_error("BAM open", thread_id));
        } else {
            if (htslib_threads > 1 &&
                hts_set_threads(bam_file, htslib_threads) < 0) {
                bool expected = false;
                if (hts_thread_warning_emitted.compare_exchange_strong(
                        expected, true)) {
                    fprintf(stderr,
                        "WARNING: HTSlib helper-thread setup failed; continuing with synchronous BAM I/O\n");
                }
            }
            header = sam_hdr_read(bam_file);
            if (!header)
                operation_status.fail(format_worker_error("BAM header read", thread_id));
            index = sam_index_load(bam_file, bamfile.c_str());
            if (!index)
                operation_status.fail(format_worker_error("BAM index load", thread_id));
            record = bam_init1();
            if (!record)
                operation_status.fail(format_worker_error("BAM record allocation", thread_id));
        }

        #pragma omp for schedule(dynamic, 1)
        for (size_t unit_index = 0; unit_index < work_units.size(); ++unit_index) {
            if (!operation_status.ok() || !bam_file || !header || !index || !record)
                continue;
            const FusedWorkUnit& unit = work_units[unit_index];
            FusedWorkUnitProfile unit_perf;
            const FusedPerfClock::time_point unit_start = FusedPerfClock::now();
            if (unit.tid < 0 || unit.tid >= header->n_targets) {
                operation_status.fail(
                    format_worker_error("invalid contig", thread_id, unit.tid));
                continue;
            }
            hts_itr_t* iterator = sam_itr_queryi(
                index, unit.tid, unit.owner_start, unit.owner_end);
            if (!iterator) {
                operation_status.fail(
                    format_worker_error("iterator creation", thread_id, unit.tid));
                continue;
            }

            const ChromSNPs* main_panel = find_chrom_panel(main_snpdat, unit.tid);
            const ChromSNPs* species_panel =
                find_chrom_panel(species_snpdat, unit.tid);
            size_t main_start_cursor = 0;
            size_t species_start_cursor = 0;
            if (main_panel) {
                main_start_cursor = (size_t)std::distance(
                    main_panel->snps.begin(),
                    std::lower_bound(
                        main_panel->snps.begin(), main_panel->snps.end(),
                        unit.owner_start,
                        [](const SNPData& snp, int position) {
                            return snp.pos < position;
                        }));
            }
            if (species_panel) {
                species_start_cursor = (size_t)std::distance(
                    species_panel->snps.begin(),
                    std::lower_bound(
                        species_panel->snps.begin(), species_panel->snps.end(),
                        unit.owner_start,
                        [](const SNPData& snp, int position) {
                            return snp.pos < position;
                        }));
            }

            int iterator_result = 0;
            while (operation_status.ok() &&
                    (iterator_result = sam_itr_next(
                        bam_file, iterator, record)) >= 0) {
                ++local_perf.iterator_records;
                ++unit_perf.iterator_records;
                if (!read_passes_filter(
                        record, default_production_read_filter())) continue;
                ++local_perf.records_passing_read_policy;
                ++unit_perf.records_passing_read_policy;

                // Indexed region iterators may return a long alignment on both
                // sides of a boundary. Alignment start uniquely selects its owner.
                if (record->core.pos < unit.owner_start ||
                    (unit.owner_end != INT_MAX &&
                     record->core.pos >= unit.owner_end)) continue;

                uint8_t* barcode_tag = bam_aux_get(record, "CB");
                if (!barcode_tag) continue;
                ++local_perf.cb_tagged_records;
                ++unit_perf.cb_tagged_records;
                const char* barcode_string = bam_aux2Z(barcode_tag);
                if (!barcode_string) continue;
                bc barcode_bits;
                str2bc(barcode_string, barcode_bits);
                const unsigned long barcode = barcode_bits.to_ulong();
                auto filtered_species_found =
                    filtered_species_lookup.find(barcode);
                const bool is_filtered =
                    filtered_species_found != filtered_species_lookup.end();
                if (is_filtered) {
                    ++local_perf.filtered_cell_records;
                    ++unit_perf.filtered_cell_records;
                } else {
                    ++local_perf.raw_only_records;
                    ++unit_perf.raw_only_records;
                }

                if (main_panel) {
                    while (main_start_cursor < main_panel->snps.size() &&
                            main_panel->snps[main_start_cursor].pos <
                                record->core.pos) {
                        ++main_start_cursor;
                    }
                }
                if (species_panel) {
                    while (species_start_cursor < species_panel->snps.size() &&
                            species_panel->snps[species_start_cursor].pos <
                                record->core.pos) {
                        ++species_start_cursor;
                    }
                }
                const int64_t probability_scaled =
                    mapq_probability_scaled(record->core.qual);
                collect_fused_alignment_observations(
                    record, main_panel, species_panel,
                    main_start_cursor, species_start_cursor,
                    probability_scaled,
                    main_observations, species_observations);
                local_perf.main_allele_observations += main_observations.size();
                local_perf.species_allele_observations += species_observations.size();
                unit_perf.main_allele_observations += main_observations.size();
                unit_perf.species_allele_observations += species_observations.size();
                if (main_observations.empty() && species_observations.empty()) continue;

                if (is_filtered) {
                    for (const FusedAlignedObservation& observation :
                            main_observations) {
                        if (accepted_site_weights_main) {
                            local_main_weights[accepted_site_weight_key(
                                unit.tid, observation.snp->pos)] +=
                                observation.ref_scaled + observation.alt_scaled;
                        }
                        if (!append_compact_observation(
                                filtered_observation_buffer, barcode, unit.tid,
                                FUSED_MAIN_PANEL, observation)) {
                            operation_status.fail(format_worker_error(
                                "filtered observation encoding", thread_id,
                                unit.tid));
                            break;
                        }
                    }
                    for (const FusedAlignedObservation& observation :
                            species_observations) {
                        if (accepted_site_weights_species) {
                            local_species_weights[accepted_site_weight_key(
                                unit.tid, observation.snp->pos)] +=
                                observation.ref_scaled + observation.alt_scaled;
                        }
                        if (!append_compact_observation(
                                filtered_observation_buffer, barcode, unit.tid,
                                FUSED_SPECIES_PANEL, observation)) {
                            operation_status.fail(format_worker_error(
                                "filtered species observation encoding",
                                thread_id, unit.tid));
                            break;
                        }
                    }

                    if (operation_status.ok() &&
                        filtered_observation_buffer.size() >=
                            FUSED_OBSERVATION_BUFFER_RECORDS) {
                        string error;
                        if (!flush_filtered_observations(
                                filtered_observation_buffer, main_snpdat,
                                species_native_targets, filtered_main_lookup,
                                filtered_species_lookup, local_perf, error)) {
                            operation_status.fail(format_worker_error(
                                error.c_str(), thread_id, unit.tid));
                        }
                    }

                    if (dump_pileup && !main_observations.empty()) {
                        const std::pair<uint64_t, uint8_t> molecule =
                            pileup_molecule_key(record);
                        FILE* spool = spool_files[(size_t)thread_id];
                        for (const FusedAlignedObservation& observation :
                                main_observations) {
                            FusedPileupSpoolRecord spool_record;
                            spool_record.barcode = (uint64_t)barcode;
                            spool_record.molecule_hash = molecule.first;
                            spool_record.site =
                                ((uint64_t)(uint32_t)unit.tid << 32) |
                                (uint64_t)(uint32_t)observation.snp->pos;
                            spool_record.ref_scaled = observation.ref_scaled;
                            spool_record.alt_scaled = observation.alt_scaled;
                            spool_record.basis = molecule.second;
                            if (fwrite(&spool_record, sizeof(spool_record), 1,
                                    spool) != 1) {
                                operation_status.fail(format_worker_error(
                                    "pileup spool write", thread_id, unit.tid));
                                break;
                            }
                        }
                    }
                } else {
                    if (!write_raw_outputs) continue;
                    for (const FusedAlignedObservation& observation :
                            main_observations) {
                        if (!append_compact_observation(
                                raw_observation_buffer, barcode, unit.tid,
                                FUSED_MAIN_PANEL, observation)) {
                            operation_status.fail(format_worker_error(
                                "raw main observation encoding", thread_id,
                                unit.tid));
                            break;
                        }
                    }
                    for (const FusedAlignedObservation& observation :
                            species_observations) {
                        if (!append_compact_observation(
                                raw_observation_buffer, barcode, unit.tid,
                                FUSED_SPECIES_PANEL, observation)) {
                            operation_status.fail(format_worker_error(
                                "raw species observation encoding", thread_id,
                                unit.tid));
                            break;
                        }
                    }
                    if (operation_status.ok() &&
                        raw_observation_buffer.size() >=
                            FUSED_OBSERVATION_BUFFER_RECORDS) {
                        string error;
                        if (!spill_raw_observations(
                                raw_observation_buffer,
                                raw_partition_order_buffer, raw_partitions,
                                local_perf, error)) {
                            operation_status.fail(format_worker_error(
                                error.c_str(), thread_id, unit.tid));
                        }
                    }
                }
            }
            if (iterator_result < -1) {
                operation_status.fail(
                    format_worker_error("iterator read", thread_id, unit.tid));
            }
            hts_itr_destroy(iterator);
            unit_perf.elapsed_nanoseconds = fused_perf_nanoseconds(
                unit_start, FusedPerfClock::now());
            work_unit_profiles[unit_index] = unit_perf;
        }

        if (operation_status.ok() && !filtered_observation_buffer.empty()) {
            string error;
            if (!flush_filtered_observations(
                    filtered_observation_buffer, main_snpdat,
                    species_native_targets, filtered_main_lookup,
                    filtered_species_lookup, local_perf, error)) {
                operation_status.fail(format_worker_error(
                    error.c_str(), thread_id));
            }
        }
        if (operation_status.ok() && !raw_observation_buffer.empty()) {
            string error;
            if (!spill_raw_observations(
                    raw_observation_buffer, raw_partition_order_buffer,
                    raw_partitions,
                    local_perf, error)) {
                operation_status.fail(format_worker_error(
                    error.c_str(), thread_id));
            }
        }

        if (record) bam_destroy1(record);
        if (index) hts_idx_destroy(index);
        if (header) bam_hdr_destroy(header);
        if (bam_file) hts_close(bam_file);
        thread_perf_counters[(size_t)thread_id] = local_perf;
    }
    print_fused_perf_phase(
        "parallel_bam_traversal", traversal_start, FusedPerfClock::now());

    FusedThreadPerfCounters total_perf;
    for (const FusedThreadPerfCounters& thread_perf : thread_perf_counters) {
        total_perf.iterator_records += thread_perf.iterator_records;
        total_perf.records_passing_read_policy +=
            thread_perf.records_passing_read_policy;
        total_perf.cb_tagged_records += thread_perf.cb_tagged_records;
        total_perf.filtered_cell_records += thread_perf.filtered_cell_records;
        total_perf.raw_only_records += thread_perf.raw_only_records;
        total_perf.main_allele_observations +=
            thread_perf.main_allele_observations;
        total_perf.species_allele_observations +=
            thread_perf.species_allele_observations;
        total_perf.main_total_target_updates +=
            thread_perf.main_total_target_updates;
        total_perf.main_pair_target_updates +=
            thread_perf.main_pair_target_updates;
        total_perf.species_total_target_updates +=
            thread_perf.species_total_target_updates;
        total_perf.species_pair_target_updates +=
            thread_perf.species_pair_target_updates;
        total_perf.raw_shard_acquisitions +=
            thread_perf.raw_shard_acquisitions;
        total_perf.raw_partition_bulk_writes +=
            thread_perf.raw_partition_bulk_writes;
        total_perf.filtered_cell_lock.applicable_records +=
            thread_perf.filtered_cell_lock.applicable_records;
        total_perf.filtered_cell_lock.samples +=
            thread_perf.filtered_cell_lock.samples;
        total_perf.filtered_cell_lock.wait_nanoseconds +=
            thread_perf.filtered_cell_lock.wait_nanoseconds;
        total_perf.filtered_cell_lock.held_nanoseconds +=
            thread_perf.filtered_cell_lock.held_nanoseconds;
        total_perf.raw_shard_lock.applicable_records +=
            thread_perf.raw_shard_lock.applicable_records;
        total_perf.raw_shard_lock.samples +=
            thread_perf.raw_shard_lock.samples;
        total_perf.raw_shard_lock.wait_nanoseconds +=
            thread_perf.raw_shard_lock.wait_nanoseconds;
        total_perf.raw_shard_lock.held_nanoseconds +=
            thread_perf.raw_shard_lock.held_nanoseconds;
        total_perf.raw_partition_lock.applicable_records +=
            thread_perf.raw_partition_lock.applicable_records;
        total_perf.raw_partition_lock.samples +=
            thread_perf.raw_partition_lock.samples;
        total_perf.raw_partition_lock.wait_nanoseconds +=
            thread_perf.raw_partition_lock.wait_nanoseconds;
        total_perf.raw_partition_lock.held_nanoseconds +=
            thread_perf.raw_partition_lock.held_nanoseconds;
    }

    uint64_t unique_raw_only_barcodes = 0;
    uint64_t raw_temporary_bytes = 0;
    if (write_raw_outputs) {
        string error;
        if (!close_raw_partitions(raw_partitions, error)) {
            operation_status.fail(error);
        }
        for (const auto& partition : raw_partitions) {
            raw_temporary_bytes += partition->bytes;
        }
    }

    if (dump_pileup && !close_fused_spools(spool_files)) {
        operation_status.fail("failed closing fused pileup spool files");
    }
    if (!operation_status.ok()) {
        fprintf(stderr, "\nERROR: fused parallel counting failed: %s\n",
            operation_status.message().c_str());
        remove_files(spool_paths);
        return false;
    }
    fprintf(stderr,
        "PERF_COUNTER name=raw_observation_partitions partitions=%llu temporary_bytes=%llu bulk_writes=%llu\n",
        (unsigned long long)(write_raw_outputs ? FUSED_RAW_PARTITIONS : 0),
        (unsigned long long)raw_temporary_bytes,
        (unsigned long long)total_perf.raw_partition_bulk_writes);

    vector<uint64_t> work_unit_elapsed;
    work_unit_elapsed.reserve(work_unit_profiles.size());
    vector<size_t> slowest_work_units;
    slowest_work_units.reserve(work_unit_profiles.size());
    for (size_t unit_index = 0; unit_index < work_unit_profiles.size(); ++unit_index) {
        work_unit_elapsed.push_back(
            work_unit_profiles[unit_index].elapsed_nanoseconds);
        slowest_work_units.push_back(unit_index);
    }
    std::sort(work_unit_elapsed.begin(), work_unit_elapsed.end());
    std::sort(slowest_work_units.begin(), slowest_work_units.end(),
        [&work_unit_profiles](size_t left, size_t right) {
            if (work_unit_profiles[left].elapsed_nanoseconds !=
                    work_unit_profiles[right].elapsed_nanoseconds) {
                return work_unit_profiles[left].elapsed_nanoseconds >
                    work_unit_profiles[right].elapsed_nanoseconds;
            }
            return left < right;
        });
    long double median_nanoseconds = 0.0L;
    if (!work_unit_elapsed.empty()) {
        const size_t middle = work_unit_elapsed.size() / 2;
        if (work_unit_elapsed.size() % 2 == 0) {
            median_nanoseconds =
                ((long double)work_unit_elapsed[middle - 1] +
                 (long double)work_unit_elapsed[middle]) / 2.0L;
        } else {
            median_nanoseconds = (long double)work_unit_elapsed[middle];
        }
    }
    const uint64_t maximum_nanoseconds = work_unit_elapsed.empty()
        ? 0 : work_unit_elapsed.back();
    fprintf(stderr,
        "PERF_WORKUNIT summary total=%llu median_seconds=%.6f max_seconds=%.6f\n",
        (unsigned long long)work_unit_profiles.size(),
        (double)(median_nanoseconds / 1000000000.0L),
        (double)maximum_nanoseconds / 1000000000.0);
    const size_t slowest_count = std::min((size_t)20, slowest_work_units.size());
    for (size_t rank = 0; rank < slowest_count; ++rank) {
        const size_t unit_index = slowest_work_units[rank];
        const FusedWorkUnit& unit = work_units[unit_index];
        const FusedWorkUnitProfile& profile = work_unit_profiles[unit_index];
        fprintf(stderr,
            "PERF_WORKUNIT rank=%llu contig=%s tid=%d start=%d end=%d elapsed_seconds=%.6f iterator_records=%llu passing_read_policy=%llu cb_tagged_records=%llu filtered_records=%llu raw_only_records=%llu main_observations=%llu species_observations=%llu\n",
            (unsigned long long)(rank + 1),
            chromosome_names[(size_t)unit.tid].c_str(), unit.tid,
            unit.owner_start, unit.owner_end,
            (double)profile.elapsed_nanoseconds / 1000000000.0,
            (unsigned long long)profile.iterator_records,
            (unsigned long long)profile.records_passing_read_policy,
            (unsigned long long)profile.cb_tagged_records,
            (unsigned long long)profile.filtered_cell_records,
            (unsigned long long)profile.raw_only_records,
            (unsigned long long)profile.main_allele_observations,
            (unsigned long long)profile.species_allele_observations);
    }
    fprintf(stderr,
        "PERF_LOCK_SAMPLE store=filtered_cell interval=%llu applicable_records=%llu samples=%llu wait_seconds=%.6f held_seconds=%.6f\n",
        (unsigned long long)FUSED_LOCK_SAMPLE_INTERVAL,
        (unsigned long long)total_perf.filtered_cell_lock.applicable_records,
        (unsigned long long)total_perf.filtered_cell_lock.samples,
        (double)total_perf.filtered_cell_lock.wait_nanoseconds / 1000000000.0,
        (double)total_perf.filtered_cell_lock.held_nanoseconds / 1000000000.0);
    fprintf(stderr,
        "PERF_LOCK_SAMPLE store=raw_shard interval=%llu applicable_records=%llu samples=%llu wait_seconds=%.6f held_seconds=%.6f\n",
        (unsigned long long)FUSED_LOCK_SAMPLE_INTERVAL,
        (unsigned long long)total_perf.raw_shard_lock.applicable_records,
        (unsigned long long)total_perf.raw_shard_lock.samples,
        (double)total_perf.raw_shard_lock.wait_nanoseconds / 1000000000.0,
        (double)total_perf.raw_shard_lock.held_nanoseconds / 1000000000.0);
    fprintf(stderr,
        "PERF_LOCK_SAMPLE store=raw_partition interval=%llu applicable_records=%llu samples=%llu wait_seconds=%.6f held_seconds=%.6f\n",
        (unsigned long long)FUSED_LOCK_SAMPLE_INTERVAL,
        (unsigned long long)total_perf.raw_partition_lock.applicable_records,
        (unsigned long long)total_perf.raw_partition_lock.samples,
        (double)total_perf.raw_partition_lock.wait_nanoseconds / 1000000000.0,
        (double)total_perf.raw_partition_lock.held_nanoseconds / 1000000000.0);

    auto merge_site_weights = [n_threads](
            vector<AcceptedSiteWeightMap>& per_thread,
            AcceptedSiteWeightMap* destination) {
        if (!destination) return;
        destination->clear();
        for (int thread = 0; thread < n_threads; ++thread) {
            for (const auto& item : per_thread[(size_t)thread]) {
                (*destination)[item.first] += item.second;
            }
        }
    };
    const FusedPerfClock::time_point accepted_merge_start = FusedPerfClock::now();
    merge_site_weights(thread_main_site_weights, accepted_site_weights_main);
    merge_site_weights(thread_species_site_weights, accepted_site_weights_species);
    print_fused_perf_phase(
        "accepted_site_weight_merging", accepted_merge_start,
        FusedPerfClock::now());

    if (write_raw_outputs) {
        FusedThreadPerfCounters reduction_perf;
        if (!write_fused_raw_counts(
                raw_counts_path, raw_species_counts_path,
                filtered_main_counts, filtered_species_counts, raw_partitions,
                raw_temporary_files, raw_temporary_base, main_snpdat,
                species_native_targets, n_samples, species_native_n_samples,
                n_threads, reduction_perf, unique_raw_only_barcodes)) {
            fprintf(stderr, "\nERROR: failed serializing fused raw count products\n");
            remove_files(spool_paths);
            return false;
        }
        total_perf.main_total_target_updates +=
            reduction_perf.main_total_target_updates;
        total_perf.main_pair_target_updates +=
            reduction_perf.main_pair_target_updates;
        total_perf.species_total_target_updates +=
            reduction_perf.species_total_target_updates;
        total_perf.species_pair_target_updates +=
            reduction_perf.species_pair_target_updates;
    } else {
        fprintf(stderr,
            "PERF_PHASE name=filtered_main_raw_serialization seconds=0.000000 skipped=1\n");
        fprintf(stderr,
            "PERF_PHASE name=filtered_species_raw_serialization seconds=0.000000 skipped=1\n");
        fprintf(stderr,
            "PERF_PHASE name=sparse_raw_only_main_species_serialization seconds=0.000000 skipped=1\n");
        fprintf(stderr,
            "PERF_PHASE name=gzip_flush_close seconds=0.000000 skipped=1\n");
    }

    if (dump_pileup) {
        const string observation_path = pileup_prefix + ".pileup_obs.tsv.gz";
        const string molecule_path = pileup_prefix + ".pileup_molecules.tsv.gz";
        if (!derive_fused_pileup_outputs(
                spool_paths, observation_path, molecule_path)) {
            fprintf(stderr, "ERROR: failed deriving public pileup views from spool\n");
            return false;
        }
    }

    fprintf(stderr,
        "PERF_COUNTER name=fused_summary iterator_records=%llu passing_read_policy=%llu cb_tagged_records=%llu filtered_cell_records=%llu raw_only_records=%llu main_allele_observations=%llu species_allele_observations=%llu main_total_target_updates=%llu main_pair_target_updates=%llu species_total_target_updates=%llu species_pair_target_updates=%llu raw_shard_acquisitions=%llu unique_raw_only_barcodes=%llu\n",
        (unsigned long long)total_perf.iterator_records,
        (unsigned long long)total_perf.records_passing_read_policy,
        (unsigned long long)total_perf.cb_tagged_records,
        (unsigned long long)total_perf.filtered_cell_records,
        (unsigned long long)total_perf.raw_only_records,
        (unsigned long long)total_perf.main_allele_observations,
        (unsigned long long)total_perf.species_allele_observations,
        (unsigned long long)total_perf.main_total_target_updates,
        (unsigned long long)total_perf.main_pair_target_updates,
        (unsigned long long)total_perf.species_total_target_updates,
        (unsigned long long)total_perf.species_pair_target_updates,
        (unsigned long long)total_perf.raw_shard_acquisitions,
        (unsigned long long)unique_raw_only_barcodes);
    fprintf(stderr, "Fused counting complete\n");
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
    
    // Filtered-barcode runs update one preallocated matrix per cell under the
    // existing per-cell mutex. Raw-barcode discovery uses independently locked
    // shards so unknown keys remain unique across workers.
    omp_set_num_threads(n_threads);
    vector<std::unique_ptr<RawCountShard>> raw_count_shards =
        make_raw_count_shards(!has_bc_list);

    std::unordered_map<unsigned long, AlignedCellCounts*> shared_count_lookup_p0;
    std::unordered_map<unsigned long, AlignedCellCounts*> shared_count_lookup_p1;
    std::unordered_map<unsigned long, AlignedCellCounts*> shared_species_native_lookup;
    if (has_bc_list) {
        shared_count_lookup_p0.reserve(valid_barcodes.size());
        shared_count_lookup_p1.reserve(valid_barcodes.size());
        if (collect_species_native) {
            shared_species_native_lookup.reserve(valid_barcodes.size());
        }
        for (unsigned long barcode : valid_barcodes) {
            auto p0_it = counts_panel0.find(barcode);
            auto p1_it = counts_panel1.find(barcode);
            if (p0_it == counts_panel0.end() || p1_it == counts_panel1.end()) {
                fprintf(stderr, "ERROR: internal dual-panel pre-allocation failure\n");
                return false;
            }
            shared_count_lookup_p0.emplace(barcode, &p0_it->second);
            shared_count_lookup_p1.emplace(barcode, &p1_it->second);
            if (collect_species_native) {
                auto native_it = species_native_counts->find(barcode);
                if (native_it == species_native_counts->end()) {
                    fprintf(stderr, "ERROR: internal native-species pre-allocation failure\n");
                    return false;
                }
                shared_species_native_lookup.emplace(barcode, &native_it->second);
            }
        }
        fprintf(stderr,
            "Filtered dual-panel accumulation uses one shared count matrix per "
            "cell/panel; per-thread dense matrices are disabled.\n");
    }
    else {
        fprintf(stderr,
            "Raw dual-panel accumulation uses %lu shared count shards; "
            "per-thread dense matrices are disabled.\n",
            (unsigned long)RAW_COUNT_SHARDS);
    }

    // --dump_pileup: per-thread per-(cell,SNP) allele evidence (interindividual
    // only).  Inner key packs (tid<<32 | pos); value is (ref_scaled, alt_scaled).
    // Empty and untouched unless dump_pileup is set.  Identical structure and
    // downstream contract to the single-panel path.
    vector<PileupObservationMap> thread_pileup(n_threads);
    vector<vector<PileupMoleculeObservation> > thread_pileup_molecules(n_threads);
    vector<size_t> thread_pileup_entries(n_threads, 0);
    vector<long> thread_pileup_rows(n_threads, 0);
    vector<long> thread_molecule_rows(n_threads, 0);
    vector<string> pileup_observation_part_paths;
    vector<string> pileup_molecule_part_paths;
    vector<gzFile> pileup_observation_part_files;
    vector<gzFile> pileup_molecule_part_files;
    string pileup_stream_error;
    if (dump_pileup && !open_parallel_pileup_parts(
            pileup_prefix, n_threads,
            pileup_observation_part_paths, pileup_observation_part_files,
            pileup_molecule_part_paths, pileup_molecule_part_files,
            pileup_stream_error)) {
        fprintf(stderr, "ERROR: %s\n", pileup_stream_error.c_str());
        return false;
    }
    if (dump_pileup) {
        fprintf(stderr,
            "Pileup memory mode: bounded worker chunks (%lu cell/site, %lu molecule/site).\n",
            (unsigned long)PILEUP_SITE_CHUNK_ENTRIES,
            (unsigned long)PILEUP_MOLECULE_CHUNK_ENTRIES);
    }
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

                        std::pair<uint64_t, uint8_t> molecule_key;
                        if (dump_pileup) {
                            molecule_key = pileup_molecule_key(record);
                        }

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
                        
                        int64_t prob_scaled =
                            mapq_probability_scaled(record->core.qual);
                        
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
                                // Precomputed targets: linear traversal, no branches
                                const auto& ttargets = snp_check->total_targets;
                                const auto& ptargets = snp_check->pair_targets;
                                auto add_precomputed = [&](CellCounts& counts) {
                                    const bool is_ref = ref_add != 0;
                                    for (const auto& t : ttargets){
                                        if (is_ref) counts.total_ref[t.total_idx] += ref_add;
                                        else counts.total_alt[t.total_idx] += alt_add;
                                    }
                                    for (const auto& p : ptargets){
                                        if (is_ref) counts.ref_counts[p.pair_idx] += ref_add;
                                        else counts.alt_counts[p.pair_idx] += alt_add;
                                    }
                                };

                                if (has_bc_list) {
                                    const auto& target_lookup = panel_id == 0 ?
                                        shared_count_lookup_p0 : shared_count_lookup_p1;
                                    auto shared_it = target_lookup.find(bc_key);
                                    if (shared_it == target_lookup.end()) {
                                        operation_status.fail(format_worker_error(
                                            "filtered dual-panel count lookup", thread_id, tid));
                                        continue;
                                    }
                                    AlignedCellCounts& shared = *shared_it->second;
                                    std::lock_guard<std::mutex> guard(shared.lock);
                                    add_precomputed(shared.counts);
                                }
                                else {
                                    RawCountShard& shard = *raw_count_shards[
                                        raw_count_shard_index(bc_key)];
                                    std::lock_guard<std::mutex> guard(shard.lock);
                                    auto* target_counts = panel_id == 0 ?
                                        &shard.panel0 : &shard.panel1;
                                    auto it = target_counts->find(bc_key);
                                    if (it == target_counts->end()){
                                        target_counts->emplace(bc_key, CellCounts(n_samples));
                                        it = target_counts->find(bc_key);
                                    }
                                    add_precomputed(it->second);
                                }

                                // --dump_pileup: record per-(cell,SNP) evidence for
                                // interindividual SNPs.  Summed within a thread; the
                                // producer sums any cross-thread duplicates.
                                if (dump_pileup && snp_check->panel_id == 0){
                                    int64_t pkey = ((int64_t)tid << 32) |
                                        (int64_t)(uint32_t)snp_check->pos;
                                    auto& cell_sites = thread_pileup[thread_id][bc_key];
                                    auto site_it = cell_sites.find(pkey);
                                    if (site_it == cell_sites.end()) {
                                        cell_sites.emplace(
                                            pkey, std::make_pair(ref_add, alt_add));
                                        ++thread_pileup_entries[thread_id];
                                    }
                                    else {
                                        site_it->second.first += ref_add;
                                        site_it->second.second += alt_add;
                                    }
                                    PileupMoleculeObservation observation;
                                    observation.barcode = bc_key;
                                    observation.molecule_hash = molecule_key.first;
                                    observation.site = pkey;
                                    observation.ref_scaled = ref_add;
                                    observation.alt_scaled = alt_add;
                                    observation.basis = molecule_key.second;
                                    thread_pileup_molecules[thread_id].push_back(
                                        observation);
                                    if (thread_pileup_entries[thread_id] >=
                                            PILEUP_SITE_CHUNK_ENTRIES) {
                                        const long written = write_collapsed_pileup_observations(
                                            pileup_observation_part_files[thread_id],
                                            thread_pileup[thread_id]);
                                        thread_pileup_entries[thread_id] = 0;
                                        if (written < 0) {
                                            operation_status.fail(format_worker_error(
                                                "pileup observation write", thread_id, tid));
                                        }
                                        else thread_pileup_rows[thread_id] += written;
                                    }
                                    if (thread_pileup_molecules[thread_id].size() >=
                                            PILEUP_MOLECULE_CHUNK_ENTRIES) {
                                        const long written = write_collapsed_pileup_molecules(
                                            pileup_molecule_part_files[thread_id],
                                            thread_pileup_molecules[thread_id]);
                                        if (written < 0) {
                                            operation_status.fail(format_worker_error(
                                                "pileup molecule write", thread_id, tid));
                                        }
                                        else thread_molecule_rows[thread_id] += written;
                                    }
                                }

                                if (collect_species_native && panel_id == 1 &&
                                    native_chrom_targets != nullptr){
                                    const size_t snp_index = (size_t)(
                                        snp_check - chrom_snps.snps.begin());
                                    if (snp_index < native_chrom_targets->site_offsets.size() &&
                                        native_chrom_targets->site_offsets[snp_index] != UINT64_MAX){
                                        if (has_bc_list) {
                                            auto native_it = shared_species_native_lookup.find(bc_key);
                                            if (native_it == shared_species_native_lookup.end()) {
                                                operation_status.fail(format_worker_error(
                                                    "filtered native-species count lookup",
                                                    thread_id, tid));
                                                continue;
                                            }
                                            AlignedCellCounts& shared_native = *native_it->second;
                                            std::lock_guard<std::mutex> guard(shared_native.lock);
                                            accumulate_species_native_targets(
                                                shared_native.counts, *native_chrom_targets,
                                                snp_index, ref_add, alt_add);
                                        }
                                        else {
                                            RawCountShard& shard = *raw_count_shards[
                                                raw_count_shard_index(bc_key)];
                                            std::lock_guard<std::mutex> guard(shard.lock);
                                            auto native_it = shard.native.find(bc_key);
                                            if (native_it == shard.native.end()){
                                                shard.native.emplace(
                                                    bc_key, CellCounts(species_native_n_samples));
                                                native_it = shard.native.find(bc_key);
                                            }
                                            accumulate_species_native_targets(
                                                native_it->second, *native_chrom_targets, snp_index,
                                                ref_add, alt_add);
                                        }
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

        if (dump_pileup && operation_status.ok()) {
            const long obs_written = write_collapsed_pileup_observations(
                pileup_observation_part_files[thread_id], thread_pileup[thread_id]);
            thread_pileup_entries[thread_id] = 0;
            const long molecule_written = write_collapsed_pileup_molecules(
                pileup_molecule_part_files[thread_id],
                thread_pileup_molecules[thread_id]);
            if (obs_written < 0 || molecule_written < 0) {
                operation_status.fail(format_worker_error(
                    "final bounded pileup write", thread_id));
            }
            else {
                thread_pileup_rows[thread_id] += obs_written;
                thread_molecule_rows[thread_id] += molecule_written;
            }
        }
    }

    if (dump_pileup) {
        string close_error;
        const bool obs_closed = close_parallel_pileup_parts(
            pileup_observation_part_files, close_error);
        const bool molecule_closed = close_parallel_pileup_parts(
            pileup_molecule_part_files, close_error);
        if (!obs_closed || !molecule_closed) operation_status.fail(close_error);
    }

    if (!operation_status.ok()){
        fprintf(stderr, "ERROR: dual-panel allele counting failed: %s\n",
            operation_status.message().c_str());
        remove_files(pileup_observation_part_paths);
        remove_files(pileup_molecule_part_paths);
        return false;
    }
    
    if (!has_bc_list) {
        fprintf(stderr, "\nMoving sharded raw-barcode counts (panel 0)...\n");
        move_sharded_counts(raw_count_shards, 0, counts_panel0);
        fprintf(stderr, "Moving sharded raw-barcode counts (panel 1)...\n");
        move_sharded_counts(raw_count_shards, 1, counts_panel1);
    }
    
    if (collect_species_native && !has_bc_list){
        fprintf(stderr, "Moving sharded native species counts (panel 1)...\n");
        move_sharded_counts(raw_count_shards, 2, *species_native_counts);
        fprintf(stderr, "Native species counts (panel 1): %lu cells, %d species\n",
            species_native_counts->size(), species_native_n_samples);
    }

    raw_count_shards.clear();
    raw_count_shards.shrink_to_fit();

    // Publish bounded gzip members after the count completes successfully.
    if (dump_pileup){
        string obs_path = pileup_prefix + ".pileup_obs.tsv.gz";
        if (!publish_concatenated_gzip_members(
                pileup_observation_part_paths, obs_path, pileup_stream_error)) {
            fprintf(stderr, "ERROR: %s\n", pileup_stream_error.c_str());
            remove_files(pileup_observation_part_paths);
            remove_files(pileup_molecule_part_paths);
            return false;
        }

        string molecule_path = pileup_prefix + ".pileup_molecules.tsv.gz";
        if (!publish_concatenated_gzip_members(
                pileup_molecule_part_paths, molecule_path, pileup_stream_error)) {
            fprintf(stderr, "ERROR: %s\n", pileup_stream_error.c_str());
            remove_files(pileup_observation_part_paths);
            remove_files(pileup_molecule_part_paths);
            return false;
        }
        long n_obs_written = 0;
        long n_molecule_rows = 0;
        for (int t = 0; t < n_threads; ++t) {
            n_obs_written += thread_pileup_rows[t];
            n_molecule_rows += thread_molecule_rows[t];
        }
        remove_files(pileup_observation_part_paths);
        remove_files(pileup_molecule_part_paths);
        fprintf(stderr, "Wrote %ld bounded pileup observation rows to %s\n",
            n_obs_written, obs_path.c_str());
        fprintf(stderr, "Wrote %ld bounded molecule/SNP rows to %s\n",
            n_molecule_rows, molecule_path.c_str());
    }
    thread_pileup.clear();
    thread_pileup.shrink_to_fit();
    thread_pileup_molecules.clear();
    thread_pileup_molecules.shrink_to_fit();

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
    const size_t mapped_size = (size_t)sb.st_size;
    if (mapped_size < sizeof(SharedVCFHeader) || header->total_size > mapped_size ||
        header->n_chromosomes < 0 || header->n_chromosomes > 8192 ||
        header->n_samples < 0 || header->n_samples > 512) {
        fprintf(stderr, "ERROR: malformed shared VCF header\n");
        munmap(ptr, mapped_size);
        close(shm_fd);
        return false;
    }
    
    fprintf(stderr, "Attached to shared VCF: %s (%d SNPs, %d chromosomes)\n",
        shm_name.c_str(), header->n_snps_total, header->n_chromosomes);
    
    // Retrieve sample names
    samples.clear();
    for (int i = 0; i < header->n_samples && i < 512; i++){
        samples.push_back(string(header->sample_names[i]));
    }
    fprintf(stderr, "Loaded %lu sample names from shared VCF\n", samples.size());
    
    // Validate all source ranges before constructing any destination object.
    set<int> seen_tids;
    for (int i = 0; i < header->n_chromosomes; ++i) {
        const size_t offset = header->chrom_offsets[i];
        const size_t count = header->chrom_snp_counts[i];
        if (!seen_tids.insert(header->chrom_tids[i]).second ||
            offset > mapped_size ||
            count > (mapped_size - offset) / sizeof(SNPData)) {
            fprintf(stderr,
                "ERROR: malformed shared VCF chromosome range at index %d\n", i);
            munmap(ptr, mapped_size);
            close(shm_fd);
            return false;
        }
    }

    // Pre-create every map destination serially, size each vector exactly once,
    // then copy independent contigs in parallel without mutating the map.
    snpdat_all.clear();
    struct SharedCopyTask {
        const SNPData* source;
        SNPData* destination;
        size_t count;
    };
    vector<SharedCopyTask> copy_tasks;
    copy_tasks.reserve((size_t)header->n_chromosomes);
    try {
        snpdat_all.reserve((size_t)header->n_chromosomes * 2 + 1);
        for (int i = 0; i < header->n_chromosomes; ++i) {
            const int tid = header->chrom_tids[i];
            const size_t n_snps = header->chrom_snp_counts[i];
            snpdat_all.emplace(tid, ChromSNPs());
            ChromSNPs& chromosome = snpdat_all.find(tid)->second;
            chromosome.snps.resize(n_snps);
            SharedCopyTask task;
            task.source = reinterpret_cast<const SNPData*>(
                (const char*)ptr + header->chrom_offsets[i]);
            task.destination = chromosome.snps.data();
            task.count = n_snps;
            copy_tasks.push_back(task);
        }
    } catch (const std::exception& error) {
        fprintf(stderr, "ERROR: shared VCF destination allocation failed: %s\n",
            error.what());
        snpdat_all.clear();
        munmap(ptr, mapped_size);
        close(shm_fd);
        return false;
    }

    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t task_index = 0; task_index < copy_tasks.size(); ++task_index) {
        const SharedCopyTask& task = copy_tasks[task_index];
        for (size_t j = 0; j < task.count; ++j) {
            // Cannot use push_back(snp_ptr[j]) because the SNPData in shared
            // memory was written via memcpy, and the var struct contains a
            // std::vector<float> gqs whose internal heap pointer is from the
            // daemon's address space. Copying it would segfault.
            // Instead, construct each SNPData from the safe POD fields only.
            const SNPData& source = task.source[j];
            SNPData& destination = task.destination[j];
            destination.pos = source.pos;
            destination.panel_id = source.panel_id;
            destination.data.ref = source.data.ref;
            destination.data.alt = source.data.alt;
            destination.data.haps1 = source.data.haps1;
            destination.data.haps2 = source.data.haps2;
            destination.data.haps_covered = source.data.haps_covered;
            destination.data.vq = source.data.vq;
            // gqs is left empty (default-constructed); it's not needed
            // for counting or conditional match fraction computation.
        }
    }

    if (munmap(ptr, mapped_size) != 0) {
        perror("munmap shared VCF");
        snpdat_all.clear();
        close(shm_fd);
        return false;
    }
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

#endif  // CELLBOUNCER_VCF_HTS_INTERFACE_REVISION == 21901
