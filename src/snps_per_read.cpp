// ============================================================================
// snps_per_read
// ----------------------------------------------------------------------------
// Streams one BAM once and, for every aligned read that the alignment places
// over a VCF SNP locus, records how many such loci that read covers. Counts are
// aggregated per (library, CB) and per SNP class, then emitted as a per-cell
// long-format table, a rolled-up histogram, and summary stats. Optional raw
// per-read dump behind --dump_reads.
//
// One BAM = one library per invocation. Run once per library; concatenate the
// per-cell tables downstream and join (library, CB) to cell-type labels.
//
// Standalone: htslib + zlib + pthread + OpenMP only. No cellbouncer object, no
// htswrapper library link (robin_hood is header-only).
//
// Provenance: V1_R1
//
// Revision history:
//   V1_R1: Initial implementation of snps_per_read_SPEC_V1_R1. Lean
//          position-only VCF load with per-locus class mask; contiguous,
//          non-overlapping position partitions with read ownership by alignment
//          start; CIGAR coverage walk (M/=/X cover, D/N/I/S/H/P do not); per-cell
//          long-format output with independent class columns plus union ANY;
//          rolled-up histogram and summary; optional raw per-read dump.
// ============================================================================

#include <algorithm>
#include <atomic>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>
#include <omp.h>

#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htswrapper/robin_hood/robin_hood.h>

using namespace std;

// ----------------------------------------------------------------------------
// Read flag filter (spec section 2.2).
//
// A read is rejected if any of these bits are set:
//   UNMAP | SECONDARY | SUPPLEMENTARY | QCFAIL | DUP   == 0xF04
//
// This is the demux scan set (0x704 = UNMAP|SECONDARY|QCFAIL|DUP) plus
// SUPPLEMENTARY, so supplementary segments do not double-count coverage for a
// single read. It differs by one bit from the synth CB push-down mask
// (0xD04 = UNMAP|SECONDARY|SUPPLEMENTARY|DUP, which keeps QCFAIL). To match the
// synth push-down exactly instead, drop BAM_FQCFAIL from the line below.
// ----------------------------------------------------------------------------
static const uint32_t FLAG_FILTER =
    BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP;

static const uint64_t READ_THRESHOLD = 5000000ULL; // reads per work unit before splitting a contig
static const int      MAX_UNITS_PER_CONTIG = 64;
static const int      MAX_CLASSES = 32;            // 32-bit per-locus class mask

// Lean position-only SNP record. Do not reuse SNPData/ChromSNPs; this tool needs
// no genotypes, ref/alt alleles, or per-sample arrays.
struct PosClass {
    int      pos;   // 0-based reference position
    uint32_t mask;  // bit i set if the locus is present in input VCF i
};

// Per-cell histogram: one count map per class, with ANY at index n_classes.
struct CellHist {
    vector<robin_hood::unordered_map<int, uint64_t>> by_class; // size n_classes + 1
};

struct WorkUnit {
    int      tid;
    int      start;      // half-open window start (inclusive)
    int      end;        // half-open window end (exclusive); INT_MAX on final window
    uint64_t est_reads;  // for largest-first sorting
};

static void die(const string& msg){
    fprintf(stderr, "ERROR: %s\n", msg.c_str());
    exit(1);
}

static void usage(){
    fprintf(stderr,
        "snps_per_read: per-read SNP-locus coverage distribution, per cell and per SNP class.\n"
        "\n"
        "Usage:\n"
        "  snps_per_read \\\n"
        "    --bam PATH \\\n"
        "    --library LABEL \\\n"
        "    --vcf LABEL:PATH        (repeatable, one per SNP class) \\\n"
        "    --output_prefix PATH \\\n"
        "    [--barcodes PATH]       (optional CB whitelist, one barcode per line) \\\n"
        "    [--min_vq INT]          (SNP QUAL floor, default 50) \\\n"
        "    [--threads INT]         (default: all available) \\\n"
        "    [--dump_reads]          (optional raw per-read output, off by default)\n");
}

// Derive a library label from the BAM basename when --library is omitted.
static string derive_library(const string& bam_path){
    size_t slash = bam_path.find_last_of('/');
    string base = (slash == string::npos) ? bam_path : bam_path.substr(slash + 1);
    if (base.size() > 4 && base.compare(base.size() - 4, 4, ".bam") == 0){
        base = base.substr(0, base.size() - 4);
    }
    return base;
}

// ----------------------------------------------------------------------------
// VCF load (spec section 6). Custom lean reader, positions only.
//
// Builds, per BAM tid, a pos -> class-mask map (OR-ing bits across inputs and
// across duplicate records at the same position), then flattens and sorts into
// a vector<PosClass> per contig.
// ----------------------------------------------------------------------------
static void load_vcfs(const vector<string>& vcf_paths,
                      const vector<string>& class_labels,
                      bam_hdr_t* bam_hdr,
                      int min_vq,
                      robin_hood::unordered_map<int, vector<PosClass>>& merged){

    int n_classes = (int)vcf_paths.size();
    // Per-tid pos -> mask, accumulated across all inputs.
    robin_hood::unordered_map<int, robin_hood::unordered_map<int, uint32_t>> build;

    for (int i = 0; i < n_classes; i++){
        const string& path = vcf_paths[i];
        uint32_t bit = (uint32_t)1 << i;

        htsFile* vfp = bcf_open(path.c_str(), "r");
        if (!vfp) die("could not open VCF/BCF: " + path);
        bcf_hdr_t* vhdr = bcf_hdr_read(vfp);
        if (!vhdr){ bcf_close(vfp); die("could not read VCF/BCF header: " + path); }
        bcf1_t* vrec = bcf_init();

        uint64_t n_read = 0, n_vq = 0, n_mapped = 0, n_offcontig = 0;
        int last_rid = -1;
        int last_tid = -1; // cache CHROM->tid mapping across consecutive records

        while (bcf_read(vfp, vhdr, vrec) >= 0){
            n_read++;
            // QUAL floor. Missing QUAL is NaN; (NaN < min_vq) is false, so missing
            // QUAL passes, matching the demux genotype-path QUAL handling.
            if (vrec->qual < (float)min_vq) continue;
            n_vq++;

            if (vrec->rid != last_rid){
                last_rid = vrec->rid;
                const char* chrom = bcf_hdr_id2name(vhdr, vrec->rid);
                last_tid = (chrom != NULL) ? sam_hdr_name2tid(bam_hdr, chrom) : -1;
            }
            if (last_tid < 0){ n_offcontig++; continue; }
            n_mapped++;

            int pos = (int)vrec->pos;
            build[last_tid][pos] |= bit; // OR across inputs and duplicate records
        }

        fprintf(stderr,
            "  [%s] %s: %llu records read, %llu pass VQ>=%d, %llu mapped to BAM contigs, %llu on off-BAM contigs\n",
            class_labels[i].c_str(), path.c_str(),
            (unsigned long long)n_read, (unsigned long long)n_vq, min_vq,
            (unsigned long long)n_mapped, (unsigned long long)n_offcontig);

        bcf_destroy(vrec);
        bcf_hdr_destroy(vhdr);
        bcf_close(vfp);
    }

    // Flatten each contig's pos->mask map into a sorted vector<PosClass>.
    uint64_t total_loci = 0, shared_loci = 0;
    for (auto& kv : build){
        int tid = kv.first;
        auto& posmap = kv.second;
        vector<PosClass>& vec = merged[tid];
        vec.reserve(posmap.size());
        for (auto& pm : posmap){
            vec.push_back(PosClass{pm.first, pm.second});
        }
        sort(vec.begin(), vec.end(),
             [](const PosClass& a, const PosClass& b){ return a.pos < b.pos; });
        total_loci += vec.size();
        for (auto& pc : vec){
            // popcount of mask > 1 means the locus is shared across classes.
            if ((pc.mask & (pc.mask - 1)) != 0) shared_loci++;
        }
    }
    fprintf(stderr,
        "  Loaded %llu distinct loci across %zu contigs (%llu shared across >1 class).\n",
        (unsigned long long)total_loci, merged.size(), (unsigned long long)shared_loci);
}

// ----------------------------------------------------------------------------
// Work-unit pool (spec section 7.1). Contiguous, non-overlapping position
// partitions so read ownership by alignment start is unambiguous.
// ----------------------------------------------------------------------------
static void build_work_units(bam_hdr_t* hdr,
                             hts_idx_t* idx,
                             vector<WorkUnit>& units){
    int n_chroms = hdr->n_targets;
    for (int tid = 0; tid < n_chroms; tid++){
        uint64_t mapped = 0, unmapped = 0;
        if (hts_idx_get_stat(idx, tid, &mapped, &unmapped) < 0){
            mapped = 0; // no index stat for this contig; one unit owns [0, INT_MAX)
        }
        int64_t clen = (int64_t)hdr->target_len[tid];

        if (mapped <= READ_THRESHOLD){
            units.push_back(WorkUnit{tid, 0, INT_MAX, mapped});
            continue;
        }

        int n_windows = (int)((mapped + READ_THRESHOLD - 1) / READ_THRESHOLD);
        if (n_windows > MAX_UNITS_PER_CONTIG) n_windows = MAX_UNITS_PER_CONTIG;
        if (clen <= 0){ // defensive: degenerate length, fall back to one unit
            units.push_back(WorkUnit{tid, 0, INT_MAX, mapped});
            continue;
        }
        int64_t wlen = (clen + n_windows - 1) / n_windows;
        uint64_t per = mapped / (uint64_t)n_windows;
        for (int w = 0; w < n_windows; w++){
            int start = (int)((int64_t)w * wlen);
            int end   = (w == n_windows - 1) ? INT_MAX : (int)((int64_t)(w + 1) * wlen);
            units.push_back(WorkUnit{tid, start, end, per});
        }
    }
    // Largest-first so heavy units are claimed early and do not form a tail.
    sort(units.begin(), units.end(),
         [](const WorkUnit& a, const WorkUnit& b){ return a.est_reads > b.est_reads; });
}

// Fetch-or-create a CellHist for a CB inside a thread-local map.
static inline CellHist& get_cell(robin_hood::unordered_map<string, CellHist>& m,
                                 const string& cb, int n_classes){
    auto it = m.find(cb);
    if (it != m.end()) return it->second;
    CellHist& ch = m[cb];
    ch.by_class.resize((size_t)n_classes + 1);
    return ch;
}

int main(int argc, char** argv){
    string bam_path, library, output_prefix, barcodes_path;
    vector<string> class_labels;
    vector<string> vcf_paths;
    int  min_vq = 50;       // matches demux / vcfhost / vcf_loader_daemon default
    int  n_threads = 0;     // 0 => all available
    bool dump_reads = false;
    bool library_set = false;

    for (int i = 1; i < argc; i++){
        string a = argv[i];
        auto need = [&](string& dest){
            if (i + 1 >= argc) die("missing value after " + a);
            dest = argv[++i];
        };
        if      (a == "--bam") need(bam_path);
        else if (a == "--library"){ need(library); library_set = true; }
        else if (a == "--output_prefix") need(output_prefix);
        else if (a == "--barcodes") need(barcodes_path);
        else if (a == "--vcf"){
            string spec; need(spec);
            size_t colon = spec.find(':');
            if (colon == string::npos || colon == 0 || colon + 1 >= spec.size())
                die("--vcf must be LABEL:PATH, got: " + spec);
            class_labels.push_back(spec.substr(0, colon));
            vcf_paths.push_back(spec.substr(colon + 1));
        }
        else if (a == "--min_vq"){ string t; need(t); min_vq = atoi(t.c_str()); }
        else if (a == "--threads"){ string t; need(t); n_threads = atoi(t.c_str()); }
        else if (a == "--dump_reads") dump_reads = true;
        else if (a == "--help" || a == "-h"){ usage(); return 0; }
        else die("unknown argument: " + a);
    }

    if (bam_path.empty() || output_prefix.empty() || vcf_paths.empty()){
        usage();
        return 1;
    }
    if ((int)vcf_paths.size() > MAX_CLASSES)
        die("at most 32 SNP classes are supported (32-bit mask)");
    // Reject duplicate class labels: they would collide in the output columns.
    for (size_t i = 0; i < class_labels.size(); i++){
        if (class_labels[i] == "ANY")
            die("class label 'ANY' is reserved");
        for (size_t j = i + 1; j < class_labels.size(); j++){
            if (class_labels[i] == class_labels[j])
                die("duplicate class label: " + class_labels[i]);
        }
    }

    if (!library_set){
        library = derive_library(bam_path);
        fprintf(stderr, "WARNING: --library not given; derived '%s' from BAM basename.\n",
                library.c_str());
    }

    int n_classes = (int)vcf_paths.size();
    int ANY = n_classes; // index of the ANY slot in CellHist.by_class

    if (n_threads <= 0) n_threads = omp_get_max_threads();
    omp_set_num_threads(n_threads);

    // Optional CB whitelist (verbatim match, consistent with synth_bam_filter.py).
    bool has_whitelist = !barcodes_path.empty();
    robin_hood::unordered_set<string> whitelist;
    if (has_whitelist){
        ifstream bf(barcodes_path.c_str());
        if (!bf) die("could not open --barcodes file: " + barcodes_path);
        string line;
        while (getline(bf, line)){
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (!line.empty()) whitelist.insert(line);
        }
        fprintf(stderr, "Loaded %zu whitelist barcodes from %s\n",
                whitelist.size(), barcodes_path.c_str());
    }

    // Open the BAM once to read the header and index for VCF mapping + work units.
    htsFile* bam0 = hts_open(bam_path.c_str(), "r");
    if (!bam0) die("could not open BAM: " + bam_path);
    bam_hdr_t* hdr0 = sam_hdr_read(bam0);
    if (!hdr0) die("could not read BAM header: " + bam_path);
    hts_idx_t* idx0 = sam_index_load(bam0, bam_path.c_str());
    if (!idx0) die("could not load BAM index (.bai/.csi) for: " + bam_path);

    fprintf(stderr, "Loading %d SNP class panel(s) (min_vq=%d)...\n", n_classes, min_vq);
    robin_hood::unordered_map<int, vector<PosClass>> merged;
    load_vcfs(vcf_paths, class_labels, hdr0, min_vq, merged);

    vector<WorkUnit> units;
    build_work_units(hdr0, idx0, units);
    fprintf(stderr, "Built %zu work units over %d contigs using %d threads.\n",
            units.size(), hdr0->n_targets, n_threads);

    hts_idx_destroy(idx0);
    bam_hdr_destroy(hdr0);
    hts_close(bam0);

    // Per-thread accumulation (no locks in the scan).
    vector<robin_hood::unordered_map<string, CellHist>> thread_maps(n_threads);
    // Per-thread raw dump temp files (only used with --dump_reads).
    vector<string> dump_paths(n_threads);
    atomic<long> units_done(0);
    atomic<long> total_eligible(0);

    #pragma omp parallel
    {
        int tid_thread = omp_get_thread_num();
        auto& local = thread_maps[tid_thread];

        htsFile*    bam_fp = hts_open(bam_path.c_str(), "r");
        bam_hdr_t*  header = bam_fp ? sam_hdr_read(bam_fp) : NULL;
        hts_idx_t*  idx    = bam_fp ? sam_index_load(bam_fp, bam_path.c_str()) : NULL;
        bam1_t*     rec    = bam_init1();

        gzFile dump = Z_NULL;
        if (dump_reads){
            string dp = output_prefix + ".per_read.t" + to_string(tid_thread) + ".tsv.gz";
            dump_paths[tid_thread] = dp;
            dump = gzopen(dp.c_str(), "wb");
        }

        if (!bam_fp || !header || !idx){
            fprintf(stderr, "ERROR: thread %d could not open BAM/header/index\n", tid_thread);
        } else {
            long local_eligible = 0;
            // Per-read class counters (reused per read; size MAX_CLASSES).
            uint32_t kcls[MAX_CLASSES];

            #pragma omp for schedule(dynamic, 1)
            for (size_t u = 0; u < units.size(); u++){
                const WorkUnit& wu = units[u];
                int tid = wu.tid;

                auto mit = merged.find(tid);
                const vector<PosClass>* vec = (mit != merged.end()) ? &mit->second : NULL;
                auto vbeg = vec ? vec->begin() : vector<PosClass>::const_iterator();
                auto vend = vec ? vec->end()   : vector<PosClass>::const_iterator();

                hts_itr_t* iter = sam_itr_queryi(idx, tid, wu.start, wu.end);
                if (!iter){
                    long d = ++units_done;
                    if (d % 200 == 0)
                        fprintf(stderr, "\rProgress: %ld/%zu units, %ld eligible reads ",
                                d, units.size(), total_eligible.load());
                    continue;
                }

                while (sam_itr_next(bam_fp, iter, rec) >= 0){
                    // Flag filter.
                    if (rec->core.flag & FLAG_FILTER) continue;
                    // Ownership: this unit owns the read iff its alignment start
                    // falls in [start, end). Contiguous half-open windows tile the
                    // contig, so every read is owned by exactly one unit.
                    int rpos = (int)rec->core.pos;
                    if (rpos < wu.start || rpos >= wu.end) continue;
                    // CB required for cell aggregation.
                    uint8_t* cb_tag = bam_aux_get(rec, "CB");
                    if (!cb_tag) continue;
                    const char* cb_c = bam_aux2Z(cb_tag);
                    if (!cb_c) continue;
                    string cb(cb_c);
                    if (has_whitelist && whitelist.find(cb) == whitelist.end()) continue;

                    local_eligible++;

                    // Per-read counters.
                    for (int c = 0; c < n_classes; c++) kcls[c] = 0;
                    uint32_t kany = 0;

                    int read_start = rpos;
                    int read_end   = (int)bam_endpos(rec); // exclusive

                    // Fast path: if no locus lies in [read_start, read_end), all
                    // class counts are 0. One lower_bound test, skip the walk.
                    bool any_locus_in_span = false;
                    vector<PosClass>::const_iterator it = vend;
                    if (vec){
                        it = lower_bound(vbeg, vend, read_start,
                            [](const PosClass& p, int v){ return p.pos < v; });
                        if (it != vend && it->pos < read_end) any_locus_in_span = true;
                    }

                    if (any_locus_in_span){
                        // CIGAR coverage walk. M/=/X place aligned bases (cover);
                        // D/N consume reference but place no base (no coverage);
                        // I/S consume read only; H/P consume neither.
                        uint32_t* cig = bam_get_cigar(rec);
                        int ref_pos = read_start;
                        for (uint32_t ci = 0; ci < rec->core.n_cigar; ci++){
                            int op  = bam_cigar_op(cig[ci]);
                            int len = bam_cigar_oplen(cig[ci]);
                            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF){
                                int bend = ref_pos + len;
                                // Skip loci left behind in a preceding D/N gap.
                                while (it != vend && it->pos < ref_pos) ++it;
                                while (it != vend && it->pos < bend){
                                    kany++;
                                    uint32_t m = it->mask;
                                    while (m){
                                        int b = __builtin_ctz(m);
                                        kcls[b]++;
                                        m &= (m - 1);
                                    }
                                    ++it;
                                }
                                ref_pos = bend;
                            } else if (op == BAM_CDEL || op == BAM_CREF_SKIP){
                                ref_pos += len; // advance reference, count nothing
                            }
                            // BAM_CINS, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD:
                            // do not consume reference; nothing to do.
                        }
                    }

                    // Tally one count per class (and ANY) into this CB's histogram.
                    CellHist& ch = get_cell(local, cb, n_classes);
                    for (int c = 0; c < n_classes; c++)
                        ch.by_class[c][(int)kcls[c]] += 1;
                    ch.by_class[ANY][(int)kany] += 1;

                    if (dump_reads && dump != Z_NULL){
                        // CB  tid  ref_start  k_ANY  k_<label0>  k_<label1> ...
                        gzprintf(dump, "%s\t%d\t%d\t%u", cb.c_str(), tid, read_start, kany);
                        for (int c = 0; c < n_classes; c++)
                            gzprintf(dump, "\t%u", kcls[c]);
                        gzprintf(dump, "\n");
                    }
                }
                hts_itr_destroy(iter);

                long d = ++units_done;
                if (d % 200 == 0)
                    fprintf(stderr, "\rProgress: %ld/%zu units, %ld eligible reads ",
                            d, units.size(), total_eligible.load());
            }
            total_eligible += local_eligible;
        }

        if (dump != Z_NULL) gzclose(dump);
        if (rec)    bam_destroy1(rec);
        if (idx)    hts_idx_destroy(idx);
        if (header) bam_hdr_destroy(header);
        if (bam_fp) hts_close(bam_fp);
    }
    fprintf(stderr, "\rScan complete: %zu/%zu units, %ld eligible reads.            \n",
            units.size(), units.size(), total_eligible.load());

    // ------------------------------------------------------------------------
    // Merge per-thread maps into one global map (spec section 7.5).
    // ------------------------------------------------------------------------
    robin_hood::unordered_map<string, CellHist> global;
    for (int t = 0; t < n_threads; t++){
        auto& tm = thread_maps[t];
        for (auto& kv : tm){
            const string& cb = kv.first;
            CellHist& src = kv.second;
            auto git = global.find(cb);
            if (git == global.end()){
                global[cb] = std::move(src);
            } else {
                CellHist& dst = git->second;
                for (int c = 0; c <= n_classes; c++){
                    for (auto& nc : src.by_class[c])
                        dst.by_class[c][nc.first] += nc.second;
                }
            }
        }
        tm = robin_hood::unordered_map<string, CellHist>(); // free as merged
    }
    fprintf(stderr, "Merged %zu cells.\n", global.size());

    // Deterministic CB order so single-thread and multi-thread outputs match.
    vector<string> cbs;
    cbs.reserve(global.size());
    for (auto& kv : global) cbs.push_back(kv.first);
    sort(cbs.begin(), cbs.end());

    // ------------------------------------------------------------------------
    // 4.1 Per-cell distribution (long format), always written.
    // Also build the rolled-up histogram (4.2) and summary inputs in one pass.
    // ------------------------------------------------------------------------
    string per_cell_path = output_prefix + ".per_cell.tsv.gz";
    gzFile pc = gzopen(per_cell_path.c_str(), "wb");
    if (!pc) die("could not open output: " + per_cell_path);
    gzprintf(pc, "library\tCB\tclass\tn_snps\tn_reads\n");

    // Rolled-up histogram: per class index (0..n_classes, ANY last), n_snps -> count.
    vector<robin_hood::unordered_map<int, uint64_t>> hist((size_t)n_classes + 1);

    auto class_name = [&](int c) -> const string& {
        static const string any_label = "ANY";
        return (c == ANY) ? any_label : class_labels[c];
    };

    for (const string& cb : cbs){
        CellHist& ch = global[cb];
        // Eligible read total for this cell (identical across classes).
        uint64_t cell_total = 0;
        for (auto& nc : ch.by_class[ANY]) cell_total += nc.second;
        if (cell_total == 0) continue;

        for (int c = 0; c <= n_classes; c++){
            auto& cm = ch.by_class[c];
            // Emit observed (count>0) bins plus a guaranteed n_snps=0 row.
            vector<int> keys;
            keys.reserve(cm.size() + 1);
            bool has_zero = false;
            for (auto& nc : cm){
                keys.push_back(nc.first);
                if (nc.first == 0) has_zero = true;
            }
            if (!has_zero) keys.push_back(0);
            sort(keys.begin(), keys.end());
            for (int k : keys){
                auto f = cm.find(k);
                uint64_t cnt = (f != cm.end()) ? f->second : 0;
                gzprintf(pc, "%s\t%s\t%s\t%d\t%llu\n",
                         library.c_str(), cb.c_str(), class_name(c).c_str(),
                         k, (unsigned long long)cnt);
                if (cnt > 0) hist[c][k] += cnt;
            }
        }
    }
    gzclose(pc);

    // ------------------------------------------------------------------------
    // 4.2 Rolled-up histogram, always written.
    // ------------------------------------------------------------------------
    string hist_path = output_prefix + ".histogram.tsv";
    FILE* hf = fopen(hist_path.c_str(), "w");
    if (!hf) die("could not open output: " + hist_path);
    fprintf(hf, "class\tn_snps\tn_reads\n");
    for (int c = 0; c <= n_classes; c++){
        vector<int> keys;
        keys.reserve(hist[c].size());
        for (auto& nc : hist[c]) keys.push_back(nc.first);
        sort(keys.begin(), keys.end());
        for (int k : keys)
            fprintf(hf, "%s\t%d\t%llu\n",
                    class_name(c).c_str(), k, (unsigned long long)hist[c][k]);
    }
    fclose(hf);

    // ------------------------------------------------------------------------
    // 4.3 Summary stats, always written and echoed to stderr.
    // Mean/median/max computed from the rolled-up histogram.
    // ------------------------------------------------------------------------
    string sum_path = output_prefix + ".summary.tsv";
    FILE* sf = fopen(sum_path.c_str(), "w");
    if (!sf) die("could not open output: " + sum_path);
    const char* sum_hdr = "class\ttotal_reads\tmean_snps_per_read\tmedian_snps_per_read\tmax_snps_per_read\n";
    fprintf(sf, "%s", sum_hdr);
    fprintf(stderr, "\n%s", sum_hdr);

    for (int c = 0; c <= n_classes; c++){
        vector<pair<int, uint64_t>> bins;
        bins.reserve(hist[c].size());
        for (auto& nc : hist[c]) bins.emplace_back(nc.first, nc.second);
        sort(bins.begin(), bins.end(),
             [](const pair<int,uint64_t>& a, const pair<int,uint64_t>& b){
                 return a.first < b.first;
             });

        uint64_t total = 0;
        double   sum_k = 0.0;
        int      maxk  = 0;
        for (auto& b : bins){
            total += b.second;
            sum_k += (double)b.first * (double)b.second;
            if (b.second > 0 && b.first > maxk) maxk = b.first;
        }
        double mean = (total > 0) ? (sum_k / (double)total) : 0.0;

        // Median from the cumulative histogram (lower/upper middle averaged).
        double median = 0.0;
        if (total > 0){
            uint64_t r1 = (total + 1) / 2; // lower middle rank (1-based)
            uint64_t r2 = total / 2 + 1;   // upper middle rank (1-based)
            int v1 = 0, v2 = 0;
            uint64_t cum = 0;
            bool got1 = false, got2 = false;
            for (auto& b : bins){
                cum += b.second;
                if (!got1 && cum >= r1){ v1 = b.first; got1 = true; }
                if (!got2 && cum >= r2){ v2 = b.first; got2 = true; }
                if (got1 && got2) break;
            }
            median = ((double)v1 + (double)v2) / 2.0;
        }

        fprintf(sf, "%s\t%llu\t%.6f\t%.1f\t%d\n",
                class_name(c).c_str(), (unsigned long long)total, mean, median, maxk);
        fprintf(stderr, "%s\t%llu\t%.6f\t%.1f\t%d\n",
                class_name(c).c_str(), (unsigned long long)total, mean, median, maxk);
    }
    fclose(sf);

    // ------------------------------------------------------------------------
    // 4.4 Raw per-read dump (optional). Concatenate per-thread gz members into
    // one multi-member gzip with a header member first.
    // ------------------------------------------------------------------------
    if (dump_reads){
        string dump_final = output_prefix + ".per_read.tsv.gz";
        // Header member.
        gzFile hd = gzopen(dump_final.c_str(), "wb");
        if (!hd) die("could not open output: " + dump_final);
        gzprintf(hd, "CB\ttid\tref_start\tk_ANY");
        for (int c = 0; c < n_classes; c++)
            gzprintf(hd, "\tk_%s", class_labels[c].c_str());
        gzprintf(hd, "\n");
        gzclose(hd);
        // Append each thread's raw gz bytes, then remove the temp.
        FILE* out = fopen(dump_final.c_str(), "ab");
        if (!out) die("could not reopen output for append: " + dump_final);
        char buf[1 << 16];
        for (int t = 0; t < n_threads; t++){
            const string& dp = dump_paths[t];
            if (dp.empty()) continue;
            FILE* in = fopen(dp.c_str(), "rb");
            if (!in) continue;
            size_t r;
            while ((r = fread(buf, 1, sizeof(buf), in)) > 0) fwrite(buf, 1, r, out);
            fclose(in);
            remove(dp.c_str());
        }
        fclose(out);
        fprintf(stderr, "Wrote raw per-read dump: %s\n", dump_final.c_str());
    }

    fprintf(stderr, "Wrote %s, %s, %s\n",
            per_cell_path.c_str(), hist_path.c_str(), sum_path.c_str());
    return 0;
}
