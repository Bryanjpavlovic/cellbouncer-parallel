// ============================================================================
// genotype_scrub_bam
// ----------------------------------------------------------------------------
// Corrected genotype-conditioned source-BAM scrubbing for the Tet2025 ambient
// benchmark (dosage-projection model). Streams one library/cell-cache BAM twice
// and, for every assigned cell, rebalances the observed ref:alt read count at
// each covered panel site to the home identity's expected alt fraction under the
// same genotype, GQ, error-rate, and dosage basis the scoring counter uses. Only
// the minimum number of reads is flipped to hit the target; read count is
// preserved exactly (base-only edits, plus MD/NM removal on rewritten records).
//
// This replaces the earlier homozygous-endpoint rule, which rewrote reads only
// where the home genotype predicted a single allele (singlet hom-ref/hom-alt and
// fusion homozygous-agree endpoints) and left every heterozygous-home and fusion
// disagreement site untouched. Those skipped sites are exactly where demux and
// the contamination estimators read cross-individual and cross-species identity,
// and for tetraploid fusion cells they are the majority of the informative
// panel. Rebalancing the per-site ref:alt ratio to the home expectation erases
// everything about a contaminant that any consumer downstream can detect,
// without needing to identify which reads were foreign.
//
// Target alt fraction by home genotype (nalt = per-genome dosage 0/1/2):
//   Singlet:  base_p = nalt_i / 2          (hom-ref 0, het 0.5, hom-alt 1)
//   Fusion:   base_p = (nalt_i + nalt_j)/4  (subsumes both-hom-ref 0,
//             both-hom-alt 1, and every disagreement/het-constituent case at
//             the balanced dosage expectation r = 0.5)
//   f_target = adjust_p_err(base_p, e_r, e_a) = base_p - base_p*e_a + (1-base_p)*e_r
//
// Both the interindividual panel and the species panel are scrubbed in one pass,
// over the union of their assayed sites. The species VCF shares the demux sample
// set exactly (same individuals, same order), so a cell's home genotype columns
// are identical on both panels; only the site set differs. Per-panel input
// signatures and a merged-union signature are recorded in the summary so a
// species-mode consumer can reject a BAM lacking a species-panel scrub.
//
// Wobble-preserving noise model: the target alt count is drawn from
// Binomial(n, f_target) with a deterministic per-cell-site seed (a hash of
// library, barcode, contig, position, panel id, and policy hash), applied
// uniformly including homozygous sites, so the whole BAM carries one consistent
// model-consistent noise structure and reruns are byte-identical. A hard_snap
// model (round(n * f_target)) is available for validation.
//
// GQ basis defaults to cellbouncer_permissive, matching the effective behavior
// of the demux/counting VCF loader (a present genotype is used regardless of the
// GQ floor); a low or missing GQ makes a sample uncovered only where the loader
// would treat it as uncovered. This keeps the scrub and the scoring counter on
// one het/hom/covered/uncovered basis.
//
// Standalone: htslib + zlib + pthread only. No cellbouncer object, no
// htswrapper link. Reuses the logic (not the code) of demux_vcf_hts
// (read_vcf_chroms genotype load), demux_parallel_hts (get_base_at_pos CIGAR
// walk and effective GQ basis), snps_per_read (coverage walk / flag mask), and
// bam_cb_cache_extract (BAM read/write/CSI skeleton).
//
// Provenance: V2_R3 (filename stable; consumed by the CellBouncer Makefile
// target genotype_scrub_bam, the SCRUB_SOURCES stage, and the cell-cache scrub
// stage, so it is not version suffixed).
//
// Revision history:
//   V2_R3: B34 fusion-sender completeness (Tet2025 NewBench remaining-closure
//          spec, Task 2). Added --fusion-r-table PATH: a header-labeled TSV keyed
//          by (library, barcode) giving each fusion cell's endogenous genome-A
//          fraction r, produced by the B35 allele typing. A fusion cell listed in
//          the table is rebalanced at its disagreement sites toward its own typed
//          r via fusion_r*pA+(1-fusion_r)*pB with r taken per cell (new
//          CellRec.fusion_r), so a fusion SENDER is cleaned of third-party leakage
//          without erasing its real, imbalanced A/B. Cells absent from the table,
//          singlets, and every cell in a run with no table fall back to the global
//          --fusion-r (default 0.5), so a receiver/singlet-sender scrub is
//          byte-identical to V2_R2. The run seed and per-cell seed are unchanged,
//          so only a table-listed fusion cell's target moves; the wobble draw stays
//          deterministic and a rerun is byte-identical. Summary gains
//          fusion_r_source, fusion_r_table_signature, fusion_r_table_rows, and
//          fusion_r_table_cells so the Python resolver can gate a fusion-sender
//          scrub on a matching r-table signature the way it gates the species panel.
//          An out-of-range or non-numeric table value fails loud rather than
//          silently reverting to 0.5.
//   V2_R2: Integration pass, reconciling two parallel builds of this rewrite
//          against the now-available coding doc (Tet2025_NewBench_Scrub_Coding_
//          Doc_V1_R1.md). Kept the one-pass union architecture (both panels in a
//          single invocation via --species-vcf), which the coding doc Section 2.6
//          mandates, over the alternative one-invocation-per-panel design. Folded
//          in the two additions the coding doc and the scrub audit endorse:
//          (1) --fusion-r (default 0.5) generalizes the fusion target to
//          r*p_A+(1-r)*p_B; the default is the balanced dosage expectation the
//          benchmark forces (spec decision 6), and the knob exists only for a
//          future fusion-sender run. (2) --mask emits a per-cell-site
//          rebalance-target TSV (barcode, identity, panel, contig, pos, ref, alt,
//          home_kind, raw_p, target_frac, n_ref, n_alt, n_other, n_refalt,
//          target_alt, dir, flips) for the het-site deviation validation. Added a
//          scrubbable-cell-site category census to the summary (the sites the old
//          homozygous-only rule skipped), n_other tracking for covered
//          third-base reads (counted for the mask, never rebalanced), and a panel
//          bitmask so a position assayed by both panels reports 'U'. Default
//          fusion behavior is byte-identical to V2_R1 (fusion_r=0.5 gives the same
//          (nalt_A+nalt_B)/4 target).
//   V2_R1: Corrected scrub model for the fresh SYNTH_AMBIENT_BENCH_V3 run.
//          Replaced the per-read homozygous-endpoint rewrite with a two-pass
//          per-cell-site dosage projection: pass 1 counts observed ref:alt per
//          (cell, site); a wobble-preserving Binomial(n, f_target) draw with a
//          deterministic per-cell-site seed sets the target alt count; pass 2
//          flips the minimum number of reads (in coordinate order) to hit it.
//          f_target uses the singlet nalt/2 and fusion (nalt_A+nalt_B)/4 dosage
//          basis with adjust_p_err error adjustment, so het-home and fusion
//          disagreement sites are rebalanced, not skipped. Added the species
//          panel via --species-vcf, loaded into the same union site index with
//          the same home columns; per-panel and merged-union input signatures
//          recorded. GQ policy default changed to cellbouncer_permissive to
//          match the scoring counter. Site-load prune widened from
//          homozygous-only to any-covered so het-home sites are retained.
//          fusion_policy reported as dosage_projection; noise_model, error_ref,
//          error_alt recorded in the summary.
//   V1_R4..V1_R1: prior homozygous-endpoint scrubber (see git history). The
//          atomic-publish, sample-order guard, 10x -1 barcode normalization,
//          tetra_refine-over-demux precedence, read-count invariant, and CSI
//          skeleton are carried forward unchanged.
// ============================================================================

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/thread_pool.h>

using namespace std;

// Synth CB push-down mask: UNMAP|SECONDARY|SUPPLEMENTARY|DUP (keeps QCFAIL).
// Reads matching this mask are written unmodified and never counted or scrubbed,
// so the scrubbed read set matches exactly the reads the synth pipeline consumes
// and pass-1 counts line up with pass-2 flips.
static const uint32_t SCRUB_SKIP_MASK = 0xD04;
static const int CSI_MIN_SHIFT = 14;

// Panel membership is a bitmask so a position assayed by both panels can be
// reported as 'U' in the audit mask. The deterministic seed uses a stable
// primary-panel id (interindividual wins when a site is in both), matching the
// dedup that keeps the interindividual copy.
static const uint8_t PANEL_BIT_INTER   = 0x1;
static const uint8_t PANEL_BIT_SPECIES = 0x2;
static inline uint8_t panel_primary_seed_id(uint8_t mask){
    return (mask & PANEL_BIT_INTER) ? 0 : 1;
}

// ----------------------------------------------------------------------------
// Small helpers
// ----------------------------------------------------------------------------
static void die(const string& msg){
    fprintf(stderr, "ERROR: %s\n", msg.c_str());
    exit(1);
}

static void usage(){
    fprintf(stderr,
        "genotype_scrub_bam: rebalance per-cell ref:alt read counts to the home\n"
        "genotype's expected alt fraction (dosage projection) at every covered\n"
        "interindividual and species panel site, for one library/cell-cache BAM.\n"
        "\n"
        "Usage:\n"
        "  genotype_scrub_bam \\\n"
        "    --bam INPUT.gex.bam \\\n"
        "    --vcf tet.vars.downsampled_20M.bcf \\\n"
        "    --assignments libN_demuxed.assignments \\\n"
        "    --samples libN_demuxed.samples \\\n"
        "    --output OUTPUT.scrubbed.bam \\\n"
        "    --report OUTPUT.scrubbed.report.tsv[.gz] \\\n"
        "    --summary OUTPUT.scrubbed.summary.json \\\n"
        "    [--species-vcf tet.vars.species_20M.bcf] \\\n"
        "    [--refined libN.refined_assignments] \\\n"
        "    [--barcodes filtered/barcodes.tsv.gz] \\\n"
        "    [--library LABEL] \\\n"
        "    [--threads INT]        (htslib BGZF codec threads, default 8) \\\n"
        "    [--min-vq INT]         (SNP QUAL floor, default 50) \\\n"
        "    [--min-gq INT]         (per-sample GQ floor, default 30) \\\n"
        "    [--tag-policy strip_md_nm|keep]   (default strip_md_nm) \\\n"
        "    [--gq-policy cellbouncer_permissive|min_gq_enforced] (default cellbouncer_permissive) \\\n"
        "    [--noise-model wobble|hard_snap]  (default wobble) \\\n"
        "    [--fusion-r FLOAT]     (endogenous genome-A fraction for fusion sites, default 0.5) \\\n"
        "    [--fusion-r-table PATH] (per-cell fusion-r TSV keyed by library,barcode; overrides --fusion-r per fusion cell) \\\n"
        "    [--mask PATH[.gz]]     (per-cell-site rebalance-target audit mask) \\\n"
        "    [--error-ref FLOAT]    (e_r for adjust_p_err, default 0.005) \\\n"
        "    [--error-alt FLOAT]    (e_a for adjust_p_err, default 0.005) \\\n"
        "    [--tmp-token STR]      (job-unique temp suffix; default pid) \\\n"
        "    [--policy-version STR] (recorded in summary; default builtin) \\\n"
        "    [--policy-hash STR]    (recorded in summary for stale-set detection)\n");
}

// Bare 16bp barcode: everything before the first '-' (drops the 10x -1 suffix).
static inline string normalize_bc(const char* s){
    const char* dash = strchr(s, '-');
    return dash ? string(s, (size_t)(dash - s)) : string(s);
}
static inline string normalize_bc(const string& s){
    size_t dash = s.find('-');
    return (dash == string::npos) ? s : s.substr(0, dash);
}

// size:mtime_ns signature, or "MISSING" when the file cannot be stat'd. Byte
// identical to source_bams.file_signature on the Python side.
static string file_sig(const string& path){
    if (path.empty()) return "MISSING";
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return "MISSING";
    char buf[128];
    snprintf(buf, sizeof(buf), "%lld:%lld.%09ld",
             (long long)st.st_size, (long long)st.st_mtime,
             (long)st.st_mtim.tv_nsec);
    return string(buf);
}

static string json_escape(const string& s){
    string out;
    out.reserve(s.size() + 8);
    for (char c : s){
        if (c == '"' || c == '\\'){ out.push_back('\\'); out.push_back(c); }
        else if (c == '\n'){ out += "\\n"; }
        else if (c == '\t'){ out += "\\t"; }
        else out.push_back(c);
    }
    return out;
}

// adjust_p_err, identical to demux_parallel_llr.cpp / demux_vcf_llr.cpp so the
// scrub targets the same expected alt fraction the scoring path computes.
static inline double adjust_p_err(double p, double e_r, double e_a){
    return p - p * e_a + (1.0 - p) * e_r;
}

// ----------------------------------------------------------------------------
// Deterministic hashing / RNG for the wobble draw. Fully reproducible across
// compilers and machines (no std::binomial_distribution, whose stream is
// implementation-defined).
// ----------------------------------------------------------------------------
static inline uint64_t splitmix64(uint64_t& x){
    uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static inline uint64_t mix64(uint64_t a, uint64_t b){
    uint64_t x = a ^ (b + 0x9E3779B97F4A7C15ULL + (a << 6) + (a >> 2));
    uint64_t s = x;
    return splitmix64(s);
}
static inline uint64_t fnv1a64(const string& s){
    uint64_t h = 1469598103934665603ULL;
    for (unsigned char c : s){ h ^= (uint64_t)c; h *= 1099511628211ULL; }
    return h;
}

// Deterministic Binomial(n, p) via a seeded Bernoulli sum. Exact and O(n);
// total work across all (cell, site) equals the pass-1 ref/alt incidence count.
static inline uint32_t wobble_alt_count(uint64_t seed, uint32_t n, double p){
    if (p <= 0.0) return 0;
    if (p >= 1.0) return n;
    uint64_t st = seed ? seed : 0x123456789ABCDEFULL;
    uint32_t k = 0;
    for (uint32_t i = 0; i < n; i++){
        uint64_t r = splitmix64(st);
        double u = (double)(r >> 11) * (1.0 / 9007199254740992.0); // 53-bit [0,1)
        if (u < p) k++;
    }
    return k;
}
static inline uint32_t hardsnap_alt_count(uint32_t n, double p){
    if (p <= 0.0) return 0;
    if (p >= 1.0) return n;
    double v = (double)n * p;
    long k = (long)llround(v);
    if (k < 0) k = 0;
    if (k > (long)n) k = (long)n;
    return (uint32_t)k;
}

// ----------------------------------------------------------------------------
// Data structures
// ----------------------------------------------------------------------------
// One panel site. geno for needed samples lives in a flat pool; geno_off is the
// base offset for this site (n_needed int8 values, dosage 0/1/2 or -1 missing).
struct Site {
    int32_t  tid;
    int32_t  pos;      // 0-based reference position
    char     ref;
    char     alt;
    uint8_t  panel;    // bitmask: PANEL_BIT_INTER | PANEL_BIT_SPECIES
    uint32_t geno_off; // index into geno_pool
};

// Home identity resolved to compact genotype columns, one record per assigned
// cell. cell_seed folds library, barcode, and policy hash once so the per-site
// seed only mixes in contig/pos/panel.
struct CellRec {
    uint8_t  kind;     // 1 = singlet, 2 = fusion
    int32_t  col_i;    // compact genotype column (needed-sample index)
    int32_t  col_j;    // second column for fusion; -1 for singlet
    uint32_t ident_id; // index into ident_names
    double   fusion_r; // per-cell endogenous genome-A fraction for a fusion cell.
                       // Set from the --fusion-r-table entry for (library, barcode)
                       // when present, otherwise the global --fusion-r. Unused for
                       // singlets (kind==1). Preserving it per cell is what lets a
                       // fusion SENDER be cleaned of third-party leakage at its
                       // disagreement sites toward its real, typed r instead of the
                       // balanced 0.5, without erasing its true A/B imbalance.
    uint64_t cell_seed;
    uint64_t reads_seen;
    uint64_t reads_rewritten;
    uint64_t bases_scrubbable;
    uint64_t bases_rewritten;
};

// Per (cell, site) aggregate. Pass 1 fills n_ref/n_alt; finalize sets the flip
// budget consumed in coordinate order during pass 2. Exactly one of
// flip_ref_to_alt / flip_alt_to_ref is nonzero.
struct Agg {
    uint32_t n_ref;
    uint32_t n_alt;
    uint32_t n_other;  // covered reads whose base is neither ref nor alt (errors)
    uint32_t flip_ref_to_alt;
    uint32_t flip_alt_to_ref;
};

static inline uint64_t agg_key(uint32_t cell_idx, uint32_t site_idx){
    return ((uint64_t)cell_idx << 32) | (uint64_t)site_idx;
}

// ----------------------------------------------------------------------------
// Assignment file parse (demux 4+ column: barcode<TAB>identity<TAB>...).
// Keys are bare barcodes. Header rows (first field barcode/cell/cb) and comment
// rows are skipped, matching synth_common.load_assignments.
// ----------------------------------------------------------------------------
static void load_assignments(const string& path,
                             unordered_map<string, string>& bc2ident){
    ifstream f(path.c_str());
    if (!f) die("could not open assignments file: " + path);
    string line;
    while (getline(f, line)){
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        size_t t1 = line.find('\t');
        if (t1 == string::npos) continue;
        string bc = line.substr(0, t1);
        string rest = line.substr(t1 + 1);
        size_t t2 = rest.find('\t');
        string ident = (t2 == string::npos) ? rest : rest.substr(0, t2);
        if (bc.empty() || ident.empty()) continue;
        string bcl = bc;
        transform(bcl.begin(), bcl.end(), bcl.begin(), ::tolower);
        if (bcl == "barcode" || bcl == "cell" || bcl == "cb") continue;
        bc2ident[normalize_bc(bc)] = ident; // last write wins (refined overrides)
    }
}

// ----------------------------------------------------------------------------
// Per-cell fusion-r table parse (Task 2, B34 fusion-sender completeness).
//
// A header-labeled TSV in the same family as the --mask output: it must carry a
// 'barcode' column and a fusion-r column (one of 'fusion_r', 'r', or
// 'endogenous_r'), and may carry a 'library' column. Keys are bare barcodes
// (normalize_bc), filtered to this run's --library when a 'library' column is
// present (so one shared table can hold every library and each scrub picks its
// own rows). The value is the cell's real endogenous genome-A fraction, produced
// by the B35 allele-typing (source_bams.py resolves it for a fusion-sender
// scrub). Rows with an out-of-range r are rejected loudly so a malformed table
// cannot silently fall back to 0.5. Returns the number of usable rows loaded for
// this library.
// ----------------------------------------------------------------------------
static size_t load_fusion_r_table(const string& path, const string& library,
                                  unordered_map<string, double>& bc2r){
    ifstream f(path.c_str());
    if (!f) die("could not open --fusion-r-table: " + path);
    string header;
    if (!getline(f, header)) die("--fusion-r-table is empty: " + path);
    if (!header.empty() && header.back() == '\r') header.pop_back();
    // Split the header on tabs and resolve column indices by name.
    vector<string> cols;
    {
        size_t start = 0;
        while (true){
            size_t tab = header.find('\t', start);
            string c = (tab == string::npos) ? header.substr(start) : header.substr(start, tab - start);
            string cl = c;
            transform(cl.begin(), cl.end(), cl.begin(), ::tolower);
            cols.push_back(cl);
            if (tab == string::npos) break;
            start = tab + 1;
        }
    }
    int idx_bc = -1, idx_r = -1, idx_lib = -1;
    for (int i = 0; i < (int)cols.size(); i++){
        const string& c = cols[i];
        if (c == "barcode" || c == "cb" || c == "cell") idx_bc = i;
        else if (c == "fusion_r" || c == "r" || c == "endogenous_r") { if (idx_r < 0) idx_r = i; }
        else if (c == "library") idx_lib = i;
    }
    if (idx_bc < 0)
        die("--fusion-r-table has no 'barcode' column: " + path);
    if (idx_r < 0)
        die("--fusion-r-table has no fusion-r column (expected 'fusion_r', 'r', or 'endogenous_r'): " + path);
    size_t n_loaded = 0, n_rows = 0;
    string line;
    while (getline(f, line)){
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        // Field split.
        vector<string> fields;
        size_t start = 0;
        while (true){
            size_t tab = line.find('\t', start);
            fields.push_back(tab == string::npos ? line.substr(start) : line.substr(start, tab - start));
            if (tab == string::npos) break;
            start = tab + 1;
        }
        n_rows++;
        if ((int)fields.size() <= idx_bc || (int)fields.size() <= idx_r) continue;
        if (idx_lib >= 0 && (int)fields.size() > idx_lib){
            if (fields[idx_lib] != library) continue; // not this library's row
        }
        const string& bc = fields[idx_bc];
        if (bc.empty()) continue;
        errno = 0;
        char* endp = nullptr;
        double r = strtod(fields[idx_r].c_str(), &endp);
        if (endp == fields[idx_r].c_str() || *endp != '\0')
            die("--fusion-r-table has a non-numeric fusion-r value: '" + fields[idx_r] +
                "' for barcode " + bc);
        if (!(r >= 0.0 && r <= 1.0))
            die("--fusion-r-table has an out-of-range fusion-r value " + to_string(r) +
                " (must be in [0,1]) for barcode " + bc);
        bc2r[normalize_bc(bc)] = r; // last write wins
        n_loaded++;
    }
    (void)n_rows;
    return n_loaded;
}

// ----------------------------------------------------------------------------
// Sample-order guard: a panel VCF header sample order must equal .samples
// exactly, or every genotype column offsets to the wrong individual.
// ----------------------------------------------------------------------------
static void sample_order_guard(const string& vcf_path, const vector<string>& samples,
                               const char* label){
    htsFile* vfp = bcf_open(vcf_path.c_str(), "r");
    if (!vfp) die(string("could not open ") + label + " VCF: " + vcf_path);
    bcf_hdr_t* vhdr = bcf_hdr_read(vfp);
    if (!vhdr){ bcf_close(vfp); die(string("could not read ") + label + " VCF header: " + vcf_path); }
    int nv = bcf_hdr_nsamples(vhdr);
    if (nv != (int)samples.size()){
        char b[256];
        snprintf(b, sizeof(b),
            "sample-order guard failed (%s): VCF has %d samples, .samples has %zu",
            label, nv, samples.size());
        bcf_hdr_destroy(vhdr); bcf_close(vfp); die(b);
    }
    for (int i = 0; i < nv; i++){
        if (samples[i] != string(vhdr->samples[i])){
            char b[512];
            snprintf(b, sizeof(b),
                "sample-order guard failed (%s) at index %d: .samples='%s' VCF='%s'",
                label, i, samples[i].c_str(), vhdr->samples[i]);
            bcf_hdr_destroy(vhdr); bcf_close(vfp); die(b);
        }
    }
    bcf_hdr_destroy(vhdr);
    bcf_close(vfp);
    fprintf(stderr, "Sample-order guard passed for %s panel (%d samples).\n", label, nv);
}

// ----------------------------------------------------------------------------
// Load one panel into the shared sites/geno_pool. Keeps a site when at least one
// needed sample has a covered genotype (0/1/2), so het-home sites are retained
// (the rebalance now targets them too), not just homozygous ones.
// ----------------------------------------------------------------------------
static void load_panel(const string& vcf_path, uint8_t panel_bit,
                       const vector<int>& needed_sorted, int n_needed,
                       const unordered_map<string, int>& seq2tid,
                       int min_vq, int min_gq, bool enforce_gq,
                       vector<Site>& sites, vector<int8_t>& geno_pool,
                       uint64_t& n_pass_out, uint64_t& n_kept_out){
    htsFile* vfp = bcf_open(vcf_path.c_str(), "r");
    if (!vfp) die("could not open panel VCF: " + vcf_path);
    bcf_hdr_t* vhdr = bcf_hdr_read(vfp);
    if (!vhdr){ bcf_close(vfp); die("could not read panel VCF header: " + vcf_path); }
    bcf1_t* rec = bcf_init();
    int num_samples = bcf_hdr_nsamples(vhdr);

    int32_t* gts = NULL; int n_gts = 0;
    float*   gqs = NULL; int n_gqs = 0;
    int last_rid = -1, last_tid = -2;
    uint64_t n_read = 0, n_pass = 0, n_kept = 0;
    vector<int8_t> tmp_geno(n_needed);

    while (bcf_read(vfp, vhdr, rec) >= 0){
        n_read++;
        if (rec->n_allele != 2) continue;

        if (rec->rid != last_rid){
            last_rid = rec->rid;
            const char* chrom = bcf_hdr_id2name(vhdr, rec->rid);
            auto it = (chrom != NULL) ? seq2tid.find(chrom) : seq2tid.end();
            last_tid = (it != seq2tid.end()) ? it->second : -1;
        }
        if (last_tid < 0) continue;

        bcf_unpack(rec, BCF_UN_STR);
        const char* a0 = rec->d.allele[0];
        const char* a1 = rec->d.allele[1];
        if (strlen(a0) != 1 || strlen(a1) != 1) continue;
        char ref = a0[0], alt = a1[0];
        if ((ref!='A'&&ref!='C'&&ref!='G'&&ref!='T') ||
            (alt!='A'&&alt!='C'&&alt!='G'&&alt!='T')) continue;
        if (ref == alt) continue;
        if (rec->qual < (float)min_vq) continue;
        n_pass++;

        int nl = bcf_get_genotypes(vhdr, rec, &gts, &n_gts);
        if (nl <= 0) continue;
        int ploidy = nl / num_samples;
        if (ploidy < 2) continue;

        int num_gq = bcf_get_format_float(vhdr, rec, "GQ", &gqs, &n_gqs);

        bool any_cov = false;
        for (int c = 0; c < n_needed; c++){
            int i = needed_sorted[c];
            int8_t g = -1;
            int32_t* gp = gts + (size_t)i * ploidy;
            bool gq_ok = true;
            // min_gq_enforced makes a valid-but-low GQ sample uncovered.
            // cellbouncer_permissive uses a present genotype regardless of the
            // GQ floor, matching the effective behavior of the demux/counting
            // VCF loader, so the scrub and the scoring counter agree on which
            // sites are het / hom / covered / uncovered.
            if (enforce_gq && num_gq >= num_samples && !bcf_float_is_missing(gqs[i]) &&
                !isnan(gqs[i]) && gqs[i] < (float)min_gq){
                gq_ok = false;
            }
            if (gq_ok && !bcf_gt_is_missing(gp[0]) && gp[0] != bcf_int32_vector_end &&
                !bcf_gt_is_missing(gp[1]) && gp[1] != bcf_int32_vector_end){
                int a = (bcf_gt_allele(gp[0]) == 1) ? 1 : 0;
                int b = (bcf_gt_allele(gp[1]) == 1) ? 1 : 0;
                g = (int8_t)(a + b);
            }
            tmp_geno[c] = g;
            if (g >= 0) any_cov = true;
        }
        if (!any_cov) continue;

        Site s;
        s.tid = last_tid;
        s.pos = (int32_t)rec->pos;
        s.ref = ref;
        s.alt = alt;
        s.panel = panel_bit;
        s.geno_off = (uint32_t)geno_pool.size();
        for (int c = 0; c < n_needed; c++) geno_pool.push_back(tmp_geno[c]);
        sites.push_back(s);
        n_kept++;
        if (n_kept % 1000000ULL == 0)
            fprintf(stderr, "  loaded %llu covered sites\r", (unsigned long long)n_kept);
    }
    if (gts) free(gts);
    if (gqs) free(gqs);
    bcf_destroy(rec);
    bcf_hdr_destroy(vhdr);
    bcf_close(vfp);
    n_pass_out = n_pass;
    n_kept_out = n_kept;
    fprintf(stderr,
        "Panel %s: %llu records, %llu passed site filter, %llu kept as covered "
        "(min_vq=%d, min_gq=%d, gq=%s).\n",
        (panel_bit == PANEL_BIT_SPECIES ? "species" : "interindividual"),
        (unsigned long long)n_read, (unsigned long long)n_pass,
        (unsigned long long)n_kept, min_vq, min_gq,
        (enforce_gq ? "min_gq_enforced" : "cellbouncer_permissive"));
}

int main(int argc, char** argv){
    string bam_path, vcf_path, species_vcf_path, assign_path, samples_path, refined_path;
    string output_path, report_path, summary_path, barcodes_path, library;
    string tag_policy = "strip_md_nm";
    string gq_policy = "cellbouncer_permissive";
    string noise_model = "wobble";
    string policy_version = "genotype_scrub_bam_v2_dosage_projection_wobble";
    string policy_hash;
    string tmp_token;
    string mask_path;
    int min_vq = 50;
    int min_gq = 30;
    int n_threads = 8;
    double error_ref = 0.005;
    double error_alt = 0.005;
    // Endogenous fraction of genome A in a heterotypic fusion. Default 0.5 is the
    // balanced dosage expectation the benchmark forces (spec decision 6); the knob
    // exists only for a future fusion-sender run that supplies a known injected r.
    double fusion_r = 0.5;
    // Optional per-cell fusion-r table (Task 2, B34 fusion-sender completeness).
    // Empty by default, so a receiver/singlet-sender scrub runs the global-r
    // fallback and is byte-identical to the prior corrected scrubber. When set,
    // a fusion cell listed in the table is rebalanced at its disagreement sites
    // toward its own typed r instead of the balanced 0.5.
    string fusion_r_table_path;

    for (int i = 1; i < argc; i++){
        string a = argv[i];
        auto need = [&](string& dst){
            if (i + 1 >= argc) die("missing value after " + a);
            dst = argv[++i];
        };
        if      (a == "--bam") need(bam_path);
        else if (a == "--vcf") need(vcf_path);
        else if (a == "--species-vcf") need(species_vcf_path);
        else if (a == "--assignments") need(assign_path);
        else if (a == "--samples") need(samples_path);
        else if (a == "--refined") need(refined_path);
        else if (a == "--output") need(output_path);
        else if (a == "--report") need(report_path);
        else if (a == "--summary") need(summary_path);
        else if (a == "--barcodes") need(barcodes_path);
        else if (a == "--library") need(library);
        else if (a == "--tag-policy") need(tag_policy);
        else if (a == "--gq-policy") need(gq_policy);
        else if (a == "--noise-model") need(noise_model);
        else if (a == "--mask") need(mask_path);
        else if (a == "--tmp-token") need(tmp_token);
        else if (a == "--policy-version") need(policy_version);
        else if (a == "--policy-hash") need(policy_hash);
        else if (a == "--min-vq"){ string t; need(t); min_vq = atoi(t.c_str()); }
        else if (a == "--min-gq"){ string t; need(t); min_gq = atoi(t.c_str()); }
        else if (a == "--error-ref"){ string t; need(t); error_ref = atof(t.c_str()); }
        else if (a == "--error-alt"){ string t; need(t); error_alt = atof(t.c_str()); }
        else if (a == "--fusion-r"){ string t; need(t); fusion_r = atof(t.c_str()); }
        else if (a == "--fusion-r-table") need(fusion_r_table_path);
        else if (a == "--threads"){ string t; need(t); n_threads = atoi(t.c_str()); }
        else if (a == "--help" || a == "-h"){ usage(); return 0; }
        else die("unknown argument: " + a);
    }

    if (bam_path.empty() || vcf_path.empty() || assign_path.empty() ||
        samples_path.empty() || output_path.empty()){
        usage();
        return 1;
    }
    if (tag_policy != "strip_md_nm" && tag_policy != "keep")
        die("--tag-policy must be strip_md_nm or keep, got: " + tag_policy);
    if (gq_policy != "min_gq_enforced" && gq_policy != "cellbouncer_permissive")
        die("--gq-policy must be cellbouncer_permissive or min_gq_enforced, got: " + gq_policy);
    if (noise_model != "wobble" && noise_model != "hard_snap")
        die("--noise-model must be wobble or hard_snap, got: " + noise_model);
    if (!(fusion_r >= 0.0 && fusion_r <= 1.0))
        die("--fusion-r must be in [0,1], got: " + to_string(fusion_r));
    if (n_threads < 1) n_threads = 1;
    bool strip_tags = (tag_policy == "strip_md_nm");
    bool enforce_gq = (gq_policy == "min_gq_enforced");
    bool wobble = (noise_model == "wobble");
    if (tmp_token.empty()) tmp_token = to_string((long)getpid());

    if (library.empty()){
        size_t slash = bam_path.find_last_of('/');
        library = (slash == string::npos) ? bam_path : bam_path.substr(slash + 1);
    }

    // Per-cell fusion-r table (optional). Loaded once, filtered to this run's
    // library. When absent, fusion_r_table.empty() drives the byte-identical
    // global-r fallback.
    unordered_map<string, double> fusion_r_table;
    size_t fusion_r_table_rows = 0;
    if (!fusion_r_table_path.empty()){
        fusion_r_table_rows = load_fusion_r_table(fusion_r_table_path, library, fusion_r_table);
        fprintf(stderr, "Loaded %zu per-cell fusion-r entries for %s from %s.\n",
                fusion_r_table_rows, library.c_str(), fusion_r_table_path.c_str());
    }

    auto t_start = chrono::steady_clock::now();

    // ------------------------------------------------------------------------
    // 1. Load demux .samples (VCF sample order) and build name -> full index.
    // ------------------------------------------------------------------------
    vector<string> samples;
    {
        ifstream sf(samples_path.c_str());
        if (!sf) die("could not open --samples file: " + samples_path);
        string line;
        while (getline(sf, line)){
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (!line.empty()) samples.push_back(line);
        }
    }
    if (samples.empty()) die("no sample names parsed from: " + samples_path);
    unordered_map<string, int> name2full;
    for (int i = 0; i < (int)samples.size(); i++) name2full[samples[i]] = i;
    fprintf(stderr, "Loaded %zu sample names from %s\n",
            samples.size(), samples_path.c_str());

    // ------------------------------------------------------------------------
    // 2. Sample-order guard for each panel VCF.
    // ------------------------------------------------------------------------
    sample_order_guard(vcf_path, samples, "interindividual");
    if (!species_vcf_path.empty())
        sample_order_guard(species_vcf_path, samples, "species");

    // ------------------------------------------------------------------------
    // 3. Load assignments (demux, then refined override), then optional
    // barcode whitelist.
    // ------------------------------------------------------------------------
    unordered_map<string, string> bc2ident;
    load_assignments(assign_path, bc2ident);
    fprintf(stderr, "Loaded %zu demux assignments.\n", bc2ident.size());
    if (!refined_path.empty()){
        size_t before = bc2ident.size();
        load_assignments(refined_path, bc2ident);
        fprintf(stderr, "Applied refined override from %s (map now %zu; was %zu).\n",
                refined_path.c_str(), bc2ident.size(), before);
    }

    bool has_whitelist = !barcodes_path.empty();
    unordered_set<string> whitelist;
    if (has_whitelist){
        gzFile bf = gzopen(barcodes_path.c_str(), "rb");
        if (!bf) die("could not open --barcodes: " + barcodes_path);
        char buf[1 << 16];
        while (gzgets(bf, buf, sizeof(buf)) != Z_NULL){
            string line(buf);
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
                line.pop_back();
            if (!line.empty()) whitelist.insert(normalize_bc(line));
        }
        gzclose(bf);
        fprintf(stderr, "Loaded %zu whitelist barcodes.\n", whitelist.size());
    }

    // ------------------------------------------------------------------------
    // 4. Resolve identities -> needed sample set + compact columns, then build
    // the per-cell record vector (indexed by compact cell id). Unresolvable
    // identities (species labels, doublets with a part not in .samples, >2
    // constituents) are left out, so their reads pass through unmodified.
    // ------------------------------------------------------------------------
    auto split_plus = [](const string& s, vector<string>& out){
        out.clear();
        size_t start = 0;
        while (true){
            size_t p = s.find('+', start);
            if (p == string::npos){ out.push_back(s.substr(start)); break; }
            out.push_back(s.substr(start, p - start));
            start = p + 1;
        }
    };

    // Pass A: collect needed full sample indices.
    unordered_set<int> needed_full;
    {
        vector<string> parts;
        for (auto& kv : bc2ident){
            if (has_whitelist && whitelist.find(kv.first) == whitelist.end()) continue;
            split_plus(kv.second, parts);
            if (parts.size() < 1 || parts.size() > 2) continue;
            bool ok = true;
            int idx[2] = {-1, -1};
            for (size_t p = 0; p < parts.size(); p++){
                auto it = name2full.find(parts[p]);
                if (it == name2full.end()){ ok = false; break; }
                idx[p] = it->second;
            }
            if (!ok) continue;
            needed_full.insert(idx[0]);
            if (parts.size() == 2) needed_full.insert(idx[1]);
        }
    }
    if (needed_full.empty())
        die("no assigned barcode resolved to a genotype column; nothing to scrub. "
            "Check --assignments/--samples pairing and barcode format.");

    vector<int> needed_sorted(needed_full.begin(), needed_full.end());
    sort(needed_sorted.begin(), needed_sorted.end());
    unordered_map<int, int> full2compact;
    for (int c = 0; c < (int)needed_sorted.size(); c++)
        full2compact[needed_sorted[c]] = c;
    int n_needed = (int)needed_sorted.size();
    fprintf(stderr, "%d distinct home samples needed across assigned cells.\n", n_needed);

    // Per-cell records + barcode -> compact cell index. cell_seed folds
    // library, barcode, and policy hash once.
    uint64_t run_seed = mix64(fnv1a64(library), fnv1a64(policy_hash));
    vector<string> ident_names;
    unordered_map<string, uint32_t> ident2id;
    vector<CellRec> cells;
    unordered_map<string, uint32_t> bc2cell;
    bc2cell.reserve(bc2ident.size());
    size_t fusion_r_table_cells = 0; // fusion cells that received a per-cell r
    {
        vector<string> parts;
        for (auto& kv : bc2ident){
            const string& bc = kv.first;
            if (has_whitelist && whitelist.find(bc) == whitelist.end()) continue;
            split_plus(kv.second, parts);
            if (parts.size() < 1 || parts.size() > 2) continue;
            int idx[2] = {-1, -1};
            bool ok = true;
            for (size_t p = 0; p < parts.size(); p++){
                auto it = name2full.find(parts[p]);
                if (it == name2full.end()){ ok = false; break; }
                idx[p] = it->second;
            }
            if (!ok) continue;
            auto iid = ident2id.find(kv.second);
            uint32_t ident_id;
            if (iid == ident2id.end()){
                ident_id = (uint32_t)ident_names.size();
                ident_names.push_back(kv.second);
                ident2id[kv.second] = ident_id;
            } else {
                ident_id = iid->second;
            }
            CellRec c;
            c.kind = (parts.size() == 2) ? 2 : 1;
            c.col_i = full2compact[idx[0]];
            c.col_j = (parts.size() == 2) ? full2compact[idx[1]] : -1;
            c.ident_id = ident_id;
            // Per-cell fusion-r: for a fusion cell (kind==2) present in the table,
            // use its typed r; otherwise the global --fusion-r. Singlets keep the
            // global value (unused). The fallback path (empty table) leaves every
            // cell at the global r, so the scrub is byte-identical to the prior
            // corrected scrubber.
            c.fusion_r = fusion_r;
            if (c.kind == 2 && !fusion_r_table.empty()){
                auto rit = fusion_r_table.find(bc);
                if (rit != fusion_r_table.end()){
                    c.fusion_r = rit->second;
                    fusion_r_table_cells++;
                }
            }
            c.cell_seed = mix64(run_seed, fnv1a64(bc));
            c.reads_seen = c.reads_rewritten = c.bases_scrubbable = c.bases_rewritten = 0;
            uint32_t cell_idx = (uint32_t)cells.size();
            cells.push_back(c);
            bc2cell[bc] = cell_idx;
        }
    }
    fprintf(stderr, "%zu cells have a scrubbable home identity.\n", cells.size());

    // ------------------------------------------------------------------------
    // 5. Open the BAM header to map VCF chrom -> tid. One shared BGZF codec
    // thread pool for every handle (input passes and output).
    // ------------------------------------------------------------------------
    htsThreadPool thread_pool = {NULL, 0};
    if (n_threads > 1){
        thread_pool.pool = hts_tpool_init(n_threads);
        if (!thread_pool.pool) die("failed to init htslib thread pool");
    }
    bam_hdr_t* hdr = NULL;
    unordered_map<string, int> seq2tid;
    {
        htsFile* bam_probe = hts_open(bam_path.c_str(), "r");
        if (!bam_probe) die("could not open --bam: " + bam_path);
        if (thread_pool.pool && hts_set_opt(bam_probe, HTS_OPT_THREAD_POOL, &thread_pool) != 0)
            die("failed to attach thread pool to input BAM");
        hdr = sam_hdr_read(bam_probe);
        if (!hdr) die("could not read BAM header: " + bam_path);
        for (int t = 0; t < hdr->n_targets; t++) seq2tid[hdr->target_name[t]] = t;
        if (sam_close(bam_probe) != 0) die("failed to close probe BAM: " + bam_path);
    }
    // Per-tid contig-name hash for the per-cell-site seed (stable across BAMs
    // that share a reference; independent of tid ordering).
    vector<uint64_t> tid_name_hash(hdr->n_targets);
    for (int t = 0; t < hdr->n_targets; t++)
        tid_name_hash[t] = fnv1a64(string(hdr->target_name[t]));

    // ------------------------------------------------------------------------
    // 6. Load both panels into the union site index. Any-covered prune keeps
    // het-home sites. Species VCF shares the demux sample columns exactly.
    // ------------------------------------------------------------------------
    vector<Site> sites;
    vector<int8_t> geno_pool;
    uint64_t inter_pass = 0, inter_kept = 0, species_pass = 0, species_kept = 0;
    load_panel(vcf_path, PANEL_BIT_INTER, needed_sorted, n_needed, seq2tid,
               min_vq, min_gq, enforce_gq, sites, geno_pool, inter_pass, inter_kept);
    if (!species_vcf_path.empty())
        load_panel(species_vcf_path, PANEL_BIT_SPECIES, needed_sorted, n_needed, seq2tid,
                   min_vq, min_gq, enforce_gq, sites, geno_pool, species_pass, species_kept);

    // Sort by (tid, pos, panel) so the interindividual copy sorts first among
    // duplicates, then dedup by (tid, pos): a position assayed by both panels is
    // kept once, from the interindividual panel, matching the demux
    // dedup-species-against-demux behavior. The kept site's panel bit is OR'd with
    // the dropped duplicate's bit so a shared position is marked 'U' in the mask,
    // and the seed uses the interindividual (primary) panel id.
    sort(sites.begin(), sites.end(), [](const Site& a, const Site& b){
        if (a.tid != b.tid) return a.tid < b.tid;
        if (a.pos != b.pos) return a.pos < b.pos;
        return a.panel < b.panel;
    });
    uint64_t overlap_dropped = 0, dup_within_panel = 0;
    {
        vector<Site> dedup;
        dedup.reserve(sites.size());
        size_t i = 0, n = sites.size();
        while (i < n){
            size_t j = i + 1;
            while (j < n && sites[j].tid == sites[i].tid && sites[j].pos == sites[i].pos) j++;
            Site kept = sites[i];
            for (size_t k = i + 1; k < j; k++){
                if ((sites[k].panel & kept.panel) == 0) overlap_dropped++;
                else dup_within_panel++;
                kept.panel |= sites[k].panel; // mark shared membership
            }
            dedup.push_back(kept);
            i = j;
        }
        sites.swap(dedup);
    }
    uint64_t kept_inter = 0, kept_species = 0, kept_union = 0;
    for (const Site& s : sites){
        if (s.panel & PANEL_BIT_SPECIES) kept_species++;
        if (s.panel & PANEL_BIT_INTER)   kept_inter++;
        if ((s.panel & PANEL_BIT_INTER) && (s.panel & PANEL_BIT_SPECIES)) kept_union++;
    }
    if (overlap_dropped || dup_within_panel)
        fprintf(stderr, "Deduped union: %llu cross-panel overlaps merged, "
                "%llu within-panel duplicate positions dropped.\n",
                (unsigned long long)overlap_dropped, (unsigned long long)dup_within_panel);

    // Per-tid [lo,hi) ranges into the flat sorted sites vector.
    unordered_map<int, pair<size_t, size_t>> tid_range;
    {
        size_t i = 0, n = sites.size();
        while (i < n){
            int tid = sites[i].tid;
            size_t j = i + 1;
            while (j < n && sites[j].tid == tid) j++;
            tid_range[tid] = make_pair(i, j);
            i = j;
        }
    }
    fprintf(stderr, "Indexed %zu union sites across %zu contigs "
            "(interindividual=%llu, species=%llu).\n",
            sites.size(), tid_range.size(),
            (unsigned long long)kept_inter, (unsigned long long)kept_species);

    // ------------------------------------------------------------------------
    // 7. PASS 1: count observed ref:alt per (cell, site) over covered sites with
    // a defined home target. Reads matching the skip mask are ignored, matching
    // pass 2.
    // ------------------------------------------------------------------------
    unordered_map<uint64_t, Agg> agg;
    uint64_t reads_seen = 0, reads_with_cb = 0, reads_with_assignment = 0;
    uint64_t bases_examined = 0, covered_incidences = 0;

    {
        htsFile* bam_in = hts_open(bam_path.c_str(), "r");
        if (!bam_in) die("could not open --bam for pass 1: " + bam_path);
        if (thread_pool.pool && hts_set_opt(bam_in, HTS_OPT_THREAD_POOL, &thread_pool) != 0)
            die("failed to attach thread pool to input BAM (pass 1)");
        bam_hdr_t* hdr1 = sam_hdr_read(bam_in);
        if (!hdr1) die("could not read BAM header (pass 1): " + bam_path);
        bam1_t* rec = bam_init1();
        int ret;
        while ((ret = sam_read1(bam_in, hdr1, rec)) >= 0){
            reads_seen++;
            if ((rec->core.flag & SCRUB_SKIP_MASK) != 0) continue;
            uint8_t* cbtag = bam_aux_get(rec, "CB");
            const char* cbz = cbtag ? bam_aux2Z(cbtag) : NULL;
            if (!cbz) continue;
            reads_with_cb++;
            string bc = normalize_bc(cbz);
            if (has_whitelist && whitelist.find(bc) == whitelist.end()) continue;
            auto hit = bc2cell.find(bc);
            if (hit == bc2cell.end()) continue;
            reads_with_assignment++;
            uint32_t cell_idx = hit->second;
            const CellRec& home = cells[cell_idx];

            int tid = rec->core.tid;
            auto rit = tid_range.find(tid);
            if (rit == tid_range.end()) continue;
            size_t lo = rit->second.first, hi = rit->second.second;
            int read_start = (int)rec->core.pos;
            int read_end   = (int)bam_endpos(rec);
            size_t sidx = lower_bound(
                sites.begin() + lo, sites.begin() + hi, read_start,
                [](const Site& s, int v){ return s.pos < v; }) - sites.begin();
            if (sidx >= hi || sites[sidx].pos >= read_end) continue;

            uint8_t* seq = bam_get_seq(rec);
            uint32_t* cig = bam_get_cigar(rec);
            int ref_pos = read_start;
            int qpos = 0;
            size_t it = sidx;
            for (uint32_t ci = 0; ci < rec->core.n_cigar && it < hi; ci++){
                int op  = bam_cigar_op(cig[ci]);
                int len = bam_cigar_oplen(cig[ci]);
                if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF){
                    int bend = ref_pos + len;
                    while (it < hi && sites[it].pos < ref_pos) ++it;
                    while (it < hi && sites[it].pos < bend){
                        const Site& s = sites[it];
                        // Home genotype defined at this site?
                        int8_t gi = geno_pool[s.geno_off + home.col_i];
                        bool covered = (home.kind == 1)
                            ? (gi >= 0)
                            : (gi >= 0 && geno_pool[s.geno_off + home.col_j] >= 0);
                        if (covered){
                            int in_read_off = qpos + (s.pos - ref_pos);
                            int base_code = bam_seqi(seq, in_read_off);
                            char obs = seq_nt16_str[base_code];
                            if (obs == s.ref || obs == s.alt){
                                bases_examined++;
                                covered_incidences++;
                                uint64_t key = agg_key(cell_idx, (uint32_t)it);
                                Agg& a = agg[key];
                                if (obs == s.alt) a.n_alt++;
                                else              a.n_ref++;
                            } else if (obs=='A'||obs=='C'||obs=='G'||obs=='T'){
                                // A covered read whose base is neither ref nor alt
                                // is a sequencing error. It carries no biallelic
                                // identity and the counter ignores it, so it is not
                                // rebalanced; it is tallied for the audit mask only.
                                uint64_t key = agg_key(cell_idx, (uint32_t)it);
                                agg[key].n_other++;
                            }
                        }
                        ++it;
                    }
                    ref_pos = bend; qpos += len;
                } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP){
                    qpos += len;
                } else if (op == BAM_CDEL || op == BAM_CREF_SKIP){
                    ref_pos += len;
                }
            }
        }
        if (ret < -1) die("error reading source BAM in pass 1: " + bam_path);
        bam_destroy1(rec);
        bam_hdr_destroy(hdr1);
        if (sam_close(bam_in) != 0) die("failed to close input BAM after pass 1: " + bam_path);
    }
    fprintf(stderr,
        "Pass 1: reads_seen=%llu reads_with_cb=%llu reads_with_assignment=%llu "
        "covered_incidences=%llu distinct_cell_sites=%zu\n",
        (unsigned long long)reads_seen, (unsigned long long)reads_with_cb,
        (unsigned long long)reads_with_assignment,
        (unsigned long long)covered_incidences, agg.size());

    // ------------------------------------------------------------------------
    // 8. FINALIZE: per (cell, site), compute the target alt fraction from the
    // home dosage with error adjustment, draw the target alt count (wobble or
    // hard snap), set the minimum flip budget in one direction, accumulate the
    // category census, and (optionally) emit the per-cell-site audit mask.
    // ------------------------------------------------------------------------
    // Recover barcode strings per cell for the mask and report.
    vector<string> cell_bc(cells.size());
    for (auto& kv : bc2cell) cell_bc[kv.second] = kv.first;

    // Category census of scrubbable cell-sites: what the old homozygous-only rule
    // could and could not touch. het/disagreement are the sites it skipped.
    uint64_t cs_sng_homref = 0, cs_sng_het = 0, cs_sng_homalt = 0;
    uint64_t cs_fus_ref = 0, cs_fus_alt = 0, cs_fus_disagree = 0, cs_fus_het = 0;

    // Optional audit mask.
    bool mask_gz = false;
    gzFile mf_gz = NULL;
    FILE*  mf = NULL;
    if (!mask_path.empty()){
        mask_gz = mask_path.size() > 3 &&
                  mask_path.compare(mask_path.size() - 3, 3, ".gz") == 0;
        const char* mhdr =
            "library\tbarcode\tidentity\tpanel\tcontig\tpos\tref\talt\thome_kind\t"
            "raw_p\ttarget_frac\tn_ref\tn_alt\tn_other\tn_refalt\ttarget_alt\tdir\tflips\n";
        if (mask_gz){
            mf_gz = gzopen(mask_path.c_str(), "wb");
            if (!mf_gz) die("could not open --mask: " + mask_path);
            gzprintf(mf_gz, "%s", mhdr);
        } else {
            mf = fopen(mask_path.c_str(), "w");
            if (!mf) die("could not open --mask: " + mask_path);
            fprintf(mf, "%s", mhdr);
        }
    }

    uint64_t sites_target_ge = 0; // (cell,site) with n>=1 that got a target
    for (auto& kv : agg){
        uint32_t cell_idx = (uint32_t)(kv.first >> 32);
        uint32_t site_idx = (uint32_t)(kv.first & 0xFFFFFFFFULL);
        Agg& a = kv.second;
        uint32_t n = a.n_ref + a.n_alt;
        a.flip_ref_to_alt = 0;
        a.flip_alt_to_ref = 0;
        if (n == 0) continue;
        const CellRec& home = cells[cell_idx];
        const Site& s = sites[site_idx];
        int8_t gi = geno_pool[s.geno_off + home.col_i];
        double base_p;
        if (home.kind == 1){
            base_p = (double)gi / 2.0;                              // singlet nalt/2
            if (gi == 0) cs_sng_homref++;
            else if (gi == 1) cs_sng_het++;
            else cs_sng_homalt++;
        } else {
            int8_t gj = geno_pool[s.geno_off + home.col_j];
            double pA = (double)gi / 2.0, pB = (double)gj / 2.0;
            base_p = home.fusion_r * pA + (1.0 - home.fusion_r) * pB;  // r*pA+(1-r)*pB, r per cell
            if (gi == 0 && gj == 0) cs_fus_ref++;
            else if (gi == 2 && gj == 2) cs_fus_alt++;
            else if ((gi == 0 && gj == 2) || (gi == 2 && gj == 0)) cs_fus_disagree++;
            if (gi == 1 || gj == 1) cs_fus_het++;
        }
        double f_target = adjust_p_err(base_p, error_ref, error_alt);
        if (f_target < 0.0) f_target = 0.0;
        if (f_target > 1.0) f_target = 1.0;

        uint32_t k;
        if (wobble){
            uint64_t site_seed = mix64(mix64(tid_name_hash[s.tid], (uint64_t)(uint32_t)s.pos),
                                       (uint64_t)panel_primary_seed_id(s.panel));
            uint64_t seed = mix64(home.cell_seed, site_seed);
            k = wobble_alt_count(seed, n, f_target);
        } else {
            k = hardsnap_alt_count(n, f_target);
        }
        const char* dir = "none";
        if (k > a.n_alt){      a.flip_ref_to_alt = k - a.n_alt; dir = "ref2alt"; } // add alt reads
        else if (k < a.n_alt){ a.flip_alt_to_ref = a.n_alt - k; dir = "alt2ref"; } // remove alt reads
        sites_target_ge++;

        if (!mask_path.empty()){
            char panelc = ((s.panel & PANEL_BIT_INTER) && (s.panel & PANEL_BIT_SPECIES)) ? 'U'
                        : (s.panel & PANEL_BIT_INTER) ? 'I' : 'S';
            uint32_t flips = a.flip_ref_to_alt + a.flip_alt_to_ref;
            const char* contig = hdr->target_name[s.tid];
            const char* bcs = cell_bc[cell_idx].c_str();
            const char* idn = ident_names[home.ident_id].c_str();
            if (mask_gz){
                gzprintf(mf_gz, "%s\t%s\t%s\t%c\t%s\t%d\t%c\t%c\t%d\t%.6f\t%.6f\t%u\t%u\t%u\t%u\t%u\t%s\t%u\n",
                    library.c_str(), bcs, idn, panelc, contig, s.pos, s.ref, s.alt,
                    (int)home.kind, base_p, f_target, a.n_ref, a.n_alt, a.n_other, n, k, dir, flips);
            } else {
                fprintf(mf, "%s\t%s\t%s\t%c\t%s\t%d\t%c\t%c\t%d\t%.6f\t%.6f\t%u\t%u\t%u\t%u\t%u\t%s\t%u\n",
                    library.c_str(), bcs, idn, panelc, contig, s.pos, s.ref, s.alt,
                    (int)home.kind, base_p, f_target, a.n_ref, a.n_alt, a.n_other, n, k, dir, flips);
            }
        }
    }
    if (mf_gz){ gzclose(mf_gz); }
    if (mf){ fclose(mf); }

    // ------------------------------------------------------------------------
    // 9. PASS 2: stream again, flip the minimum reads in coordinate order, write
    // every record exactly once. Coordinate order is identical to pass 1, so the
    // flip selection is deterministic and byte-reproducible.
    // ------------------------------------------------------------------------
    string tmp_out = output_path + ".tmp." + tmp_token;
    string tmp_csi = tmp_out + ".csi";
    htsFile* bam_out = hts_open(tmp_out.c_str(), "wb");
    if (!bam_out) die("could not open temporary output BAM: " + tmp_out);
    if (thread_pool.pool && hts_set_opt(bam_out, HTS_OPT_THREAD_POOL, &thread_pool) != 0)
        die("failed to attach thread pool to output BAM");

    uint64_t reads_written = 0, reads_rewritten = 0;
    uint64_t bases_scrubbable = 0;
    uint64_t bases_rewritten_to_ref = 0, bases_rewritten_to_alt = 0;
    uint64_t md_nm_tags_stripped = 0;

    {
        htsFile* bam_in = hts_open(bam_path.c_str(), "r");
        if (!bam_in) die("could not open --bam for pass 2: " + bam_path);
        if (thread_pool.pool && hts_set_opt(bam_in, HTS_OPT_THREAD_POOL, &thread_pool) != 0)
            die("failed to attach thread pool to input BAM (pass 2)");
        bam_hdr_t* hdr2 = sam_hdr_read(bam_in);
        if (!hdr2) die("could not read BAM header (pass 2): " + bam_path);
        if (sam_hdr_write(bam_out, hdr2) < 0) die("failed to write BAM header to: " + tmp_out);
        bam1_t* rec = bam_init1();
        int ret;
        while ((ret = sam_read1(bam_in, hdr2, rec)) >= 0){
            bool rewrote_any = false;
            if ((rec->core.flag & SCRUB_SKIP_MASK) == 0){
                uint8_t* cbtag = bam_aux_get(rec, "CB");
                const char* cbz = cbtag ? bam_aux2Z(cbtag) : NULL;
                if (cbz){
                    string bc = normalize_bc(cbz);
                    if (!has_whitelist || whitelist.find(bc) != whitelist.end()){
                        auto hit = bc2cell.find(bc);
                        if (hit != bc2cell.end()){
                            uint32_t cell_idx = hit->second;
                            CellRec& home = cells[cell_idx];
                            home.reads_seen++;

                            int tid = rec->core.tid;
                            auto rit = tid_range.find(tid);
                            if (rit != tid_range.end()){
                                size_t lo = rit->second.first, hi = rit->second.second;
                                int read_start = (int)rec->core.pos;
                                int read_end   = (int)bam_endpos(rec);
                                size_t sidx = lower_bound(
                                    sites.begin() + lo, sites.begin() + hi, read_start,
                                    [](const Site& s, int v){ return s.pos < v; }) - sites.begin();
                                if (sidx < hi && sites[sidx].pos < read_end){
                                    uint8_t* seq = bam_get_seq(rec);
                                    uint32_t* cig = bam_get_cigar(rec);
                                    int ref_pos = read_start;
                                    int qpos = 0;
                                    size_t it = sidx;
                                    for (uint32_t ci = 0; ci < rec->core.n_cigar && it < hi; ci++){
                                        int op  = bam_cigar_op(cig[ci]);
                                        int len = bam_cigar_oplen(cig[ci]);
                                        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF){
                                            int bend = ref_pos + len;
                                            while (it < hi && sites[it].pos < ref_pos) ++it;
                                            while (it < hi && sites[it].pos < bend){
                                                const Site& s = sites[it];
                                                int8_t gi = geno_pool[s.geno_off + home.col_i];
                                                bool covered = (home.kind == 1)
                                                    ? (gi >= 0)
                                                    : (gi >= 0 && geno_pool[s.geno_off + home.col_j] >= 0);
                                                if (covered){
                                                    int in_read_off = qpos + (s.pos - ref_pos);
                                                    int base_code = bam_seqi(seq, in_read_off);
                                                    char obs = seq_nt16_str[base_code];
                                                    if (obs == s.ref || obs == s.alt){
                                                        bases_scrubbable++;
                                                        home.bases_scrubbable++;
                                                        auto ai = agg.find(agg_key(cell_idx, (uint32_t)it));
                                                        if (ai != agg.end()){
                                                            Agg& a = ai->second;
                                                            char newbase = 0;
                                                            if (obs == s.ref && a.flip_ref_to_alt > 0){
                                                                newbase = s.alt; a.flip_ref_to_alt--;
                                                            } else if (obs == s.alt && a.flip_alt_to_ref > 0){
                                                                newbase = s.ref; a.flip_alt_to_ref--;
                                                            }
                                                            if (newbase != 0){
                                                                uint8_t enc = seq_nt16_table[(unsigned char)newbase];
                                                                if ((in_read_off & 1) == 0)
                                                                    seq[in_read_off >> 1] =
                                                                        (seq[in_read_off >> 1] & 0x0F) | (enc << 4);
                                                                else
                                                                    seq[in_read_off >> 1] =
                                                                        (seq[in_read_off >> 1] & 0xF0) | enc;
                                                                rewrote_any = true;
                                                                home.bases_rewritten++;
                                                                if (newbase == s.ref) bases_rewritten_to_ref++;
                                                                else                  bases_rewritten_to_alt++;
                                                            }
                                                        }
                                                    }
                                                }
                                                ++it;
                                            }
                                            ref_pos = bend; qpos += len;
                                        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP){
                                            qpos += len;
                                        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP){
                                            ref_pos += len;
                                        }
                                    }
                                }
                            }
                            if (rewrote_any){
                                reads_rewritten++;
                                home.reads_rewritten++;
                                if (strip_tags){
                                    uint8_t* p;
                                    if ((p = bam_aux_get(rec, "MD")) != NULL){ bam_aux_del(rec, p); md_nm_tags_stripped++; }
                                    if ((p = bam_aux_get(rec, "NM")) != NULL){ bam_aux_del(rec, p); md_nm_tags_stripped++; }
                                }
                            }
                        }
                    }
                }
            }
            if (sam_write1(bam_out, hdr2, rec) < 0) die("failed to write record to: " + tmp_out);
            reads_written++;
        }
        if (ret < -1) die("error reading source BAM in pass 2: " + bam_path);
        bam_destroy1(rec);
        bam_hdr_destroy(hdr2);
        if (sam_close(bam_in) != 0) die("failed to close input BAM after pass 2: " + bam_path);
    }
    if (sam_close(bam_out) != 0) die("failed to close output BAM cleanly: " + tmp_out);

    if (reads_written != reads_seen){
        char b[128];
        snprintf(b, sizeof(b), "read-count invariant violated: seen=%llu written=%llu",
                 (unsigned long long)reads_seen, (unsigned long long)reads_written);
        die(b);
    }

    // Sanity: every flip budget must be fully consumed (each direction is
    // bounded by the pass-1 ref/alt count, so pass 2 always satisfies it).
    uint64_t flips_unconsumed = 0;
    for (auto& kv : agg){
        flips_unconsumed += kv.second.flip_ref_to_alt + kv.second.flip_alt_to_ref;
    }
    if (flips_unconsumed != 0){
        char b[160];
        snprintf(b, sizeof(b),
                 "flip budget not fully consumed (%llu remaining); pass 1/pass 2 read "
                 "sets disagree", (unsigned long long)flips_unconsumed);
        die(b);
    }

    // ------------------------------------------------------------------------
    // 10. Index the temp BAM (validates it is readable end to end).
    // ------------------------------------------------------------------------
    string csi_path = output_path + ".csi";
    if (sam_index_build3(tmp_out.c_str(), tmp_csi.c_str(), CSI_MIN_SHIFT, n_threads) != 0){
        remove(tmp_out.c_str());
        die("failed to build CSI index for temp BAM: " + tmp_csi);
    }

    auto t_end = chrono::steady_clock::now();
    double elapsed = chrono::duration<double>(t_end - t_start).count();
    double thru = (elapsed > 0.0) ? (double)reads_seen / elapsed : 0.0;

    // ------------------------------------------------------------------------
    // 11. Per-cell report TSV to a temp path (plain or gzip by suffix).
    // ------------------------------------------------------------------------
    string tmp_report;
    if (!report_path.empty()){
        bool gz = report_path.size() > 3 &&
                  report_path.compare(report_path.size() - 3, 3, ".gz") == 0;
        tmp_report = report_path + ".tmp." + tmp_token;
        const char* hdr_line =
            "library\tbarcode\tidentity\treads_seen\treads_rewritten\tbases_scrubbable\tbases_rewritten\n";
        // cell_bc was built in the finalize step above and is reused here.
        if (gz){
            gzFile rf = gzopen(tmp_report.c_str(), "wb");
            if (!rf) die("could not open temp --report: " + tmp_report);
            gzprintf(rf, "%s", hdr_line);
            for (uint32_t c = 0; c < cells.size(); c++){
                CellRec& h = cells[c];
                gzprintf(rf, "%s\t%s\t%s\t%llu\t%llu\t%llu\t%llu\n",
                         library.c_str(), cell_bc[c].c_str(), ident_names[h.ident_id].c_str(),
                         (unsigned long long)h.reads_seen, (unsigned long long)h.reads_rewritten,
                         (unsigned long long)h.bases_scrubbable, (unsigned long long)h.bases_rewritten);
            }
            gzclose(rf);
        } else {
            FILE* rf = fopen(tmp_report.c_str(), "w");
            if (!rf) die("could not open temp --report: " + tmp_report);
            fprintf(rf, "%s", hdr_line);
            for (uint32_t c = 0; c < cells.size(); c++){
                CellRec& h = cells[c];
                fprintf(rf, "%s\t%s\t%s\t%llu\t%llu\t%llu\t%llu\n",
                        library.c_str(), cell_bc[c].c_str(), ident_names[h.ident_id].c_str(),
                        (unsigned long long)h.reads_seen, (unsigned long long)h.reads_rewritten,
                        (unsigned long long)h.bases_scrubbable, (unsigned long long)h.bases_rewritten);
            }
            fclose(rf);
        }
    }

    // ------------------------------------------------------------------------
    // 12. Summary JSON to a temp path. Field names the Python resolver reads
    // (reads_seen/reads_written/policy_hash/*_signature/reads_rewritten/
    // bases_rewritten_to_*) are preserved; species_vcf_signature and
    // panels_union_signature are added for species-mode provenance.
    // ------------------------------------------------------------------------
    string tmp_summary;
    if (!summary_path.empty()){
        tmp_summary = summary_path + ".tmp." + tmp_token;
        FILE* jf = fopen(tmp_summary.c_str(), "w");
        if (!jf) die("could not open temp --summary: " + tmp_summary);
        string inter_sig = file_sig(vcf_path);
        string species_sig = file_sig(species_vcf_path); // "MISSING" if absent
        string union_sig = inter_sig + "|" + species_sig;
        fprintf(jf, "{\n");
        fprintf(jf, "  \"library\": \"%s\",\n", json_escape(library).c_str());
        fprintf(jf, "  \"input_bam\": \"%s\",\n", json_escape(bam_path).c_str());
        fprintf(jf, "  \"output_bam\": \"%s\",\n", json_escape(output_path).c_str());
        fprintf(jf, "  \"input_bam_signature\": \"%s\",\n", file_sig(bam_path).c_str());
        fprintf(jf, "  \"vcf_signature\": \"%s\",\n", inter_sig.c_str());
        fprintf(jf, "  \"species_vcf_signature\": \"%s\",\n", species_sig.c_str());
        fprintf(jf, "  \"panels_union_signature\": \"%s\",\n", json_escape(union_sig).c_str());
        fprintf(jf, "  \"assignments_signature\": \"%s\",\n", file_sig(assign_path).c_str());
        fprintf(jf, "  \"refined_signature\": \"%s\",\n",
                refined_path.empty() ? "MISSING" : file_sig(refined_path).c_str());
        fprintf(jf, "  \"samples_signature\": \"%s\",\n", file_sig(samples_path).c_str());
        fprintf(jf, "  \"policy_version\": \"%s\",\n", json_escape(policy_version).c_str());
        fprintf(jf, "  \"policy_hash\": \"%s\",\n", json_escape(policy_hash).c_str());
        fprintf(jf, "  \"min_vq\": %d,\n", min_vq);
        fprintf(jf, "  \"min_gq\": %d,\n", min_gq);
        fprintf(jf, "  \"tag_policy\": \"%s\",\n", json_escape(tag_policy).c_str());
        fprintf(jf, "  \"gq_policy\": \"%s\",\n", json_escape(gq_policy).c_str());
        fprintf(jf, "  \"fusion_policy\": \"dosage_projection\",\n");
        fprintf(jf, "  \"fusion_r\": %.6f,\n", fusion_r);
        fprintf(jf, "  \"fusion_r_source\": \"%s\",\n",
                fusion_r_table_path.empty() ? "global_default" : "per_cell_table");
        fprintf(jf, "  \"fusion_r_table_signature\": \"%s\",\n", file_sig(fusion_r_table_path).c_str());
        fprintf(jf, "  \"fusion_r_table_rows\": %zu,\n", fusion_r_table_rows);
        fprintf(jf, "  \"fusion_r_table_cells\": %zu,\n", fusion_r_table_cells);
        fprintf(jf, "  \"noise_model\": \"%s\",\n", json_escape(noise_model).c_str());
        fprintf(jf, "  \"error_ref\": %.6g,\n", error_ref);
        fprintf(jf, "  \"error_alt\": %.6g,\n", error_alt);
        fprintf(jf, "  \"species_panel\": %s,\n", species_vcf_path.empty() ? "false" : "true");
        fprintf(jf, "  \"panels_scrubbed\": \"%s\",\n",
                species_vcf_path.empty() ? "interindividual" : "interindividual+species");
        fprintf(jf, "  \"mask_emitted\": %s,\n", mask_path.empty() ? "false" : "true");
        fprintf(jf, "  \"needed_home_samples\": %d,\n", n_needed);
        fprintf(jf, "  \"interindividual_sites_pass\": %llu,\n", (unsigned long long)inter_pass);
        fprintf(jf, "  \"species_sites_pass\": %llu,\n", (unsigned long long)species_pass);
        fprintf(jf, "  \"interindividual_sites_kept\": %llu,\n", (unsigned long long)kept_inter);
        fprintf(jf, "  \"species_sites_kept\": %llu,\n", (unsigned long long)kept_species);
        fprintf(jf, "  \"union_overlap_sites\": %llu,\n", (unsigned long long)kept_union);
        fprintf(jf, "  \"union_sites\": %zu,\n", sites.size());
        fprintf(jf, "  \"overlap_dropped\": %llu,\n", (unsigned long long)overlap_dropped);
        fprintf(jf, "  \"cells_scrubbable\": %zu,\n", cells.size());
        fprintf(jf, "  \"census_singlet_homref\": %llu,\n", (unsigned long long)cs_sng_homref);
        fprintf(jf, "  \"census_singlet_het\": %llu,\n", (unsigned long long)cs_sng_het);
        fprintf(jf, "  \"census_singlet_homalt\": %llu,\n", (unsigned long long)cs_sng_homalt);
        fprintf(jf, "  \"census_fusion_agree_ref\": %llu,\n", (unsigned long long)cs_fus_ref);
        fprintf(jf, "  \"census_fusion_agree_alt\": %llu,\n", (unsigned long long)cs_fus_alt);
        fprintf(jf, "  \"census_fusion_disagree\": %llu,\n", (unsigned long long)cs_fus_disagree);
        fprintf(jf, "  \"census_fusion_has_het\": %llu,\n", (unsigned long long)cs_fus_het);
        fprintf(jf, "  \"reads_seen\": %llu,\n", (unsigned long long)reads_seen);
        fprintf(jf, "  \"reads_written\": %llu,\n", (unsigned long long)reads_written);
        fprintf(jf, "  \"reads_with_cb\": %llu,\n", (unsigned long long)reads_with_cb);
        fprintf(jf, "  \"reads_with_assignment\": %llu,\n", (unsigned long long)reads_with_assignment);
        fprintf(jf, "  \"reads_rewritten\": %llu,\n", (unsigned long long)reads_rewritten);
        fprintf(jf, "  \"cell_site_targets\": %llu,\n", (unsigned long long)sites_target_ge);
        fprintf(jf, "  \"bases_examined\": %llu,\n", (unsigned long long)bases_examined);
        fprintf(jf, "  \"covered_incidences\": %llu,\n", (unsigned long long)covered_incidences);
        fprintf(jf, "  \"bases_scrubbable\": %llu,\n", (unsigned long long)bases_scrubbable);
        fprintf(jf, "  \"bases_rewritten_to_ref\": %llu,\n", (unsigned long long)bases_rewritten_to_ref);
        fprintf(jf, "  \"bases_rewritten_to_alt\": %llu,\n", (unsigned long long)bases_rewritten_to_alt);
        fprintf(jf, "  \"md_nm_tags_stripped\": %llu,\n", (unsigned long long)md_nm_tags_stripped);
        fprintf(jf, "  \"elapsed_seconds\": %.3f,\n", elapsed);
        fprintf(jf, "  \"throughput_reads_per_second\": %.1f\n", thru);
        fprintf(jf, "}\n");
        fclose(jf);
    }

    // ------------------------------------------------------------------------
    // 13. Publish. All temps are written before any final is touched, and each
    // stale final is removed immediately before its rename. Order: BAM, CSI,
    // report, summary; the summary is the final ready marker. Consumers must
    // require summary + BAM + CSI together and validate policy_hash and input
    // signatures (including species_vcf_signature) before trusting the set.
    // ------------------------------------------------------------------------
    remove(csi_path.c_str());
    if (rename(tmp_out.c_str(), output_path.c_str()) != 0)
        die("failed to publish BAM " + tmp_out + " -> " + output_path);
    if (rename(tmp_csi.c_str(), csi_path.c_str()) != 0)
        die("failed to publish CSI " + tmp_csi + " -> " + csi_path);
    if (!report_path.empty()){
        remove(report_path.c_str());
        if (rename(tmp_report.c_str(), report_path.c_str()) != 0)
            die("failed to publish report " + tmp_report + " -> " + report_path);
    }
    if (!summary_path.empty()){
        remove(summary_path.c_str());
        if (rename(tmp_summary.c_str(), summary_path.c_str()) != 0)
            die("failed to publish summary " + tmp_summary + " -> " + summary_path);
    }

    bam_hdr_destroy(hdr);
    if (thread_pool.pool){ hts_tpool_destroy(thread_pool.pool); thread_pool.pool = NULL; }

    fprintf(stderr,
        "Done. reads_seen=%llu reads_written=%llu reads_rewritten=%llu "
        "bases_rewritten=%llu (to_ref=%llu to_alt=%llu) in %.1fs (%.0f reads/s).\n",
        (unsigned long long)reads_seen, (unsigned long long)reads_written,
        (unsigned long long)reads_rewritten,
        (unsigned long long)(bases_rewritten_to_ref + bases_rewritten_to_alt),
        (unsigned long long)bases_rewritten_to_ref, (unsigned long long)bases_rewritten_to_alt,
        elapsed, thru);
    fprintf(stderr, "Wrote %s and %s\n", output_path.c_str(), csi_path.c_str());
    return 0;
}
