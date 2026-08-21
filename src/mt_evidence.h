#ifndef CELLBOUNCER_MT_EVIDENCE_H
#define CELLBOUNCER_MT_EVIDENCE_H

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace mt_evidence {

struct Site {
    int rid = -1;
    int64_t pos = -1; // 0-based
    char ref = 'N';
    char alt = 'N';
    std::vector<int8_t> state; // -1 unavailable, 0 ref-homoplasmic, 1 alt-homoplasmic
};

struct Panel {
    std::string chrom;
    std::vector<std::string> samples;
    std::unordered_map<std::string,int> sample_index;
    std::vector<Site> sites;
};

struct AlleleCount {
    uint64_t ref = 0;
    uint64_t alt = 0;
    uint64_t total() const { return ref + alt; }
};

using CellCounts = std::unordered_map<std::string, std::unordered_map<int, AlleleCount>>;

struct CountStats {
    uint64_t reads_seen = 0;
    uint64_t reads_accepted = 0;
    uint64_t observations = 0;
    uint64_t molecules = 0;
    uint64_t conflicting_molecules = 0;
};

struct Score {
    double log_likelihood = 0.0;
    uint64_t molecules = 0;
    int sites = 0;
    bool scoreable = false;
};

std::string canonical_pair_key(const std::string& a, const std::string& b);
std::vector<std::string> split_genotype(const std::string& genotype);

Panel load_panel(const std::string& vcf_path,
                 const std::string& mito_chrom,
                 int min_depth,
                 double homoplasmy_af,
                 const std::unordered_set<int64_t>* allowed_positions = nullptr);

CellCounts count_bam_once(const std::string& bam_path,
                          const Panel& panel,
                          const std::unordered_set<std::string>& target_barcodes,
                          const std::string& barcode_tag,
                          const std::string& umi_tag,
                          int min_mapq,
                          int min_baseq,
                          int threads,
                          CountStats& stats);

double binomial_log_pmf(uint64_t alt, uint64_t ref, double alt_probability);
double beta_binomial_log_pmf(uint64_t alt, uint64_t ref, double alt_probability, double rho);

Score score_genotype(const Panel& panel,
                     const std::unordered_map<int,AlleleCount>& counts,
                     const std::string& genotype,
                     double error_rate,
                     double overdispersion_rho = 0.0);

double fit_overdispersion_rho(const Panel& panel,
                              const std::unordered_map<int,AlleleCount>& counts,
                              const std::string& genotype,
                              double error_rate,
                              double rho_max = 0.25);

double fusion_parent2_fraction(const Panel& panel,
                               const std::unordered_map<int,AlleleCount>& counts,
                               const std::string& parent1,
                               const std::string& parent2,
                               uint64_t* informative_molecules = nullptr,
                               int* informative_sites = nullptr);

} // namespace mt_evidence

#endif
