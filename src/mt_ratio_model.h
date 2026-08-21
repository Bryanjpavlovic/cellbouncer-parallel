#ifndef CELLBOUNCER_MT_RATIO_MODEL_H
#define CELLBOUNCER_MT_RATIO_MODEL_H

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

struct MtRatioObservation {
    int site_index = -1;
    uint64_t ref = 0;
    uint64_t alt = 0;
    double parent1_alt_probability = 0.0;
    double parent2_alt_probability = 0.0;
    double ambient_alt_probability = 0.0;
    bool ratio_informative = false;
    bool ambient_anchor = false;
};

struct MtRatioFitConfig {
    bool use_beta_binomial = true;
    bool ambient_enabled = true;
    bool ambient_fixed = false;
    double fixed_ambient = 0.0;
    double overdispersion_initial = 0.02;
    double overdispersion_max = 0.25;
    bool overdispersion_fixed = false;
    double fixed_overdispersion = 0.02;
    bool overdispersion_penalized = false;
    double overdispersion_prior_mean = 0.02;
    double overdispersion_prior_strength = 0.0;
    int max_iterations = 500;
    double tolerance = 1e-8;
};

struct MtRatioFitResult {
    double parent1_weight = std::numeric_limits<double>::quiet_NaN();
    double parent2_weight = std::numeric_limits<double>::quiet_NaN();
    double parent2_fraction = std::numeric_limits<double>::quiet_NaN();
    double ambient_fraction = std::numeric_limits<double>::quiet_NaN();
    double overdispersion_rho = std::numeric_limits<double>::quiet_NaN();
    double log_likelihood = std::numeric_limits<double>::quiet_NaN();
    double objective_value = std::numeric_limits<double>::quiet_NaN();
    double ratio_se = std::numeric_limits<double>::quiet_NaN();
    double ratio_se_robust = std::numeric_limits<double>::quiet_NaN();
    double ambient_se = std::numeric_limits<double>::quiet_NaN();
    double ambient_se_robust = std::numeric_limits<double>::quiet_NaN();
    double overdispersion_se = std::numeric_limits<double>::quiet_NaN();
    double information_condition = std::numeric_limits<double>::quiet_NaN();
    double min_information_eigenvalue = std::numeric_limits<double>::quiet_NaN();
    int iterations = 0;
    bool converged = false;
    bool joint_ambient = false;
};

struct MtRatioFragmentSite {
    int site_index = -1;
    uint8_t allele = 0;
    double parent1_alt_probability = 0.0;
    double parent2_alt_probability = 0.0;
    double ambient_alt_probability = 0.0;
    bool ratio_informative = false;
    bool ambient_anchor = false;
};

struct MtRatioFragmentObservation {
    std::vector<MtRatioFragmentSite> sites;
};

struct MtRatioProfileResult {
    double parent2_ci_low = std::numeric_limits<double>::quiet_NaN();
    double parent2_ci_high = std::numeric_limits<double>::quiet_NaN();
    double parent1_ci_low = std::numeric_limits<double>::quiet_NaN();
    double parent1_ci_high = std::numeric_limits<double>::quiet_NaN();
    double profile_ci_level = 0.95;
    double single_parent_epsilon = 0.02;
    double delta_ll_parent1_only = std::numeric_limits<double>::quiet_NaN();
    double delta_ll_both = std::numeric_limits<double>::quiet_NaN();
    double delta_ll_parent2_only = std::numeric_limits<double>::quiet_NaN();
    int evaluations = 0;
    int failed_evaluations = 0;
    std::string inheritance_class = "AMBIGUOUS";
    std::string inheritance_class_reason = "PROFILE_NOT_AVAILABLE";
    std::string profile_status = "FAILED";
    std::vector<double> grid_r;
    std::vector<double> grid_delta_log_likelihood;
};

double mt_ratio_log_likelihood(
    const std::vector<MtRatioObservation>& observations,
    double parent2_fraction,
    double ambient_fraction,
    double rho,
    bool use_beta_binomial);

double mt_ratio_predicted_alt_probability(
    const MtRatioObservation& observation,
    double parent2_fraction,
    double ambient_fraction);

double mt_ratio_count_log_likelihood(
    uint64_t ref,
    uint64_t alt,
    double alt_probability,
    double rho,
    bool use_beta_binomial);

double mt_ratio_site_log_likelihood(
    const MtRatioObservation& observation,
    double parent2_fraction,
    double ambient_fraction,
    double rho,
    bool use_beta_binomial);

MtRatioFitResult fit_mt_ratio(
    const std::vector<MtRatioObservation>& observations,
    const MtRatioFitConfig& config);

MtRatioFitResult fit_mt_ratio_at_fixed_fraction(
    const std::vector<MtRatioObservation>& observations,
    const MtRatioFitConfig& config,
    double fixed_parent2_fraction,
    const MtRatioFitResult* warm_start = nullptr);

MtRatioProfileResult profile_mt_ratio(
    const std::vector<MtRatioObservation>& observations,
    const MtRatioFitConfig& config,
    const MtRatioFitResult& unconstrained_fit,
    double single_parent_epsilon,
    double grid_step,
    bool retain_grid);

double mt_fragment_ratio_log_likelihood(
    const std::vector<MtRatioFragmentObservation>& fragments,
    double parent2_fraction,
    double ambient_fraction);

MtRatioFitResult fit_mt_fragment_ratio(
    const std::vector<MtRatioFragmentObservation>& fragments,
    const MtRatioFitConfig& config);

MtRatioFitResult fit_mt_fragment_ratio_at_fixed_fraction(
    const std::vector<MtRatioFragmentObservation>& fragments,
    const MtRatioFitConfig& config,
    double fixed_parent2_fraction,
    const MtRatioFitResult* warm_start = nullptr);

MtRatioProfileResult profile_mt_fragment_ratio(
    const std::vector<MtRatioFragmentObservation>& fragments,
    const MtRatioFitConfig& config,
    const MtRatioFitResult& unconstrained_fit,
    double single_parent_epsilon,
    double grid_step,
    bool retain_grid);

#endif
