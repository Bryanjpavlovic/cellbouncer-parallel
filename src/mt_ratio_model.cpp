#include "mt_ratio_model.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace {

constexpr double kTiny = 1e-12;
constexpr double kRhoLower = 1e-8;
constexpr double kProfileDrop95 = 1.920729410347062;

inline double clamp_probability(double value) {
    return std::max(kTiny, std::min(1.0 - kTiny, value));
}

inline double predicted_alt_probability(const MtRatioObservation& obs,
                                        double parent2_fraction,
                                        double ambient_fraction) {
    const double r = std::max(0.0, std::min(1.0, parent2_fraction));
    const double c = std::max(0.0, std::min(0.99, ambient_fraction));
    const double parental = (1.0 - r) * obs.parent1_alt_probability +
                            r * obs.parent2_alt_probability;
    return clamp_probability((1.0 - c) * parental + c * obs.ambient_alt_probability);
}

double log_combinatorial_term(uint64_t alt, uint64_t ref) {
    const double n = static_cast<double>(alt + ref);
    const double k = static_cast<double>(alt);
    return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) -
           std::lgamma(n - k + 1.0);
}

double binomial_log_pmf_cached(uint64_t alt, uint64_t ref, double q,
                               double log_combinatorial) {
    const double n = static_cast<double>(alt + ref);
    const double k = static_cast<double>(alt);
    q = clamp_probability(q);
    return log_combinatorial + k * std::log(q) +
           (n - k) * std::log1p(-q);
}

double binomial_log_pmf(uint64_t alt, uint64_t ref, double q) {
    return binomial_log_pmf_cached(alt, ref, q, log_combinatorial_term(alt, ref));
}

double beta_binomial_log_pmf_cached(uint64_t alt, uint64_t ref, double q, double rho,
                                    double log_combinatorial) {
    if (rho <= 1e-9 || alt + ref <= 1) {
        return binomial_log_pmf_cached(alt, ref, q, log_combinatorial);
    }
    q = clamp_probability(q);
    rho = std::max(1e-9, std::min(0.999999, rho));
    const double concentration = 1.0 / rho - 1.0;
    const double alpha = std::max(kTiny, q * concentration);
    const double beta = std::max(kTiny, (1.0 - q) * concentration);
    const double n = static_cast<double>(alt + ref);
    const double k = static_cast<double>(alt);
    return log_combinatorial +
           std::lgamma(k + alpha) + std::lgamma(n - k + beta) -
           std::lgamma(n + alpha + beta) -
           std::lgamma(alpha) - std::lgamma(beta) +
           std::lgamma(alpha + beta);
}

double beta_binomial_log_pmf(uint64_t alt, uint64_t ref, double q, double rho) {
    return beta_binomial_log_pmf_cached(
        alt, ref, q, rho, log_combinatorial_term(alt, ref));
}

double site_log_likelihood(const MtRatioObservation& obs,
                           double parent2_fraction,
                           double ambient_fraction,
                           double rho,
                           bool use_beta_binomial) {
    const double q = predicted_alt_probability(obs, parent2_fraction, ambient_fraction);
    return use_beta_binomial
        ? beta_binomial_log_pmf(obs.alt, obs.ref, q, rho)
        : binomial_log_pmf(obs.alt, obs.ref, q);
}

std::vector<double> prepare_log_combinatorial_terms(
        const std::vector<MtRatioObservation>& observations) {
    std::vector<double> terms;
    terms.reserve(observations.size());
    for (const MtRatioObservation& obs : observations) {
        terms.push_back(log_combinatorial_term(obs.alt, obs.ref));
    }
    return terms;
}

double mt_ratio_log_likelihood_cached(
        const std::vector<MtRatioObservation>& observations,
        const std::vector<double>& log_combinatorial_terms,
        double parent2_fraction,
        double ambient_fraction,
        double rho,
        bool use_beta_binomial) {
    double value = 0.0;
    for (size_t i = 0; i < observations.size(); ++i) {
        const MtRatioObservation& obs = observations[i];
        const double q = predicted_alt_probability(obs, parent2_fraction, ambient_fraction);
        value += use_beta_binomial
            ? beta_binomial_log_pmf_cached(obs.alt, obs.ref, q, rho,
                                           log_combinatorial_terms[i])
            : binomial_log_pmf_cached(obs.alt, obs.ref, q,
                                      log_combinatorial_terms[i]);
    }
    return value;
}

double rho_penalty(const MtRatioFitConfig& config, double rho) {
    if (!config.use_beta_binomial || !config.overdispersion_penalized ||
        config.overdispersion_prior_strength <= 0.0) return 0.0;
    const double scale = std::max(config.overdispersion_max, 1e-6);
    const double z = (rho - config.overdispersion_prior_mean) / scale;
    return 0.5 * config.overdispersion_prior_strength * z * z;
}

double ratio_fit_objective(
        const std::vector<MtRatioObservation>& observations,
        const MtRatioFitConfig& config,
        double parent2_fraction,
        double ambient_fraction,
        double rho,
        const std::vector<double>* log_combinatorial_terms = nullptr) {
    const double ll = log_combinatorial_terms
        ? mt_ratio_log_likelihood_cached(observations, *log_combinatorial_terms,
                                         parent2_fraction, ambient_fraction, rho,
                                         config.use_beta_binomial)
        : mt_ratio_log_likelihood(observations, parent2_fraction, ambient_fraction, rho,
                                  config.use_beta_binomial);
    return ll - rho_penalty(config, rho);
}

template <typename Function>
double golden_maximize(Function function,
                       double lower,
                       double upper,
                       int max_iterations,
                       double tolerance,
                       int& iterations) {
    if (!(upper > lower)) {
        iterations = 0;
        return lower;
    }
    const double phi = (std::sqrt(5.0) - 1.0) / 2.0;
    double left = lower;
    double right = upper;
    double x1 = right - phi * (right - left);
    double x2 = left + phi * (right - left);
    double f1 = function(x1);
    double f2 = function(x2);
    iterations = 0;
    while (iterations < max_iterations && right - left > tolerance) {
        if (f1 < f2) {
            left = x1;
            x1 = x2;
            f1 = f2;
            x2 = left + phi * (right - left);
            f2 = function(x2);
        } else {
            right = x2;
            x2 = x1;
            f2 = f1;
            x1 = right - phi * (right - left);
            f1 = function(x1);
        }
        ++iterations;
    }
    const double interior = 0.5 * (left + right);
    const double f_interior = function(interior);
    const double f_lower = function(lower);
    const double f_upper = function(upper);
    if (f_lower >= f_interior && f_lower >= f_upper) return lower;
    if (f_upper >= f_interior && f_upper >= f_lower) return upper;
    return interior;
}

struct FitStart {
    double ratio;
    double ambient;
    double rho;
};

std::vector<FitStart> deterministic_starts(const MtRatioFitConfig& config) {
    const double fixed_c = !config.ambient_enabled ? 0.0 : config.fixed_ambient;
    const bool ambient_free = config.ambient_enabled && !config.ambient_fixed;
    std::vector<FitStart> starts;
    starts.push_back({0.50, ambient_free ? 0.05 : fixed_c, config.overdispersion_initial});
    starts.push_back({0.20, ambient_free ? 0.10 : fixed_c, config.overdispersion_initial});
    starts.push_back({0.80, ambient_free ? 0.10 : fixed_c, config.overdispersion_initial});
    starts.push_back({0.50, ambient_free ? 0.30 : fixed_c,
                      std::min(0.10, config.overdispersion_max)});
    starts.push_back({0.10, ambient_free ? 0.50 : fixed_c,
                      std::min(0.05, config.overdispersion_max)});
    starts.push_back({0.90, ambient_free ? 0.50 : fixed_c,
                      std::min(0.05, config.overdispersion_max)});
    return starts;
}

MtRatioFitResult fit_internal_legacy(const std::vector<MtRatioObservation>& observations,
                                     const MtRatioFitConfig& config,
                                     const std::optional<double>& fixed_r,
                                     const MtRatioFitResult* warm_start,
                                     bool continuation_only = false) {
    const bool use_beta = config.use_beta_binomial;
    const std::vector<double> log_combinatorial_terms =
        prepare_log_combinatorial_terms(observations);
    const bool ambient_free = config.ambient_enabled && !config.ambient_fixed;
    const double fixed_c = !config.ambient_enabled ? 0.0 :
                           (config.ambient_fixed ? config.fixed_ambient : 0.0);
    std::vector<FitStart> starts;
    if (warm_start && std::isfinite(warm_start->parent2_fraction) &&
        std::isfinite(warm_start->ambient_fraction)) {
        starts.push_back({fixed_r ? *fixed_r : warm_start->parent2_fraction,
                          ambient_free ? warm_start->ambient_fraction : fixed_c,
                          use_beta && std::isfinite(warm_start->overdispersion_rho)
                              ? warm_start->overdispersion_rho
                              : config.overdispersion_initial});
    }
    if (!continuation_only) {
        const std::vector<FitStart> defaults = deterministic_starts(config);
        starts.insert(starts.end(), defaults.begin(), defaults.end());
    }
    MtRatioFitResult best;
    best.joint_ambient = ambient_free;
    double best_ll = -std::numeric_limits<double>::infinity();
    bool found_converged = false;
    const bool prefer_converged = fixed_r.has_value();
    for (const FitStart& start : starts) {
        double r = fixed_r ? std::max(0.0, std::min(1.0, *fixed_r))
                           : std::max(0.0, std::min(1.0, start.ratio));
        double c = ambient_free ? std::max(0.0, std::min(0.99, start.ambient)) : fixed_c;
        double rho = use_beta
            ? std::max(kRhoLower, std::min(config.overdispersion_max, start.rho))
            : 0.0;
        bool converged = false;
        int iteration = 0;
        double previous_ll = mt_ratio_log_likelihood_cached(observations, log_combinatorial_terms, r, c, rho, use_beta);
        if (!std::isfinite(previous_ll)) continue;
        for (; iteration < config.max_iterations; ++iteration) {
            const double old_r = r, old_c = c, old_rho = rho;
            int inner = 0;
            if (!fixed_r) {
                r = golden_maximize([&](double candidate) {
                    return mt_ratio_log_likelihood_cached(observations, log_combinatorial_terms, candidate, c, rho, use_beta);
                }, 0.0, 1.0, 100, config.tolerance, inner);
            }
            if (ambient_free) {
                c = golden_maximize([&](double candidate) {
                    return mt_ratio_log_likelihood_cached(observations, log_combinatorial_terms, r, candidate, rho, use_beta);
                }, 0.0, 0.99, 100, config.tolerance, inner);
            }
            if (use_beta) {
                rho = golden_maximize([&](double candidate) {
                    return mt_ratio_log_likelihood_cached(observations, log_combinatorial_terms, r, c, candidate, true);
                }, kRhoLower, config.overdispersion_max, 100, config.tolerance, inner);
                const double boundary_snap = std::max(10.0 * config.tolerance, 1e-7);
                if (rho <= kRhoLower + boundary_snap) rho = kRhoLower;
                if (rho >= config.overdispersion_max - boundary_snap) rho = config.overdispersion_max;
            }
            const double current_ll = mt_ratio_log_likelihood_cached(observations, log_combinatorial_terms, r, c, rho, use_beta);
            if (!std::isfinite(current_ll)) break;
            const double parameter_delta = std::max({std::fabs(r-old_r),std::fabs(c-old_c),std::fabs(rho-old_rho)});
            const double relative_ll_delta = std::fabs(current_ll-previous_ll)/(1.0+std::fabs(previous_ll));
            previous_ll=current_ll;
            const double parameter_tolerance=std::max(10.0*config.tolerance,1e-7);
            const double likelihood_tolerance=std::max(config.tolerance,1e-10);
            if (parameter_delta<=parameter_tolerance && relative_ll_delta<=likelihood_tolerance) {
                converged=true; ++iteration; break;
            }
        }
        const double ll=mt_ratio_log_likelihood_cached(observations, log_combinatorial_terms, r, c, rho, use_beta);
        if (!std::isfinite(ll)) continue;
        const bool take = prefer_converged
            ? ((converged && !found_converged) || (converged == found_converged && ll > best_ll))
            : (ll > best_ll);
        if (take) {
            found_converged=found_converged||converged; best_ll=ll;
            best.parent2_fraction=r; best.ambient_fraction=c;
            best.parent1_weight=(1.0-c)*(1.0-r); best.parent2_weight=(1.0-c)*r;
            best.overdispersion_rho=use_beta?rho:0.0; best.log_likelihood=ll;
            best.objective_value=ll; best.iterations=iteration; best.converged=converged;
            best.joint_ambient=ambient_free;
        }
    }
    return best;
}

MtRatioFitResult fit_internal(const std::vector<MtRatioObservation>& observations,
                              const MtRatioFitConfig& config,
                              const std::optional<double>& fixed_r,
                              const MtRatioFitResult* warm_start,
                              bool continuation_only = false) {
    if (!config.overdispersion_fixed && !config.overdispersion_penalized) {
        return fit_internal_legacy(
            observations, config, fixed_r, warm_start, continuation_only);
    }
    const std::vector<double> log_combinatorial_terms =
        prepare_log_combinatorial_terms(observations);
    const bool use_beta = config.use_beta_binomial;
    const bool rho_fixed = use_beta && config.overdispersion_fixed;
    const bool ambient_free = config.ambient_enabled && !config.ambient_fixed;
    const double fixed_c = !config.ambient_enabled ? 0.0 :
                           (config.ambient_fixed ? config.fixed_ambient : 0.0);

    std::vector<FitStart> starts;
    if (warm_start && std::isfinite(warm_start->parent2_fraction) &&
        std::isfinite(warm_start->ambient_fraction)) {
        starts.push_back({fixed_r ? *fixed_r : warm_start->parent2_fraction,
                          ambient_free ? warm_start->ambient_fraction : fixed_c,
                          use_beta && std::isfinite(warm_start->overdispersion_rho)
                              ? warm_start->overdispersion_rho
                              : config.overdispersion_initial});
    }
    if (!continuation_only) {
        const std::vector<FitStart> defaults = deterministic_starts(config);
        starts.insert(starts.end(), defaults.begin(), defaults.end());
    }

    MtRatioFitResult best;
    best.joint_ambient = ambient_free;
    double best_objective = -std::numeric_limits<double>::infinity();
    bool found_converged = false;
    const bool prefer_converged = fixed_r.has_value();

    for (const FitStart& start : starts) {
        double r = fixed_r ? std::max(0.0, std::min(1.0, *fixed_r))
                           : std::max(0.0, std::min(1.0, start.ratio));
        double c = ambient_free ? std::max(0.0, std::min(0.99, start.ambient)) : fixed_c;
        double rho = use_beta
            ? (rho_fixed
                ? std::max(kRhoLower, std::min(config.overdispersion_max, config.fixed_overdispersion))
                : std::max(kRhoLower, std::min(config.overdispersion_max, start.rho)))
            : 0.0;
        bool converged = false;
        int iteration = 0;
        double previous_objective = ratio_fit_objective(observations, config, r, c, rho, &log_combinatorial_terms);
        if (!std::isfinite(previous_objective)) continue;

        for (; iteration < config.max_iterations; ++iteration) {
            const double old_r = r;
            const double old_c = c;
            const double old_rho = rho;
            int inner = 0;
            if (!fixed_r) {
                r = golden_maximize(
                    [&](double candidate) {
                        return ratio_fit_objective(observations, config, candidate, c, rho, &log_combinatorial_terms);
                    }, 0.0, 1.0, 100, config.tolerance, inner);
            }

            if (ambient_free) {
                c = golden_maximize(
                    [&](double candidate) {
                        return ratio_fit_objective(observations, config, r, candidate, rho, &log_combinatorial_terms);
                    }, 0.0, 0.99, 100, config.tolerance, inner);
            }

            if (use_beta && !rho_fixed) {
                rho = golden_maximize(
                    [&](double candidate) {
                        return ratio_fit_objective(observations, config, r, c, candidate, &log_combinatorial_terms);
                    }, kRhoLower, config.overdispersion_max, 100, config.tolerance, inner);
                const double boundary_snap = std::max(10.0 * config.tolerance, 1e-7);
                if (rho <= kRhoLower + boundary_snap) rho = kRhoLower;
                if (rho >= config.overdispersion_max - boundary_snap) {
                    rho = config.overdispersion_max;
                }
            }

            const double current_objective = ratio_fit_objective(observations, config, r, c, rho, &log_combinatorial_terms);
            if (!std::isfinite(current_objective)) break;
            const double parameter_delta = std::max({std::fabs(r - old_r),
                                                     std::fabs(c - old_c),
                                                     std::fabs(rho - old_rho)});
            const double relative_ll_delta =
                std::fabs(current_objective - previous_objective) /
                (1.0 + std::fabs(previous_objective));
            previous_objective = current_objective;
            const double parameter_tolerance = std::max(10.0 * config.tolerance, 1e-7);
            const double likelihood_tolerance = std::max(config.tolerance, 1e-10);
            if (parameter_delta <= parameter_tolerance &&
                relative_ll_delta <= likelihood_tolerance) {
                converged = true;
                ++iteration;
                break;
            }
        }

        const double ll = mt_ratio_log_likelihood_cached(observations, log_combinatorial_terms, r, c, rho, use_beta);
        const double objective = ratio_fit_objective(observations, config, r, c, rho, &log_combinatorial_terms);
        if (!std::isfinite(ll) || !std::isfinite(objective)) continue;
        // The unconstrained path deliberately preserves the historical behavior:
        // select the highest likelihood regardless of its convergence flag.  For
        // constrained profile points, prefer a converged retry when one exists.
        const bool take = prefer_converged
            ? ((converged && !found_converged) ||
               (converged == found_converged && objective > best_objective))
            : (objective > best_objective);
        if (take) {
            found_converged = found_converged || converged;
            best_objective = objective;
            best.parent2_fraction = r;
            best.ambient_fraction = c;
            best.parent1_weight = (1.0 - c) * (1.0 - r);
            best.parent2_weight = (1.0 - c) * r;
            best.overdispersion_rho = use_beta ? rho : 0.0;
            best.log_likelihood = ll;
            best.objective_value = objective;
            best.iterations = iteration;
            best.converged = converged;
            best.joint_ambient = ambient_free;
        }
    }
    return best;
}

double finite_step(double value, double lower, double upper) {
    const double distance = std::min(value - lower, upper - value);
    if (distance <= 1e-7) return 0.0;
    return std::min(1e-4, 0.25 * distance);
}

void compute_uncertainty(const std::vector<MtRatioObservation>& observations,
                         const MtRatioFitConfig& config,
                         MtRatioFitResult& fit) {
    const bool use_beta = config.use_beta_binomial;
    const bool ambient_parameter_free = config.ambient_enabled && !config.ambient_fixed;
    const double r = fit.parent2_fraction;
    const double c = fit.ambient_fraction;
    const double rho = use_beta ? fit.overdispersion_rho : 0.0;
    const double hr = finite_step(r, 0.0, 1.0);
    const double hc = ambient_parameter_free ? finite_step(c, 0.0, 0.99) : 0.0;
    if (hr <= 0.0) return;

    const auto total_ll = [&](double rr, double cc) {
        return mt_ratio_log_likelihood(observations, rr, cc, rho, use_beta);
    };
    const double center = total_ll(r, c);
    const double i_rr = -(total_ll(r + hr, c) - 2.0 * center +
                          total_ll(r - hr, c)) / (hr * hr);

    if (!ambient_parameter_free) {
        if (i_rr > kTiny && std::isfinite(i_rr)) {
            fit.ratio_se = std::sqrt(1.0 / i_rr);
            double b_rr = 0.0;
            for (const MtRatioObservation& obs : observations) {
                const double score = (site_log_likelihood(obs, r + hr, c, rho, use_beta) -
                                      site_log_likelihood(obs, r - hr, c, rho, use_beta)) /
                                     (2.0 * hr);
                b_rr += score * score;
            }
            fit.ratio_se_robust = std::sqrt(std::max(0.0, b_rr / (i_rr * i_rr)));
            fit.information_condition = 1.0;
            fit.min_information_eigenvalue = i_rr;
        }
    } else if (hc > 0.0) {
        const double i_cc = -(total_ll(r, c + hc) - 2.0 * center +
                              total_ll(r, c - hc)) / (hc * hc);
        const double i_rc = -(total_ll(r + hr, c + hc) - total_ll(r + hr, c - hc) -
                              total_ll(r - hr, c + hc) + total_ll(r - hr, c - hc)) /
                            (4.0 * hr * hc);
        const double trace = i_rr + i_cc;
        const double discriminant = std::sqrt(std::max(0.0,
            (i_rr - i_cc) * (i_rr - i_cc) + 4.0 * i_rc * i_rc));
        const double lambda_max = 0.5 * (trace + discriminant);
        const double lambda_min = 0.5 * (trace - discriminant);
        fit.min_information_eigenvalue = lambda_min;
        fit.information_condition = lambda_max / std::max(lambda_min, kTiny);
        const double determinant = i_rr * i_cc - i_rc * i_rc;
        if (determinant > kTiny && i_rr > 0.0 && i_cc > 0.0 &&
            std::isfinite(determinant)) {
            const double inv_rr = i_cc / determinant;
            const double inv_rc = -i_rc / determinant;
            const double inv_cc = i_rr / determinant;
            fit.ratio_se = std::sqrt(std::max(0.0, inv_rr));
            fit.ambient_se = std::sqrt(std::max(0.0, inv_cc));

            double b_rr = 0.0;
            double b_rc = 0.0;
            double b_cc = 0.0;
            for (const MtRatioObservation& obs : observations) {
                const double score_r =
                    (site_log_likelihood(obs, r + hr, c, rho, use_beta) -
                     site_log_likelihood(obs, r - hr, c, rho, use_beta)) / (2.0 * hr);
                const double score_c =
                    (site_log_likelihood(obs, r, c + hc, rho, use_beta) -
                     site_log_likelihood(obs, r, c - hc, rho, use_beta)) / (2.0 * hc);
                b_rr += score_r * score_r;
                b_rc += score_r * score_c;
                b_cc += score_c * score_c;
            }
            const double v_rr = inv_rr * (b_rr * inv_rr + b_rc * inv_rc) +
                                inv_rc * (b_rc * inv_rr + b_cc * inv_rc);
            const double v_cc = inv_rc * (b_rr * inv_rc + b_rc * inv_cc) +
                                inv_cc * (b_rc * inv_rc + b_cc * inv_cc);
            fit.ratio_se_robust = std::sqrt(std::max(0.0, v_rr));
            fit.ambient_se_robust = std::sqrt(std::max(0.0, v_cc));
        }
    }

    if (use_beta && !config.overdispersion_fixed) {
        const double h_rho = finite_step(rho, kRhoLower, config.overdispersion_max);
        if (h_rho > 0.0) {
            const double ll_plus = mt_ratio_log_likelihood(observations, r, c, rho + h_rho, true);
            const double ll_minus = mt_ratio_log_likelihood(observations, r, c, rho - h_rho, true);
            const double info = -(ll_plus - 2.0 * center + ll_minus) /
                                (h_rho * h_rho);
            if (info > kTiny && std::isfinite(info)) {
                fit.overdispersion_se = std::sqrt(1.0 / info);
            }
        }
    }
}

bool profile_fit_valid(const MtRatioFitResult& fit) {
    return fit.converged && std::isfinite(fit.log_likelihood) &&
           std::isfinite(fit.objective_value) &&
           std::isfinite(fit.parent2_fraction) && std::isfinite(fit.ambient_fraction);
}

struct ProfilePoint {
    double r = std::numeric_limits<double>::quiet_NaN();
    double delta = std::numeric_limits<double>::quiet_NaN();
    MtRatioFitResult fit;
    bool ok = false;
};

ProfilePoint evaluate_profile_point(const std::vector<MtRatioObservation>& observations,
                                    const MtRatioFitConfig& config,
                                    const MtRatioFitResult& unconstrained_fit,
                                    double r,
                                    const MtRatioFitResult* warm_start) {
    ProfilePoint point;
    point.r = r;
    point.fit = fit_mt_ratio_at_fixed_fraction(observations, config, r, warm_start);
    point.ok = profile_fit_valid(point.fit);
    if (point.ok) {
        point.delta = point.fit.objective_value - unconstrained_fit.objective_value;
        // A fixed-r fit cannot legitimately improve on the unconstrained optimum;
        // absorb small optimization noise and keep the profile relative scale valid.
        if (point.delta > 0.0) point.delta = 0.0;
    }
    return point;
}

double refine_crossing(const std::vector<MtRatioObservation>& observations,
                       const MtRatioFitConfig& config,
                       const MtRatioFitResult& unconstrained_fit,
                       double inside_r,
                       double outside_r,
                       const MtRatioFitResult* warm_start,
                       int& evaluations,
                       int& failed_evaluations,
                       bool& failed) {
    double inside = inside_r;
    double outside = outside_r;
    MtRatioFitResult warm = warm_start ? *warm_start : unconstrained_fit;
    const double target = -kProfileDrop95;
    for (int iter = 0; iter < 80; ++iter) {
        const double mid = 0.5 * (inside + outside);
        ProfilePoint p = evaluate_profile_point(observations, config, unconstrained_fit, mid, &warm);
        ++evaluations;
        if (!p.ok) {
            ++failed_evaluations;
            failed = true;
            return inside;
        }
        warm = p.fit;
        if (p.delta >= target) inside = mid;
        else outside = mid;
        if (std::fabs(outside - inside) <= 1e-5 ||
            std::fabs(p.delta - target) <= std::max(config.tolerance, 1e-8)) {
            break;
        }
    }
    return 0.5 * (inside + outside);
}

}  // namespace

double mt_ratio_log_likelihood(const std::vector<MtRatioObservation>& observations,
                               double parent2_fraction,
                               double ambient_fraction,
                               double rho,
                               bool use_beta_binomial) {
    double value = 0.0;
    for (const MtRatioObservation& obs : observations) {
        value += site_log_likelihood(obs, parent2_fraction, ambient_fraction,
                                     rho, use_beta_binomial);
    }
    return value;
}

double mt_ratio_predicted_alt_probability(const MtRatioObservation& observation,
                                          double parent2_fraction,
                                          double ambient_fraction) {
    return predicted_alt_probability(observation, parent2_fraction, ambient_fraction);
}

double mt_ratio_count_log_likelihood(uint64_t ref, uint64_t alt,
                                     double alt_probability, double rho,
                                     bool use_beta_binomial) {
    return use_beta_binomial
        ? beta_binomial_log_pmf(alt, ref, alt_probability, rho)
        : binomial_log_pmf(alt, ref, alt_probability);
}

double mt_ratio_site_log_likelihood(const MtRatioObservation& observation,
                                    double parent2_fraction,
                                    double ambient_fraction,
                                    double rho,
                                    bool use_beta_binomial) {
    return site_log_likelihood(observation, parent2_fraction, ambient_fraction, rho,
                               use_beta_binomial);
}

MtRatioFitResult fit_mt_ratio(const std::vector<MtRatioObservation>& observations,
                              const MtRatioFitConfig& config) {
    MtRatioFitResult fit = fit_internal(observations, config, std::nullopt, nullptr);
    if (std::isfinite(fit.log_likelihood)) compute_uncertainty(observations, config, fit);
    return fit;
}

MtRatioFitResult fit_mt_ratio_at_fixed_fraction(
    const std::vector<MtRatioObservation>& observations,
    const MtRatioFitConfig& config,
    double fixed_parent2_fraction,
    const MtRatioFitResult* warm_start) {
    const double r = std::max(0.0, std::min(1.0, fixed_parent2_fraction));
    const bool use_beta = config.use_beta_binomial;
    const bool ambient_free = config.ambient_enabled && !config.ambient_fixed;
    const bool rho_free = use_beta && !config.overdispersion_fixed;
    const std::vector<double> log_combinatorial_terms =
        prepare_log_combinatorial_terms(observations);

    // Once r is fixed, zero- and one-dimensional nuisance problems do not
    // benefit from the generic multi-start coordinate optimizer.  Solve the
    // remaining bounded dimension exactly with the same objective, bounds,
    // tolerance, and boundary handling used by the generic path.
    if (static_cast<int>(ambient_free) + static_cast<int>(rho_free) <= 1) {
        const double fixed_c = !config.ambient_enabled ? 0.0 :
                               (config.ambient_fixed ? config.fixed_ambient : 0.0);
        double c = fixed_c;
        double rho = use_beta
            ? (config.overdispersion_fixed
                ? std::max(kRhoLower, std::min(config.overdispersion_max,
                                               config.fixed_overdispersion))
                : std::max(kRhoLower, std::min(config.overdispersion_max,
                                               config.overdispersion_initial)))
            : 0.0;
        int inner_iterations = 0;

        if (ambient_free) {
            c = golden_maximize(
                [&](double candidate) {
                    return ratio_fit_objective(observations, config, r, candidate, rho, &log_combinatorial_terms);
                }, 0.0, 0.99, 100, config.tolerance, inner_iterations);
        } else if (rho_free) {
            rho = golden_maximize(
                [&](double candidate) {
                    return ratio_fit_objective(observations, config, r, c, candidate, &log_combinatorial_terms);
                }, kRhoLower, config.overdispersion_max, 100, config.tolerance,
                inner_iterations);
            const double boundary_snap = std::max(10.0 * config.tolerance, 1e-7);
            if (rho <= kRhoLower + boundary_snap) rho = kRhoLower;
            if (rho >= config.overdispersion_max - boundary_snap) {
                rho = config.overdispersion_max;
            }
        }

        MtRatioFitResult fit;
        fit.joint_ambient = ambient_free;
        const double ll = mt_ratio_log_likelihood_cached(
            observations, log_combinatorial_terms, r, c, rho, use_beta);
        const double objective = ratio_fit_objective(
            observations, config, r, c, rho, &log_combinatorial_terms);
        if (!std::isfinite(ll) || !std::isfinite(objective)) return fit;

        fit.parent2_fraction = r;
        fit.ambient_fraction = c;
        fit.parent1_weight = (1.0 - c) * (1.0 - r);
        fit.parent2_weight = (1.0 - c) * r;
        fit.overdispersion_rho = use_beta ? rho : 0.0;
        fit.log_likelihood = ll;
        fit.objective_value = objective;
        fit.iterations = 1;
        fit.converged = true;
        return fit;
    }

    // With two free nuisance parameters, use the adjacent profile point as a
    // true continuation attempt.  Pay the six global deterministic starts
    // only when that ordinary finite/converged attempt fails.
    if (warm_start) {
        MtRatioFitResult continuation =
            fit_internal(observations, config, r, warm_start, true);
        if (continuation.converged &&
            std::isfinite(continuation.log_likelihood) &&
            std::isfinite(continuation.objective_value)) {
            return continuation;
        }
    }
    return fit_internal(observations, config, r, nullptr, false);
}

MtRatioProfileResult profile_mt_ratio(const std::vector<MtRatioObservation>& observations,
                                      const MtRatioFitConfig& config,
                                      const MtRatioFitResult& unconstrained_fit,
                                      double single_parent_epsilon,
                                      double grid_step,
                                      bool retain_grid) {
    MtRatioProfileResult out;
    out.single_parent_epsilon = single_parent_epsilon;
    if (!profile_fit_valid(unconstrained_fit) || observations.empty() ||
        !(single_parent_epsilon > 0.0 && single_parent_epsilon < 0.5) ||
        !(grid_step > 0.0 && grid_step <= 0.05)) {
        out.inheritance_class_reason = "UNCONSTRAINED_FIT_NOT_PROFILEABLE";
        return out;
    }

    const int n_intervals = std::max(1, static_cast<int>(std::ceil(1.0 / grid_step - 1e-12)));
    const double actual_grid_step = 1.0 / static_cast<double>(n_intervals);
    const int n_points = n_intervals + 1;
    std::vector<ProfilePoint> points(static_cast<size_t>(n_points));
    const int center = std::max(0, std::min(n_intervals,
        static_cast<int>(std::llround(unconstrained_fit.parent2_fraction / actual_grid_step))));

    auto eval_index = [&](int idx, const MtRatioFitResult* warm) {
        const double r = idx == n_intervals ? 1.0 : idx * actual_grid_step;
        points[static_cast<size_t>(idx)] =
            evaluate_profile_point(observations, config, unconstrained_fit, r, warm);
        ++out.evaluations;
        if (!points[static_cast<size_t>(idx)].ok) ++out.failed_evaluations;
    };

    eval_index(center, &unconstrained_fit);
    const MtRatioFitResult* warm_left = points[static_cast<size_t>(center)].ok
        ? &points[static_cast<size_t>(center)].fit : &unconstrained_fit;
    for (int idx = center - 1; idx >= 0; --idx) {
        eval_index(idx, warm_left);
        if (points[static_cast<size_t>(idx)].ok) warm_left = &points[static_cast<size_t>(idx)].fit;
    }
    const MtRatioFitResult* warm_right = points[static_cast<size_t>(center)].ok
        ? &points[static_cast<size_t>(center)].fit : &unconstrained_fit;
    for (int idx = center + 1; idx < n_points; ++idx) {
        eval_index(idx, warm_right);
        if (points[static_cast<size_t>(idx)].ok) warm_right = &points[static_cast<size_t>(idx)].fit;
    }

    if (retain_grid) {
        out.grid_r.reserve(static_cast<size_t>(n_points));
        out.grid_delta_log_likelihood.reserve(static_cast<size_t>(n_points));
        for (int idx = 0; idx < n_points; ++idx) {
            const double r = idx == n_intervals ? 1.0 : idx * actual_grid_step;
            out.grid_r.push_back(r);
            out.grid_delta_log_likelihood.push_back(points[static_cast<size_t>(idx)].ok
                ? points[static_cast<size_t>(idx)].delta
                : std::numeric_limits<double>::quiet_NaN());
        }
    }

    const double target = -kProfileDrop95;
    std::vector<int> accepted;
    for (int idx = 0; idx < n_points; ++idx) {
        const ProfilePoint& p = points[static_cast<size_t>(idx)];
        if (p.ok && p.delta >= target) accepted.push_back(idx);
    }

    bool separated_resolved_modes = false;
    for (size_t k = 1; k < accepted.size(); ++k) {
        if (accepted[k] == accepted[k - 1] + 1) continue;
        bool gap_fully_resolved = true;
        bool resolved_rejection = false;
        for (int idx = accepted[k - 1] + 1; idx < accepted[k]; ++idx) {
            const ProfilePoint& gap = points[static_cast<size_t>(idx)];
            if (!gap.ok) {
                gap_fully_resolved = false;
                break;
            }
            if (gap.delta < target) resolved_rejection = true;
        }
        if (gap_fully_resolved && resolved_rejection) {
            separated_resolved_modes = true;
            break;
        }
    }

    // The unconstrained MLE is the guaranteed interior point of its own profile
    // interval.  Anchor both searches there rather than requiring a coarse grid
    // point to land inside a potentially much narrower interval.
    const double mle_r = std::max(0.0, std::min(1.0, unconstrained_fit.parent2_fraction));
    double ci_low = mle_r;
    double ci_high = mle_r;
    bool lower_refinement_failed = false;
    bool upper_refinement_failed = false;
    bool lower_resolved = mle_r <= 1e-12;
    bool upper_resolved = mle_r >= 1.0 - 1e-12;
    if (lower_resolved) ci_low = 0.0;
    if (upper_resolved) ci_high = 1.0;

    if (!lower_resolved) {
        double inside_r = mle_r;
        const MtRatioFitResult* inside_fit = &unconstrained_fit;
        for (int idx = n_intervals; idx >= 0; --idx) {
            const ProfilePoint& p = points[static_cast<size_t>(idx)];
            if (!p.ok || p.r >= mle_r - 1e-12) continue;
            if (p.delta >= target) {
                inside_r = p.r;
                inside_fit = &p.fit;
                if (p.r <= 1e-12) {
                    ci_low = 0.0;
                    lower_resolved = true;
                    break;
                }
                continue;
            }
            ci_low = refine_crossing(observations, config, unconstrained_fit,
                                     inside_r, p.r, inside_fit, out.evaluations,
                                     out.failed_evaluations, lower_refinement_failed);
            lower_resolved = !lower_refinement_failed;
            break;
        }
    }

    if (!upper_resolved) {
        double inside_r = mle_r;
        const MtRatioFitResult* inside_fit = &unconstrained_fit;
        for (int idx = 0; idx <= n_intervals; ++idx) {
            const ProfilePoint& p = points[static_cast<size_t>(idx)];
            if (!p.ok || p.r <= mle_r + 1e-12) continue;
            if (p.delta >= target) {
                inside_r = p.r;
                inside_fit = &p.fit;
                if (p.r >= 1.0 - 1e-12) {
                    ci_high = 1.0;
                    upper_resolved = true;
                    break;
                }
                continue;
            }
            ci_high = refine_crossing(observations, config, unconstrained_fit,
                                      inside_r, p.r, inside_fit, out.evaluations,
                                      out.failed_evaluations, upper_refinement_failed);
            upper_resolved = !upper_refinement_failed;
            break;
        }
    }

    if (!lower_resolved && !upper_resolved) {
        out.profile_status = "FAILED";
        out.inheritance_class = "AMBIGUOUS";
        out.inheritance_class_reason = "PROFILE_INTERVAL_BOUNDARIES_UNRESOLVED";
        return out;
    }
    out.parent2_ci_low = std::max(0.0, std::min(1.0, ci_low));
    out.parent2_ci_high = std::max(0.0, std::min(1.0, ci_high));
    out.parent1_ci_low = 1.0 - out.parent2_ci_high;
    out.parent1_ci_high = 1.0 - out.parent2_ci_low;

    std::vector<ProfilePoint> boundary_cache;
    auto boundary_point = [&](double r) {
        const int grid_index = std::max(0, std::min(n_intervals,
            static_cast<int>(std::llround(r / actual_grid_step))));
        const ProfilePoint& grid_point = points[static_cast<size_t>(grid_index)];
        if (std::fabs(grid_point.r - r) <= 1e-12 && grid_point.ok) return grid_point;

        for (const ProfilePoint& cached : boundary_cache) {
            if (std::fabs(cached.r - r) <= 1e-12) return cached;
        }

        ProfilePoint evaluated = evaluate_profile_point(
            observations, config, unconstrained_fit, r, &unconstrained_fit);
        ++out.evaluations;
        if (!evaluated.ok) ++out.failed_evaluations;
        boundary_cache.push_back(evaluated);
        return evaluated;
    };

    auto best_region_delta = [&](double lower, double upper,
                                 bool include_lower, bool include_upper) {
        const double mle_r = unconstrained_fit.parent2_fraction;
        const bool mle_ge = include_lower ? mle_r >= lower - 1e-12 : mle_r > lower + 1e-12;
        const bool mle_le = include_upper ? mle_r <= upper + 1e-12 : mle_r < upper - 1e-12;
        double best = (mle_ge && mle_le) ? 0.0 : -std::numeric_limits<double>::infinity();
        for (const ProfilePoint& p : points) {
            if (!p.ok) continue;
            const bool ge = include_lower ? p.r >= lower - 1e-12 : p.r > lower + 1e-12;
            const bool le = include_upper ? p.r <= upper + 1e-12 : p.r < upper - 1e-12;
            if (ge && le) best = std::max(best, p.delta);
        }
        // Reuse resolved grid boundaries and evaluate each off-grid boundary at
        // most once.  A failed grid point still gets the historical retry with
        // the unconstrained fit as the warm start.
        for (double r : {lower, upper}) {
            if (r < 0.0 || r > 1.0) continue;
            const ProfilePoint p = boundary_point(r);
            if (!p.ok) continue;
            const bool ge = include_lower ? p.r >= lower - 1e-12 : p.r > lower + 1e-12;
            const bool le = include_upper ? p.r <= upper + 1e-12 : p.r < upper - 1e-12;
            if (ge && le) best = std::max(best, p.delta);
        }
        if (!std::isfinite(best)) return std::numeric_limits<double>::quiet_NaN();
        return std::max(0.0, -best);
    };

    out.delta_ll_parent1_only = best_region_delta(0.0, single_parent_epsilon, true, true);
    out.delta_ll_both = best_region_delta(single_parent_epsilon,
                                          1.0 - single_parent_epsilon,
                                          false, false);
    out.delta_ll_parent2_only = best_region_delta(1.0 - single_parent_epsilon, 1.0,
                                                  true, true);

    if (separated_resolved_modes) {
        out.profile_status = "MULTIMODAL";
    } else if (!lower_resolved || !upper_resolved ||
               lower_refinement_failed || upper_refinement_failed) {
        out.profile_status = "PARTIAL";
    } else {
        out.profile_status = "PASS";
    }

    if (out.profile_status != "PASS") {
        out.inheritance_class = "AMBIGUOUS";
        out.inheritance_class_reason = out.profile_status == "MULTIMODAL"
            ? "PROFILE_MULTIMODAL" : "PROFILE_PARTIALLY_RESOLVED";
    } else if (out.parent2_ci_high <= single_parent_epsilon) {
        out.inheritance_class = "ONLY_PARENT1";
        out.inheritance_class_reason = "PROFILE_INTERVAL_WITHIN_PARENT1_REGION";
    } else if (out.parent2_ci_low >= 1.0 - single_parent_epsilon) {
        out.inheritance_class = "ONLY_PARENT2";
        out.inheritance_class_reason = "PROFILE_INTERVAL_WITHIN_PARENT2_REGION";
    } else if (out.parent2_ci_low > single_parent_epsilon &&
               out.parent2_ci_high < 1.0 - single_parent_epsilon) {
        out.inheritance_class = "BOTH";
        out.inheritance_class_reason = "PROFILE_INTERVAL_WITHIN_BOTH_REGION";
    } else {
        out.inheritance_class = "AMBIGUOUS";
        out.inheritance_class_reason = "PROFILE_INTERVAL_CROSSES_INHERITANCE_REGION";
    }
    return out;
}

namespace {

double log_weighted_pair(double log_a, double weight_a,
                         double log_b, double weight_b) {
    if (weight_a <= 0.0) return log_b + std::log(std::max(weight_b, kTiny));
    if (weight_b <= 0.0) return log_a + std::log(std::max(weight_a, kTiny));
    const double a = log_a + std::log(weight_a);
    const double b = log_b + std::log(weight_b);
    const double m = std::max(a, b);
    return m + std::log(std::exp(a - m) + std::exp(b - m));
}

double fragment_site_log_emission(const MtRatioFragmentSite& site,
                                  double alt_probability) {
    const double p = clamp_probability(alt_probability);
    return site.allele == 1 ? std::log(p) : std::log1p(-p);
}

double fragment_observation_log_likelihood(const MtRatioFragmentObservation& fragment,
                                           double parent2_fraction,
                                           double ambient_fraction) {
    if (fragment.sites.empty()) return 0.0;
    double lp1 = 0.0;
    double lp2 = 0.0;
    double la = 0.0;
    for (const MtRatioFragmentSite& site : fragment.sites) {
        lp1 += fragment_site_log_emission(site, site.parent1_alt_probability);
        lp2 += fragment_site_log_emission(site, site.parent2_alt_probability);
        la += fragment_site_log_emission(site, site.ambient_alt_probability);
    }
    const double r = std::max(0.0, std::min(1.0, parent2_fraction));
    const double c = std::max(0.0, std::min(0.99, ambient_fraction));
    const double parental = log_weighted_pair(lp1, 1.0 - r, lp2, r);
    if (c <= 0.0) return parental;
    return log_weighted_pair(parental, 1.0 - c, la, c);
}

MtRatioFitResult fit_fragment_internal(
        const std::vector<MtRatioFragmentObservation>& fragments,
        const MtRatioFitConfig& config,
        const std::optional<double>& fixed_r,
        const MtRatioFitResult* warm_start) {
    const bool ambient_free = config.ambient_enabled && !config.ambient_fixed;
    const double fixed_c = !config.ambient_enabled ? 0.0 :
                           (config.ambient_fixed ? config.fixed_ambient : 0.0);
    std::vector<FitStart> starts;
    if (warm_start && std::isfinite(warm_start->parent2_fraction) &&
        std::isfinite(warm_start->ambient_fraction)) {
        starts.push_back({fixed_r ? *fixed_r : warm_start->parent2_fraction,
                          ambient_free ? warm_start->ambient_fraction : fixed_c, 0.0});
    }
    const std::vector<FitStart> defaults = deterministic_starts(config);
    starts.insert(starts.end(), defaults.begin(), defaults.end());

    MtRatioFitResult best;
    best.joint_ambient = ambient_free;
    double best_ll = -std::numeric_limits<double>::infinity();
    bool found_converged = false;
    const bool prefer_converged = fixed_r.has_value();

    for (const FitStart& start : starts) {
        double r = fixed_r ? std::max(0.0, std::min(1.0, *fixed_r))
                           : std::max(0.0, std::min(1.0, start.ratio));
        double c = ambient_free ? std::max(0.0, std::min(0.99, start.ambient)) : fixed_c;
        bool converged = false;
        int iteration = 0;
        double previous_ll = mt_fragment_ratio_log_likelihood(fragments, r, c);
        if (!std::isfinite(previous_ll)) continue;
        for (; iteration < config.max_iterations; ++iteration) {
            const double old_r = r;
            const double old_c = c;
            int inner = 0;
            if (!fixed_r) {
                r = golden_maximize(
                    [&](double candidate) {
                        return mt_fragment_ratio_log_likelihood(fragments, candidate, c);
                    }, 0.0, 1.0, 100, config.tolerance, inner);
            }
            if (ambient_free) {
                c = golden_maximize(
                    [&](double candidate) {
                        return mt_fragment_ratio_log_likelihood(fragments, r, candidate);
                    }, 0.0, 0.99, 100, config.tolerance, inner);
            }
            const double ll = mt_fragment_ratio_log_likelihood(fragments, r, c);
            if (!std::isfinite(ll)) break;
            const double parameter_delta = std::max(std::fabs(r - old_r), std::fabs(c - old_c));
            const double relative_ll_delta = std::fabs(ll - previous_ll) /
                                             (1.0 + std::fabs(previous_ll));
            previous_ll = ll;
            if (parameter_delta <= std::max(10.0 * config.tolerance, 1e-7) &&
                relative_ll_delta <= std::max(config.tolerance, 1e-10)) {
                converged = true;
                ++iteration;
                break;
            }
        }
        const double ll = mt_fragment_ratio_log_likelihood(fragments, r, c);
        if (!std::isfinite(ll)) continue;
        const bool take = prefer_converged
            ? ((converged && !found_converged) ||
               (converged == found_converged && ll > best_ll))
            : (ll > best_ll);
        if (!take) continue;
        found_converged = found_converged || converged;
        best_ll = ll;
        best.parent2_fraction = r;
        best.ambient_fraction = c;
        best.parent1_weight = (1.0 - c) * (1.0 - r);
        best.parent2_weight = (1.0 - c) * r;
        best.overdispersion_rho = 0.0;
        best.log_likelihood = ll;
        best.objective_value = ll;
        best.iterations = iteration;
        best.converged = converged;
        best.joint_ambient = ambient_free;
    }
    return best;
}

void compute_fragment_uncertainty(
        const std::vector<MtRatioFragmentObservation>& fragments,
        const MtRatioFitConfig& config,
        MtRatioFitResult& fit) {
    const bool ambient_free = config.ambient_enabled && !config.ambient_fixed;
    const double r = fit.parent2_fraction;
    const double c = fit.ambient_fraction;
    const double hr = finite_step(r, 0.0, 1.0);
    const double hc = ambient_free ? finite_step(c, 0.0, 0.99) : 0.0;
    if (hr <= 0.0) return;
    const double center = mt_fragment_ratio_log_likelihood(fragments, r, c);
    const auto ll = [&](double rr, double cc) {
        return mt_fragment_ratio_log_likelihood(fragments, rr, cc);
    };
    const double i_rr = -(ll(r + hr, c) - 2.0 * center + ll(r - hr, c)) / (hr * hr);
    if (!ambient_free) {
        if (i_rr > kTiny && std::isfinite(i_rr)) {
            fit.ratio_se = std::sqrt(1.0 / i_rr);
            double b_rr = 0.0;
            for (const MtRatioFragmentObservation& fragment : fragments) {
                const double score =
                    (fragment_observation_log_likelihood(fragment, r + hr, c) -
                     fragment_observation_log_likelihood(fragment, r - hr, c)) / (2.0 * hr);
                b_rr += score * score;
            }
            fit.ratio_se_robust = std::sqrt(std::max(0.0, b_rr / (i_rr * i_rr)));
            fit.information_condition = 1.0;
            fit.min_information_eigenvalue = i_rr;
        }
        return;
    }
    if (hc <= 0.0) return;
    const double i_cc = -(ll(r, c + hc) - 2.0 * center + ll(r, c - hc)) / (hc * hc);
    const double i_rc = -(ll(r + hr, c + hc) - ll(r + hr, c - hc) -
                          ll(r - hr, c + hc) + ll(r - hr, c - hc)) / (4.0 * hr * hc);
    const double trace = i_rr + i_cc;
    const double disc = std::sqrt(std::max(0.0,
        (i_rr - i_cc) * (i_rr - i_cc) + 4.0 * i_rc * i_rc));
    const double lambda_max = 0.5 * (trace + disc);
    const double lambda_min = 0.5 * (trace - disc);
    fit.min_information_eigenvalue = lambda_min;
    fit.information_condition = lambda_max / std::max(lambda_min, kTiny);
    const double determinant = i_rr * i_cc - i_rc * i_rc;
    if (!(determinant > kTiny && i_rr > 0.0 && i_cc > 0.0 && std::isfinite(determinant))) {
        return;
    }
    const double inv_rr = i_cc / determinant;
    const double inv_rc = -i_rc / determinant;
    const double inv_cc = i_rr / determinant;
    fit.ratio_se = std::sqrt(std::max(0.0, inv_rr));
    fit.ambient_se = std::sqrt(std::max(0.0, inv_cc));
    double b_rr = 0.0, b_rc = 0.0, b_cc = 0.0;
    for (const MtRatioFragmentObservation& fragment : fragments) {
        const double sr =
            (fragment_observation_log_likelihood(fragment, r + hr, c) -
             fragment_observation_log_likelihood(fragment, r - hr, c)) / (2.0 * hr);
        const double sc =
            (fragment_observation_log_likelihood(fragment, r, c + hc) -
             fragment_observation_log_likelihood(fragment, r, c - hc)) / (2.0 * hc);
        b_rr += sr * sr;
        b_rc += sr * sc;
        b_cc += sc * sc;
    }
    const double v_rr = inv_rr * (b_rr * inv_rr + b_rc * inv_rc) +
                        inv_rc * (b_rc * inv_rr + b_cc * inv_rc);
    const double v_cc = inv_rc * (b_rr * inv_rc + b_rc * inv_cc) +
                        inv_cc * (b_rc * inv_rc + b_cc * inv_cc);
    fit.ratio_se_robust = std::sqrt(std::max(0.0, v_rr));
    fit.ambient_se_robust = std::sqrt(std::max(0.0, v_cc));
}

struct FragmentProfilePoint {
    double r = std::numeric_limits<double>::quiet_NaN();
    double delta = std::numeric_limits<double>::quiet_NaN();
    MtRatioFitResult fit;
    bool ok = false;
};

FragmentProfilePoint evaluate_fragment_profile_point(
        const std::vector<MtRatioFragmentObservation>& fragments,
        const MtRatioFitConfig& config,
        const MtRatioFitResult& unconstrained_fit,
        double r,
        const MtRatioFitResult* warm_start) {
    FragmentProfilePoint point;
    point.r = r;
    point.fit = fit_mt_fragment_ratio_at_fixed_fraction(fragments, config, r, warm_start);
    point.ok = profile_fit_valid(point.fit);
    if (point.ok) {
        point.delta = point.fit.log_likelihood - unconstrained_fit.log_likelihood;
        if (point.delta > 0.0) point.delta = 0.0;
    }
    return point;
}

double refine_fragment_crossing(
        const std::vector<MtRatioFragmentObservation>& fragments,
        const MtRatioFitConfig& config,
        const MtRatioFitResult& unconstrained_fit,
        double inside_r,
        double outside_r,
        const MtRatioFitResult* warm_start,
        int& evaluations,
        int& failed_evaluations,
        bool& failed) {
    double inside = inside_r;
    double outside = outside_r;
    MtRatioFitResult warm = warm_start ? *warm_start : unconstrained_fit;
    const double target = -kProfileDrop95;
    for (int iter = 0; iter < 80; ++iter) {
        const double mid = 0.5 * (inside + outside);
        FragmentProfilePoint p = evaluate_fragment_profile_point(
            fragments, config, unconstrained_fit, mid, &warm);
        ++evaluations;
        if (!p.ok) {
            ++failed_evaluations;
            failed = true;
            return inside;
        }
        warm = p.fit;
        if (p.delta >= target) inside = mid;
        else outside = mid;
        if (std::fabs(outside - inside) <= 1e-5 ||
            std::fabs(p.delta - target) <= std::max(config.tolerance, 1e-8)) break;
    }
    return 0.5 * (inside + outside);
}

}  // namespace

double mt_fragment_ratio_log_likelihood(
        const std::vector<MtRatioFragmentObservation>& fragments,
        double parent2_fraction,
        double ambient_fraction) {
    double total = 0.0;
    for (const MtRatioFragmentObservation& fragment : fragments) {
        total += fragment_observation_log_likelihood(fragment, parent2_fraction, ambient_fraction);
    }
    return total;
}

MtRatioFitResult fit_mt_fragment_ratio(
        const std::vector<MtRatioFragmentObservation>& fragments,
        const MtRatioFitConfig& config) {
    MtRatioFitResult fit = fit_fragment_internal(fragments, config, std::nullopt, nullptr);
    if (std::isfinite(fit.log_likelihood)) compute_fragment_uncertainty(fragments, config, fit);
    return fit;
}

MtRatioFitResult fit_mt_fragment_ratio_at_fixed_fraction(
        const std::vector<MtRatioFragmentObservation>& fragments,
        const MtRatioFitConfig& config,
        double fixed_parent2_fraction,
        const MtRatioFitResult* warm_start) {
    (void)warm_start;
    const double r = std::max(0.0, std::min(1.0, fixed_parent2_fraction));
    const bool ambient_free = config.ambient_enabled && !config.ambient_fixed;

    // A fixed-r fragment fit has at most one nuisance parameter (ambient).
    // Avoid repeating the same bounded search for every deterministic start.
    if (!ambient_free) {
        const double c = !config.ambient_enabled ? 0.0 : config.fixed_ambient;
        MtRatioFitResult fit;
        fit.joint_ambient = false;
        const double ll = mt_fragment_ratio_log_likelihood(fragments, r, c);
        if (!std::isfinite(ll)) return fit;
        fit.parent2_fraction = r;
        fit.ambient_fraction = c;
        fit.parent1_weight = (1.0 - c) * (1.0 - r);
        fit.parent2_weight = (1.0 - c) * r;
        fit.overdispersion_rho = 0.0;
        fit.log_likelihood = ll;
        fit.objective_value = ll;
        fit.iterations = 1;
        fit.converged = true;
        return fit;
    }

    int inner_iterations = 0;
    const double c = golden_maximize(
        [&](double candidate) {
            return mt_fragment_ratio_log_likelihood(fragments, r, candidate);
        }, 0.0, 0.99, 100, config.tolerance, inner_iterations);
    MtRatioFitResult fit;
    fit.joint_ambient = true;
    const double ll = mt_fragment_ratio_log_likelihood(fragments, r, c);
    if (!std::isfinite(ll)) return fit;
    fit.parent2_fraction = r;
    fit.ambient_fraction = c;
    fit.parent1_weight = (1.0 - c) * (1.0 - r);
    fit.parent2_weight = (1.0 - c) * r;
    fit.overdispersion_rho = 0.0;
    fit.log_likelihood = ll;
    fit.objective_value = ll;
    fit.iterations = 1;
    fit.converged = true;
    return fit;
}

MtRatioProfileResult profile_mt_fragment_ratio(
        const std::vector<MtRatioFragmentObservation>& fragments,
        const MtRatioFitConfig& config,
        const MtRatioFitResult& unconstrained_fit,
        double single_parent_epsilon,
        double grid_step,
        bool retain_grid) {
    MtRatioProfileResult out;
    out.single_parent_epsilon = single_parent_epsilon;
    if (!profile_fit_valid(unconstrained_fit) || fragments.empty() ||
        !(single_parent_epsilon > 0.0 && single_parent_epsilon < 0.5) ||
        !(grid_step > 0.0 && grid_step <= 0.05)) {
        out.inheritance_class_reason = "UNCONSTRAINED_FIT_NOT_PROFILEABLE";
        return out;
    }
    const int n_intervals = std::max(1, static_cast<int>(std::ceil(1.0 / grid_step - 1e-12)));
    const double actual_grid_step = 1.0 / static_cast<double>(n_intervals);
    const int n_points = n_intervals + 1;
    std::vector<FragmentProfilePoint> points(static_cast<size_t>(n_points));
    const int center = std::max(0, std::min(n_intervals,
        static_cast<int>(std::llround(unconstrained_fit.parent2_fraction / actual_grid_step))));
    auto eval_index = [&](int idx, const MtRatioFitResult* warm) {
        const double r = idx == n_intervals ? 1.0 : idx * actual_grid_step;
        points[static_cast<size_t>(idx)] = evaluate_fragment_profile_point(
            fragments, config, unconstrained_fit, r, warm);
        ++out.evaluations;
        if (!points[static_cast<size_t>(idx)].ok) ++out.failed_evaluations;
    };
    eval_index(center, &unconstrained_fit);
    const MtRatioFitResult* warm_left = points[static_cast<size_t>(center)].ok
        ? &points[static_cast<size_t>(center)].fit : &unconstrained_fit;
    for (int idx = center - 1; idx >= 0; --idx) {
        eval_index(idx, warm_left);
        if (points[static_cast<size_t>(idx)].ok) warm_left = &points[static_cast<size_t>(idx)].fit;
    }
    const MtRatioFitResult* warm_right = points[static_cast<size_t>(center)].ok
        ? &points[static_cast<size_t>(center)].fit : &unconstrained_fit;
    for (int idx = center + 1; idx < n_points; ++idx) {
        eval_index(idx, warm_right);
        if (points[static_cast<size_t>(idx)].ok) warm_right = &points[static_cast<size_t>(idx)].fit;
    }
    if (retain_grid) {
        out.grid_r.reserve(static_cast<size_t>(n_points));
        out.grid_delta_log_likelihood.reserve(static_cast<size_t>(n_points));
        for (int idx = 0; idx < n_points; ++idx) {
            out.grid_r.push_back(idx == n_intervals ? 1.0 : idx * actual_grid_step);
            out.grid_delta_log_likelihood.push_back(points[static_cast<size_t>(idx)].ok
                ? points[static_cast<size_t>(idx)].delta
                : std::numeric_limits<double>::quiet_NaN());
        }
    }
    const double target = -kProfileDrop95;
    std::vector<int> accepted;
    for (int idx = 0; idx < n_points; ++idx) {
        if (points[static_cast<size_t>(idx)].ok && points[static_cast<size_t>(idx)].delta >= target) {
            accepted.push_back(idx);
        }
    }
    bool separated_resolved_modes = false;
    for (size_t k = 1; k < accepted.size(); ++k) {
        if (accepted[k] == accepted[k - 1] + 1) continue;
        bool gap_fully_resolved = true;
        bool resolved_rejection = false;
        for (int idx = accepted[k - 1] + 1; idx < accepted[k]; ++idx) {
            const FragmentProfilePoint& gap = points[static_cast<size_t>(idx)];
            if (!gap.ok) { gap_fully_resolved = false; break; }
            if (gap.delta < target) resolved_rejection = true;
        }
        if (gap_fully_resolved && resolved_rejection) {
            separated_resolved_modes = true;
            break;
        }
    }

    const double mle_r = std::max(0.0, std::min(1.0, unconstrained_fit.parent2_fraction));
    double ci_low = mle_r;
    double ci_high = mle_r;
    bool lower_failed = false, upper_failed = false;
    bool lower_resolved = mle_r <= 1e-12;
    bool upper_resolved = mle_r >= 1.0 - 1e-12;
    if (lower_resolved) ci_low = 0.0;
    if (upper_resolved) ci_high = 1.0;

    if (!lower_resolved) {
        double inside_r = mle_r;
        const MtRatioFitResult* inside_fit = &unconstrained_fit;
        for (int idx = n_intervals; idx >= 0; --idx) {
            const FragmentProfilePoint& p = points[static_cast<size_t>(idx)];
            if (!p.ok || p.r >= mle_r - 1e-12) continue;
            if (p.delta >= target) {
                inside_r = p.r;
                inside_fit = &p.fit;
                if (p.r <= 1e-12) {
                    ci_low = 0.0;
                    lower_resolved = true;
                    break;
                }
                continue;
            }
            ci_low = refine_fragment_crossing(fragments, config, unconstrained_fit,
                inside_r, p.r, inside_fit, out.evaluations,
                out.failed_evaluations, lower_failed);
            lower_resolved = !lower_failed;
            break;
        }
    }

    if (!upper_resolved) {
        double inside_r = mle_r;
        const MtRatioFitResult* inside_fit = &unconstrained_fit;
        for (int idx = 0; idx <= n_intervals; ++idx) {
            const FragmentProfilePoint& p = points[static_cast<size_t>(idx)];
            if (!p.ok || p.r <= mle_r + 1e-12) continue;
            if (p.delta >= target) {
                inside_r = p.r;
                inside_fit = &p.fit;
                if (p.r >= 1.0 - 1e-12) {
                    ci_high = 1.0;
                    upper_resolved = true;
                    break;
                }
                continue;
            }
            ci_high = refine_fragment_crossing(fragments, config, unconstrained_fit,
                inside_r, p.r, inside_fit, out.evaluations,
                out.failed_evaluations, upper_failed);
            upper_resolved = !upper_failed;
            break;
        }
    }

    if (!lower_resolved && !upper_resolved) {
        out.profile_status = "FAILED";
        out.inheritance_class = "AMBIGUOUS";
        out.inheritance_class_reason = "PROFILE_INTERVAL_BOUNDARIES_UNRESOLVED";
        return out;
    }
    out.parent2_ci_low = std::max(0.0, std::min(1.0, ci_low));
    out.parent2_ci_high = std::max(0.0, std::min(1.0, ci_high));
    out.parent1_ci_low = 1.0 - out.parent2_ci_high;
    out.parent1_ci_high = 1.0 - out.parent2_ci_low;

    std::vector<FragmentProfilePoint> boundary_cache;
    auto boundary_point = [&](double r) {
        const int grid_index = std::max(0, std::min(n_intervals,
            static_cast<int>(std::llround(r / actual_grid_step))));
        const FragmentProfilePoint& grid_point = points[static_cast<size_t>(grid_index)];
        if (std::fabs(grid_point.r - r) <= 1e-12 && grid_point.ok) return grid_point;

        for (const FragmentProfilePoint& cached : boundary_cache) {
            if (std::fabs(cached.r - r) <= 1e-12) return cached;
        }

        FragmentProfilePoint evaluated = evaluate_fragment_profile_point(
            fragments, config, unconstrained_fit, r, &unconstrained_fit);
        ++out.evaluations;
        if (!evaluated.ok) ++out.failed_evaluations;
        boundary_cache.push_back(evaluated);
        return evaluated;
    };

    auto best_region_delta = [&](double lower, double upper,
                                 bool include_lower, bool include_upper) {
        const double mle_r = unconstrained_fit.parent2_fraction;
        const bool mle_ge = include_lower ? mle_r >= lower - 1e-12 : mle_r > lower + 1e-12;
        const bool mle_le = include_upper ? mle_r <= upper + 1e-12 : mle_r < upper - 1e-12;
        double best = (mle_ge && mle_le) ? 0.0 : -std::numeric_limits<double>::infinity();
        for (const FragmentProfilePoint& p : points) {
            if (!p.ok) continue;
            const bool ge = include_lower ? p.r >= lower - 1e-12 : p.r > lower + 1e-12;
            const bool le = include_upper ? p.r <= upper + 1e-12 : p.r < upper - 1e-12;
            if (ge && le) best = std::max(best, p.delta);
        }
        for (double r : {lower, upper}) {
            const FragmentProfilePoint p = boundary_point(r);
            if (!p.ok) continue;
            const bool ge = include_lower ? p.r >= lower - 1e-12 : p.r > lower + 1e-12;
            const bool le = include_upper ? p.r <= upper + 1e-12 : p.r < upper - 1e-12;
            if (ge && le) best = std::max(best, p.delta);
        }
        return std::isfinite(best) ? std::max(0.0, -best)
                                   : std::numeric_limits<double>::quiet_NaN();
    };
    out.delta_ll_parent1_only = best_region_delta(0.0, single_parent_epsilon, true, true);
    out.delta_ll_both = best_region_delta(single_parent_epsilon, 1.0 - single_parent_epsilon,
                                          false, false);
    out.delta_ll_parent2_only = best_region_delta(1.0 - single_parent_epsilon, 1.0, true, true);
    if (separated_resolved_modes) out.profile_status = "MULTIMODAL";
    else if (!lower_resolved || !upper_resolved ||
             lower_failed || upper_failed) out.profile_status = "PARTIAL";
    else out.profile_status = "PASS";

    if (out.profile_status != "PASS") {
        out.inheritance_class = "AMBIGUOUS";
        out.inheritance_class_reason = out.profile_status == "MULTIMODAL"
            ? "PROFILE_MULTIMODAL" : "PROFILE_PARTIALLY_RESOLVED";
    } else if (out.parent2_ci_high <= single_parent_epsilon) {
        out.inheritance_class = "ONLY_PARENT1";
        out.inheritance_class_reason = "PROFILE_INTERVAL_WITHIN_PARENT1_REGION";
    } else if (out.parent2_ci_low >= 1.0 - single_parent_epsilon) {
        out.inheritance_class = "ONLY_PARENT2";
        out.inheritance_class_reason = "PROFILE_INTERVAL_WITHIN_PARENT2_REGION";
    } else if (out.parent2_ci_low > single_parent_epsilon &&
               out.parent2_ci_high < 1.0 - single_parent_epsilon) {
        out.inheritance_class = "BOTH";
        out.inheritance_class_reason = "PROFILE_INTERVAL_WITHIN_BOTH_REGION";
    } else {
        out.inheritance_class = "AMBIGUOUS";
        out.inheritance_class_reason = "PROFILE_INTERVAL_CROSSES_INHERITANCE_REGION";
    }
    return out;
}
