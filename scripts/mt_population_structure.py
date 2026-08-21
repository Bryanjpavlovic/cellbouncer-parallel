#!/usr/bin/env python3
"""Likelihood-profile population structure for mt_fusion_ratio outputs.

First-wave models:
  M0: one common canonical parent-2 ratio.
  M1: one truncated-normal distribution of latent cell ratios.
  M2: two discrete latent ratio subpopulations.

The program never modifies its per-cell inputs. Site-calibration provenance is
propagated from the C++ outputs, and an optional comparison file can report
changes in population-structure calls between calibrated and reference runs.
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import platform
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.special import logsumexp, ndtr

PROFILE_DROP_95 = 1.920729410347062
EPS = 1e-300
POPULATION_IMPLEMENTATION_VERSION = "mt_population_structure_v2_correctness_20260818"
M2_NEAR_BEST_LL_TOLERANCE = 0.1
M2_PROPORTION_DISAGREEMENT_TOLERANCE = 0.05


@dataclass
class M0Result:
    ratio: float
    ci_low: float
    ci_high: float
    log_likelihood: float
    bic: float
    grid_index: int


@dataclass
class M1Result:
    location: float
    tau: float
    distribution_mean: float
    distribution_sd: float
    sd_upper_95: float
    log_likelihood: float
    bic: float
    fit_status: str


@dataclass
class M2Result:
    ratio_low: float
    ratio_high: float
    proportion_low: float
    proportion_high: float
    separation: float
    confident_membership_fraction: float
    effective_low: float
    effective_high: float
    log_likelihood: float
    bic: float
    p_low: np.ndarray
    p_high: np.ndarray
    fit_status: str


def normalize_profiles(log_likelihoods: np.ndarray) -> np.ndarray:
    arr = np.asarray(log_likelihoods, dtype=float)
    if arr.ndim != 2:
        raise ValueError("profile likelihood array must be two-dimensional")
    if not np.all(np.isfinite(arr)):
        raise ValueError("profile likelihood array contains nonfinite values")
    return arr - np.max(arr, axis=1, keepdims=True)


def _refined_m0_curve(grid: np.ndarray, curve: np.ndarray, idx: int):
    """Return one coherent local likelihood representation around the M0 peak."""
    coefficients = None
    local_low = local_high = math.nan
    ratio = float(grid[idx])

    if 0 < idx < len(grid) - 1:
        x = grid[idx - 1 : idx + 2]
        y = curve[idx - 1 : idx + 2]
        try:
            candidate = np.polyfit(x, y, 2)
        except Exception:
            candidate = None
        if candidate is not None:
            a, b, _ = candidate
            if np.isfinite(a) and a < 0.0 and np.isfinite(b):
                vertex = float(-b / (2.0 * a))
                if float(x[0]) <= vertex <= float(x[-1]):
                    coefficients = candidate
                    local_low = float(x[0])
                    local_high = float(x[-1])
                    ratio = vertex

    def evaluate(x_value: float) -> float:
        x_value = float(x_value)
        if coefficients is not None and local_low <= x_value <= local_high:
            return float(np.polyval(coefficients, x_value))
        return float(np.interp(x_value, grid, curve))

    return evaluate, ratio, evaluate(ratio)


def _profile_interval_from_refined_curve(
        grid: np.ndarray, evaluate, ratio: float, peak: float,
        drop: float = PROFILE_DROP_95) -> Tuple[float, float]:
    threshold = float(peak) - float(drop)

    def offset(x_value: float) -> float:
        return evaluate(float(x_value)) - threshold

    knots = sorted(set(float(x) for x in grid) | {float(ratio)})
    ratio_pos = min(range(len(knots)), key=lambda i: abs(knots[i] - ratio))

    if offset(knots[0]) >= 0.0:
        lo = knots[0]
    else:
        lo = math.nan
        inside = float(ratio)
        for j in range(ratio_pos - 1, -1, -1):
            outside = knots[j]
            f_out = offset(outside)
            if f_out <= 0.0:
                if f_out == 0.0:
                    lo = outside
                else:
                    lo = float(brentq(lambda x: offset(x), outside, inside, xtol=1e-12, rtol=1e-12))
                break
            inside = outside

    if offset(knots[-1]) >= 0.0:
        hi = knots[-1]
    else:
        hi = math.nan
        inside = float(ratio)
        for j in range(ratio_pos + 1, len(knots)):
            outside = knots[j]
            f_out = offset(outside)
            if f_out <= 0.0:
                if f_out == 0.0:
                    hi = outside
                else:
                    hi = float(brentq(lambda x: offset(x), inside, outside, xtol=1e-12, rtol=1e-12))
                break
            inside = outside

    lo = max(0.0, min(1.0, lo)) if np.isfinite(lo) else math.nan
    hi = max(0.0, min(1.0, hi)) if np.isfinite(hi) else math.nan
    return lo, hi


def fit_m0(profiles: np.ndarray, grid: np.ndarray) -> M0Result:
    profiles = normalize_profiles(profiles)
    curve = np.sum(profiles, axis=0)
    idx = int(np.argmax(curve))
    evaluate, ratio, ll = _refined_m0_curve(grid, curve, idx)
    lo, hi = _profile_interval_from_refined_curve(grid, evaluate, ratio, ll)
    n = profiles.shape[0]
    bic = -2.0 * ll + math.log(n)
    return M0Result(ratio, lo, hi, ll, bic, idx)


def _bin_edges(grid: np.ndarray) -> np.ndarray:
    edges = np.empty(len(grid) + 1, dtype=float)
    edges[0] = 0.0
    edges[-1] = 1.0
    edges[1:-1] = 0.5 * (grid[:-1] + grid[1:])
    return edges


def truncated_normal_weights(grid: np.ndarray, mu: float, tau: float) -> np.ndarray:
    tau = max(float(tau), 1e-12)
    edges = _bin_edges(grid)
    z0 = (0.0 - mu) / tau
    z1 = (1.0 - mu) / tau
    norm = max(float(ndtr(z1) - ndtr(z0)), EPS)
    masses = np.diff(ndtr((edges - mu) / tau)) / norm
    masses = np.maximum(masses, EPS)
    masses /= masses.sum()
    return masses


def m1_log_likelihood(profiles: np.ndarray, grid: np.ndarray, mu: float, tau: float) -> float:
    weights = truncated_normal_weights(grid, mu, tau)
    return float(np.sum(logsumexp(profiles + np.log(weights)[None, :], axis=1)))


def _distribution_moments(grid: np.ndarray, mu: float, tau: float) -> Tuple[float, float]:
    w = truncated_normal_weights(grid, mu, tau)
    mean = float(np.sum(w * grid))
    var = float(np.sum(w * (grid - mean) ** 2))
    return mean, math.sqrt(max(0.0, var))


def _optimize_mu_for_tau(profiles: np.ndarray, grid: np.ndarray, tau: float, starts: Sequence[float]) -> Tuple[float, float]:
    candidates: List[Tuple[float, float]] = []
    for start in starts:
        span = max(0.1, min(0.5, 4.0 * tau))
        lo = max(0.0, float(start) - span)
        hi = min(1.0, float(start) + span)
        if hi <= lo:
            lo, hi = 0.0, 1.0
        res = minimize_scalar(
            lambda mu: -m1_log_likelihood(profiles, grid, float(mu), tau),
            bounds=(lo, hi), method="bounded", options={"xatol": 1e-7, "maxiter": 300}
        )
        if res.success and np.isfinite(res.fun):
            candidates.append((float(res.x), float(-res.fun)))
    if not candidates:
        return math.nan, -math.inf
    return max(candidates, key=lambda x: x[1])


def fit_m1(profiles: np.ndarray, grid: np.ndarray, m0_ratio: float) -> M1Result:
    profiles = normalize_profiles(profiles)
    step = float(np.min(np.diff(grid)))
    tau_min = max(step / 2.0, 1e-5)
    tau_max = 0.5
    starts = [
        (m0_ratio, tau_min),
        (m0_ratio, max(0.02, tau_min)),
        (m0_ratio, 0.05),
        (m0_ratio, 0.10),
        (0.25, 0.10),
        (0.50, 0.10),
        (0.75, 0.10),
        (0.50, 0.25),
    ]
    bounds = [(0.0, 1.0), (math.log(tau_min), math.log(tau_max))]
    solutions = []
    for mu0, tau0 in starts:
        x0 = np.array([np.clip(mu0, 0.0, 1.0), math.log(np.clip(tau0, tau_min, tau_max))])
        res = minimize(
            lambda x: -m1_log_likelihood(profiles, grid, float(x[0]), math.exp(float(x[1]))),
            x0, method="L-BFGS-B", bounds=bounds,
            options={"ftol": 1e-11, "gtol": 1e-8, "maxiter": 1000},
        )
        if np.isfinite(res.fun):
            solutions.append((float(-res.fun), float(res.x[0]), math.exp(float(res.x[1])), bool(res.success)))
    if not solutions:
        return M1Result(math.nan, math.nan, math.nan, math.nan, math.nan,
                        -math.inf, math.inf, "FAILED")
    best_ll, best_mu, best_tau, best_success = max(solutions, key=lambda x: x[0])
    mean, sd = _distribution_moments(grid, best_mu, best_tau)
    n = profiles.shape[0]
    bic = -2.0 * best_ll + 2.0 * math.log(n)

    # Profile the spread parameter. Each tau re-optimizes location; the reported
    # upper bound is on the actual SD of the bounded discretized distribution.
    tau_grid = np.unique(np.concatenate([
        np.geomspace(tau_min, tau_max, 60),
        np.array([best_tau]),
    ]))
    accepted_sds: List[float] = []
    profile_starts = [best_mu, m0_ratio, 0.25, 0.5, 0.75]
    for tau in tau_grid:
        mu, ll = _optimize_mu_for_tau(profiles, grid, float(tau), profile_starts)
        if np.isfinite(ll) and ll >= best_ll - PROFILE_DROP_95 - 1e-8:
            _, this_sd = _distribution_moments(grid, mu, float(tau))
            accepted_sds.append(this_sd)
    sd_upper = max(accepted_sds) if accepted_sds else math.nan

    successful_lls = [x[0] for x in solutions if x[3]]
    fit_status = "PASS" if best_success or successful_lls else "OPTIMIZER_WARNING"
    if successful_lls and max(successful_lls) - min(successful_lls) > 1e-4:
        # Different starts can land on equivalent boundary parameterizations; only
        # flag material observed-likelihood disagreement.
        near_best = [ll for ll in successful_lls if best_ll - ll < 1.0]
        if near_best and best_ll - min(near_best) > 0.1:
            fit_status = "MULTISTART_DISAGREEMENT"
    return M1Result(best_mu, best_tau, mean, sd, sd_upper, best_ll, bic, fit_status)


def _nearest_idx(grid: np.ndarray, value: float) -> int:
    return int(np.argmin(np.abs(grid - float(value))))


def _m2_observed_ll(profiles: np.ndarray, low_idx: int, high_idx: int, pi_low: float) -> Tuple[float, np.ndarray, np.ndarray]:
    pi_low = float(np.clip(pi_low, 1e-9, 1.0 - 1e-9))
    a = math.log(pi_low) + profiles[:, low_idx]
    b = math.log1p(-pi_low) + profiles[:, high_idx]
    den = np.logaddexp(a, b)
    p_low = np.exp(a - den)
    p_high = 1.0 - p_low
    return float(np.sum(den)), p_low, p_high


def fit_m2(profiles: np.ndarray, grid: np.ndarray, cell_mles: np.ndarray, membership_threshold: float) -> M2Result:
    profiles = normalize_profiles(profiles)
    mles = np.asarray(cell_mles, dtype=float)
    q = np.nanquantile(mles, [0.15, 0.25, 0.50, 0.75, 0.85])
    start_pairs = [
        (q[0], q[4]), (q[1], q[3]), (0.10, 0.90), (0.20, 0.80),
        (0.05, 0.50), (0.50, 0.95), (0.25, 0.75), (q[0], q[2]), (q[2], q[4]),
    ]
    solutions = []
    for low0, high0 in start_pairs:
        low_idx, high_idx = sorted((_nearest_idx(grid, low0), _nearest_idx(grid, high0)))
        if low_idx == high_idx:
            if high_idx < len(grid) - 1:
                high_idx += 1
            elif low_idx > 0:
                low_idx -= 1
            else:
                continue
        pi = 0.5
        prev_ll = -math.inf
        converged = False
        for _ in range(300):
            ll, p_low, p_high = _m2_observed_ll(profiles, low_idx, high_idx, pi)
            pi = float(np.mean(p_low))
            score_low = p_low @ profiles
            score_high = p_high @ profiles
            new_low = int(np.argmax(score_low))
            new_high = int(np.argmax(score_high))
            if new_low > new_high:
                new_low, new_high = new_high, new_low
                pi = 1.0 - pi
            if new_low == new_high:
                # A collapsed point-mass mixture is not an M2 solution.
                break
            if new_low == low_idx and new_high == high_idx and abs(ll - prev_ll) <= 1e-9 * (1.0 + abs(prev_ll)):
                converged = True
                low_idx, high_idx = new_low, new_high
                break
            low_idx, high_idx = new_low, new_high
            prev_ll = ll
        if low_idx == high_idx:
            continue
        ll, p_low, p_high = _m2_observed_ll(profiles, low_idx, high_idx, pi)
        solutions.append((ll, low_idx, high_idx, pi, p_low, p_high, converged))

    if not solutions:
        n = profiles.shape[0]
        return M2Result(math.nan, math.nan, math.nan, math.nan, math.nan, math.nan,
                        math.nan, math.nan, -math.inf, math.inf,
                        np.full(n, np.nan), np.full(n, np.nan), "FAILED")
    converged_solutions = [x for x in solutions if x[-1]]
    candidate_solutions = converged_solutions if converged_solutions else solutions
    ll, low_idx, high_idx, pi, p_low, p_high, converged = max(
        candidate_solutions, key=lambda x: x[0])
    low = float(grid[low_idx])
    high = float(grid[high_idx])
    confident = float(np.mean(np.maximum(p_low, p_high) >= membership_threshold))
    n = profiles.shape[0]
    bic = -2.0 * ll + 3.0 * math.log(n)
    fit_status = "PASS" if converged_solutions else "OPTIMIZER_WARNING"

    if converged_solutions:
        near_best = [
            x for x in converged_solutions
            if ll - x[0] <= M2_NEAR_BEST_LL_TOLERANCE + 1e-12
        ]
        ratio_tolerance = max(0.01, 2.0 * float(np.min(np.diff(grid))))
        near_best_params = []
        for solution in near_best:
            this_low = float(grid[solution[1]])
            this_high = float(grid[solution[2]])
            this_pi = float(solution[3])
            near_best_params.append((
                this_low, this_high, this_pi, 1.0 - this_pi,
                this_high - this_low,
            ))
        for i in range(len(near_best_params)):
            for j in range(i + 1, len(near_best_params)):
                left = near_best_params[i]
                right = near_best_params[j]
                ratio_disagreement = (
                    abs(left[0] - right[0]) > ratio_tolerance or
                    abs(left[1] - right[1]) > ratio_tolerance or
                    abs(left[4] - right[4]) > ratio_tolerance
                )
                proportion_disagreement = (
                    abs(left[2] - right[2]) > M2_PROPORTION_DISAGREEMENT_TOLERANCE or
                    abs(left[3] - right[3]) > M2_PROPORTION_DISAGREEMENT_TOLERANCE
                )
                if ratio_disagreement or proportion_disagreement:
                    fit_status = "MULTISTART_DISAGREEMENT"
                    break
            if fit_status == "MULTISTART_DISAGREEMENT":
                break

    return M2Result(
        low, high, float(pi), float(1.0 - pi), high - low, confident,
        float(np.sum(p_low)), float(np.sum(p_high)), float(ll), float(bic),
        p_low, p_high, fit_status,
    )


def classify_population(m0: M0Result, m1: M1Result, m2: M2Result, n: int,
                        meaningful_sd: float, min_cells: int, delta_bic: float,
                        min_component_fraction: float, min_component_cells: int,
                        min_component_separation: float, min_confident_fraction: float) -> Tuple[str, str, int]:
    if n < min_cells:
        return "INSUFFICIENT_CELLS", f"n_profiled={n}<min_cells={min_cells}", 0

    model_bics = {"M0": m0.bic}
    if m1.fit_status != "FAILED" and np.isfinite(m1.bic):
        model_bics["M1"] = m1.bic
    if m2.fit_status != "FAILED" and np.isfinite(m2.bic):
        model_bics["M2"] = m2.bic
    finite_bics = {name: bic for name, bic in model_bics.items() if np.isfinite(bic)}
    if not finite_bics:
        return "HETEROGENEITY_UNRESOLVED", "no_valid_population_model", 0
    best_bic = min(finite_bics.values())
    competitive = {name for name, bic in finite_bics.items() if bic <= best_bic + delta_bic}
    if "M1" in competitive and "DISAGREEMENT" in m1.fit_status:
        return "HETEROGENEITY_UNRESOLVED", "competitive_M1_multistart_optimization_disagreement", 0
    if "M2" in competitive and "DISAGREEMENT" in m2.fit_status:
        return "HETEROGENEITY_UNRESOLVED", "competitive_M2_multistart_optimization_disagreement", 0

    m2_decisive = (m2.bic <= m0.bic - delta_bic and m2.bic <= m1.bic - delta_bic)
    m2_practical = (
        m2.separation >= min_component_separation and
        min(m2.proportion_low, m2.proportion_high) >= min_component_fraction and
        min(m2.effective_low, m2.effective_high) >= min_component_cells and
        m2.confident_membership_fraction >= min_confident_fraction
    )
    if m2_decisive and m2_practical:
        return "DISCRETE_SUBPOPULATIONS", "M2_decisive_and_components_meet_practical_thresholds", 0

    m1_decisive = m1.bic <= m0.bic - delta_bic
    if m1_decisive and not (m2_decisive and m2_practical) and m1.distribution_sd >= meaningful_sd:
        return "CONTINUOUS_HETEROGENEITY", "M1_decisive_with_meaningful_between_cell_sd", 0

    m1_not_decisive = m1.bic > m0.bic - delta_bic
    m2_not_decisive = m2.bic > m0.bic - delta_bic
    if m1_not_decisive and m2_not_decisive and np.isfinite(m1.sd_upper_95) and m1.sd_upper_95 <= meaningful_sd:
        return "COMMON_RATIO_SUPPORTED", "heterogeneity_not_decisive_and_sd_upper_within_meaningful_threshold", 1

    if m2_decisive and not m2_practical:
        return "HETEROGENEITY_UNRESOLVED", "M2_fit_improves_but_component_requirements_fail", 0
    if m0.bic <= m1.bic and np.isfinite(m1.sd_upper_95) and m1.sd_upper_95 > meaningful_sd:
        return "HETEROGENEITY_UNRESOLVED", "M0_not_rejected_but_meaningful_heterogeneity_not_excluded", 0
    return "HETEROGENEITY_UNRESOLVED", "model_evidence_or_practical_effect_not_decisive", 0


def analyze_profiles(profiles: np.ndarray, grid: np.ndarray, cell_mles: np.ndarray,
                     meaningful_sd: float = 0.05, min_cells: int = 20,
                     delta_bic: float = 10.0, min_component_fraction: float = 0.10,
                     min_component_cells: int = 10, min_component_separation: float = 0.10,
                     membership_threshold: float = 0.80,
                     min_confident_fraction: float = 0.70) -> Dict[str, object]:
    profiles = normalize_profiles(profiles)
    m0 = fit_m0(profiles, grid)
    m1 = fit_m1(profiles, grid, m0.ratio)
    m2 = fit_m2(profiles, grid, cell_mles, membership_threshold)
    structure, reason, pooling = classify_population(
        m0, m1, m2, profiles.shape[0], meaningful_sd, min_cells, delta_bic,
        min_component_fraction, min_component_cells, min_component_separation,
        min_confident_fraction,
    )
    return {"m0": m0, "m1": m1, "m2": m2,
            "population_structure": structure,
            "population_structure_reason": reason,
            "pooling_supported": pooling}


def _read_tables(paths: Sequence[str]) -> pd.DataFrame:
    frames = []
    for path in paths:
        df = pd.read_csv(path, sep="\t", dtype={"library_id": str, "barcode": str})
        df["_source_path"] = path
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _read_reconciled_cells(paths: Sequence[str]) -> pd.DataFrame:
    frames = []
    for path in paths:
        df = pd.read_csv(path, sep="\t", dtype=str)
        if "library_id" not in df.columns and "library" in df.columns:
            df = df.rename(columns={"library": "library_id"})
        df["_source_path"] = path
        frames.append(df)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def _canonical_donor_genotype(value: object) -> str:
    text = "" if value is None or pd.isna(value) else str(value).strip()
    if not text or text.upper() in {"NA", "NAN", "NONE", "."}:
        return ""
    parts = [x.strip() for x in text.split("+") if x.strip()]
    return "+".join(sorted(parts))


def _read_profile_rows(paths: Sequence[str]) -> pd.DataFrame:
    frames = []
    for path in paths:
        df = pd.read_csv(path, sep="\t", dtype={"library_id": str, "barcode": str})
        df["_source_path"] = path
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _write_exclusions(prefix: str, exclusions: List[Dict[str, str]]) -> None:
    path = prefix + ".mt_population_exclusions.tsv"
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    cols = [
        "library_id", "barcode", "identity", "canonical_parent1",
        "canonical_parent2", "group_id", "exclusion_reason", "detail",
    ]
    pd.DataFrame(exclusions, columns=cols).to_csv(path, sep="\t", index=False, na_rep="NA")


def _enrich_exclusions(exclusions: List[Dict[str, str]], ratio: pd.DataFrame,
                       ratio_meta: pd.DataFrame, group_cols: Sequence[str]) -> None:
    context: Dict[Tuple[str, str], Dict[str, str]] = {}

    def add_frame(frame: pd.DataFrame, include_group: bool) -> None:
        if frame.empty or "library_id" not in frame.columns or "barcode" not in frame.columns:
            return
        can_group = include_group and all(col in frame.columns for col in group_cols)
        for _, row in frame.iterrows():
            key = (str(row.get("library_id", "NA")), str(row.get("barcode", "NA")))
            entry = context.setdefault(key, {})
            for col in ("identity", "canonical_parent1", "canonical_parent2"):
                value = row.get(col, "NA")
                if pd.notna(value) and str(value).strip() not in {"", "NA", "nan"}:
                    entry[col] = str(value)
            if can_group:
                values = [row.get(col, "NA") for col in group_cols]
                entry["group_id"] = _group_id(group_cols, values)

    add_frame(ratio, False)
    add_frame(ratio_meta, True)
    for exclusion in exclusions:
        key = (str(exclusion.get("library_id", "NA")), str(exclusion.get("barcode", "NA")))
        entry = context.get(key, {})
        for col in ("identity", "canonical_parent1", "canonical_parent2", "group_id"):
            exclusion.setdefault(col, entry.get(col, "NA"))


def _check_duplicates(df: pd.DataFrame, label: str, exclusions: List[Dict[str, str]]) -> bool:
    if df.empty:
        return False
    dup = df.duplicated(["library_id", "barcode"], keep=False)
    if not dup.any():
        return False
    for _, row in df.loc[dup, ["library_id", "barcode"]].drop_duplicates().iterrows():
        exclusions.append({"library_id": str(row.library_id), "barcode": str(row.barcode),
                           "exclusion_reason": "DUPLICATE_INPUT_KEY", "detail": label})
    return True


def _parse_profile_vector(row: pd.Series) -> Tuple[np.ndarray, np.ndarray]:
    start = float(row["grid_start"])
    step = float(row["grid_step"])
    points = int(row["grid_points"])
    tokens = str(row["delta_log_likelihood_csv"]).split(",")
    if len(tokens) != points:
        raise ValueError(f"grid_points={points} but vector has {len(tokens)} values")
    vals = np.array([math.nan if x == "NA" else float(x) for x in tokens], dtype=float)
    grid = start + step * np.arange(points, dtype=float)
    if points:
        if abs(grid[0]) < 1e-8:
            grid[0] = 0.0
        if abs(grid[-1] - 1.0) < max(1e-8, step * 1e-4):
            grid[-1] = 1.0
    return grid, vals


def _interpolate_profiles(profile_df: pd.DataFrame, exclusions: List[Dict[str, str]]) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    parsed = []
    steps = []
    for idx, row in profile_df.iterrows():
        key = (str(row.library_id), str(row.barcode))
        try:
            grid, vals = _parse_profile_vector(row)
        except Exception as exc:
            exclusions.append({"library_id": key[0], "barcode": key[1],
                               "exclusion_reason": "NONFINITE_PROFILE", "detail": str(exc)})
            continue
        finite = np.isfinite(vals)
        if len(grid) < 2 or not np.isclose(grid[0], 0.0) or not np.isclose(grid[-1], 1.0):
            exclusions.append({"library_id": key[0], "barcode": key[1],
                               "exclusion_reason": "NONFINITE_PROFILE", "detail": "profile grid does not span [0,1]"})
            continue
        if finite.sum() < 2 or not finite[0] or not finite[-1]:
            exclusions.append({"library_id": key[0], "barcode": key[1],
                               "exclusion_reason": "NONFINITE_PROFILE", "detail": "profile lacks finite boundary support"})
            continue
        if not finite.all():
            vals = np.interp(grid, grid[finite], vals[finite])
        vals = vals - np.max(vals)
        parsed.append((idx, grid, vals))
        steps.append(float(np.min(np.diff(grid))))
    if not parsed:
        return profile_df.iloc[0:0].copy(), np.array([]), np.empty((0, 0))
    finest = min(steps)
    n_intervals = max(1, int(round(1.0 / finest)))
    common_grid = np.linspace(0.0, 1.0, n_intervals + 1)
    rows = []
    matrices = []
    for idx, grid, vals in parsed:
        interp = np.interp(common_grid, grid, vals)
        interp -= np.max(interp)
        rows.append(profile_df.loc[idx])
        matrices.append(interp)
    return pd.DataFrame(rows).reset_index(drop=True), common_grid, np.vstack(matrices)


def _canonical_raw_support(row: pd.Series) -> Tuple[float, float]:
    p1 = pd.to_numeric(row.get("raw_parent1_support"), errors="coerce")
    p2 = pd.to_numeric(row.get("raw_parent2_support"), errors="coerce")
    match = str(row.get("assignment_matches_canonical", "NA"))
    if not (np.isfinite(p1) and np.isfinite(p2)):
        return math.nan, math.nan
    if match in {"1", "1.0", "True", "true"}:
        return float(p1), float(p2)
    if match in {"0", "0.0", "False", "false"}:
        return float(p2), float(p1)
    return math.nan, math.nan


def _group_id(group_cols: Sequence[str], values: Sequence[object]) -> str:
    return "|".join(f"{col}={str(value)}" for col, value in zip(group_cols, values))


def _selected_model_summary(structure: str, m0: M0Result, m1: M1Result, m2: M2Result) -> Tuple[float, float]:
    if structure == "DISCRETE_SUBPOPULATIONS":
        mean = m2.proportion_low * m2.ratio_low + m2.proportion_high * m2.ratio_high
        var = (m2.proportion_low * (m2.ratio_low - mean) ** 2 +
               m2.proportion_high * (m2.ratio_high - mean) ** 2)
        return float(mean), math.sqrt(max(0.0, var))
    if structure == "CONTINUOUS_HETEROGENEITY":
        return m1.distribution_mean, m1.distribution_sd
    if structure == "COMMON_RATIO_SUPPORTED":
        return m0.ratio, 0.0
    # No population model is selected for unresolved or insufficient groups.
    # Keep the explicit M0/M1/M2 and descriptive summaries available, but do
    # not turn the lowest BIC into an implied biological selection.
    return math.nan, math.nan


def _group_calibration_status(group: pd.DataFrame) -> str:
    if "site_calibration_status" not in group.columns:
        return "UNCALIBRATED"
    vals = {str(x).strip().upper() for x in group["site_calibration_status"].dropna()}
    vals -= {"", "NA", ".", "NAN", "NOT_REQUESTED", "NO_USED_SITES"}
    if not vals:
        return "UNCALIBRATED"
    if vals == {"FULL"}:
        return "CALIBRATED"
    if vals == {"FALLBACK"}:
        return "CALIBRATION_REQUESTED_NO_MATCH"
    return "PARTIAL_CALIBRATION"


def _write_calibration_comparison(args: argparse.Namespace, groups_path: str, group_cols: Sequence[str]) -> None:
    if not getattr(args, "compare_to_groups", None):
        return
    current = pd.read_csv(groups_path, sep="\t", dtype=str)
    reference = pd.read_csv(getattr(args, "compare_to_groups"), sep="\t", dtype=str)
    join_cols = [c for c in group_cols if c in current.columns and c in reference.columns]
    if len(join_cols) != len(group_cols):
        missing = [c for c in group_cols if c not in current.columns or c not in reference.columns]
        raise ValueError("calibration comparison missing group column(s): " + ", ".join(missing))
    keep = join_cols + ["population_structure", "calibration_status", "common_ratio_mle", "selected_model_population_mean"]
    left = reference[[c for c in keep if c in reference.columns]].copy()
    right = current[[c for c in keep if c in current.columns]].copy()
    comp = left.merge(right, on=join_cols, how="outer", suffixes=("_reference", "_current"), indicator=True)
    comp["structure_changed"] = np.where(
        comp["_merge"] == "both",
        (comp.get("population_structure_reference") != comp.get("population_structure_current")).astype(int),
        np.nan,
    )
    comp["comparison_status"] = comp["_merge"].map({"both": "MATCHED", "left_only": "REFERENCE_ONLY", "right_only": "CURRENT_ONLY"})
    comp.drop(columns=["_merge"]).to_csv(args.output_prefix + ".mt_population_calibration_comparison.tsv", sep="\t", index=False, na_rep="NA")


def run_population(args: argparse.Namespace) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(args.output_prefix)), exist_ok=True)
    exclusions: List[Dict[str, str]] = []
    ratio = _read_tables(args.ratio_tsv)
    profile = _read_profile_rows(args.profile_tsv)
    reconciled_mode = bool(getattr(args, "reconciled_cells", None))
    if reconciled_mode:
        metadata = _read_reconciled_cells(args.reconciled_cells)
        metadata_label = "reconciled-cell input"
    else:
        metadata = pd.read_csv(
            args.metadata_tsv, sep="\t", dtype={"library_id": str, "barcode": str})
        metadata_label = "metadata input"

    required_ratio = [
        "library_id", "barcode", "canonical_parent1", "canonical_parent2",
        "canonical_parent2_fraction", "parent2_profile_ci_low",
        "parent2_profile_ci_high", "ratio_molecules", "ratio_sites_used",
        "status", "profile_status",
        "assignment_matches_canonical", "raw_parent1_support", "raw_parent2_support",
    ]
    required_profile = [
        "library_id", "barcode", "canonical_parent1", "canonical_parent2",
        "canonical_parent2_fraction_mle", "grid_start", "grid_step",
        "grid_points", "delta_log_likelihood_csv", "profile_status",
    ]
    missing_ratio = [x for x in required_ratio if x not in ratio.columns]
    missing_profile = [x for x in required_profile if x not in profile.columns]
    if missing_ratio or missing_profile:
        details = []
        if missing_ratio:
            details.append("ratio input columns missing: " + ", ".join(missing_ratio))
        if missing_profile:
            details.append("profile input columns missing: " + ", ".join(missing_profile))
        _write_exclusions(args.output_prefix, exclusions)
        raise ValueError("; ".join(details))

    group_cols = [x.strip() for x in args.group_by.split(",") if x.strip()]
    if reconciled_mode:
        required_meta = [
            "library_id", "barcode", "reconciled_donor_genotype",
            "reconciled_droplet_state", "reconciled_uid",
            "uid_candidate_count", "uid_resolution_status",
        ]
    else:
        required_meta = ["library_id", "barcode", "parent_pair_id", "fusion_replicate_id"]
    missing_required = [x for x in required_meta if x not in metadata.columns]
    if missing_required:
        _write_exclusions(args.output_prefix, exclusions)
        raise ValueError(
            f"required {metadata_label} columns missing: " + ", ".join(missing_required))

    duplicate_flags = [
        _check_duplicates(ratio, "ratio input", exclusions),
        _check_duplicates(profile, "profile input", exclusions),
        _check_duplicates(metadata, metadata_label, exclusions),
    ]
    if any(duplicate_flags):
        _write_exclusions(args.output_prefix, exclusions)
        raise ValueError("duplicate library_id+barcode keys detected; see exclusions table")

    for frame in (ratio, profile, metadata):
        frame["library_id"] = frame["library_id"].astype(str)
        frame["barcode"] = frame["barcode"].astype(str)

    ratio_meta = ratio.merge(
        metadata, on=["library_id", "barcode"], how="left", indicator=True,
        suffixes=("", "_meta"))
    missing_meta = ratio_meta["_merge"] == "left_only"
    for _, row in ratio_meta.loc[missing_meta, ["library_id", "barcode"]].iterrows():
        exclusions.append({
            "library_id": str(row.library_id), "barcode": str(row.barcode),
            "exclusion_reason": "MISSING_RECONCILED_CELL" if reconciled_mode else "MISSING_METADATA",
            "detail": "no matching reconciled-cell row" if reconciled_mode else "no metadata row",
        })
    ratio_meta = ratio_meta.loc[~missing_meta].drop(columns="_merge").copy()

    if reconciled_mode:
        # Parent-pair identity belongs to the mitochondrial ratio output itself.
        # Physical fusion replicate identity belongs to the reconciliation layer.
        ratio_meta["parent_pair_id"] = (
            ratio_meta["canonical_parent1"].astype(str) + "+" +
            ratio_meta["canonical_parent2"].astype(str))
        uid_count = pd.to_numeric(ratio_meta["uid_candidate_count"], errors="coerce")
        uid_text = ratio_meta["reconciled_uid"].fillna("").astype(str).str.strip()
        status_text = ratio_meta["uid_resolution_status"].fillna("").astype(str).str.strip()
        droplet_text = ratio_meta["reconciled_droplet_state"].fillna("").astype(str).str.strip()

        keep = np.ones(len(ratio_meta), dtype=bool)
        for pos, (_, row) in enumerate(ratio_meta.iterrows()):
            lib = str(row["library_id"])
            bc = str(row["barcode"])
            n_uid = pd.to_numeric(pd.Series([row.get("uid_candidate_count")]), errors="coerce").iloc[0]
            uid = "" if pd.isna(row.get("reconciled_uid")) else str(row.get("reconciled_uid")).strip()
            uid_status = "" if pd.isna(row.get("uid_resolution_status")) else str(row.get("uid_resolution_status")).strip()
            droplet = "" if pd.isna(row.get("reconciled_droplet_state")) else str(row.get("reconciled_droplet_state")).strip()
            expected_pair = _canonical_donor_genotype(
                f"{row.get('canonical_parent1', '')}+{row.get('canonical_parent2', '')}")
            reconciled_pair = _canonical_donor_genotype(row.get("reconciled_donor_genotype"))
            if droplet != "SINGLE_CELL":
                keep[pos] = False
                exclusions.append({"library_id": lib, "barcode": bc,
                    "exclusion_reason": "RECONCILED_NOT_SINGLE_CELL",
                    "detail": f"reconciled_droplet_state={droplet or 'NA'}"})
            elif reconciled_pair != expected_pair:
                keep[pos] = False
                exclusions.append({"library_id": lib, "barcode": bc,
                    "exclusion_reason": "RECONCILED_IDENTITY_MISMATCH",
                    "detail": f"mt_pair={expected_pair}; reconciled_donor_genotype={reconciled_pair or 'NA'}"})
            elif not np.isfinite(n_uid) or int(n_uid) < 1 or not uid:
                keep[pos] = False
                exclusions.append({"library_id": lib, "barcode": bc,
                    "exclusion_reason": "FUSION_REPLICATE_UNRESOLVED",
                    "detail": f"uid_resolution_status={uid_status or 'NA'}"})
            elif int(n_uid) != 1:
                # A set-valued UID is valid reconciliation metadata, but it does
                # not identify one independent fusion replicate for population pooling.
                keep[pos] = False
                exclusions.append({"library_id": lib, "barcode": bc,
                    "exclusion_reason": "MULTIPLE_FUSION_REPLICATES_UNRESOLVED",
                    "detail": f"reconciled_uid={uid}; uid_candidate_count={int(n_uid)}"})
        ratio_meta = ratio_meta.loc[keep].copy()
        ratio_meta["fusion_replicate_id"] = ratio_meta["reconciled_uid"].astype(str)
        # assay_mode remains owned by the C++ result; do not overwrite it from metadata.
    else:
        # Legacy external metadata remains supported for standalone use.
        for group_col in group_cols:
            metadata_col = group_col + "_meta"
            if metadata_col in ratio_meta.columns:
                ratio_meta[group_col] = ratio_meta[metadata_col]

    missing_group = [x for x in group_cols if x not in ratio_meta.columns]
    if missing_group:
        _write_exclusions(args.output_prefix, exclusions)
        raise ValueError("requested group columns missing after input join: " + ", ".join(missing_group))

    # Validate profile orientation against ratio output before interpolation.
    ratio_orientation = ratio[["library_id", "barcode", "canonical_parent1", "canonical_parent2"]].copy()
    profile_columns = list(profile.columns)
    prof_join = profile.merge(ratio_orientation, on=["library_id", "barcode"], how="left",
                              suffixes=("_profile", "_ratio"), indicator=True)
    keep_profile = np.ones(len(prof_join), dtype=bool)
    for i, row in prof_join.iterrows():
        if row["_merge"] == "left_only":
            keep_profile[i] = False
            continue
        if (str(row.get("canonical_parent1_profile")) != str(row.get("canonical_parent1_ratio")) or
                str(row.get("canonical_parent2_profile")) != str(row.get("canonical_parent2_ratio"))):
            keep_profile[i] = False
            exclusions.append({"library_id": str(row.library_id), "barcode": str(row.barcode),
                               "exclusion_reason": "PARENT_ORIENTATION_MISMATCH",
                               "detail": "ratio/profile canonical parents differ"})
        if str(row.get("profile_status", "")) == "FAILED":
            keep_profile[i] = False
            exclusions.append({"library_id": str(row.library_id), "barcode": str(row.barcode),
                               "exclusion_reason": "PROFILE_FAILED", "detail": "C++ profile_status=FAILED"})
    kept = prof_join.loc[keep_profile].copy()
    reconstructed = {}
    for col in profile_columns:
        source_col = col
        if col in {"canonical_parent1", "canonical_parent2"}:
            source_col = col + "_profile"
        reconstructed[col] = kept[source_col].to_numpy()
    profile = pd.DataFrame(reconstructed).reset_index(drop=True)

    valid_profile, common_grid, matrix = _interpolate_profiles(profile, exclusions)
    matrix_by_key = {(str(r.library_id), str(r.barcode)): matrix[i]
                     for i, r in valid_profile.iterrows()}

    ratio_keys = set(zip(ratio_meta.library_id.astype(str), ratio_meta.barcode.astype(str)))
    profile_keys = set(matrix_by_key)
    for lib, bc in sorted(ratio_keys - profile_keys):
        if not any(e["library_id"] == lib and e["barcode"] == bc for e in exclusions):
            exclusions.append({"library_id": lib, "barcode": bc,
                               "exclusion_reason": "MISSING_PROFILE", "detail": "no valid profile row"})

    raw_support = ratio_meta.apply(_canonical_raw_support, axis=1)
    ratio_meta["_canonical_raw_parent1"] = [x[0] for x in raw_support]
    ratio_meta["_canonical_raw_parent2"] = [x[1] for x in raw_support]
    ratio_meta["canonical_parent2_fraction"] = pd.to_numeric(
        ratio_meta.get("canonical_parent2_fraction"), errors="coerce")

    group_rows = []
    cell_rows = []
    grouped = ratio_meta.groupby(group_cols, dropna=False, sort=True)
    for group_values, group in grouped:
        if not isinstance(group_values, tuple):
            group_values = (group_values,)
        gid = _group_id(group_cols, group_values)
        keys = list(zip(group.library_id.astype(str), group.barcode.astype(str)))
        modeled_positions = [i for i, key in enumerate(keys) if key in matrix_by_key]
        profiled = group.iloc[modeled_positions].copy()
        profiles = np.vstack([matrix_by_key[keys[i]] for i in modeled_positions]) if modeled_positions else np.empty((0, len(common_grid)))
        mles = pd.to_numeric(profiled.get("canonical_parent2_fraction"), errors="coerce").to_numpy(dtype=float) if len(profiled) else np.array([])
        finite_mle = np.isfinite(mles)
        if len(profiled) and not finite_mle.all():
            # A profile can still be modeled without an MLE for M0/M1, but M2 starts
            # need finite cell locations. Use each profile's grid maximum instead.
            mles = np.where(finite_mle, mles, common_grid[np.argmax(profiles, axis=1)])

        if len(profiled):
            result = analyze_profiles(
                profiles, common_grid, mles,
                meaningful_sd=args.meaningful_sd,
                min_cells=args.min_cells,
                delta_bic=args.delta_bic,
                min_component_fraction=args.min_component_fraction,
                min_component_cells=args.min_component_cells,
                min_component_separation=args.min_component_separation,
                membership_threshold=args.membership_threshold,
                min_confident_fraction=args.min_confident_fraction,
            )
            m0, m1, m2 = result["m0"], result["m1"], result["m2"]
            structure = result["population_structure"]
            reason = result["population_structure_reason"]
            pooling = result["pooling_supported"]
        else:
            m0 = M0Result(math.nan, math.nan, math.nan, -math.inf, math.inf, -1)
            m1 = M1Result(math.nan, math.nan, math.nan, math.nan, math.nan, -math.inf, math.inf, "FAILED")
            m2 = M2Result(math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan,
                          -math.inf, math.inf, np.array([]), np.array([]), "FAILED")
            structure, reason, pooling = "INSUFFICIENT_CELLS", "no_valid_profiles", 0

        finite_ratio = pd.to_numeric(group["canonical_parent2_fraction"], errors="coerce").dropna().to_numpy(dtype=float)
        raw1 = pd.to_numeric(group["_canonical_raw_parent1"], errors="coerce").to_numpy(dtype=float)
        raw2 = pd.to_numeric(group["_canonical_raw_parent2"], errors="coerce").to_numpy(dtype=float)
        raw_ok = np.isfinite(raw1) & np.isfinite(raw2)
        denom = float(np.sum(raw1[raw_ok] + raw2[raw_ok])) if raw_ok.any() else 0.0
        raw_fraction = float(np.sum(raw2[raw_ok]) / denom) if denom > 0 else math.nan
        selected_mean, selected_sd = _selected_model_summary(structure, m0, m1, m2)
        model_fit_status = f"M1:{m1.fit_status};M2:{m2.fit_status}"

        row = {
            "group_id": gid,
            "population_implementation_version": POPULATION_IMPLEMENTATION_VERSION,
            **{c: v for c, v in zip(group_cols, group_values)},
            "n_cells_input": int(len(group)),
            "n_cells_with_ratio": int(np.isfinite(pd.to_numeric(group["canonical_parent2_fraction"], errors="coerce")).sum()),
            "n_cells_profiled": int(len(profiled)),
            "n_cells_excluded": int(len(group) - len(profiled)),
            "population_structure": structure,
            "population_structure_reason": reason,
            "pooling_supported": int(pooling),
            "meaningful_sd_threshold": args.meaningful_sd,
            "sd_upper_95": m1.sd_upper_95,
            "common_ratio_mle": m0.ratio,
            "common_ratio_ci_low": m0.ci_low,
            "common_ratio_ci_high": m0.ci_high,
            "selected_model_population_mean": selected_mean,
            "selected_model_population_sd": selected_sd,
            "cell_weighted_mean_mle": float(np.mean(finite_ratio)) if len(finite_ratio) else math.nan,
            "cell_weighted_median_mle": float(np.median(finite_ratio)) if len(finite_ratio) else math.nan,
            "cell_weighted_sd_mle": float(np.std(finite_ratio, ddof=1)) if len(finite_ratio) > 1 else (0.0 if len(finite_ratio) == 1 else math.nan),
            "molecule_weighted_raw_parent2_fraction": raw_fraction,
            "m0_log_likelihood": m0.log_likelihood,
            "m0_bic": m0.bic,
            "m1_location": m1.location,
            "m1_distribution_mean": m1.distribution_mean,
            "m1_distribution_sd": m1.distribution_sd,
            "m1_log_likelihood": m1.log_likelihood,
            "m1_bic": m1.bic,
            "m2_ratio_low": m2.ratio_low,
            "m2_ratio_high": m2.ratio_high,
            "m2_proportion_low": m2.proportion_low,
            "m2_proportion_high": m2.proportion_high,
            "m2_separation": m2.separation,
            "m2_confident_membership_fraction": m2.confident_membership_fraction,
            "m2_log_likelihood": m2.log_likelihood,
            "m2_bic": m2.bic,
            "delta_bic_m1_vs_m0": m0.bic - m1.bic,
            "delta_bic_m2_vs_m0": m0.bic - m2.bic,
            "delta_bic_m2_vs_m1": m1.bic - m2.bic,
            "model_fit_status": model_fit_status,
            "calibration_status": _group_calibration_status(profiled if len(profiled) else group),
            "delta_bic_threshold": args.delta_bic,
            "min_component_fraction": args.min_component_fraction,
            "min_component_cells": args.min_component_cells,
            "min_component_separation": args.min_component_separation,
            "membership_threshold": args.membership_threshold,
            "min_confident_fraction": args.min_confident_fraction,
            "min_cells_threshold": args.min_cells,
        }
        group_rows.append(row)

        if len(profiled):
            for local_i, (_, cell) in enumerate(profiled.iterrows()):
                mle = float(pd.to_numeric(cell.get("canonical_parent2_fraction"), errors="coerce"))
                p_low = float(m2.p_low[local_i]) if len(m2.p_low) else math.nan
                p_high = float(m2.p_high[local_i]) if len(m2.p_high) else math.nan
                assignment = "NA"
                assignment_status = "NOT_SELECTED"
                if structure == "DISCRETE_SUBPOPULATIONS":
                    if max(p_low, p_high) >= args.membership_threshold:
                        assignment = "LOW" if p_low >= p_high else "HIGH"
                        assignment_status = "CONFIDENT"
                    else:
                        assignment = "AMBIGUOUS"
                        assignment_status = "BELOW_MEMBERSHIP_THRESHOLD"
                cell_rows.append({
                    "library_id": str(cell.library_id),
                    "barcode": str(cell.barcode),
                    "group_id": gid,
                    "canonical_parent1": cell.get("canonical_parent1", "NA"),
                    "canonical_parent2": cell.get("canonical_parent2", "NA"),
                    "canonical_parent2_fraction_mle": mle,
                    "parent2_profile_ci_low": pd.to_numeric(cell.get("parent2_profile_ci_low"), errors="coerce"),
                    "parent2_profile_ci_high": pd.to_numeric(cell.get("parent2_profile_ci_high"), errors="coerce"),
                    "ratio_molecules": pd.to_numeric(cell.get("ratio_molecules"), errors="coerce"),
                    "ratio_sites_used": pd.to_numeric(cell.get("ratio_sites_used"), errors="coerce"),
                    "primary_fit_status": cell.get("status", "NA"),
                    "profile_status": cell.get("profile_status", "NA"),
                    "population_structure": structure,
                    "common_ratio_residual": mle - m0.ratio if np.isfinite(mle) and np.isfinite(m0.ratio) else math.nan,
                    "p_component_low": p_low,
                    "p_component_high": p_high,
                    "component_assignment": assignment,
                    "component_assignment_status": assignment_status,
                })

    groups_path = args.output_prefix + ".mt_population_groups.tsv"
    pd.DataFrame(group_rows).to_csv(groups_path, sep="\t", index=False, na_rep="NA")
    cell_columns = [
        "library_id", "barcode", "group_id", "canonical_parent1", "canonical_parent2",
        "canonical_parent2_fraction_mle", "parent2_profile_ci_low", "parent2_profile_ci_high",
        "ratio_molecules", "ratio_sites_used", "primary_fit_status", "profile_status",
        "population_structure", "common_ratio_residual", "p_component_low", "p_component_high",
        "component_assignment", "component_assignment_status",
    ]
    pd.DataFrame(cell_rows, columns=cell_columns).to_csv(
        args.output_prefix + ".mt_population_cells.tsv", sep="\t", index=False, na_rep="NA")
    _enrich_exclusions(exclusions, ratio, ratio_meta, group_cols)
    _write_exclusions(args.output_prefix, exclusions)

    revisions = sorted(set(str(x) for x in ratio.get("source_revision", pd.Series(dtype=str)).dropna()))
    versions = sorted(set(str(x) for x in ratio.get("model_version", pd.Series(dtype=str)).dropna()))
    params = [
        ("ratio_tsv", ";".join(args.ratio_tsv)),
        ("profile_tsv", ";".join(args.profile_tsv)),
        ("metadata_source", "reconciled_cells" if reconciled_mode else "metadata_tsv"),
        ("reconciled_cells", ";".join(args.reconciled_cells) if reconciled_mode else "NA"),
        ("metadata_tsv", args.metadata_tsv if not reconciled_mode else "NA"),
        ("output_prefix", args.output_prefix),
        ("group_by", ",".join(group_cols)),
        ("population_implementation_version", POPULATION_IMPLEMENTATION_VERSION),
        ("source_revisions", ",".join(revisions)),
        ("model_versions", ",".join(versions)),
        ("profile_grid_step", float(common_grid[1] - common_grid[0]) if len(common_grid) > 1 else math.nan),
        ("meaningful_sd", args.meaningful_sd),
        ("min_cells", args.min_cells),
        ("delta_bic", args.delta_bic),
        ("min_component_fraction", args.min_component_fraction),
        ("min_component_cells", args.min_component_cells),
        ("min_component_separation", args.min_component_separation),
        ("membership_threshold", args.membership_threshold),
        ("min_confident_fraction", args.min_confident_fraction),
        ("calibration_status", ";".join(sorted(set(str(x) for x in pd.DataFrame(group_rows).get("calibration_status", pd.Series(dtype=str)).dropna()))) or "UNCALIBRATED"),
        ("compare_to_groups", getattr(args, "compare_to_groups", None) or "NA"),
        ("python_version", platform.python_version()),
        ("numpy_version", np.__version__),
        ("pandas_version", pd.__version__),
        ("scipy_version", scipy.__version__),
        ("run_timestamp_utc", datetime.now(timezone.utc).isoformat()),
    ]
    pd.DataFrame(params, columns=["parameter", "value"]).to_csv(
        args.output_prefix + ".mt_population_run_parameters.tsv", sep="\t", index=False)
    _write_calibration_comparison(args, groups_path, group_cols)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Infer within-replicate mitochondrial population structure from per-cell profile likelihoods")
    p.add_argument("--ratio-tsv", nargs="+", required=True)
    p.add_argument("--profile-tsv", nargs="+", required=True)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--reconciled-cells", nargs="+",
        help=("Rich identity-reconciliation cell table(s). The MT population layer derives "
              "parent_pair_id from the ratio output and uses a single resolved reconciled_uid "
              "as fusion_replicate_id."))
    source.add_argument(
        "--metadata-tsv",
        help="Legacy standalone cell metadata TSV with library_id/barcode/parent_pair_id/fusion_replicate_id")
    p.add_argument("--group-by", default="assay_mode,parent_pair_id,fusion_replicate_id")
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--meaningful-sd", type=float, default=0.05)
    p.add_argument("--min-cells", type=int, default=20)
    p.add_argument("--delta-bic", type=float, default=10.0)
    p.add_argument("--min-component-fraction", type=float, default=0.10)
    p.add_argument("--min-component-cells", type=int, default=10)
    p.add_argument("--min-component-separation", type=float, default=0.10)
    p.add_argument("--membership-threshold", type=float, default=0.80)
    p.add_argument("--min-confident-fraction", type=float, default=0.70)
    p.add_argument("--compare-to-groups", default=None,
                   help="Optional previous .mt_population_groups.tsv for calibrated-vs-reference call comparison")
    return p


def _validate_args(args: argparse.Namespace) -> None:
    if args.meaningful_sd <= 0 or args.meaningful_sd >= 0.5:
        raise ValueError("--meaningful-sd must be between 0 and 0.5")
    if args.min_cells < 1 or args.min_component_cells < 1:
        raise ValueError("cell-count thresholds must be positive")
    if args.delta_bic < 0:
        raise ValueError("--delta-bic must be nonnegative")
    for name in ("min_component_fraction", "membership_threshold", "min_confident_fraction"):
        value = getattr(args, name)
        if not 0 < value <= 1:
            raise ValueError(f"--{name.replace('_', '-')} must be in (0,1]")
    if not 0 < args.min_component_separation <= 1:
        raise ValueError("--min-component-separation must be in (0,1]")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        _validate_args(args)
        run_population(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
