#!/usr/bin/env python3
"""Pure helpers and contracts for the CellBouncer arm-CNV hybrid v1 model.

Hybrid v1 keeps expression and ASE as independently calibrated branches.  This
module intentionally contains no BAM, pileup, or demultiplexing code and never
forms a multiplied expression-plus-ASE posterior.
"""

from __future__ import annotations

import bisect
import hashlib
import math
import statistics
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import numpy as np

from tetra_arm_common import (
    EXPRESSION_MODEL_SCHEMA,
    HYBRID_CALIBRATION_SCHEMA,
    HYBRID_CALL_SCHEMA,
    HYBRID_CELL_MANIFEST_SCHEMA,
    HYBRID_COMPONENT_SCHEMA,
    HYBRID_CONTRACT_SCHEMA,
    HYBRID_PAIR_SCHEMA,
    HYBRID_QC_SCHEMA,
    HYBRID_SCORE_SCHEMA,
    HYBRID_SHARD_SCHEMA,
    HYBRID_UID_SCHEMA,
    clean,
    finite_float,
    natural_key,
    validate_tsv_schema,
)


PROGRAM_VERSION = "1.0.11"

EXACT_STATES = (
    "DONOR_A_LOSS",
    "DONOR_B_LOSS",
    "DONOR_A_GAIN",
    "DONOR_B_GAIN",
)
ASE_DIRECTIONS = ("DONOR_A_DEPLETED", "DONOR_A_ENRICHED")
COPY_STATES = ("LOSS", "GAIN")

STATE_LOG_ODDS_SHIFT = {
    "BALANCED": 0.0,
    "DONOR_A_LOSS": math.log(1.0 / 2.0),
    "DONOR_B_LOSS": math.log(2.0),
    "DONOR_A_GAIN": math.log(3.0 / 2.0),
    "DONOR_B_GAIN": math.log(2.0 / 3.0),
}
STATE_AXES = {
    "DONOR_A_LOSS": ("LOSS", "DONOR_A_DEPLETED", "DONOR_A"),
    "DONOR_B_LOSS": ("LOSS", "DONOR_A_ENRICHED", "DONOR_B"),
    "DONOR_A_GAIN": ("GAIN", "DONOR_A_ENRICHED", "DONOR_A"),
    "DONOR_B_GAIN": ("GAIN", "DONOR_A_DEPLETED", "DONOR_B"),
}
AXES_TO_STATE = {
    (copy_state, ase_direction): state
    for state, (copy_state, ase_direction, _donor) in STATE_AXES.items()
}
ASE_DIRECTION_STATES = {
    "DONOR_A_DEPLETED": ("DONOR_A_LOSS", "DONOR_B_GAIN"),
    "DONOR_A_ENRICHED": ("DONOR_B_LOSS", "DONOR_A_GAIN"),
}

LOSS_SHIFT = math.log(0.75)
GAIN_SHIFT = math.log(1.25)
EXPRESSION_PATTERN_SHIFTS = {
    "BALANCED": (0.0, 0.0),
    "WHOLE_LOSS": (LOSS_SHIFT, LOSS_SHIFT),
    "WHOLE_GAIN": (GAIN_SHIFT, GAIN_SHIFT),
    "P_LOSS": (LOSS_SHIFT, 0.0),
    "P_GAIN": (GAIN_SHIFT, 0.0),
    "Q_LOSS": (0.0, LOSS_SHIFT),
    "Q_GAIN": (0.0, GAIN_SHIFT),
    "RECIPROCAL_P_GAIN_Q_LOSS": (GAIN_SHIFT, LOSS_SHIFT),
    "RECIPROCAL_P_LOSS_Q_GAIN": (LOSS_SHIFT, GAIN_SHIFT),
    "OUTLIER": (0.0, 0.0),
}
EXPRESSION_PATTERNS = tuple(EXPRESSION_PATTERN_SHIFTS)

EXPRESSION_DIRECTION_PATTERNS = {
    ("p", "LOSS"): (
        "WHOLE_LOSS", "P_LOSS", "RECIPROCAL_P_LOSS_Q_GAIN"),
    ("p", "GAIN"): (
        "WHOLE_GAIN", "P_GAIN", "RECIPROCAL_P_GAIN_Q_LOSS"),
    ("q", "LOSS"): (
        "WHOLE_LOSS", "Q_LOSS", "RECIPROCAL_P_GAIN_Q_LOSS"),
    ("q", "GAIN"): (
        "WHOLE_GAIN", "Q_GAIN", "RECIPROCAL_P_LOSS_Q_GAIN"),
}

# For a test of one arm, the sister arm is a nuisance state under both the
# event and null hypotheses.  Comparing an arm event only with BALANCED makes
# a real sister-arm event look like evidence on the tested arm because the
# reciprocal component can fit the sister arm better.  These composite nulls
# explicitly profile/marginalize over every non-outlier sister-arm state.
EXPRESSION_DIRECTION_NULL_PATTERNS = {
    "p": ("BALANCED", "Q_LOSS", "Q_GAIN"),
    "q": ("BALANCED", "P_LOSS", "P_GAIN"),
}

# For each tested arm, represent the nine non-outlier components as a complete
# tested-state x sister-state grid.  Population weights are reduced to a
# *sister-arm marginal* and that same nuisance distribution is used under the
# tested-arm null and both alternatives.  This is essential: separately
# normalizing learned event/null weights lets population correlations make a
# q-only event look like a p event (and can make a true p event look null).
EXPRESSION_TESTED_SISTER_GRID = {
    "p": {
        "BALANCED": {
            "BALANCED": "BALANCED", "LOSS": "P_LOSS", "GAIN": "P_GAIN"},
        "LOSS": {
            "BALANCED": "Q_LOSS", "LOSS": "WHOLE_LOSS",
            "GAIN": "RECIPROCAL_P_GAIN_Q_LOSS"},
        "GAIN": {
            "BALANCED": "Q_GAIN", "LOSS": "RECIPROCAL_P_LOSS_Q_GAIN",
            "GAIN": "WHOLE_GAIN"},
    },
    "q": {
        "BALANCED": {
            "BALANCED": "BALANCED", "LOSS": "Q_LOSS", "GAIN": "Q_GAIN"},
        "LOSS": {
            "BALANCED": "P_LOSS", "LOSS": "WHOLE_LOSS",
            "GAIN": "RECIPROCAL_P_LOSS_Q_GAIN"},
        "GAIN": {
            "BALANCED": "P_GAIN", "LOSS": "RECIPROCAL_P_GAIN_Q_LOSS",
            "GAIN": "WHOLE_GAIN"},
    },
}

BASELINE_LEVELS = (
    "LIBRARY_VALID_GROUP_EXTERNAL_PAIR",
    "LIBRARY_EXTERNAL_PAIR",
    "COHORT_VALID_GROUP_EXTERNAL_PAIR",
    "COHORT_EXTERNAL_PAIR",
)

EVIDENCE_CLASSES = (
    "CONCORDANT_BOTH",
    "EXPRESSION_ONLY",
    "ASE_ONLY",
    "DISCORDANT",
    "EXPRESSION_OUTLIER",
    "BALANCED",
    "INSUFFICIENT_EVIDENCE",
)

CONFIDENCE_TIERS = (
    "HIGH_CONFIDENCE_RNA_CANDIDATE",
    "SINGLE_MODALITY_CANDIDATE",
    "REVIEW_DISCORDANT",
    "NO_CALL",
)

MISSING_IDENTIFIERS = {"", ".", "NA", "N/A", "NONE", "NULL", "UNKNOWN"}

FATAL_HIGH_CONFIDENCE_FLAGS = {
    "PAIRWIDE_EVENT_NOT_IDENTIFIABLE_FROM_RNA",
    "UNADJUSTED_CELL_STATE",
    "HIGH_AMBIENT_EXPRESSION_SENSITIVITY",
    "INSUFFICIENT_GENE_BREADTH",
}


HYBRID_REQUIRED_FIELDS = {
    HYBRID_CELL_MANIFEST_SCHEMA: {
        "library", "barcode", "donor_a", "donor_b", "donor_pair",
        "hybrid_target_eligible", "hybrid_target_reasons",
        "expression_reference_eligible", "expression_reference_reasons",
        "cell_group_source", "cell_group_target_chromosome_excluded",
        "schema_version",
    },
    EXPRESSION_MODEL_SCHEMA: {
        "library", "barcode", "chromosome", "p_arm", "q_arm",
        "donor_a", "donor_b", "donor_pair", "uid", "calibration_group",
        "hybrid_target_eligible", "expression_reference_eligible",
        "cell_group_source", "cell_group_target_chromosome_excluded",
        "p_score", "q_score", "p_fold0_score", "p_fold1_score",
        "q_fold0_score", "q_fold1_score", "mapped_autosomal_library_size",
        "expression_input_state", "normalization_method",
        "gene_fold_method", "schema_version",
    },
    HYBRID_SHARD_SCHEMA: {
        "library", "barcode", "chromosome", "arm",
        "hybrid_target_eligible", "expression_reference_eligible",
        "biological_block", "biological_block_source", "crossfit_fold",
        "ase_present", "expression_present", "schema_version",
    },
    HYBRID_SCORE_SCHEMA: {
        "library", "barcode", "calibration_group", "uid", "donor_a",
        "donor_b", "donor_pair", "chromosome", "arm",
        "hybrid_target_eligible", "ase_eligible", "expression_eligible",
        "ase_status", "expression_status", "ase_calibration_level",
        "expression_calibration_level", "ase_calibration_cells",
        "expression_calibration_cells", "ase_log_bf_DONOR_A_LOSS",
        "ase_log_bf_DONOR_B_LOSS", "ase_log_bf_DONOR_A_GAIN",
        "ase_log_bf_DONOR_B_GAIN", "ase_best_exact_state", "ase_direction",
        "ase_p_DONOR_A_DEPLETED", "ase_p_DONOR_A_ENRICHED",
        "ase_p_floor", "expression_component", "expression_best_copy_state",
        "expression_log_bf_LOSS", "expression_log_bf_GAIN",
        "expression_p_LOSS", "expression_p_GAIN", "expression_fold0_state",
        "expression_fold1_state", "expression_fold_replication_status",
        "expression_top1_fraction", "expression_top5_fraction",
        "expression_top10_fraction", "expression_nonzero_genes",
        "expression_p_floor", "expression_input_state",
        "expression_depth_bin", "expression_breadth_bin",
        "hybrid_conjunction_p_DONOR_A_LOSS",
        "hybrid_conjunction_p_DONOR_B_LOSS",
        "hybrid_conjunction_p_DONOR_A_GAIN",
        "hybrid_conjunction_p_DONOR_B_GAIN", "provisional_evidence_class",
        "provisional_resolved_state", "discordance_reason",
        "confounding_flags", "qc_flags", "schema_version",
    },
    HYBRID_COMPONENT_SCHEMA: {
        "chromosome", "fold", "calibration_level", "component",
        "weight", "status", "schema_version",
    },
    HYBRID_CALIBRATION_SCHEMA: {
        "chromosome", "calibration_key", "branch", "calibration_level",
        "calibration_cells", "status", "schema_version",
    },
    HYBRID_CALL_SCHEMA: {
        "library", "barcode", "chromosome", "arm", "evidence_class",
        "copy_state", "donor_origin", "resolved_state",
        "confidence_tier", "call_state", "call_status", "schema_version",
    },
    HYBRID_UID_SCHEMA: {
        "uid", "donor_pair", "chromosome", "distinct_uid_blocks",
        "summary_status", "schema_version",
    },
    HYBRID_PAIR_SCHEMA: {
        "donor_pair", "chromosome", "arm", "distinct_uid_blocks",
        "summary_status", "schema_version",
    },
}


def validate_hybrid_table(path: str, schema: str,
                          allow_header_only: bool = True) -> tuple[list[str], int]:
    required = HYBRID_REQUIRED_FIELDS.get(schema)
    if required is None:
        raise ValueError(f"no hybrid validator is registered for {schema}")
    return validate_tsv_schema(
        path, schema, required, allow_header_only=allow_header_only)


def present_identifier(value: object) -> str:
    value = clean(value)
    return "" if value.upper() in MISSING_IDENTIFIERS else value


def canonical_chromosome(value: object) -> str:
    text = clean(value)
    if text.lower().startswith("chr"):
        text = text[3:]
    if text.isdigit():
        return str(int(text))
    return text.upper()


def is_autosomal_chromosome(value: object) -> bool:
    text = canonical_chromosome(value)
    return text.isdigit() and 1 <= int(text) <= 22


def chromosome_from_arm(value: object) -> str:
    text = clean(value)
    if text.lower().endswith(("p", "q")):
        text = text[:-1]
    return canonical_chromosome(text)


def arm_side(value: object) -> str:
    text = clean(value).lower()
    if text.endswith("p"):
        return "p"
    if text.endswith("q"):
        return "q"
    return ""


def stable_fold(value: object, folds: int, namespace: bytes = b"cellbouncer") -> int:
    if folds < 1:
        raise ValueError("fold count must be positive")
    digest = hashlib.blake2b(
        str(value).encode("utf-8"), digest_size=16, person=namespace[:16]).digest()
    return int.from_bytes(digest[:8], "big") % folds


def stable_gene_fold(gene_identifier: object, folds: int = 2) -> int:
    identifier = clean(gene_identifier)
    if not identifier:
        raise ValueError("canonical gene identifier is empty")
    return stable_fold(identifier, folds, b"tetra-gene-fold")


def biological_block(row: Mapping[str, object]) -> tuple[str, str, list[str]]:
    uid = present_identifier(row.get("uid"))
    if uid:
        return f"UID:{uid}", "UID", []
    pair = present_identifier(row.get("donor_pair"))
    if pair:
        return f"PAIR:{pair}", "DONOR_PAIR", []
    library = clean(row.get("library"))
    barcode = clean(row.get("barcode"))
    if not library or not barcode:
        raise ValueError("biological block fallback requires library and barcode")
    return (
        f"CELL:{library}:{barcode}",
        "LIBRARY_BARCODE_FALLBACK",
        ["BIOLOGICAL_BLOCK_LIBRARY_BARCODE_FALLBACK"],
    )


def biological_crossfit_fold(row: Mapping[str, object], folds: int) -> int:
    block, _source, _flags = biological_block(row)
    return stable_fold(block, folds, b"tetra-cell-fold")


def bounded(value: float, lower: float, upper: float) -> float:
    return min(upper, max(lower, value))


def logistic(value: float) -> float:
    if value >= 0.0:
        return 1.0 / (1.0 + math.exp(-value))
    exponential = math.exp(value)
    return exponential / (1.0 + exponential)


def logit(value: float) -> float:
    clipped = bounded(value, 1e-12, 1.0 - 1e-12)
    return math.log(clipped / (1.0 - clipped))


def logsumexp(values: Sequence[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return -math.inf
    maximum = max(finite)
    return maximum + math.log(sum(math.exp(value - maximum) for value in finite))


def format_number(value: object) -> str:
    number = finite_float(value)
    return "NA" if not math.isfinite(number) else f"{number:.17g}"


def split_flags(value: object) -> set[str]:
    return {
        clean(token) for token in str(value or "").split(";")
        if clean(token) and clean(token).upper() not in {"PASS", "NONE"}
    }


def join_flags(values: Iterable[object]) -> str:
    flags: set[str] = set()
    for value in values:
        if isinstance(value, (set, list, tuple)):
            for child in value:
                flags.update(split_flags(child))
        else:
            flags.update(split_flags(value))
    return ";".join(sorted(flags, key=natural_key)) or "PASS"


def weighted_median(values: Sequence[float], weights: Sequence[float]) -> float:
    pairs = sorted(
        (float(value), max(0.0, float(weight)))
        for value, weight in zip(values, weights)
        if math.isfinite(value) and math.isfinite(weight) and weight > 0.0)
    if not pairs:
        return math.nan
    total = sum(weight for _value, weight in pairs)
    cumulative = 0.0
    for value, weight in pairs:
        cumulative += weight
        if cumulative >= total / 2.0:
            return value
    return pairs[-1][0]


@dataclass(frozen=True)
class RobustBivariateFit:
    center: tuple[float, float]
    covariance: tuple[tuple[float, float], tuple[float, float]]
    dispersion: float
    observations: int
    complete_observations: int
    retained_observations: int
    status: str


def _clamp_covariance(covariance: np.ndarray, minimum_variance: float,
                      maximum_variance: float,
                      maximum_abs_correlation: float) -> np.ndarray:
    covariance = np.asarray(covariance, dtype=float)
    covariance = (covariance + covariance.T) / 2.0
    variances = np.clip(np.diag(covariance), minimum_variance, maximum_variance)
    correlation = 0.0
    denominator = math.sqrt(float(variances[0] * variances[1]))
    if denominator > 0.0 and math.isfinite(float(covariance[0, 1])):
        correlation = bounded(
            float(covariance[0, 1]) / denominator,
            -maximum_abs_correlation, maximum_abs_correlation)
    result = np.array([
        [variances[0], correlation * denominator],
        [correlation * denominator, variances[1]],
    ], dtype=float)
    eigenvalues, eigenvectors = np.linalg.eigh(result)
    eigenvalues = np.clip(eigenvalues, minimum_variance, maximum_variance)
    return eigenvectors @ np.diag(eigenvalues) @ eigenvectors.T


def robust_bivariate_fit(
        points: Sequence[Sequence[float]], minimum_sigma: float = 0.10,
        maximum_sigma: float = 2.0, maximum_abs_correlation: float = 0.95,
        degrees_freedom: float = 4.0, iterations: int = 8,
        center_hint: Sequence[float] | None = None) -> RobustBivariateFit:
    """Fit a deterministic robust p/q center and covariance.

    Missing p or q values contribute to their marginal center/scale but not to
    covariance.  The event components remain fixed relative to this balanced
    fit; this routine never learns event-specific centers or widths.
    """
    if minimum_sigma <= 0.0 or maximum_sigma < minimum_sigma:
        raise ValueError("invalid robust covariance bounds")
    array = np.asarray(points, dtype=float)
    if array.size == 0:
        array = np.empty((0, 2), dtype=float)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError("bivariate points must have exactly two columns")
    finite = np.isfinite(array)
    observations = int(np.sum(np.any(finite, axis=1)))
    complete = array[np.all(finite, axis=1)]
    minimum_variance = minimum_sigma ** 2
    maximum_variance = maximum_sigma ** 2

    center_values: list[float] = []
    for dimension in range(2):
        marginal = array[finite[:, dimension], dimension]
        if center_hint is not None and math.isfinite(float(center_hint[dimension])):
            center_values.append(float(center_hint[dimension]))
        elif marginal.size:
            center_values.append(float(np.median(marginal)))
        else:
            center_values.append(0.0)
    center = np.asarray(center_values, dtype=float)

    marginal_sigmas = []
    for dimension in range(2):
        marginal = array[finite[:, dimension], dimension]
        if marginal.size:
            mad = 1.4826 * float(np.median(np.abs(marginal - center[dimension])))
            marginal_sigmas.append(bounded(
                mad if math.isfinite(mad) and mad > 0.0 else minimum_sigma,
                minimum_sigma, maximum_sigma))
        else:
            marginal_sigmas.append(minimum_sigma)
    covariance = np.diag(np.square(marginal_sigmas))
    retained = len(complete)

    if len(complete) >= 3:
        covariance = np.cov(complete, rowvar=False, ddof=1)
        covariance = _clamp_covariance(
            covariance, minimum_variance, maximum_variance,
            maximum_abs_correlation)
        weights = np.ones(len(complete), dtype=float)
        for _ in range(max(1, iterations)):
            inverse = np.linalg.inv(covariance)
            residual = complete - center
            distance = np.einsum("ij,jk,ik->i", residual, inverse, residual)
            weights = (degrees_freedom + 2.0) / (degrees_freedom + distance)
            total = float(np.sum(weights))
            if total <= 0.0:
                break
            center = np.sum(complete * weights[:, None], axis=0) / total
            residual = complete - center
            covariance = (
                residual.T @ (residual * weights[:, None]) / max(total, 1.0))
            covariance = _clamp_covariance(
                covariance, minimum_variance, maximum_variance,
                maximum_abs_correlation)
        retained = int(np.sum(weights >= 0.25))

    correlation_denominator = math.sqrt(
        float(covariance[0, 0] * covariance[1, 1]))
    correlation = (
        float(covariance[0, 1]) / correlation_denominator
        if correlation_denominator > 0.0 else 0.0)
    status = "PASS" if observations >= 3 else "WEAK_COMPONENT_LOW_SUPPORT"
    return RobustBivariateFit(
        (float(center[0]), float(center[1])),
        ((float(covariance[0, 0]), float(covariance[0, 1])),
         (float(covariance[1, 0]), float(covariance[1, 1]))),
        correlation, observations, len(complete), retained, status)


def multivariate_student_log_density(
        point: Sequence[float], mean: Sequence[float],
        covariance: Sequence[Sequence[float]], degrees_freedom: float = 4.0,
        covariance_scale: float = 1.0) -> float:
    observed = [index for index, value in enumerate(point)
                if math.isfinite(float(value))]
    if not observed:
        return math.nan
    scale = max(float(covariance_scale), 1e-12)
    if len(observed) == 1:
        index = observed[0]
        variance = float(covariance[index][index]) * scale
        if variance <= 0.0 or not math.isfinite(variance):
            return -math.inf
        residual = float(point[index]) - float(mean[index])
        log_determinant = math.log(variance)
        distance = residual * residual / variance
    else:
        residual_p = float(point[0]) - float(mean[0])
        residual_q = float(point[1]) - float(mean[1])
        variance_p = float(covariance[0][0]) * scale
        variance_q = float(covariance[1][1]) * scale
        cross = (float(covariance[0][1]) + float(covariance[1][0])) \
            * 0.5 * scale
        determinant = variance_p * variance_q - cross * cross
        if determinant <= 0.0 or not math.isfinite(determinant):
            return -math.inf
        log_determinant = math.log(determinant)
        distance = (
            variance_q * residual_p * residual_p
            - 2.0 * cross * residual_p * residual_q
            + variance_p * residual_q * residual_q) / determinant
    dimension = len(observed)
    df = max(float(degrees_freedom), 1e-6)
    return (
        math.lgamma((df + dimension) / 2.0) - math.lgamma(df / 2.0)
        - 0.5 * (dimension * math.log(df * math.pi) + log_determinant)
        - 0.5 * (df + dimension) * math.log1p(distance / df)
    )


@dataclass(frozen=True)
class FixedMixtureFit:
    weights: Mapping[str, float]
    iterations: int
    log_likelihood: float
    status: str


@dataclass(frozen=True)
class AseBiasFit:
    observations: int
    total_weight: float
    delta: float
    rho: float
    observed_fraction: float
    retained_observations: int


def fit_ase_delta(
        observations: Sequence[tuple[float, float, float]],
        robust_weights: Sequence[float] | None = None,
        maximum_weight: float = 50.0) -> float:
    """Native-compatible mapping-bias score-equation fit."""
    if not observations:
        return 0.0

    values = np.asarray(observations, dtype=float)
    fraction = values[:, 0]
    weight = np.minimum(values[:, 1], maximum_weight)
    if robust_weights is not None:
        weight = weight * np.asarray(robust_weights, dtype=float)
    base = np.clip(values[:, 2], 1e-12, 1.0 - 1e-12)
    base_logit = np.log(base / (1.0 - base))

    def score(delta: float) -> float:
        eta = base_logit + delta
        expected = np.empty_like(eta)
        positive = eta >= 0.0
        expected[positive] = 1.0 / (1.0 + np.exp(-eta[positive]))
        exponential = np.exp(eta[~positive])
        expected[~positive] = exponential / (1.0 + exponential)
        return float(np.sum(weight * (fraction - expected)))

    lower, upper = -6.0, 6.0
    if score(lower) <= 0.0:
        return lower
    if score(upper) >= 0.0:
        return upper
    for _ in range(80):
        midpoint = (lower + upper) / 2.0
        if score(midpoint) > 0.0:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0


def fit_ase_bias_component(
        observations: Sequence[tuple[float, float, float]],
        rho_observations: Sequence[tuple[float, float, float]] | None = None,
        n_cells: int | None = None,
        maximum_weight: float = 50.0, huber_z: float = 2.5,
        minimum_robust_weight: float = 0.25, default_rho: float = 0.02,
        minimum_rho: float = 0.001, maximum_rho: float = 0.25,
        ) -> AseBiasFit:
    """Port the legacy robust ASE bias/rho calibration for hybrid use."""
    if not observations:
        return AseBiasFit(0, 0.0, 0.0, default_rho, math.nan, 0)
    rho_scale = (list(rho_observations) if rho_observations is not None
                 else list(observations))
    if n_cells is None:
        n_cells = len(observations)

    def expected_fractions(subset: np.ndarray, delta: float) -> np.ndarray:
        base = np.clip(subset[:, 2], 1e-12, 1.0 - 1e-12)
        eta = np.log(base / (1.0 - base)) + delta
        result = np.empty_like(eta)
        positive = eta >= 0.0
        result[positive] = 1.0 / (1.0 + np.exp(-eta[positive]))
        exponential = np.exp(eta[~positive])
        result[~positive] = exponential / (1.0 + exponential)
        return result

    def estimate_rho(delta: float, weights: Sequence[float],
                     subset: Sequence[tuple[float, float, float]]) -> float:
        if not subset:
            return default_rho
        values = np.asarray(subset, dtype=float)
        selected = values[:, 1] > 1.0
        if not np.any(selected):
            return default_rho
        values = values[selected]
        weight = np.asarray(weights, dtype=float)[selected]
        expected = expected_fractions(values, delta)
        binomial = expected * (1.0 - expected) / values[:, 1]
        numerator = float(np.sum(weight * (
            (values[:, 0] - expected) ** 2 - binomial)))
        denominator = float(np.sum(weight * expected * (1.0 - expected)
                                   * (1.0 - 1.0 / values[:, 1])))
        if denominator <= 0.0:
            return default_rho
        return bounded(
            max(0.0, numerator) / denominator, minimum_rho, maximum_rho)

    baseline_robust = [1.0] * len(observations)
    rho_robust = [1.0] * len(rho_scale)
    delta, rho = 0.0, default_rho
    for _ in range(5):
        delta = fit_ase_delta(observations, baseline_robust, maximum_weight)
        rho = estimate_rho(delta, rho_robust, rho_scale)

        def robust_weights(
                subset: Sequence[tuple[float, float, float]]) -> list[float]:
            if not subset:
                return []
            values = np.asarray(subset, dtype=float)
            expected = expected_fractions(values, delta)
            variance = np.maximum(expected * (1.0 - expected) * (
                rho + (1.0 - rho) / np.maximum(values[:, 1], 1e-8)), 1e-6)
            z_score = np.abs(values[:, 0] - expected) / np.sqrt(variance)
            return np.minimum(1.0, huber_z / np.maximum(
                z_score, 1e-12)).tolist()

        baseline_robust = robust_weights(observations)
        rho_robust = robust_weights(rho_scale)

    retained_mask = [weight >= minimum_robust_weight
                     for weight in baseline_robust]
    if sum(retained_mask) >= 3:
        retained = [observation for observation, keep in zip(
            observations, retained_mask) if keep]
        delta = fit_ase_delta(retained, None, maximum_weight)
        working = [(1.0, observation) for observation in retained]
    else:
        working = list(zip(baseline_robust, observations))
    retained_rho_mask = [weight >= minimum_robust_weight
                         for weight in rho_robust]
    if sum(retained_rho_mask) >= 3:
        retained_rho = [observation for observation, keep in zip(
            rho_scale, retained_rho_mask) if keep]
        rho = estimate_rho(delta, [1.0] * len(retained_rho), retained_rho)
    else:
        rho = estimate_rho(delta, rho_robust, rho_scale)
    if len(rho_scale) < 3:
        rho = default_rho
    rho = bounded(rho, minimum_rho, maximum_rho)
    total_weight = sum(
        weight * min(count, maximum_weight)
        for weight, (_fraction, count, _base) in working)
    observed_weight = sum(count for _fraction, count, _base in observations)
    observed = (
        sum(fraction * count for fraction, count, _base in observations)
        / observed_weight if observed_weight else math.nan)
    return AseBiasFit(
        n_cells, total_weight, delta, rho, observed,
        sum(retained_mask))


def component_log_likelihoods(
        point: Sequence[float], center: Sequence[float],
        covariance: Sequence[Sequence[float]], degrees_freedom: float = 4.0,
        outlier_scale: float = 16.0) -> dict[str, float]:
    result: dict[str, float] = {}
    for component, shift in EXPRESSION_PATTERN_SHIFTS.items():
        mean = (center[0] + shift[0], center[1] + shift[1])
        result[component] = multivariate_student_log_density(
            point, mean, covariance, degrees_freedom,
            outlier_scale if component == "OUTLIER" else 1.0)
    return result


def _component_density_matrix(
        points: np.ndarray, center: Sequence[float],
        covariance: Sequence[Sequence[float]], degrees_freedom: float,
        outlier_scale: float) -> np.ndarray:
    """Evaluate fixed Student components in batches by observed arm pattern."""
    finite = np.isfinite(points)
    both = finite[:, 0] & finite[:, 1]
    p_only = finite[:, 0] & ~finite[:, 1]
    q_only = ~finite[:, 0] & finite[:, 1]
    densities = np.full((len(points), len(EXPRESSION_PATTERNS)), -math.inf)
    df = max(float(degrees_freedom), 1e-6)
    for column, component in enumerate(EXPRESSION_PATTERNS):
        shift_p, shift_q = EXPRESSION_PATTERN_SHIFTS[component]
        mean_p = center[0] + shift_p
        mean_q = center[1] + shift_q
        scale = max(float(outlier_scale if component == "OUTLIER" else 1.0),
                    1e-12)
        variance_p = float(covariance[0][0]) * scale
        variance_q = float(covariance[1][1]) * scale
        if p_only.any() and variance_p > 0.0 and math.isfinite(variance_p):
            residual = points[p_only, 0] - mean_p
            constant = (math.lgamma((df + 1.0) / 2.0)
                        - math.lgamma(df / 2.0)
                        - 0.5 * (math.log(df * math.pi)
                                 + math.log(variance_p)))
            densities[p_only, column] = (
                constant - 0.5 * (df + 1.0)
                * np.log1p((residual * residual / variance_p) / df))
        if q_only.any() and variance_q > 0.0 and math.isfinite(variance_q):
            residual = points[q_only, 1] - mean_q
            constant = (math.lgamma((df + 1.0) / 2.0)
                        - math.lgamma(df / 2.0)
                        - 0.5 * (math.log(df * math.pi)
                                 + math.log(variance_q)))
            densities[q_only, column] = (
                constant - 0.5 * (df + 1.0)
                * np.log1p((residual * residual / variance_q) / df))
        cross = (float(covariance[0][1]) + float(covariance[1][0])) \
            * 0.5 * scale
        determinant = variance_p * variance_q - cross * cross
        if both.any() and determinant > 0.0 and math.isfinite(determinant):
            residual_p = points[both, 0] - mean_p
            residual_q = points[both, 1] - mean_q
            distance = (
                variance_q * residual_p * residual_p
                - 2.0 * cross * residual_p * residual_q
                + variance_p * residual_q * residual_q) / determinant
            constant = (math.lgamma((df + 2.0) / 2.0)
                        - math.lgamma(df / 2.0)
                        - 0.5 * (2.0 * math.log(df * math.pi)
                                 + math.log(determinant)))
            densities[both, column] = (
                constant - 0.5 * (df + 2.0) * np.log1p(distance / df))
    return densities


def fit_fixed_center_mixture(
        points: Sequence[Sequence[float]], center: Sequence[float],
        covariance: Sequence[Sequence[float]], degrees_freedom: float = 4.0,
        outlier_scale: float = 16.0, minimum_weight: float = 1e-6,
        maximum_iterations: int = 100, tolerance: float = 1e-10,
        ) -> FixedMixtureFit:
    point_array = np.asarray(points, dtype=float)
    usable = (point_array[np.isfinite(point_array).any(axis=1)]
              if len(points) else np.empty((0, 2), dtype=float))
    count = len(EXPRESSION_PATTERNS)
    if len(usable) == 0:
        weights = {component: 1.0 / count for component in EXPRESSION_PATTERNS}
        return FixedMixtureFit(weights, 0, math.nan, "NO_DATA")
    initial = {component: 0.2 / (count - 1) for component in EXPRESSION_PATTERNS}
    initial["BALANCED"] = 0.8
    weights = initial
    # Component centers/covariance do not move during EM. Evaluate their
    # densities and stabilized exponentials once per discovery point.
    components = EXPRESSION_PATTERNS
    densities = _component_density_matrix(
        usable, center, covariance, degrees_freedom, outlier_scale)
    row_maximum = np.max(densities, axis=1)
    fixed_densities = np.exp(densities - row_maximum[:, None])
    previous = -math.inf
    final_log_likelihood = -math.inf
    completed = 0
    for iteration in range(1, maximum_iterations + 1):
        weight_array = np.asarray([
            max(weights[component], minimum_weight)
            for component in components])
        denominators = np.einsum(
            "ij,j->i", fixed_densities, weight_array, optimize=False)
        final_log_likelihood = float(np.sum(
            row_maximum + np.log(denominators)))
        # Sum responsibilities without allocating two full point-by-component
        # arrays on every iteration. The component densities remain fixed.
        totals = weight_array * np.einsum(
            "ij,i->j", fixed_densities, np.reciprocal(denominators),
            optimize=False)
        raw = {
            component: max(minimum_weight, float(totals[index]) / len(usable))
            for index, component in enumerate(components)
        }
        denominator = sum(raw.values())
        weights = {component: raw[component] / denominator
                   for component in EXPRESSION_PATTERNS}
        completed = iteration
        if math.isfinite(previous) and abs(final_log_likelihood - previous) <= (
                tolerance * max(1.0, abs(previous))):
            break
        previous = final_log_likelihood
    status = "PASS" if len(usable) >= 3 else "WEAK_COMPONENT_LOW_SUPPORT"
    return FixedMixtureFit(weights, completed, final_log_likelihood, status)


def best_expression_component(
        point: Sequence[float], center: Sequence[float],
        covariance: Sequence[Sequence[float]], weights: Mapping[str, float],
        degrees_freedom: float = 4.0, outlier_scale: float = 16.0) -> str:
    likelihoods = component_log_likelihoods(
        point, center, covariance, degrees_freedom, outlier_scale)
    return max(
        EXPRESSION_PATTERNS,
        key=lambda component: (
            math.log(max(float(weights.get(component, 0.0)), 1e-300))
            + likelihoods[component],
            -EXPRESSION_PATTERNS.index(component),
        ))


def expression_direction_log_bfs(
        point: Sequence[float], side: str, center: Sequence[float],
        covariance: Sequence[Sequence[float]], degrees_freedom: float = 4.0,
        weights: Mapping[str, float] | None = None,
        outlier_scale: float = 16.0,
        ) -> dict[str, float]:
    """Return sister-arm-marginalized, population-weighted log BFs.

    Learned weights contribute only the common sister-arm nuisance marginal.
    The identical marginal is used for BALANCED, LOSS, and GAIN on the tested
    arm, so population correlations cannot supply evidence for (or against)
    the tested-arm state.
    """
    if side not in {"p", "q"}:
        return {"LOSS": math.nan, "GAIN": math.nan}
    likelihoods = component_log_likelihoods(
        point, center, covariance, degrees_freedom, outlier_scale)
    if weights is None:
        weights = {
            component: 1.0 / len(EXPRESSION_PATTERNS)
            for component in EXPRESSION_PATTERNS
        }

    grid = EXPRESSION_TESTED_SISTER_GRID[side]
    sister_weights = {
        sister_state: sum(max(float(weights.get(
            grid[sister_state][tested_state], 0.0)), 0.0)
            for tested_state in ("BALANCED", "LOSS", "GAIN"))
        for sister_state in ("BALANCED", "LOSS", "GAIN")
    }
    normalizer = sum(sister_weights.values())
    if not math.isfinite(normalizer) or normalizer <= 0.0:
        sister_weights = {state: 1.0 / 3.0
                          for state in ("BALANCED", "LOSS", "GAIN")}
    else:
        sister_weights = {
            state: max(value / normalizer, 1e-300)
            for state, value in sister_weights.items()
        }

    def tested_state_log_likelihood(tested_state: str) -> float:
        return logsumexp([
            math.log(sister_weights[sister_state])
            + likelihoods[grid[sister_state][tested_state]]
            for sister_state in ("BALANCED", "LOSS", "GAIN")
        ])

    null = tested_state_log_likelihood("BALANCED")
    return {
        direction: tested_state_log_likelihood(direction) - null
        for direction in COPY_STATES
    }


def direction_from_log_bfs(loss: float, gain: float,
                           minimum_log_bf: float = 0.0) -> str:
    if not math.isfinite(loss) and not math.isfinite(gain):
        return "NO_DATA"
    loss_value = loss if math.isfinite(loss) else -math.inf
    gain_value = gain if math.isfinite(gain) else -math.inf
    if max(loss_value, gain_value) <= minimum_log_bf:
        return "BALANCED"
    return "LOSS" if loss_value >= gain_value else "GAIN"


def empirical_upper_tail(statistic: float,
                         sorted_null: Sequence[float]) -> float:
    if not math.isfinite(statistic) or not sorted_null:
        return math.nan
    index = bisect.bisect_left(sorted_null, statistic - 1e-12)
    return (len(sorted_null) - index + 1.0) / (len(sorted_null) + 1.0)


def empirical_p_floor(null_count: int) -> float:
    return 1.0 / (null_count + 1.0) if null_count >= 0 else math.nan


def bh_adjust(values: Sequence[float]) -> list[float]:
    result = [math.nan] * len(values)
    finite_indices = [index for index, value in enumerate(values)
                      if math.isfinite(value)]
    ordered = sorted(finite_indices, key=lambda index: (values[index], index))
    running = 1.0
    total = len(ordered)
    for reverse_rank, index in enumerate(reversed(ordered), start=1):
        rank = total - reverse_rank + 1
        running = min(running, min(1.0, values[index] * total / rank))
        result[index] = running
    return result


def by_adjust(values: Sequence[float]) -> list[float]:
    finite_count = sum(math.isfinite(value) for value in values)
    if finite_count == 0:
        return [math.nan] * len(values)
    harmonic = sum(1.0 / rank for rank in range(1, finite_count + 1))
    adjusted_input = [
        min(1.0, value * harmonic) if math.isfinite(value) else math.nan
        for value in values
    ]
    return bh_adjust(adjusted_input)


def conjunction_pvalues(expression_p: Mapping[str, float],
                        ase_p: Mapping[str, float]) -> dict[str, float]:
    result: dict[str, float] = {}
    for state in EXACT_STATES:
        copy_state, ase_direction, _donor = STATE_AXES[state]
        left = finite_float(expression_p.get(copy_state))
        right = finite_float(ase_p.get(ase_direction))
        result[state] = (
            max(left, right)
            if math.isfinite(left) and math.isfinite(right) else math.nan)
    return result


def student_t_log_density(value: float, mean: float, sigma: float,
                          degrees_freedom: float = 4.0) -> float:
    sigma = max(float(sigma), 1e-6)
    standardized = (float(value) - float(mean)) / sigma
    return -math.log(sigma) - 0.5 * (degrees_freedom + 1.0) * math.log1p(
        standardized * standardized / degrees_freedom)


def quasi_ase_log_likelihood(effective_a: float, effective_weight: float,
                             probability: float, rho: float) -> float:
    if effective_weight <= 0.0:
        return 0.0
    fraction = bounded(effective_a / effective_weight, 0.0, 1.0)
    probability = bounded(probability, 1e-8, 1.0 - 1e-8)
    overdispersion = bounded(rho, 0.0, 1.0)
    variance = probability * (1.0 - probability) * (
        overdispersion + (1.0 - overdispersion) / effective_weight)
    return student_t_log_density(
        fraction, probability, math.sqrt(max(variance, 1e-8)))


def ambient_adjusted_fraction(
        state: str, orientation: str, baseline_logit: float,
        contamination: float, ambient_a: float, mapping_delta: float) -> float:
    if state not in STATE_LOG_ODDS_SHIFT:
        raise ValueError(f"unknown ASE state: {state}")
    cellular = logistic(baseline_logit + STATE_LOG_ODDS_SHIFT[state])
    mixed = ((1.0 - bounded(contamination, 0.0, 0.999999)) * cellular
             + bounded(contamination, 0.0, 0.999999)
             * bounded(ambient_a, 0.0, 1.0))
    return bounded(logistic(logit(mixed) + mapping_delta), 1e-8, 1.0 - 1e-8)


def contamination_quadrature(contamination: float, standard_error: float,
                              strata: int = 15) -> list[tuple[float, float]]:
    center = bounded(contamination, 0.0, 0.999999)
    if not math.isfinite(standard_error) or standard_error <= 0.0:
        return [(center, 1.0)]
    normal = statistics.NormalDist()
    lower_probability = normal.cdf((0.0 - center) / standard_error)
    upper_probability = normal.cdf((1.0 - center) / standard_error)
    interval = upper_probability - lower_probability
    if interval <= 1e-15:
        return [(center, 1.0)]
    result = []
    for index in range(strata):
        probability = lower_probability + (index + 0.5) * interval / strata
        probability = bounded(probability, 1e-15, 1.0 - 1e-15)
        value = center + standard_error * normal.inv_cdf(probability)
        result.append((bounded(value, 0.0, 0.999999), 1.0 / strata))
    return result


def ase_state_log_bfs(
        evidence: Mapping[str, object], baseline_logit: float,
        mapping_offsets: Mapping[str, float], rhos: Mapping[str, float],
        site_fallback_weight: float = 0.50) -> dict[str, float]:
    contamination = bounded(finite_float(evidence.get("ambient_c"), 0.0),
                            0.0, 0.999999)
    standard_error = finite_float(evidence.get("ambient_c_se"))
    nodes = contamination_quadrature(contamination, standard_error)
    integrands: dict[str, list[float]] = {
        state: [] for state in STATE_LOG_ODDS_SHIFT}
    evidence_weight = (
        site_fallback_weight
        if clean(evidence.get("evidence_status")) == "PASS_SITE_FALLBACK"
        else 1.0)
    orientation_evidence = tuple((
        orientation,
        finite_float(evidence.get(f"effective_a_{orientation}"), 0.0),
        finite_float(evidence.get(f"effective_weight_{orientation}"), 0.0),
        finite_float(evidence.get(f"ambient_a_{orientation}"), 0.5),
        float(mapping_offsets.get(orientation, 0.0)),
        float(rhos.get(orientation, 0.02)),
    ) for orientation in ("ref", "alt", "mixed"))
    for node_contamination, node_weight in nodes:
        likelihoods: dict[str, float] = {}
        for state in STATE_LOG_ODDS_SHIFT:
            value = 0.0
            for (orientation, effective_a, effective_weight, ambient_a,
                 mapping_delta, rho) in orientation_evidence:
                probability = ambient_adjusted_fraction(
                    state, orientation, baseline_logit, node_contamination,
                    ambient_a, mapping_delta)
                value += quasi_ase_log_likelihood(
                    effective_a, effective_weight, probability, rho)
            likelihoods[state] = evidence_weight * value
        for state in STATE_LOG_ODDS_SHIFT:
            integrands[state].append(math.log(node_weight) + likelihoods[state])
    integrated = {state: logsumexp(values) for state, values in integrands.items()}
    balanced = integrated["BALANCED"]
    return {state: integrated[state] - balanced for state in EXACT_STATES}


def ase_direction_statistics(log_bfs: Mapping[str, float]) -> dict[str, float]:
    return {
        direction: max(0.0, max(finite_float(log_bfs.get(state), -math.inf)
                                for state in states))
        for direction, states in ASE_DIRECTION_STATES.items()
    }


def inverse_ambient_logit(observed_fraction: float, contamination: float,
                          ambient_fraction: float,
                          mapping_delta: float = 0.0) -> float:
    unbiased = logistic(logit(observed_fraction) - mapping_delta)
    contamination = bounded(contamination, 0.0, 0.999999)
    cellular = ((unbiased - contamination * bounded(ambient_fraction, 0.0, 1.0))
                / max(1.0 - contamination, 1e-8))
    return logit(bounded(cellular, 1e-6, 1.0 - 1e-6))


def copy_and_donor_for_state(state: str) -> tuple[str, str]:
    if state not in STATE_AXES:
        return "UNRESOLVED", "UNRESOLVED"
    copy_state, _direction, donor = STATE_AXES[state]
    return copy_state, donor


def classify_pq_relationship(p_state: str, q_state: str) -> str:
    p_event = p_state in EXACT_STATES
    q_event = q_state in EXACT_STATES
    if not p_event and not q_event:
        return "NO_RESOLVED_EVENT"
    if p_event and not q_event:
        return "P_ARM_ONLY"
    if q_event and not p_event:
        return "Q_ARM_ONLY"
    p_copy, p_donor = copy_and_donor_for_state(p_state)
    q_copy, q_donor = copy_and_donor_for_state(q_state)
    if p_copy == q_copy and p_donor == q_donor:
        return "WHOLE_CHROMOSOME_CONCORDANT"
    if p_copy != q_copy:
        return "RECIPROCAL_OR_ISOCHROMOSOME_LIKE"
    return "MULTI_ARM_MIXED_DONOR_ORIGIN"


def resolve_evidence(
        expression_state: str, expression_significant: bool,
        ase_direction: str, ase_significant: bool,
        conjunction_state: str, conjunction_significant: bool,
        expression_fold_status: str, branch_data_available: bool,
        confounding_flags: Iterable[str] = (),
        ) -> dict[str, str]:
    """Resolve final semantics without letting one branch supply both axes."""
    flags = {clean(value) for value in confounding_flags if clean(value)}
    if expression_fold_status in {"DISAGREE", "OUTLIER_COMPONENT"}:
        evidence_class = "EXPRESSION_OUTLIER"
    elif expression_significant and ase_significant:
        evidence_class = (
            "CONCORDANT_BOTH" if conjunction_significant
            and conjunction_state in EXACT_STATES else "DISCORDANT")
    elif expression_significant:
        evidence_class = "EXPRESSION_ONLY"
        flags.add("EXPRESSION_ONLY_DONOR_UNRESOLVED")
    elif ase_significant:
        evidence_class = "ASE_ONLY"
        flags.add("ASE_ONLY_COPY_UNRESOLVED")
    elif branch_data_available:
        evidence_class = "BALANCED"
    else:
        evidence_class = "INSUFFICIENT_EVIDENCE"

    copy_state = "UNRESOLVED"
    donor_origin = "UNRESOLVED"
    resolved_state = "NO_CALL"
    if evidence_class == "CONCORDANT_BOTH":
        resolved_state = conjunction_state
        copy_state, donor_origin = copy_and_donor_for_state(conjunction_state)
    elif evidence_class == "EXPRESSION_ONLY":
        copy_state = expression_state if expression_state in COPY_STATES else "UNRESOLVED"
    elif evidence_class == "ASE_ONLY":
        donor_origin = (
            "DONOR_A_DEPLETED" if ase_direction == "DONOR_A_DEPLETED"
            else "DONOR_A_ENRICHED" if ase_direction == "DONOR_A_ENRICHED"
            else "UNRESOLVED")

    fatal = bool(flags & FATAL_HIGH_CONFIDENCE_FLAGS)
    if evidence_class == "CONCORDANT_BOTH" and not fatal:
        confidence = "HIGH_CONFIDENCE_RNA_CANDIDATE"
        call_state = resolved_state
        call_status = "PASS_EVENT"
        interpretation = "CONCORDANT_RNA_INFERRED_DOSAGE_AND_DONOR_DIRECTION"
    elif evidence_class in {"EXPRESSION_ONLY", "ASE_ONLY"}:
        confidence = "SINGLE_MODALITY_CANDIDATE"
        call_state = evidence_class
        call_status = "REVIEW_SINGLE_MODALITY"
        interpretation = (
            "UNORIENTED_RNA_DOSAGE_CANDIDATE"
            if evidence_class == "EXPRESSION_ONLY"
            else "COPY_UNRESOLVED_DONOR_DIRECTION_CANDIDATE")
    elif evidence_class in {"DISCORDANT", "EXPRESSION_OUTLIER"} or (
            evidence_class == "CONCORDANT_BOTH" and fatal):
        confidence = "REVIEW_DISCORDANT"
        call_state = evidence_class
        call_status = "REVIEW"
        interpretation = "RNA_BRANCH_DISCORDANCE_OR_CONFOUNDING_REVIEW"
    elif evidence_class == "BALANCED":
        confidence = "NO_CALL"
        call_state = "BALANCED"
        call_status = "BALANCED"
        interpretation = "NO_SIGNIFICANT_RNA_ARM_EVENT"
    else:
        confidence = "NO_CALL"
        call_state = "NO_CALL"
        call_status = "INSUFFICIENT_POWER"
        interpretation = "INSUFFICIENT_RNA_EVIDENCE"

    return {
        "evidence_class": evidence_class,
        "copy_state": copy_state,
        "donor_origin": donor_origin,
        "resolved_state": resolved_state,
        "confidence_tier": confidence,
        "call_state": call_state,
        "call_status": call_status,
        "rna_dosage_interpretation": interpretation,
        "confounding_flags": join_flags(flags),
    }


__all__ = [
    "ASE_DIRECTIONS", "ASE_DIRECTION_STATES", "AXES_TO_STATE",
    "BASELINE_LEVELS", "CONFIDENCE_TIERS", "COPY_STATES", "EVIDENCE_CLASSES",
    "EXACT_STATES", "EXPRESSION_DIRECTION_NULL_PATTERNS",
    "EXPRESSION_PATTERNS", "EXPRESSION_PATTERN_SHIFTS",
    "AseBiasFit", "FixedMixtureFit", "GAIN_SHIFT", "LOSS_SHIFT", "PROGRAM_VERSION",
    "RobustBivariateFit", "STATE_AXES", "STATE_LOG_ODDS_SHIFT",
    "ambient_adjusted_fraction", "arm_side", "ase_direction_statistics",
    "ase_state_log_bfs", "best_expression_component", "bh_adjust",
    "biological_block", "biological_crossfit_fold", "bounded", "by_adjust",
    "canonical_chromosome", "chromosome_from_arm", "classify_pq_relationship",
    "component_log_likelihoods", "conjunction_pvalues",
    "contamination_quadrature", "copy_and_donor_for_state",
    "direction_from_log_bfs", "empirical_p_floor", "empirical_upper_tail",
    "expression_direction_log_bfs", "fit_fixed_center_mixture",
    "fit_ase_bias_component", "fit_ase_delta",
    "format_number", "inverse_ambient_logit", "is_autosomal_chromosome",
    "join_flags", "logistic", "logit", "logsumexp",
    "multivariate_student_log_density", "quasi_ase_log_likelihood",
    "present_identifier", "resolve_evidence", "robust_bivariate_fit",
    "split_flags", "stable_fold",
    "stable_gene_fold", "student_t_log_density", "validate_hybrid_table",
    "weighted_median",
]
