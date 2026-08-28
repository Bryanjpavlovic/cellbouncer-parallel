#!/usr/bin/env python3
"""Aggregate probabilities only for final supported reconciliation events.

The input contract is intentionally narrow.  Every scored cell has exactly two
biological hypotheses:

1. the original top assignment from the allowed/constrained demux roster; and
2. the exact alternative identity nominated by identity reconciliation.

The supported-event report is restricted to final ``REASSIGN_GENOTYPE`` cells
in the established v10 population-supported sample-swap classes and attached
``REVIEW_CELLULAR_ORIGIN`` context cells.  Applied, recommended-not-applied,
and held cells are summarized separately; held cells never vote in the
authoritative reassignment verdict.  Unexpected-identity review rows are
written separately.  Generic constrained
runners, unconstrained winners/runners, unmade donor combinations, technical
multiplets, keep/unresolved/conflicted decisions, occupancy/ploidy routes, and
standalone review cells are excluded.  For the initial frozen library 19
evaluation, only site-primary scorer rows with site-aligned QC may vote;
molecule-primary rows remain diagnostic and unevaluable.  Legacy total-
likelihood percentages are retained for audit compatibility but are explicitly
non-headline because they saturate with accumulated evidence.  Experimental
relative-fit scores normalize the site log-likelihood difference by common
discriminating evidence depth; they compare average fit between the two fixed
candidates and are not calibrated real-world identity probabilities.  This
program writes tables only; plotting remains a separate step.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import hashlib
import math
import os
import shutil
import statistics
from collections import Counter, defaultdict
from pathlib import Path


SCHEMA = (
    "identity_probability_aggregate_v6_3_experimental_relative_fit_"
    "scoped_sample_swap_events"
)
PAIR_SCHEMA = "identity_pair_probability_v3_reconciliation_targeted"
PAIR_MANIFEST_SCHEMA = "identity_reconciliation_score_pair_manifest_v2"
LEGACY_PAIR_MANIFEST_SCHEMA = "identity_reconciliation_score_pair_manifest_v1"
ACCEPTED_PAIR_MANIFEST_SCHEMAS = {
    PAIR_MANIFEST_SCHEMA,
    LEGACY_PAIR_MANIFEST_SCHEMA,
}
PAIR_CONTRACT = (
    "ORIGINAL_ALLOWED_DEMUX_VS_RECONCILIATION_NOMINATED_SWAP_ONLY"
)
ORIGINAL_ROLE = "ORIGINAL_ALLOWED_DEMUX"
PROPOSED_ROLE = "RECONCILIATION_NOMINATED_SWAP"
PAIR_COMPARISON = "reconciliation_swap"
AGGREGATE_SCOPE_CONTRACT = (
    "FINAL_V10_SAMPLE_SWAP_EVENTS_SEPARATE_APPLIED_RECOMMENDED_AND_HELD_"
    "FROZEN_SITE_ONLY_QC_PASS"
)
APPLIED_REASSIGNMENT = "APPLIED_REASSIGNMENT"
RECOMMENDED_NOT_APPLIED = "RECOMMENDED_NOT_APPLIED"
SUPPORTED_EVENT_HELD_CELL = "SUPPORTED_EVENT_HELD_CELL"
REVIEW_ONLY_UNEXPECTED_IDENTITY = "REVIEW_ONLY_UNEXPECTED_IDENTITY"
REASSIGNMENT_SCOPES = {
    APPLIED_REASSIGNMENT,
    RECOMMENDED_NOT_APPLIED,
}
HEADLINE_POPULATION_SCOPES = {
    *REASSIGNMENT_SCOPES,
    SUPPORTED_EVENT_HELD_CELL,
}
SUPPORTED_SWAP_EVENT_CLASSES = {
    "LIKELY_UNEXPECTED_INTACT_BIOLOGICAL_LINE",
    "LIKELY_UNEXPECTED_SINGLET_POPULATION",
}
PREFERRED_AMBIENT_CONDITION = "IND_CK_RF_SX0_GATED_RFREE_PFIT"
PROBABILITY_SUM_TOLERANCE_PP = 1e-4
TIE_TOLERANCE_PP = 1e-9
RELATIVE_FIT_BASIS = (
    "SITE_DELTA_LOG_LIKELIHOOD_A_MINUS_B_DIVIDED_BY_"
    "DISCRIMINATING_EVIDENCE_DEPTH"
)
RELATIVE_FIT_INTERPRETATION = (
    "EXPERIMENTAL_AVERAGE_PAIRWISE_FIT_BETWEEN_FIXED_CANDIDATES_"
    "NOT_PROBABILITY_IDENTITY_IS_CORRECT"
)
BIOLOGICAL_ADMISSIBILITY = {
    "SINGLET_IDENTITY_CANDIDATE",
    "BIOLOGICAL_SINGLE_CELL_ALLOWED",
}

PAIR_EVIDENCE_FIELDS = [
    "comparison_status",
    "probability_basis",
    "evidence_basis",
    "delta_log_likelihood_a_minus_b",
    "site_delta_log_likelihood_a_minus_b",
    "site_candidate_a_probability_pct",
    "site_candidate_b_probability_pct",
    "molecule_balanced_delta_log_likelihood_a_minus_b",
    "molecule_balanced_candidate_a_probability_pct",
    "molecule_balanced_candidate_b_probability_pct",
    "shared_donor_components",
    "n_common_observed_snps",
    "n_discriminating_snps",
    "n_nondiscriminating_snps",
    "common_evidence_depth",
    "discriminating_evidence_depth",
    "n_snps_favor_preferred",
    "n_snps_favor_alternative",
    "n_snps_neutral",
    "genotype_model_similarity_pct",
    "chromosomes_covered",
    "effective_independent_snps",
    "maximum_site_contribution_fraction",
    "top_five_site_contribution_fraction",
    "n_independent_molecules",
    "effective_independent_molecules",
    "maximum_molecule_contribution_fraction",
    "molecule_umi_gene_fraction",
    "molecule_evidence_status",
    "preferred_residual_mismatch",
    "alternative_residual_mismatch",
    "absolute_fit_status",
    "probability_without_top_site_pct",
    "probability_without_top_five_sites_pct",
    "probability_without_top_molecule_pct",
    "minimum_leave_one_chromosome_out_probability_pct",
    "winner_changed_after_influence_removal",
    "site_bootstrap_win_fraction",
    "downsample_50pct_win_fraction",
    "downsample_basis",
    "minimum_error_sensitivity_probability_pct",
    "error_sensitivity_stable",
    "ambient_sensitivity_status",
    "error_ref",
    "error_alt",
    "warnings",
]

CELL_FIELDS = [
    "library", "barcode", "score_pair_id",
    "score_population_scope", "supported_event_key",
    "original_allowed_demux_assignment", "reconciliation_nominated_swap",
    "headline_probability_status", "headline_probability_exclusion_reason",
    "primary_probability_basis", "headline_qc_basis_alignment",
    "site_top_site_removal_status", "site_top_five_removal_status",
    "combined_influence_removal_flag_role", "comparison_outcome",
    "scorer_candidate_a_probability_pct", "scorer_candidate_b_probability_pct",
    "scorer_preferred_assignment",
    "winning_assignment", "winning_probability_pct", "losing_assignment",
    "losing_probability_pct", "original_assignment_probability_pct",
    "proposed_swap_probability_pct", "proposed_swap_preferred",
    "probability_gap_percentage_points", "alternative_closeness_out_of_100",
    "original_assignment_relative_fit_score_out_of_100",
    "proposed_swap_relative_fit_score_out_of_100",
    "proposed_minus_original_relative_fit_score_points",
    "relative_fit_score_status", "relative_fit_score_basis",
    "relative_fit_score_interpretation",
    "probability_interpretation", "reconciliation_event_id",
    "reconciliation_event_class", "reconciliation_event_confidence",
    "reconciliation_final_action", "reconciliation_decision_confidence",
    "reconciliation_reassignment_applied",
    "reconciliation_current_refined_assignment", "source_identity",
    "replaced_component", "replacement_component", "preserved_partner",
    "original_project_genotype_status", "original_biological_admissibility",
    "proposed_project_genotype_status", "proposed_biological_admissibility",
    "snps_favor_original", "snps_favor_proposed_swap",
    "original_residual_mismatch", "proposed_swap_residual_mismatch",
    "similarity_explanation", "evidence_summary",
    *PAIR_EVIDENCE_FIELDS,
    "ambient_original_c", "ambient_reconciled_c",
    "ambient_reconciled_minus_original_c", "ambient_stratum",
    "ambient_assignment_diagnostic", "score_scope_contract",
    "aggregate_scope_contract", "schema_version",
]

EVENT_METRIC_FIELDS = [
    "reconciliation_event_class",
    "reconciliation_event_confidence", "n_cells_scored",
    "n_cells_in_authoritative_verdict", "authoritative_reassignment_scope",
    "n_applied_reassignment_cells", "n_recommended_not_applied_cells",
    "n_supported_reassignment_cells", "n_supported_event_held_cells",
    "n_review_only_cells", "score_population_scope_counts",
    "n_cells_probability_available", "n_cells_proposed_swap_preferred",
    "n_cells_original_preferred", "n_cells_tied",
    "n_cells_probability_unavailable", "probability_exclusion_reason_counts",
    "fraction_cells_proposed_swap_preferred",
    "median_proposed_swap_probability_pct",
    "proposed_swap_probability_q10_pct",
    "proposed_swap_probability_q25_pct",
    "proposed_swap_probability_q75_pct",
    "proposed_swap_probability_q90_pct",
    "median_original_probability_pct", "median_winning_probability_pct",
    "median_alternative_closeness_out_of_100",
    "n_cells_relative_fit_available",
    "median_proposed_swap_relative_fit_score_out_of_100",
    "proposed_swap_relative_fit_score_q10_out_of_100",
    "proposed_swap_relative_fit_score_q25_out_of_100",
    "proposed_swap_relative_fit_score_q75_out_of_100",
    "proposed_swap_relative_fit_score_q90_out_of_100",
    "median_original_assignment_relative_fit_score_out_of_100",
    "median_proposed_minus_original_relative_fit_score_points",
    "median_discriminating_snps", "median_effective_independent_snps",
    "median_independent_molecules", "fraction_molecule_evidence_available",
    "median_maximum_site_contribution_fraction",
    "median_top_five_site_contribution_fraction",
    "site_top_site_removal_status_counts",
    "site_top_five_removal_status_counts",
    "median_downsample_50pct_win_fraction",
    "median_probability_without_top_site_pct",
    "median_probability_without_top_five_sites_pct",
    "median_minimum_error_sensitivity_probability_pct",
    "fraction_error_sensitivity_stable", "absolute_fit_status_counts",
    "warning_counts", "final_action_counts", "reassignment_applied_counts",
    "n_cells_with_controlled_ambient", "median_ambient_change_c",
    "event_result", "event_evidence_summary", "score_scope_contract",
    "aggregate_scope_contract", "schema_version",
]

EVENT_FIELDS = [
    "library", "reconciliation_event_id", "reconciliation_nominated_swap",
    "original_assignments", "n_original_assignment_strata",
    *EVENT_METRIC_FIELDS,
]

EVENT_STRATUM_FIELDS = [
    "library", "reconciliation_event_id", "original_allowed_demux_assignment",
    "reconciliation_nominated_swap", *EVENT_METRIC_FIELDS,
]

EVENT_SCOPE_FIELDS = [
    "library", "reconciliation_event_id", "original_allowed_demux_assignment",
    "reconciliation_nominated_swap", "score_population_scope",
    *EVENT_METRIC_FIELDS,
]

CROSS_LIBRARY_FIELDS = [
    "reconciliation_nominated_swap", "n_libraries_nominated",
    "libraries_nominated", "n_reconciliation_events", "event_ids",
    "original_assignments", "n_cells_scored",
    "n_cells_probability_available", "fraction_cells_proposed_swap_preferred",
    "median_proposed_swap_probability_pct",
    "median_proposed_swap_relative_fit_score_out_of_100",
    "median_original_assignment_relative_fit_score_out_of_100",
    "event_class_counts",
    "libraries_where_identity_is_expected", "other_expected_replicate_libraries",
    "cross_library_status", "score_scope_contract",
    "aggregate_scope_contract", "schema_version",
]

DISTRIBUTION_FIELDS = [
    "library", "reconciliation_event_id", "original_allowed_demux_assignment",
    "reconciliation_nominated_swap", "score_population_scope",
    "proposed_probability_bin_start_pct",
    "proposed_probability_bin_end_pct", "n_cells", "schema_version",
]

SCOPE_SUMMARY_FIELDS = [
    "library", "n_source_pair_scores", "n_supported_event_keys",
    "n_applied_reassignment_cells", "n_recommended_not_applied_cells",
    "n_supported_reassignment_cells", "n_supported_event_held_cells",
    "n_headline_supported_cells", "n_review_only_unexpected_identity_cells",
    "n_excluded_after_scoring", "n_headline_supported_events",
    "n_review_only_events", "aggregate_scope_contract", "schema_version",
]

PROVENANCE_FIELDS = [
    "library", "validation_summary", "validation_summary_sha256",
    "decision_table", "decision_table_sha256", "pair_manifest",
    "pair_manifest_sha256", "pair_summary", "pair_summary_sha256",
    "probability_output", "probability_output_sha256",
    "probability_provenance", "probability_provenance_sha256",
    "pair_manifest_schema", "pair_builder_sha256",
    "metadata_expected_genotypes_sha256",
    "metadata_global_biological_lines_sha256", "metadata_global_donors_sha256",
    "aggregator_sha256", "scorer_binary", "scorer_binary_sha256",
    "software_hash_status", "source_pair_scores_reused", "schema_version",
]


def open_text(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return open(path, mode, newline="")


def clean(value, default="NA"):
    text = "" if value is None else str(value).strip()
    return default if text in {"", ".", "NA", "nan", "NaN", "None"} else text


def fnum(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return math.nan
    return result if math.isfinite(result) else math.nan


def fint(value, default=0):
    number = fnum(value)
    return int(number) if math.isfinite(number) else default


def is_true(value):
    return clean(value, "FALSE").upper() in {"TRUE", "1", "YES", "Y"}


def fmt(value, digits=6):
    number = fnum(value)
    return "NA" if not math.isfinite(number) else f"{number:.{digits}g}"


def stable_logistic_percent(value):
    """Return 100 * logistic(value) without exponential overflow."""
    number = fnum(value)
    if not math.isfinite(number):
        return math.nan
    if number >= 0.0:
        return 100.0 / (1.0 + math.exp(-number))
    exp_number = math.exp(number)
    return 100.0 * exp_number / (1.0 + exp_number)


def relative_fit_scores(row):
    """Compute depth-normalized site relative fit for original A vs proposal B.

    Candidate A is contractually the original allowed demux assignment and
    candidate B is the reconciliation-nominated swap.  The scorer reports
    ``site_delta_log_likelihood_a_minus_b``; negating it orients positive
    values toward the proposed candidate.  Dividing by discriminating evidence
    depth removes the direct accumulation with evidence quantity.  The two
    returned scores sum to 100 and describe average pairwise fit only.
    """
    delta_a_minus_b = fnum(row.get("site_delta_log_likelihood_a_minus_b"))
    evidence_depth = fnum(row.get("discriminating_evidence_depth"))
    if not math.isfinite(delta_a_minus_b):
        return {
            "original": math.nan,
            "proposed": math.nan,
            "gap": math.nan,
            "status": "UNAVAILABLE_SITE_DELTA_LOG_LIKELIHOOD",
        }
    if not math.isfinite(evidence_depth) or evidence_depth <= 0.0:
        return {
            "original": math.nan,
            "proposed": math.nan,
            "gap": math.nan,
            "status": "UNAVAILABLE_NONPOSITIVE_DISCRIMINATING_EVIDENCE_DEPTH",
        }
    proposed_mean_advantage = -delta_a_minus_b / evidence_depth
    proposed_score = stable_logistic_percent(proposed_mean_advantage)
    original_score = 100.0 - proposed_score
    return {
        "original": original_score,
        "proposed": proposed_score,
        "gap": proposed_score - original_score,
        "status": "AVAILABLE",
    }


def relative_fit_available(row):
    return (
        clean(row.get("relative_fit_score_status"), "").upper() == "AVAILABLE"
        and math.isfinite(fnum(
            row.get("original_assignment_relative_fit_score_out_of_100")
        ))
        and math.isfinite(fnum(
            row.get("proposed_swap_relative_fit_score_out_of_100")
        ))
    )


def canonical_genotype(value):
    text = clean(value, "")
    if not text:
        return ""
    for prefix in ("D[", "T[", "UNKNOWN_SINGLE_CELL["):
        if text.startswith(prefix) and text.endswith("]"):
            text = text[len(prefix):-1]
            break
    if text.startswith("M{"):
        return text
    parts = [part.strip() for part in text.replace("x", "+").split("+")]
    return "+".join(sorted(part for part in parts if part))


def donor_components(value):
    genotype = canonical_genotype(value)
    if not genotype or genotype.startswith("M{"):
        return []
    return [part for part in genotype.split("+") if part]


def parse_libraries(tokens):
    values = set()
    for token in tokens:
        for part in str(token).split(","):
            part = part.strip().lower()
            if part.startswith("lib"):
                part = part[3:]
            if not part:
                continue
            if "-" in part:
                start, end = (int(value) for value in part.split("-", 1))
                values.update(range(min(start, end), max(start, end) + 1))
            else:
                values.add(int(part))
    return sorted(values)


def read_table(path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"missing TSV header: {path}")
        return list(reader.fieldnames), [dict(row) for row in reader]


def require_columns(header, required, path):
    missing = sorted(set(required) - set(header))
    if missing:
        raise ValueError(
            f"{path}: required columns missing: " + ", ".join(missing)
        )


def write_table(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(
        path.name + f".tmp.{os.getpid()}" +
        (".gz" if path.name.endswith(".gz") else "")
    )
    with open_text(temporary, "wt") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "NA") for field in fields})
    os.replace(temporary, path)


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def acquire_output_lock(output_root):
    """Hold a nonblocking process lock for the complete multi-file write."""
    lock_path = Path(output_root) / ".identity_probability_aggregate.lock"
    handle = open(lock_path, "a+")
    try:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        handle.close()
        raise RuntimeError(
            f"another identity probability aggregate is writing {output_root}"
        ) from exc
    handle.seek(0)
    handle.truncate()
    handle.write(f"pid={os.getpid()}\n")
    handle.flush()
    return handle


def optional_bool(value):
    text = clean(value, "").upper()
    if text in {"TRUE", "1", "YES", "Y"}:
        return True
    if text in {"FALSE", "0", "NO", "N"}:
        return False
    return None


def validate_reconciliation_summary(identity_root):
    path = Path(identity_root) / "validation" / "validation_summary.tsv"
    if not path.is_file():
        raise FileNotFoundError(
            f"validated reconciliation summary is missing: {path}"
        )
    header, rows = read_table(path)
    require_columns(header, {"check", "status", "n_failures"}, path)
    if not rows:
        raise ValueError(f"{path}: validation summary is empty")
    failures = [
        row for row in rows
        if clean(row.get("status"), "").upper() != "PASS"
        or fint(row.get("n_failures"), -1) != 0
    ]
    if failures:
        preview = ", ".join(
            f"{clean(row.get('check'))}:{clean(row.get('status'))}"
            for row in failures[:5]
        )
        raise ValueError(
            f"{path}: reconciliation validation is not PASS ({preview})"
        )
    return path, sha256_file(path)


def median(values):
    finite = [fnum(value) for value in values if math.isfinite(fnum(value))]
    return statistics.median(finite) if finite else math.nan


def quantile(values, fraction):
    finite = sorted(fnum(value) for value in values if math.isfinite(fnum(value)))
    if not finite:
        return math.nan
    if len(finite) == 1:
        return finite[0]
    position = fraction * (len(finite) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return finite[lower]
    weight = position - lower
    return finite[lower] * (1.0 - weight) + finite[upper] * weight


def fraction(numerator, denominator):
    return "NA" if denominator == 0 else fmt(numerator / denominator)


def counter_text(values):
    counts = Counter(clean(value) for value in values)
    return ",".join(f"{key}:{counts[key]}" for key in sorted(counts)) or "NONE"


def interpretive_warning_tokens(row):
    warnings = clean(row.get("warnings"), "NONE")
    if warnings == "NONE":
        return []
    return [
        warning.strip() for warning in warnings.split(",")
        if warning.strip() and warning.strip() != "INFLUENCE_REMOVAL_UNSTABLE"
    ]


def warning_counter_text(rows):
    """Summarize only warnings that are valid for the headline QC basis.

    The scorer's combined influence warning includes leave-one-contig behavior,
    so it is retained in cell-level provenance but cannot enter event-level
    interpretation for this approximately 4,000-sequence reference.
    """
    counts = Counter()
    for row in rows:
        for warning in interpretive_warning_tokens(row):
            counts[warning] += 1
    return ",".join(f"{key}:{counts[key]}" for key in sorted(counts)) or "NONE"


def newest_named_file(root, basename):
    if not root or not Path(root).is_dir():
        return None
    matches = [path for path in Path(root).rglob(basename) if path.is_file()]
    return max(matches, key=lambda path: (path.stat().st_mtime_ns, str(path))) if matches else None


def load_ambient(args, libraries):
    selected = set(libraries)
    contrast_path = (
        Path(args.ambient_contrast_cells)
        if args.ambient_contrast_cells else
        newest_named_file(args.ambient_root, "reconciliation_planned_contrast_cells.tsv")
    )
    if not contrast_path or not contrast_path.is_file():
        return {}, [], {
            "status": "UNAVAILABLE", "condition": "NA",
            "contrast_path": "NA", "donor_burden_path": "NA",
        }
    header, rows = read_table(contrast_path)
    require_columns(
        header,
        {"library", "condition", "contrast", "barcode", "left_rate", "right_rate"},
        contrast_path,
    )
    conditions = sorted({
        clean(row.get("condition")) for row in rows
        if fint(row.get("library"), -1) in selected
    })
    if args.ambient_condition:
        condition = args.ambient_condition
        if condition not in conditions:
            raise ValueError(
                f"ambient condition {condition!r} absent from {contrast_path}"
            )
    elif PREFERRED_AMBIENT_CONDITION in conditions:
        condition = PREFERRED_AMBIENT_CONDITION
    else:
        condition = conditions[0] if conditions else "NA"

    cell_values = {}
    for row in rows:
        library = fint(row.get("library"), -1)
        if library not in selected or clean(row.get("condition")) != condition:
            continue
        if clean(row.get("contrast")) != "combined_production_change":
            continue
        barcode = clean(row.get("barcode"))
        left = fnum(row.get("left_rate"))
        right = fnum(row.get("right_rate"))
        if barcode != "NA" and math.isfinite(left) and math.isfinite(right):
            cell_values[(library, barcode)] = (
                left, right, right - left, clean(row.get("stratum"), "UNAVAILABLE")
            )

    burden_path = contrast_path.parent / "reconciliation_exact_donor_burden.tsv"
    if not burden_path.is_file():
        burden_path = newest_named_file(
            args.ambient_root, "reconciliation_exact_donor_burden.tsv"
        )
    burden_rows = []
    if burden_path and burden_path.is_file():
        _, source_rows = read_table(burden_path)
        burden_rows = [
            {**row, "schema_version": SCHEMA}
            for row in source_rows
            if fint(row.get("library"), -1) in selected
            and clean(row.get("condition")) == condition
        ]
    return cell_values, burden_rows, {
        "status": "AVAILABLE",
        "condition": condition,
        "contrast_path": str(contrast_path),
        "donor_burden_path": str(burden_path) if burden_path else "NA",
    }


def load_expected_identity_libraries(identity_root):
    path = Path(identity_root) / "metadata" / "library_expected_genotypes.tsv"
    if not path.is_file():
        return defaultdict(set)
    header, rows = read_table(path)
    require_columns(header, {"library", "canonical_genotype"}, path)
    result = defaultdict(set)
    for row in rows:
        genotype = canonical_genotype(row.get("canonical_genotype"))
        library = clean(row.get("library"))
        if genotype and library != "NA":
            result[genotype].add(library)
    return result


def load_project_biological_universe(identity_root):
    metadata_root = Path(identity_root) / "metadata"
    line_path = metadata_root / "global_biological_lines.tsv"
    donor_path = metadata_root / "global_donors.tsv"
    if not line_path.is_file() or not donor_path.is_file():
        raise FileNotFoundError(
            "targeted score validation requires global_biological_lines.tsv "
            "and global_donors.tsv"
        )
    line_header, line_rows = read_table(line_path)
    donor_header, donor_rows = read_table(donor_path)
    require_columns(line_header, {"canonical_genotype"}, line_path)
    require_columns(donor_header, {"donor_id"}, donor_path)
    lines = {
        canonical_genotype(row.get("canonical_genotype"))
        for row in line_rows
        if canonical_genotype(row.get("canonical_genotype"))
    }
    donors = {
        clean(row.get("donor_id"), "")
        for row in donor_rows
        if clean(row.get("donor_id"), "")
    }
    return lines, donors


def is_real_biological_identity(genotype, global_lines, global_donors):
    components = donor_components(genotype)
    if len(components) == 1:
        return components[0] in global_donors
    if len(components) == 2:
        return canonical_genotype(genotype) in global_lines
    return False


def validate_probability_provenance(
        path, manifest_path, pair_summary_path, probability_path, scorer_binary,
        allow_legacy_missing):
    if not path.is_file():
        if allow_legacy_missing:
            return {}, "LEGACY_PROVENANCE_NOT_RECORDED", "NA"
        raise FileNotFoundError(
            f"probability provenance sidecar is missing: {path}"
        )
    header, rows = read_table(path)
    require_columns(
        header,
        {
            "pair_manifest_sha256", "probability_output_sha256",
            "pair_summary_sha256",
            "scorer_binary_sha256", "counts_path", "counts_sha256",
            "samples_path", "samples_sha256", "assignments_path",
            "assignments_sha256", "pileup_sites_path",
            "pileup_sites_sha256", "pileup_observations_path",
            "pileup_observations_sha256", "pileup_molecules_path",
            "pileup_molecules_sha256",
        },
        path,
    )
    if len(rows) != 1:
        raise ValueError(f"{path}: expected one probability provenance row")
    row = rows[0]
    observed_manifest = sha256_file(manifest_path)
    observed_probability = sha256_file(probability_path)
    if clean(row.get("pair_manifest_sha256"), "") != observed_manifest:
        raise ValueError(f"{path}: probability was not scored from this pair manifest")
    if clean(row.get("pair_summary_sha256"), "") != sha256_file(pair_summary_path):
        raise ValueError(f"{path}: probability pair-summary hash mismatch")
    if clean(row.get("probability_output_sha256"), "") != observed_probability:
        raise ValueError(f"{path}: probability output hash mismatch")
    for prefix in (
        "counts", "samples", "assignments", "pileup_sites",
        "pileup_observations", "pileup_molecules",
    ):
        input_path = clean(row.get(prefix + "_path"), "")
        recorded_hash = clean(row.get(prefix + "_sha256"), "")
        if not input_path:
            continue
        input_file = Path(input_path)
        observed_hash = (
            sha256_file(input_file)
            if input_file.is_file() and input_file.stat().st_size > 0 else "NA"
        )
        if recorded_hash != observed_hash:
            raise ValueError(
                f"{path}: scored input changed after probability generation: "
                f"{input_path}"
            )
    recorded_scorer = clean(row.get("scorer_binary_sha256"), "")
    if scorer_binary and Path(scorer_binary).is_file():
        current_scorer = sha256_file(scorer_binary)
        software_status = (
            "MATCH" if recorded_scorer == current_scorer
            else "WARNING_CURRENT_SCORER_DIFFERS_FROM_FROZEN_RUN"
        )
    else:
        software_status = "CURRENT_SCORER_UNAVAILABLE_HASH_NOT_BLOCKING"
    return row, software_status, sha256_file(path)


def same_nuclear_donor_identity(left, right):
    left_components = donor_components(left)
    right_components = donor_components(right)
    return (
        bool(left_components)
        and bool(right_components)
        and len(set(left_components)) == 1
        and len(set(right_components)) == 1
        and left_components[0] == right_components[0]
    )


def load_pair_manifest(
        path, library, expected_libraries, global_lines, global_donors):
    header, rows = read_table(path)
    required = {
        "library", "barcode", "hypothesis_id", "donor_genotype",
        "current_donor_genotype", "candidate_origin",
        "project_genotype_status", "biological_admissibility",
        "expected_genotype_status", "score_pair_id", "score_pair_role",
        "reconciliation_event_id", "reconciliation_event_class",
        "reconciliation_event_confidence", "reconciliation_final_action",
        "reconciliation_decision_confidence",
        "reconciliation_reassignment_applied",
        "reconciliation_current_refined_assignment",
        "original_demux_assignment", "proposed_donor_genotype",
        "score_scope_contract", "schema_version",
    }
    require_columns(header, required, path)
    by_barcode = defaultdict(dict)
    pair_id_to_barcode = {}
    for row in rows:
        barcode = clean(row.get("barcode"))
        role = clean(row.get("score_pair_role"))
        if clean(row.get("library")) != f"lib{library}":
            raise ValueError(f"{path}: row has the wrong library for {barcode}")
        if barcode == "NA":
            raise ValueError(f"{path}: blank barcode")
        if role not in {ORIGINAL_ROLE, PROPOSED_ROLE}:
            raise ValueError(f"{path}: forbidden score_pair_role {role!r}")
        if role in by_barcode[barcode]:
            raise ValueError(f"{path}: duplicate role {role} for {barcode}")
        if clean(row.get("score_scope_contract")) != PAIR_CONTRACT:
            raise ValueError(f"{path}: score-scope contract mismatch for {barcode}")
        row_schema = clean(row.get("schema_version"))
        if row_schema not in ACCEPTED_PAIR_MANIFEST_SCHEMAS:
            raise ValueError(
                f"{path}: unsupported pair-manifest schema {row_schema!r}"
            )
        pair_id = clean(row.get("score_pair_id"))
        previous_barcode = pair_id_to_barcode.get(pair_id)
        if previous_barcode is not None and previous_barcode != barcode:
            raise ValueError(f"{path}: duplicate score_pair_id {pair_id}")
        pair_id_to_barcode[pair_id] = barcode
        by_barcode[barcode][role] = row

    result = {}
    for barcode, role_rows in by_barcode.items():
        if set(role_rows) != {ORIGINAL_ROLE, PROPOSED_ROLE}:
            raise ValueError(f"{path}: incomplete targeted pair for {barcode}")
        original = role_rows[ORIGINAL_ROLE]
        proposed = role_rows[PROPOSED_ROLE]
        pair_id = clean(original.get("score_pair_id"))
        if pair_id != clean(proposed.get("score_pair_id")):
            raise ValueError(f"{path}: inconsistent score_pair_id for {barcode}")
        if clean(original.get("schema_version")) != clean(
                proposed.get("schema_version")):
            raise ValueError(f"{path}: inconsistent manifest schema for {barcode}")
        original_identity = canonical_genotype(original.get("donor_genotype"))
        proposed_identity = canonical_genotype(proposed.get("donor_genotype"))
        if original_identity != canonical_genotype(original.get("original_demux_assignment")):
            raise ValueError(f"{path}: candidate A is not the original demux identity for {barcode}")
        if proposed_identity != canonical_genotype(proposed.get("proposed_donor_genotype")):
            raise ValueError(f"{path}: candidate B is not the reconciliation proposal for {barcode}")
        if original_identity == proposed_identity:
            raise ValueError(f"{path}: identical original/proposed identities for {barcode}")
        if same_nuclear_donor_identity(original_identity, proposed_identity):
            raise ValueError(
                f"{path}: A versus A+A is not nuclear-scoreable for {barcode}"
            )
        library_name = f"lib{library}"
        if library_name not in expected_libraries.get(original_identity, set()):
            raise ValueError(
                f"{path}: original identity is absent from the allowed library "
                f"roster for {barcode}: {original_identity}"
            )
        if library_name in expected_libraries.get(proposed_identity, set()):
            raise ValueError(
                f"{path}: proposed identity is already allowed in the library "
                f"roster for {barcode}: {proposed_identity}"
            )
        if not is_real_biological_identity(
                original_identity, global_lines, global_donors):
            raise ValueError(
                f"{path}: original identity is not a real project biological "
                f"identity for {barcode}: {original_identity}"
            )
        if not is_real_biological_identity(
                proposed_identity, global_lines, global_donors):
            raise ValueError(
                f"{path}: proposed swap is not a real project biological line "
                f"for {barcode}: {proposed_identity}"
            )
        for role, row in role_rows.items():
            if clean(row.get("biological_admissibility")) not in BIOLOGICAL_ADMISSIBILITY:
                raise ValueError(
                    f"{path}: non-biological {role} candidate escaped the pair gate for {barcode}"
                )
            if canonical_genotype(row.get("donor_genotype")).startswith("M{"):
                raise ValueError(f"{path}: technical multiplet escaped the pair gate for {barcode}")
        if clean(proposed.get("expected_genotype_status")) == "EXPECTED":
            raise ValueError(f"{path}: expected-roster alternative escaped swap-only gate for {barcode}")
        if not clean(proposed.get("reconciliation_event_id"), ""):
            raise ValueError(f"{path}: proposed swap has no reconciliation event for {barcode}")
        result[barcode] = {
            "score_pair_id": pair_id,
            "original": original,
            "proposed": proposed,
            "original_identity": original_identity,
            "proposed_identity": proposed_identity,
            "manifest_schema": clean(original.get("schema_version")),
        }
    return result


def pair_event_key(pair):
    event_id = clean(
        pair["proposed"].get("reconciliation_event_id"), ""
    )
    proposed = pair["proposed_identity"]
    return event_id, proposed


def pair_event_key_text(library, event_key):
    event_id, proposed = event_key
    return f"lib{library}|{event_id}|{proposed}"


def aggregate_nonidentity_event_reason(pair):
    proposed = pair["proposed"]
    event_class = clean(
        proposed.get("reconciliation_event_class"), ""
    ).upper()
    action = clean(
        proposed.get("reconciliation_final_action"), ""
    ).upper()
    if "HOMOTET" in event_class or "OCCUPANCY" in event_class:
        return "HOMOTET_OCCUPANCY_REVIEW_NOT_IDENTITY_SWAP"
    if "PLOIDY" in event_class or action == "RECLASSIFY_PLOIDY":
        return "PLOIDY_RECLASSIFICATION_NOT_IDENTITY_SWAP"
    if event_class in {
        "NO_METADATA_PROBLEM",
        "EXPECTED_COMPOSITE_COMPONENT_SINGLET_POPULATION",
    }:
        return "NOT_A_SAMPLE_SWAP_RECONCILIATION_EVENT"
    return ""


def supported_event_keys_for_manifest(manifest):
    return {
        pair_event_key(pair)
        for pair in manifest.values()
        if clean(
            pair["proposed"].get("reconciliation_final_action"), ""
        ).upper() == "REASSIGN_GENOTYPE"
        and not aggregate_nonidentity_event_reason(pair)
        and clean(
            pair["proposed"].get("reconciliation_event_class"), ""
        ).upper() in SUPPORTED_SWAP_EVENT_CLASSES
        and all(pair_event_key(pair))
    }


def reassignment_population_scope(value, context):
    applied = optional_bool(value)
    if applied is True:
        return APPLIED_REASSIGNMENT
    if applied is False:
        return RECOMMENDED_NOT_APPLIED
    raise ValueError(
        f"{context}: reconciliation_reassignment_applied must be an explicit "
        "TRUE or FALSE value"
    )


def classify_pair_population_scope(pair, supported_event_keys):
    if aggregate_nonidentity_event_reason(pair):
        return ""
    action = clean(
        pair["proposed"].get("reconciliation_final_action"), ""
    ).upper()
    event_class = clean(
        pair["proposed"].get("reconciliation_event_class"), ""
    ).upper()
    event_key = pair_event_key(pair)
    if (
        action == "REASSIGN_GENOTYPE"
        and event_class in SUPPORTED_SWAP_EVENT_CLASSES
    ):
        return reassignment_population_scope(
            pair["proposed"].get("reconciliation_reassignment_applied"),
            f"{clean(pair['proposed'].get('library'))} "
            f"{clean(pair['proposed'].get('barcode'))}",
        )
    if (
        action == "REVIEW_CELLULAR_ORIGIN"
        and event_key in supported_event_keys
    ):
        return SUPPORTED_EVENT_HELD_CELL
    if action == "REVIEW_UNEXPECTED_IDENTITY":
        return REVIEW_ONLY_UNEXPECTED_IDENTITY
    return ""


def aggregate_scope_exclusion_reason(pair, supported_event_keys):
    proposed = pair["proposed"]
    action = clean(proposed.get("reconciliation_final_action"), "").upper()
    event_key = pair_event_key(pair)
    nonidentity_reason = aggregate_nonidentity_event_reason(pair)
    if nonidentity_reason:
        return nonidentity_reason
    event_class = clean(
        proposed.get("reconciliation_event_class"), ""
    ).upper()
    if (
        action == "REASSIGN_GENOTYPE"
        and event_class not in SUPPORTED_SWAP_EVENT_CLASSES
    ):
        return "REASSIGNMENT_NOT_A_SUPPORTED_SAMPLE_SWAP_EVENT"
    if action == "KEEP":
        return "FINAL_ACTION_KEEP_NOT_A_SUPPORTED_SWAP"
    if action == "REVIEW_CELLULAR_ORIGIN":
        if event_key not in supported_event_keys:
            return "STANDALONE_CELLULAR_ORIGIN_REVIEW_NOT_SUPPORTED_SWAP"
    if action == "KEEP_CURRENT_CONFLICTED":
        return "CONFLICTED_PROPOSAL_NOT_SUPPORTED_SWAP"
    if action == "UNRESOLVED_INSUFFICIENT_EVIDENCE":
        return "UNRESOLVED_PROPOSAL_NOT_SUPPORTED_SWAP"
    return "FINAL_ACTION_OUTSIDE_SUPPORTED_SWAP_SCOPE"


def validate_optional_scope_annotations(
        path, library, pair, population_scope, supported_event_key):
    proposed = pair["proposed"]
    annotated_scope = clean(proposed.get("score_population_scope"), "")
    annotated_key = clean(proposed.get("supported_event_key"), "")
    if annotated_scope == "SUPPORTED_REASSIGNMENT":
        annotated_scope = reassignment_population_scope(
            proposed.get("reconciliation_reassignment_applied"),
            f"{path} lib{library} {clean(proposed.get('barcode'))}",
        )
    # Legacy broad v1 manifests annotated excluded final actions before v6
    # narrowed the supported sample-swap population.  Their annotation is
    # retained as provenance but is not authoritative for an excluded row.
    if population_scope and annotated_scope and annotated_scope != population_scope:
        raise ValueError(
            f"{path}: score_population_scope mismatch for "
            f"{clean(proposed.get('barcode'))}: annotated={annotated_scope!r} "
            f"computed={population_scope!r}"
        )
    if population_scope and annotated_key and annotated_key != supported_event_key:
        raise ValueError(
            f"{path}: supported_event_key mismatch for "
            f"{clean(proposed.get('barcode'))}: annotated={annotated_key!r} "
            f"computed={supported_event_key!r}"
        )


def primary_site_probabilities(row):
    """Return site probabilities for the frozen site-only evaluation arm."""
    site_a = fnum(row.get("site_candidate_a_probability_pct"))
    site_b = fnum(row.get("site_candidate_b_probability_pct"))
    if math.isfinite(site_a) and math.isfinite(site_b):
        return site_a, site_b, "SITE_LIKELIHOOD", ""
    raw_a = fnum(row.get("candidate_a_probability_pct"))
    raw_b = fnum(row.get("candidate_b_probability_pct"))
    probability_basis = clean(row.get("probability_basis"), "").upper()
    if (
        "SITE" in probability_basis
        and math.isfinite(raw_a)
        and math.isfinite(raw_b)
    ):
        return raw_a, raw_b, "SITE_LIKELIHOOD_LEGACY_COLUMNS", ""
    return (
        math.nan, math.nan, "UNAVAILABLE",
        "SITE_PRIMARY_PROBABILITY_UNAVAILABLE",
    )


def probability_pair_well_formed(left, right):
    return (
        math.isfinite(left)
        and math.isfinite(right)
        and 0.0 <= left <= 100.0
        and 0.0 <= right <= 100.0
        and abs(math.fsum((left, right)) - 100.0)
        <= PROBABILITY_SUM_TOLERANCE_PP
    )


def scorer_primary_is_molecule(row):
    """V3 reorients fit/stability fields when molecule evidence is primary."""
    basis = clean(row.get("probability_basis"), "").upper()
    return "MOLECULE" in basis or molecule_evidence_available(row)


def site_removal_status(row, field):
    probability = fnum(row.get(field))
    if not math.isfinite(probability):
        return "UNAVAILABLE"
    if probability < 0.0 or probability > 100.0:
        return "INVALID_PROBABILITY"
    if probability <= 50.0:
        return "SITE_WINNER_NOT_PRESERVED"
    return "SITE_WINNER_PRESERVED"


def headline_probability_eligibility(row, p_original, p_proposed, basis):
    reasons = []
    if clean(row.get("comparison_status"), "").upper() != "PASS":
        reasons.append(
            "COMPARISON_STATUS_" + clean(row.get("comparison_status"), "MISSING").upper()
        )
    if not (math.isfinite(p_original) and math.isfinite(p_proposed)):
        reasons.append("SITE_PRIMARY_PROBABILITY_UNAVAILABLE")
    else:
        if not (0.0 <= p_original <= 100.0):
            reasons.append("SITE_CANDIDATE_A_PROBABILITY_OUT_OF_RANGE")
        if not (0.0 <= p_proposed <= 100.0):
            reasons.append("SITE_CANDIDATE_B_PROBABILITY_OUT_OF_RANGE")
        if abs(math.fsum((p_original, p_proposed)) - 100.0) > PROBABILITY_SUM_TOLERANCE_PP:
            reasons.append("SITE_PROBABILITIES_DO_NOT_SUM_TO_100")
    if basis == "UNAVAILABLE":
        reasons.append("SITE_PRIMARY_BASIS_UNAVAILABLE")
    if scorer_primary_is_molecule(row):
        reasons.append("SCORER_MOLECULE_PRIMARY_QC_NOT_SITE_ALIGNED")
    else:
        if clean(row.get("absolute_fit_status"), "").upper() != "PASS":
            reasons.append(
                "ABSOLUTE_FIT_STATUS_"
                + clean(row.get("absolute_fit_status"), "MISSING").upper()
            )
        top_site_status = site_removal_status(
            row, "probability_without_top_site_pct"
        )
        top_five_status = site_removal_status(
            row, "probability_without_top_five_sites_pct"
        )
        if top_site_status != "SITE_WINNER_PRESERVED":
            reasons.append("TOP_SITE_REMOVAL_" + top_site_status)
        if top_five_status != "SITE_WINNER_PRESERVED":
            reasons.append("TOP_FIVE_SITE_REMOVAL_" + top_five_status)
        error_stable = optional_bool(row.get("error_sensitivity_stable"))
        if error_stable is False:
            reasons.append("ERROR_SENSITIVITY_UNSTABLE")
    reasons = sorted(set(reasons))
    return (not reasons), (";".join(reasons) if reasons else "PASS")


def probability_available(row):
    return (
        clean(row.get("headline_probability_status"), "").upper() == "ELIGIBLE"
        and math.isfinite(fnum(row.get("original_assignment_probability_pct")))
        and math.isfinite(fnum(row.get("proposed_swap_probability_pct")))
    )


def molecule_evidence_available(row):
    status = clean(row.get("molecule_evidence_status"), "")
    return (
        status not in {
            "", "MOLECULE_SIDECAR_UNAVAILABLE",
            "MOLECULE_INDEPENDENCE_UNAVAILABLE",
        }
        and fint(row.get("n_independent_molecules")) > 0
    )


def molecule_evidence_phrase(row):
    if molecule_evidence_available(row):
        return (
            f"{clean(row.get('n_independent_molecules'))} "
            "independent molecules"
        )
    return "independent molecules not measured"


def similarity_explanation(row):
    parts = []
    shared = clean(row.get("shared_donor_components"), "NONE")
    if shared != "NONE":
        parts.append(f"shared donor component(s): {shared}")
    similarity = fnum(row.get("genotype_model_similarity_pct"))
    if math.isfinite(similarity):
        parts.append(f"observed genotype-model similarity {fmt(similarity)}%")
    parts.append(
        f"{clean(row.get('n_discriminating_snps'), '0')} discriminating SNPs"
    )
    parts.append(molecule_evidence_phrase(row))
    warnings = interpretive_warning_tokens(row)
    if warnings:
        parts.append("stability/concentration flags: " + ",".join(warnings))
    return "; ".join(parts)


def build_cell_row(
        library, barcode, manifest, probability, ambient,
        score_population_scope, supported_event_key):
    original_row = manifest["original"]
    proposed_row = manifest["proposed"]
    original = manifest["original_identity"]
    proposed = manifest["proposed_identity"]
    raw_original = fnum(probability.get("candidate_a_probability_pct"))
    raw_proposed = fnum(probability.get("candidate_b_probability_pct"))
    p_original, p_proposed, primary_basis, primary_reason = (
        primary_site_probabilities(probability)
    )
    eligible, eligibility_reason = headline_probability_eligibility(
        probability, p_original, p_proposed, primary_basis
    )
    if primary_reason and primary_reason not in eligibility_reason:
        eligibility_reason = (
            primary_reason if eligibility_reason == "PASS"
            else eligibility_reason + ";" + primary_reason
        )
        eligible = False
    molecule_primary = scorer_primary_is_molecule(probability)
    available = probability_pair_well_formed(p_original, p_proposed)
    tied = available and abs(p_original - p_proposed) <= TIE_TOLERANCE_PP
    if tied:
        winner = loser = "TIE"
        winning_probability = losing_probability = p_original
        comparison_outcome = "TIE"
    elif available and p_original > p_proposed:
        winner, loser = original, proposed
        winning_probability, losing_probability = p_original, p_proposed
        comparison_outcome = "ORIGINAL_PREFERRED"
    elif available:
        winner, loser = proposed, original
        winning_probability, losing_probability = p_proposed, p_original
        comparison_outcome = "PROPOSED_SWAP_PREFERRED"
    else:
        winner = loser = "NA"
        winning_probability = losing_probability = math.nan
        comparison_outcome = "UNAVAILABLE"

    relative_fit = relative_fit_scores(probability)
    if relative_fit["status"] == "AVAILABLE" and available:
        relative_gap = relative_fit["gap"]
        probability_direction = (
            1 if p_proposed > p_original else
            -1 if p_proposed < p_original else 0
        )
        relative_direction = (
            1 if relative_gap > TIE_TOLERANCE_PP else
            -1 if relative_gap < -TIE_TOLERANCE_PP else 0
        )
        if probability_direction != relative_direction:
            raise ValueError(
                f"lib{library} {barcode}: normalized relative-fit direction "
                "does not match the site-probability direction"
            )

    preferred = clean(probability.get("preferred_assignment"))
    preferred_identity = canonical_genotype(preferred)
    preferred_sites = clean(probability.get("n_snps_favor_preferred"), "0")
    alternative_sites = clean(probability.get("n_snps_favor_alternative"), "0")
    preferred_residual = clean(probability.get("preferred_residual_mismatch"))
    alternative_residual = clean(probability.get("alternative_residual_mismatch"))
    if preferred_identity == original:
        original_sites, proposed_sites = preferred_sites, alternative_sites
        original_residual, proposed_residual = preferred_residual, alternative_residual
    elif preferred_identity == proposed:
        original_sites, proposed_sites = alternative_sites, preferred_sites
        original_residual, proposed_residual = alternative_residual, preferred_residual
    else:
        original_sites = proposed_sites = "NA"
        original_residual = proposed_residual = "NA"

    left, right, delta, ambient_stratum = ambient or (
        math.nan, math.nan, math.nan, "UNAVAILABLE"
    )
    if math.isfinite(delta):
        ambient_diagnostic = (
            "LOWER_AFTER_RECONCILIATION" if delta < 0 else
            "HIGHER_AFTER_RECONCILIATION" if delta > 0 else "NO_CHANGE"
        )
    else:
        ambient_diagnostic = "UNAVAILABLE"

    row = {
        "library": library,
        "barcode": barcode,
        "score_pair_id": manifest["score_pair_id"],
        "score_population_scope": score_population_scope,
        "supported_event_key": supported_event_key,
        "original_allowed_demux_assignment": original,
        "reconciliation_nominated_swap": proposed,
        "headline_probability_status": "ELIGIBLE" if eligible else "UNEVALUABLE",
        "headline_probability_exclusion_reason": eligibility_reason,
        "primary_probability_basis": primary_basis,
        "headline_qc_basis_alignment": (
            "SCORER_MOLECULE_PRIMARY_QC_NOT_SITE_ALIGNED"
            if molecule_primary else "FROZEN_SITE_ONLY_QC_ALIGNED"
        ),
        "site_top_site_removal_status": (
            "NOT_APPLICABLE_SCORER_MOLECULE_PRIMARY"
            if molecule_primary else site_removal_status(
                probability, "probability_without_top_site_pct"
            )
        ),
        "site_top_five_removal_status": (
            "NOT_APPLICABLE_SCORER_MOLECULE_PRIMARY"
            if molecule_primary else site_removal_status(
                probability, "probability_without_top_five_sites_pct"
            )
        ),
        "combined_influence_removal_flag_role": (
            "NON_AUTHORITATIVE_INCLUDES_LEAVE_ONE_CONTIG"
        ),
        "comparison_outcome": comparison_outcome,
        "scorer_candidate_a_probability_pct": fmt(raw_original),
        "scorer_candidate_b_probability_pct": fmt(raw_proposed),
        "scorer_preferred_assignment": preferred,
        "winning_assignment": winner,
        "winning_probability_pct": fmt(winning_probability),
        "losing_assignment": loser,
        "losing_probability_pct": fmt(losing_probability),
        "original_assignment_probability_pct": fmt(p_original),
        "proposed_swap_probability_pct": fmt(p_proposed),
        "proposed_swap_preferred": (
            "TIE" if tied else
            "TRUE" if available and p_proposed > p_original else
            "FALSE" if available else "UNAVAILABLE"
        ),
        "probability_gap_percentage_points": fmt(
            winning_probability - losing_probability
        ),
        "alternative_closeness_out_of_100": fmt(2.0 * losing_probability),
        "original_assignment_relative_fit_score_out_of_100": fmt(
            relative_fit["original"]
        ),
        "proposed_swap_relative_fit_score_out_of_100": fmt(
            relative_fit["proposed"]
        ),
        "proposed_minus_original_relative_fit_score_points": fmt(
            relative_fit["gap"]
        ),
        "relative_fit_score_status": relative_fit["status"],
        "relative_fit_score_basis": RELATIVE_FIT_BASIS,
        "relative_fit_score_interpretation": RELATIVE_FIT_INTERPRETATION,
        "probability_interpretation": (
            "LEGACY_TOTAL_SITE_LIKELIHOOD_PERCENTAGE_RETAINED_FOR_AUDIT_"
            "SATURATED_NOT_HEADLINE_NOT_REAL_WORLD_ERROR_RATE"
            if not molecule_primary else
            "DIAGNOSTIC_SITE_PROBABILITY_EXCLUDED_BECAUSE_SCORER_QC_IS_MOLECULE_ORIENTED"
        ),
        "reconciliation_event_id": clean(proposed_row.get("reconciliation_event_id")),
        "reconciliation_event_class": clean(proposed_row.get("reconciliation_event_class")),
        "reconciliation_event_confidence": clean(proposed_row.get("reconciliation_event_confidence")),
        "reconciliation_final_action": clean(proposed_row.get("reconciliation_final_action")),
        "reconciliation_decision_confidence": clean(proposed_row.get("reconciliation_decision_confidence")),
        "reconciliation_reassignment_applied": clean(proposed_row.get("reconciliation_reassignment_applied")),
        "reconciliation_current_refined_assignment": clean(proposed_row.get("reconciliation_current_refined_assignment")),
        "source_identity": clean(proposed_row.get("source_identity")),
        "replaced_component": clean(proposed_row.get("replaced_component")),
        "replacement_component": clean(proposed_row.get("replacement_component")),
        "preserved_partner": clean(proposed_row.get("preserved_partner")),
        "original_project_genotype_status": clean(original_row.get("project_genotype_status")),
        "original_biological_admissibility": clean(original_row.get("biological_admissibility")),
        "proposed_project_genotype_status": clean(proposed_row.get("project_genotype_status")),
        "proposed_biological_admissibility": clean(proposed_row.get("biological_admissibility")),
        "snps_favor_original": original_sites,
        "snps_favor_proposed_swap": proposed_sites,
        "original_residual_mismatch": original_residual,
        "proposed_swap_residual_mismatch": proposed_residual,
        "similarity_explanation": similarity_explanation(probability),
        "ambient_original_c": fmt(left),
        "ambient_reconciled_c": fmt(right),
        "ambient_reconciled_minus_original_c": fmt(delta),
        "ambient_stratum": ambient_stratum,
        "ambient_assignment_diagnostic": ambient_diagnostic,
        "score_scope_contract": PAIR_CONTRACT,
        "aggregate_scope_contract": AGGREGATE_SCOPE_CONTRACT,
        "schema_version": SCHEMA,
    }
    for field in PAIR_EVIDENCE_FIELDS:
        row[field] = clean(probability.get(field))
    if not molecule_evidence_available(probability):
        for field in (
            "n_independent_molecules", "effective_independent_molecules",
            "maximum_molecule_contribution_fraction", "molecule_umi_gene_fraction",
            "probability_without_top_molecule_pct",
        ):
            row[field] = "NA"
    if eligible and available and not tied and relative_fit["status"] == "AVAILABLE":
        row["evidence_summary"] = (
            f"{proposed} relative fit {fmt(relative_fit['proposed'])}/100 versus "
            f"original {original} {fmt(relative_fit['original'])}/100; "
            f"{clean(probability.get('n_discriminating_snps'), '0')} "
            "discriminating SNPs; "
            f"{molecule_evidence_phrase(probability)}. The relative-fit score "
            "compares average fit between these two candidates and is not the "
            "probability that the identity is correct."
        )
    elif eligible and tied and relative_fit["status"] == "AVAILABLE":
        row["evidence_summary"] = (
            f"{original} and {proposed} have equal relative fit at "
            f"{fmt(relative_fit['original'])}/100 each under the site-based "
            "two-model nuclear comparison."
        )
    elif available and relative_fit["status"] == "AVAILABLE":
        row["evidence_summary"] = (
            "Diagnostic relative fit excluded from headline inference; "
            f"reason={eligibility_reason}; original={fmt(relative_fit['original'])}/100; "
            f"proposed={fmt(relative_fit['proposed'])}/100."
        )
    else:
        row["evidence_summary"] = (
            "No valid two-candidate nuclear probability; status="
            f"{clean(probability.get('comparison_status'))}; "
            f"{clean(probability.get('n_discriminating_snps'), '0')} "
            "discriminating SNPs."
        )
    return row


def event_summary(key, rows):
    library, event_id, original, proposed = key
    evaluable = [row for row in rows if probability_available(row)]
    proposed_wins = sum(row.get("proposed_swap_preferred") == "TRUE" for row in evaluable)
    original_wins = sum(row.get("proposed_swap_preferred") == "FALSE" for row in evaluable)
    ties = sum(row.get("proposed_swap_preferred") == "TIE" for row in evaluable)
    unavailable = len(rows) - len(evaluable)
    proposed_probabilities = [row.get("proposed_swap_probability_pct") for row in evaluable]
    original_probabilities = [row.get("original_assignment_probability_pct") for row in evaluable]
    relative_fit_rows = [row for row in evaluable if relative_fit_available(row)]
    proposed_relative_fit_scores = [
        row.get("proposed_swap_relative_fit_score_out_of_100")
        for row in relative_fit_rows
    ]
    original_relative_fit_scores = [
        row.get("original_assignment_relative_fit_score_out_of_100")
        for row in relative_fit_rows
    ]
    relative_fit_gaps = [
        row.get("proposed_minus_original_relative_fit_score_points")
        for row in relative_fit_rows
    ]
    fraction_proposed = proposed_wins / len(evaluable) if evaluable else math.nan
    if not evaluable:
        result = "NO_ELIGIBLE_NUCLEAR_PROBABILITIES"
    elif ties == len(evaluable):
        result = "ALL_ELIGIBLE_CELLS_TIED"
    elif proposed_wins > original_wins:
        result = "PROPOSED_SWAP_PREFERRED_IN_MAJORITY_OF_ELIGIBLE_CELLS"
    elif original_wins > proposed_wins:
        result = "ORIGINAL_ASSIGNMENT_PREFERRED_IN_MAJORITY_OF_ELIGIBLE_CELLS"
    else:
        result = "ORIGINAL_AND_PROPOSED_SPLIT_EVENLY_AMONG_ELIGIBLE_CELLS"
    ambient_deltas = [
        row.get("ambient_reconciled_minus_original_c") for row in evaluable
        if math.isfinite(fnum(row.get("ambient_reconciled_minus_original_c")))
    ]
    molecule_rows = [
        row for row in evaluable if molecule_evidence_available(row)
    ]
    molecule_available = len(molecule_rows)
    population_scope_counts = Counter(
        clean(row.get("score_population_scope")) for row in rows
    )
    reassignment_scopes = sorted(
        scope for scope in REASSIGNMENT_SCOPES
        if population_scope_counts.get(scope, 0)
    )
    eligible_reassignment_rows = [
        row for row in evaluable
        if clean(row.get("score_population_scope")) in REASSIGNMENT_SCOPES
    ]
    n_authoritative = (
        len(eligible_reassignment_rows) if len(reassignment_scopes) == 1 else 0
    )
    authoritative_scope = (
        reassignment_scopes[0] if len(reassignment_scopes) == 1
        else "MIXED_REASSIGNMENT_APPLICATION_STATE" if reassignment_scopes
        else "NONE"
    )
    if len(reassignment_scopes) > 1:
        result = "NO_AUTHORITATIVE_VERDICT_MIXED_APPLICATION_STATE"
    error_stable = sum(
        is_true(row.get("error_sensitivity_stable")) for row in evaluable
    )
    row = {
        "library": library,
        "reconciliation_event_id": event_id,
        "original_allowed_demux_assignment": original,
        "reconciliation_nominated_swap": proposed,
        "reconciliation_event_class": counter_text(
            row.get("reconciliation_event_class") for row in rows
        ),
        "reconciliation_event_confidence": counter_text(
            row.get("reconciliation_event_confidence") for row in rows
        ),
        "n_cells_scored": len(rows),
        "n_cells_in_authoritative_verdict": n_authoritative,
        "authoritative_reassignment_scope": authoritative_scope,
        "n_applied_reassignment_cells": population_scope_counts.get(
            APPLIED_REASSIGNMENT, 0
        ),
        "n_recommended_not_applied_cells": population_scope_counts.get(
            RECOMMENDED_NOT_APPLIED, 0
        ),
        "n_supported_reassignment_cells": sum(
            population_scope_counts.get(scope, 0)
            for scope in REASSIGNMENT_SCOPES
        ),
        "n_supported_event_held_cells": population_scope_counts.get(
            SUPPORTED_EVENT_HELD_CELL, 0
        ),
        "n_review_only_cells": population_scope_counts.get(
            REVIEW_ONLY_UNEXPECTED_IDENTITY, 0
        ),
        "score_population_scope_counts": ",".join(
            f"{scope}:{population_scope_counts[scope]}"
            for scope in sorted(population_scope_counts)
        ) or "NONE",
        "n_cells_probability_available": len(evaluable),
        "n_cells_proposed_swap_preferred": proposed_wins,
        "n_cells_original_preferred": original_wins,
        "n_cells_tied": ties,
        "n_cells_probability_unavailable": unavailable,
        "probability_exclusion_reason_counts": counter_text(
            row.get("headline_probability_exclusion_reason")
            for row in rows if not probability_available(row)
        ),
        "fraction_cells_proposed_swap_preferred": fmt(fraction_proposed),
        "median_proposed_swap_probability_pct": fmt(median(proposed_probabilities)),
        "proposed_swap_probability_q10_pct": fmt(quantile(proposed_probabilities, 0.10)),
        "proposed_swap_probability_q25_pct": fmt(quantile(proposed_probabilities, 0.25)),
        "proposed_swap_probability_q75_pct": fmt(quantile(proposed_probabilities, 0.75)),
        "proposed_swap_probability_q90_pct": fmt(quantile(proposed_probabilities, 0.90)),
        "median_original_probability_pct": fmt(median(original_probabilities)),
        "median_winning_probability_pct": fmt(median(
            row.get("winning_probability_pct") for row in evaluable
        )),
        "median_alternative_closeness_out_of_100": fmt(median(
            row.get("alternative_closeness_out_of_100") for row in evaluable
        )),
        "n_cells_relative_fit_available": len(relative_fit_rows),
        "median_proposed_swap_relative_fit_score_out_of_100": fmt(median(
            proposed_relative_fit_scores
        )),
        "proposed_swap_relative_fit_score_q10_out_of_100": fmt(quantile(
            proposed_relative_fit_scores, 0.10
        )),
        "proposed_swap_relative_fit_score_q25_out_of_100": fmt(quantile(
            proposed_relative_fit_scores, 0.25
        )),
        "proposed_swap_relative_fit_score_q75_out_of_100": fmt(quantile(
            proposed_relative_fit_scores, 0.75
        )),
        "proposed_swap_relative_fit_score_q90_out_of_100": fmt(quantile(
            proposed_relative_fit_scores, 0.90
        )),
        "median_original_assignment_relative_fit_score_out_of_100": fmt(median(
            original_relative_fit_scores
        )),
        "median_proposed_minus_original_relative_fit_score_points": fmt(median(
            relative_fit_gaps
        )),
        "median_discriminating_snps": fmt(median(
            row.get("n_discriminating_snps") for row in evaluable
        )),
        "median_effective_independent_snps": fmt(median(
            row.get("effective_independent_snps") for row in evaluable
        )),
        "median_independent_molecules": fmt(median(
            row.get("n_independent_molecules") for row in molecule_rows
        )),
        "fraction_molecule_evidence_available": fraction(
            molecule_available, len(evaluable)
        ),
        "median_maximum_site_contribution_fraction": fmt(median(
            row.get("maximum_site_contribution_fraction") for row in evaluable
        )),
        "median_top_five_site_contribution_fraction": fmt(median(
            row.get("top_five_site_contribution_fraction") for row in evaluable
        )),
        "site_top_site_removal_status_counts": counter_text(
            row.get("site_top_site_removal_status") for row in rows
        ),
        "site_top_five_removal_status_counts": counter_text(
            row.get("site_top_five_removal_status") for row in rows
        ),
        "median_downsample_50pct_win_fraction": fmt(median(
            row.get("downsample_50pct_win_fraction") for row in evaluable
        )),
        "median_probability_without_top_site_pct": fmt(median(
            row.get("probability_without_top_site_pct") for row in evaluable
        )),
        "median_probability_without_top_five_sites_pct": fmt(median(
            row.get("probability_without_top_five_sites_pct")
            for row in evaluable
        )),
        "median_minimum_error_sensitivity_probability_pct": fmt(median(
            row.get("minimum_error_sensitivity_probability_pct")
            for row in evaluable
        )),
        "fraction_error_sensitivity_stable": fraction(
            error_stable, len(evaluable)
        ),
        "absolute_fit_status_counts": counter_text(
            row.get("absolute_fit_status") for row in rows
        ),
        "warning_counts": warning_counter_text(rows),
        "final_action_counts": counter_text(
            row.get("reconciliation_final_action") for row in rows
        ),
        "reassignment_applied_counts": counter_text(
            row.get("reconciliation_reassignment_applied") for row in rows
        ),
        "n_cells_with_controlled_ambient": len(ambient_deltas),
        "median_ambient_change_c": fmt(median(ambient_deltas)),
        "event_result": result,
        "score_scope_contract": PAIR_CONTRACT,
        "aggregate_scope_contract": AGGREGATE_SCOPE_CONTRACT,
        "schema_version": SCHEMA,
    }
    molecule_summary = (
        f"median independent molecules {row['median_independent_molecules']}"
        if molecule_rows else "independent molecules not measured"
    )
    if len(reassignment_scopes) > 1:
        row["event_evidence_summary"] = (
            "Applied and recommended-not-applied cells coexist in this scope; "
            "zero cells entered an authoritative combined verdict."
        )
    else:
        row["event_evidence_summary"] = (
            f"{proposed} versus original {original}: {proposed_wins}/{len(evaluable)} "
            f"eligible cells prefer the proposed swap, {original_wins} prefer the "
            f"original, and {ties} are tied; median relative fit proposed "
            f"{row['median_proposed_swap_relative_fit_score_out_of_100']}/100 "
            f"versus original "
            f"{row['median_original_assignment_relative_fit_score_out_of_100']}/100; "
            f"median discriminating SNPs {row['median_discriminating_snps']}; "
            f"{molecule_summary}. Relative fit compares average fit between the "
            "two fixed candidates and is not identity correctness probability."
        )
    return row


def consolidated_event_summary(key, rows):
    library, event_id, proposed = key
    originals = sorted({
        row["original_allowed_demux_assignment"] for row in rows
    })
    decision_rows = [
        row for row in rows
        if row.get("score_population_scope") in REASSIGNMENT_SCOPES
    ]
    if not decision_rows:
        raise ValueError(
            f"lib{library} {event_id}: supported event has no reassignment cells"
        )
    summary = event_summary(
        (library, event_id, ",".join(originals), proposed), decision_rows
    )
    all_scope_counts = Counter(
        clean(row.get("score_population_scope")) for row in rows
    )
    summary["n_cells_scored"] = len(rows)
    summary["n_applied_reassignment_cells"] = all_scope_counts.get(
        APPLIED_REASSIGNMENT, 0
    )
    summary["n_recommended_not_applied_cells"] = all_scope_counts.get(
        RECOMMENDED_NOT_APPLIED, 0
    )
    summary["n_supported_reassignment_cells"] = len(decision_rows)
    summary["n_supported_event_held_cells"] = all_scope_counts.get(
        SUPPORTED_EVENT_HELD_CELL, 0
    )
    summary["score_population_scope_counts"] = ",".join(
        f"{scope}:{all_scope_counts[scope]}" for scope in sorted(all_scope_counts)
    ) or "NONE"
    decision_scope_set = {
        row.get("score_population_scope") for row in decision_rows
    }
    summary["n_cells_in_authoritative_verdict"] = (
        sum(probability_available(row) for row in decision_rows)
        if len(decision_scope_set) == 1 else 0
    )
    if len(decision_scope_set) != 1:
        summary["authoritative_reassignment_scope"] = (
            "MIXED_REASSIGNMENT_APPLICATION_STATE"
        )
        summary["event_result"] = (
            "NO_AUTHORITATIVE_VERDICT_MIXED_APPLICATION_STATE"
        )
        summary["event_evidence_summary"] = (
            "Applied and recommended-not-applied cells coexist in one event; "
            "no combined verdict was calculated."
        )
    else:
        summary["event_evidence_summary"] += (
            f" {summary['n_supported_event_held_cells']} attached held cells "
            "were excluded from the reassignment verdict."
        )
    summary["original_assignments"] = ",".join(originals)
    summary["n_original_assignment_strata"] = len(originals)
    return summary


def consolidated_review_event_summary(key, rows):
    library, event_id, proposed = key
    originals = sorted({
        row["original_allowed_demux_assignment"] for row in rows
    })
    summary = event_summary(
        (library, event_id, ",".join(originals), proposed), rows
    )
    summary["original_assignments"] = ",".join(originals)
    summary["n_original_assignment_strata"] = len(originals)
    summary["n_cells_in_authoritative_verdict"] = 0
    summary["authoritative_reassignment_scope"] = "NONE"
    summary["event_result"] = "REVIEW_ONLY_NO_REASSIGNMENT_VERDICT"
    summary["event_evidence_summary"] = (
        "Unexpected-identity review probabilities are diagnostic only; no "
        "cells entered a reassignment verdict."
    )
    return summary


def build_distribution_rows(cell_rows):
    counts = Counter()
    for row in cell_rows:
        probability = fnum(row.get("proposed_swap_probability_pct"))
        if not probability_available(row):
            continue
        lower = min(95, max(0, int(probability // 5) * 5))
        key = (
            row["library"], row["reconciliation_event_id"],
            row["original_allowed_demux_assignment"],
            row["reconciliation_nominated_swap"],
            row["score_population_scope"], lower,
        )
        counts[key] += 1
    return [{
        "library": key[0],
        "reconciliation_event_id": key[1],
        "original_allowed_demux_assignment": key[2],
        "reconciliation_nominated_swap": key[3],
        "score_population_scope": key[4],
        "proposed_probability_bin_start_pct": key[5],
        "proposed_probability_bin_end_pct": key[5] + 5,
        "n_cells": count,
        "schema_version": SCHEMA,
    } for key, count in sorted(counts.items())]


def build_cross_library_rows(cell_rows, event_rows, expected_libraries):
    cells_by_proposed = defaultdict(list)
    events_by_proposed = defaultdict(list)
    for row in cell_rows:
        if row.get("score_population_scope") not in REASSIGNMENT_SCOPES:
            continue
        cells_by_proposed[row["reconciliation_nominated_swap"]].append(row)
    for row in event_rows:
        events_by_proposed[row["reconciliation_nominated_swap"]].append(row)
    output = []
    for proposed, rows in sorted(cells_by_proposed.items()):
        events = events_by_proposed[proposed]
        nominated_libraries = sorted({f"lib{int(row['library'])}" for row in rows})
        expected = sorted(expected_libraries.get(proposed, set()))
        other_expected = sorted(set(expected) - set(nominated_libraries))
        evaluable = [row for row in rows if probability_available(row)]
        relative_fit_rows = [
            row for row in evaluable if relative_fit_available(row)
        ]
        proposed_wins = sum(
            row.get("proposed_swap_preferred") == "TRUE" for row in evaluable
        )
        if len(nominated_libraries) > 1:
            status = "RECONCILIATION_NOMINATION_RECURS_ACROSS_LIBRARIES"
        elif other_expected:
            status = "NOMINATED_ONCE_WITH_EXPECTED_REPLICATE_LIBRARY_CONTEXT"
        else:
            status = "SINGLE_LIBRARY_NOMINATION"
        output.append({
            "reconciliation_nominated_swap": proposed,
            "n_libraries_nominated": len(nominated_libraries),
            "libraries_nominated": ",".join(nominated_libraries),
            "n_reconciliation_events": len(events),
            "event_ids": ",".join(sorted({row["reconciliation_event_id"] for row in events})),
            "original_assignments": ",".join(sorted({row["original_allowed_demux_assignment"] for row in rows})),
            "n_cells_scored": len(rows),
            "n_cells_probability_available": len(evaluable),
            "fraction_cells_proposed_swap_preferred": fraction(proposed_wins, len(evaluable)),
            "median_proposed_swap_probability_pct": fmt(median(
                row.get("proposed_swap_probability_pct") for row in evaluable
            )),
            "median_proposed_swap_relative_fit_score_out_of_100": fmt(median(
                row.get("proposed_swap_relative_fit_score_out_of_100")
                for row in relative_fit_rows
            )),
            "median_original_assignment_relative_fit_score_out_of_100": fmt(median(
                row.get("original_assignment_relative_fit_score_out_of_100")
                for row in relative_fit_rows
            )),
            "event_class_counts": counter_text(
                row.get("reconciliation_event_class") for row in rows
            ),
            "libraries_where_identity_is_expected": ",".join(expected) or "NONE",
            "other_expected_replicate_libraries": ",".join(other_expected) or "NONE",
            "cross_library_status": status,
            "score_scope_contract": PAIR_CONTRACT,
            "aggregate_scope_contract": AGGREGATE_SCOPE_CONTRACT,
            "schema_version": SCHEMA,
        })
    return output


def write_plot_handoff(plot_root, run_summary):
    handoff = f"""# Supported reconciliation-event score plotting handoff

## Authoritative population boundary

The headline plotting population is already finalized in these tables.  It
contains only:

1. cells whose final v10 action is `REASSIGN_GENOTYPE` within
   `LIKELY_UNEXPECTED_INTACT_BIOLOGICAL_LINE` or
   `LIKELY_UNEXPECTED_SINGLET_POPULATION`; and
2. `REVIEW_CELLULAR_ORIGIN` cells attached to the same exact
   `(library, reconciliation_event_id, proposed identity)` as at least one
   `REASSIGN_GENOTYPE` cell.

The populations remain distinct through `score_population_scope`:
`APPLIED_REASSIGNMENT`, `RECOMMENDED_NOT_APPLIED`, and
`SUPPORTED_EVENT_HELD_CELL`.  The authoritative event verdict uses only one
reassignment scope.  Held cells never vote.  If applied and recommended cells
coexist in one event, no combined verdict is emitted.

`REVIEW_UNEXPECTED_IDENTITY` rows are retained only in the three files whose
names begin `review_only_`.  They are not part of any headline swap result.
Never plot `KEEP`, conflicted, unresolved, occupancy/ploidy, standalone cellular
origin review, below-event-mass cell-level cleanup, expected-roster,
technical-multiplet, or invented identities.

## Initial library 19 evaluation: tables only

Do not generate or interpret plots during the initial frozen library 19
evaluation.  Inspect the cell, event, scope, exclusion, and provenance tables
first.  If plotting is authorized after that review, keep the report to three
operational views:

1. a ranked supported-event table, with applied, recommended, and attached-held counts;
2. per-event original-versus-proposed relative-fit distributions derived from
   the cell-table relative-fit fields, visibly split into applied,
   recommended-not-applied, and held scopes; and
3. evidence/stability plus controlled ambient and replicate context for those
   same supported events.

Do not add all-cell PCA, generic identity landscapes, open/closed/bridge plots,
unconstrained transitions, or generic candidate discovery.  Review-only rows
may appear on one explicitly optional page and must never affect headline ranks.

## Authoritative files

- `cells_scored_for_reconciliation_swaps.tsv.gz`: headline supported cells only.
- `swap_events_within_library.tsv`: one consolidated row per supported event.
- `swap_event_original_assignment_strata.tsv`: optional original-assignment detail.
- `swap_event_scope_summary.tsv`: independent applied/recommended/held summaries.
- `swap_candidates_across_libraries.tsv`: supported proposed identities only.
- `proposed_swap_probability_distributions.tsv`: legacy saturated total-likelihood
  percentage bins retained for audit compatibility; never use as headline score.
- `review_only_*`: separate unexpected-identity review outputs.
- `score_population_scope_summary.tsv`: audit from source scores to final scopes.
- `score_pair_exclusions.tsv`: exclusions from pair construction and aggregation.
- `ambient_donor_burden.tsv`: optional controlled ambient evidence.
- `identity_score_provenance.tsv`: hashes linking validation, decisions, pairs,
  probabilities, metadata, scorer, and aggregator.

The plot-data manifest is authoritative.  Ignore stale unlisted files that may
exist in the directory from an earlier aggregate.

## Relative-fit and legacy-probability language

Use `proposed_swap_relative_fit_score_out_of_100` and
`original_assignment_relative_fit_score_out_of_100` as the experimental
pairwise-fit display.  They sum to 100 and compare average site fit after
normalizing the site log-likelihood difference by discriminating evidence
depth.  They are not probabilities that either biological identity is correct.
The legacy `*_probability_pct` fields use total accumulated likelihood, are
saturated in Library 19, and are retained only for audit compatibility.  Do
not use them as headline scores.  Missing molecule evidence means “not
measured”; never display it as zero or substitute SNP counts.

For the frozen site-only library 19 arm, eligibility uses the site likelihood,
PASS comparison/fit status, and explicit top-site and top-five-site removal
checks.  Rows for which the v3 scorer selected molecule-primary QC are
unevaluable because their fit and stability fields are not aligned to the site
winner.  Neither site nor molecule support is declared universally
authoritative.  Ties are neutral.

The raw combined influence flag and leave-one-chromosome value remain only for
source provenance.  Never use either in eligibility, event summaries, or
plots: the combined flag includes the invalid approximately 4,000-sequence
leave-one-contig diagnostic.  Top-contig, top-five-contig, 1/5/10% evidence
removal, and evidence-balanced genomic-block diagnostics are not present in
the frozen scorer output and are explicitly deferred.

- Run schema: `{SCHEMA}`
- Aggregate scope: `{AGGREGATE_SCOPE_CONTRACT}`
- Headline supported cells: {run_summary['n_cells_scored']}
- Applied reassignments: {run_summary['n_applied_reassignment_cells']}
- Recommended, not applied: {run_summary['n_recommended_not_applied_cells']}
- Attached held cells: {run_summary['n_supported_event_held_cells']}
- Separate review-only cells: {run_summary['n_review_only_cells']}
- Supported within-library events: {run_summary['n_library_swap_events']}
"""
    (plot_root / "plotting_handoff.md").write_text(handoff)
    dictionary = f"""# Supported reconciliation-event score data dictionary

The experimental display fields are
`original_assignment_relative_fit_score_out_of_100` and
`proposed_swap_relative_fit_score_out_of_100`.  They sum to 100 when available.
They are calculated by dividing the oriented site log-likelihood difference by
`discriminating_evidence_depth` and applying a two-candidate logistic
transformation.  They measure average comparative fit between the two fixed
candidates; they are neither accumulated decision confidence nor calibrated
identity-correctness probabilities.

The legacy `original_assignment_probability_pct` and
`proposed_swap_probability_pct` fields also sum to 100 but are based on total
accumulated site likelihood and saturate in Library 19.  They and
`alternative_closeness_out_of_100` are retained only for provenance and
backward-compatible diagnostics, not as headline scores.

Eligibility is limited to scorer rows whose primary QC orientation is also
site based.  `site_top_site_removal_status` and
`site_top_five_removal_status` are the only influence-removal gates.  The raw
combined influence and leave-one-chromosome fields are non-authoritative source
provenance.  Ambient values are separate controlled annotations and do not
modify nuclear support.  `score_scope_contract` must always equal
`{PAIR_CONTRACT}`.  `aggregate_scope_contract` must always equal
`{AGGREGATE_SCOPE_CONTRACT}`.

`score_population_scope` has four possible retained values:

- `{APPLIED_REASSIGNMENT}`: the reconciliation actually rewrote the assignment;
- `{RECOMMENDED_NOT_APPLIED}`: decisive recommendation recorded without rewrite;
- `{SUPPORTED_EVENT_HELD_CELL}`: occupancy-held cell attached to an exact
  applied or recommended event; it never votes in that event's verdict;
- `{REVIEW_ONLY_UNEXPECTED_IDENTITY}`: separate review-only output, never headline.
"""
    (plot_root / "data_dictionary.md").write_text(dictionary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--identity-root", required=True)
    parser.add_argument("--libraries", nargs="+", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--ambient-root", default=None)
    parser.add_argument("--ambient-contrast-cells", default=None)
    parser.add_argument("--ambient-condition", default=None)
    parser.add_argument(
        "--reuse-frozen-probabilities", action="store_true",
        help=("Aggregate already-scored v1/v2 pair manifests and probability "
              "outputs without regenerating them; missing legacy scorer "
              "provenance is reported, not fabricated"),
    )
    parser.add_argument(
        "--scorer-binary", default=None,
        help=("Current tetra_score_calls path used only to record a warning-only "
              "software-hash comparison against scored-run provenance"),
    )
    args = parser.parse_args()

    libraries = parse_libraries(args.libraries)
    identity_root = Path(args.identity_root).resolve()
    output_root = Path(args.output_root).resolve()
    plot_root = output_root / "plot_data"
    output_root.mkdir(parents=True, exist_ok=True)
    output_lock = acquire_output_lock(output_root)
    plot_root.mkdir(parents=True, exist_ok=True)

    validation_path, validation_sha256 = validate_reconciliation_summary(
        identity_root
    )
    aggregator_sha256 = sha256_file(Path(__file__).resolve())
    scorer_binary = (
        str(Path(args.scorer_binary).resolve())
        if args.scorer_binary and Path(args.scorer_binary).is_file()
        else "NA"
    )
    current_scorer_sha256 = (
        sha256_file(scorer_binary) if scorer_binary != "NA" else "NA"
    )

    ambient_values, ambient_burden_rows, ambient_provenance = load_ambient(
        args, libraries
    )
    expected_libraries = load_expected_identity_libraries(identity_root)
    global_lines, global_donors = load_project_biological_universe(
        identity_root
    )
    cell_rows = []
    review_only_cell_rows = []
    manifest_summary_rows = []
    provenance_rows = []
    exclusion_counts = Counter()
    scope_counts_by_library = defaultdict(Counter)
    supported_event_keys_by_library = {}
    review_event_keys_by_library = defaultdict(set)
    libraries_with_source_pairs = []
    libraries_with_pairs = []
    libraries_with_review_only = []

    for library in libraries:
        pair_root = identity_root / "score_pairs"
        manifest_path = pair_root / f"lib{library}.reconciliation_score_pairs.tsv.gz"
        summary_path = pair_root / f"lib{library}.reconciliation_score_pair_summary.tsv"
        probability_path = (
            identity_root / "nuclear" /
            f"lib{library}.identity_pair_probabilities.tsv.gz"
        )
        if not manifest_path.is_file() or not summary_path.is_file():
            raise FileNotFoundError(
                f"lib{library}: targeted score-pair manifest/summary missing"
            )
        manifest = load_pair_manifest(
            manifest_path, library, expected_libraries,
            global_lines, global_donors,
        )
        summary_header, summaries = read_table(summary_path)
        require_columns(
            summary_header,
            {"library", "n_score_pairs", "n_manifest_rows",
             "n_candidate_rows_excluded", "exclusion_reason_counts",
             "source_decisions", "source_decisions_sha256",
             "score_scope_contract", "schema_version"},
            summary_path,
        )
        if len(summaries) != 1:
            raise ValueError(f"{summary_path}: expected exactly one summary row")
        summary = summaries[0]
        summary_schema = clean(summary.get("schema_version"))
        if summary_schema not in ACCEPTED_PAIR_MANIFEST_SCHEMAS:
            raise ValueError(
                f"{summary_path}: unsupported pair-summary schema {summary_schema!r}"
            )
        if clean(summary.get("score_scope_contract")) != PAIR_CONTRACT:
            raise ValueError(f"{summary_path}: score-scope contract mismatch")
        if fint(summary.get("n_score_pairs"), -1) != len(manifest):
            raise ValueError(f"{summary_path}: pair count does not match manifest")
        if fint(summary.get("n_manifest_rows"), -1) != 2 * len(manifest):
            raise ValueError(f"{summary_path}: manifest row count is not two per pair")
        decisions_path = (
            identity_root / "decisions" / f"lib{library}.reconciled_cells.tsv.gz"
        )
        if not decisions_path.is_file():
            raise FileNotFoundError(
                f"lib{library}: current decision table is missing: {decisions_path}"
            )
        decision_sha256 = sha256_file(decisions_path)
        if clean(summary.get("source_decisions_sha256"), "") != decision_sha256:
            raise ValueError(
                f"{summary_path}: score pairs are stale relative to {decisions_path}"
            )
        if (
            summary_schema == PAIR_MANIFEST_SCHEMA
            and validation_path.stat().st_mtime_ns < decisions_path.stat().st_mtime_ns
        ):
            raise ValueError(
                f"{validation_path}: validation predates {decisions_path}"
            )
        if (
            summary_schema == PAIR_MANIFEST_SCHEMA
            and clean(summary.get("validation_summary_sha256"), "")
            != validation_sha256
        ):
            raise ValueError(
                f"{summary_path}: score pairs do not reference the current "
                "validation summary"
            )
        metadata_paths = {
            "metadata_expected_genotypes_sha256": (
                identity_root / "metadata" / "library_expected_genotypes.tsv"
            ),
            "metadata_global_biological_lines_sha256": (
                identity_root / "metadata" / "global_biological_lines.tsv"
            ),
            "metadata_global_donors_sha256": (
                identity_root / "metadata" / "global_donors.tsv"
            ),
        }
        metadata_hashes = {
            field: sha256_file(path) for field, path in metadata_paths.items()
        }
        if summary_schema == PAIR_MANIFEST_SCHEMA:
            for field, observed in metadata_hashes.items():
                if clean(summary.get(field), "") != observed:
                    raise ValueError(
                        f"{summary_path}: {field} does not match current metadata"
                    )
        manifest_summary_rows.append(summary)
        reason_text = clean(summary.get("exclusion_reason_counts"), "NONE")
        if reason_text != "NONE":
            for item in reason_text.split(","):
                reason, count = item.rsplit(":", 1)
                exclusion_counts[(library, "PAIR_BUILD", reason)] += int(count)

        if not manifest:
            if probability_path.is_file():
                exclusion_counts[(
                    library, "AGGREGATE_SCOPE",
                    "STALE_PROBABILITY_FILE_IGNORED_FOR_ZERO_PAIR_MANIFEST",
                )] += 1
            provenance_rows.append({
                "library": library,
                "validation_summary": str(validation_path),
                "validation_summary_sha256": validation_sha256,
                "decision_table": str(decisions_path),
                "decision_table_sha256": decision_sha256,
                "pair_manifest": str(manifest_path),
                "pair_manifest_sha256": sha256_file(manifest_path),
                "pair_summary": str(summary_path),
                "pair_summary_sha256": sha256_file(summary_path),
                "probability_output": "NA", "probability_output_sha256": "NA",
                "probability_provenance": "NA", "probability_provenance_sha256": "NA",
                "pair_manifest_schema": summary_schema,
                "pair_builder_sha256": clean(summary.get("score_pair_builder_sha256")),
                **metadata_hashes,
                "aggregator_sha256": aggregator_sha256,
                "scorer_binary": scorer_binary,
                "scorer_binary_sha256": current_scorer_sha256,
                "software_hash_status": "NOT_APPLICABLE_ZERO_PAIRS",
                "source_pair_scores_reused": str(args.reuse_frozen_probabilities).upper(),
                "schema_version": SCHEMA,
            })
            continue
        libraries_with_source_pairs.append(library)
        if not probability_path.is_file():
            raise FileNotFoundError(
                f"lib{library}: targeted probability output missing: {probability_path}"
            )
        manifest_schemas = {
            pair["manifest_schema"] for pair in manifest.values()
        }
        if manifest_schemas != {summary_schema}:
            raise ValueError(
                f"{manifest_path}: manifest/summary schema mismatch "
                f"({sorted(manifest_schemas)} versus {summary_schema})"
            )
        probability_provenance_path = (
            identity_root / "nuclear" /
            f"lib{library}.identity_pair_probability_provenance.tsv"
        )
        probability_provenance, software_hash_status, probability_provenance_sha = (
            validate_probability_provenance(
                probability_provenance_path, manifest_path, summary_path,
                probability_path,
                None if scorer_binary == "NA" else scorer_binary,
                allow_legacy_missing=(
                    args.reuse_frozen_probabilities
                    and summary_schema == LEGACY_PAIR_MANIFEST_SCHEMA
                ),
            )
        )
        provenance_rows.append({
            "library": library,
            "validation_summary": str(validation_path),
            "validation_summary_sha256": validation_sha256,
            "decision_table": str(decisions_path),
            "decision_table_sha256": decision_sha256,
            "pair_manifest": str(manifest_path),
            "pair_manifest_sha256": sha256_file(manifest_path),
            "pair_summary": str(summary_path),
            "pair_summary_sha256": sha256_file(summary_path),
            "probability_output": str(probability_path),
            "probability_output_sha256": sha256_file(probability_path),
            "probability_provenance": (
                str(probability_provenance_path)
                if probability_provenance_path.is_file() else "NA"
            ),
            "probability_provenance_sha256": probability_provenance_sha,
            "pair_manifest_schema": summary_schema,
            "pair_builder_sha256": clean(summary.get("score_pair_builder_sha256")),
            **metadata_hashes,
            "aggregator_sha256": aggregator_sha256,
            "scorer_binary": scorer_binary,
            "scorer_binary_sha256": current_scorer_sha256,
            "software_hash_status": software_hash_status,
            "source_pair_scores_reused": str(args.reuse_frozen_probabilities).upper(),
            "schema_version": SCHEMA,
        })
        probability_header, probabilities = read_table(probability_path)
        require_columns(
            probability_header,
            {
                "library", "barcode", "score_pair_id", "comparison",
                "comparison_status", "candidate_a", "candidate_b",
                "candidate_a_role", "candidate_b_role",
                "preferred_assignment", "preferred_probability_pct",
                "alternative_assignment", "alternative_probability_pct",
                "candidate_a_probability_pct", "candidate_b_probability_pct",
                "schema_version", *PAIR_EVIDENCE_FIELDS,
            },
            probability_path,
        )
        by_barcode = {}
        for probability in probabilities:
            barcode = clean(probability.get("barcode"))
            if barcode in by_barcode:
                raise ValueError(f"{probability_path}: duplicate barcode {barcode}")
            if clean(probability.get("comparison")) != PAIR_COMPARISON:
                raise ValueError(
                    f"{probability_path}: forbidden comparison "
                    f"{clean(probability.get('comparison'))!r}"
                )
            if clean(probability.get("schema_version")) != PAIR_SCHEMA:
                raise ValueError(f"{probability_path}: unsupported probability schema")
            by_barcode[barcode] = probability
        if set(by_barcode) != set(manifest):
            missing = sorted(set(manifest) - set(by_barcode))[:5]
            extra = sorted(set(by_barcode) - set(manifest))[:5]
            raise ValueError(
                f"{probability_path}: probability/manifest barcode mismatch; "
                f"missing={missing} extra={extra}"
            )

        supported_event_keys = supported_event_keys_for_manifest(manifest)
        supported_event_keys_by_library[library] = supported_event_keys
        for barcode in sorted(manifest):
            pair = manifest[barcode]
            probability = by_barcode[barcode]
            if clean(probability.get("score_pair_id")) != pair["score_pair_id"]:
                raise ValueError(f"lib{library} {barcode}: score_pair_id mismatch")
            if clean(probability.get("candidate_a_role")) != ORIGINAL_ROLE:
                raise ValueError(f"lib{library} {barcode}: candidate A role is not original")
            if clean(probability.get("candidate_b_role")) != PROPOSED_ROLE:
                raise ValueError(f"lib{library} {barcode}: candidate B role is not proposed swap")
            if canonical_genotype(probability.get("candidate_a")) != pair["original_identity"]:
                raise ValueError(f"lib{library} {barcode}: candidate A identity mismatch")
            if canonical_genotype(probability.get("candidate_b")) != pair["proposed_identity"]:
                raise ValueError(f"lib{library} {barcode}: candidate B identity mismatch")
            p_original = fnum(probability.get("candidate_a_probability_pct"))
            p_proposed = fnum(probability.get("candidate_b_probability_pct"))
            if math.isfinite(p_original) != math.isfinite(p_proposed):
                raise ValueError(f"lib{library} {barcode}: only one candidate probability is finite")
            if math.isfinite(p_original):
                probability_sum = math.fsum((p_original, p_proposed))
                if (
                    p_original < 0.0 or p_original > 100.0
                    or p_proposed < 0.0 or p_proposed > 100.0
                    or abs(probability_sum - 100.0) > PROBABILITY_SUM_TOLERANCE_PP
                ):
                    raise ValueError(
                        f"lib{library} {barcode}: candidate probabilities do not "
                        f"sum to 100 ({probability_sum})"
                    )
            population_scope = classify_pair_population_scope(
                pair, supported_event_keys
            )
            event_key = pair_event_key(pair)
            supported_event_key = (
                pair_event_key_text(library, event_key)
                if population_scope in HEADLINE_POPULATION_SCOPES else "NA"
            )
            validate_optional_scope_annotations(
                manifest_path, library, pair, population_scope,
                supported_event_key,
            )
            cell_row = build_cell_row(
                library, barcode, pair, probability,
                ambient_values.get((library, barcode)),
                population_scope or "EXCLUDED_FROM_SUPPORTED_EVENT_SCOPE",
                supported_event_key,
            )
            if population_scope in HEADLINE_POPULATION_SCOPES:
                cell_rows.append(cell_row)
                scope_counts_by_library[library][population_scope] += 1
            elif population_scope == REVIEW_ONLY_UNEXPECTED_IDENTITY:
                review_only_cell_rows.append(cell_row)
                review_event_keys_by_library[library].add(event_key)
                scope_counts_by_library[library][population_scope] += 1
            else:
                reason = aggregate_scope_exclusion_reason(
                    pair, supported_event_keys
                )
                exclusion_counts[(library, "AGGREGATE_SCOPE", reason)] += 1
                scope_counts_by_library[library]["EXCLUDED"] += 1

        if any(
            scope_counts_by_library[library].get(scope, 0)
            for scope in HEADLINE_POPULATION_SCOPES
        ):
            libraries_with_pairs.append(library)
        if scope_counts_by_library[library].get(
                REVIEW_ONLY_UNEXPECTED_IDENTITY, 0):
            libraries_with_review_only.append(library)

    event_groups = defaultdict(list)
    event_stratum_groups = defaultdict(list)
    event_scope_groups = defaultdict(list)
    for row in cell_rows:
        event_key = (
            int(row["library"]), row["reconciliation_event_id"],
            row["reconciliation_nominated_swap"],
        )
        stratum_key = (
            int(row["library"]), row["reconciliation_event_id"],
            row["original_allowed_demux_assignment"],
            row["reconciliation_nominated_swap"],
        )
        event_groups[event_key].append(row)
        scope_key = (*stratum_key, row["score_population_scope"])
        event_scope_groups[scope_key].append(row)
        if row["score_population_scope"] in REASSIGNMENT_SCOPES:
            event_stratum_groups[stratum_key].append(row)
    event_rows = [
        consolidated_event_summary(key, rows)
        for key, rows in sorted(event_groups.items())
    ]
    event_stratum_rows = [
        event_summary(key, rows)
        for key, rows in sorted(event_stratum_groups.items())
    ]
    event_scope_rows = []
    for key, rows in sorted(event_scope_groups.items()):
        library, event_id, original, proposed, scope = key
        row = event_summary(
            (library, event_id, original, proposed), rows
        )
        row["score_population_scope"] = scope
        if scope == SUPPORTED_EVENT_HELD_CELL:
            row["n_cells_in_authoritative_verdict"] = 0
            row["authoritative_reassignment_scope"] = "NONE"
            row["event_result"] = "HELD_CELL_DIAGNOSTIC_ONLY_NO_REASSIGNMENT_VERDICT"
            row["event_evidence_summary"] = (
                "Attached held-cell probabilities are diagnostic only; no "
                "cells entered a reassignment verdict."
            )
        event_scope_rows.append(row)
    review_event_groups = defaultdict(list)
    for row in review_only_cell_rows:
        key = (
            int(row["library"]), row["reconciliation_event_id"],
            row["reconciliation_nominated_swap"],
        )
        review_event_groups[key].append(row)
    review_event_rows = [
        consolidated_review_event_summary(key, rows)
        for key, rows in sorted(review_event_groups.items())
    ]
    expected_supported_events = sum(
        len(keys) for keys in supported_event_keys_by_library.values()
    )
    if len(event_rows) != expected_supported_events:
        raise ValueError(
            "supported-event accounting mismatch: "
            f"event summaries={len(event_rows)} "
            f"supported keys={expected_supported_events}"
        )
    cross_library_rows = build_cross_library_rows(
        cell_rows, event_rows, expected_libraries
    )
    distribution_rows = build_distribution_rows(cell_rows)
    review_distribution_rows = build_distribution_rows(review_only_cell_rows)
    exclusion_rows = [{
        "library": library,
        "exclusion_stage": stage,
        "exclusion_reason": reason,
        "n_cells": count,
        "schema_version": SCHEMA,
    } for (library, stage, reason), count in sorted(exclusion_counts.items())]
    scope_summary_rows = []
    for library in libraries:
        counts = scope_counts_by_library[library]
        supported_count = sum(
            counts.get(scope, 0) for scope in HEADLINE_POPULATION_SCOPES
        )
        scope_summary_rows.append({
            "library": library,
            "n_source_pair_scores": sum(counts.values()),
            "n_supported_event_keys": len(
                supported_event_keys_by_library.get(library, set())
            ),
            "n_supported_reassignment_cells": counts.get(
                APPLIED_REASSIGNMENT, 0
            ) + counts.get(
                RECOMMENDED_NOT_APPLIED, 0
            ),
            "n_applied_reassignment_cells": counts.get(
                APPLIED_REASSIGNMENT, 0
            ),
            "n_recommended_not_applied_cells": counts.get(
                RECOMMENDED_NOT_APPLIED, 0
            ),
            "n_supported_event_held_cells": counts.get(
                SUPPORTED_EVENT_HELD_CELL, 0
            ),
            "n_headline_supported_cells": supported_count,
            "n_review_only_unexpected_identity_cells": counts.get(
                REVIEW_ONLY_UNEXPECTED_IDENTITY, 0
            ),
            "n_excluded_after_scoring": counts.get("EXCLUDED", 0),
            "n_headline_supported_events": len(
                supported_event_keys_by_library.get(library, set())
            ),
            "n_review_only_events": len(
                review_event_keys_by_library.get(library, set())
            ),
            "aggregate_scope_contract": AGGREGATE_SCOPE_CONTRACT,
            "schema_version": SCHEMA,
        })

    write_table(
        output_root / "cell_swap_identity_scores.tsv.gz",
        cell_rows, CELL_FIELDS,
    )
    write_table(
        output_root / "library_swap_event_summary.tsv",
        event_rows, EVENT_FIELDS,
    )
    write_table(
        output_root / "library_swap_event_strata.tsv",
        event_stratum_rows, EVENT_STRATUM_FIELDS,
    )
    write_table(
        output_root / "library_swap_event_scope_summary.tsv",
        event_scope_rows, EVENT_SCOPE_FIELDS,
    )
    write_table(
        output_root / "cross_library_swap_summary.tsv",
        cross_library_rows, CROSS_LIBRARY_FIELDS,
    )
    write_table(
        output_root / "proposed_swap_probability_distributions.tsv",
        distribution_rows, DISTRIBUTION_FIELDS,
    )
    write_table(
        output_root / "review_only_identity_scores.tsv.gz",
        review_only_cell_rows, CELL_FIELDS,
    )
    write_table(
        output_root / "review_only_event_summary.tsv",
        review_event_rows, EVENT_FIELDS,
    )
    write_table(
        output_root / "review_only_probability_distributions.tsv",
        review_distribution_rows, DISTRIBUTION_FIELDS,
    )
    write_table(
        output_root / "score_population_scope_summary.tsv",
        scope_summary_rows, SCOPE_SUMMARY_FIELDS,
    )
    write_table(
        output_root / "score_pair_exclusion_summary.tsv",
        exclusion_rows,
        ["library", "exclusion_stage", "exclusion_reason", "n_cells",
         "schema_version"],
    )
    write_table(
        output_root / "identity_score_provenance.tsv",
        provenance_rows, PROVENANCE_FIELDS,
    )
    ambient_fields = (
        list(ambient_burden_rows[0]) if ambient_burden_rows else
        ["library", "condition", "schema_version"]
    )
    write_table(
        output_root / "ambient_donor_burden.tsv",
        ambient_burden_rows, ambient_fields,
    )

    run_summary = {
        "schema_version": SCHEMA,
        "identity_root": str(identity_root),
        "output_root": str(output_root),
        "libraries_requested": ",".join(str(value) for value in libraries),
        "libraries_with_source_pair_scores": (
            ",".join(str(value) for value in libraries_with_source_pairs)
            or "NONE"
        ),
        "libraries_with_scored_pairs": (
            ",".join(str(value) for value in libraries_with_pairs) or "NONE"
        ),
        "libraries_with_review_only_candidates": (
            ",".join(str(value) for value in libraries_with_review_only)
            or "NONE"
        ),
        "n_source_pair_scores": sum(
            fint(row.get("n_source_pair_scores"))
            for row in scope_summary_rows
        ),
        "n_cells_scored": len(cell_rows),
        "n_applied_reassignment_cells": sum(
            row.get("score_population_scope") == APPLIED_REASSIGNMENT
            for row in cell_rows
        ),
        "n_recommended_not_applied_cells": sum(
            row.get("score_population_scope") == RECOMMENDED_NOT_APPLIED
            for row in cell_rows
        ),
        "n_supported_reassignment_cells": sum(
            row.get("score_population_scope") in REASSIGNMENT_SCOPES
            for row in cell_rows
        ),
        "n_supported_event_held_cells": sum(
            row.get("score_population_scope") == SUPPORTED_EVENT_HELD_CELL
            for row in cell_rows
        ),
        "n_review_only_cells": len(review_only_cell_rows),
        "n_postscore_scope_excluded": sum(
            counts.get("EXCLUDED", 0)
            for counts in scope_counts_by_library.values()
        ),
        "n_supported_event_keys": expected_supported_events,
        "n_library_swap_events": len(event_rows),
        "n_library_swap_event_strata": len(event_stratum_rows),
        "n_library_swap_event_scope_rows": len(event_scope_rows),
        "n_review_only_events": len(review_event_rows),
        "n_distinct_nominated_identities": len(cross_library_rows),
        "n_candidate_rows_excluded": sum(
            fint(row.get("n_candidate_rows_excluded"))
            for row in manifest_summary_rows
        ),
        "ambient_status": ambient_provenance["status"],
        "ambient_condition": ambient_provenance["condition"],
        "ambient_contrast_path": ambient_provenance["contrast_path"],
        "ambient_donor_burden_path": ambient_provenance["donor_burden_path"],
        "score_scope_contract": PAIR_CONTRACT,
        "aggregate_scope_contract": AGGREGATE_SCOPE_CONTRACT,
        "population_supported_swap_event_classes": ",".join(
            sorted(SUPPORTED_SWAP_EVENT_CLASSES)
        ),
        "source_pair_scores_reused": str(
            args.reuse_frozen_probabilities
        ).upper(),
        "headline_probability_basis": (
            "LEGACY_FROZEN_SITE_LIKELIHOOD_FOR_QC_AND_DIRECTION_ONLY_"
            "NOT_HEADLINE_SCORE"
        ),
        "headline_display_score": "EXPERIMENTAL_DEPTH_NORMALIZED_RELATIVE_FIT",
        "relative_fit_score_status": "EXPERIMENTAL_LIBRARY19_EVALUATION",
        "relative_fit_score_basis": RELATIVE_FIT_BASIS,
        "relative_fit_score_interpretation": RELATIVE_FIT_INTERPRETATION,
        "legacy_total_likelihood_probability_role": (
            "AUDIT_COMPATIBILITY_ONLY_SATURATED_NOT_CALIBRATED_CERTAINTY"
        ),
        "molecule_probability_role": (
            "DIAGNOSTIC_UNTIL_LIBRARY19_EVALUATION;MOLECULE_PRIMARY_ROWS_UNEVALUABLE"
        ),
        "headline_qc_contract": (
            "COMPARISON_PASS;SITE_ALIGNED_ABSOLUTE_FIT_PASS;BOUNDED_SITE_"
            "PROBABILITIES;TOP_SITE_AND_TOP_FIVE_SITE_WINNER_PRESERVED;"
            "ERROR_SENSITIVITY_NOT_FAILED"
        ),
        "held_cells_vote_in_reassignment_verdict": "FALSE",
        "ties_vote_for_original": "FALSE",
        "leave_one_chromosome_metric_authoritative": "FALSE",
        "combined_influence_removal_flag_authoritative": "FALSE",
        "replacement_contig_and_block_diagnostics_status": (
            "NOT_AVAILABLE_IN_FROZEN_CPP_OUTPUT;DEFERRED_REQUIRES_NEW_SCORING"
        ),
        "validation_summary_sha256": validation_sha256,
        "aggregator_sha256": aggregator_sha256,
        "standalone_keep_scored": "FALSE",
        "unresolved_or_conflicted_scored": "FALSE",
        "standalone_cellular_origin_review_scored": "FALSE",
        "below_event_mass_reassignments_scored": "FALSE",
        "review_only_in_headline": "FALSE",
        "generic_discovery_candidates_scored": "FALSE",
        "technical_multiplets_scored_as_identities": "FALSE",
        "plot_code_generated": "FALSE",
    }
    run_summary_fields = list(run_summary)
    write_table(
        output_root / "identity_score_run_summary.tsv",
        [run_summary], run_summary_fields,
    )

    plot_files = {
        "cells_scored_for_reconciliation_swaps.tsv.gz":
            output_root / "cell_swap_identity_scores.tsv.gz",
        "swap_events_within_library.tsv":
            output_root / "library_swap_event_summary.tsv",
        "swap_event_original_assignment_strata.tsv":
            output_root / "library_swap_event_strata.tsv",
        "swap_event_scope_summary.tsv":
            output_root / "library_swap_event_scope_summary.tsv",
        "swap_candidates_across_libraries.tsv":
            output_root / "cross_library_swap_summary.tsv",
        "proposed_swap_probability_distributions.tsv":
            output_root / "proposed_swap_probability_distributions.tsv",
        "review_only_identity_scores.tsv.gz":
            output_root / "review_only_identity_scores.tsv.gz",
        "review_only_event_summary.tsv":
            output_root / "review_only_event_summary.tsv",
        "review_only_probability_distributions.tsv":
            output_root / "review_only_probability_distributions.tsv",
        "score_population_scope_summary.tsv":
            output_root / "score_population_scope_summary.tsv",
        "score_pair_exclusions.tsv":
            output_root / "score_pair_exclusion_summary.tsv",
        "ambient_donor_burden.tsv":
            output_root / "ambient_donor_burden.tsv",
        "identity_score_run_summary.tsv":
            output_root / "identity_score_run_summary.tsv",
        "identity_score_provenance.tsv":
            output_root / "identity_score_provenance.tsv",
    }
    for name, source in plot_files.items():
        shutil.copyfile(source, plot_root / name)
    write_plot_handoff(plot_root, run_summary)

    manifest_rows = []
    grains = {
        "cells_scored_for_reconciliation_swaps.tsv.gz": "one scored reconciliation-nominated cell",
        "swap_events_within_library.tsv": "one library/event/proposed identity",
        "swap_event_original_assignment_strata.tsv": "one library/event/original/proposed identity stratum",
        "swap_event_scope_summary.tsv": "one library/event/original/proposed identity/population scope",
        "swap_candidates_across_libraries.tsv": "one reconciliation-nominated biological identity",
        "proposed_swap_probability_distributions.tsv": "one event and five-point proposed-probability bin",
        "review_only_identity_scores.tsv.gz": "one separately retained unexpected-identity review cell",
        "review_only_event_summary.tsv": "one separate unexpected-identity review event",
        "review_only_probability_distributions.tsv": "one review-only event and five-point proposed-probability bin",
        "score_population_scope_summary.tsv": "one library with source-to-final scope counts",
        "score_pair_exclusions.tsv": "one library and exclusion reason",
        "ambient_donor_burden.tsv": "controlled ambient donor-burden row when available",
        "identity_score_run_summary.tsv": "one run",
        "identity_score_provenance.tsv": "one selected library and frozen input-hash contract",
        "plotting_handoff.md": "instructions",
        "data_dictionary.md": "definitions",
    }
    for name in sorted(grains):
        path = plot_root / name
        manifest_rows.append({
            "file": name,
            "grain": grains[name],
            "required": "TRUE" if name not in {"ambient_donor_burden.tsv"} else "OPTIONAL_HEADER_ONLY",
            "bytes": path.stat().st_size,
            "schema_version": SCHEMA,
        })
    write_table(
        plot_root / "plot_data_manifest.tsv",
        manifest_rows,
        ["file", "grain", "required", "bytes", "schema_version"],
    )
    print(
        f"Wrote {len(cell_rows)} supported-event cell scores across "
        f"{len(event_rows)} within-library events; retained "
        f"{len(review_only_cell_rows)} unexpected-identity review cells "
        "in separate outputs"
    )


if __name__ == "__main__":
    main()
