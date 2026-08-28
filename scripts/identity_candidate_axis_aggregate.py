#!/usr/bin/env python3
"""Aggregate fixed-pair candidate-axis results across requested libraries.

The candidate-axis position is an uncalibrated geometric location between two
fixed genotype expectations.  It is not confidence, certainty, accuracy, FDR,
a posterior probability, or a probability that either biological identity is
correct.  This helper never changes reconciliation decisions.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import importlib.util
import math
import os
import shutil
import statistics
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


PAIR_SCHEMA = "identity_candidate_axis_pair_manifest_v1"
RAW_SCHEMA = "identity_candidate_axis_pair_score_v2_sampling_adjusted_fit_diagnostic"
CELL_SCHEMA = "identity_candidate_axis_cell_v2_sampling_adjusted_fit_diagnostic"
AGGREGATE_SCHEMA = "identity_candidate_axis_aggregate_v3_multilibrary_sampling_adjusted_fit"
PAIR_AUDIT_SCHEMA = "identity_candidate_axis_pair_source_audit_v1"
SCORE_PROVENANCE_SCHEMA = "identity_candidate_axis_score_provenance_v1"
SAMPLING_ADJUSTED_BRIER_FIELDS = (
    "candidate_a_observed_brier_mean",
    "candidate_a_expected_sampling_brier_mean",
    "candidate_a_excess_brier_mean",
    "candidate_b_observed_brier_mean",
    "candidate_b_expected_sampling_brier_mean",
    "candidate_b_excess_brier_mean",
)
RAW_RESIDUAL_THRESHOLD_FLAGS = {
    "NEITHER_ABOVE_LEGACY_THRESHOLD",
    "CANDIDATE_A_ONLY_ABOVE_LEGACY_THRESHOLD",
    "CANDIDATE_B_ONLY_ABOVE_LEGACY_THRESHOLD",
    "BOTH_ABOVE_LEGACY_THRESHOLD",
    "UNAVAILABLE",
}
SAMPLING_ADJUSTED_SUMMARY_FIELDS = [
    "median_candidate_a_excess_brier_mean",
    "median_candidate_b_excess_brier_mean",
]
RELATIVE_FIT_BASIS = (
    "SITE_DELTA_LOG_LIKELIHOOD_A_MINUS_B_DIVIDED_BY_"
    "DISCRIMINATING_EVIDENCE_DEPTH"
)
RELATIVE_FIT_INTERPRETATION = (
    "EXPERIMENTAL_AVERAGE_PAIRWISE_FIT_BETWEEN_FIXED_CANDIDATES_"
    "NOT_PROBABILITY_IDENTITY_IS_CORRECT"
)
OPERATIONAL_CONTRACT = (
    "ORIGINAL_ALLOWED_DEMUX_VS_RECONCILIATION_NOMINATED_SWAP_ONLY"
)
RETAINED_CONTRACT = (
    "ORIGINAL_ALLOWED_DEMUX_VS_FROZEN_SUPPORTED_EVENT_PROPOSAL_CONTRAST_ONLY"
)
POPULATIONS = (
    "APPLIED_REASSIGNMENT",
    "RECOMMENDED_NOT_APPLIED",
    "SUPPORTED_EVENT_HELD_CELL",
    "REVIEW_ONLY_UNEXPECTED_IDENTITY",
    "RETAINED_ORIGINAL_CONTRAST_ONLY",
)
REASSIGNMENT_SCOPES = {"APPLIED_REASSIGNMENT", "RECOMMENDED_NOT_APPLIED"}
PAIR_EXCLUSION_FIELDS = [
    "exclusion_stage", "library", "barcode", "original_demux_assignment",
    "proposed_donor_genotype", "source_reconciliation_event_id",
    "source_reconciliation_proposed_identity", "reconciliation_event_id",
    "reconciliation_event_class", "reconciliation_final_action",
    "exclusion_reason", "detail", "selected_supported_event_id",
    "selected_supported_event_proposal", "score_scope_contract",
    "schema_version",
]

BASELINE_FIELDS = [
    "original_assignment_relative_fit_score_out_of_100",
    "proposed_swap_relative_fit_score_out_of_100",
    "proposed_minus_original_relative_fit_score_points",
    "relative_fit_score_status",
    "relative_fit_score_basis",
    "relative_fit_score_interpretation",
]
ESTABLISHED_FIELDS = [
    "headline_probability_status", "headline_probability_exclusion_reason",
    "primary_probability_basis", "headline_qc_basis_alignment",
    "absolute_fit_status", "original_residual_mismatch",
    "proposed_swap_residual_mismatch", "snps_favor_original",
    "snps_favor_proposed_swap", "minimum_error_sensitivity_probability_pct",
    "error_sensitivity_stable", "ambient_sensitivity_status",
    "ambient_original_c", "ambient_reconciled_c",
    "ambient_reconciled_minus_original_c", "ambient_stratum",
    "ambient_assignment_diagnostic", "reconciliation_final_action",
    "reconciliation_reassignment_applied", "comparison_outcome",
    "proposed_swap_preferred", "original_assignment_probability_pct",
    "proposed_swap_probability_pct", "probability_without_top_site_pct",
    "probability_without_top_five_sites_pct", "comparison_status",
    "probability_basis", "error_ref", "error_alt",
]

SUMMARY_METRIC_FIELDS = [
    "n_cells_total", "n_cells_candidate_axis_available",
    "n_cells_legacy_qc_eligible", "n_cells_headline_eligible",
    "n_cells_original_side", "n_cells_proposal_side", "n_cells_tied",
    "n_cells_unavailable", "candidate_axis_status_counts",
    "candidate_axis_evidence_basis_counts", "candidate_axis_position_raw_q10",
    "candidate_axis_position_raw_q25", "median_candidate_axis_position_raw",
    "candidate_axis_position_raw_q75", "candidate_axis_position_raw_q90",
    "n_available_direction_fraction_denominator",
    "n_original_side_direction_fraction_numerator",
    "n_proposal_side_direction_fraction_numerator",
    "n_tied_direction_fraction_numerator",
    "fraction_available_cells_original_side",
    "fraction_available_cells_proposal_side",
    "fraction_available_cells_tied",
    "median_observed_design_candidate_separation_rms_out_of_100",
    "median_observed_design_candidate_similarity_complement_out_of_100",
    "median_n_candidate_axis_discriminating_sites",
    "median_n_primary_evidence_units", "median_candidate_axis_design_mass",
    "median_minimum_leave_one_fold_out_candidate_axis_position",
    "median_fraction_leave_one_fold_out_positions_proposal_side",
    "candidate_axis_fold_direction_stability_status_counts",
    "n_direction_changed_after_top_primary_unit_removal",
    "n_direction_changed_after_top_five_primary_units_removal",
    "top_primary_unit_removal_status_counts",
    "top_five_primary_units_removal_status_counts",
    "raw_residual_threshold_flag_counts",
    "absolute_fit_status_legacy_counts",
    "candidate_axis_vs_v6_3_direction_status_counts",
    "v6_3_baseline_origin_counts", "n_cells_v6_3_available",
    "n_cells_v6_3_original_side", "n_cells_v6_3_proposal_side",
    "n_cells_v6_3_tied", "n_cells_v6_3_unavailable",
    "v6_3_proposed_relative_fit_score_q10_out_of_100",
    "v6_3_proposed_relative_fit_score_q25_out_of_100",
    "median_v6_3_proposed_relative_fit_score_out_of_100",
    "v6_3_proposed_relative_fit_score_q75_out_of_100",
    "v6_3_proposed_relative_fit_score_q90_out_of_100",
    "median_v6_3_original_relative_fit_score_out_of_100",
    "eligible_n_cells_candidate_axis_available",
    "eligible_candidate_axis_position_raw_q10",
    "eligible_candidate_axis_position_raw_q25",
    "eligible_median_candidate_axis_position_raw",
    "eligible_candidate_axis_position_raw_q75",
    "eligible_candidate_axis_position_raw_q90", "voting_semantics",
]


def open_text(path, mode="rt"):
    return gzip.open(path, mode, newline="") if str(path).endswith(".gz") else open(
        path, mode, newline=""
    )


def clean(value, default="NA"):
    text = "" if value is None else str(value).strip()
    return default if text in {"", ".", "NA", "nan", "NaN", "None"} else text


def finite(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def fmt(value, digits=12):
    number = finite(value)
    return "NA" if not math.isfinite(number) else f"{number:.{digits}g}"


def file_identity(path):
    target = Path(path).resolve()
    if not target.is_file():
        raise FileNotFoundError(f"required file is missing: {target}")
    stat_result = target.stat()
    return {
        "path": str(target), "size_bytes": stat_result.st_size,
        "mtime_ns": stat_result.st_mtime_ns,
    }


def read_tsv(path):
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is missing a header: {path}")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"TSV has duplicate header fields: {path}")
        return [dict(row) for row in reader]


def read_tsv_fields(path):
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV is missing a header: {path}")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"TSV has duplicate header fields: {path}")
        return list(reader.fieldnames)


def write_tsv(path, rows, fields):
    path = Path(path)
    temporary = path.with_name(
        path.stem + f".tmp.{os.getpid()}" + path.suffix
    )
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    try:
        with open_text(temporary, "wt") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
                extrasaction="raise",
            )
            writer.writeheader()
            for row in rows:
                rendered = {}
                for field in fields:
                    value = clean(row.get(field))
                    if value.lower() in {"nan", "inf", "+inf", "-inf"}:
                        raise ValueError(
                            f"nonfinite token in {path}: field={field} value={value}"
                        )
                    rendered[field] = value
                writer.writerow(rendered)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def unique_rows(rows, key_fields, label):
    result = {}
    for row in rows:
        key = tuple(clean(row.get(field), "") for field in key_fields)
        if any(not value for value in key):
            raise ValueError(f"{label} has a blank key: {key}")
        if key in result:
            raise ValueError(f"{label} has a duplicate key: {key}")
        result[key] = row
    return result


def canonical_genotype(value):
    raw = clean(value, "")
    if not raw:
        return ""
    if raw.startswith("M{"):
        return raw
    for prefix in ("D[", "T[", "UNKNOWN_SINGLE_CELL["):
        if raw.startswith(prefix) and raw.endswith("]"):
            raw = raw[len(prefix):-1]
            break
    return "+".join(sorted(part.strip() for part in raw.replace("x", "+").split("+") if part.strip()))


def import_v6(path):
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("candidate_axis_frozen_v6_3", path)
    if spec is None or spec.loader is None:
        raise ValueError(f"could not load V6.3 source: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    required = (
        "relative_fit_scores", "primary_site_probabilities",
        "headline_probability_eligibility", "build_cell_row", "fmt",
        "RELATIVE_FIT_BASIS", "RELATIVE_FIT_INTERPRETATION",
        "TIE_TOLERANCE_PP",
    )
    missing = [name for name in required if not hasattr(module, name)]
    if missing:
        raise ValueError("V6.3 source is missing required symbols: " + ", ".join(missing))
    return module


def parse_input_specs(values):
    result = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"input identity must be LABEL=PATH, saw {value!r}")
        label, path = value.split("=", 1)
        if not label or label in result:
            raise ValueError(f"blank/duplicate input identity label: {label!r}")
        target = Path(path).resolve()
        if not target.is_file():
            raise FileNotFoundError(f"required scorer input is missing: {target}")
        result[label] = str(target)
    return result


def score_provenance_main(argv):
    parser = argparse.ArgumentParser(description="Write simple scorer run metadata")
    parser.add_argument("--output", required=True)
    parser.add_argument("--library", required=True)
    parser.add_argument("--command", required=True)
    parser.add_argument("--input", action="append", default=[])
    parser.add_argument("--raw-output", required=True)
    parser.add_argument("--scorer-binary", required=True)
    parser.add_argument("--temp-root", required=True)
    parser.add_argument("--error-ref", required=True)
    parser.add_argument("--error-alt", required=True)
    parser.add_argument("--min-evidence", required=True)
    parser.add_argument("--poor-fit-residual", required=True)
    args = parser.parse_args(argv)
    inputs = parse_input_specs(args.input)
    scorer_binary = file_identity(args.scorer_binary)
    raw = file_identity(args.raw_output)
    pair_manifest_rows = data_row_count(inputs["pair_manifest"])
    pair_summary_rows = read_tsv(inputs["pair_summary"])
    if len(pair_summary_rows) != 1:
        raise ValueError("candidate-axis pair summary must contain exactly one row")
    rows = [{
        "library": args.library, "command": args.command,
        "pair_manifest_path": inputs["pair_manifest"],
        "pair_summary_path": inputs["pair_summary"],
        "samples_path": inputs["samples"],
        "pileup_sites_path": inputs["pileup_sites"],
        "pileup_observations_path": inputs["pileup_observations"],
        "raw_output_path": raw["path"], "raw_output_size_bytes": raw["size_bytes"],
        "raw_output_mtime_ns": raw["mtime_ns"],
        "raw_output_row_count": data_row_count(args.raw_output),
        "pair_manifest_row_count": pair_manifest_rows,
        "pair_summary_row_count": len(pair_summary_rows),
        "pair_summary_n_score_pairs": clean(
            pair_summary_rows[0].get("n_score_pairs"), "0"
        ),
        "scorer_binary_path": scorer_binary["path"],
        "scorer_binary_size_bytes": scorer_binary["size_bytes"],
        "scorer_binary_mtime_ns": scorer_binary["mtime_ns"],
        "candidate_axis_temp_root": str(Path(args.temp_root).resolve()),
        "error_ref": args.error_ref, "error_alt": args.error_alt,
        "min_evidence": args.min_evidence,
        "poor_fit_residual": args.poor_fit_residual,
        "schema_version": SCORE_PROVENANCE_SCHEMA,
    }]
    fields = list(rows[0])
    write_tsv(Path(args.output).resolve(), rows, fields)
    return 0


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--identity-root", required=True)
    parser.add_argument("--libraries", nargs="+", required=True)
    parser.add_argument("--event-id", default="")
    parser.add_argument("--proposed-identity", default="")
    parser.add_argument("--pair-root", required=True)
    parser.add_argument("--score-root", required=True)
    parser.add_argument("--v6-3-root", default="")
    parser.add_argument("--ambient-path", default="")
    parser.add_argument("--ambient-condition", default="")
    parser.add_argument("--output-root", required=True)
    return parser.parse_args(argv)


def parse_library(values):
    numbers = set()
    for value in values:
        for token in str(value).replace(",", " ").split():
            if "-" in token:
                left, right = token.split("-", 1)
                numbers.update(range(int(left), int(right) + 1))
            else:
                token = token[3:] if token.lower().startswith("lib") else token
                numbers.add(int(token))
    if not numbers:
        raise ValueError("candidate-axis aggregate requires at least one library")
    return [f"lib{number}" for number in sorted(numbers)]


def load_manifest_pairs(path, library, event_id="", proposal=""):
    required_fields = {
        "schema_version", "library", "barcode", "score_pair_id",
        "score_pair_role", "score_pair_source", "score_population_scope",
        "population_votes_in_authoritative_event", "supported_event_key",
        "selected_supported_event_id", "selected_supported_event_proposal",
        "reconciliation_event_id", "original_demux_assignment",
        "proposed_donor_genotype", "candidate_b_fixed_identity",
        "score_scope_contract", "candidate_origin", "donor_genotype",
    }
    missing_fields = required_fields - set(read_tsv_fields(path))
    if missing_fields:
        raise ValueError(
            f"candidate-axis pair manifest is missing fields {sorted(missing_fields)}: {path}"
        )
    rows = read_tsv(path)
    pairs = defaultdict(list)
    for row in rows:
        if clean(row.get("schema_version"), "") != PAIR_SCHEMA:
            raise ValueError(f"candidate-axis pair manifest schema mismatch: {path}")
        if clean(row.get("library"), "") != library:
            raise ValueError(f"candidate-axis pair manifest has a cross-library row: {path}")
        row_event = clean(row.get("selected_supported_event_id"), "")
        row_proposal = canonical_genotype(
            row.get("selected_supported_event_proposal")
        )
        if not row_event or not row_proposal:
            raise ValueError("candidate-axis pair manifest has a blank event key")
        if event_id and (row_event != event_id or row_proposal != proposal):
            raise ValueError("candidate-axis pair manifest targeted event mismatch")
        pairs[(library, clean(row.get("barcode"), ""), clean(row.get("score_pair_id"), ""))].append(row)
    result = {}
    shared = [
        "library", "barcode", "score_pair_id", "score_pair_source",
        "score_population_scope", "population_votes_in_authoritative_event",
        "supported_event_key", "selected_supported_event_id",
        "selected_supported_event_proposal", "reconciliation_event_id",
        "reconciliation_event_class", "reconciliation_event_confidence",
        "reconciliation_final_action", "reconciliation_decision_confidence",
        "reconciliation_reassignment_applied", "original_demux_assignment",
        "proposed_donor_genotype", "reconciliation_nominated_swap",
        "candidate_b_fixed_identity", "pair_construction_mode",
        "score_scope_contract", "schema_version",
    ]
    for key, pair_rows in pairs.items():
        if len(pair_rows) != 2:
            raise ValueError(f"candidate-axis pair must contain exactly two rows: {key}")
        by_role = {clean(row.get("score_pair_role"), ""): row for row in pair_rows}
        if len(by_role) != 2 or "ORIGINAL_ALLOWED_DEMUX" not in by_role:
            raise ValueError(f"candidate-axis pair has duplicate/missing roles: {key}")
        original = by_role["ORIGINAL_ALLOWED_DEMUX"]
        b_roles = set(by_role) - {"ORIGINAL_ALLOWED_DEMUX"}
        if b_roles not in (
            {"RECONCILIATION_NOMINATED_SWAP"},
            {"FROZEN_SUPPORTED_EVENT_PROPOSAL_CONTRAST"},
        ):
            raise ValueError(f"candidate-axis pair has unsupported B role: {key}")
        proposed = by_role[next(iter(b_roles))]
        for field in shared:
            if clean(original.get(field)) != clean(proposed.get(field)):
                raise ValueError(f"candidate-axis pair metadata mismatch: {key} {field}")
        if canonical_genotype(original.get("donor_genotype")) != canonical_genotype(
            original.get("original_demux_assignment")
        ):
            raise ValueError(f"candidate A is not the original allowed demux identity: {key}")
        pair_event = clean(original.get("selected_supported_event_id"), "")
        pair_proposal = canonical_genotype(
            original.get("selected_supported_event_proposal")
        )
        if canonical_genotype(proposed.get("donor_genotype")) != pair_proposal:
            raise ValueError(f"candidate B is not the selected fixed proposal: {key}")
        expected_event_key = f"{library}|{pair_event}|{pair_proposal}"
        if clean(original.get("supported_event_key"), "") != expected_event_key:
            raise ValueError(f"candidate-axis supported event key mismatch: {key}")
        if clean(original.get("reconciliation_event_id"), "") != pair_event:
            raise ValueError(f"candidate-axis reconciliation event mismatch: {key}")
        population = clean(original.get("score_population_scope"), "")
        if population not in POPULATIONS:
            raise ValueError(f"unsupported candidate-axis population: {population}")
        retained = population == "RETAINED_ORIGINAL_CONTRAST_ONLY"
        expected_contract = RETAINED_CONTRACT if retained else OPERATIONAL_CONTRACT
        if clean(original.get("score_scope_contract"), "") != expected_contract:
            raise ValueError(f"candidate-axis pair contract mismatch: {key}")
        votes = clean(original.get("population_votes_in_authoritative_event"), "").upper()
        if votes != ("TRUE" if population in REASSIGNMENT_SCOPES else "FALSE"):
            raise ValueError(f"candidate-axis pair voting annotation mismatch: {key}")
        result[key] = {
            "original": original, "proposed": proposed, "retained": retained,
            "event_id": pair_event, "proposal": pair_proposal,
        }
    return result, rows


def validate_raw_join(manifest, raw, key):
    if clean(raw.get("schema_version"), "") != RAW_SCHEMA:
        raise ValueError(f"raw candidate-axis schema mismatch: {key}")
    required_fit_fields = set(SAMPLING_ADJUSTED_BRIER_FIELDS) | {
        "raw_residual_threshold_flag",
        "candidate_a_residual_mismatch_legacy",
        "candidate_b_residual_mismatch_legacy",
        "absolute_fit_status_legacy",
        "comparison_status_legacy",
        "poor_fit_residual_threshold_legacy",
        "n_candidate_axis_discriminating_sites",
    }
    missing_fit_fields = required_fit_fields - set(raw)
    if missing_fit_fields:
        raise ValueError(
            f"raw candidate-axis sampling-adjusted fit fields are missing: "
            f"{key} {sorted(missing_fit_fields)}"
        )
    discriminating_value = finite(raw.get("n_candidate_axis_discriminating_sites"))
    if not math.isfinite(discriminating_value) or discriminating_value < 0 or \
            discriminating_value != int(discriminating_value):
        raise ValueError(
            f"raw candidate-axis discriminating-site count is invalid: {key}"
        )
    discriminating_sites = int(discriminating_value)
    brier_values = [finite(raw.get(field)) for field in SAMPLING_ADJUSTED_BRIER_FIELDS]
    if discriminating_sites > 0 and not all(
            math.isfinite(value) for value in brier_values):
        raise ValueError(
            f"raw candidate-axis Brier diagnostics are unavailable despite "
            f"discriminating sites: {key}"
        )
    if discriminating_sites == 0 and any(
            math.isfinite(value) for value in brier_values):
        raise ValueError(
            f"raw candidate-axis Brier diagnostics must be unavailable without "
            f"discriminating sites: {key}"
        )
    residual_a = finite(raw.get("candidate_a_residual_mismatch_legacy"))
    residual_b = finite(raw.get("candidate_b_residual_mismatch_legacy"))
    threshold = finite(raw.get("poor_fit_residual_threshold_legacy"))
    flag = clean(raw.get("raw_residual_threshold_flag"), "")
    if flag not in RAW_RESIDUAL_THRESHOLD_FLAGS:
        raise ValueError(f"raw candidate-axis residual threshold flag is invalid: {key}")
    if not all(math.isfinite(value) for value in (residual_a, residual_b, threshold)):
        expected_flag = "UNAVAILABLE"
    elif residual_a > threshold and residual_b > threshold:
        expected_flag = "BOTH_ABOVE_LEGACY_THRESHOLD"
    elif residual_a > threshold:
        expected_flag = "CANDIDATE_A_ONLY_ABOVE_LEGACY_THRESHOLD"
    elif residual_b > threshold:
        expected_flag = "CANDIDATE_B_ONLY_ABOVE_LEGACY_THRESHOLD"
    else:
        expected_flag = "NEITHER_ABOVE_LEGACY_THRESHOLD"
    if flag != expected_flag:
        raise ValueError(
            f"raw candidate-axis residual threshold flag disagrees with legacy "
            f"residuals: {key} expected={expected_flag} observed={flag}"
        )
    original, proposed = manifest["original"], manifest["proposed"]
    comparisons = {
        "library": (original.get("library"), raw.get("library")),
        "barcode": (original.get("barcode"), raw.get("barcode")),
        "score_pair_id": (original.get("score_pair_id"), raw.get("score_pair_id")),
        "supported_event_key": (original.get("supported_event_key"), raw.get("supported_event_key")),
        "candidate_a": (original.get("donor_genotype"), raw.get("candidate_a")),
        "candidate_b": (proposed.get("donor_genotype"), raw.get("candidate_b")),
        "candidate_a_role": (original.get("score_pair_role"), raw.get("candidate_a_role")),
        "candidate_b_role": (proposed.get("score_pair_role"), raw.get("candidate_b_role")),
        "candidate_a_origin": (original.get("candidate_origin"), raw.get("candidate_a_origin")),
        "candidate_b_origin": (proposed.get("candidate_origin"), raw.get("candidate_b_origin")),
        "score_pair_source": (original.get("score_pair_source"), raw.get("score_pair_source")),
        "score_population_scope": (original.get("score_population_scope"), raw.get("score_population_scope")),
        "population_votes_in_authoritative_event": (
            original.get("population_votes_in_authoritative_event"),
            raw.get("population_votes_in_authoritative_event"),
        ),
        "score_scope_contract": (original.get("score_scope_contract"), raw.get("score_scope_contract")),
    }
    for field, (expected, observed) in comparisons.items():
        mismatch = (
            canonical_genotype(expected) != canonical_genotype(observed)
            if field in {"candidate_a", "candidate_b"}
            else clean(expected) != clean(observed)
        )
        if mismatch:
            raise ValueError(
                f"manifest/raw candidate-axis join mismatch: {key} {field} "
                f"expected={expected!r} observed={observed!r}"
            )


def load_optional_frozen_v6(root, library):
    if not root:
        return {}, {}, {"cell": "NA", "review": "NA"}
    root = Path(root).resolve()
    cell_path = root / "cell_swap_identity_scores.tsv.gz"
    review_path = root / "review_only_identity_scores.tsv.gz"
    if not cell_path.is_file() or not review_path.is_file():
        return {}, {}, {"cell": "NA", "review": "NA"}
    def index(path, label):
        indexed = {}
        for row in read_tsv(path):
            key = (normalize_library(row.get("library")), clean(row.get("barcode"), ""))
            if key[0] != library:
                continue
            if not key[1] or key in indexed:
                raise ValueError(f"{label} has a blank or duplicate key: {key}")
            indexed[key] = row
        return indexed
    cell = index(cell_path, "frozen V6.3 cell")
    review = index(review_path, "frozen V6.3 review-only")
    if set(cell) & set(review):
        raise ValueError("frozen V6.3 cell and review-only tables are not disjoint")
    return cell, review, {
        "cell": file_identity(cell_path), "review": file_identity(review_path)
    }


def v6_direction_from_row(row, tie_tolerance):
    preferred = clean(row.get("proposed_swap_preferred"), "").upper()
    outcome = clean(row.get("comparison_outcome"), "").upper()
    if preferred == "TRUE":
        if outcome and outcome != "PROPOSED_SWAP_PREFERRED":
            raise ValueError("frozen V6.3 proposed_swap_preferred/outcome mismatch")
        return "PROPOSAL_SIDE"
    if preferred == "FALSE":
        if outcome and outcome != "ORIGINAL_PREFERRED":
            raise ValueError("frozen V6.3 proposed_swap_preferred/outcome mismatch")
        return "ORIGINAL_SIDE"
    if preferred == "TIE" or outcome == "TIE":
        return "TIE"
    proposed = finite(row.get("proposed_swap_probability_pct"))
    original = finite(row.get("original_assignment_probability_pct"))
    if not (math.isfinite(proposed) and math.isfinite(original)):
        return "UNAVAILABLE"
    gap = proposed - original
    return "PROPOSAL_SIDE" if gap > tie_tolerance else (
        "ORIGINAL_SIDE" if gap < -tie_tolerance else "TIE"
    )


def compare_directions(axis_direction, v6_direction):
    if axis_direction == "UNAVAILABLE":
        return "CANDIDATE_AXIS_UNAVAILABLE"
    if v6_direction == "UNAVAILABLE":
        return "V6_3_UNAVAILABLE"
    if axis_direction == "TIE" and v6_direction == "TIE":
        return "BOTH_TIED"
    if axis_direction == "TIE":
        return "CANDIDATE_AXIS_TIE"
    if v6_direction == "TIE":
        return "V6_3_TIE"
    return "CONCORDANT" if axis_direction == v6_direction else "DISCORDANT"


def retained_probability_row(raw):
    delta = finite(raw.get(
        "v6_3_compatible_site_delta_log_likelihood_a_minus_b"
    ))
    preferred = "NA" if not math.isfinite(delta) else (
        clean(raw.get("candidate_a")) if delta >= 0
        else clean(raw.get("candidate_b"))
    )
    return {
        "comparison_status": clean(raw.get("comparison_status_legacy")),
        "probability_basis": clean(raw.get("probability_basis_legacy")),
        "site_delta_log_likelihood_a_minus_b": clean(
            raw.get("v6_3_compatible_site_delta_log_likelihood_a_minus_b")
        ),
        "discriminating_evidence_depth": clean(
            raw.get("v6_3_compatible_discriminating_evidence_depth")
        ),
        "site_candidate_a_probability_pct": clean(
            raw.get("v6_3_compatible_site_candidate_a_probability_pct")
        ),
        "site_candidate_b_probability_pct": clean(
            raw.get("v6_3_compatible_site_candidate_b_probability_pct")
        ),
        "candidate_a_probability_pct": clean(raw.get("site_candidate_a_probability_pct")),
        "candidate_b_probability_pct": clean(raw.get("site_candidate_b_probability_pct")),
        "molecule_evidence_status": "MOLECULE_SIDECAR_UNAVAILABLE",
        "n_independent_molecules": "NA",
        "absolute_fit_status": clean(raw.get("absolute_fit_status_legacy")),
        "probability_without_top_site_pct": clean(
            raw.get("probability_without_top_site_pct_legacy")
        ),
        "probability_without_top_five_sites_pct": clean(
            raw.get("probability_without_top_five_sites_pct_legacy")
        ),
        "minimum_error_sensitivity_probability_pct": clean(
            raw.get("minimum_error_sensitivity_probability_pct_legacy")
        ),
        "error_sensitivity_stable": clean(raw.get("error_sensitivity_stable_legacy")),
        "preferred_assignment": preferred,
        "n_snps_favor_preferred": (
            "NA" if not math.isfinite(delta) else (
                clean(raw.get("n_sites_favor_candidate_a"))
                if delta >= 0 else clean(raw.get("n_sites_favor_candidate_b"))
            )
        ),
        "n_snps_favor_alternative": (
            "NA" if not math.isfinite(delta) else (
                clean(raw.get("n_sites_favor_candidate_b"))
                if delta >= 0 else clean(raw.get("n_sites_favor_candidate_a"))
            )
        ),
        "preferred_residual_mismatch": (
            "NA" if not math.isfinite(delta) else (
                clean(raw.get("candidate_a_residual_mismatch_legacy"))
                if delta >= 0 else clean(raw.get("candidate_b_residual_mismatch_legacy"))
            )
        ),
        "alternative_residual_mismatch": (
            "NA" if not math.isfinite(delta) else (
                clean(raw.get("candidate_b_residual_mismatch_legacy"))
                if delta >= 0 else clean(raw.get("candidate_a_residual_mismatch_legacy"))
            )
        ),
        "n_discriminating_snps": clean(raw.get("n_candidate_axis_discriminating_sites")),
        "n_common_observed_snps": clean(raw.get("n_common_observed_nuclear_sites")),
        "common_evidence_depth": clean(raw.get("total_common_observed_evidence_depth")),
        "error_ref": clean(raw.get("error_ref")),
        "error_alt": clean(raw.get("error_alt")),
        "ambient_sensitivity_status": "NOT_EVALUATED_NO_FROZEN_ALLELE_PROFILE",
        "warnings": clean(raw.get("warnings"), "NONE"),
    }


def frozen_v3_error_pair(rows):
    values = {}
    for field in ("error_ref", "error_alt"):
        observed = {finite(row.get(field)) for row in rows}
        if any(not math.isfinite(value) for value in observed):
            raise ValueError(f"frozen V3 {field} contains a nonfinite value")
        if len(observed) != 1:
            raise ValueError(
                f"frozen V3 {field} must be one unanimous run-level value"
            )
        values[field] = next(iter(observed))
    if not (
        0 <= values["error_ref"] <= 1
        and 0 <= values["error_alt"] <= 1
        and values["error_ref"] + values["error_alt"] < 1
    ):
        raise ValueError("frozen V3 error pair is outside the transform domain")
    return values["error_ref"], values["error_alt"]


def v6_fields_from_probability(module, probability):
    result = module.relative_fit_scores(probability)
    return {
        "original_assignment_relative_fit_score_out_of_100": module.fmt(
            result["original"]
        ),
        "proposed_swap_relative_fit_score_out_of_100": module.fmt(
            result["proposed"]
        ),
        "proposed_minus_original_relative_fit_score_points": module.fmt(
            result["gap"]
        ),
        "relative_fit_score_status": result["status"],
        "relative_fit_score_basis": module.RELATIVE_FIT_BASIS,
        "relative_fit_score_interpretation": module.RELATIVE_FIT_INTERPRETATION,
    }


def v6_direction_from_probability(module, probability):
    candidate_a, candidate_b, _, _ = module.primary_site_probabilities(probability)
    if not (math.isfinite(candidate_a) and math.isfinite(candidate_b)):
        return "UNAVAILABLE"
    gap = candidate_b - candidate_a
    return "PROPOSAL_SIDE" if gap > module.TIE_TOLERANCE_PP else (
        "ORIGINAL_SIDE" if gap < -module.TIE_TOLERANCE_PP else "TIE"
    )


def reconstruct_established(module, manifest, probability, library, retained):
    bundle = {
        "original": manifest["original"],
        "proposed": manifest["proposed"],
        "original_identity": canonical_genotype(
            manifest["original"].get("donor_genotype")
        ),
        "proposed_identity": canonical_genotype(
            manifest["proposed"].get("donor_genotype")
        ),
        "score_pair_id": clean(manifest["original"].get("score_pair_id")),
    }
    row = module.build_cell_row(
        int(library[3:]),
        clean(manifest["original"].get("barcode")),
        bundle,
        probability,
        None,
        clean(manifest["original"].get("score_population_scope")),
        clean(manifest["original"].get("supported_event_key")),
    )
    row["library"] = library
    if retained:
        eligible, reason = module.headline_probability_eligibility(
            probability,
            *module.primary_site_probabilities(probability)[:2],
            module.primary_site_probabilities(probability)[2],
        )
        row["legacy_pair_qc_status"] = "ELIGIBLE" if eligible else "UNEVALUABLE"
        row["legacy_pair_qc_exclusion_reason"] = reason
        row["computed_component_headline_probability_status"] = (
            "ELIGIBLE" if eligible else "UNEVALUABLE"
        )
        row["computed_component_headline_probability_exclusion_reason"] = reason
        row["headline_probability_status"] = "UNEVALUABLE"
        row["headline_probability_exclusion_reason"] = (
            "NONHEADLINE_RETAINED_ORIGINAL_CONTRAST"
        )
        row["reconciliation_nominated_swap"] = "NA"
    else:
        row["legacy_pair_qc_status"] = row["headline_probability_status"]
        row["legacy_pair_qc_exclusion_reason"] = row[
            "headline_probability_exclusion_reason"
        ]
        row["computed_component_headline_probability_status"] = row[
            "headline_probability_status"
        ]
        row["computed_component_headline_probability_exclusion_reason"] = row[
            "headline_probability_exclusion_reason"
        ]
    return row


def build_cell(
    module,
    library,
    key,
    manifest,
    raw,
    decision,
    frozen_v3,
    frozen_v6,
):
    validate_raw_join(manifest, raw, key)
    retained = manifest["retained"]
    original = canonical_genotype(manifest["original"].get("donor_genotype"))
    proposed = canonical_genotype(manifest["proposed"].get("donor_genotype"))
    if decision is None:
        raise ValueError(f"manifest row is missing from the decision table: {key}")
    if canonical_genotype(decision.get("original_demux_assignment")) != original:
        raise ValueError(f"decision/manifest original identity mismatch: {key}")

    if retained:
        probability = retained_probability_row(raw)
        established = reconstruct_established(
            module, manifest, probability, library, retained=True
        )
        relative = v6_fields_from_probability(module, probability)
        baseline_origin = "RECOMPUTED_FOR_EVENT_CONTRAST_WITH_V6_3"
        v6_direction = v6_direction_from_probability(module, probability)
    else:
        expected_old_pair_id = f"{library}:{clean(raw.get('barcode'))}:RECONCILIATION_SWAP"
        if frozen_v3 is not None:
            if canonical_genotype(frozen_v3.get("candidate_a")) != original or \
                    canonical_genotype(frozen_v3.get("candidate_b")) != proposed:
                raise ValueError(f"frozen V3 candidate orientation mismatch: {key}")
            if clean(frozen_v3.get("score_pair_id")) != expected_old_pair_id:
                raise ValueError(f"frozen V3 score_pair_id mismatch: {key}")
            new_delta = finite(raw.get("site_delta_log_likelihood_a_minus_b"))
            old_delta = finite(frozen_v3.get("site_delta_log_likelihood_a_minus_b"))
            if math.isfinite(new_delta) and math.isfinite(old_delta):
                if (new_delta > 0) != (old_delta > 0) and \
                        new_delta != 0 and old_delta != 0:
                    raise ValueError(
                        f"new/frozen legacy site-likelihood direction mismatch: {key}"
                    )
        if frozen_v6 is not None:
            for field in (
                "library", "barcode", "score_pair_id", "score_population_scope",
                "supported_event_key", "reconciliation_event_id",
                "original_allowed_demux_assignment", "reconciliation_nominated_swap",
            ):
                expected = {
                    "library": library,
                    "barcode": clean(raw.get("barcode")),
                    "score_pair_id": expected_old_pair_id,
                    "score_population_scope": clean(
                        manifest["original"].get("score_population_scope")
                    ),
                    "supported_event_key": clean(
                        manifest["original"].get("supported_event_key")
                    ),
                    "reconciliation_event_id": clean(
                        manifest["original"].get("selected_supported_event_id")
                    ),
                    "original_allowed_demux_assignment": original,
                    "reconciliation_nominated_swap": proposed,
                }[field]
                observed = frozen_v6.get(field)
                mismatch = (
                    canonical_genotype(expected) != canonical_genotype(observed)
                    if field in {
                        "original_allowed_demux_assignment",
                        "reconciliation_nominated_swap",
                    }
                    else clean(expected) != clean(observed)
                )
                if mismatch:
                    raise ValueError(f"frozen V6.3 join mismatch: {key} {field}")
            established = dict(frozen_v6)
            relative = {field: clean(frozen_v6.get(field)) for field in BASELINE_FIELDS}
            baseline_origin = "FROZEN_EXISTING"
            v6_direction = v6_direction_from_row(
                frozen_v6, module.TIE_TOLERANCE_PP
            )
        elif frozen_v3 is not None:
            established = reconstruct_established(
                module, manifest, frozen_v3, library, retained=False
            )
            relative = v6_fields_from_probability(module, frozen_v3)
            baseline_origin = (
                "RECOMPUTED_FROM_FROZEN_V3_FIELDS_WITH_V6_3"
            )
            v6_direction = v6_direction_from_probability(module, frozen_v3)
        else:
            probability = retained_probability_row(raw)
            probability["error_ref"] = "0.001"
            probability["error_alt"] = "0.001"
            established = reconstruct_established(
                module, manifest, probability, library, retained=False
            )
            relative = v6_fields_from_probability(module, probability)
            baseline_origin = (
                "RECOMPUTED_FROM_CANDIDATE_AXIS_SCORE_WITH_V6_3"
            )
            v6_direction = v6_direction_from_probability(module, probability)

    axis_direction = clean(raw.get("candidate_axis_direction"), "UNAVAILABLE")
    cell = dict(established)
    for field, value in raw.items():
        if field == "schema_version":
            cell["raw_score_schema_version"] = value
        elif field not in ESTABLISHED_FIELDS and field not in BASELINE_FIELDS:
            cell[field] = value
    for legacy_alias in (
        "absolute_fit_status",
        "comparison_status",
        "original_residual_mismatch",
        "proposed_swap_residual_mismatch",
        "preferred_residual_mismatch",
        "alternative_residual_mismatch",
    ):
        cell.pop(legacy_alias, None)
    cell.update(relative)
    cell.update({
        "library": library,
        "barcode": clean(raw.get("barcode")),
        "score_pair_id": clean(raw.get("score_pair_id")),
        "selected_supported_event_id": clean(
            manifest["original"].get("selected_supported_event_id")
        ),
        "selected_supported_event_proposal": proposed,
        "candidate_b_fixed_identity": proposed,
        "source_reconciliation_event_id": clean(
            manifest["original"].get("source_reconciliation_event_id")
        ),
        "source_reconciliation_proposed_identity": canonical_genotype(
            manifest["original"].get("source_reconciliation_proposed_identity")
        ) or "NA",
        "population_votes_in_authoritative_event": clean(
            manifest["original"].get("population_votes_in_authoritative_event")
        ),
        "v6_3_baseline_origin": baseline_origin,
        "candidate_axis_vs_v6_3_direction_status": compare_directions(
            axis_direction, v6_direction
        ),
        "v6_3_categorical_direction": v6_direction,
        "candidate_axis_probability_interpretation": (
            "FIXED_PAIR_GEOMETRIC_POSITION_AND_BRIER_DIRECTION_NOT_"
            "PROBABILITY_IDENTITY_IS_CORRECT"
        ),
        "reconciliation_nominated_swap": proposed if not retained else "NA",
        "reconciliation_final_action": clean(decision.get("final_action")),
        "reconciliation_reassignment_applied": clean(
            decision.get("reassignment_applied")
        ),
        "original_demux_assignment": original,
        "schema_version": CELL_SCHEMA,
    })
    if "legacy_pair_qc_status" not in cell:
        eligible = clean(cell.get("headline_probability_status"), "") == "ELIGIBLE"
        cell["legacy_pair_qc_status"] = "ELIGIBLE" if eligible else "UNEVALUABLE"
        cell["legacy_pair_qc_exclusion_reason"] = (
            "NONE" if eligible else clean(
                cell.get("headline_probability_exclusion_reason"), "UNEVALUABLE"
            )
        )
    if retained:
        cell["snps_favor_original"] = clean(raw.get("n_sites_favor_candidate_a"))
        cell["snps_favor_proposed_swap"] = clean(
            raw.get("n_sites_favor_candidate_b")
        )
    return cell


def quantile(values, q):
    numbers = sorted(value for value in (finite(item) for item in values) if math.isfinite(value))
    if not numbers:
        return math.nan
    if len(numbers) == 1:
        return numbers[0]
    position = q * (len(numbers) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return numbers[lower]
    fraction = position - lower
    return numbers[lower] * (1.0 - fraction) + numbers[upper] * fraction


def counter_text(values):
    counts = Counter(clean(value) for value in values)
    return ",".join(f"{key}:{counts[key]}" for key in sorted(counts)) or "NONE"


def is_available(row):
    return clean(row.get("candidate_axis_status"), "") == "AVAILABLE"


def is_legacy_eligible(row):
    return clean(row.get("legacy_pair_qc_status"), "") == "ELIGIBLE"


def summary_metrics(rows, include_sampling_adjusted_fit=False):
    available = [row for row in rows if is_available(row)]
    eligible = [row for row in available if is_legacy_eligible(row)]
    directions = Counter(clean(row.get("candidate_axis_direction")) for row in available)
    positions = [row.get("candidate_axis_position_raw") for row in available]
    eligible_positions = [row.get("candidate_axis_position_raw") for row in eligible]
    v6_available = [
        row for row in rows
        if clean(row.get("relative_fit_score_status"), "") == "AVAILABLE"
    ]
    v6_directions = Counter(clean(row.get("v6_3_categorical_direction")) for row in rows)
    result = {
        "n_cells_total": len(rows),
        "n_cells_candidate_axis_available": len(available),
        "n_cells_legacy_qc_eligible": sum(is_legacy_eligible(row) for row in rows),
        "n_cells_headline_eligible": sum(
            clean(row.get("headline_probability_status"), "") == "ELIGIBLE"
            for row in rows
        ),
        "n_cells_original_side": directions["ORIGINAL_SIDE"],
        "n_cells_proposal_side": directions["PROPOSAL_SIDE"],
        "n_cells_tied": directions["TIE"],
        "n_cells_unavailable": len(rows) - len(available),
        "candidate_axis_status_counts": counter_text(
            row.get("candidate_axis_status") for row in rows
        ),
        "candidate_axis_evidence_basis_counts": counter_text(
            row.get("candidate_axis_evidence_basis") for row in rows
        ),
        "candidate_axis_position_raw_q10": fmt(quantile(positions, 0.10)),
        "candidate_axis_position_raw_q25": fmt(quantile(positions, 0.25)),
        "median_candidate_axis_position_raw": fmt(quantile(positions, 0.50)),
        "candidate_axis_position_raw_q75": fmt(quantile(positions, 0.75)),
        "candidate_axis_position_raw_q90": fmt(quantile(positions, 0.90)),
        "n_available_direction_fraction_denominator": len(available),
        "n_original_side_direction_fraction_numerator": directions["ORIGINAL_SIDE"],
        "n_proposal_side_direction_fraction_numerator": directions["PROPOSAL_SIDE"],
        "n_tied_direction_fraction_numerator": directions["TIE"],
        "fraction_available_cells_original_side": fmt(
            directions["ORIGINAL_SIDE"] / len(available) if available else math.nan
        ),
        "fraction_available_cells_proposal_side": fmt(
            directions["PROPOSAL_SIDE"] / len(available) if available else math.nan
        ),
        "fraction_available_cells_tied": fmt(
            directions["TIE"] / len(available) if available else math.nan
        ),
        "median_observed_design_candidate_separation_rms_out_of_100": fmt(quantile(
            [row.get("observed_design_candidate_separation_rms_out_of_100") for row in available], 0.5
        )),
        "median_observed_design_candidate_similarity_complement_out_of_100": fmt(quantile(
            [row.get("observed_design_candidate_similarity_complement_out_of_100") for row in available], 0.5
        )),
        "median_n_candidate_axis_discriminating_sites": fmt(quantile(
            [row.get("n_candidate_axis_discriminating_sites") for row in available], 0.5
        )),
        "median_n_primary_evidence_units": fmt(quantile(
            [row.get("n_primary_evidence_units") for row in available], 0.5
        )),
        "median_candidate_axis_design_mass": fmt(quantile(
            [row.get("candidate_axis_design_mass") for row in available], 0.5
        )),
        "median_minimum_leave_one_fold_out_candidate_axis_position": fmt(quantile(
            [row.get("minimum_leave_one_fold_out_candidate_axis_position") for row in available], 0.5
        )),
        "median_fraction_leave_one_fold_out_positions_proposal_side": fmt(quantile(
            [row.get("fraction_leave_one_fold_out_positions_proposal_side") for row in available], 0.5
        )),
        "candidate_axis_fold_direction_stability_status_counts": counter_text(
            row.get("candidate_axis_fold_direction_stability_status") for row in rows
        ),
        "n_direction_changed_after_top_primary_unit_removal": sum(
            clean(row.get("candidate_axis_direction_preserved_without_top_primary_unit"), "") == "FALSE"
            for row in rows
        ),
        "n_direction_changed_after_top_five_primary_units_removal": sum(
            clean(row.get("candidate_axis_direction_preserved_without_top_five_primary_units"), "") == "FALSE"
            for row in rows
        ),
        "top_primary_unit_removal_status_counts": counter_text(
            row.get("top_primary_unit_removal_status") for row in rows
        ),
        "top_five_primary_units_removal_status_counts": counter_text(
            row.get("top_five_primary_units_removal_status") for row in rows
        ),
        "raw_residual_threshold_flag_counts": counter_text(
            row.get("raw_residual_threshold_flag") for row in rows
        ),
        "absolute_fit_status_legacy_counts": counter_text(
            row.get("absolute_fit_status_legacy") for row in rows
        ),
        "candidate_axis_vs_v6_3_direction_status_counts": counter_text(
            row.get("candidate_axis_vs_v6_3_direction_status") for row in rows
        ),
        "v6_3_baseline_origin_counts": counter_text(
            row.get("v6_3_baseline_origin") for row in rows
        ),
        "n_cells_v6_3_available": len(v6_available),
        "n_cells_v6_3_original_side": v6_directions["ORIGINAL_SIDE"],
        "n_cells_v6_3_proposal_side": v6_directions["PROPOSAL_SIDE"],
        "n_cells_v6_3_tied": v6_directions["TIE"],
        "n_cells_v6_3_unavailable": v6_directions["UNAVAILABLE"],
        "v6_3_proposed_relative_fit_score_q10_out_of_100": fmt(quantile(
            [row.get("proposed_swap_relative_fit_score_out_of_100") for row in v6_available], 0.10
        )),
        "v6_3_proposed_relative_fit_score_q25_out_of_100": fmt(quantile(
            [row.get("proposed_swap_relative_fit_score_out_of_100") for row in v6_available], 0.25
        )),
        "median_v6_3_proposed_relative_fit_score_out_of_100": fmt(quantile(
            [row.get("proposed_swap_relative_fit_score_out_of_100") for row in v6_available], 0.50
        )),
        "v6_3_proposed_relative_fit_score_q75_out_of_100": fmt(quantile(
            [row.get("proposed_swap_relative_fit_score_out_of_100") for row in v6_available], 0.75
        )),
        "v6_3_proposed_relative_fit_score_q90_out_of_100": fmt(quantile(
            [row.get("proposed_swap_relative_fit_score_out_of_100") for row in v6_available], 0.90
        )),
        "median_v6_3_original_relative_fit_score_out_of_100": fmt(quantile(
            [row.get("original_assignment_relative_fit_score_out_of_100") for row in v6_available], 0.50
        )),
        "eligible_n_cells_candidate_axis_available": len(eligible),
        "eligible_candidate_axis_position_raw_q10": fmt(quantile(eligible_positions, 0.10)),
        "eligible_candidate_axis_position_raw_q25": fmt(quantile(eligible_positions, 0.25)),
        "eligible_median_candidate_axis_position_raw": fmt(quantile(eligible_positions, 0.50)),
        "eligible_candidate_axis_position_raw_q75": fmt(quantile(eligible_positions, 0.75)),
        "eligible_candidate_axis_position_raw_q90": fmt(quantile(eligible_positions, 0.90)),
        "voting_semantics": (
            "CANDIDATE_AXIS_NEVER_VOTES;ONLY_EXISTING_REASSIGNMENT_SCOPES_"
            "DESCRIBE_PRIOR_AUTHORITATIVE_EVENT_VOTING"
        ),
    }
    if include_sampling_adjusted_fit:
        result.update({
            "median_candidate_a_excess_brier_mean": fmt(quantile(
                [row.get("candidate_a_excess_brier_mean") for row in rows], 0.5
            )),
            "median_candidate_b_excess_brier_mean": fmt(quantile(
                [row.get("candidate_b_excess_brier_mean") for row in rows], 0.5
            )),
        })
    return result


def make_summary_rows(cells, group_fields, include_sampling_adjusted_fit=False):
    grouped = defaultdict(list)
    for cell in cells:
        grouped[tuple(clean(cell.get(field)) for field in group_fields)].append(cell)
    rows = []
    for key in sorted(grouped):
        row = dict(zip(group_fields, key))
        row.update(summary_metrics(
            grouped[key], include_sampling_adjusted_fit
        ))
        row["schema_version"] = AGGREGATE_SCHEMA
        rows.append(row)
    return rows


def primary_event_summary(
    cells, group_fields, include_sampling_adjusted_fit=False
):
    grouped = defaultdict(list)
    for cell in cells:
        grouped[tuple(clean(cell.get(field)) for field in group_fields)].append(cell)
    rows = []
    for key in sorted(grouped):
        group = grouped[key]
        scopes = {clean(row.get("score_population_scope")) for row in group} & REASSIGNMENT_SCOPES
        row = dict(zip(group_fields, key))
        if len(scopes) == 1:
            scope = next(iter(scopes))
            primary = [item for item in group if clean(item.get("score_population_scope")) == scope]
            row.update(summary_metrics(primary, include_sampling_adjusted_fit))
            row["authoritative_summary_status"] = "AVAILABLE"
            row["authoritative_reassignment_scope"] = scope
        elif len(scopes) > 1:
            row.update({
                field: "NA" for field in (
                    SUMMARY_METRIC_FIELDS + (
                        SAMPLING_ADJUSTED_SUMMARY_FIELDS
                        if include_sampling_adjusted_fit else []
                    )
                )
            })
            row["n_cells_total"] = len(group)
            row["authoritative_summary_status"] = (
                "NO_AUTHORITATIVE_VERDICT_MIXED_APPLICATION_STATE"
            )
            row["authoritative_reassignment_scope"] = "MIXED"
        else:
            row.update({
                field: "NA" for field in (
                    SUMMARY_METRIC_FIELDS + (
                        SAMPLING_ADJUSTED_SUMMARY_FIELDS
                        if include_sampling_adjusted_fit else []
                    )
                )
            })
            row["n_cells_total"] = len(group)
            row["authoritative_summary_status"] = "NO_REASSIGNMENT_SCOPE"
            row["authoritative_reassignment_scope"] = "NONE"
        counts = Counter(clean(item.get("score_population_scope")) for item in group)
        for scope in POPULATIONS:
            row["n_" + scope.lower()] = counts[scope]
        row["schema_version"] = AGGREGATE_SCHEMA
        rows.append(row)
    return rows


def contrast_rows(cells):
    fixed_fields = [
        "row_type", "library", "selected_supported_event_id",
        "selected_supported_event_proposal",
        "original_allowed_demux_assignment", "n_applied",
        "n_retained_contrast", "median_applied_candidate_axis_position",
        "median_retained_contrast_candidate_axis_position",
        "applied_minus_retained_median_candidate_axis_position",
        "fraction_applied_proposal_side",
        "fraction_retained_contrast_original_side",
        "n_applied_legacy_qc_eligible", "n_retained_legacy_qc_eligible",
        "legacy_qc_eligible_median_applied_candidate_axis_position",
        "legacy_qc_eligible_median_retained_contrast_candidate_axis_position",
        "legacy_qc_eligible_applied_minus_retained_median_candidate_axis_position",
        "n_original_strata_with_applied", "n_original_strata_with_retained",
        "n_original_strata_with_both",
        "median_of_stratum_applied_minus_retained_differences",
        "pooled_contrast_interpretation", "schema_version",
    ]
    if not cells:
        return [], fixed_fields
    by_event = defaultdict(list)
    for cell in cells:
        event_key = (
            clean(cell.get("library")),
            clean(cell.get("selected_supported_event_id")),
            canonical_genotype(cell.get("selected_supported_event_proposal")),
        )
        by_event[event_key].append(cell)
    rows = []
    for (library, event_id, proposal), event_cells in sorted(by_event.items()):
        by_stratum = defaultdict(lambda: defaultdict(list))
        for cell in event_cells:
            scope = clean(cell.get("score_population_scope"))
            if scope in {
                "APPLIED_REASSIGNMENT", "RETAINED_ORIGINAL_CONTRAST_ONLY"
            }:
                stratum = clean(cell.get("original_allowed_demux_assignment"))
                by_stratum[stratum][scope].append(cell)
        differences = []
        for stratum in sorted(by_stratum):
            applied = by_stratum[stratum]["APPLIED_REASSIGNMENT"]
            retained = by_stratum[stratum]["RETAINED_ORIGINAL_CONTRAST_ONLY"]
            if not applied or not retained:
                continue
            applied_available = [row for row in applied if is_available(row)]
            retained_available = [row for row in retained if is_available(row)]
            applied_median = quantile(
                [row.get("candidate_axis_position_raw") for row in applied_available], 0.5
            )
            retained_median = quantile(
                [row.get("candidate_axis_position_raw") for row in retained_available], 0.5
            )
            difference = applied_median - retained_median if (
                math.isfinite(applied_median) and math.isfinite(retained_median)
            ) else math.nan
            if math.isfinite(difference):
                differences.append(difference)
            applied_eligible = [row for row in applied_available if is_legacy_eligible(row)]
            retained_eligible = [row for row in retained_available if is_legacy_eligible(row)]
            ae = quantile([row.get("candidate_axis_position_raw") for row in applied_eligible], 0.5)
            re = quantile([row.get("candidate_axis_position_raw") for row in retained_eligible], 0.5)
            rows.append({
                "row_type": "MATCHED_ORIGINAL_ASSIGNMENT_STRATUM",
                "library": library, "selected_supported_event_id": event_id,
                "selected_supported_event_proposal": proposal,
                "original_allowed_demux_assignment": stratum,
                "n_applied": len(applied), "n_retained_contrast": len(retained),
                "median_applied_candidate_axis_position": fmt(applied_median),
                "median_retained_contrast_candidate_axis_position": fmt(retained_median),
                "applied_minus_retained_median_candidate_axis_position": fmt(difference),
                "fraction_applied_proposal_side": fmt(
                    sum(clean(row.get("candidate_axis_direction")) == "PROPOSAL_SIDE" for row in applied_available) /
                    len(applied_available) if applied_available else math.nan
                ),
                "fraction_retained_contrast_original_side": fmt(
                    sum(clean(row.get("candidate_axis_direction")) == "ORIGINAL_SIDE" for row in retained_available) /
                    len(retained_available) if retained_available else math.nan
                ),
                "n_applied_legacy_qc_eligible": len(applied_eligible),
                "n_retained_legacy_qc_eligible": len(retained_eligible),
                "legacy_qc_eligible_median_applied_candidate_axis_position": fmt(ae),
                "legacy_qc_eligible_median_retained_contrast_candidate_axis_position": fmt(re),
                "legacy_qc_eligible_applied_minus_retained_median_candidate_axis_position": fmt(
                    ae - re if math.isfinite(ae) and math.isfinite(re) else math.nan
                ),
                "pooled_contrast_interpretation": "NOT_POOLED_MATCHED_STRATUM",
                "schema_version": AGGREGATE_SCHEMA,
            })
        strata_applied = {
            key for key, value in by_stratum.items()
            if value["APPLIED_REASSIGNMENT"]
        }
        strata_retained = {
            key for key, value in by_stratum.items()
            if value["RETAINED_ORIGINAL_CONTRAST_ONLY"]
        }
        rows.append({
            "row_type": "EVENT_LEVEL_MATCHED_STRATA_SUMMARY",
            "library": library, "selected_supported_event_id": event_id,
            "selected_supported_event_proposal": proposal,
            "original_allowed_demux_assignment": "ALL_MATCHED_STRATA",
            "n_original_strata_with_applied": len(strata_applied),
            "n_original_strata_with_retained": len(strata_retained),
            "n_original_strata_with_both": len(strata_applied & strata_retained),
            "median_of_stratum_applied_minus_retained_differences": fmt(
                quantile(differences, 0.5)
            ),
            "pooled_contrast_interpretation": "DESCRIPTIVE_COMPOSITION_CONFOUNDED",
            "schema_version": AGGREGATE_SCHEMA,
        })
    fields = list(fixed_fields)
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    return rows, fields


def normalize_library(value):
    text = clean(value, "").lower()
    if text.startswith("lib"):
        text = text[3:]
    return f"lib{int(text)}" if text.isdigit() else clean(value, "")


def require_pass_audit(path, expected_schema):
    rows = read_tsv(path)
    if not rows:
        raise ValueError(f"required audit is empty: {path}")
    if any(clean(row.get("schema_version"), "") != expected_schema for row in rows):
        raise ValueError(f"audit schema mismatch: {path}")
    if any(clean(row.get("status"), "") != "PASS" for row in rows):
        stops = [row for row in rows if clean(row.get("status"), "") != "PASS"]
        raise ValueError(f"audit contains STOP/non-PASS rows: {path}: {stops[:3]}")
    return rows


def data_row_count(path):
    with open_text(path, "rt") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def aggregate_main(argv):
    args = parse_args(argv)
    libraries = parse_library(args.libraries)
    event_id = clean(args.event_id, "")
    proposal = canonical_genotype(args.proposed_identity)
    if bool(event_id) != bool(proposal):
        raise ValueError("event ID and proposal must be supplied together")
    if event_id and len(libraries) != 1:
        raise ValueError("targeted candidate-axis aggregation requires one library")
    identity_root = Path(args.identity_root).resolve()
    pair_root = Path(args.pair_root).resolve()
    score_root = Path(args.score_root).resolve()
    output_root = Path(args.output_root).resolve()
    parent = output_root.parent
    parent.mkdir(parents=True, exist_ok=True)
    if output_root.exists():
        raise FileExistsError(
            f"candidate-axis aggregate destination already exists; refuses overwrite: {output_root}"
        )
    lock_path = parent / ".candidate_axis_aggregate.lock"
    lock_handle = open(lock_path, "a+")
    try:
        try:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(
                f"another candidate-axis aggregate holds the nonblocking lock: {lock_path}"
            ) from exc
        staging = Path(tempfile.mkdtemp(
            prefix=".candidate_axis_aggregate.staging.", dir=parent
        )).resolve()
        try:
            all_cells = []
            all_exclusions = []
            provenance_rows = []
            run_summary_rows = []
            raw_field_order = []
            for library in libraries:
                pair_manifest_path = pair_root / f"{library}.candidate_axis_pairs.tsv.gz"
                pair_summary_path = pair_root / f"{library}.candidate_axis_pair_summary.tsv"
                pair_exclusions_path = pair_root / f"{library}.candidate_axis_pair_exclusions.tsv.gz"
                pair_audit_path = pair_root / f"{library}.candidate_axis_pair_source_audit.tsv"
                raw_path = score_root / f"{library}.candidate_axis_scores.tsv.gz"
                score_provenance_path = score_root / f"{library}.candidate_axis_score_provenance.tsv"
                decisions_path = identity_root / "decisions" / f"{library}.reconciled_cells.tsv.gz"
                frozen_v3_path = identity_root / "nuclear" / f"{library}.identity_pair_probabilities.tsv.gz"
                frozen_v3_provenance_path = identity_root / "nuclear" / f"{library}.identity_pair_probability_provenance.tsv"
                for path in (
                    pair_manifest_path, pair_summary_path, pair_exclusions_path,
                    pair_audit_path, raw_path, score_provenance_path,
                    decisions_path,
                ):
                    if not path.is_file():
                        raise FileNotFoundError(
                            f"{library}: required candidate-axis aggregate input is missing: {path}"
                        )

                require_pass_audit(pair_audit_path, PAIR_AUDIT_SCHEMA)
                pair_summary_rows = read_tsv(pair_summary_path)
                if len(pair_summary_rows) != 1 or clean(
                    pair_summary_rows[0].get("schema_version"), ""
                ) != PAIR_SCHEMA:
                    raise ValueError(f"{library}: pair summary must be one V1 row")
                pair_summary = pair_summary_rows[0]
                if clean(pair_summary.get("library"), "") != library:
                    raise ValueError(f"{library}: pair summary library mismatch")
                if event_id and (
                    clean(pair_summary.get("selected_supported_event_id"), "") != event_id
                    or canonical_genotype(
                        pair_summary.get("selected_supported_event_proposal")
                    ) != proposal
                ):
                    raise ValueError(f"{library}: targeted pair summary event mismatch")

                score_provenance_rows = read_tsv(score_provenance_path)
                if len(score_provenance_rows) != 1 or clean(
                    score_provenance_rows[0].get("schema_version"), ""
                ) != SCORE_PROVENANCE_SCHEMA:
                    raise ValueError(f"{library}: scorer provenance must be one V1 row")
                score_provenance = score_provenance_rows[0]
                if normalize_library(score_provenance.get("library")) != library:
                    raise ValueError(f"{library}: scorer provenance library mismatch")

                manifest_pairs, manifest_rows = load_manifest_pairs(
                    pair_manifest_path, library, event_id, proposal
                )
                if int(clean(score_provenance.get("pair_manifest_row_count"), "-1")) != len(manifest_rows):
                    raise ValueError(f"{library}: scorer provenance manifest-row count mismatch")
                if int(clean(score_provenance.get("pair_summary_row_count"), "-1")) != 1:
                    raise ValueError(f"{library}: scorer provenance pair-summary row count mismatch")
                if int(clean(score_provenance.get("pair_summary_n_score_pairs"), "-1")) != len(manifest_pairs):
                    raise ValueError(f"{library}: scorer provenance pair-count mismatch")
                for field in read_tsv_fields(raw_path):
                    if field not in raw_field_order:
                        raw_field_order.append(field)
                raw_rows = read_tsv(raw_path)
                raw_by_key = unique_rows(
                    raw_rows, ["library", "barcode", "score_pair_id"],
                    f"{library} candidate-axis raw scores",
                )
                if set(raw_by_key) != set(manifest_pairs):
                    raise ValueError(
                        f"{library}: candidate-axis manifest/raw one-to-one accounting failed; "
                        f"missing={sorted(set(manifest_pairs)-set(raw_by_key))[:5]} "
                        f"extra={sorted(set(raw_by_key)-set(manifest_pairs))[:5]}"
                    )
                if int(clean(score_provenance.get("raw_output_row_count"), "-1")) != len(raw_rows):
                    raise ValueError(f"{library}: scorer provenance raw-row count mismatch")

                decision_rows = read_tsv(decisions_path)
                decisions = {}
                for row in decision_rows:
                    if normalize_library(row.get("library")) != library:
                        raise ValueError(f"{library}: decision table contains a cross-library row")
                    barcode = clean(row.get("barcode"), "")
                    if not barcode or barcode in decisions:
                        raise ValueError(f"{library}: decision table has blank/duplicate barcode: {barcode!r}")
                    decisions[barcode] = row

                module = None
                frozen_v3 = {}
                v6_source = None
                error_ref = error_alt = 0.001
                if manifest_pairs:
                    v6_text = clean(pair_summary.get("v6_3_aggregator_path"), "")
                    if not v6_text:
                        raise ValueError(f"{library}: pair-bearing summary lacks V6.3 source")
                    v6_source = Path(v6_text).resolve()
                    module = import_v6(v6_source)
                    if frozen_v3_path.is_file():
                        frozen_v3_rows = read_tsv(frozen_v3_path)
                        for row in frozen_v3_rows:
                            if normalize_library(row.get("library")) != library:
                                raise ValueError(
                                    f"{library}: frozen V3 table contains a cross-library row"
                                )
                            barcode = clean(row.get("barcode"), "")
                            if not barcode or barcode in frozen_v3:
                                raise ValueError(
                                    f"{library}: frozen V3 table has blank/duplicate barcode: {barcode!r}"
                                )
                            frozen_v3[barcode] = row
                        error_ref, error_alt = frozen_v3_error_pair(
                            frozen_v3_rows
                        )

                frozen_cell, frozen_review, frozen_v6_records = load_optional_frozen_v6(
                    args.v6_3_root, library
                )
                frozen_v6_by_barcode = {}
                for source in (frozen_cell, frozen_review):
                    for (_, barcode), row in source.items():
                        if barcode in frozen_v6_by_barcode:
                            raise ValueError(f"{library}: duplicate frozen V6.3 barcode: {barcode}")
                        frozen_v6_by_barcode[barcode] = row

                ambient_by_barcode = {}
                ambient_record = "NA"
                if args.ambient_path:
                    ambient_record = file_identity(args.ambient_path)
                    for row in read_tsv(args.ambient_path):
                        if args.ambient_condition and clean(row.get("condition"), "") != args.ambient_condition:
                            continue
                        if normalize_library(row.get("library")) != library:
                            continue
                        barcode = clean(row.get("barcode"), "")
                        if not barcode or barcode in ambient_by_barcode:
                            raise ValueError(f"{library}: ambient annotations have blank/duplicate key: {barcode!r}")
                        ambient_by_barcode[barcode] = row

                library_cells = []
                for key in sorted(manifest_pairs):
                    barcode = key[1]
                    manifest = manifest_pairs[key]
                    retained = manifest["retained"]
                    decision = decisions.get(barcode)
                    if decision is None:
                        raise ValueError(f"{library}: manifest row is missing from the decision table: {key}")
                    if retained:
                        if clean(decision.get("final_action"), "").upper() != "KEEP" or \
                                clean(decision.get("reassignment_applied"), "").upper() != "FALSE" or \
                                canonical_genotype(decision.get("reconciled_donor_genotype", "")) != \
                                canonical_genotype(decision.get("original_demux_assignment", "")):
                            raise ValueError(f"{library}: retained decision no longer satisfies KEEP/original contract: {key}")
                    else:
                        if clean(decision.get("event_id"), "") != manifest["event_id"] or \
                                canonical_genotype(decision.get("proposed_donor_genotype", "")) != manifest["proposal"]:
                            raise ValueError(f"{library}: operational decision no longer matches its pair event: {key}")
                    v3 = None if retained else frozen_v3.get(barcode)
                    frozen = None if retained else frozen_v6_by_barcode.get(barcode)
                    cell = build_cell(
                        module, library, key, manifest, raw_by_key[key],
                        decision, v3, frozen,
                    )
                    if barcode in ambient_by_barcode:
                        annotation = ambient_by_barcode[barcode]
                        for field in (
                            "ambient_sensitivity_status", "ambient_original_c",
                            "ambient_reconciled_c", "ambient_reconciled_minus_original_c",
                            "ambient_stratum", "ambient_assignment_diagnostic",
                        ):
                            if field in annotation:
                                cell[field] = clean(annotation.get(field))
                    library_cells.append(cell)

                exclusion_rows = read_tsv(pair_exclusions_path)
                if any(clean(row.get("exclusion_stage"), "") != "PAIR_CONSTRUCTION" for row in exclusion_rows):
                    raise ValueError(f"{library}: exclusions may only originate during pair construction")
                if int(float(clean(pair_summary.get("n_manifest_rows"), "0"))) != len(manifest_rows):
                    raise ValueError(f"{library}: pair summary manifest-row count mismatch")
                if int(float(clean(pair_summary.get("n_score_pairs"), "0"))) != len(library_cells):
                    raise ValueError(f"{library}: pair summary pair count mismatch")
                if int(float(clean(pair_summary.get("n_pair_construction_exclusions"), "0"))) != len(exclusion_rows):
                    raise ValueError(f"{library}: pair summary exclusion count mismatch")

                all_cells.extend(library_cells)
                all_exclusions.extend(exclusion_rows)
                pair_count = len(library_cells)
                provenance_rows.append({
                    "library": library,
                    "selected_event_key": clean(pair_summary.get("supported_event_key")),
                    "pair_source_audit_path": str(pair_audit_path),
                    "pair_manifest_path": str(pair_manifest_path),
                    "pair_summary_path": str(pair_summary_path),
                    "candidate_axis_raw_score_path": str(raw_path),
                    "candidate_axis_score_provenance_path": str(score_provenance_path),
                    "decision_path": str(decisions_path),
                    "validation_path": clean(pair_summary.get("validation_summary")),
                    "samples_path": clean(pair_summary.get("samples_path")),
                    "pileup_sites_path": clean(pair_summary.get("pileup_sites_path")),
                    "pileup_observations_path": clean(pair_summary.get("pileup_observations_path")),
                    "frozen_v3_probability_path": (
                        str(frozen_v3_path)
                        if pair_count and frozen_v3_path.is_file() else "NA"
                    ),
                    "frozen_v3_provenance_path": (
                        str(frozen_v3_provenance_path)
                        if pair_count and frozen_v3_provenance_path.is_file() else "NA"
                    ),
                    "frozen_v6_3_cell_path": frozen_v6_records["cell"]["path"] if isinstance(frozen_v6_records["cell"], dict) else "NA",
                    "frozen_v6_3_review_only_path": frozen_v6_records["review"]["path"] if isinstance(frozen_v6_records["review"], dict) else "NA",
                    "ambient_path": ambient_record["path"] if isinstance(ambient_record, dict) else "NA",
                    "scorer_binary_path": clean(score_provenance.get("scorer_binary_path")),
                    "pair_builder_path": clean(pair_summary.get("pair_builder_path")),
                    "candidate_axis_aggregator_path": str(Path(__file__).resolve()),
                    "v6_3_aggregator_path": str(v6_source) if v6_source else "NA",
                    "candidate_axis_formula_version": clean(library_cells[0].get("candidate_axis_formula_version")) if library_cells else "WEIGHTED_BRIER_FIXED_PAIR_PROJECTION_UNCLIPPED_V1",
                    "candidate_prediction_transform_version": clean(library_cells[0].get("candidate_prediction_transform_version")) if library_cells else "HARD_GT_DOSAGE_EXISTING_ASYMMETRIC_ERROR_TRANSFORM_V1",
                    "primary_evidence_basis": clean(library_cells[0].get("candidate_axis_evidence_basis")) if library_cells else "NUCLEAR_SITE_BALANCED_FIXED_PRIMARY",
                    "error_ref": (
                        format(error_ref, ".17g") if pair_count else "NA"
                    ),
                    "error_alt": (
                        format(error_alt, ".17g") if pair_count else "NA"
                    ),
                    "fold_definition_version": clean(library_cells[0].get("candidate_axis_fold_definition_version")) if library_cells else "CANDIDATE_AXIS_GREEDY_DESIGN_MASS_SITE_GROUPS_PROJECT_FNV1A64_COMPAT_V1",
                    "numerical_tolerance_version": clean(library_cells[0].get("numerical_tolerance_version")) if library_cells else "IEEE754_LONG_DOUBLE_SCALE64_V1",
                    "bucket_count": clean(library_cells[0].get("candidate_axis_bucket_count")) if library_cells else "0",
                    "warnings": "NONE" if pair_count else "NO_FINALIZED_SUPPORTED_EVENTS",
                    "schema_version": AGGREGATE_SCHEMA,
                })
                population_counts = Counter(
                    clean(cell.get("score_population_scope"))
                    for cell in library_cells
                )
                zero_event = (
                    pair_count == 0 and
                    clean(pair_summary.get("selected_supported_event_id"), "") == "NONE"
                )
                run_summary_rows.append({
                    "library": library,
                    "selected_supported_event_id": clean(pair_summary.get("selected_supported_event_id")),
                    "selected_supported_event_proposal": clean(pair_summary.get("selected_supported_event_proposal")),
                    "n_manifest_pairs": len(manifest_pairs),
                    "n_raw_score_rows": len(raw_rows), "n_cell_rows": len(library_cells),
                    "n_pair_construction_exclusions": len(exclusion_rows),
                    "population_scope_counts": counter_text(
                        cell.get("score_population_scope") for cell in library_cells
                    ),
                    "n_applied_reassignment": population_counts["APPLIED_REASSIGNMENT"],
                    "n_recommended_not_applied": population_counts["RECOMMENDED_NOT_APPLIED"],
                    "n_supported_event_held": population_counts["SUPPORTED_EVENT_HELD_CELL"],
                    "n_review_only": population_counts["REVIEW_ONLY_UNEXPECTED_IDENTITY"],
                    "n_retained_original_contrast": population_counts["RETAINED_ORIGINAL_CONTRAST_ONLY"],
                    "candidate_axis_status_counts": counter_text(
                        cell.get("candidate_axis_status") for cell in library_cells
                    ),
                    "accounting_status": "PASS",
                    "scientific_pilot_status": (
                        "NO_FINALIZED_SUPPORTED_EVENTS" if zero_event
                        else "READY_FOR_SCIENTIFIC_REVIEW"
                    ),
                    "assignment_outputs_changed": "FALSE", "plots_created": "FALSE",
                    "resamples_run": "0", "schema_version": AGGREGATE_SCHEMA,
                })

            library_rank = {library: index for index, library in enumerate(libraries)}
            all_cells.sort(key=lambda row: (
                library_rank[clean(row.get("library"))],
                clean(row.get("selected_supported_event_id")),
                canonical_genotype(row.get("selected_supported_event_proposal")),
                clean(row.get("barcode")), clean(row.get("score_pair_id")),
            ))
            all_exclusions.sort(key=lambda row: (
                library_rank[clean(row.get("library"))],
                clean(row.get("selected_supported_event_id")),
                clean(row.get("barcode")), clean(row.get("exclusion_reason")),
            ))
            cell_fields = ["schema_version"]
            for field in raw_field_order:
                if field not in cell_fields:
                    cell_fields.append(field)
            for row in all_cells:
                for field in row:
                    if field not in cell_fields:
                        cell_fields.append(field)
            cell_path = staging / "candidate_axis_cell_scores.tsv.gz"
            write_tsv(cell_path, all_cells, cell_fields)

            event_fields = [
                "library", "selected_supported_event_id",
                "selected_supported_event_proposal",
            ]
            strata_fields = event_fields + ["original_allowed_demux_assignment"]
            scope_fields = event_fields + ["score_population_scope"]
            event_rows = primary_event_summary(
                all_cells, event_fields, include_sampling_adjusted_fit=True
            )
            strata_rows = primary_event_summary(all_cells, strata_fields)
            scope_rows = make_summary_rows(
                all_cells, scope_fields, include_sampling_adjusted_fit=True
            )
            event_summary_path = staging / "candidate_axis_event_summary.tsv"
            strata_summary_path = staging / "candidate_axis_event_strata_summary.tsv"
            scope_summary_path = staging / "candidate_axis_event_scope_summary.tsv"
            event_summary_fields = event_fields + [
                "authoritative_summary_status", "authoritative_reassignment_scope",
                *["n_" + scope.lower() for scope in POPULATIONS],
                *SUMMARY_METRIC_FIELDS, *SAMPLING_ADJUSTED_SUMMARY_FIELDS,
                "schema_version",
            ]
            strata_summary_fields = strata_fields + [
                "authoritative_summary_status", "authoritative_reassignment_scope",
                *["n_" + scope.lower() for scope in POPULATIONS],
                *SUMMARY_METRIC_FIELDS, "schema_version",
            ]
            scope_summary_fields = (
                scope_fields + SUMMARY_METRIC_FIELDS +
                SAMPLING_ADJUSTED_SUMMARY_FIELDS + ["schema_version"]
            )
            write_tsv(event_summary_path, event_rows, event_summary_fields)
            write_tsv(strata_summary_path, strata_rows, strata_summary_fields)
            write_tsv(scope_summary_path, scope_rows, scope_summary_fields)

            contrast, contrast_fields = contrast_rows(all_cells)
            contrast_path = staging / "candidate_axis_applied_vs_retained_contrast.tsv"
            write_tsv(contrast_path, contrast, contrast_fields)
            exclusions_path = staging / "candidate_axis_pair_exclusions.tsv.gz"
            exclusion_fields = (
                list(all_exclusions[0]) if all_exclusions
                else PAIR_EXCLUSION_FIELDS
            )
            write_tsv(exclusions_path, all_exclusions, exclusion_fields)

            provenance_path = staging / "candidate_axis_provenance.tsv"
            write_tsv(provenance_path, provenance_rows, list(provenance_rows[0]))

            run_summary_path = staging / "candidate_axis_run_summary.tsv"
            write_tsv(run_summary_path, run_summary_rows, list(run_summary_rows[0]))

            outputs = [
                (cell_path, CELL_SCHEMA, "one row per admitted fixed candidate pair"),
                (event_summary_path, AGGREGATE_SCHEMA, "selected event authoritative reassignment scope"),
                (strata_summary_path, AGGREGATE_SCHEMA, "selected event by original-assignment stratum"),
                (scope_summary_path, AGGREGATE_SCHEMA, "selected event by population scope"),
                (contrast_path, AGGREGATE_SCHEMA, "matched applied-versus-retained stratum contrast"),
                (exclusions_path, PAIR_SCHEMA, "pre-manifest pair-construction exclusion"),
                (provenance_path, AGGREGATE_SCHEMA, "one aggregate provenance row"),
                (run_summary_path, AGGREGATE_SCHEMA, "one run accounting row"),
            ]
            manifest_output_rows = []
            for path, schema, grain in outputs:
                manifest_output_rows.append({
                    "relative_path": path.name, "schema_version": schema,
                    "row_count": data_row_count(path),
                    "semantic_grain": grain,
                })
            output_manifest_path = staging / "candidate_axis_output_manifest.tsv"
            write_tsv(
                output_manifest_path, manifest_output_rows,
                ["relative_path", "schema_version", "row_count", "semantic_grain"],
            )
            if output_root.exists():
                raise FileExistsError(f"candidate-axis aggregate destination appeared during staging: {output_root}")
            os.rename(staging, output_root)
            print(
                f"PASS: candidate-axis aggregate wrote {len(all_cells)} cell rows "
                f"for {len(libraries)} requested libraries"
            )
            return 0
        except Exception:
            if staging.exists() and staging.parent == parent:
                shutil.rmtree(staging)
            raise
    finally:
        try:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
        finally:
            lock_handle.close()


def self_test():
    def require(condition, message):
        if not condition:
            raise AssertionError(message)

    require(quantile([0, 100], 0.5) == 50, "linear quantile invariant")
    require(compare_directions("PROPOSAL_SIDE", "PROPOSAL_SIDE") == "CONCORDANT",
            "direction concordance invariant")
    require(compare_directions("ORIGINAL_SIDE", "PROPOSAL_SIDE") == "DISCORDANT",
            "direction discordance invariant")
    try:
        unique_rows([{"k": "x"}, {"k": "x"}], ["k"], "synthetic")
    except ValueError:
        pass
    else:
        raise AssertionError("duplicate key rejection invariant")

    original = {
        "library": "lib19", "barcode": "BC", "score_pair_id": "PAIR",
        "supported_event_key": "lib19|EVENT|A+B", "donor_genotype": "A",
        "score_pair_role": "ORIGINAL_ALLOWED_DEMUX", "candidate_origin": "ORIGINAL",
        "score_pair_source": "IDENTITY_RECONCILIATION", "score_population_scope": "APPLIED_REASSIGNMENT",
        "population_votes_in_authoritative_event": "TRUE", "score_scope_contract": OPERATIONAL_CONTRACT,
    }
    proposed = dict(original, donor_genotype="A+B",
                    score_pair_role="RECONCILIATION_NOMINATED_SWAP",
                    candidate_origin="RECONCILIATION_NOMINATION")
    raw = {
        "schema_version": RAW_SCHEMA, "library": "lib19", "barcode": "BC",
        "score_pair_id": "PAIR", "supported_event_key": "lib19|EVENT|A+B",
        "candidate_a": "A", "candidate_b": "A+B",
        "candidate_a_role": "ORIGINAL_ALLOWED_DEMUX",
        "candidate_b_role": "RECONCILIATION_NOMINATED_SWAP",
        "candidate_a_origin": "ORIGINAL", "candidate_b_origin": "RECONCILIATION_NOMINATION",
        "score_pair_source": "IDENTITY_RECONCILIATION",
        "score_population_scope": "APPLIED_REASSIGNMENT",
        "population_votes_in_authoritative_event": "TRUE",
        "score_scope_contract": OPERATIONAL_CONTRACT,
        "n_candidate_axis_discriminating_sites": "1",
        "candidate_a_observed_brier_mean": "0.04",
        "candidate_a_expected_sampling_brier_mean": "0.02",
        "candidate_a_excess_brier_mean": "0.02",
        "candidate_b_observed_brier_mean": "0.16",
        "candidate_b_expected_sampling_brier_mean": "0.03",
        "candidate_b_excess_brier_mean": "0.13",
        "candidate_a_residual_mismatch_legacy": "0.2",
        "candidate_b_residual_mismatch_legacy": "0.4",
        "absolute_fit_status_legacy": "PASS",
        "comparison_status_legacy": "PASS",
        "poor_fit_residual_threshold_legacy": "0.3",
        "raw_residual_threshold_flag":
            "CANDIDATE_B_ONLY_ABOVE_LEGACY_THRESHOLD",
    }
    validate_raw_join({"original": original, "proposed": proposed}, raw,
                      ("lib19", "BC", "PAIR"))
    bad_raw = dict(raw, candidate_b="THIRD")
    try:
        validate_raw_join({"original": original, "proposed": proposed}, bad_raw,
                          ("lib19", "BC", "PAIR"))
    except ValueError:
        pass
    else:
        raise AssertionError("candidate orientation join invariant")

    synthetic_cells = [
        {
            "library": "lib19", "selected_supported_event_id": "EVENT",
            "selected_supported_event_proposal": "A+B",
            "original_allowed_demux_assignment": "A",
            "score_population_scope": "APPLIED_REASSIGNMENT",
            "candidate_axis_status": "AVAILABLE", "candidate_axis_position_raw": "75",
            "candidate_axis_direction": "PROPOSAL_SIDE", "legacy_pair_qc_status": "ELIGIBLE",
            "population_votes_in_authoritative_event": "TRUE",
        },
        {
            "library": "lib19", "selected_supported_event_id": "EVENT",
            "selected_supported_event_proposal": "A+B",
            "original_allowed_demux_assignment": "A",
            "score_population_scope": "RETAINED_ORIGINAL_CONTRAST_ONLY",
            "candidate_axis_status": "AVAILABLE", "candidate_axis_position_raw": "25",
            "candidate_axis_direction": "ORIGINAL_SIDE", "legacy_pair_qc_status": "ELIGIBLE",
            "population_votes_in_authoritative_event": "FALSE",
        },
        {
            "library": "lib20", "selected_supported_event_id": "EVENT_2",
            "selected_supported_event_proposal": "C+D",
            "original_allowed_demux_assignment": "C",
            "score_population_scope": "SUPPORTED_EVENT_HELD_CELL",
            "candidate_axis_status": "AVAILABLE", "candidate_axis_position_raw": "60",
            "candidate_axis_direction": "PROPOSAL_SIDE", "legacy_pair_qc_status": "ELIGIBLE",
            "population_votes_in_authoritative_event": "FALSE",
        },
    ]
    scope = make_summary_rows(
        synthetic_cells,
        ["library", "selected_supported_event_id", "selected_supported_event_proposal",
         "score_population_scope"],
    )
    require(sum(int(row["n_cells_total"]) for row in scope) == len(synthetic_cells),
            "scope accounting invariant")
    event_summaries = make_summary_rows(
        synthetic_cells,
        ["library", "selected_supported_event_id", "selected_supported_event_proposal"],
    )
    require(len(event_summaries) == 2,
            "multi-library and multi-event grouping invariant")
    sampling_cells = [dict(
        synthetic_cells[0],
        candidate_a_excess_brier_mean="0.01",
        candidate_b_excess_brier_mean="0.03",
        raw_residual_threshold_flag="NEITHER_ABOVE_LEGACY_THRESHOLD",
        absolute_fit_status_legacy="PASS",
    )]
    sampling_summary = make_summary_rows(
        sampling_cells,
        ["library", "selected_supported_event_id", "score_population_scope"],
        include_sampling_adjusted_fit=True,
    )
    require(
        sampling_summary[0]["median_candidate_a_excess_brier_mean"] == "0.01"
        and sampling_summary[0]["median_candidate_b_excess_brier_mean"] == "0.03"
        and sampling_summary[0]["raw_residual_threshold_flag_counts"] ==
            "NEITHER_ABOVE_LEGACY_THRESHOLD:1",
        "sampling-adjusted fit summary invariant",
    )
    contrasts, _ = contrast_rows(synthetic_cells)
    require({row["library"] for row in contrasts} == {"lib19", "lib20"},
            "contrasts remain library/event local")
    empty_contrast, empty_fields = contrast_rows([])
    require(not empty_contrast and "schema_version" in empty_fields,
            "zero-pair header-only contrast invariant")
    require(parse_library(["20", "19", "19"]) == ["lib19", "lib20"],
            "library parsing is sorted and unique")

    with tempfile.TemporaryDirectory(prefix="candidate_axis_aggregate_selftest.") as directory:
        target = Path(directory) / "atomic.tsv.gz"
        write_tsv(target, [{"value": "1", "schema_version": "synthetic"}],
                  ["value", "schema_version"])
        require(data_row_count(target) == 1,
                "atomic output accounting invariant")
        partial = Path(directory) / "partial.final"
        staging = Path(directory) / ".staging"
        staging.mkdir()
        try:
            raise RuntimeError("injected")
        except RuntimeError:
            shutil.rmtree(staging)
        require(not partial.exists() and not staging.exists(),
                "injected-failure cleanup invariant")
        empty_manifest = Path(directory) / "empty_pairs.tsv.gz"
        write_tsv(empty_manifest, [], sorted({
            "schema_version", "library", "barcode", "score_pair_id",
            "score_pair_role", "score_pair_source", "score_population_scope",
            "population_votes_in_authoritative_event", "supported_event_key",
            "selected_supported_event_id", "selected_supported_event_proposal",
            "reconciliation_event_id", "original_demux_assignment",
            "proposed_donor_genotype", "candidate_b_fixed_identity",
            "score_scope_contract", "candidate_origin", "donor_genotype",
        }))
        empty_pairs, empty_rows = load_manifest_pairs(
            empty_manifest, "lib20"
        )
        require(not empty_pairs and not empty_rows,
                "zero-pair manifest acceptance invariant")
    print("PASS: identity candidate-axis aggregate self-test")
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv == ["--self-test"]:
        return self_test()
    if argv and argv[0] == "score-provenance":
        return score_provenance_main(argv[1:])
    return aggregate_main(argv)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
