#!/usr/bin/env python3
"""Create a self-contained audit report for a Tetraploid Arm CNV run."""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import os
import sys
from collections import Counter
from pathlib import Path

from tetra_arm_common import (
    atomic_text,
    clean,
    file_record,
    finite_float,
    natural_key,
    open_text,
    require_file,
    require_outputs_absent,
    write_json_atomic,
    write_tsv_atomic,
)


PROGRAM_VERSION = "2.6.0"
REPORT_SCHEMA = "tetra_arm_cnv_report_v3"
CALL_SCHEMA = "tetra_arm_cnv_calls_v2"
CALIBRATION_SCHEMA = "tetra_arm_calibration_v2"
UID_SCHEMA = "tetra_arm_uid_chromosome_flags_v3"
PAIR_SCHEMA = "tetra_arm_donor_pair_arm_summary_v3"
CALL_QC_SCHEMA = "tetra_arm_call_qc_v2"
CALL_CONTRACT_SCHEMA = "tetra_arm_call_contract_v2"
EVENT_STATES = (
    "DONOR_A_LOSS", "DONOR_B_LOSS", "DONOR_A_GAIN", "DONOR_B_GAIN")
TERMINAL_STATES = {
    "PASS_NO_HETEROTYPIC_TARGETS",
    "PASS_NO_OBSERVED_ASE",
    "PASS_NO_CALLABLE_ASE",
}


def first(row: dict[str, str], *names: str, default: str = "") -> str:
    for name in names:
        value = clean(row.get(name))
        if value:
            return value
    return default


def is_autosomal_chromosome(value: object) -> bool:
    chromosome = clean(value).upper().removeprefix("CHR")
    return chromosome.isdigit() and 1 <= int(chromosome) <= 22


def selected_state_q(row: dict[str, str], state: str) -> float:
    """Use the selected event state's multiplicity-adjusted q when available."""
    if state in EVENT_STATES:
        value = finite_float(row.get(f"empirical_q_{state}"))
        if math.isfinite(value):
            return value
    return finite_float(first(row, "empirical_q_value", "q_value", "fdr"))


def parse_qc(path: str) -> list[tuple[str, str]]:
    """Accept both metric/value QC and one-row headered QC tables."""
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = [row for row in reader if row]
    if not rows:
        raise ValueError(f"empty QC file: {path}")
    if len(rows[0]) == 2 and [cell.lower() for cell in rows[0]] == ["metric", "value"]:
        values = []
        seen = set()
        for line_number, row in enumerate(rows[1:], start=2):
            metric = clean(row[0]) if len(row) == 2 else ""
            if len(row) != 2 or not metric or metric in seen:
                raise ValueError(
                    f"malformed or duplicate QC metric at {path}:{line_number}")
            seen.add(metric)
            values.append((metric, row[1]))
        return values
    header = rows[0]
    if (len(rows) != 2 or len(header) != len(set(header)) or
            any(not clean(value) for value in header) or
            len(rows[1]) != len(header)):
        raise ValueError(f"unrecognized QC table: {path}")
    return list(zip(header, rows[1]))


def small_table(headers: list[str], rows: list[list[object]], limit: int | None = None) -> str:
    if limit is not None:
        rows = rows[:limit]
    head = "".join(f"<th>{html.escape(str(value))}</th>" for value in headers)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(
            f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def bar_table(counter: Counter, total: int, limit: int = 30) -> str:
    rows = []
    for label, count in sorted(counter.items(), key=lambda item: (-item[1], natural_key(item[0])))[:limit]:
        fraction = count / total if total else 0.0
        bar = (
            '<div class="bar-track"><div class="bar-fill" style="width:'
            f'{100.0 * fraction:.3f}%"></div></div>')
        rows.append(
            "<tr><td>" + html.escape(str(label)) + "</td><td>" + str(count) +
            "</td><td>" + f"{100.0 * fraction:.2f}%" + "</td><td>" + bar + "</td></tr>")
    return (
        "<table><thead><tr><th>Category</th><th>Rows</th><th>Fraction</th>"
        "<th>Distribution</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>")


def stream_calls(path: str, top_limit: int) -> dict:
    totals = Counter()
    states = Counter()
    best_states = Counter()
    statuses = Counter()
    evidence = Counter()
    libraries = Counter()
    arms = Counter()
    donor_pairs = Counter()
    cells = set()
    pass_event_cells = set()
    top: list[tuple[float, dict[str, str]]] = []
    header: list[str] = []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        if not header or not {"barcode", "arm"} <= set(header):
            raise ValueError(f"call table lacks barcode/arm columns: {path}")
        required = {
            "library", "barcode", "arm", "chromosome", "best_state",
            "best_state_posterior", "empirical_q_resolution_floor",
            "call_state", "call_status", "call_schema_version",
        }
        missing = sorted(required - set(header))
        if missing:
            raise ValueError(f"call table lacks v2 columns {missing}: {path}")
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"malformed call row {path}:{line_number}")
            if clean(row.get("call_schema_version")) != CALL_SCHEMA:
                raise ValueError(
                    f"incompatible call schema at {path}:{line_number}; "
                    f"expected {CALL_SCHEMA}")
            library = first(row, "library", default="NA")
            barcode = first(row, "barcode", default="NA")
            arm = first(row, "arm", default="NA")
            state = first(row, "best_state", "cnv_state", "call", default="UNKNOWN")
            call_state = first(row, "call_state", default="NO_CALL")
            status = first(row, "call_status", "model_status", "status", default="UNKNOWN")
            basis = first(row, "evidence_status", "evidence_basis", default="UNKNOWN")
            pair = first(row, "donor_pair", default="NA")
            qvalue = selected_state_q(row, state)
            q_resolution = finite_float(row.get("empirical_q_resolution_floor"))
            posterior = finite_float(first(row, "best_state_posterior", "posterior", "call_probability"))
            autosomal = is_autosomal_chromosome(row.get("chromosome"))
            totals["rows"] += 1
            cells.add((library, barcode))
            states[call_state] += 1
            best_states[state] += 1
            statuses[status] += 1
            evidence[basis] += 1
            libraries[library] += 1
            arms[arm] += 1
            donor_pairs[pair] += 1
            if not autosomal:
                totals["exploratory_non_autosomal_rows"] += 1
            is_balanced = state.upper() in {"BALANCED", "NO_CALL", "UNKNOWN"}
            is_pass_event = status == "PASS_EVENT" and call_state not in {
                "BALANCED", "NO_CALL", "UNKNOWN"}
            is_q05_candidate = math.isfinite(qvalue) and qvalue <= 0.05 and not is_balanced
            if is_q05_candidate:
                totals["nonbalanced_candidate_rows_q05"] += 1
            if is_pass_event:
                totals["pass_event_rows"] += 1
                pass_event_cells.add((library, barcode))
            if not is_balanced:
                totals["nonbalanced_best_state_rows"] += 1
            rank = ((-math.log10(max(qvalue, 1e-300))) if math.isfinite(qvalue)
                    else (posterior if math.isfinite(posterior) else -math.inf))
            if not is_balanced and math.isfinite(rank):
                selected = {
                    "library": library,
                    "barcode": barcode,
                    "arm": arm,
                    "donor_pair": pair,
                    "best_state": state,
                    "posterior": "NA" if not math.isfinite(posterior) else f"{posterior:.4g}",
                    "q_value": "NA" if not math.isfinite(qvalue) else f"{qvalue:.4g}",
                    "q_resolution": (
                        "NA" if not math.isfinite(q_resolution)
                        else f"{q_resolution:.4g}"),
                    "inference_scope": (
                        "AUTOSOMAL" if autosomal else "EXPLORATORY_NON_AUTOSOMAL"),
                    "status": status,
                }
                top.append((rank, selected))
                if len(top) > 4 * top_limit:
                    top = sorted(top, key=lambda item: item[0], reverse=True)[:top_limit]
    top = [row for _rank, row in sorted(top, key=lambda item: item[0], reverse=True)[:top_limit]]
    totals["cells"] = len(cells)
    totals["pass_event_cells"] = len(pass_event_cells)
    return {
        "header": header,
        "totals": dict(totals),
        "states": states,
        "best_states": best_states,
        "statuses": statuses,
        "evidence": evidence,
        "libraries": libraries,
        "arms": arms,
        "donor_pairs": donor_pairs,
        "top": top,
    }


def load_optional_rows(
        path: str, maximum: int = 10000, expected_schema: str = "",
        match_field: str = "", match_value: str = "",
        ) -> tuple[int, list[str], list[dict[str, str]], int]:
    if not path:
        return 0, [], [], 0
    target = require_file(path)
    count = 0
    matching = 0
    retained: list[dict[str, str]] = []
    with open_text(target) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        if not header or len(header) != len(set(header)):
            raise ValueError(f"headerless or duplicate-column table: {target}")
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"malformed row {target}:{line_number}")
            if (expected_schema and
                    clean(row.get("schema_version")) != expected_schema):
                raise ValueError(
                    f"incompatible schema at {target}:{line_number}; "
                    f"expected {expected_schema}")
            count += 1
            if match_field and clean(row.get(match_field)).upper() == match_value.upper():
                matching += 1
            if len(retained) < maximum:
                retained.append({str(key): str(value) for key, value in row.items()})
    return count, header, retained, matching


def require_table_columns(header: list[str], required: set[str], label: str) -> None:
    missing = sorted(required - set(header))
    if missing:
        raise ValueError(f"{label} lacks required v2 columns: {missing}")


def companion_contract_path(calls_path: str, explicit: str) -> str:
    if explicit:
        return require_file(explicit, "caller contract")
    suffix = ".arm_calls.tsv.gz"
    if not calls_path.endswith(suffix):
        raise ValueError(
            "--call-contract is required when --calls does not end in " + suffix)
    return require_file(
        calls_path[:-len(suffix)] + ".contract.json", "caller contract")


def load_call_contract(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"caller contract is not a JSON object: {path}")
    if clean(payload.get("schema_version")) != CALL_CONTRACT_SCHEMA:
        raise ValueError(
            f"caller contract schema is not {CALL_CONTRACT_SCHEMA}: {path}")
    expected_outputs = {
        "calls": CALL_SCHEMA,
        "calibration": CALIBRATION_SCHEMA,
        "uid_chromosome_flags": UID_SCHEMA,
        "donor_pair_arm_summary": PAIR_SCHEMA,
        "qc": CALL_QC_SCHEMA,
        "contract": CALL_CONTRACT_SCHEMA,
    }
    output_schemas = payload.get("output_schemas")
    if not isinstance(output_schemas, dict) or any(
            clean(output_schemas.get(key)) != value
            for key, value in expected_outputs.items()):
        raise ValueError("caller contract output_schemas do not match report v3")
    return payload


def main_impl(args) -> int:
    calls_path = require_file(args.calls, "arm call table")
    calibration_path = require_file(args.calibration, "calibration table") if args.calibration else ""
    uid_path = require_file(args.uid_summary, "UID summary") if args.uid_summary else ""
    pair_path = (require_file(args.donor_pair_summary, "donor-pair arm summary")
                 if args.donor_pair_summary else "")
    contract_path = companion_contract_path(calls_path, args.call_contract)
    call_contract = load_call_contract(contract_path)
    result = stream_calls(calls_path, args.top_calls)
    calibration_count, calibration_header, calibration_rows, _ = load_optional_rows(
        calibration_path, expected_schema=CALIBRATION_SCHEMA)
    uid_count, uid_header, uid_rows, _ = load_optional_rows(
        uid_path, expected_schema=UID_SCHEMA)
    pair_count, pair_header, pair_rows, recurrent_pair_rows = load_optional_rows(
        pair_path, expected_schema=PAIR_SCHEMA,
        match_field="recurrence_flag", match_value="TRUE")
    if calibration_path:
        require_table_columns(
            calibration_header,
            {"library", "barcode", "crossfit_fold", "excluded_chromosome",
             "cell_loo_baseline_status", "ref_source_level",
             "ref_quasi_overdispersion_rho", "expression_source_level",
             "calibration_status", "schema_version"},
            "calibration table")
    if uid_path:
        require_table_columns(
            uid_header,
            {"uid", "donor_pair", "chromosome", "p_state", "q_state",
             "directional_support_basis", "min_directional_log_bf",
             "whole_chromosome_flag", "summary_status", "schema_version"},
            "UID table")
    if pair_path:
        require_table_columns(
            pair_header,
            {"library", "calibration_group", "donor_pair", "arm",
             "tested_cells", "tested_uid_blocks", "supporting_fraction",
             "directional_support_basis", "min_directional_log_bf",
             "partial_conjunction_method",
             "dependence_assumption", "fdr_interpretation",
             "partial_conjunction_q_value",
             "partial_conjunction_q_resolution_floor", "recurrence_flag",
             "summary_status", "schema_version"},
            "donor-pair arm summary")

    qc_paths = [require_file(path, "QC input") for path in args.qc_files]
    derived_call_qc = (calls_path[:-len(".arm_calls.tsv.gz")] + ".qc.tsv"
                       if calls_path.endswith(".arm_calls.tsv.gz") else "")
    if (derived_call_qc and os.path.isfile(derived_call_qc) and
            os.path.abspath(derived_call_qc) not in qc_paths):
        qc_paths.append(require_file(derived_call_qc, "caller QC"))
    qc_records = []
    qc_failures = []
    call_qc_records = []
    for path in qc_paths:
        values = parse_qc(path)
        record = {key: value for key, value in values}
        record["file"] = os.path.basename(path)
        qc_records.append(record)
        status = clean(record.get("status")).upper()
        if status and not status.startswith("PASS"):
            qc_failures.append({"file": path, "status": status})
        if clean(record.get("schema_version")) == CALL_QC_SCHEMA:
            call_qc_records.append(record)

    if len(call_qc_records) != 1:
        raise ValueError(
            "report requires exactly one tetra_arm_call_qc_v2 input or "
            "auto-associated caller QC file")
    call_qc = call_qc_records[0]
    terminal_state = clean(call_contract.get("terminal_state")).upper() or "NONE"
    contract_status = clean(call_contract.get("status")).upper()
    qc_terminal = clean(call_qc.get("terminal_state")).upper() or "NONE"
    qc_status = clean(call_qc.get("status")).upper()
    if terminal_state not in {"NONE", *TERMINAL_STATES}:
        raise ValueError(f"unknown caller terminal state: {terminal_state}")
    if qc_terminal != terminal_state or qc_status != contract_status:
        raise ValueError("caller QC and contract terminal/status fields disagree")
    expected_call_status = "PASS" if terminal_state == "NONE" else terminal_state
    if contract_status != expected_call_status:
        raise ValueError(
            f"caller status {contract_status!r} disagrees with terminal state "
            f"{terminal_state!r}")
    call_rows = result["totals"].get("rows", 0)
    if terminal_state == "PASS_NO_HETEROTYPIC_TARGETS" and call_rows != 0:
        raise ValueError("PASS_NO_HETEROTYPIC_TARGETS requires header-only calls")
    if terminal_state == "PASS_NO_CALLABLE_ASE" and call_rows == 0:
        raise ValueError("PASS_NO_CALLABLE_ASE requires explicit non-callable rows")
    if terminal_state == "NONE" and call_rows == 0:
        raise ValueError("header-only calls require an explicit caller terminal state")

    totals = result["totals"]
    overall = expected_call_status if not qc_failures else "REVIEW"
    summary_rows = [
        {"metric": "schema_version", "value": REPORT_SCHEMA,
         "definition": "Report contract"},
        {"metric": "status", "value": overall,
         "definition": "Caller PASS/terminal outcome unless an input QC bundle requests review"},
        {"metric": "caller_terminal_state", "value": terminal_state,
         "definition": "Exact v2 caller scientific terminal state; NONE for a regular call"},
        {"metric": "call_rows", "value": totals.get("rows", 0),
         "definition": "Cell by chromosome-arm result rows"},
        {"metric": "cells", "value": totals.get("cells", 0),
         "definition": "Distinct library/barcode keys represented"},
        {"metric": "nonbalanced_best_state_rows",
         "value": totals.get("nonbalanced_best_state_rows", 0),
         "definition": "Rows whose maximum-posterior state is not balanced"},
        {"metric": "pass_event_rows", "value": totals.get("pass_event_rows", 0),
         "definition": "Arm rows passing posterior, empirical FDR, and evidence QC"},
        {"metric": "pass_event_cells", "value": totals.get("pass_event_cells", 0),
         "definition": "Distinct cells with at least one PASS_EVENT arm"},
        {"metric": "nonbalanced_candidate_rows_q05",
         "value": totals.get("nonbalanced_candidate_rows_q05", 0),
         "definition": "Nonbalanced best states with empirical q at most 0.05, including review candidates"},
        {"metric": "exploratory_non_autosomal_rows",
         "value": totals.get("exploratory_non_autosomal_rows", 0),
         "definition": "Explicitly exploratory sex/non-autosomal rows; never production PASS"},
        {"metric": "calibration_rows", "value": calibration_count,
         "definition": "Query-level calibration records written by the caller"},
        {"metric": "uid_summary_rows", "value": uid_count,
         "definition": "Physical-fusion UID aggregate rows"},
        {"metric": "donor_pair_arm_summary_rows", "value": pair_count,
         "definition": "Library/group/donor-pair arm aggregate rows"},
        {"metric": "recurrent_donor_pair_arm_events", "value": recurrent_pair_rows,
         "definition": "Primary recurrence rows passing the caller's working-FDR rules"},
        {"metric": "qc_files", "value": len(qc_paths),
         "definition": "Worker QC tables included"},
        {"metric": "qc_review_files", "value": len(qc_failures),
         "definition": "Worker QC tables whose status is not PASS"},
    ]

    output_dir = Path(os.path.abspath(args.output_dir))
    summary_path = output_dir / "tetra_arm_cnv_summary.tsv"
    json_path = output_dir / "tetra_arm_cnv_summary.json"
    html_path = output_dir / "tetra_arm_cnv_report.html"
    require_outputs_absent((summary_path, json_path, html_path))
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv_atomic(str(summary_path), summary_rows, ["metric", "value", "definition"])

    payload = {
        "schema_version": REPORT_SCHEMA,
        "release": PROGRAM_VERSION,
        "status": overall,
        "caller_terminal_state": terminal_state,
        "inputs": {
            "calls": file_record(calls_path),
            "calibration": file_record(calibration_path) if calibration_path else None,
            "uid_summary": file_record(uid_path) if uid_path else None,
            "donor_pair_summary": file_record(pair_path) if pair_path else None,
            "call_contract": file_record(contract_path),
            "qc": [file_record(path) for path in qc_paths],
        },
        "totals": totals,
        "state_counts": dict(result["states"]),
        "best_state_counts": dict(result["best_states"]),
        "status_counts": dict(result["statuses"]),
        "evidence_counts": dict(result["evidence"]),
        "library_counts": dict(result["libraries"]),
        "arm_counts": dict(result["arms"]),
        "donor_pair_counts": dict(result["donor_pairs"]),
        "calibration_rows": calibration_count,
        "uid_summary_rows": uid_count,
        "donor_pair_arm_summary_rows": pair_count,
        "recurrent_donor_pair_arm_events": recurrent_pair_rows,
        "qc_review": qc_failures,
        "top_nonbalanced_calls": result["top"],
    }
    write_json_atomic(str(json_path), payload)

    top_headers = [
        "Library", "Barcode", "Arm", "Donor pair", "Best state",
        "Working pseudo-posterior", "Selected-state empirical q",
        "q-resolution bound", "Inference scope", "Status"]
    top_rows = [[row[key] for key in (
        "library", "barcode", "arm", "donor_pair", "best_state", "posterior",
        "q_value", "q_resolution", "inference_scope", "status")]
        for row in result["top"]]
    calibration_preview = []
    if calibration_rows:
        preferred = [field for field in (
            "library", "barcode", "uid", "calibration_group", "donor_pair",
            "chromosome", "arm", "crossfit_fold", "excluded_chromosome",
            "target_effective_ase_weight", "group_genomewide_baseline_logit",
            "cell_loo_baseline_donor_a_fraction", "cell_loo_baseline_arms",
            "cell_loo_baseline_status", "ref_source_level",
            "ref_calibration_cells", "ref_orientation_mapping_logit_offset",
            "ref_quasi_overdispersion_rho", "alt_source_level",
            "alt_orientation_mapping_logit_offset",
            "alt_quasi_overdispersion_rho", "mixed_source_level",
            "mixed_quasi_overdispersion_rho", "expression_source_level",
            "expression_metric", "expression_calibration_cells",
            "calibration_status") if field in calibration_rows[0]]
        calibration_preview = [[row.get(field, "") for field in preferred]
                               for row in calibration_rows[:50]]
    else:
        preferred = []
    uid_preview = []
    uid_fields = []
    if uid_rows:
        uid_fields = [field for field in (
            "uid", "donor_pair", "chromosome", "n_cells", "libraries",
            "p_arm", "p_evaluable_cells", "p_concordant_cells", "p_state",
            "p_best_state_posterior", "p_empirical_q_value",
            "p_empirical_q_resolution_floor", "q_arm", "q_evaluable_cells",
            "q_concordant_cells", "q_state", "q_best_state_posterior",
            "q_empirical_q_value", "q_empirical_q_resolution_floor",
            "directional_support_basis", "min_directional_log_bf",
            "paired_pq_concordant_cells",
            "whole_chromosome_flag", "whole_chromosome_state", "summary_status")
            if field in uid_rows[0]]
        uid_preview = [[row.get(field, "") for field in uid_fields] for row in uid_rows[:50]]
    pair_preview = []
    pair_fields = []
    if pair_rows:
        pair_fields = [field for field in (
            "library", "calibration_group", "donor_pair", "arm", "libraries",
            "total_cells", "eligible_uid_blocks", "tested_cells",
            "tested_uid_blocks", "qc_blocked_cells", "supporting_cells",
            "supporting_fraction", "concordant_cells", "concordance",
            "best_state",
            "best_state_posterior", "pooled_effective_units",
            "partial_conjunction_r", "partial_conjunction_method",
            "fisher_terms", "fisher_df", "dependence_assumption",
            "directional_support_basis", "min_directional_log_bf",
            "partial_conjunction_p_value", "partial_conjunction_q_value",
            "partial_conjunction_q_resolution_floor", "fdr_interpretation",
            "recurrence_flag", "summary_status") if field in pair_rows[0]]
        pair_preview = [[row.get(field, "") for field in pair_fields]
                        for row in pair_rows[:50]]

    css = """
    :root { color-scheme: light; --ink:#17212b; --muted:#596875; --line:#d9e0e6;
      --accent:#176b87; --soft:#eef6f8; --warn:#9a5d00; }
    body { font-family: system-ui,-apple-system,Segoe UI,sans-serif; color:var(--ink);
      max-width:1280px; margin:0 auto; padding:28px; line-height:1.42; }
    h1 { margin-bottom:4px; } h2 { margin-top:32px; border-bottom:1px solid var(--line); padding-bottom:6px; }
    .subtle { color:var(--muted); } .cards { display:grid; grid-template-columns:repeat(auto-fit,minmax(175px,1fr)); gap:12px; }
    .card { background:var(--soft); border:1px solid var(--line); border-radius:8px; padding:14px; }
    .number { font-size:1.7rem; font-weight:700; } table { border-collapse:collapse; width:100%; font-size:.9rem; }
    th,td { border:1px solid var(--line); padding:6px 8px; text-align:left; vertical-align:top; }
    th { background:#f5f7f9; position:sticky; top:0; } .bar-track { width:100%; min-width:90px; background:#e8edf0; height:10px; }
    .bar-fill { height:10px; background:var(--accent); } code { background:#f3f5f6; padding:1px 4px; }
    .review { color:var(--warn); font-weight:700; }
    """
    cards = [
        ("Run status", overall),
        ("Caller terminal", terminal_state),
        ("Cells", f"{totals.get('cells', 0):,}"),
        ("Cell-arm rows", f"{totals.get('rows', 0):,}"),
        ("PASS_EVENT cells", f"{totals.get('pass_event_cells', 0):,}"),
        ("Calibration rows", f"{calibration_count:,}"),
        ("UID aggregate rows", f"{uid_count:,}"),
        ("Recurrent pair-arms", f"{recurrent_pair_rows:,}"),
    ]
    card_html = "".join(
        f'<div class="card"><div class="subtle">{html.escape(label)}</div>'
        f'<div class="number">{html.escape(value)}</div></div>' for label, value in cards)
    terminal_descriptions = {
        "NONE": "Regular caller completion with callable ASE rows.",
        "PASS_NO_HETEROTYPIC_TARGETS": (
            "No heterotypic two-donor target cells were present; the call table "
            "is intentionally header-only."),
        "PASS_NO_OBSERVED_ASE": (
            "Heterotypic targets were present but no ASE row was observed. "
            "Expression-only NO_DATA rows may still be present."),
        "PASS_NO_CALLABLE_ASE": (
            "ASE rows were present, but none had eligible positive directional "
            "evidence for calling."),
    }
    sections = [
        "<!doctype html><html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
        f"<title>{html.escape(args.title)}</title><style>{css}</style></head><body>",
        f"<h1>{html.escape(args.title)}</h1>",
        f'<p class="subtle">Schema <code>{REPORT_SCHEMA}</code>. State weights are working pseudo-posteriors; production events require selected-state held-out empirical q support. Aggregate donor-direction concordance uses prior-independent ASE event-vs-balanced evidence, while individual-cell calls retain the conservative event prior. Downstream identity and modality fields are QC context, not duplicated likelihood evidence.</p>',
        f'<div class="cards">{card_html}</div>',
        "<h2>Caller outcome</h2>",
        f"<p><code>{html.escape(terminal_state)}</code>: "
        f"{html.escape(terminal_descriptions[terminal_state])}</p>",
        "<h2>Final call-state distribution</h2>", bar_table(result["states"], totals.get("rows", 0)),
        "<h2>Maximum-posterior candidate state</h2>", bar_table(result["best_states"], totals.get("rows", 0)),
        "<h2>Call status</h2>", bar_table(result["statuses"], totals.get("rows", 0)),
        "<h2>Evidence status</h2>", bar_table(result["evidence"], totals.get("rows", 0)),
        "<h2>Top nonbalanced calls</h2>",
        small_table(top_headers, top_rows) if top_rows else "<p>No rankable nonbalanced calls.</p>",
    ]
    if preferred:
        sections.extend(["<h2>Calibration preview</h2>", small_table(preferred, calibration_preview)])
    if uid_fields:
        sections.extend(["<h2>Physical-fusion UID preview</h2>", small_table(uid_fields, uid_preview)])
    if pair_fields:
        sections.extend([
            "<h2>Donor-pair arm recurrence preview</h2>",
            "<p class=\"subtle\">This is the primary library/group recurrence view. Its Fisher partial-conjunction q-values are working-FDR evidence under the reported dependence assumption.</p>",
            small_table(pair_fields, pair_preview)])
    if qc_records:
        qc_fields = sorted({key for row in qc_records for key in row}, key=natural_key)
        qc_table = [[row.get(field, "") for field in qc_fields] for row in qc_records]
        sections.extend(["<h2>Worker QC</h2>", small_table(qc_fields, qc_table)])
    sections.extend([
        "<h2>Interpretation</h2>",
        "<p>A gain or loss is reportable only when donor-directed ASE has usable interindividual sites and the selected aggregate state passes its held-out empirical q, resolution, recurrence, and donor-direction concordance thresholds. Expression supports total-copy direction when present. Individual-cell event calls remain prior-regularized and are not prerequisites for aggregate support. A nonbalanced maximum working pseudo-posterior without matching state-specific empirical support remains a candidate. Sex-chromosome rows are exploratory only.</p>",
        "</body></html>",
    ])
    with atomic_text(str(html_path), gzip_output=False) as handle:
        handle.write("\n".join(sections))
        handle.write("\n")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a self-contained audit report for calibrated chromosome-arm calls.")
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {PROGRAM_VERSION}")
    parser.add_argument("--calls", required=True)
    parser.add_argument("--calibration", default="")
    parser.add_argument("--uid-summary", default="")
    parser.add_argument("--donor-pair-summary", default="")
    parser.add_argument(
        "--call-contract", default="",
        help=("Caller v2 contract; defaults to the .contract.json companion "
              "of --calls"))
    parser.add_argument("--qc-files", nargs="*", default=[])
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--title", default="Tetraploid chromosome-arm ASE/CNV report")
    parser.add_argument("--top-calls", type=int, default=100)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.top_calls < 1 or args.top_calls > 10000:
        print("ERROR: --top-calls must be between 1 and 10000", file=sys.stderr)
        return 2
    try:
        return main_impl(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
