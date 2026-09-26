#!/usr/bin/env python3
"""Summarize completed tetra-arm CALL results without exporting cell-level data.

Uses only Python's standard library. The existing REPORT outputs are left intact.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import html
import json
import os
from pathlib import Path
import re
import sys
import tempfile


VERSION = "1.0.0"
PREFIX = "tetra_arm_cnv_hybrid"
EVIDENCE = (
    "CONCORDANT_BOTH", "EXPRESSION_ONLY", "ASE_ONLY", "DISCORDANT",
    "EXPRESSION_OUTLIER", "BALANCED", "INSUFFICIENT_EVIDENCE",
)
EXACT_STATES = {"DONOR_A_LOSS", "DONOR_B_LOSS", "DONOR_A_GAIN", "DONOR_B_GAIN"}
CANDIDATES = EVIDENCE[:3]
REQUIRED_CALL = {
    "library", "barcode", "arm", "evidence_class", "resolved_state",
    "confidence_tier", "call_status", "schema_version",
}
REQUIRED_MODEL = {
    "library", "barcode", "arm", "ase_eligible", "expression_eligible",
    "ase_status", "expression_status", "expression_fold_replication_status",
    "discordance_reason",
}
REQUIRED_PAIR = {
    "donor_pair", "arm", "distinct_uid_blocks", "summary_status", "total_cells",
    "expression_support_uid_blocks", "ase_support_uid_blocks", "concordant_uid_blocks",
    "pairwide_identifiability_flag", "best_resolved_state",
}
REQUIRED_UID = {
    "uid", "donor_pair", "chromosome", "summary_status", "n_cells",
    "distinct_uid_blocks", "p_resolved_state", "q_resolved_state",
}


def tsv_rows(path: Path, required: set[str]):
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {', '.join(sorted(missing))}")
        for row in reader:
            yield row


def tsv_header(path: Path) -> set[str]:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return set(next(csv.reader(handle, delimiter="\t"), []))


def json_read(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def yes(value: str) -> bool:
    return str(value).strip().upper() in {"1", "TRUE", "YES"}


def number(value: str, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def known(value: str) -> bool:
    return str(value or "").strip().upper() not in {"", ".", "NA", "N/A", "NONE", "NULL"}


def percentage(part: int, total: int) -> str:
    return f"{100.0 * part / total:.1f}%" if total else "n/a"


def arm_order(arm: str):
    match = re.search(r"(?:chr)?(\d+)([pq]?)$", arm, re.I)
    return (int(match.group(1)), {"p": 0, "q": 1, "": 2}[match.group(2).lower()], arm) if match else (9999, 3, arm)


def write_atomic(path: Path, write):
    descriptor, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            write(handle)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def table_atomic(path: Path, fields: list[str], rows: list[dict]):
    def emit(handle):
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    write_atomic(path, emit)


def rollup(run: Path) -> dict:
    call_prefix = run / "call" / PREFIX
    calls_file = Path(f"{call_prefix}.arm_calls.tsv.gz")
    pair_file = Path(f"{call_prefix}.donor_pair_arm_summary.tsv.gz")
    uid_file = Path(f"{call_prefix}.uid_chromosome_flags.tsv.gz")
    call_contract_file = Path(f"{call_prefix}.contract.json")
    for path in (calls_file, pair_file, uid_file, call_contract_file):
        if not path.is_file():
            raise FileNotFoundError(f"Required completed CALL output missing: {path}")

    by_arm: dict[str, Counter] = defaultdict(Counter)
    by_cell: dict[tuple[str, str], Counter] = defaultdict(Counter)
    evidence, tiers, statuses, flags, schema_versions = Counter(), Counter(), Counter(), Counter(), Counter()
    libraries: set[str] = set()
    # This index is bounded by observed CALL rows, and contains no expression arrays.
    call_index: dict[tuple[str, str, str], str] = {}
    q_coverage = Counter()
    call_fields = tsv_header(calls_file)
    rows = 0
    for row in tsv_rows(calls_file, REQUIRED_CALL):
        lib, barcode, arm, cls = (row[name].strip() for name in ("library", "barcode", "arm", "evidence_class"))
        if not lib or not barcode or not arm or not cls:
            raise ValueError(f"{calls_file}: empty library/barcode/arm/evidence_class on row {rows + 2}")
        if cls not in EVIDENCE:
            raise ValueError(f"{calls_file}: unknown evidence_class {cls!r} on row {rows + 2}")
        if not row["schema_version"].strip():
            raise ValueError(f"{calls_file}: blank schema_version on row {rows + 2}")
        schema_versions[row["schema_version"].strip()] += 1
        key = (lib, barcode, arm)
        if key in call_index:
            raise ValueError(f"{calls_file}: duplicate library/barcode/arm on row {rows + 2}")
        call_index[key] = cls
        libraries.add(lib)
        rows += 1
        cell = by_cell[(lib, barcode)]
        arm_stats = by_arm[arm]
        arm_stats["observed"] += 1
        arm_stats[cls] += 1
        cell["observed"] += 1
        cell[cls] += 1
        evidence[cls] += 1
        tiers[row.get("confidence_tier", "") or "UNREPORTED"] += 1
        statuses[row.get("call_status", "") or "UNREPORTED"] += 1
        for flag in (row.get("confounding_flags") or "").split(";"):
            if known(flag):
                flags[flag.strip()] += 1
        for field in ("expression_q", "ase_q", "hybrid_q"):
            raw = row.get(field, "")
            if known(raw):
                try:
                    value = float(raw)
                except ValueError:
                    raise ValueError(f"{calls_file}: invalid {field} on row {rows + 1}: {raw}") from None
                if not 0 <= value <= 1:
                    raise ValueError(f"{calls_file}: {field} outside [0,1] on row {rows + 1}")
                q_coverage[f"{field}_available"] += 1
                if value <= 0.05:
                    q_coverage[f"{field}_at_most_0_05"] += 1
    if not rows:
        raise ValueError(f"{calls_file}: no CALL rows; no report generated")
    if len(schema_versions) != 1:
        raise ValueError(f"{calls_file}: mixed CALL schema versions {dict(schema_versions)}")

    score_files = sorted((run / "model").glob("*.scores.tsv.gz"))
    score_coverage = Counter()
    seen_scores: set[tuple[str, str, str]] = set()
    ase_status, expression_status, fold_status, reasons = Counter(), Counter(), Counter(), Counter()
    model_flags = Counter()
    if not score_files:
        score_coverage["files_missing"] = 1
    for path in score_files:
        for score in tsv_rows(path, REQUIRED_MODEL):
            key = (score["library"].strip(), score["barcode"].strip(), score["arm"].strip())
            if key in seen_scores:
                score_coverage["duplicate_model_keys"] += 1
                continue
            seen_scores.add(key)
            if key not in call_index:
                score_coverage["model_rows_without_call"] += 1
                continue
            score_coverage["matched"] += 1
            arm_stats = by_arm[key[2]]
            cell = by_cell[key[:2]]
            for name in ("ase_eligible", "expression_eligible"):
                if yes(score.get(name, "")):
                    arm_stats[name] += 1
                    cell[name] += 1
                    score_coverage[name] += 1
            if yes(score.get("ase_eligible", "")) and yes(score.get("expression_eligible", "")):
                arm_stats["both_eligible"] += 1
                cell["both_eligible"] += 1
                score_coverage["both_eligible"] += 1
            for name, dest in (("ase_status", ase_status), ("expression_status", expression_status),
                               ("expression_fold_replication_status", fold_status), ("discordance_reason", reasons)):
                dest[score.get(name) or "UNREPORTED"] += 1
            if call_index[key] == "EXPRESSION_OUTLIER":
                reason = score.get("discordance_reason") or "UNREPORTED"
                arm_stats[f"outlier_reason:{reason}"] += 1
                model_flags[reason] += 1
    # A duplicate/misaligned score inventory would otherwise inflate power rates.
    if score_coverage["matched"] > rows or score_coverage["model_rows_without_call"] or score_coverage["duplicate_model_keys"]:
        score_coverage["alignment_warning"] = 1
    score_coverage["unmatched_call_rows"] = max(0, rows - score_coverage["matched"])
    complete_model_diagnostics = score_coverage["matched"] == rows and not score_coverage["alignment_warning"]

    pair_entries = []
    by_pair: dict[str, dict] = {}
    for row in tsv_rows(pair_file, REQUIRED_PAIR):
        pair = row["donor_pair"].strip() or "UNREPORTED"
        entry = by_pair.setdefault(pair, {
            "donor_pair": pair, "logical_arm_summaries": 0, "distinct_uid_blocks_max": 0,
            "cells_max": 0, "resolved_arm_summaries": 0, "concordant_uid_block_summaries": 0,
            "pairwide_identifiability_flag_summaries": 0, "expression_support_uid_block_summaries": 0,
            "ase_support_uid_block_summaries": 0,
        })
        entry["logical_arm_summaries"] += 1
        entry["distinct_uid_blocks_max"] = max(entry["distinct_uid_blocks_max"], number(row.get("distinct_uid_blocks")))
        entry["cells_max"] = max(entry["cells_max"], number(row.get("total_cells")))
        entry["expression_support_uid_block_summaries"] += int(number(row.get("expression_support_uid_blocks")) > 0)
        entry["ase_support_uid_block_summaries"] += int(number(row.get("ase_support_uid_blocks")) > 0)
        entry["concordant_uid_block_summaries"] += int(number(row.get("concordant_uid_blocks")) > 0)
        entry["pairwide_identifiability_flag_summaries"] += int(yes(row.get("pairwide_identifiability_flag", "")))
        entry["resolved_arm_summaries"] += int(
            row.get("summary_status") == "PASS_EVENT" or row.get("best_resolved_state") in EXACT_STATES
        )
        pair_entries.append(row)

    by_uid: dict[tuple[str, str], dict] = {}
    uid_statuses = Counter()
    for row in tsv_rows(uid_file, REQUIRED_UID):
        uid_statuses[row.get("summary_status") or "UNREPORTED"] += 1
        uid, pair = row["uid"].strip(), row["donor_pair"].strip()
        entry = by_uid.setdefault((uid, pair), {
            "uid": uid, "donor_pair": pair, "logical_group_summaries": 0,
            "cells_max": 0, "resolved_group_summaries": 0,
            "distinct_uid_blocks_max": 0,
        })
        entry["logical_group_summaries"] += 1
        entry["cells_max"] = max(entry["cells_max"], number(row.get("n_cells")))
        entry["distinct_uid_blocks_max"] = max(entry["distinct_uid_blocks_max"], number(row.get("distinct_uid_blocks")))
        entry["resolved_group_summaries"] += int(
            row.get("summary_status") == "PASS_EVENT" or
            any(row.get(field) in EXACT_STATES
                for field in ("p_resolved_state", "q_resolved_state"))
        )

    burden = []
    for category in ("CONCORDANT_BOTH", "EXPRESSION_ONLY", "ASE_ONLY", "DISCORDANT", "EXPRESSION_OUTLIER"):
        bins = Counter()
        for cell in by_cell.values():
            value = cell[category]
            label = "0" if value == 0 else "1" if value == 1 else "2" if value == 2 else "3-4" if value <= 4 else "5-9" if value <= 9 else "10+"
            bins[label] += 1
        for label in ("0", "1", "2", "3-4", "5-9", "10+"):
            burden.append({"evidence_class": category, "arms_per_cell": label, "cells": bins[label],
                           "percent_of_observed_cells": percentage(bins[label], len(by_cell))})

    arm_rows = []
    for arm, counts in sorted(by_arm.items(), key=lambda item: arm_order(item[0])):
        record = {
            "logical_arm": arm, "observed_cell_arm_rows": counts["observed"],
            "ase_eligible": counts["ase_eligible"] if complete_model_diagnostics else "",
            "expression_eligible": counts["expression_eligible"] if complete_model_diagnostics else "",
            "both_eligible": counts["both_eligible"] if complete_model_diagnostics else "",
            "both_eligible_percent": percentage(counts["both_eligible"], counts["observed"]) if complete_model_diagnostics else "n/a",
            "candidate_rows": sum(counts[name] for name in CANDIDATES),
            "expression_outlier_percent": percentage(counts["EXPRESSION_OUTLIER"], counts["observed"]),
            "insufficient_percent": percentage(counts["INSUFFICIENT_EVIDENCE"], counts["observed"]),
        }
        record.update({name.lower(): counts[name] for name in EVIDENCE})
        arm_rows.append(record)

    call_contract = json_read(call_contract_file)
    model_contracts = sorted((run / "model").glob("*.contract.json"))
    model_states = Counter(json_read(path).get("status", "UNREPORTED") for path in model_contracts)
    overview = {
        "observed_cell_arm_rows": rows, "unique_observed_cells": len(by_cell),
        "call_schema_version": next(iter(schema_versions)),
        "distinct_logical_arms": len(by_arm), "libraries": sorted(libraries, key=lambda v: (not v.isdigit(), number(v), v)),
        "evidence_class": dict(evidence), "confidence_tier": dict(tiers), "call_status": dict(statuses),
        "candidate_rows": sum(evidence[name] for name in CANDIDATES),
        "resolved_calls": statuses["PASS_EVENT"],
        "expression_outlier_cells": sum(cell["EXPRESSION_OUTLIER"] > 0 for cell in by_cell.values()),
        "call_q_fields_present": {name: name in call_fields for name in ("expression_q", "ase_q", "hybrid_q")},
        "call_confounding_flags_field_present": "confounding_flags" in call_fields,
        "q_value_coverage_at_report_threshold_0_05": dict(q_coverage),
        "confounding_flags": dict(flags),
        "model_score_coverage": dict(score_coverage),
        "complete_model_diagnostics": complete_model_diagnostics,
        "ase_status": dict(ase_status), "expression_status": dict(expression_status),
        "expression_fold_replication_status": dict(fold_status),
        "discordance_reason": dict(reasons),
        "expression_outlier_reasons": dict(model_flags),
        "pair_arm_summary_rows": len(pair_entries), "uid_logical_group_summary_rows": sum(v["logical_group_summaries"] for v in by_uid.values()),
        "unique_uids": len(by_uid), "uid_group_status": dict(uid_statuses),
        "execution": {
            "call_contract_status": call_contract.get("status", "UNREPORTED"),
            "call_contract_release": call_contract.get("release", "UNREPORTED"),
            "model_contracts": len(model_contracts), "model_contract_statuses": dict(model_states),
        },
    }
    return {
        "schema_version": "tetra_arm_interpretive_summary_v1",
        "run_root": str(run),
        "input_files": [str(calls_file), str(pair_file), str(uid_file), str(call_contract_file)] + [str(p) for p in score_files],
        "overview": overview, "arms": arm_rows,
        "cell_burden": burden, "donor_pairs": sorted(by_pair.values(), key=lambda entry: entry["donor_pair"]),
        "uid_support": sorted(by_uid.values(), key=lambda entry: (entry["donor_pair"], entry["uid"])),
    }


def html_report(result: dict) -> str:
    d, arms = result["overview"], result["arms"]
    n = d["observed_cell_arm_rows"]
    ev = d["evidence_class"]
    get = lambda name: ev.get(name, 0)
    esc = lambda value: html.escape(str(value), quote=True)
    def fmt(value):
        return f"{value:,}" if isinstance(value, int) else esc(value)
    def table(fields, rows):
        return ("<div class='scroll'><table><thead><tr>" + "".join(f"<th>{esc(title)}</th>" for _, title in fields) +
                "</tr></thead><tbody>" + "".join("<tr>" + "".join(f"<td>{fmt(row.get(key, ''))}</td>" for key, _ in fields) + "</tr>" for row in rows) +
                "</tbody></table></div>")
    cards = [
        ("Observed cells", d["unique_observed_cells"]), ("Observed cell-arm rows", n),
        ("Logical arms", d["distinct_logical_arms"]), ("Exact-state PASS_EVENT rows", d["resolved_calls"]),
        ("Single-modality candidate rows", get("EXPRESSION_ONLY") + get("ASE_ONLY")),
        ("Expression outlier rows", get("EXPRESSION_OUTLIER")),
        ("Cells with an outlier", d["expression_outlier_cells"]),
    ]
    if d["resolved_calls"] == 0 and d["candidate_rows"] == 0:
        headline = "No resolved RNA-inferred arm events were identified at this run's thresholds."
        explanation = (f"{get('EXPRESSION_OUTLIER'):,} expression outlier rows are review flags, not CNV calls; "
                       f"{get('INSUFFICIENT_EVIDENCE'):,} rows lack sufficient evidence, and {get('BALANCED'):,} were classified balanced. "
                       "This outcome does not establish that all cells are copy-number balanced.")
    elif d["resolved_calls"] == 0:
        headline = ("No exact-state PASS_EVENT rows; "
                    f"{d['candidate_rows']:,} partial or confounded candidate rows require review.")
        explanation = ("Expression-only rows lack donor orientation and ASE-only rows cannot distinguish loss from gain. "
                       "These are cell-arm observations, not validated DNA CNV events.")
    else:
        headline = f"{d['resolved_calls']:,} exact-state RNA-inferred PASS_EVENT rows require biological review."
        explanation = ("Concordant, expression-only, and ASE-only candidate classes have different interpretations. "
                       "Candidate rows are cell-arm observations, not independent clones or validated DNA CNVs.")
    denominator_note = ("All displayed percentages use observed CALL cell-arm rows for the corresponding arm. "
                        "Cells are counted once in the burden summary. Missing cell-arm combinations are not assumed balanced or powered.")
    if d["complete_model_diagnostics"]:
        score = d["model_score_coverage"]
        power = (f"Among {n:,} matched model score rows, {score.get('both_eligible', 0):,} "
                 f"({percentage(score.get('both_eligible', 0), n)}) were eligible for both ASE and expression. "
                 "Eligibility reports whether a branch could be evaluated; it is not an event call.")
    else:
        score = d["model_score_coverage"]
        power = (f"Model score diagnostics are incomplete or misaligned ({score.get('matched', 0):,}/{n:,} matched rows). "
                 "Eligibility fractions are suppressed until all CALL rows have exactly one matching model score.")
    execution = d["execution"]
    model_states = execution["model_contract_statuses"]
    execution_text = (f"CALL contract: {esc(execution['call_contract_status'])}; "
                      f"MODEL contracts: {execution['model_contracts']:,}, status counts {esc(model_states)}. "
                      "PASS is execution/QC status and does not mean that any biological event was detected.")
    distribution = "".join(
        f"<div class='distrow'><span>{esc(label.replace('_', ' ').title())}</span>"
        f"<div class='track'><div class='fill {esc(label.lower())}' style='width:{min(100, 100*get(label)/n):.2f}%'></div></div>"
        f"<strong>{get(label):,}</strong><span>{percentage(get(label), n)}</span></div>"
        for label in EVIDENCE
    )
    arm_fields = [("logical_arm", "Logical arm"), ("observed_cell_arm_rows", "Observed rows"),
                  ("both_eligible_percent", "Both eligible"), ("concordant_both", "Concordant"),
                  ("expression_only", "Expression only"), ("ase_only", "ASE only"),
                  ("discordant", "Branch discordant"), ("expression_outlier", "Expression outliers"),
                  ("expression_outlier_percent", "Outlier rate"), ("balanced", "Balanced"),
                  ("insufficient_evidence", "Insufficient")]
    burden_rows = [r for r in result["cell_burden"] if r["evidence_class"] == "EXPRESSION_OUTLIER"]
    reason_rows = [{"reason": key, "rows": value} for key, value in sorted(d["expression_outlier_reasons"].items(), key=lambda item: (-item[1], item[0]))]
    if not d["complete_model_diagnostics"]:
        reason_rows = [{"reason": "Incomplete or misaligned MODEL score inventory", "rows": "n/a"}]
    elif not reason_rows:
        reason_rows = [{"reason": "Model score diagnostics unavailable", "rows": "n/a"}]
    status_rows = [{"status": key, "rows": value, "percent": percentage(value, n)} for key, value in sorted(d["call_status"].items())]
    q_rows = [{"branch": name.replace("_", " "),
               "q_available": d["q_value_coverage_at_report_threshold_0_05"].get(f"{name}_available", 0) if d["call_q_fields_present"][name] else "field absent",
               "q_at_most_0_05": d["q_value_coverage_at_report_threshold_0_05"].get(f"{name}_at_most_0_05", 0) if d["call_q_fields_present"][name] else "n/a"}
              for name in ("expression_q", "ase_q", "hybrid_q")]
    branch_rows = [{"branch": name, "status": key, "rows": count} for name, mapping in
                   (("ASE", d["ase_status"]), ("Expression", d["expression_status"]),
                    ("Expression fold replication", d["expression_fold_replication_status"]))
                   for key, count in sorted(mapping.items(), key=lambda item: (-item[1], item[0]))]
    if not d["complete_model_diagnostics"]:
        branch_rows = [{"branch": "MODEL score diagnostics", "status": "Incomplete or misaligned", "rows": "n/a"}]
    flag_rows = [{"flag": key, "rows": value} for key, value in sorted(d["confounding_flags"].items(), key=lambda x: (-x[1], x[0]))]
    pair_rows = result["donor_pairs"]
    uid_status_rows = [{"status": key, "logical_group_summaries": value} for key, value in
                       sorted(d["uid_group_status"].items(), key=lambda item: (-item[1], item[0]))]
    def section(title, content):
        return f"<section><h2>{esc(title)}</h2>{content}</section>"
    styles = """
    :root{font-family:system-ui,sans-serif;color:#193047;background:#f5f7fa}body{max-width:1180px;margin:0 auto;padding:30px 22px 80px}
    h1{font-size:2.1rem;margin:0 0 12px}h2{font-size:1.25rem;margin:0 0 14px}p{line-height:1.55}
    .kicker{color:#526980;font-weight:700;letter-spacing:.1em;text-transform:uppercase;font-size:.75rem}
    .lede{font-size:1.24rem;font-weight:700;color:#123d59}.caution{padding:15px 17px;background:#fff5e5;border-left:4px solid #c87d15}
    .cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(168px,1fr));gap:10px;margin:20px 0}
    .card,section{border:1px solid #dbe4eb;background:white;border-radius:10px;padding:20px;margin:17px 0}
    .card{margin:0;padding:15px}.card strong{display:block;font-size:1.55rem;color:#145574}.card span{color:#526980;font-size:.85rem}
    .muted{color:#536878;font-size:.93rem}.distrow{display:grid;grid-template-columns:minmax(150px,1.8fr) minmax(100px,4fr) 85px 65px;gap:10px;align-items:center;margin:10px 0;font-size:.92rem}
    .track{height:15px;border-radius:10px;background:#e9eef2;overflow:hidden}.fill{height:100%;background:#467d96}
    .expression_outlier{background:#e79a38}.balanced{background:#2d9c79}.insufficient_evidence{background:#a6b6c4}
    .scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:.84rem}th,td{text-align:right;padding:8px 10px;border-bottom:1px solid #e4eaf0;white-space:nowrap}
    th:first-child,td:first-child{text-align:left}th{position:sticky;top:0;background:#eaf2f6;color:#17374c}
    tbody tr:nth-child(even){background:#f8fafb}details{margin-top:12px}summary{cursor:pointer;color:#145574;font-weight:600}
    """
    sections = [
        section("What the run found", f"<p>{esc(explanation)}</p><p class='muted'>{esc(denominator_note)}</p>{distribution}"),
        section("Which logical arms carried the review signal", "<p class='muted'>Outlier rates use observed rows for each arm. These ancestral-reference logical arms are labels assigned by the explicit BED; they are not physical contigs.</p>" + table(arm_fields, arms)),
        section("How outliers are distributed across cells", "<p class='muted'>This measures concentration of expression-outlier review flags across distinct observed cells. One cell can have multiple flagged arms.</p>" + table([("arms_per_cell", "Outlier arms per cell"), ("cells", "Cells"), ("percent_of_observed_cells", "% of cells")], burden_rows)),
        section("Power and review reasons", "<p>" + esc(power) + "</p>" +
                "<h3>Reported adjusted q-values</h3><p class='muted'>Expression and ASE q-values summarize the separate branches; hybrid q-values concern exact-state conjunction from the conservative maximum of branch p-values. A small q-value alone does not bypass fold and quality gates. 'Field absent' means the completed CALL table lacks that optional summary column.</p>" +
                table([("branch", "Branch"), ("q_available", "Nonmissing q"), ("q_at_most_0_05", "q ≤ 0.05")], q_rows) +
                "<h3>Outlier subtypes</h3>" + table([("reason", "MODEL discordance reason"), ("rows", "Outlier rows")], reason_rows) +
                "<details><summary>Branch and fold statuses</summary>" + table([("branch", "Branch"), ("status", "Status"), ("rows", "Rows")], branch_rows) + "</details>" +
                "<h3>CALL disposition</h3>" + table([("status", "Status"), ("rows", "Rows"), ("percent", "Observed rows")], status_rows) +
                "<details><summary>Confounding flags by row</summary>" +
                ("<p class='muted'>Confounding flags were absent from the CALL schema.</p>" if not d["call_confounding_flags_field_present"] else
                 "<p class='muted'>Flags can overlap. Their counts are not an additional event tally.</p>" + table([("flag", "Flag"), ("rows", "Rows")], flag_rows)) + "</details>"),
        section("Support across donor pairs and UIDs", "<p class='muted'>The source summary tables reduce cells within a UID to one biological block. Counts below report arm/group summaries with at least one support block; they are not sums of independent cells. Maximum cells and UID blocks are used because counts recur across arms.</p>" +
                table([("donor_pair", "Donor pair"), ("cells_max", "Cells (maximum)"), ("distinct_uid_blocks_max", "UID blocks (maximum)"),
                       ("logical_arm_summaries", "Arm summaries"), ("expression_support_uid_block_summaries", "Arms with expression support"),
                       ("ase_support_uid_block_summaries", "Arms with ASE support"), ("concordant_uid_block_summaries", "Arms with concordant support"),
                       ("resolved_arm_summaries", "Resolved arm summaries")], pair_rows) +
                f"<h3>UID summary: {d['unique_uids']:,} distinct UIDs</h3>" +
                table([("status", "UID logical-group status"), ("logical_group_summaries", "Summaries")], uid_status_rows)),
        section("Execution and interpretation", f"<p>CALL schema: {esc(d['call_schema_version'])}. {execution_text}</p><p><strong>BALANCED</strong> means no significant RNA arm event with available branch data; it does not establish DNA balance. Expression outliers are REVIEW / NO_CALL, even if the confidence tier says REVIEW_DISCORDANT. RNA dosage is not a direct DNA copy-number measurement; ASE alone does not resolve loss versus gain, and expression alone does not establish donor origin. Cell state, ambient RNA and donor-pair effects can alter apparent signals. Orthogonal DNA or matched ATAC dosage validation is outside this report.</p>"),
    ]
    cards_html = "".join(f"<div class='card'><strong>{fmt(value)}</strong><span>{esc(label)}</span></div>" for label, value in cards)
    return ("<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
            "<title>Hybrid arm findings | Interpretive report</title><style>" + styles + "</style></head><body>"
            "<header><div class='kicker'>Tetraploid hybrid arm analysis</div><h1>Interpretive run summary</h1>"
            f"<p class='lede'>{esc(headline)}</p><p class='caution'>{esc(explanation)}</p>"
            f"<div class='cards'>{cards_html}</div></header>" + "".join(sections) + "</body></html>\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True, help="Absolute completed hybrid run root")
    parser.add_argument("--output-dir", type=Path, help="Absolute output directory; defaults to RUN/report_interpretive")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = parser.parse_args(argv)
    if not args.run_root.is_absolute() or (args.output_dir and not args.output_dir.is_absolute()):
        parser.error("--run-root and --output-dir must be absolute paths")
    run = args.run_root
    output = args.output_dir or run / "report_interpretive"
    try:
        result = rollup(run)
        output.mkdir(parents=True, exist_ok=True)
        if list(output.glob(f"{PREFIX}.*")):
            raise FileExistsError(f"Interpretive outputs already exist in {output}; use a new --output-dir to preserve them")
        table_atomic(output / f"{PREFIX}.arm_summary.tsv", list(result["arms"][0]), result["arms"])
        table_atomic(output / f"{PREFIX}.cell_burden.tsv", list(result["cell_burden"][0]), result["cell_burden"])
        if result["donor_pairs"]:
            table_atomic(output / f"{PREFIX}.donor_pair_summary.tsv", list(result["donor_pairs"][0]), result["donor_pairs"])
        if result["uid_support"]:
            table_atomic(output / f"{PREFIX}.uid_support_summary.tsv", list(result["uid_support"][0]), result["uid_support"])
        write_atomic(output / f"{PREFIX}.interpretive_report.html", lambda handle: handle.write(html_report(result)))
        # Written last: presence of this file signals that generation finished successfully.
        write_atomic(output / f"{PREFIX}.interpretive_summary.json", lambda handle: json.dump(result, handle, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, csv.Error, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    d = result["overview"]
    print(f"Wrote {output / (PREFIX + '.interpretive_report.html')}")
    print(f"Observed cells={d['unique_observed_cells']:,}; rows={d['observed_cell_arm_rows']:,}; "
          f"candidate rows={d['candidate_rows']:,}; expression outliers={d['evidence_class'].get('EXPRESSION_OUTLIER', 0):,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
