#!/usr/bin/env python3
"""Cross-check POSTHOC unexpected genotype components against NN ploidy calls.

This is analysis-only. It never changes demux, NN, tetra_refine, or POSTHOC calls.

For each library, the script reads the current swap report to identify the top
supported unexpected component, then selects the exact POSTHOC cells whose
unconstrained genotype call contains that component and whose call-QC is not
LOW_EVIDENCE. Those barcodes are joined to the existing NN ploidy calls.

The key biological check is whether an unexpected component that appears as a
single-genotype audit call is expression-NN diploid-like (low P(tet)), while
unexpected-component fusion calls are allowed to be tetraploid-like.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
from collections import Counter
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

NA = "NA"


def open_text(path: Path):
    return gzip.open(path, "rt", newline="") if str(path).endswith(".gz") else open(path, "r", newline="")


def read_tsv(path: Path) -> List[dict]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_one_row(path: Path) -> dict:
    rows = read_tsv(path)
    if not rows:
        raise ValueError(f"No data rows in {path}")
    return rows[0]


def finite_float(value) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.upper() in {"NA", "NAN", "NONE", "."}:
        return None
    try:
        x = float(text)
    except ValueError:
        return None
    return x if math.isfinite(x) else None


def pct(n: int, d: int) -> float:
    return n / d if d else float("nan")


def fmt(value):
    if value is None:
        return NA
    if isinstance(value, float):
        if not math.isfinite(value):
            return NA
        return f"{value:.10g}"
    return str(value)


def quantile(values: Sequence[float], q: float) -> float:
    vals = sorted(float(x) for x in values if x is not None and math.isfinite(float(x)))
    if not vals:
        return float("nan")
    if len(vals) == 1:
        return vals[0]
    pos = q * (len(vals) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    w = pos - lo
    return vals[lo] * (1.0 - w) + vals[hi] * w


def identity_components(identity: str) -> List[str]:
    if identity is None:
        return []
    text = str(identity).strip()
    if not text or text.upper() in {"NA", "NONE", "."}:
        return []
    return [x.strip() for x in text.split("+") if x.strip()]


def distinct_components(identity: str) -> List[str]:
    return sorted(set(identity_components(identity)))


def classify_audit_identity(identity: str, unexpected_component: str) -> str:
    comps = identity_components(identity)
    uniq = sorted(set(comps))
    if not uniq:
        return "missing_audit_identity"
    if len(uniq) == 1:
        if unexpected_component in uniq:
            return "unexpected_single_genotype"
        return "single_genotype_other"
    if unexpected_component in uniq:
        return "unexpected_heterotypic_fusion"
    return "heterotypic_other"


def parse_libs(text: str) -> List[int]:
    out = set()
    for part in str(text).replace(" ", ",").split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            lo, hi = int(a), int(b)
            if hi < lo:
                lo, hi = hi, lo
            out.update(range(lo, hi + 1))
        else:
            out.add(int(part))
    return sorted(out)


def load_call_qc(path: Path) -> Dict[str, dict]:
    out = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"barcode", "call_qc_flags"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing call-QC columns {sorted(missing)}")
        for row in reader:
            out[row["barcode"]] = row
    return out


def load_nn(path: Path) -> Dict[str, dict]:
    out = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"barcode", "prob_tetraploid"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing NN columns {sorted(missing)}")
        for row in reader:
            out[row["barcode"]] = row
    return out


def summarize_probs(rows: Sequence[dict], prefix: str) -> dict:
    probs = [r["nn_prob_tetraploid"] for r in rows if r.get("nn_prob_tetraploid") is not None]
    n = len(rows)
    m = len(probs)
    result = {
        f"{prefix}_n": n,
        f"{prefix}_nn_matched": m,
        f"{prefix}_nn_match_fraction": pct(m, n),
        f"{prefix}_ptet_q05": quantile(probs, 0.05),
        f"{prefix}_ptet_q25": quantile(probs, 0.25),
        f"{prefix}_ptet_median": quantile(probs, 0.50),
        f"{prefix}_ptet_q75": quantile(probs, 0.75),
        f"{prefix}_ptet_q95": quantile(probs, 0.95),
        f"{prefix}_ptet_lt_0_10_fraction": pct(sum(p < 0.10 for p in probs), m),
        f"{prefix}_ptet_lt_0_50_fraction": pct(sum(p < 0.50 for p in probs), m),
        f"{prefix}_ptet_ge_0_90_fraction": pct(sum(p >= 0.90 for p in probs), m),
    }
    qc_values = [str(r.get("nn_qc_pass", "")).strip().lower() for r in rows if r.get("nn_qc_pass") not in (None, "")]
    result[f"{prefix}_nn_qc_pass_fraction"] = pct(sum(x in {"1", "true", "t", "yes"} for x in qc_values), len(qc_values))
    return result


def write_tsv(path: Path, rows: Sequence[dict], fields: Sequence[str], gzip_output: bool = False):
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if gzip_output else open
    mode = "wt" if gzip_output else "w"
    with opener(path, mode, newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({f: fmt(row.get(f)) for f in fields})


def analyze_library(lib: str, audit_root: Path, nn_root: Path, out_root: Path) -> Tuple[dict, List[dict]]:
    libdir = audit_root / lib
    report_path = libdir / f"{lib}.swap_report.tsv"
    scores_path = libdir / f"{lib}.swap_scores.tsv"
    qc_path = libdir / f"{lib}.call_qc.tsv.gz"
    nn_path = nn_root / f"{lib}.ploidy_calls_nn.tsv"

    for path in (report_path, scores_path, qc_path, nn_path):
        if not path.is_file():
            raise FileNotFoundError(f"{lib}: required input missing: {path}")

    report = read_one_row(report_path)
    component = str(report.get("top_unexpected_component_supported", "")).strip()
    verdict = str(report.get("audit_verdict", "")).strip()
    reported_fraction = finite_float(report.get("top_unexpected_component_fraction_supported"))

    call_qc = load_call_qc(qc_path)
    nn = load_nn(nn_path)

    # Reconstruct the full set of individuals expected anywhere in this library
    # from the swap report.  The current swap_scores.tsv schema reserves an
    # unexpected_pool_components column, but older/current summary outputs may
    # leave that per-cell field blank even though the library-level counter was
    # computed correctly.  Recomputing from audit_unconstrained_best exactly
    # matches the library-level definition and avoids relying on that blank field.
    expected_pool_components = set()
    expected_identity_text = str(report.get("expected_identity", "") or "").strip()
    for expected_identity in expected_identity_text.split(","):
        expected_identity = expected_identity.strip()
        if expected_identity:
            expected_pool_components.update(identity_components(expected_identity))

    score_rows = read_tsv(scores_path)
    implicated: List[dict] = []
    total_supported = 0

    for row in score_rows:
        bc = row.get("barcode", "")
        qc = call_qc.get(bc, {})
        flags = str(qc.get("call_qc_flags", ""))
        supported = "LOW_EVIDENCE" not in flags
        if supported:
            total_supported += 1

        audit_best = row.get("audit_unconstrained_best", "")
        pool_unexpected = set(identity_components(audit_best)) - expected_pool_components
        if not component or component.upper() == "NA" or component not in pool_unexpected or not supported:
            continue

        nnrow = nn.get(bc, {})
        p_tet = finite_float(nnrow.get("prob_tetraploid"))
        category = classify_audit_identity(audit_best, component)
        implicated.append({
            "library": lib,
            "barcode": bc,
            "unexpected_component": component,
            "audit_unconstrained_best": audit_best,
            "audit_identity_class": category,
            "raw_demux_assignment": row.get("raw_demux_assignment", row.get("original_assignment", "")),
            "refined_biological_assignment": row.get("refined_biological_assignment", ""),
            "refined_assignment_source": row.get("refined_assignment_source", ""),
            "dosage_concordance_assigned": finite_float(row.get("dosage_concordance_assigned")),
            "delta_ll_best_vs_expected": finite_float(row.get("delta_ll_best_vs_expected")),
            "call_qc_flags": flags,
            "nn_call": nnrow.get("ploidy_call", ""),
            "nn_ploidy_probability": finite_float(nnrow.get("ploidy_probability")),
            "nn_prob_tetraploid": p_tet,
            "nn_qc_pass": nnrow.get("qc_pass", ""),
        })

    # Same-library background = all NN calls other than the implicated set.
    implicated_bcs = {r["barcode"] for r in implicated}
    background = []
    for bc, nnrow in nn.items():
        if bc in implicated_bcs:
            continue
        background.append({
            "nn_prob_tetraploid": finite_float(nnrow.get("prob_tetraploid")),
            "nn_qc_pass": nnrow.get("qc_pass", ""),
        })

    single = [r for r in implicated if r["audit_identity_class"] == "unexpected_single_genotype"]
    fusion = [r for r in implicated if r["audit_identity_class"] == "unexpected_heterotypic_fusion"]

    summary = {
        "library": lib,
        "audit_verdict": verdict,
        "unexpected_component": component or NA,
        "reported_top_supported_fraction": reported_fraction,
        "n_swap_score_cells": len(score_rows),
        "n_supported_swap_cells": total_supported,
        "n_implicated_supported_cells": len(implicated),
        "recomputed_implicated_fraction_all_audit_cells": pct(len(implicated), len(score_rows)),
        "recomputed_implicated_fraction_supported_audit_cells": pct(len(implicated), total_supported),
        "audit_single_genotype_n": len(single),
        "audit_heterotypic_fusion_n": len(fusion),
        "audit_other_n": len(implicated) - len(single) - len(fusion),
    }
    summary.update(summarize_probs(implicated, "implicated"))
    summary.update(summarize_probs(single, "single_genotype"))
    summary.update(summarize_probs(fusion, "heterotypic_fusion"))
    summary.update(summarize_probs(background, "library_background"))

    # A compact interpretation field. This is descriptive, not a hard gate.
    single_probs = [r["nn_prob_tetraploid"] for r in single if r.get("nn_prob_tetraploid") is not None]
    if not component or component.upper() == "NA":
        interpretation = "NO_UNEXPECTED_COMPONENT"
    elif len(single_probs) >= 10 and pct(sum(p < 0.50 for p in single_probs), len(single_probs)) >= 0.80:
        interpretation = "FOREIGN_SINGLE_GENOTYPE_IS_NN_DIPLOID_LIKE"
    elif len(single_probs) >= 10 and pct(sum(p >= 0.90 for p in single_probs), len(single_probs)) >= 0.80:
        interpretation = "FOREIGN_SINGLE_GENOTYPE_IS_NN_TET_LIKE"
    elif len(single_probs) >= 10:
        interpretation = "FOREIGN_SINGLE_GENOTYPE_HAS_MIXED_NN_PLOIDY"
    elif fusion:
        interpretation = "SIGNAL_DOMINATED_BY_FUSION_CALLS"
    else:
        interpretation = "TOO_FEW_SINGLE_GENOTYPE_CELLS"
    summary["nn_crosscheck_interpretation"] = interpretation

    outdir = out_root / lib
    cell_fields = [
        "library", "barcode", "unexpected_component", "audit_unconstrained_best",
        "audit_identity_class", "raw_demux_assignment", "refined_biological_assignment",
        "refined_assignment_source", "dosage_concordance_assigned",
        "delta_ll_best_vs_expected", "call_qc_flags", "nn_call",
        "nn_ploidy_probability", "nn_prob_tetraploid", "nn_qc_pass",
    ]
    write_tsv(outdir / f"{lib}.unexpected_component_nn_cells.tsv.gz", implicated, cell_fields, gzip_output=True)
    return summary, implicated


SUMMARY_FIELDS = [
    "library", "audit_verdict", "unexpected_component", "reported_top_supported_fraction",
    "n_swap_score_cells", "n_supported_swap_cells", "n_implicated_supported_cells",
    "recomputed_implicated_fraction_all_audit_cells", "recomputed_implicated_fraction_supported_audit_cells",
    "audit_single_genotype_n", "audit_heterotypic_fusion_n", "audit_other_n",
    "implicated_n", "implicated_nn_matched", "implicated_nn_match_fraction",
    "implicated_ptet_q05", "implicated_ptet_q25", "implicated_ptet_median", "implicated_ptet_q75", "implicated_ptet_q95",
    "implicated_ptet_lt_0_10_fraction", "implicated_ptet_lt_0_50_fraction", "implicated_ptet_ge_0_90_fraction", "implicated_nn_qc_pass_fraction",
    "single_genotype_n", "single_genotype_nn_matched", "single_genotype_nn_match_fraction",
    "single_genotype_ptet_q05", "single_genotype_ptet_q25", "single_genotype_ptet_median", "single_genotype_ptet_q75", "single_genotype_ptet_q95",
    "single_genotype_ptet_lt_0_10_fraction", "single_genotype_ptet_lt_0_50_fraction", "single_genotype_ptet_ge_0_90_fraction", "single_genotype_nn_qc_pass_fraction",
    "heterotypic_fusion_n", "heterotypic_fusion_nn_matched", "heterotypic_fusion_nn_match_fraction",
    "heterotypic_fusion_ptet_q05", "heterotypic_fusion_ptet_q25", "heterotypic_fusion_ptet_median", "heterotypic_fusion_ptet_q75", "heterotypic_fusion_ptet_q95",
    "heterotypic_fusion_ptet_lt_0_10_fraction", "heterotypic_fusion_ptet_lt_0_50_fraction", "heterotypic_fusion_ptet_ge_0_90_fraction", "heterotypic_fusion_nn_qc_pass_fraction",
    "library_background_n", "library_background_nn_matched", "library_background_nn_match_fraction",
    "library_background_ptet_q05", "library_background_ptet_q25", "library_background_ptet_median", "library_background_ptet_q75", "library_background_ptet_q95",
    "library_background_ptet_lt_0_10_fraction", "library_background_ptet_lt_0_50_fraction", "library_background_ptet_ge_0_90_fraction", "library_background_nn_qc_pass_fraction",
    "nn_crosscheck_interpretation",
]


def write_library_summary(out_root: Path, lib: str, summary: dict):
    outdir = out_root / lib
    write_tsv(outdir / f"{lib}.unexpected_component_nn_summary.tsv", [summary], SUMMARY_FIELDS)


def aggregate(out_root: Path, libs: Sequence[int]):
    rows = []
    all_cells: List[dict] = []
    cell_fields = [
        "library", "barcode", "unexpected_component", "audit_unconstrained_best",
        "audit_identity_class", "raw_demux_assignment", "refined_biological_assignment",
        "refined_assignment_source", "dosage_concordance_assigned",
        "delta_ll_best_vs_expected", "call_qc_flags", "nn_call",
        "nn_ploidy_probability", "nn_prob_tetraploid", "nn_qc_pass",
    ]
    for n in libs:
        lib = f"lib{n}"
        summary_path = out_root / lib / f"{lib}.unexpected_component_nn_summary.tsv"
        cells_path = out_root / lib / f"{lib}.unexpected_component_nn_cells.tsv.gz"
        if not summary_path.is_file():
            continue
        rows.extend(read_tsv(summary_path))
        if cells_path.is_file():
            all_cells.extend(read_tsv(cells_path))

    write_tsv(out_root / "all_libraries.unexpected_component_nn_summary.tsv", rows, SUMMARY_FIELDS)
    write_tsv(out_root / "all_libraries.unexpected_component_nn_cells.tsv.gz", all_cells, cell_fields, gzip_output=True)

    flagged = [r for r in rows if r.get("unexpected_component") not in (None, "", NA) and r.get("audit_verdict") == "UNEXPECTED_COMPONENT_SIGNAL"]
    with open(out_root / "summary.txt", "w") as out:
        out.write("UNEXPECTED COMPONENT x NN PLOIDY CROSS-CHECK\n")
        out.write("=============================================\n\n")
        out.write(f"Libraries summarized: {len(rows)}\n")
        out.write(f"Libraries with UNEXPECTED_COMPONENT_SIGNAL: {len(flagged)}\n\n")
        for r in flagged:
            out.write(
                f"{r.get('library')}  {r.get('unexpected_component')}  "
                f"supported_fraction={r.get('reported_top_supported_fraction')}  "
                f"n={r.get('n_implicated_supported_cells')}  "
                f"single_n={r.get('audit_single_genotype_n')}  "
                f"single_median_Ptet={r.get('single_genotype_ptet_median')}  "
                f"single_Ptet<0.5={r.get('single_genotype_ptet_lt_0_50_fraction')}  "
                f"fusion_n={r.get('audit_heterotypic_fusion_n')}  "
                f"fusion_median_Ptet={r.get('heterotypic_fusion_ptet_median')}  "
                f"interpretation={r.get('nn_crosscheck_interpretation')}\n"
            )

    print(f"Wrote aggregate NN cross-check: {out_root / 'all_libraries.unexpected_component_nn_summary.tsv'}")
    print(f"Wrote cell-level implicated set: {out_root / 'all_libraries.unexpected_component_nn_cells.tsv.gz'}")
    print(f"Wrote compact summary: {out_root / 'summary.txt'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lib", help="Library, e.g. lib11")
    ap.add_argument("--audit-root", type=Path, required=True)
    ap.add_argument("--nn-root", type=Path, required=True)
    ap.add_argument("--out-root", type=Path, required=True)
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("--libraries", default="1-40")
    args = ap.parse_args()

    args.out_root.mkdir(parents=True, exist_ok=True)

    if args.aggregate:
        aggregate(args.out_root, parse_libs(args.libraries))
        return

    if not args.lib:
        ap.error("--lib is required unless --aggregate is used")

    lib = args.lib if str(args.lib).startswith("lib") else f"lib{args.lib}"
    print(f"[{lib}] joining POSTHOC unexpected-component cells to NN ploidy calls", flush=True)
    summary, implicated = analyze_library(lib, args.audit_root, args.nn_root, args.out_root)
    write_library_summary(args.out_root, lib, summary)
    print(
        f"[{lib}] component={summary['unexpected_component']} verdict={summary['audit_verdict']} "
        f"implicated={len(implicated)} single={summary['audit_single_genotype_n']} "
        f"fusion={summary['audit_heterotypic_fusion_n']} "
        f"single_median_Ptet={fmt(summary.get('single_genotype_ptet_median'))}",
        flush=True,
    )


if __name__ == "__main__":
    main()
