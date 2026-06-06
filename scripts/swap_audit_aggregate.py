#!/usr/bin/env python3
"""Aggregate per-library swap-audit reports."""
from __future__ import annotations

import argparse
import csv
import glob
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List

from swap_audit_lib import NA, identity_type, open_text, read_refined_assignments, refined_assignment_counts, read_tsv, write_tsv


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--audit_root", required=True)
    p.add_argument("--refined_assignments_root", default=None)
    p.add_argument("--out_prefix", default=None, help="Default: <audit_root>/all_libraries")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--skip-existing", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def union_fields(rows: List[Dict[str, str]]) -> List[str]:
    preferred = [
        "library", "preflight_status", "audit_verdict", "audit_flags", "feature_mode",
        "has_call_qc", "has_species_qc", "has_atac_qc", "has_calibrated_thresholds",
        "missing_optional_features", "warnings", "thresholds_source",
    ]
    seen = []
    for f in preferred:
        if any(f in r for r in rows): seen.append(f)
    for r in rows:
        for f in r:
            if f not in seen: seen.append(f)
    return seen


def load_refined_counts(path: Path) -> Dict[str, int] | None:
    """Load tetra_refine output and count refined_assignment, not original_assignment."""
    if not path.exists():
        return None
    refined = read_refined_assignments(path)
    if not refined:
        return Counter()
    return refined_assignment_counts(refined)


def main() -> int:
    args = parse_args()
    root = Path(args.audit_root).resolve()
    out_prefix = Path(args.out_prefix).resolve() if args.out_prefix else root / "all_libraries"
    paths = sorted(root.glob("lib*/lib*.swap_report.tsv"))
    rows: List[Dict[str, str]] = []
    for p in paths:
        try:
            rows.extend(read_tsv(p))
        except Exception as e:
            rows.append({"library": p.parent.name, "preflight_status": "FAIL_READ_REPORT", "audit_verdict": "NOT_RUN", "warnings": f"FAIL_READ_REPORT:{e}"})
    if not rows:
        raise SystemExit(f"ERROR: no per-library reports found under {root}/lib*/lib*.swap_report.tsv")
    fields = union_fields(rows)
    write_tsv(str(out_prefix) + ".swap_audit.tsv", rows, fields)
    flagged = [r for r in rows if r.get("audit_verdict") != "OK" or (r.get("audit_flags") or "")]
    write_tsv(str(out_prefix) + ".flagged.tsv", flagged, fields)

    cap_fields = ["library", "feature_mode", "has_call_qc", "has_species_qc", "has_atac_qc", "has_refined_assignments", "has_calibrated_thresholds", "missing_optional_features", "warnings"]
    write_tsv(str(out_prefix) + ".capability_summary.tsv", rows, cap_fields)

    # Expected vs observed/audit-best matrix across libraries.
    all_cols = sorted({r.get("audit_best_identity_unconstrained", NA) for r in rows if r.get("audit_best_identity_unconstrained")})
    mat = []
    for r in rows:
        mr = {"library": r.get("library", NA), "expected_identity": r.get("expected_identity", NA)}
        for c in all_cols:
            mr[c] = 1 if r.get("audit_best_identity_unconstrained") == c else 0
        mat.append(mr)
    write_tsv(str(out_prefix) + ".expected_vs_observed_matrix.tsv", mat, ["library", "expected_identity"] + all_cols)

    summary_rows = []
    for r in rows:
        lib = r.get("library", NA)
        warn = [x for x in (r.get("warnings", "") or "").split(",") if x]
        ref_counts = None
        if args.refined_assignments_root and lib != NA:
            base = Path(args.refined_assignments_root)
            candidates = [base / lib / f"{lib}.refined_assignments", base / f"{lib}.refined_assignments", base / f"{lib}_demuxed.refined_assignments"]
            for c in candidates:
                ref_counts = load_refined_counts(c)
                if ref_counts is not None: break
        if ref_counts is None:
            warn.append("REFINED_ASSIGNMENTS_MISSING_COMPOSITION_SKIPPED")
        summary_rows.append({
            "library": lib,
            "assignment_source": "refined" if ref_counts is not None else "original_demux",
            "n_cells": r.get("n_cells", NA),
            "median_dosage_concordance": r.get("median_dosage_concordance_assigned", NA),
            "median_dosage_gap_constrained": r.get("median_dosage_gap_constrained", NA),
            "frac_cells_low_evidence": r.get("frac_cells_low_evidence", NA),
            "frac_cells_neg_gap_constrained": r.get("frac_cells_neg_gap_constrained", NA),
            "frac_cells_species_conflict": r.get("frac_cells_species_conflict", NA),
            "frac_cells_species_exact_match": r.get("frac_cells_species_exact_match", NA),
            "frac_cells_species_component_missing": r.get("frac_cells_species_component_missing", NA),
            "frac_cells_species_unexpected_extra": r.get("frac_cells_species_unexpected_extra", NA),
            "frac_cells_species_partial_overlap_extra": r.get("frac_cells_species_partial_overlap_extra", NA),
            "frac_cells_species_disjoint_wrong": r.get("frac_cells_species_disjoint_wrong", NA),
            "frac_cells_species_unexpected_or_disjoint": r.get("frac_cells_species_unexpected_or_disjoint", NA),
            "median_species_support_expected": r.get("median_species_support_expected", NA),
            "median_total_depth": r.get("median_total_depth", NA),
            "n_refined_cells": r.get("refined_n_cells", ref_counts.get("total", NA) if ref_counts else NA),
            "n_refined_total_available": r.get("refined_n_total_available", NA),
            "n_refined_overlap_with_audit": r.get("refined_n_overlap_with_audit", NA),
            "refined_overlap_fraction": r.get("refined_overlap_fraction", NA),
            "n_refined_singlet": r.get("refined_n_singlet", ref_counts.get("singlet", NA) if ref_counts else NA),
            "n_refined_heterotypic": r.get("refined_n_heterotypic", ref_counts.get("heterotypic", NA) if ref_counts else NA),
            "n_refined_homotypic": r.get("refined_n_homotypic", ref_counts.get("homotypic", NA) if ref_counts else NA),
            "n_refined_plus_identity": r.get("refined_n_plus_identity", ref_counts.get("plus_identity", NA) if ref_counts else NA),
            "n_refined_biological_fusion": r.get("refined_n_biological_fusion", ref_counts.get("biological_fusion", NA) if ref_counts else NA),
            "n_refined_doublet": r.get("refined_n_doublet", ref_counts.get("doublet", NA) if ref_counts else NA),
            "n_refined_changed": r.get("refined_n_changed", ref_counts.get("changed", NA) if ref_counts else NA),
            "n_refined_droplet_flagged": r.get("refined_n_droplet_flagged", ref_counts.get("droplet_flagged", NA) if ref_counts else NA),
            "warnings": sorted(set(warn)),
        })
    write_tsv(str(out_prefix) + ".call_qc_summary.tsv", summary_rows, [
        "library", "assignment_source", "n_cells", "median_dosage_concordance", "median_dosage_gap_constrained",
        "frac_cells_low_evidence", "frac_cells_neg_gap_constrained", "frac_cells_species_conflict",
        "frac_cells_species_exact_match", "frac_cells_species_component_missing",
        "frac_cells_species_unexpected_extra", "frac_cells_species_partial_overlap_extra",
        "frac_cells_species_disjoint_wrong", "frac_cells_species_unexpected_or_disjoint",
        "median_species_support_expected", "median_total_depth",
        "n_refined_cells", "n_refined_total_available", "n_refined_overlap_with_audit", "refined_overlap_fraction",
        "n_refined_singlet", "n_refined_heterotypic", "n_refined_homotypic",
        "n_refined_plus_identity", "n_refined_biological_fusion", "n_refined_doublet",
        "n_refined_changed", "n_refined_droplet_flagged", "warnings"
    ])
    print(f"Wrote {out_prefix}.swap_audit.tsv")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
