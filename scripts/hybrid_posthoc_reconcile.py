#!/usr/bin/env python3
"""Post-hoc hybrid reconciliation over independent native individual/species outputs.

This script intentionally does NOT implement a joint likelihood and does NOT reuse
legacy species_regularize.  It reads already-completed native individual outputs,
native species outputs, and post-hoc QC/swap-audit sidecars, then writes a small
single-row reconciliation summary for the library.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
from pathlib import Path
from typing import Dict, List, Tuple

NA = "NA"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--lib", required=True, help="Library, e.g. lib19 or 19")
    p.add_argument("--audit_root", required=True)
    p.add_argument("--quant_root", required=True)
    p.add_argument("--individual_condition", default="IND_RF")
    p.add_argument("--species_condition", default="SP_RF")
    p.add_argument("--fixed_species_condition", default="SP_FIXED_EMPTY")
    p.add_argument("--out", required=True)
    return p.parse_args()


def norm_lib(x: str) -> str:
    x = str(x).strip()
    if x.lower().startswith("lib"):
        return "lib" + str(int(x[3:]))
    return "lib" + str(int(x))


def lib_num(lib: str) -> int:
    return int(norm_lib(lib)[3:])


def read_tsv_first(path: Path) -> Dict[str, str]:
    if not path.exists():
        return {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            return dict(row)
    return {}


def read_tsv_rows(path: Path, max_rows: int | None = None) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    opener = gzip.open if str(path).endswith(".gz") else open
    rows: List[Dict[str, str]] = []
    with opener(path, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
            if max_rows is not None and len(rows) >= max_rows:
                break
    return rows


def file_status(path: Path) -> str:
    if not path.exists():
        return "missing"
    if path.stat().st_size == 0:
        return "empty"
    return "present"


def count_rows(path: Path) -> str:
    if not path.exists():
        return "0"
    opener = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with opener(path, "rt", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i == 0 and ("\t" in line):
                continue
            if line.strip():
                n += 1
    return str(n)


def _expand_species_token(tok: str) -> list[str]:
    tok = str(tok).strip()
    if tok in {"Hy", "Chinobo-mCherry"}:
        return ["B", "C"]
    return [tok] if tok and tok != NA else []


def fold_species_set_value(value: str) -> str:
    raw = str(value or "").strip()
    if not raw or raw == NA:
        return NA
    entries = []
    for entry in raw.split(";"):
        vals = set()
        for part in entry.replace("+", ",").split(","):
            for sp in _expand_species_token(part):
                vals.add(sp)
        if vals:
            entries.append(",".join(sorted(vals)))
    entries = sorted(set(entries))
    return ";".join(entries) if entries else NA


def profile_top(path: Path, n: int = 4) -> str:
    """Best-effort profile summary: return up to n name:value entries.

    Supports common profile layouts with either explicit name/value columns or a
    first two-column TSV.  If parsing fails, returns NA instead of hard-failing.
    """
    if not path.exists():
        return NA
    rows = read_tsv_rows(path)
    if not rows:
        return NA
    candidates: List[Tuple[str, float]] = []
    for r in rows:
        keys = list(r.keys())
        name = None
        val = None
        for k in ("sample", "species", "identity", "indiv", "name"):
            if k in r and r[k]:
                name = r[k]
                break
        for k in ("contam_prof", "profile", "pi", "fraction", "frac", "mean"):
            if k in r and r[k] not in ("", ".", "NA"):
                try:
                    val = float(r[k])
                    break
                except ValueError:
                    pass
        if name is None and len(keys) >= 1:
            name = r.get(keys[0], "")
        if val is None and len(keys) >= 2:
            try:
                val = float(r.get(keys[1], "nan"))
            except ValueError:
                val = None
        if name and val is not None:
            candidates.append((name, val))
    if not candidates:
        return NA
    candidates.sort(key=lambda x: x[1], reverse=True)
    return ",".join(f"{k}:{v:.6g}" for k, v in candidates[:n])


def main() -> int:
    args = parse_args()
    lib = norm_lib(args.lib)
    n = lib_num(lib)
    audit_root = Path(args.audit_root)
    quant_root = Path(args.quant_root)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    audit_dir = audit_root / lib
    swap_report = audit_dir / f"{lib}.swap_report.tsv"
    species_qc = audit_dir / f"{lib}.species_qc.tsv"
    call_qc = audit_dir / f"{lib}.call_qc.tsv.gz"

    indiv_prefix = quant_root / args.individual_condition / lib / f"lib{n}_demuxed"
    sp_prefix = quant_root / args.species_condition / lib / f"lib{n}_demuxed"
    fixed_prefix = quant_root / args.fixed_species_condition / lib / f"lib{n}_demuxed" if args.fixed_species_condition else None

    swap = read_tsv_first(swap_report)
    spqc = read_tsv_first(species_qc)

    row: Dict[str, str] = {
        "library": lib,
        "hybrid_method": "HYB_POSTHOC_RECONCILIATION_NOT_JOINT_LIKELIHOOD",
        "separation_note": "individual_native_and_species_native_outputs_read_separately_no_species_regularize_no_set_species_mode",
        "individual_condition": args.individual_condition,
        "species_condition": args.species_condition,
        "fixed_species_condition": args.fixed_species_condition or NA,
        "audit_verdict": swap.get("audit_verdict", NA),
        "audit_flags": swap.get("audit_flags", NA),
        "feature_mode": swap.get("feature_mode", NA),
        "n_cells_swap_report": swap.get("n_cells", NA),
        "median_dosage_concordance_assigned": swap.get("median_dosage_concordance_assigned", NA),
        "median_dosage_gap_constrained": swap.get("median_dosage_gap_constrained", NA),
        "frac_cells_low_evidence": swap.get("frac_cells_low_evidence", NA),
        "frac_cells_neg_gap_constrained": swap.get("frac_cells_neg_gap_constrained", NA),
        "median_species_support_expected": spqc.get("median_species_support_expected", swap.get("median_species_support_expected", NA)),
        "frac_cells_species_conflict": spqc.get("frac_cells_species_conflict", swap.get("frac_cells_species_conflict", NA)),
        "frac_cells_species_exact_match": spqc.get("frac_cells_species_exact_match", swap.get("frac_cells_species_exact_match", NA)),
        "frac_cells_species_component_missing": spqc.get("frac_cells_species_component_missing", swap.get("frac_cells_species_component_missing", NA)),
        "frac_cells_species_unexpected_extra": spqc.get("frac_cells_species_unexpected_extra", swap.get("frac_cells_species_unexpected_extra", NA)),
        "frac_cells_species_partial_overlap_extra": spqc.get("frac_cells_species_partial_overlap_extra", swap.get("frac_cells_species_partial_overlap_extra", NA)),
        "frac_cells_species_disjoint_wrong": spqc.get("frac_cells_species_disjoint_wrong", swap.get("frac_cells_species_disjoint_wrong", NA)),
        "frac_cells_species_unexpected_or_disjoint": spqc.get("frac_cells_species_unexpected_or_disjoint", swap.get("frac_cells_species_unexpected_or_disjoint", NA)),
        "expected_species_set": fold_species_set_value(swap.get("species_expected") or spqc.get("expected_species_set", NA)),
        "observed_species_evidence": spqc.get("observed_species_evidence", NA),
        "call_qc_status": file_status(call_qc),
        "call_qc_rows": count_rows(call_qc),
        "species_qc_status": file_status(species_qc),
        "swap_report_status": file_status(swap_report),
        "individual_contam_rate_status": file_status(Path(str(indiv_prefix) + ".contam_rate")),
        "individual_contam_prof_status": file_status(Path(str(indiv_prefix) + ".contam_prof")),
        "species_contam_rate_status": file_status(Path(str(sp_prefix) + ".contam_rate")),
        "species_prof_status": file_status(Path(str(sp_prefix) + ".species_prof")),
        "fixed_species_contam_rate_status": file_status(Path(str(fixed_prefix) + ".contam_rate")) if fixed_prefix else NA,
        "fixed_species_prof_status": file_status(Path(str(fixed_prefix) + ".species_prof")) if fixed_prefix else NA,
        "individual_profile_top": profile_top(Path(str(indiv_prefix) + ".contam_prof")),
        "species_profile_top": profile_top(Path(str(sp_prefix) + ".species_prof")),
        "fixed_species_profile_top": profile_top(Path(str(fixed_prefix) + ".species_prof")) if fixed_prefix else NA,
    }

    fields = list(row.keys())
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
