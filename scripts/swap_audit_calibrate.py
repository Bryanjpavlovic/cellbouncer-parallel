#!/usr/bin/env python3
"""Conservative placeholder threshold calibrator for swap-audit V1.

This emits an editable recommended_thresholds.tsv seeded from defaults. It is a
safe first implementation: production can run with defaults immediately, then
replace this with perturbation-based calibration once clean-library truth labels
are curated.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

DEFAULTS = {
    "threshold_min_signal_llr": (10.0, "LLR", "Minimum median original demux LLR"),
    "threshold_min_dosage_concordance": (0.7, "fraction", "Minimum median dosage concordance when call_qc is available"),
    "threshold_max_ambiguous_frac": (0.30, "fraction", "Maximum fraction of cells with n_close > 0"),
    "threshold_max_low_evidence_frac": (0.30, "fraction", "Maximum LOW_EVIDENCE fraction when call_qc is available"),
    "threshold_min_cells": (100, "cells", "Minimum cells for library-level interpretation"),
    "threshold_clean_overlap_frac": (0.90, "fraction", "Clean overlap fraction"),
    "threshold_clean_best_frac": (0.30, "fraction", "Minimum dominant audit-best fraction for clean call"),
    "threshold_full_swap_frac": (0.50, "fraction", "Fraction overlap=0 for likely full swap"),
    "threshold_one_component_swap_frac": (0.50, "fraction", "Fraction overlap=1 for likely component swap"),
    "threshold_coherent_alt_frac": (0.25, "fraction", "Dominant alternate fraction for coherent swap"),
    "threshold_ambiguous_delta_ll": (20.0, "LLR", "Small median delta threshold for ambiguous close genotype"),
    "threshold_ambiguous_gap": (0.0, "fraction", "Gap threshold for ambiguous close genotype"),
    "threshold_homotypic_frac": (0.50, "fraction", "High homotypic/unresolved fraction threshold"),
    "threshold_species_conflict_frac": (0.20, "fraction", "Species conflict fraction threshold"),
    "threshold_clean_concordance": (0.70, "fraction", "Clean concordance threshold when call_qc is available"),
}

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--audit_root", required=True)
    p.add_argument("--clean_libraries", default="")
    p.add_argument("--expected_pool_metadata", required=True)
    p.add_argument("--n_replicates", type=int, default=100)
    p.add_argument("--out", required=True)
    args = p.parse_args()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    libs = [x for x in args.clean_libraries.split(",") if x]
    with open(out, "w") as fh:
        fh.write(f"# calibration_date: {time.strftime('%Y-%m-%d')}\n")
        fh.write(f"# audit_root: {Path(args.audit_root).resolve()}\n")
        fh.write(f"# n_clean_libraries: {len(libs)}\n")
        fh.write(f"# n_replicates_per_perturbation: {args.n_replicates}\n")
        fh.write("# status: defaults_seeded; replace with empirical perturbation calibration after clean truth set is finalized\n")
        fh.write("threshold_name\tvalue\tunits\tdescription\n")
        for name, (value, units, desc) in DEFAULTS.items():
            fh.write(f"{name}\t{value}\t{units}\t{desc}\n")
    print(f"Wrote {out}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
