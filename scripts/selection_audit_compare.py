#!/usr/bin/env python3
"""
selection_audit_compare.py
Created in conversation: https://claude.ai/chat/CURRENT_CONVERSATION

Quantify how often the maximin winner-selection criterion (current
demux_parallel: pick the identity with the highest min_margin) disagrees with
the maxllr criterion (Nathan's original demux_vcf del(2) elimination survivor:
pick the identity with the highest maxllr).

Input: one or more .diagnostics.gz files written by `demux_parallel` run with
--dump_selection_audit. That flag appends these columns:
    maximin_winner   the actual assignment (highest min_margin)
    maximin_score    that winner's min_margin
    maxllr_winner    identity the maxllr criterion would pick
    maxllr_score     that identity's maxllr
    selection_agree  1 if the two criteria pick the same identity, else 0

Reports the overall disagreement rate and breaks it down by confidence so you
can see whether disagreements concentrate in low-confidence cells (expected and
harmless) or reach confident calls (would warrant attention). Confidence is
read from the `posterior` column when present, else from min_margin bands.

Usage:
    module purge
    module load miniforge/3
    module load genomics-base
    python3 selection_audit_compare.py \
        /path/to/libN_demuxed.diagnostics.gz [more.diagnostics.gz ...] \
        --out-disagreements disagreements.tsv.gz

    # Or point at a root and let it find audited diagnostics:
    python3 selection_audit_compare.py --root /mnt/beegfs/.../mapping_output

Revision history at bottom of file.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from collections import Counter


AUDIT_COLS = ["maximin_winner", "maximin_score",
              "maxllr_winner", "maxllr_score", "selection_agree"]


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def find_audited(root):
    """Walk root for *.diagnostics.gz that carry the audit columns."""
    hits = []
    for dirpath, _dirs, files in os.walk(root):
        for f in files:
            if f.endswith(".diagnostics.gz"):
                p = os.path.join(dirpath, f)
                try:
                    with _open(p) as fh:
                        header = fh.readline().rstrip("\n").split("\t")
                    if all(c in header for c in AUDIT_COLS):
                        hits.append(p)
                except OSError:
                    continue
    return sorted(hits)


def _conf_band(posterior, min_margin):
    """Coarse confidence band for breakdown. Prefer posterior if available."""
    if posterior is not None:
        if posterior >= 0.999:
            return "posterior>=0.999"
        if posterior >= 0.99:
            return "0.99-0.999"
        if posterior >= 0.9:
            return "0.9-0.99"
        return "<0.9"
    # fall back to min_margin bands
    if min_margin is None:
        return "unknown"
    if min_margin >= 100:
        return "min_margin>=100"
    if min_margin >= 10:
        return "min_margin 10-100"
    if min_margin >= 1:
        return "min_margin 1-10"
    return "min_margin<1"


def process_file(path, agg, disagree_writer=None, libname=None):
    """Update aggregate counters from one audited diagnostics file."""
    with _open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        missing = [c for c in AUDIT_COLS if c not in idx]
        if missing:
            raise SystemExit(
                f"ERROR: {path} lacks audit columns {missing}. "
                "Re-run demux_parallel with --dump_selection_audit.")
        i_assn = idx.get("assignment")
        i_post = idx.get("posterior")
        i_mm = idx.get("min_margin")
        i_mxw = idx["maxllr_winner"]
        i_agree = idx["selection_agree"]

        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= i_agree:
                continue
            agg["n_total"] += 1
            try:
                agree = int(f[i_agree])
            except ValueError:
                continue
            post = None
            if i_post is not None and i_post < len(f):
                try:
                    post = float(f[i_post])
                except ValueError:
                    post = None
            mm = None
            if i_mm is not None and i_mm < len(f):
                try:
                    mm = float(f[i_mm])
                except ValueError:
                    mm = None
            band = _conf_band(post, mm)
            agg["band_total"][band] += 1
            if agree == 0:
                agg["n_disagree"] += 1
                agg["band_disagree"][band] += 1
                if disagree_writer is not None:
                    assn = f[i_assn] if i_assn is not None and i_assn < len(f) else "?"
                    mxw = f[i_mxw] if i_mxw < len(f) else "?"
                    bc = f[0]
                    disagree_writer.write(
                        f"{libname or os.path.basename(path)}\t{bc}\t{assn}\t{mxw}\t"
                        f"{'' if post is None else f'{post:.6f}'}\t"
                        f"{'' if mm is None else f'{mm:.4f}'}\n")


def main():
    ap = argparse.ArgumentParser(
        description="Count maximin-vs-maxllr selection disagreements from "
                    "audited demux_parallel diagnostics.")
    ap.add_argument("diagnostics", nargs="*",
                    help="One or more *.diagnostics.gz files (run with "
                         "--dump_selection_audit).")
    ap.add_argument("--root",
                    help="Recursively find audited *.diagnostics.gz under this "
                         "directory instead of (or in addition to) listing them.")
    ap.add_argument("--out-disagreements",
                    help="Optional path to write the per-cell disagreement list "
                         "(.tsv or .tsv.gz).")
    args = ap.parse_args()

    files = list(args.diagnostics)
    if args.root:
        found = find_audited(args.root)
        print(f"Found {len(found)} audited diagnostics file(s) under {args.root}")
        files.extend(found)
    files = sorted(set(files))
    if not files:
        raise SystemExit("ERROR: no input files. Pass paths or --root.")

    agg = {
        "n_total": 0,
        "n_disagree": 0,
        "band_total": Counter(),
        "band_disagree": Counter(),
    }

    dw = None
    if args.out_disagreements:
        dw = (gzip.open(args.out_disagreements, "wt")
              if args.out_disagreements.endswith(".gz")
              else open(args.out_disagreements, "w"))
        dw.write("source\tbarcode\tmaximin_winner\tmaxllr_winner\tposterior\tmin_margin\n")

    try:
        for p in files:
            n0, d0 = agg["n_total"], agg["n_disagree"]
            process_file(p, agg, disagree_writer=dw,
                         libname=os.path.basename(p).split(".")[0])
            n = agg["n_total"] - n0
            d = agg["n_disagree"] - d0
            rate = (100.0 * d / n) if n else 0.0
            print(f"  {os.path.basename(p)}: {n} cells, {d} disagreements ({rate:.3f}%)")
    finally:
        if dw is not None:
            dw.close()

    total = agg["n_total"]
    dis = agg["n_disagree"]
    print("\n" + "=" * 60)
    print(f"TOTAL cells:          {total}")
    print(f"Selection disagreements: {dis} ({(100.0*dis/total) if total else 0:.4f}%)")
    print("=" * 60)
    print("\nBy confidence band (disagreements / cells in band):")
    bands = sorted(agg["band_total"], key=lambda b: -agg["band_total"][b])
    for b in bands:
        bt = agg["band_total"][b]
        bd = agg["band_disagree"][b]
        rate = (100.0 * bd / bt) if bt else 0.0
        print(f"  {b:24s} {bd:>8d} / {bt:<8d}  ({rate:.3f}%)")
    if args.out_disagreements:
        print(f"\nPer-cell disagreement list written to: {args.out_disagreements}")
    print("\nInterpretation: disagreements concentrated in low-confidence bands "
          "(posterior<0.9 / small min_margin) indicate the two criteria only "
          "differ on already-ambiguous cells. Disagreements in high-confidence "
          "bands would warrant a closer look.")


if __name__ == "__main__":
    main()


# =============================================================================
# Revision History
# =============================================================================
#   V1_R1  Initial implementation. Reads .diagnostics.gz files written with
#          demux_parallel --dump_selection_audit and reports the maximin-vs-
#          maxllr winner disagreement rate overall and by confidence band.
#          Optional per-cell disagreement dump. --root auto-discovers audited
#          files by checking for the audit columns in the header.
# =============================================================================
