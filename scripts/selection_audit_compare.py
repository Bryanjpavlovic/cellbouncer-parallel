#!/usr/bin/env python3
"""Compare production maximin selection with the audit-only max-LLR comparator.

The comparator reported by ``demux_parallel --dump_selection_audit`` is the
non-mutating statistic ``argmax(maxllr)``. It is not the historical destructive
process-of-elimination selector and this script makes no parity claim with that
historical implementation.

Required diagnostic columns:
    maximin_winner
    maximin_score
    max_llr_comparator_winner
    max_llr_comparator_score
    selection_agree

``margin_softmax_score`` is used only as a descriptive confidence summary. It is
not interpreted as a calibrated Bayesian posterior probability.
"""
from __future__ import annotations

import argparse
import gzip
import os
from collections import Counter
from typing import IO, Dict, Iterable, List, Optional


AUDIT_COLS = [
    "maximin_winner",
    "maximin_score",
    "max_llr_comparator_winner",
    "max_llr_comparator_score",
    "selection_agree",
]


def _open(path: str) -> IO[str]:
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def _parse_optional_float(value: str) -> Optional[float]:
    if value in {"", ".", "NA", "nan", "NaN"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def find_audited(root: str) -> List[str]:
    """Find diagnostics carrying the corrected audit schema under ``root``."""
    hits: List[str] = []
    for dirpath, _dirs, files in os.walk(root):
        for filename in files:
            if not filename.endswith(".diagnostics.gz"):
                continue
            path = os.path.join(dirpath, filename)
            try:
                with _open(path) as handle:
                    header = handle.readline().rstrip("\n").split("\t")
                if all(column in header for column in AUDIT_COLS):
                    hits.append(path)
            except OSError:
                continue
    return sorted(hits)


def _confidence_band(
    margin_softmax_score: Optional[float], min_margin: Optional[float]
) -> str:
    """Return a descriptive band, preferring the margin-softmax summary."""
    if margin_softmax_score is not None:
        if margin_softmax_score >= 0.999:
            return "margin_softmax>=0.999"
        if margin_softmax_score >= 0.99:
            return "margin_softmax 0.99-0.999"
        if margin_softmax_score >= 0.9:
            return "margin_softmax 0.9-0.99"
        return "margin_softmax<0.9"
    if min_margin is None:
        return "unknown"
    if min_margin >= 100:
        return "min_margin>=100"
    if min_margin >= 10:
        return "min_margin 10-100"
    if min_margin >= 1:
        return "min_margin 1-10"
    return "min_margin<1"


def process_file(
    path: str,
    aggregate: Dict[str, object],
    disagreement_writer: Optional[IO[str]] = None,
    source_name: Optional[str] = None,
) -> None:
    """Update counters from one corrected-schema diagnostics file."""
    with _open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {column: position for position, column in enumerate(header)}
        missing = [column for column in AUDIT_COLS if column not in index]
        if missing:
            raise SystemExit(
                f"ERROR: {path} lacks corrected audit columns {missing}. "
                "Re-run demux_parallel with --dump_selection_audit."
            )

        assignment_index = index.get("assignment")
        softmax_index = index.get("margin_softmax_score")
        margin_index = index.get("min_margin")
        comparator_index = index["max_llr_comparator_winner"]
        agree_index = index["selection_agree"]

        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= agree_index:
                continue
            agree_text = fields[agree_index]
            if agree_text not in {"0", "1"}:
                # NA means one of the two statistics had no evaluable winner.
                aggregate["n_not_comparable"] += 1
                continue

            agree = int(agree_text)
            aggregate["n_comparable"] += 1
            softmax = (
                _parse_optional_float(fields[softmax_index])
                if softmax_index is not None and softmax_index < len(fields)
                else None
            )
            min_margin = (
                _parse_optional_float(fields[margin_index])
                if margin_index is not None and margin_index < len(fields)
                else None
            )
            band = _confidence_band(softmax, min_margin)
            aggregate["band_total"][band] += 1

            if agree == 0:
                aggregate["n_disagree"] += 1
                aggregate["band_disagree"][band] += 1
                if disagreement_writer is not None:
                    assignment = (
                        fields[assignment_index]
                        if assignment_index is not None and assignment_index < len(fields)
                        else "NA"
                    )
                    comparator = fields[comparator_index]
                    disagreement_writer.write(
                        f"{source_name or os.path.basename(path)}\t{fields[0]}\t"
                        f"{assignment}\t{comparator}\t"
                        f"{'' if softmax is None else f'{softmax:.6f}'}\t"
                        f"{'' if min_margin is None else f'{min_margin:.6f}'}\n"
                    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compare maximin with the non-mutating max_llr_comparator from "
            "audited demux_parallel diagnostics."
        )
    )
    parser.add_argument(
        "diagnostics",
        nargs="*",
        help="One or more *.diagnostics.gz files generated with --dump_selection_audit.",
    )
    parser.add_argument(
        "--root",
        help="Recursively discover corrected-schema audited diagnostics under this directory.",
    )
    parser.add_argument(
        "--out-disagreements",
        help="Optional .tsv or .tsv.gz path for per-cell disagreements.",
    )
    args = parser.parse_args()

    files = list(args.diagnostics)
    if args.root:
        found = find_audited(args.root)
        print(f"Found {len(found)} corrected audited diagnostics file(s) under {args.root}")
        files.extend(found)
    files = sorted(set(files))
    if not files:
        raise SystemExit("ERROR: no input files. Pass paths or --root.")

    aggregate: Dict[str, object] = {
        "n_comparable": 0,
        "n_not_comparable": 0,
        "n_disagree": 0,
        "band_total": Counter(),
        "band_disagree": Counter(),
    }

    writer: Optional[IO[str]] = None
    if args.out_disagreements:
        writer = (
            gzip.open(args.out_disagreements, "wt")
            if args.out_disagreements.endswith(".gz")
            else open(args.out_disagreements, "w")
        )
        writer.write(
            "source\tbarcode\tmaximin_assignment\tmax_llr_comparator_winner\t"
            "margin_softmax_score\tmin_margin\n"
        )

    try:
        for path in files:
            before_comparable = int(aggregate["n_comparable"])
            before_disagree = int(aggregate["n_disagree"])
            before_unavailable = int(aggregate["n_not_comparable"])
            process_file(
                path,
                aggregate,
                disagreement_writer=writer,
                source_name=os.path.basename(path).split(".")[0],
            )
            comparable = int(aggregate["n_comparable"]) - before_comparable
            disagree = int(aggregate["n_disagree"]) - before_disagree
            unavailable = int(aggregate["n_not_comparable"]) - before_unavailable
            rate = 100.0 * disagree / comparable if comparable else 0.0
            print(
                f"  {os.path.basename(path)}: {comparable} comparable cells, "
                f"{disagree} disagreements ({rate:.3f}%), "
                f"{unavailable} not comparable"
            )
    finally:
        if writer is not None:
            writer.close()

    comparable = int(aggregate["n_comparable"])
    disagree = int(aggregate["n_disagree"])
    unavailable = int(aggregate["n_not_comparable"])
    print("\n" + "=" * 68)
    print(f"Comparable cells:                    {comparable}")
    print(f"Comparator unavailable/not comparable: {unavailable}")
    print(
        "Selection disagreements:             "
        f"{disagree} ({(100.0 * disagree / comparable) if comparable else 0.0:.4f}%)"
    )
    print("=" * 68)
    print("\nBy descriptive confidence band (disagreements / comparable cells):")
    band_total: Counter = aggregate["band_total"]
    band_disagree: Counter = aggregate["band_disagree"]
    for band in sorted(band_total, key=lambda value: -band_total[value]):
        total = band_total[band]
        different = band_disagree[band]
        rate = 100.0 * different / total if total else 0.0
        print(f"  {band:32s} {different:>8d} / {total:<8d} ({rate:.3f}%)")

    if args.out_disagreements:
        print(f"\nPer-cell disagreement list written to: {args.out_disagreements}")
    print(
        "\nInterpretation: this compares maximin with argmax(maxllr) only. "
        "margin_softmax_score is a margin-derived descriptive score, not a "
        "calibrated posterior probability."
    )


if __name__ == "__main__":
    main()
