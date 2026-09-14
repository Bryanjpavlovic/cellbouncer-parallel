#!/usr/bin/env python3
"""Project chromosome-arm intervals through HAL into ancestral coordinates."""

from __future__ import annotations

import argparse
import collections
import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from tetra_arm_common import (
    RELEASE,
    atomic_text,
    require_file,
    require_outputs_absent,
    write_json_atomic,
    write_tsv_atomic,
)


QC_FIELDS = (
    "source_genome", "target_genome", "source_contigs", "source_intervals",
    "source_logical_arms", "raw_lifted_intervals", "mapped_logical_arms",
    "target_contigs", "bed_intervals", "mapped_bases", "conflicting_bases",
    "status", "schema_version",
)


def chromosome_from_arm(arm: str) -> str:
    return arm[:-1] if arm.endswith(("p", "q")) else arm


def is_sex_arm(arm: str) -> bool:
    return chromosome_from_arm(arm).lower().removeprefix("chr") in {"x", "y"}


def read_source_arms(path: str, include_sex: bool):
    rows: list[tuple[str, int, int, str]] = []
    names: set[str] = set()
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected BED4")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: noninteger coordinates") from error
            chrom, name = fields[0].strip(), fields[3].strip()
            if not chrom or not name or start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid BED4 interval")
            if not include_sex and is_sex_arm(name):
                continue
            rows.append((chrom, start, end, name))
            names.add(name)
    if not rows:
        raise ValueError(f"source arms BED contains no selected intervals: {path}")
    return rows, names


def read_lifted(path: str, allowed_names: set[str]):
    intervals: dict[str, list[tuple[int, int, str]]] = collections.defaultdict(list)
    row_count = 0
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected lifted BED4")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: noninteger coordinates") from error
            contig, name = fields[0].strip(), fields[3].strip()
            if not contig or start < 0 or end <= start or name not in allowed_names:
                raise ValueError(f"{path}:{line_number}: invalid lifted arm interval")
            intervals[contig].append((start, end, name))
            row_count += 1
    return intervals, row_count


def disjoint_segments(intervals: dict[str, list[tuple[int, int, str]]]):
    output: list[tuple[str, int, int, str]] = []
    conflicting_bases = 0
    for contig in sorted(intervals):
        events: dict[int, collections.Counter[str]] = collections.defaultdict(
            collections.Counter)
        for start, end, arm in intervals[contig]:
            events[start][arm] += 1
            events[end][arm] -= 1
        active: collections.Counter[str] = collections.Counter()
        previous = None
        for position in sorted(events):
            active_arms = [arm for arm, count in active.items() if count > 0]
            if previous is not None and position > previous:
                if len(active_arms) == 1:
                    arm = active_arms[0]
                    if (output and output[-1][0] == contig and
                            output[-1][2] == previous and output[-1][3] == arm):
                        prior = output[-1]
                        output[-1] = (prior[0], prior[1], position, prior[3])
                    else:
                        output.append((contig, previous, position, arm))
                elif len(active_arms) > 1:
                    conflicting_bases += position - previous
            for arm, delta in events[position].items():
                active[arm] += delta
                if active[arm] == 0:
                    del active[arm]
            previous = position
    return output, conflicting_bases


def hal_genomes(hal_stats: str, hal_path: str) -> set[str]:
    result = subprocess.run(
        [hal_stats, hal_path, "--genomes"], check=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    genomes = set(result.stdout.split())
    if not genomes:
        raise ValueError(f"halStats returned no genomes for {hal_path}")
    return genomes


def hal_sequence_lengths(
        hal_stats: str, hal_path: str, genome: str) -> dict[str, int]:
    result = subprocess.run(
        [hal_stats, hal_path, "--sequenceStats", genome], check=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    lengths: dict[str, int] = {}
    for row in csv.reader(result.stdout.splitlines()):
        if not row or row[0] == "SequenceName":
            continue
        if len(row) < 2:
            raise ValueError(f"malformed halStats sequence row for {genome}")
        try:
            length = int(row[1])
        except ValueError as error:
            raise ValueError(
                f"invalid HAL sequence length for {genome}:{row[0]}") from error
        if not row[0] or length <= 0 or row[0] in lengths:
            raise ValueError(f"invalid duplicate HAL sequence for {genome}:{row[0]}")
        lengths[row[0]] = length
    if not lengths:
        raise ValueError(f"halStats returned no sequences for genome {genome}")
    return lengths


def main_impl(args) -> int:
    hal_path = require_file(args.hal, "HAL alignment")
    source_arms_path = require_file(args.source_arms, "source arms BED")
    output = os.path.abspath(args.output)
    qc = os.path.abspath(args.qc)
    contract = os.path.abspath(args.contract)
    require_outputs_absent((output, qc, contract))
    if not args.source_genome.strip() or not args.target_genome.strip():
        raise ValueError("HAL source and target genome names cannot be empty")

    genomes = hal_genomes(args.hal_stats, hal_path)
    missing = [name for name in (args.source_genome, args.target_genome)
               if name not in genomes]
    if missing:
        raise ValueError(
            "genome name(s) absent from HAL: " + ", ".join(missing))

    source_rows, source_names = read_source_arms(
        source_arms_path, args.include_sex_chromosomes)
    source_lengths = hal_sequence_lengths(
        args.hal_stats, hal_path, args.source_genome)
    bed_ends: dict[str, int] = collections.defaultdict(int)
    for chrom, _start, end, _name in source_rows:
        if chrom not in source_lengths:
            raise ValueError(
                f"source arms contig is absent from HAL genome "
                f"{args.source_genome}: {chrom}")
        if end > source_lengths[chrom]:
            raise ValueError(
                f"source arm exceeds HAL sequence length for "
                f"{args.source_genome}:{chrom}: {end} > {source_lengths[chrom]}")
        bed_ends[chrom] = max(bed_ends[chrom], end)
    mismatched_lengths = [
        f"{chrom}:{bed_ends[chrom]}!={source_lengths[chrom]}"
        for chrom in sorted(bed_ends)
        if bed_ends[chrom] != source_lengths[chrom]
    ]
    if mismatched_lengths:
        raise ValueError(
            "source arms BED does not terminate at the HAL sequence length; "
            "the source coordinate build may be incompatible: " +
            ", ".join(mismatched_lengths[:5]))
    with tempfile.TemporaryDirectory(prefix="tetra_arm_hal_liftover.") as directory:
        root = Path(directory)
        filtered_source = root / "source_arms.bed"
        raw_lifted = root / "lifted_arms.bed"
        with open(filtered_source, "x", encoding="utf-8", newline="") as handle:
            for chrom, start, end, name in source_rows:
                handle.write(f"{chrom}\t{start}\t{end}\t{name}\n")
        command = [
            args.hal_liftover, hal_path, args.source_genome,
            str(filtered_source), args.target_genome, str(raw_lifted),
        ]
        completed = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if completed.returncode != 0:
            detail = completed.stderr.strip() or completed.stdout.strip()
            raise ValueError(
                f"halLiftover failed with exit {completed.returncode}: {detail}")
        if not raw_lifted.is_file() or raw_lifted.stat().st_size == 0:
            raise ValueError("halLiftover produced no mapped arm intervals")
        lifted, raw_count = read_lifted(str(raw_lifted), source_names)

    segments, conflicting_bases = disjoint_segments(lifted)
    if not segments:
        raise ValueError("HAL projection contains no unambiguous arm intervals")
    with atomic_text(output) as handle:
        for contig, start, end, arm in segments:
            handle.write(f"{contig}\t{start}\t{end}\t{arm}\n")

    mapped_names = sorted({row[3] for row in segments})
    qc_row = {
        "source_genome": args.source_genome,
        "target_genome": args.target_genome,
        "source_contigs": len(bed_ends),
        "source_intervals": len(source_rows),
        "source_logical_arms": len(source_names),
        "raw_lifted_intervals": raw_count,
        "mapped_logical_arms": len(mapped_names),
        "target_contigs": len({row[0] for row in segments}),
        "bed_intervals": len(segments),
        "mapped_bases": sum(row[2] - row[1] for row in segments),
        "conflicting_bases": conflicting_bases,
        "status": "PASS",
        "schema_version": "tetra_arm_hal_liftover_qc_v1",
    }
    write_tsv_atomic(qc, (qc_row,), QC_FIELDS)
    write_json_atomic(contract, {
        "schema_version": "tetra_arm_hal_liftover_contract_v1",
        "release": RELEASE,
        "hal": hal_path,
        "source_genome": args.source_genome,
        "target_genome": args.target_genome,
        "source_arms": source_arms_path,
        "include_sex_chromosomes": bool(args.include_sex_chromosomes),
        "output": output,
        "mapped_logical_arms": mapped_names,
        "metrics": qc_row,
        "status": "PASS",
    })
    return 0


def self_test() -> None:
    source_names = {"chr1p", "chr1q"}
    with tempfile.TemporaryDirectory(prefix="tetra_arm_hal_test.") as directory:
        lifted = Path(directory) / "lifted.bed"
        lifted.write_text(
            "ancA\t0\t100\tchr1p\n"
            "ancA\t80\t120\tchr1q\n"
            "ancB\t10\t50\tchr1p\n",
            encoding="utf-8")
        intervals, count = read_lifted(str(lifted), source_names)
        rows, conflicts = disjoint_segments(intervals)
        expected = [
            ("ancA", 0, 80, "chr1p"),
            ("ancA", 100, 120, "chr1q"),
            ("ancB", 10, 50, "chr1p"),
        ]
        if count != 3 or conflicts != 20 or rows != expected:
            raise AssertionError(
                f"unexpected HAL normalization result: {rows}, conflicts={conflicts}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {RELEASE}")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--hal")
    parser.add_argument("--source-genome")
    parser.add_argument("--target-genome")
    parser.add_argument("--source-arms")
    parser.add_argument("--output")
    parser.add_argument("--qc")
    parser.add_argument("--contract")
    parser.add_argument("--hal-stats", default="halStats")
    parser.add_argument("--hal-liftover", default="halLiftover")
    parser.add_argument("--include-sex-chromosomes", action="store_true")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.self_test:
            self_test()
            print("tetra_arm_hal_liftover self-test PASS")
            return 0
        required = (
            args.hal, args.source_genome, args.target_genome, args.source_arms,
            args.output, args.qc, args.contract,
        )
        if not all(required):
            raise ValueError(
                "--hal, --source-genome, --target-genome, --source-arms, "
                "--output, --qc, and --contract are required")
        return main_impl(args)
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
