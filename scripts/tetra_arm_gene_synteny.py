#!/usr/bin/env python3
"""Build ancestral-coordinate arm intervals from annotation gene names.

The output is a BED4 file suitable for tetra_arm_ase.  Multiple disjoint
intervals may share the same biological arm label; tetra_arm_ase aggregates
all such intervals into one arm.  The default anchored-contigs mode uses the
ordered mapped genes as synteny anchors, partitions each anchored contig at
midpoints between arm transitions, and covers the complete contig.  The older
gene-spans mode limits output to annotated gene spans.
"""

from __future__ import annotations

import argparse
import collections
import os
import re
import sys
import tempfile
from pathlib import Path

from tetra_arm_common import (
    RELEASE,
    atomic_text,
    open_text,
    require_file,
    require_outputs_absent,
    write_json_atomic,
    write_tsv_atomic,
)


QC_FIELDS = (
    "projection_mode", "reference_contigs", "reference_bases",
    "annotation_rows", "feature_rows", "named_features", "mapped_features",
    "unmapped_features", "mapped_gene_names", "contigs", "logical_arms",
    "bed_intervals", "mapped_bases", "anchored_contigs",
    "single_arm_contigs", "arm_transitions", "unanchored_contigs",
    "unanchored_bases", "ambiguous_anchor_positions", "conflicting_bases",
    "status", "schema_version",
)


def chromosome_from_arm(arm: str) -> str:
    return arm[:-1] if arm.endswith(("p", "q")) else arm


def is_sex_arm(arm: str) -> bool:
    value = chromosome_from_arm(arm).lower().removeprefix("chr")
    return value in {"x", "y"}


def load_gene_arms(path: str, include_sex: bool) -> dict[str, str]:
    result: dict[str, str] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 2 or not fields[0].strip() or not fields[1].strip():
                raise ValueError(f"malformed gene-arm row {path}:{line_number}")
            gene, arm = fields[0].strip(), fields[1].strip()
            if not include_sex and is_sex_arm(arm):
                continue
            previous = result.get(gene)
            if previous is not None and previous != arm:
                raise ValueError(f"gene maps to multiple arms in {path}: {gene}")
            result[gene] = arm
    if not result:
        raise ValueError(f"gene-arm map contains no selected entries: {path}")
    return result


def parse_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for raw in text.split(";"):
        field = raw.strip()
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        else:
            parts = field.split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts
        value = value.strip().strip('"').strip("'")
        if key.strip() and value:
            result.setdefault(key.strip(), value)
    return result


def candidate_names(attributes: dict[str, str]) -> list[str]:
    names: list[str] = []
    for key in ("gene_name", "Name", "gene", "gene_id", "ID"):
        value = attributes.get(key, "").strip()
        if not value:
            continue
        for prefix in ("gene:", "gene-", "rna:"):
            if value.startswith(prefix):
                value = value[len(prefix):]
        for candidate in (value, re.sub(r"\.\d+$", "", value)):
            if candidate and candidate not in names:
                names.append(candidate)
    return names


def build_intervals(annotation: str, gene_arms: dict[str, str],
                    feature_types: set[str]):
    intervals: dict[str, list[tuple[int, int, str]]] = collections.defaultdict(list)
    metrics = collections.Counter()
    mapped_names = set()
    with open_text(annotation) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            metrics["annotation_rows"] += 1
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 9:
                raise ValueError(
                    f"malformed GTF/GFF row {annotation}:{line_number}; expected 9 columns")
            if fields[2].lower() not in feature_types:
                continue
            metrics["feature_rows"] += 1
            attributes = parse_attributes(fields[8])
            names = candidate_names(attributes)
            if not names:
                continue
            metrics["named_features"] += 1
            matches = {gene_arms[name] for name in names if name in gene_arms}
            if not matches:
                metrics["unmapped_features"] += 1
                continue
            if len(matches) != 1:
                metrics["conflicting_name_features"] += 1
                continue
            try:
                start = int(fields[3]) - 1
                end = int(fields[4])
            except ValueError as error:
                raise ValueError(
                    f"noninteger coordinates at {annotation}:{line_number}") from error
            if start < 0 or end <= start or not fields[0].strip():
                raise ValueError(f"invalid feature interval at {annotation}:{line_number}")
            arm = next(iter(matches))
            intervals[fields[0].strip()].append((start, end, arm))
            mapped_names.update(name for name in names if gene_arms.get(name) == arm)
            metrics["mapped_features"] += 1
    metrics["mapped_gene_names"] = len(mapped_names)
    return intervals, metrics


def disjoint_segments(intervals: dict[str, list[tuple[int, int, str]]]):
    output: list[tuple[str, int, int, str]] = []
    conflicting_bases = 0
    for contig in sorted(intervals):
        events: dict[int, collections.Counter[str]] = collections.defaultdict(collections.Counter)
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


def load_fai(path: str) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 2 or not fields[0].strip():
                raise ValueError(f"malformed FASTA index row {path}:{line_number}")
            try:
                length = int(fields[1])
            except ValueError as error:
                raise ValueError(
                    f"noninteger contig length at {path}:{line_number}") from error
            contig = fields[0].strip()
            if length <= 0:
                raise ValueError(f"invalid contig length at {path}:{line_number}")
            if contig in lengths:
                raise ValueError(f"duplicate contig in FASTA index: {contig}")
            lengths[contig] = length
    if not lengths:
        raise ValueError(f"FASTA index contains no contigs: {path}")
    return lengths


def anchored_contig_segments(
        intervals: dict[str, list[tuple[int, int, str]]],
        contig_lengths: dict[str, int]):
    """Partition complete anchored contigs by nearest ordered arm anchor."""
    output: list[tuple[str, int, int, str]] = []
    metrics = collections.Counter()
    for contig in sorted(intervals):
        if contig not in contig_lengths:
            raise ValueError(
                f"mapped annotation contig is absent from FASTA index: {contig}")
        contig_length = contig_lengths[contig]
        by_center: dict[int, set[str]] = collections.defaultdict(set)
        for start, end, arm in intervals[contig]:
            if end > contig_length:
                raise ValueError(
                    f"annotation interval exceeds indexed contig length: {contig}:{end}")
            by_center[(start + end) // 2].add(arm)

        anchors: list[tuple[int, str]] = []
        for center in sorted(by_center):
            labels = by_center[center]
            if len(labels) != 1:
                metrics["ambiguous_anchor_positions"] += 1
                continue
            anchors.append((center, next(iter(labels))))
        if not anchors:
            metrics["unanchored_contigs"] += 1
            metrics["unanchored_bases"] += contig_length
            continue

        # Consecutive anchors with the same arm form one run.  A boundary
        # between unlike runs is halfway between their nearest anchor centers.
        runs: list[list[object]] = []
        for center, arm in anchors:
            if runs and runs[-1][2] == arm:
                runs[-1][1] = center
            else:
                runs.append([center, center, arm])
        boundaries = [
            (int(runs[index][1]) + int(runs[index + 1][0])) // 2
            for index in range(len(runs) - 1)
        ]
        start = 0
        for index, run in enumerate(runs):
            end = boundaries[index] if index < len(boundaries) else contig_length
            if end <= start:
                raise ValueError(
                    f"nonincreasing synteny partition boundary on contig {contig}")
            output.append((contig, start, end, str(run[2])))
            start = end
        metrics["anchored_contigs"] += 1
        metrics["arm_transitions"] += max(0, len(runs) - 1)
        if len(runs) == 1:
            metrics["single_arm_contigs"] += 1

    for contig, length in contig_lengths.items():
        if contig not in intervals:
            metrics["unanchored_contigs"] += 1
            metrics["unanchored_bases"] += length
    return output, metrics


def main_impl(args) -> int:
    annotation = require_file(args.annotation, "ancestral gene annotation")
    gene_arm_path = require_file(args.gene_arms, "gene-arm map")
    output = os.path.abspath(args.output)
    qc = os.path.abspath(args.qc)
    contract = os.path.abspath(args.contract)
    require_outputs_absent((output, qc, contract))
    feature_types = {value.strip().lower() for value in args.feature_type if value.strip()}
    if not feature_types:
        raise ValueError("at least one --feature-type is required")

    projection_mode = args.projection_mode
    fai_path = ""
    contig_lengths: dict[str, int] = {}
    if projection_mode == "anchored-contigs":
        if not args.reference_fai:
            raise ValueError(
                "--reference-fai is required for --projection-mode anchored-contigs")
        fai_path = require_file(args.reference_fai, "ancestral FASTA index")
        contig_lengths = load_fai(fai_path)

    gene_arms = load_gene_arms(gene_arm_path, args.include_sex_chromosomes)
    intervals, metrics = build_intervals(annotation, gene_arms, feature_types)
    if projection_mode == "anchored-contigs":
        segments, partition_metrics = anchored_contig_segments(
            intervals, contig_lengths)
        metrics.update(partition_metrics)
        conflicting_bases = 0
    else:
        segments, conflicting_bases = disjoint_segments(intervals)
    if not segments:
        raise ValueError(
            "no unambiguous annotation intervals mapped to the gene-arm table")

    with atomic_text(output) as handle:
        for contig, start, end, arm in segments:
            handle.write(f"{contig}\t{start}\t{end}\t{arm}\n")

    logical_arms = sorted({row[3] for row in segments})
    qc_row = {
        "projection_mode": projection_mode,
        "reference_contigs": len(contig_lengths),
        "reference_bases": sum(contig_lengths.values()),
        "annotation_rows": metrics["annotation_rows"],
        "feature_rows": metrics["feature_rows"],
        "named_features": metrics["named_features"],
        "mapped_features": metrics["mapped_features"],
        "unmapped_features": metrics["unmapped_features"],
        "mapped_gene_names": metrics["mapped_gene_names"],
        "contigs": len({row[0] for row in segments}),
        "logical_arms": len(logical_arms),
        "bed_intervals": len(segments),
        "mapped_bases": sum(row[2] - row[1] for row in segments),
        "anchored_contigs": metrics["anchored_contigs"],
        "single_arm_contigs": metrics["single_arm_contigs"],
        "arm_transitions": metrics["arm_transitions"],
        "unanchored_contigs": metrics["unanchored_contigs"],
        "unanchored_bases": metrics["unanchored_bases"],
        "ambiguous_anchor_positions": metrics["ambiguous_anchor_positions"],
        "conflicting_bases": conflicting_bases,
        "status": "PASS",
        "schema_version": "tetra_arm_gene_synteny_qc_v1",
    }
    write_tsv_atomic(qc, (qc_row,), QC_FIELDS)
    write_json_atomic(contract, {
        "schema_version": "tetra_arm_gene_synteny_contract_v1",
        "release": RELEASE,
        "annotation": annotation,
        "gene_arms": gene_arm_path,
        "reference_fai": fai_path,
        "projection_mode": projection_mode,
        "feature_types": sorted(feature_types),
        "include_sex_chromosomes": bool(args.include_sex_chromosomes),
        "output": output,
        "logical_arms": logical_arms,
        "metrics": qc_row,
        "status": "PASS",
    })
    return 0


def self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="tetra_arm_gene_synteny.") as directory:
        root = Path(directory)
        annotation = root / "test.gtf"
        mapping = root / "arms.tsv"
        reference_fai = root / "test.fa.fai"
        annotation.write_text(
            'ancA\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "G1"; gene_name "G1";\n'
            'ancA\ttest\tgene\t80\t150\t.\t+\t.\tgene_id "G2"; gene_name "G2";\n'
            'ancB\ttest\tgene\t5\t25\t.\t+\t.\tgene_id "G3"; gene_name "G3";\n',
            encoding="utf-8")
        mapping.write_text("G1\tchr1p\nG2\tchr1q\nG3\tchr1p\n", encoding="utf-8")
        reference_fai.write_text("ancA\t200\t0\t0\t0\nancB\t80\t0\t0\t0\n",
                                 encoding="utf-8")
        args = argparse.Namespace(
            annotation=str(annotation), gene_arms=str(mapping),
            reference_fai=str(reference_fai), projection_mode="anchored-contigs",
            output=str(root / "out.bed"), qc=str(root / "out.qc.tsv"),
            contract=str(root / "out.contract.json"), feature_type=["gene"],
            include_sex_chromosomes=False)
        main_impl(args)
        rows = (root / "out.bed").read_text(encoding="utf-8").splitlines()
        expected = [
            "ancA\t0\t82\tchr1p", "ancA\t82\t200\tchr1q",
            "ancB\t0\t80\tchr1p",
        ]
        if rows != expected:
            raise AssertionError(f"unexpected self-test BED: {rows}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {RELEASE}")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--annotation")
    parser.add_argument("--gene-arms")
    parser.add_argument("--reference-fai")
    parser.add_argument(
        "--projection-mode", choices=("anchored-contigs", "gene-spans"),
        default="anchored-contigs")
    parser.add_argument("--output")
    parser.add_argument("--qc")
    parser.add_argument("--contract")
    parser.add_argument("--feature-type", action="append", default=["gene"])
    parser.add_argument("--include-sex-chromosomes", action="store_true")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.self_test:
            self_test()
            print("tetra_arm_gene_synteny self-test PASS")
            return 0
        required = (args.annotation, args.gene_arms, args.output, args.qc, args.contract)
        if not all(required):
            raise ValueError(
                "--annotation, --gene-arms, --output, --qc, and --contract are required")
        return main_impl(args)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
