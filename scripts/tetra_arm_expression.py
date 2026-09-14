#!/usr/bin/env python3
"""Aggregate filtered RNA MEX counts into chromosome-arm dosage evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import sys
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse

from tetra_arm_common import (
    EXPRESSION_SCHEMA,
    RELEASE,
    canonical_barcode,
    clean,
    file_record,
    natural_key,
    open_text,
    read_tsv,
    require_file,
    require_outputs_absent,
    write_json_atomic,
    write_tsv_atomic,
)


FIELDS = [
    "library", "barcode", "arm", "chromosome", "arm_counts",
    "other_autosomal_counts", "reference_autosomal_counts",
    "total_autosomal_counts", "arm_fraction", "log2_arm_to_other",
    "log2_arm_to_reference", "mapped_genes_on_arm", "nonzero_genes_on_arm",
    "matrix_value_type", "expression_input_state", "schema_version",
]


REFERENCE_DEFINITION = (
    "sum of mapped expression counts on autosomes other than the target "
    "chromosome; both target-chromosome arms are excluded"
)


def load_lines(path: str) -> list[str]:
    result = []
    with open_text(path) as handle:
        for line in handle:
            value = line.rstrip("\r\n")
            if value:
                result.append(value)
    return result


def load_gene_arms(path: str) -> dict[str, str]:
    mapping = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if not fields or not fields[0]:
                continue
            if len(fields) < 2 or not fields[1]:
                raise ValueError(f"malformed gene-arm row {path}:{line_number}")
            gene, arm = fields[0], fields[1]
            previous = mapping.get(gene)
            if previous is not None and previous != arm:
                raise ValueError(f"gene maps to multiple arms: {gene}")
            mapping[gene] = arm
    if not mapping:
        raise ValueError(f"gene-arm map contains no entries: {path}")
    return mapping


def chromosome_from_arm(arm: str) -> str:
    if arm.endswith("p") or arm.endswith("q"):
        return arm[:-1]
    return arm


def is_autosomal(arm: str) -> bool:
    chromosome = chromosome_from_arm(arm)
    normalized = chromosome.lower().removeprefix("chr")
    return normalized.isdigit() and 1 <= int(normalized) <= 22


def main_impl(args) -> int:
    library = int(args.library)
    barcodes_path = require_file(args.barcodes, "MEX barcodes")
    features_path = require_file(args.features, "MEX features")
    matrix_path = require_file(args.matrix, "MEX matrix")
    manifest_path = require_file(args.cell_manifest, "cell manifest")
    gene_arm_path = require_file(args.gene_arms, "gene-arm map")
    output_dir = Path(os.path.abspath(args.output_dir))
    evidence_path = output_dir / f"lib{library}.arm_expression.tsv.gz"
    qc_path = output_dir / f"lib{library}.expression_qc.tsv"
    contract_path = output_dir / f"lib{library}.expression_contract.json"
    require_outputs_absent((evidence_path, qc_path, contract_path))

    manifest = {}
    manifest_rows_total = 0
    for row in read_tsv(manifest_path):
        manifest_rows_total += 1
        barcode = canonical_barcode(row.get("barcode", ""))
        donor_a = clean(row.get("donor_a"))
        donor_b = clean(row.get("donor_b"))
        if not args.include_nonheterotypic and (
                not donor_a or not donor_b or donor_a == donor_b):
            continue
        if not barcode or barcode in manifest:
            raise ValueError(f"empty/duplicate barcode in {manifest_path}")
        manifest[barcode] = row
    if not manifest:
        output_dir.mkdir(parents=True, exist_ok=True)
        write_tsv_atomic(str(evidence_path), (), FIELDS)
        write_tsv_atomic(
            str(qc_path),
            [{
                "library": library, "manifest_cells": manifest_rows_total,
                "selected_manifest_cells": 0, "matrix_cells": "NA",
                "overlap_cells": 0, "missing_manifest_cells": 0,
                "matrix_genes": "NA", "mapped_genes": "NA", "arms": 0,
                "rows": 0, "status": "PASS_NO_HETEROTYPIC_TARGETS",
                "schema_version": "tetra_arm_expression_qc_v1",
            }],
            ["library", "manifest_cells", "selected_manifest_cells", "matrix_cells",
             "overlap_cells", "missing_manifest_cells", "matrix_genes",
             "mapped_genes", "arms", "rows", "status", "schema_version"])
        write_json_atomic(
            str(contract_path), {
                "schema_version": "tetra_arm_expression_contract_v1",
                "release": RELEASE, "library": library,
                "inputs": {
                    "barcodes": file_record(barcodes_path),
                    "features": file_record(features_path),
                    "matrix": file_record(matrix_path),
                    "cell_manifest": file_record(manifest_path),
                    "gene_arms": file_record(gene_arm_path),
                },
                "cells": 0, "arms": [], "rows": 0,
                "evidence_fields": FIELDS,
                "reference_autosomal_counts_definition": REFERENCE_DEFINITION,
                "terminal_state": "PASS_NO_HETEROTYPIC_TARGETS",
                "status": "PASS",
            })
        return 0

    barcode_lines = load_lines(barcodes_path)
    canonical_to_column = {}
    for column, raw in enumerate(barcode_lines):
        barcode = canonical_barcode(raw.split("\t", 1)[0])
        if barcode in canonical_to_column:
            raise ValueError(f"canonical barcode collision in {barcodes_path}: {barcode}")
        canonical_to_column[barcode] = column
    selected = sorted(set(manifest) & set(canonical_to_column), key=natural_key)
    missing = sorted(set(manifest) - set(canonical_to_column), key=natural_key)
    if missing and not args.allow_missing_barcodes:
        raise ValueError(
            f"{len(missing)} manifest barcodes are absent from filtered MEX; "
            f"first={missing[:5]}")
    if not selected:
        raise ValueError("no cell-manifest barcodes overlap the expression matrix")

    feature_lines = load_lines(features_path)
    gene_arms = load_gene_arms(gene_arm_path)
    row_arm = []
    arms = set()
    mapped_gene_counts = {}
    for line in feature_lines:
        fields = line.split("\t")
        gene_id = fields[0] if fields else ""
        gene_name = fields[1] if len(fields) > 1 else gene_id
        arm = gene_arms.get(gene_name) or gene_arms.get(gene_id) or ""
        if arm and (args.include_sex_chromosomes or is_autosomal(arm)):
            row_arm.append(arm)
            arms.add(arm)
            mapped_gene_counts[arm] = mapped_gene_counts.get(arm, 0) + 1
        else:
            row_arm.append("")
    arms = sorted(arms, key=natural_key)
    if not arms:
        raise ValueError("no expression features map to selected chromosome arms")
    arm_to_index = {arm: index for index, arm in enumerate(arms)}
    row_arm_index = np.array(
        [arm_to_index.get(arm, -1) for arm in row_arm], dtype=np.int32)

    with (gzip.open(matrix_path, "rb") if matrix_path.endswith(".gz")
          else open(matrix_path, "rb")) as handle:
        matrix = scipy.io.mmread(handle)
    if matrix.shape != (len(feature_lines), len(barcode_lines)):
        raise ValueError(
            f"MEX dimension mismatch: matrix={matrix.shape}, "
            f"features={len(feature_lines)}, barcodes={len(barcode_lines)}")
    coo = scipy.sparse.coo_matrix(matrix)
    del matrix
    coo.sum_duplicates()
    coo.eliminate_zeros()
    if (not np.all(np.isfinite(coo.data)) or
            np.any(np.asarray(coo.data, dtype=np.float64) < 0.0)):
        raise ValueError("expression matrix contains non-finite or negative values")
    matrix_is_integer = bool(np.allclose(coo.data, np.rint(coo.data)))
    if args.ambient_corrected:
        matrix_value_type = (
            "ambient_corrected_counts"
            if matrix_is_integer else "fractional_corrected_counts")
    else:
        matrix_value_type = (
            "integer_counts" if matrix_is_integer else "noninteger_values")

    old_to_selected = np.full(len(barcode_lines), -1, dtype=np.int32)
    for selected_index, barcode in enumerate(selected):
        old_to_selected[canonical_to_column[barcode]] = selected_index
    selected_columns = old_to_selected[coo.col]
    selected_arms = row_arm_index[coo.row]
    keep = (selected_columns >= 0) & (selected_arms >= 0)
    if not np.any(keep):
        raise ValueError("selected cells have no counts in mapped arm genes")
    combined = (
        selected_columns[keep].astype(np.int64) * len(arms)
        + selected_arms[keep].astype(np.int64)
    )
    sums = np.bincount(
        combined, weights=np.asarray(coo.data[keep], dtype=np.float64),
        minlength=len(selected) * len(arms),
    ).reshape(len(selected), len(arms))
    nonzero_gene_counts = np.bincount(
        combined, minlength=len(selected) * len(arms),
    ).reshape(len(selected), len(arms))
    autosomal_mask = np.array([is_autosomal(arm) for arm in arms], dtype=bool)
    autosomal_totals = sums[:, autosomal_mask].sum(axis=1)
    autosomal_chromosome_counts = {}
    for chromosome in sorted(
            {chromosome_from_arm(arm) for arm in arms if is_autosomal(arm)},
            key=natural_key):
        chromosome_mask = np.array([
            is_autosomal(arm) and chromosome_from_arm(arm) == chromosome
            for arm in arms
        ], dtype=bool)
        autosomal_chromosome_counts[chromosome] = sums[:, chromosome_mask].sum(axis=1)

    row_count = len(selected) * len(arms)

    def expression_rows():
        for cell_index, barcode in enumerate(selected):
            autosomal_total = float(autosomal_totals[cell_index])
            for arm_index, arm in enumerate(arms):
                count = float(sums[cell_index, arm_index])
                chromosome = chromosome_from_arm(arm)
                arm_is_autosomal = is_autosomal(arm)
                other = max(
                    0.0,
                    autosomal_total - count if arm_is_autosomal else autosomal_total,
                )
                same_chromosome = autosomal_chromosome_counts.get(chromosome)
                same_chromosome_count = (
                    float(same_chromosome[cell_index])
                    if same_chromosome is not None else 0.0)
                reference = max(
                    0.0, float(autosomal_totals[cell_index]) - same_chromosome_count)
                fraction = (
                    count / autosomal_total
                    if arm_is_autosomal and autosomal_total > 0 else math.nan)
                log_ratio = math.log2((count + args.pseudocount) /
                                      (other + args.pseudocount))
                reference_log_ratio = math.log2(
                    (count + args.pseudocount) /
                    (reference + args.pseudocount))
                yield {
                    "library": library,
                    "barcode": barcode,
                    "arm": arm,
                    "chromosome": chromosome,
                    "arm_counts": f"{count:.17g}",
                    "other_autosomal_counts": f"{other:.17g}",
                    "reference_autosomal_counts": f"{reference:.17g}",
                    "total_autosomal_counts": f"{autosomal_total:.17g}",
                    "arm_fraction": "NA" if not math.isfinite(fraction) else f"{fraction:.17g}",
                    "log2_arm_to_other": f"{log_ratio:.17g}",
                    "log2_arm_to_reference": f"{reference_log_ratio:.17g}",
                    "mapped_genes_on_arm": mapped_gene_counts[arm],
                    "nonzero_genes_on_arm": int(
                        nonzero_gene_counts[cell_index, arm_index]),
                    "matrix_value_type": matrix_value_type,
                    "expression_input_state": (
                        "UPSTREAM_AMBIENT_CORRECTED"
                        if args.ambient_corrected else "OBSERVED_FILTERED_COUNTS"),
                    "schema_version": EXPRESSION_SCHEMA,
                }

    output_dir.mkdir(parents=True, exist_ok=True)
    written = write_tsv_atomic(str(evidence_path), expression_rows(), FIELDS)
    if written != row_count:
        raise RuntimeError(f"expression row-count mismatch: {written} != {row_count}")
    write_tsv_atomic(
        str(qc_path),
        [{
            "library": library,
            "manifest_cells": manifest_rows_total,
            "selected_manifest_cells": len(manifest),
            "matrix_cells": len(barcode_lines),
            "overlap_cells": len(selected),
            "missing_manifest_cells": len(missing),
            "matrix_genes": len(feature_lines),
            "mapped_genes": int(np.sum(row_arm_index >= 0)),
            "arms": len(arms),
            "rows": row_count,
            "status": "PASS",
            "schema_version": "tetra_arm_expression_qc_v1",
        }],
        ["library", "manifest_cells", "selected_manifest_cells", "matrix_cells", "overlap_cells",
         "missing_manifest_cells", "matrix_genes", "mapped_genes", "arms",
         "rows", "status", "schema_version"])
    write_json_atomic(str(contract_path), {
        "schema_version": "tetra_arm_expression_contract_v1",
        "release": RELEASE,
        "library": library,
        "inputs": {
            "barcodes": file_record(barcodes_path),
            "features": file_record(features_path),
            "matrix": file_record(matrix_path),
            "cell_manifest": file_record(manifest_path),
            "gene_arms": file_record(gene_arm_path),
        },
        "cells": len(selected),
        "arms": arms,
        "rows": row_count,
        "evidence_fields": FIELDS,
        "reference_autosomal_counts_definition": REFERENCE_DEFINITION,
        "status": "PASS",
        "expression_input_state": (
            "UPSTREAM_AMBIENT_CORRECTED"
            if args.ambient_corrected else "OBSERVED_FILTERED_COUNTS"),
    })
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Aggregate a filtered 10x MEX into per-cell chromosome-arm expression evidence.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {RELEASE}")
    parser.add_argument("--library", type=int, required=True)
    parser.add_argument("--barcodes", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--cell-manifest", required=True)
    parser.add_argument("--gene-arms", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--pseudocount", type=float, default=0.5)
    parser.add_argument("--include-sex-chromosomes", action="store_true")
    parser.add_argument("--allow-missing-barcodes", action="store_true")
    parser.add_argument(
        "--include-nonheterotypic", action="store_true",
        help="also emit expression evidence for cells without two distinct donors")
    parser.add_argument(
        "--ambient-corrected", action="store_true",
        help="declare that --matrix is the upstream GEX-ambient corrected MEX")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        return main_impl(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
