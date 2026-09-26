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
    EXPRESSION_MODEL_SCHEMA,
    EXPRESSION_SCHEMA,
    HYBRID_CELL_MANIFEST_SCHEMA,
    RELEASE,
    canonical_barcode,
    clean,
    file_record,
    natural_key,
    open_text,
    read_tsv,
    require_file,
    require_outputs_absent,
    truthy,
    validate_tsv_schema,
    write_json_atomic,
    write_tsv_atomic,
)
from tetra_arm_hybrid_common import stable_gene_fold


FIELDS = [
    "library", "barcode", "arm", "chromosome", "arm_counts",
    "other_autosomal_counts", "reference_autosomal_counts",
    "total_autosomal_counts", "arm_fraction", "log2_arm_to_other",
    "log2_arm_to_reference", "mapped_genes_on_arm", "nonzero_genes_on_arm",
    "matrix_value_type", "expression_input_state", "schema_version",
]

MODEL_FIELDS = [
    "library", "barcode", "chromosome", "p_arm", "q_arm",
    "donor_a", "donor_b", "donor_pair", "uid", "calibration_group",
    "hybrid_target_eligible", "expression_reference_eligible",
    "cell_group_source", "cell_group_target_chromosome_excluded",
    "p_raw_count_total", "q_raw_count_total", "reference_raw_count_total",
    "mapped_autosomal_library_size",
    "p_mapped_genes", "q_mapped_genes", "reference_mapped_genes",
    "p_nonzero_genes", "q_nonzero_genes", "reference_nonzero_genes",
    "p_mean_log1p_cpm", "q_mean_log1p_cpm", "reference_mean_log1p_cpm",
    "p_score", "q_score",
    "p_fold0_mean_log1p_cpm", "p_fold1_mean_log1p_cpm",
    "q_fold0_mean_log1p_cpm", "q_fold1_mean_log1p_cpm",
    "reference_fold0_mean_log1p_cpm", "reference_fold1_mean_log1p_cpm",
    "p_fold0_mapped_genes", "p_fold1_mapped_genes",
    "q_fold0_mapped_genes", "q_fold1_mapped_genes",
    "reference_fold0_mapped_genes", "reference_fold1_mapped_genes",
    "p_fold0_score", "p_fold1_score", "q_fold0_score", "q_fold1_score",
    "p_top1_fraction", "p_top5_fraction", "p_top10_fraction",
    "q_top1_fraction", "q_top5_fraction", "q_top10_fraction",
    "reference_top1_fraction", "reference_top5_fraction",
    "reference_top10_fraction", "expressed_autosomal_genes",
    "ambient_c", "matrix_value_type", "expression_input_state",
    "normalization_method", "gene_fold_method", "expression_scale_factor",
    "score_formula", "reference_definition", "schema_version",
]


REFERENCE_DEFINITION = (
    "sum of mapped expression counts on autosomes other than the target "
    "chromosome; both target-chromosome arms are excluded"
)

MODEL_REFERENCE_DEFINITION = (
    "complete mapped autosomal gene universe excluding every gene on the "
    "target chromosome; zeros remain in each region denominator"
)
MODEL_SCORE_FORMULA = (
    "ln(expm1(mean_log1p_cpm_region)+1e-9)-"
    "ln(expm1(mean_log1p_cpm_chromosome_excluded_reference)+1e-9)"
)
NORMALIZATION_METHOD = "NATURAL_LOG1P_CPM_COMPLETE_MAPPED_GENE_DENOMINATOR"
GENE_FOLD_METHOD = "BLAKE2B_CANONICAL_GENE_IDENTIFIER_MOD_2"


def format_number(value: float) -> str:
    return "NA" if not math.isfinite(value) else f"{value:.17g}"


def score_from_means(region_mean: float, reference_mean: float,
                     epsilon: float = 1e-9) -> float:
    if not (math.isfinite(region_mean) and math.isfinite(reference_mean)):
        return math.nan
    return (math.log(max(math.expm1(region_mean), 0.0) + epsilon)
            - math.log(max(math.expm1(reference_mean), 0.0) + epsilon))


def top_k_sums(group_ids: np.ndarray, values: np.ndarray,
               group_count: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return exact top-1/5/10 sums per integer group without dense genes."""
    outputs = [np.zeros(group_count, dtype=np.float64) for _ in range(3)]
    if values.size == 0:
        return outputs[0], outputs[1], outputs[2]
    order = np.lexsort((-values, group_ids))
    ordered_groups = group_ids[order]
    ordered_values = values[order]
    starts = np.flatnonzero(np.r_[True, ordered_groups[1:] != ordered_groups[:-1]])
    ends = np.r_[starts[1:], len(order)]
    for start, end in zip(starts, ends):
        group = int(ordered_groups[start])
        segment = ordered_values[start:end]
        outputs[0][group] = float(np.sum(segment[:1]))
        outputs[1][group] = float(np.sum(segment[:5]))
        outputs[2][group] = float(np.sum(segment[:10]))
    return outputs[0], outputs[1], outputs[2]


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
    if args.emit_model_evidence and not args.hybrid_cell_manifest:
        raise ValueError("--emit-model-evidence requires --hybrid-cell-manifest")
    if args.gene_folds != 2:
        raise ValueError("hybrid v1 requires exactly --gene-folds 2")
    if not math.isfinite(args.expression_scale_factor) or args.expression_scale_factor <= 0:
        raise ValueError("--expression-scale-factor must be finite and positive")

    barcodes_path = require_file(args.barcodes, "MEX barcodes")
    features_path = require_file(args.features, "MEX features")
    matrix_path = require_file(args.matrix, "MEX matrix")
    legacy_manifest_path = require_file(args.cell_manifest, "cell manifest")
    gene_arm_path = require_file(args.gene_arms, "gene-arm map")
    hybrid_manifest_path = ""
    if args.emit_model_evidence:
        hybrid_manifest_path = require_file(
            args.hybrid_cell_manifest, "hybrid cell manifest")
        validate_tsv_schema(
            hybrid_manifest_path, HYBRID_CELL_MANIFEST_SCHEMA,
            {"library", "barcode", "hybrid_target_eligible",
             "expression_reference_eligible", "schema_version"})

    output_dir = Path(os.path.abspath(args.output_dir))
    evidence_path = output_dir / f"lib{library}.arm_expression.tsv.gz"
    model_path = output_dir / f"lib{library}.arm_expression_model.tsv.gz"
    qc_path = output_dir / f"lib{library}.expression_qc.tsv"
    contract_path = output_dir / f"lib{library}.expression_contract.json"
    output_paths = [evidence_path, qc_path, contract_path]
    if args.emit_model_evidence:
        output_paths.append(model_path)
    require_outputs_absent(output_paths)

    manifest: dict[str, dict[str, str]] = {}
    manifest_rows_total = 0
    selected_manifest_path = hybrid_manifest_path or legacy_manifest_path
    for row in read_tsv(selected_manifest_path):
        manifest_rows_total += 1
        raw_library = clean(row.get("library")).lower().removeprefix("lib")
        if raw_library and int(raw_library) != library:
            raise ValueError(
                f"wrong library in {selected_manifest_path}: {raw_library}")
        barcode = canonical_barcode(row.get("barcode", ""))
        if args.emit_model_evidence:
            if clean(row.get("schema_version")) != HYBRID_CELL_MANIFEST_SCHEMA:
                raise ValueError(
                    f"incompatible hybrid manifest schema: {selected_manifest_path}")
            if not (truthy(row.get("hybrid_target_eligible")) or
                    truthy(row.get("expression_reference_eligible"))):
                continue
        else:
            donor_a = clean(row.get("donor_a"))
            donor_b = clean(row.get("donor_b"))
            if not args.include_nonheterotypic and (
                    not donor_a or not donor_b or donor_a == donor_b):
                continue
        if not barcode or barcode in manifest:
            raise ValueError(f"empty/duplicate barcode in {selected_manifest_path}")
        manifest[barcode] = row

    if not manifest:
        output_dir.mkdir(parents=True, exist_ok=True)
        write_tsv_atomic(str(evidence_path), (), FIELDS)
        if args.emit_model_evidence:
            write_tsv_atomic(
                str(model_path), (), MODEL_FIELDS, deterministic_gzip=True)
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
        contract = {
            "schema_version": "tetra_arm_expression_contract_v1",
            "release": RELEASE, "library": library,
            "inputs": {
                "barcodes": file_record(barcodes_path),
                "features": file_record(features_path),
                "matrix": file_record(matrix_path),
                "cell_manifest": file_record(legacy_manifest_path),
                "gene_arms": file_record(gene_arm_path),
            },
            "cells": 0, "arms": [], "rows": 0,
            "evidence_fields": FIELDS,
            "reference_autosomal_counts_definition": REFERENCE_DEFINITION,
            "terminal_state": "PASS_NO_HETEROTYPIC_TARGETS",
            "status": "PASS",
        }
        if args.emit_model_evidence:
            contract["inputs"]["hybrid_cell_manifest"] = file_record(
                hybrid_manifest_path)
            contract["model_evidence"] = {
                "path": str(model_path), "schema_version": EXPRESSION_MODEL_SCHEMA,
                "rows": 0, "fields": MODEL_FIELDS,
                "chromosomes": [],
                "normalization_method": NORMALIZATION_METHOD,
                "gene_fold_method": GENE_FOLD_METHOD,
                "gene_folds": args.gene_folds,
                "scale_factor": args.expression_scale_factor,
                "epsilon": 1e-9,
                "score_formula": MODEL_SCORE_FORMULA,
                "reference_definition": MODEL_REFERENCE_DEFINITION,
                "raw_and_corrected_rows_mixed": False,
                "terminal_state": "PASS_NO_ELIGIBLE_HYBRID_CELLS",
            }
        write_json_atomic(str(contract_path), contract)
        return 0

    barcode_lines = load_lines(barcodes_path)
    canonical_to_column: dict[str, int] = {}
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
    row_arm: list[str] = []
    row_gene_id: list[str] = []
    mapped_gene_ids: set[str] = set()
    arms_set: set[str] = set()
    mapped_gene_counts: dict[str, int] = {}
    mapped_fold_gene_counts: dict[tuple[str, int], int] = {}
    for line_number, line in enumerate(feature_lines, start=1):
        fields = line.split("\t")
        gene_id = clean(fields[0] if fields else "")
        gene_name = clean(fields[1] if len(fields) > 1 else gene_id)
        arm = gene_arms.get(gene_name) or gene_arms.get(gene_id) or ""
        if arm and (args.include_sex_chromosomes or is_autosomal(arm)):
            if args.emit_model_evidence and not gene_id:
                raise ValueError(f"mapped feature has no canonical gene ID: {features_path}:{line_number}")
            if args.emit_model_evidence and gene_id in mapped_gene_ids:
                raise ValueError(f"duplicate mapped canonical gene ID: {gene_id}")
            mapped_gene_ids.add(gene_id)
            row_arm.append(arm)
            row_gene_id.append(gene_id)
            arms_set.add(arm)
            mapped_gene_counts[arm] = mapped_gene_counts.get(arm, 0) + 1
            fold = (stable_gene_fold(gene_id, args.gene_folds)
                    if args.emit_model_evidence else 0)
            mapped_fold_gene_counts[(arm, fold)] = (
                mapped_fold_gene_counts.get((arm, fold), 0) + 1)
        else:
            row_arm.append("")
            row_gene_id.append(gene_id)
    arms = sorted(arms_set, key=natural_key)
    if not arms:
        raise ValueError("no expression features map to selected chromosome arms")
    arm_to_index = {arm: index for index, arm in enumerate(arms)}
    row_arm_index = np.asarray(
        [arm_to_index.get(arm, -1) for arm in row_arm], dtype=np.int32)
    row_fold_index = np.asarray([
        stable_gene_fold(gene_id, args.gene_folds)
        if args.emit_model_evidence and arm else (0 if arm else -1)
        for gene_id, arm in zip(row_gene_id, row_arm)
    ], dtype=np.int8)

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
    if args.emit_model_evidence and not args.ambient_corrected and not matrix_is_integer:
        raise ValueError(
            "raw hybrid expression evidence requires integer filtered MEX counts")
    if args.ambient_corrected:
        matrix_value_type = (
            "ambient_corrected_counts"
            if matrix_is_integer else "fractional_corrected_counts")
    else:
        matrix_value_type = (
            "integer_counts" if matrix_is_integer else "noninteger_values")
    input_state = (
        "UPSTREAM_AMBIENT_CORRECTED"
        if args.ambient_corrected else "OBSERVED_FILTERED_COUNTS")

    old_to_selected = np.full(len(barcode_lines), -1, dtype=np.int32)
    for selected_index, barcode in enumerate(selected):
        old_to_selected[canonical_to_column[barcode]] = selected_index
    selected_columns_all = old_to_selected[coo.col]
    selected_arms_all = row_arm_index[coo.row]
    keep = (selected_columns_all >= 0) & (selected_arms_all >= 0)
    if not np.any(keep) and not args.emit_model_evidence:
        raise ValueError("selected cells have no counts in mapped arm genes")
    selected_columns = selected_columns_all[keep].astype(np.int64)
    selected_arms = selected_arms_all[keep].astype(np.int64)
    values = np.asarray(coo.data[keep], dtype=np.float64)
    selected_folds = row_fold_index[coo.row[keep]].astype(np.int64)
    combined = selected_columns * len(arms) + selected_arms
    sums = np.bincount(
        combined, weights=values,
        minlength=len(selected) * len(arms),
    ).reshape(len(selected), len(arms))
    nonzero_gene_counts = np.bincount(
        combined, minlength=len(selected) * len(arms),
    ).reshape(len(selected), len(arms))
    top1, top5, top10 = top_k_sums(
        combined, values, len(selected) * len(arms))
    top1 = top1.reshape(len(selected), len(arms))
    top5 = top5.reshape(len(selected), len(arms))
    top10 = top10.reshape(len(selected), len(arms))

    autosomal_mask = np.asarray([is_autosomal(arm) for arm in arms], dtype=bool)
    autosomal_totals = sums[:, autosomal_mask].sum(axis=1)
    autosomal_nonzero = nonzero_gene_counts[:, autosomal_mask].sum(axis=1)
    autosomal_chromosome_counts: dict[str, np.ndarray] = {}
    autosomal_chromosome_nonzero: dict[str, np.ndarray] = {}
    chromosomes = sorted(
        {chromosome_from_arm(arm) for arm in arms if is_autosomal(arm)},
        key=natural_key)
    chromosome_masks: dict[str, np.ndarray] = {}
    chromosome_sides: dict[str, dict[str, int]] = {
        chromosome: {} for chromosome in chromosomes}
    for arm_index, arm in enumerate(arms):
        if not is_autosomal(arm):
            continue
        chromosome = chromosome_from_arm(arm)
        side = arm[-1].lower() if arm.lower().endswith(("p", "q")) else ""
        if side:
            if side in chromosome_sides[chromosome]:
                raise ValueError(
                    f"multiple expression arms for chromosome {chromosome}{side}")
            chromosome_sides[chromosome][side] = arm_index
    for chromosome in chromosomes:
        chromosome_mask = np.asarray([
            is_autosomal(arm) and chromosome_from_arm(arm) == chromosome
            for arm in arms
        ], dtype=bool)
        chromosome_masks[chromosome] = chromosome_mask
        autosomal_chromosome_counts[chromosome] = sums[:, chromosome_mask].sum(axis=1)
        autosomal_chromosome_nonzero[chromosome] = (
            nonzero_gene_counts[:, chromosome_mask].sum(axis=1))

    z_sums = np.zeros((len(selected), len(arms)), dtype=np.float64)
    fold_z_sums = np.zeros((len(selected), len(arms), args.gene_folds),
                           dtype=np.float64)
    if values.size:
        denominators = autosomal_totals[selected_columns]
        z_values = np.zeros_like(values)
        positive = denominators > 0.0
        z_values[positive] = np.log1p(
            args.expression_scale_factor * values[positive] / denominators[positive])
        z_sums = np.bincount(
            combined, weights=z_values,
            minlength=len(selected) * len(arms)).reshape(len(selected), len(arms))
        combined_fold = (
            (selected_columns * len(arms) + selected_arms) * args.gene_folds
            + selected_folds)
        fold_z_sums = np.bincount(
            combined_fold, weights=z_values,
            minlength=len(selected) * len(arms) * args.gene_folds,
        ).reshape(len(selected), len(arms), args.gene_folds)

    arm_gene_denominators = np.asarray(
        [mapped_gene_counts[arm] for arm in arms], dtype=np.int64)
    arm_fold_denominators = np.asarray([
        [mapped_fold_gene_counts.get((arm, fold), 0)
         for fold in range(args.gene_folds)]
        for arm in arms
    ], dtype=np.int64)
    total_autosomal_genes = int(np.sum(arm_gene_denominators[autosomal_mask]))
    total_autosomal_fold_genes = np.sum(
        arm_fold_denominators[autosomal_mask, :], axis=0)
    autosomal_z_totals = z_sums[:, autosomal_mask].sum(axis=1)
    autosomal_fold_z_totals = fold_z_sums[:, autosomal_mask, :].sum(axis=1)

    reference_top = {
        1: np.zeros((len(selected), len(chromosomes)), dtype=np.float64),
        5: np.zeros((len(selected), len(chromosomes)), dtype=np.float64),
        10: np.zeros((len(selected), len(chromosomes)), dtype=np.float64),
    }
    if args.emit_model_evidence and values.size:
        value_chromosomes = np.asarray([
            chromosome_from_arm(arms[index]) for index in selected_arms],
            dtype=object)
        cell_order = np.argsort(selected_columns, kind="stable")
        ordered_cells = selected_columns[cell_order]
        ordered_values = values[cell_order]
        ordered_chromosomes = value_chromosomes[cell_order]
        starts = np.flatnonzero(np.r_[True, ordered_cells[1:] != ordered_cells[:-1]])
        ends = np.r_[starts[1:], len(cell_order)]
        chromosome_to_index = {value: index for index, value in enumerate(chromosomes)}
        for start, end in zip(starts, ends):
            cell_index = int(ordered_cells[start])
            local_order = np.argsort(-ordered_values[start:end], kind="stable")
            local_values = ordered_values[start:end][local_order]
            local_chromosomes = ordered_chromosomes[start:end][local_order]
            for chromosome, chromosome_index in chromosome_to_index.items():
                retained = local_values[local_chromosomes != chromosome][:10]
                reference_top[1][cell_index, chromosome_index] = float(
                    np.sum(retained[:1]))
                reference_top[5][cell_index, chromosome_index] = float(
                    np.sum(retained[:5]))
                reference_top[10][cell_index, chromosome_index] = float(
                    np.sum(retained[:10]))

    row_count = len(selected) * len(arms)
    model_row_count = len(selected) * len(chromosomes)

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
                reference = max(0.0, autosomal_total - same_chromosome_count)
                fraction = (
                    count / autosomal_total
                    if arm_is_autosomal and autosomal_total > 0 else math.nan)
                log_ratio = math.log2((count + args.pseudocount) /
                                      (other + args.pseudocount))
                reference_log_ratio = math.log2(
                    (count + args.pseudocount) /
                    (reference + args.pseudocount))
                yield {
                    "library": library, "barcode": barcode, "arm": arm,
                    "chromosome": chromosome, "arm_counts": f"{count:.17g}",
                    "other_autosomal_counts": f"{other:.17g}",
                    "reference_autosomal_counts": f"{reference:.17g}",
                    "total_autosomal_counts": f"{autosomal_total:.17g}",
                    "arm_fraction": format_number(fraction),
                    "log2_arm_to_other": f"{log_ratio:.17g}",
                    "log2_arm_to_reference": f"{reference_log_ratio:.17g}",
                    "mapped_genes_on_arm": mapped_gene_counts[arm],
                    "nonzero_genes_on_arm": int(
                        nonzero_gene_counts[cell_index, arm_index]),
                    "matrix_value_type": matrix_value_type,
                    "expression_input_state": input_state,
                    "schema_version": EXPRESSION_SCHEMA,
                }

    def fraction_or_na(numerator: float, denominator: float) -> str:
        return format_number(numerator / denominator) if denominator > 0.0 else "NA"

    def model_rows():
        chromosome_to_index = {value: index for index, value in enumerate(chromosomes)}
        for cell_index, barcode in enumerate(selected):
            manifest_row = manifest[barcode]
            for chromosome in chromosomes:
                chromosome_index = chromosome_to_index[chromosome]
                mask = chromosome_masks[chromosome]
                sides = chromosome_sides[chromosome]
                chromosome_gene_count = int(np.sum(arm_gene_denominators[mask]))
                reference_gene_count = total_autosomal_genes - chromosome_gene_count
                chromosome_fold_genes = np.sum(arm_fold_denominators[mask, :], axis=0)
                reference_fold_genes = total_autosomal_fold_genes - chromosome_fold_genes
                chromosome_z = float(np.sum(z_sums[cell_index, mask]))
                reference_z = float(autosomal_z_totals[cell_index] - chromosome_z)
                reference_mean = (
                    reference_z / reference_gene_count
                    if reference_gene_count > 0 else math.nan)
                chromosome_fold_z = np.sum(fold_z_sums[cell_index, mask, :], axis=0)
                reference_fold_z = autosomal_fold_z_totals[cell_index, :] - chromosome_fold_z
                reference_fold_means = [
                    float(reference_fold_z[fold] / reference_fold_genes[fold])
                    if reference_fold_genes[fold] > 0 else math.nan
                    for fold in range(args.gene_folds)
                ]
                reference_count = max(
                    0.0, float(autosomal_totals[cell_index]
                               - autosomal_chromosome_counts[chromosome][cell_index]))
                reference_nonzero = max(
                    0, int(autosomal_nonzero[cell_index]
                           - autosomal_chromosome_nonzero[chromosome][cell_index]))

                def arm_metrics(side: str) -> dict[str, object]:
                    arm_index = sides.get(side)
                    if arm_index is None:
                        return {
                            "arm": "NA", "count": math.nan, "genes": 0,
                            "nonzero": 0, "mean": math.nan, "score": math.nan,
                            "fold_means": [math.nan, math.nan],
                            "fold_genes": [0, 0],
                            "fold_scores": [math.nan, math.nan],
                            "top": [math.nan, math.nan, math.nan],
                        }
                    count = float(sums[cell_index, arm_index])
                    genes = int(arm_gene_denominators[arm_index])
                    mean = (float(z_sums[cell_index, arm_index]) / genes
                            if genes > 0 else math.nan)
                    fold_genes = [int(value) for value in arm_fold_denominators[arm_index, :]]
                    fold_means = [
                        float(fold_z_sums[cell_index, arm_index, fold]
                              / fold_genes[fold])
                        if fold_genes[fold] > 0 else math.nan
                        for fold in range(args.gene_folds)
                    ]
                    return {
                        "arm": arms[arm_index], "count": count, "genes": genes,
                        "nonzero": int(nonzero_gene_counts[cell_index, arm_index]),
                        "mean": mean, "score": score_from_means(mean, reference_mean),
                        "fold_means": fold_means, "fold_genes": fold_genes,
                        "fold_scores": [
                            score_from_means(fold_means[fold], reference_fold_means[fold])
                            for fold in range(args.gene_folds)],
                        "top": [
                            (float(values_array[cell_index, arm_index]) / count
                             if count > 0 else math.nan)
                            for values_array in (top1, top5, top10)
                        ],
                    }

                p = arm_metrics("p")
                q = arm_metrics("q")
                yield {
                    "library": library, "barcode": barcode,
                    "chromosome": chromosome, "p_arm": p["arm"], "q_arm": q["arm"],
                    "donor_a": clean(manifest_row.get("donor_a")) or "NA",
                    "donor_b": clean(manifest_row.get("donor_b")) or "NA",
                    "donor_pair": clean(manifest_row.get("donor_pair")) or "NA",
                    "uid": clean(manifest_row.get("uid")) or "NA",
                    "calibration_group": clean(
                        manifest_row.get("calibration_group")) or f"lib{library}",
                    "hybrid_target_eligible": int(truthy(
                        manifest_row.get("hybrid_target_eligible"))),
                    "expression_reference_eligible": int(truthy(
                        manifest_row.get("expression_reference_eligible"))),
                    "cell_group_source": clean(
                        manifest_row.get("cell_group_source")) or "LIBRARY_FALLBACK",
                    "cell_group_target_chromosome_excluded": int(truthy(
                        manifest_row.get("cell_group_target_chromosome_excluded"))),
                    "p_raw_count_total": format_number(p["count"]),
                    "q_raw_count_total": format_number(q["count"]),
                    "reference_raw_count_total": format_number(reference_count),
                    "mapped_autosomal_library_size": format_number(
                        float(autosomal_totals[cell_index])),
                    "p_mapped_genes": p["genes"], "q_mapped_genes": q["genes"],
                    "reference_mapped_genes": reference_gene_count,
                    "p_nonzero_genes": p["nonzero"], "q_nonzero_genes": q["nonzero"],
                    "reference_nonzero_genes": reference_nonzero,
                    "p_mean_log1p_cpm": format_number(p["mean"]),
                    "q_mean_log1p_cpm": format_number(q["mean"]),
                    "reference_mean_log1p_cpm": format_number(reference_mean),
                    "p_score": format_number(p["score"]),
                    "q_score": format_number(q["score"]),
                    "p_fold0_mean_log1p_cpm": format_number(p["fold_means"][0]),
                    "p_fold1_mean_log1p_cpm": format_number(p["fold_means"][1]),
                    "q_fold0_mean_log1p_cpm": format_number(q["fold_means"][0]),
                    "q_fold1_mean_log1p_cpm": format_number(q["fold_means"][1]),
                    "reference_fold0_mean_log1p_cpm": format_number(
                        reference_fold_means[0]),
                    "reference_fold1_mean_log1p_cpm": format_number(
                        reference_fold_means[1]),
                    "p_fold0_mapped_genes": p["fold_genes"][0],
                    "p_fold1_mapped_genes": p["fold_genes"][1],
                    "q_fold0_mapped_genes": q["fold_genes"][0],
                    "q_fold1_mapped_genes": q["fold_genes"][1],
                    "reference_fold0_mapped_genes": int(reference_fold_genes[0]),
                    "reference_fold1_mapped_genes": int(reference_fold_genes[1]),
                    "p_fold0_score": format_number(p["fold_scores"][0]),
                    "p_fold1_score": format_number(p["fold_scores"][1]),
                    "q_fold0_score": format_number(q["fold_scores"][0]),
                    "q_fold1_score": format_number(q["fold_scores"][1]),
                    "p_top1_fraction": format_number(p["top"][0]),
                    "p_top5_fraction": format_number(p["top"][1]),
                    "p_top10_fraction": format_number(p["top"][2]),
                    "q_top1_fraction": format_number(q["top"][0]),
                    "q_top5_fraction": format_number(q["top"][1]),
                    "q_top10_fraction": format_number(q["top"][2]),
                    "reference_top1_fraction": fraction_or_na(
                        reference_top[1][cell_index, chromosome_index], reference_count),
                    "reference_top5_fraction": fraction_or_na(
                        reference_top[5][cell_index, chromosome_index], reference_count),
                    "reference_top10_fraction": fraction_or_na(
                        reference_top[10][cell_index, chromosome_index], reference_count),
                    "expressed_autosomal_genes": int(autosomal_nonzero[cell_index]),
                    "ambient_c": clean(manifest_row.get("ambient_c")) or "NA",
                    "matrix_value_type": matrix_value_type,
                    "expression_input_state": input_state,
                    "normalization_method": NORMALIZATION_METHOD,
                    "gene_fold_method": GENE_FOLD_METHOD,
                    "expression_scale_factor": f"{args.expression_scale_factor:.17g}",
                    "score_formula": MODEL_SCORE_FORMULA,
                    "reference_definition": MODEL_REFERENCE_DEFINITION,
                    "schema_version": EXPRESSION_MODEL_SCHEMA,
                }

    output_dir.mkdir(parents=True, exist_ok=True)
    written = write_tsv_atomic(str(evidence_path), expression_rows(), FIELDS)
    if written != row_count:
        raise RuntimeError(f"expression row-count mismatch: {written} != {row_count}")
    model_written = 0
    if args.emit_model_evidence:
        model_written = write_tsv_atomic(
            str(model_path), model_rows(), MODEL_FIELDS,
            deterministic_gzip=True)
        if model_written != model_row_count:
            raise RuntimeError(
                f"expression-model row-count mismatch: {model_written} != {model_row_count}")
    write_tsv_atomic(
        str(qc_path),
        [{
            "library": library, "manifest_cells": manifest_rows_total,
            "selected_manifest_cells": len(manifest),
            "matrix_cells": len(barcode_lines), "overlap_cells": len(selected),
            "missing_manifest_cells": len(missing), "matrix_genes": len(feature_lines),
            "mapped_genes": int(np.sum(row_arm_index >= 0)), "arms": len(arms),
            "rows": row_count, "status": "PASS",
            "schema_version": "tetra_arm_expression_qc_v1",
        }],
        ["library", "manifest_cells", "selected_manifest_cells", "matrix_cells",
         "overlap_cells", "missing_manifest_cells", "matrix_genes", "mapped_genes",
         "arms", "rows", "status", "schema_version"])
    contract = {
        "schema_version": "tetra_arm_expression_contract_v1",
        "release": RELEASE, "library": library,
        "inputs": {
            "barcodes": file_record(barcodes_path),
            "features": file_record(features_path),
            "matrix": file_record(matrix_path),
            "cell_manifest": file_record(legacy_manifest_path),
            "gene_arms": file_record(gene_arm_path),
        },
        "cells": len(selected), "arms": arms, "rows": row_count,
        "evidence_fields": FIELDS,
        "reference_autosomal_counts_definition": REFERENCE_DEFINITION,
        "status": "PASS", "expression_input_state": input_state,
    }
    if args.emit_model_evidence:
        contract["inputs"]["hybrid_cell_manifest"] = file_record(
            hybrid_manifest_path)
        contract["model_evidence"] = {
            "path": str(model_path), "schema_version": EXPRESSION_MODEL_SCHEMA,
            "rows": model_written, "fields": MODEL_FIELDS,
            "chromosomes": chromosomes,
            "normalization_method": NORMALIZATION_METHOD,
            "gene_fold_method": GENE_FOLD_METHOD,
            "gene_folds": args.gene_folds,
            "scale_factor": args.expression_scale_factor, "epsilon": 1e-9,
            "score_formula": MODEL_SCORE_FORMULA,
            "reference_definition": MODEL_REFERENCE_DEFINITION,
            "raw_and_corrected_rows_mixed": False,
        }
    write_json_atomic(str(contract_path), contract)
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
    parser.add_argument(
        "--emit-model-evidence", action="store_true",
        help="emit independent chromosome-level p/q hybrid expression evidence")
    parser.add_argument(
        "--hybrid-cell-manifest", default="",
        help="absolute hybrid manifest used for target/reference eligibility")
    parser.add_argument("--expression-scale-factor", type=float, default=10000.0)
    parser.add_argument("--gene-folds", type=int, default=2)
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
