#!/usr/bin/env python3
"""Join one library's hybrid inputs and shard them by chromosome."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Mapping

from tetra_arm_common import (
    ASE_SCHEMA,
    EXPRESSION_MODEL_SCHEMA,
    EXPRESSION_SCHEMA,
    HYBRID_CELL_MANIFEST_SCHEMA,
    HYBRID_SHARD_SCHEMA,
    canonical_barcode,
    clean,
    file_record,
    finite_float,
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
from tetra_arm_hybrid_common import (
    PROGRAM_VERSION,
    arm_side,
    biological_block,
    canonical_chromosome,
    chromosome_from_arm,
    is_autosomal_chromosome,
    join_flags,
    present_identifier,
    stable_fold,
    validate_hybrid_table,
)


INPUT_HEADER = (
    "library", "hybrid_cell_manifest", "ase", "expression", "expression_model")

ASE_CALIBRATION_FIELDS = (
    "library", "barcode", "ase_calibration_observations_ref",
    "ase_calibration_observations_alt", "ase_calibration_observations_mixed",
    "schema_version",
)
ASE_CALIBRATION_SCHEMA = "tetra_arm_hybrid_ase_calibration_payload_v1"

ASE_REQUIRED = {
    "library", "barcode", "donor_a", "donor_b", "donor_pair", "arm",
    "chromosome", "ambient_c", "ambient_c_se", "n_sites",
    "n_molecules", "n_informative_molecules",
    "effective_a_ref", "effective_b_ref", "effective_weight_ref",
    "effective_a_alt", "effective_b_alt", "effective_weight_alt",
    "effective_a_mixed", "effective_b_mixed", "effective_weight_mixed",
    "ambient_a_ref", "ambient_a_alt", "ambient_a_mixed",
    "ambient_genotyped_mass", "qname_fallback_fraction", "evidence_basis",
    "model_eligible", "evidence_status", "schema_version",
}

EXPRESSION_REQUIRED = {
    "library", "barcode", "arm", "chromosome", "arm_counts",
    "total_autosomal_counts", "mapped_genes_on_arm", "nonzero_genes_on_arm",
    "matrix_value_type", "expression_input_state", "schema_version",
}

ASE_COPY_FIELDS = (
    "arm_start", "arm_end", "ambient_c", "ambient_c_se", "n_sites",
    "n_molecules", "n_informative_molecules", "effective_a_ref",
    "effective_b_ref", "effective_weight_ref", "effective_a_alt",
    "effective_b_alt", "effective_weight_alt", "effective_a_mixed",
    "effective_b_mixed", "effective_weight_mixed", "ambient_a_ref",
    "ambient_a_alt", "ambient_a_mixed", "ambient_genotyped_mass",
    "qname_fallback_fraction", "mean_sites_per_molecule", "evidence_basis",
    "model_eligible", "evidence_status",
)

EXPRESSION_COPY_FIELDS = (
    "arm_counts", "other_autosomal_counts", "reference_autosomal_counts",
    "total_autosomal_counts", "arm_fraction", "log2_arm_to_other",
    "log2_arm_to_reference", "mapped_genes_on_arm", "nonzero_genes_on_arm",
    "matrix_value_type", "expression_input_state",
)

MODEL_COPY_FIELDS = (
    "p_arm", "q_arm", "p_raw_count_total", "q_raw_count_total",
    "reference_raw_count_total", "mapped_autosomal_library_size",
    "p_mapped_genes", "q_mapped_genes", "reference_mapped_genes",
    "p_nonzero_genes", "q_nonzero_genes", "reference_nonzero_genes",
    "p_mean_log1p_cpm", "q_mean_log1p_cpm", "reference_mean_log1p_cpm",
    "p_score", "q_score", "p_fold0_mean_log1p_cpm",
    "p_fold1_mean_log1p_cpm", "q_fold0_mean_log1p_cpm",
    "q_fold1_mean_log1p_cpm", "reference_fold0_mean_log1p_cpm",
    "reference_fold1_mean_log1p_cpm", "p_fold0_mapped_genes",
    "p_fold1_mapped_genes", "q_fold0_mapped_genes", "q_fold1_mapped_genes",
    "reference_fold0_mapped_genes", "reference_fold1_mapped_genes",
    "p_fold0_score", "p_fold1_score", "q_fold0_score", "q_fold1_score",
    "p_top1_fraction", "p_top5_fraction", "p_top10_fraction",
    "q_top1_fraction", "q_top5_fraction", "q_top10_fraction",
    "reference_top1_fraction", "reference_top5_fraction",
    "reference_top10_fraction", "expressed_autosomal_genes", "ambient_c",
    "matrix_value_type", "expression_input_state", "normalization_method",
    "gene_fold_method", "expression_scale_factor", "score_formula",
    "reference_definition",
)

# Every field consumed below is structural for the corresponding producer
# schema. Additional columns remain allowed by the required-subset validator.
ASE_REQUIRED.update(ASE_COPY_FIELDS)
ASE_REQUIRED.update({
    "n_molecules_ref", "n_molecules_alt", "n_molecules_mixed",
    "a_ref", "b_ref", "a_alt", "b_alt", "a_mixed", "b_mixed",
    "n_ambiguous", "soft_a", "soft_b", "soft_a_ref", "soft_b_ref",
    "soft_a_alt", "soft_b_alt", "soft_a_mixed", "soft_b_mixed",
    "soft_a_sumsq_ref", "soft_a_sumsq_alt", "soft_a_sumsq_mixed",
})
EXPRESSION_REQUIRED.update(EXPRESSION_COPY_FIELDS)
MODEL_REQUIRED = {
    "library", "barcode", "chromosome", "p_arm", "q_arm", "donor_a",
    "donor_b", "donor_pair", "uid", "calibration_group",
    "hybrid_target_eligible", "expression_reference_eligible",
    "cell_group_source", "cell_group_target_chromosome_excluded",
    *MODEL_COPY_FIELDS, "schema_version",
}

SHARD_FIELDS = [
    "library", "barcode", "chromosome", "arm", "arm_side", "uid",
    "donor_a", "donor_b", "donor_pair", "calibration_group",
    "hybrid_target_eligible", "expression_reference_eligible",
    "ase_calibration_eligible", "cell_group_source",
    "cell_group_target_chromosome_excluded", "biological_block",
    "biological_block_source", "crossfit_fold", "ase_present",
    "expression_present", "expression_model_present", "ase_status",
    "expression_status", "ase_depth_bin", "ase_ambient_bin",
    "expression_depth_bin", "expression_breadth_bin",
    *(f"ase_{field}" for field in ASE_COPY_FIELDS),
    *(f"legacy_expression_{field}" for field in EXPRESSION_COPY_FIELDS),
    *(f"expression_model_{field}" for field in MODEL_COPY_FIELDS),
    "ase_loco_arms", "ase_loco_effective_a_ref",
    "ase_loco_effective_weight_ref", "ase_loco_ambient_a_ref",
    "ase_loco_effective_a_alt", "ase_loco_effective_weight_alt",
    "ase_loco_ambient_a_alt", "ase_loco_effective_a_mixed",
    "ase_loco_effective_weight_mixed", "ase_loco_ambient_a_mixed",
    "ase_calibration_observations_ref",
    "ase_calibration_observations_alt",
    "ase_calibration_observations_mixed",
    "confounding_flags", "qc_flags", "schema_version",
]

NATIVE_ASE_CALIBRATION_POLICY = {
    "require_evidence_status": "PASS",
    "require_model_eligible": True,
    "reject_site_fallback": True,
    "min_sites": 2,
    "min_effective_weight": 4.0,
    "max_qname_fallback_fraction": 0.50,
    "min_ambient_genotyped_mass": 0.50,
    "require_positive_ambient_c_se": True,
}

UNIT_INTERVAL_TOLERANCE = 1e-12


def require_absolute(path: str, label: str) -> str:
    if not os.path.isabs(path):
        raise ValueError(f"{label} must be absolute: {path}")
    return require_file(path, label)


def load_input_row(path: str, library: int) -> dict[str, str]:
    target = require_file(path, "hybrid input manifest")
    with open(target, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != INPUT_HEADER:
            raise ValueError(
                f"hybrid input manifest header must be {INPUT_HEADER}: {target}")
        matches = [row for row in reader
                   if clean(row.get("library")).lower().removeprefix("lib")
                   == str(library)]
    if len(matches) != 1:
        raise ValueError(
            f"hybrid input manifest must contain exactly one lib{library} row")
    row = matches[0]
    for field in INPUT_HEADER[1:]:
        row[field] = require_absolute(row[field], field)
    return row


def load_manifest(path: str, library: int) -> dict[str, dict[str, str]]:
    validate_hybrid_table(path, HYBRID_CELL_MANIFEST_SCHEMA)
    result: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        observed_library = clean(row.get("library")).lower().removeprefix("lib")
        barcode = canonical_barcode(row.get("barcode", ""))
        if observed_library != str(library) or not barcode or barcode in result:
            raise ValueError(f"wrong library or duplicate hybrid manifest key: {path}")
        row["library"], row["barcode"] = str(library), barcode
        if truthy(row.get("hybrid_target_eligible")):
            donor_a = present_identifier(row.get("donor_a"))
            donor_b = present_identifier(row.get("donor_b"))
            donor_pair = present_identifier(row.get("donor_pair"))
            if not donor_a or not donor_b or donor_a == donor_b or not donor_pair:
                raise ValueError(
                    f"target-eligible cell lacks an ordered heterotypic pair: "
                    f"lib{library}/{barcode}")
        result[barcode] = row
    return result


def load_keyed_rows(
        path: str, library: int, required: set[str], schema: str,
        key_kind: str) -> tuple[dict[tuple[str, ...], dict[str, str]], list[str]]:
    header, _count = validate_tsv_schema(path, schema, required)
    result: dict[tuple[str, ...], dict[str, str]] = {}
    for row in read_tsv(path):
        observed_library = clean(row.get("library")).lower().removeprefix("lib")
        barcode = canonical_barcode(row.get("barcode", ""))
        if observed_library != str(library) or not barcode:
            raise ValueError(f"wrong library or empty barcode: {path}")
        row["library"], row["barcode"] = str(library), barcode
        if key_kind == "chromosome":
            chromosome = canonical_chromosome(row.get("chromosome"))
            if not chromosome:
                raise ValueError(f"empty chromosome in {path}")
            row["chromosome"] = chromosome
            key = (barcode, chromosome)
        else:
            arm = clean(row.get("arm"))
            if not arm:
                raise ValueError(f"empty arm in {path}")
            row["chromosome"] = canonical_chromosome(
                row.get("chromosome") or chromosome_from_arm(arm))
            key = (barcode, arm)
        if key in result:
            raise ValueError(f"duplicate {key_kind} key {key}: {path}")
        result[key] = row
    return result, header


def parse_chromosomes(values: list[str]) -> list[str]:
    result: set[str] = set()
    for value in values:
        for token in str(value).split(","):
            chromosome = canonical_chromosome(token)
            if chromosome:
                result.add(chromosome)
    return sorted(result, key=natural_key)


def validate_ase_values(rows: Mapping[tuple[str, ...], dict[str, str]],
                        path: str) -> None:
    """Reject malformed soft-evidence identities before any aggregation."""
    for key, row in rows.items():
        context = f"{path}:{'/'.join(key)}"
        count_fields = (
            "n_sites", "n_molecules", "n_informative_molecules",
            "n_molecules_ref", "n_molecules_alt", "n_molecules_mixed",
            "a_ref", "b_ref", "a_alt", "b_alt", "a_mixed", "b_mixed",
            "n_ambiguous",
        )
        counts: dict[str, int] = {}
        for field in count_fields:
            value = finite_float(row.get(field))
            if not math.isfinite(value) or value < 0.0 or not value.is_integer():
                raise ValueError(f"invalid nonnegative integer {field}: {context}")
            counts[field] = int(value)
        if sum(counts[field] for field in (
                "a_ref", "b_ref", "a_alt", "b_alt", "a_mixed", "b_mixed")) \
                != counts["n_informative_molecules"]:
            raise ValueError(f"informative hard-count mismatch: {context}")
        contamination = finite_float(row.get("ambient_c"))
        if not math.isfinite(contamination) or not 0.0 <= contamination <= 1.0:
            raise ValueError(f"ambient_c is outside [0,1]: {context}")
        standard_error = finite_float(row.get("ambient_c_se"))
        if math.isfinite(standard_error) and standard_error < 0.0:
            raise ValueError(f"ambient_c_se is negative: {context}")
        qname_fraction = finite_float(row.get("qname_fallback_fraction"))
        if (not math.isfinite(qname_fraction)
                or not 0.0 <= qname_fraction <= 1.0):
            raise ValueError(
                f"qname_fallback_fraction is outside [0,1]: {context}")
        total_effective_weight = 0.0
        for orientation in ("ref", "alt", "mixed"):
            units = counts[f"n_molecules_{orientation}"]
            soft_a = finite_float(row.get(f"soft_a_{orientation}"))
            soft_b = finite_float(row.get(f"soft_b_{orientation}"))
            sumsq = finite_float(row.get(f"soft_a_sumsq_{orientation}"))
            effective_a = finite_float(row.get(f"effective_a_{orientation}"))
            effective_b = finite_float(row.get(f"effective_b_{orientation}"))
            weight = finite_float(row.get(f"effective_weight_{orientation}"))
            tolerance = max(1e-7, 1e-6 * max(units, 1))
            if (not all(math.isfinite(value) and value >= 0.0 for value in (
                    soft_a, soft_b, sumsq))
                    or not math.isclose(
                        soft_a + soft_b, units,
                        rel_tol=1e-6, abs_tol=tolerance)):
                raise ValueError(
                    f"inconsistent soft molecule evidence ({orientation}): "
                    f"{context}")
            lower_sumsq = soft_a * soft_a / units if units else 0.0
            if sumsq < lower_sumsq - tolerance or sumsq > soft_a + tolerance:
                raise ValueError(
                    f"invalid soft-a sumsq ({orientation}): {context}")
            if (not all(math.isfinite(value) and value >= 0.0 for value in (
                    effective_a, effective_b, weight))
                    or not math.isclose(
                        effective_a + effective_b, weight,
                        rel_tol=1e-6, abs_tol=tolerance)
                    or weight > units + tolerance):
                raise ValueError(
                    f"inconsistent soft effective evidence ({orientation}): "
                    f"{context}")
            identity_weight = 4.0 * sumsq - 4.0 * soft_a + units
            if not math.isclose(
                    weight, identity_weight,
                    rel_tol=1e-6, abs_tol=tolerance):
                raise ValueError(
                    f"effective-weight identity mismatch ({orientation}): "
                    f"{context}")
            total_effective_weight += weight
            ambient = finite_float(row.get(f"ambient_a_{orientation}"))
            if weight > 0.0 and (
                    not math.isfinite(ambient)
                    or ambient < -UNIT_INTERVAL_TOLERANCE
                    or ambient > 1.0 + UNIT_INTERVAL_TOLERANCE):
                raise ValueError(
                    f"invalid ambient_a_{orientation}: {context}")
            if weight > 0.0:
                row[f"ambient_a_{orientation}"] = (
                    f"{min(1.0, max(0.0, ambient)):.17g}")
        genotyped_mass = finite_float(row.get("ambient_genotyped_mass"))
        if math.isfinite(genotyped_mass):
            if (genotyped_mass < -UNIT_INTERVAL_TOLERANCE
                    or genotyped_mass > 1.0 + UNIT_INTERVAL_TOLERANCE):
                raise ValueError(
                    f"ambient_genotyped_mass is outside [0,1]: {context}")
            row["ambient_genotyped_mass"] = (
                f"{min(1.0, max(0.0, genotyped_mass)):.17g}")
        elif total_effective_weight > 0.0:
            raise ValueError(
                "ambient_genotyped_mass is unavailable with positive "
                f"effective weight: {context}")
        global_soft_a = finite_float(row.get("soft_a"))
        global_soft_b = finite_float(row.get("soft_b"))
        orientation_soft_a = sum(finite_float(
            row.get(f"soft_a_{orientation}"), 0.0)
            for orientation in ("ref", "alt", "mixed"))
        orientation_soft_b = sum(finite_float(
            row.get(f"soft_b_{orientation}"), 0.0)
            for orientation in ("ref", "alt", "mixed"))
        if (not math.isfinite(global_soft_a) or not math.isfinite(global_soft_b)
                or not math.isclose(
                    global_soft_a, orientation_soft_a,
                    rel_tol=1e-6, abs_tol=1e-6)
                or not math.isclose(
                    global_soft_b, orientation_soft_b,
                    rel_tol=1e-6, abs_tol=1e-6)):
            raise ValueError(f"global/orientation soft-count mismatch: {context}")


def ase_calibration_row_usable(row: Mapping[str, object]) -> bool:
    if (clean(row.get("evidence_status")) != "PASS"
            or not truthy(row.get("model_eligible"))
            or not is_autosomal_chromosome(row.get("chromosome"))):
        return False
    contamination = finite_float(row.get("ambient_c"))
    standard_error = finite_float(row.get("ambient_c_se"))
    total_weight = sum(max(0.0, finite_float(
        row.get(f"effective_weight_{orientation}"), 0.0))
        for orientation in ("ref", "alt", "mixed"))
    return (
        math.isfinite(contamination) and contamination < 1.0
        and finite_float(row.get("n_sites"), 0.0)
        >= NATIVE_ASE_CALIBRATION_POLICY["min_sites"]
        and total_weight
        >= NATIVE_ASE_CALIBRATION_POLICY["min_effective_weight"]
        and finite_float(row.get("qname_fallback_fraction"), 1.0)
        <= NATIVE_ASE_CALIBRATION_POLICY["max_qname_fallback_fraction"]
        and finite_float(row.get("ambient_genotyped_mass"), 0.0)
        >= NATIVE_ASE_CALIBRATION_POLICY["min_ambient_genotyped_mass"]
        and math.isfinite(standard_error) and standard_error > 0.0)


def ase_target_safety_status(row: Mapping[str, object] | None) -> str:
    """Apply native high-confidence ASE safety gates before target scoring."""
    if row is None:
        return "NO_DATA"
    evidence_status = clean(row.get("evidence_status")) or "NO_DATA"
    evidence_basis = clean(row.get("evidence_basis")).upper()
    if evidence_status == "PASS_SITE_FALLBACK" or "SITE" in evidence_basis:
        return "UNSAFE_SITE_FALLBACK"
    if evidence_status != "PASS":
        return evidence_status
    contamination = finite_float(row.get("ambient_c"))
    if math.isfinite(contamination) and contamination >= 1.0:
        return "UNSAFE_AMBIENT_C_AT_UPPER_BOUNDARY"
    if (finite_float(row.get("qname_fallback_fraction"), 1.0) >
            NATIVE_ASE_CALIBRATION_POLICY["max_qname_fallback_fraction"]):
        return "UNSAFE_HIGH_QNAME_FALLBACK"
    if (finite_float(row.get("ambient_genotyped_mass"), 0.0) <
            NATIVE_ASE_CALIBRATION_POLICY["min_ambient_genotyped_mass"]):
        return "UNSAFE_LOW_AMBIENT_GENOTYPED_MASS"
    standard_error = finite_float(row.get("ambient_c_se"))
    if not math.isfinite(standard_error) or standard_error <= 0.0:
        return "UNSAFE_AMBIENT_UNCERTAINTY_UNAVAILABLE"
    return "PASS"


def calibration_observation_payload(
        ase_by_cell: Mapping[str, list[dict[str, str]]],
        barcode: str) -> dict[str, str]:
    """Serialize native-scale cell-arm observations once per library shard.

    Mapping-bias calibration pools these observations once per cell, whereas
    rho and the cell LOCO baseline retain the individual cell-arm records.
    The model reads this payload from one carrier chromosome instead of
    repeating a genome-wide list in every chromosome shard.
    """
    retained = sorted([
        row for row in ase_by_cell.get(barcode, [])
        if ase_calibration_row_usable(row)
    ], key=lambda row: (
        natural_key(canonical_chromosome(row.get("chromosome"))),
        natural_key(row.get("arm", ""))))
    result: dict[str, str] = {}
    for orientation in ("ref", "alt", "mixed"):
        observations: list[list[str]] = []
        for row in retained:
            weight = max(0.0, finite_float(
                row.get(f"effective_weight_{orientation}"), 0.0))
            if weight <= 0.0:
                continue
            effective_a = max(0.0, finite_float(
                row.get(f"effective_a_{orientation}"), 0.0))
            contamination = min(0.999999, max(
                0.0, finite_float(row.get("ambient_c"), 0.0)))
            ambient = min(1.0, max(0.0, finite_float(
                row.get(f"ambient_a_{orientation}"), 0.5)))
            observations.append([
                canonical_chromosome(row.get("chromosome")),
                clean(row.get("arm")),
                f"{effective_a / weight:.17g}", f"{weight:.17g}",
                f"{contamination:.17g}", f"{ambient:.17g}",
            ])
        result[f"ase_calibration_observations_{orientation}"] = json.dumps(
            observations, ensure_ascii=True, separators=(",", ":"))
    return result


def loco_summary(ase_by_cell: Mapping[str, list[dict[str, str]]],
                 barcode: str, chromosome: str) -> dict[str, object]:
    retained = sorted([
        row for row in ase_by_cell.get(barcode, [])
        if canonical_chromosome(row.get("chromosome")) != chromosome
        and ase_calibration_row_usable(row)
    ], key=lambda row: (
        natural_key(canonical_chromosome(row.get("chromosome"))),
        natural_key(row.get("arm", ""))))
    result: dict[str, object] = {
        "ase_loco_arms": len({clean(row.get("arm")) for row in retained}),
    }
    for orientation in ("ref", "alt", "mixed"):
        effective_a = 0.0
        effective_weight = 0.0
        ambient_weighted = 0.0
        for row in retained:
            weight = max(0.0, finite_float(
                row.get(f"effective_weight_{orientation}"), 0.0))
            effective_a += max(0.0, finite_float(
                row.get(f"effective_a_{orientation}"), 0.0))
            effective_weight += weight
            ambient_weighted += weight * min(1.0, max(
                0.0, finite_float(row.get(f"ambient_a_{orientation}"), 0.5)))
        result[f"ase_loco_effective_a_{orientation}"] = f"{effective_a:.17g}"
        result[f"ase_loco_effective_weight_{orientation}"] = (
            f"{effective_weight:.17g}")
        result[f"ase_loco_ambient_a_{orientation}"] = (
            f"{ambient_weighted / effective_weight:.17g}"
            if effective_weight > 0.0 else "NA")
    return result


def main_impl(args) -> int:
    library = int(args.library)
    if args.crossfit_folds < 2:
        raise ValueError("--crossfit-folds must be at least two")
    input_row = load_input_row(args.input_manifest, library)
    manifest = load_manifest(input_row["hybrid_cell_manifest"], library)
    ase, _ase_header = load_keyed_rows(
        input_row["ase"], library, ASE_REQUIRED, ASE_SCHEMA, "arm")
    validate_ase_values(ase, input_row["ase"])
    expression, _expression_header = load_keyed_rows(
        input_row["expression"], library, EXPRESSION_REQUIRED,
        EXPRESSION_SCHEMA, "arm")
    expression_model, _model_header = load_keyed_rows(
        input_row["expression_model"], library,
        MODEL_REQUIRED,
        EXPRESSION_MODEL_SCHEMA, "chromosome")

    for source_name, rows in (("ASE", ase), ("expression", expression),
                              ("expression model", expression_model)):
        for key, row in rows.items():
            barcode = key[0]
            if barcode not in manifest:
                raise ValueError(
                    f"{source_name} cell is absent from hybrid manifest: "
                    f"lib{library}/{barcode}")
            if source_name == "ASE":
                manifest_row = manifest[barcode]
                if (clean(row.get("donor_a")) != clean(manifest_row.get("donor_a"))
                        or clean(row.get("donor_b"))
                        != clean(manifest_row.get("donor_b"))
                        or truthy(row.get("model_eligible")) != truthy(
                            manifest_row.get("hybrid_target_eligible"))):
                    raise ValueError(
                        f"ASE/manifest identity or eligibility mismatch for "
                        f"lib{library}/{barcode}")
                for ambient_field in ("ambient_c", "ambient_c_se"):
                    ase_value = finite_float(row.get(ambient_field))
                    manifest_value = finite_float(manifest_row.get(ambient_field))
                    if (math.isfinite(ase_value) != math.isfinite(manifest_value)
                            or (math.isfinite(ase_value) and not math.isclose(
                                ase_value, manifest_value,
                                rel_tol=1e-8, abs_tol=1e-10))):
                        raise ValueError(
                            f"{ambient_field} mismatch for "
                            f"lib{library}/{barcode}")
            elif source_name == "expression model":
                manifest_row = manifest[barcode]
                if (clean(row.get("donor_a")) != clean(
                        manifest_row.get("donor_a"))
                        or clean(row.get("donor_b")) != clean(
                            manifest_row.get("donor_b"))
                        or clean(row.get("donor_pair")) != clean(
                            manifest_row.get("donor_pair"))
                        or clean(row.get("uid")) != clean(
                            manifest_row.get("uid"))
                        or clean(row.get("calibration_group")) != clean(
                            manifest_row.get("calibration_group"))
                        or truthy(row.get("hybrid_target_eligible")) != truthy(
                            manifest_row.get("hybrid_target_eligible"))
                        or truthy(row.get("expression_reference_eligible"))
                        != truthy(manifest_row.get(
                            "expression_reference_eligible"))
                        or clean(row.get("cell_group_source")) != clean(
                            manifest_row.get("cell_group_source"))
                        or truthy(row.get(
                            "cell_group_target_chromosome_excluded")) != truthy(
                                manifest_row.get(
                                    "cell_group_target_chromosome_excluded"))):
                    raise ValueError(
                        f"expression-model manifest join mismatch for "
                        f"lib{library}/{barcode}")

    universe: set[tuple[str, str, str]] = set()
    for (barcode, arm), row in ase.items():
        universe.add((barcode, canonical_chromosome(row.get("chromosome")), arm))
    for (barcode, arm), row in expression.items():
        universe.add((barcode, canonical_chromosome(row.get("chromosome")), arm))
    for (barcode, chromosome), row in expression_model.items():
        for field in ("p_arm", "q_arm"):
            arm = present_identifier(row.get(field))
            if arm:
                observed = canonical_chromosome(chromosome_from_arm(arm))
                if observed != chromosome:
                    raise ValueError(
                        f"expression-model arm/chromosome mismatch: {arm}/{chromosome}")
                universe.add((barcode, chromosome, arm))

    discovered = {chromosome for _barcode, chromosome, _arm in universe
                  if chromosome}
    requested = parse_chromosomes(args.chromosomes)
    chromosomes = requested or sorted(discovered, key=natural_key)
    unexpected = sorted(discovered - set(chromosomes), key=natural_key)
    if unexpected:
        raise ValueError(
            "inputs contain chromosomes absent from --chromosomes: "
            + ",".join(unexpected))
    eligible_barcodes = {
        barcode for barcode, row in manifest.items()
        if truthy(row.get("hybrid_target_eligible"))
        or truthy(row.get("expression_reference_eligible"))
    }
    missing_model_keys = [
        (barcode, chromosome) for barcode in sorted(
            eligible_barcodes, key=natural_key)
        for chromosome in chromosomes
        if (barcode, chromosome) not in expression_model
    ]
    if missing_model_keys:
        raise ValueError(
            "eligible hybrid cells are missing chromosome expression-model "
            f"rows; first={missing_model_keys[:5]}")

    output_dir = Path(os.path.abspath(args.output_dir))
    qc_path = output_dir / f"lib{library}.hybrid_shard_qc.tsv"
    contract_path = output_dir / f"lib{library}.hybrid_shard_contract.json"
    ase_calibration_path = (
        output_dir / f"lib{library}.ase_calibration_payload.tsv.gz")
    shard_paths = {
        chromosome: output_dir / f"lib{library}.chr{chromosome}.hybrid_shard.tsv.gz"
        for chromosome in chromosomes
    }
    require_outputs_absent([
        qc_path, contract_path, ase_calibration_path, *shard_paths.values()])

    ase_by_cell: dict[str, list[dict[str, str]]] = defaultdict(list)
    for (barcode, _arm), row in ase.items():
        ase_by_cell[barcode].append(row)

    calibration_carrier_chromosome = chromosomes[0] if chromosomes else ""
    calibration_payloads = {
        barcode: calibration_observation_payload(ase_by_cell, barcode)
        for barcode in sorted(ase_by_cell, key=natural_key)
    }

    rows_by_chromosome: dict[str, list[dict[str, object]]] = {
        chromosome: [] for chromosome in chromosomes}
    status_counts = Counter()
    for barcode, chromosome, arm in sorted(
            universe, key=lambda value: tuple(natural_key(part) for part in value)):
        if chromosome not in rows_by_chromosome:
            continue
        manifest_row = manifest[barcode]
        ase_row = ase.get((barcode, arm))
        expression_row = expression.get((barcode, arm))
        model_row = expression_model.get((barcode, chromosome))
        if ase_row is not None and canonical_chromosome(
                ase_row.get("chromosome")) != chromosome:
            raise ValueError(f"ASE chromosome mismatch for lib{library}/{barcode}/{arm}")
        if expression_row is not None and canonical_chromosome(
                expression_row.get("chromosome")) != chromosome:
            raise ValueError(
                f"expression chromosome mismatch for lib{library}/{barcode}/{arm}")
        block, block_source, block_flags = biological_block(manifest_row)
        side = arm_side(arm)
        score = finite_float(model_row.get(f"{side}_score")) if (
            model_row is not None and side) else math.nan
        expression_status = "PASS" if math.isfinite(score) else "NO_DATA"
        ase_status = ase_target_safety_status(ase_row)
        depth = (finite_float(
            model_row.get("mapped_autosomal_library_size"), 0.0)
            if model_row is not None else 0.0)
        breadth = (finite_float(
            model_row.get(f"{side}_nonzero_genes"), 0.0)
            if model_row is not None and side else 0.0)
        ase_weight = (sum(max(0.0, finite_float(
            ase_row.get(f"effective_weight_{orientation}"), 0.0))
            for orientation in ("ref", "alt", "mixed"))
            if ase_row is not None else 0.0)
        ase_ambient = (min(0.999999, max(
            0.0, finite_float(ase_row.get("ambient_c"), 0.0)))
            if ase_row is not None else 0.0)
        row: dict[str, object] = {
            "library": library, "barcode": barcode, "chromosome": chromosome,
            "arm": arm, "arm_side": side or "NA",
            "uid": clean(manifest_row.get("uid")) or "NA",
            "donor_a": clean(manifest_row.get("donor_a")) or "NA",
            "donor_b": clean(manifest_row.get("donor_b")) or "NA",
            "donor_pair": clean(manifest_row.get("donor_pair")) or "NA",
            "calibration_group": clean(
                manifest_row.get("calibration_group")) or f"lib{library}",
            "hybrid_target_eligible": int(truthy(
                manifest_row.get("hybrid_target_eligible"))),
            "expression_reference_eligible": int(truthy(
                manifest_row.get("expression_reference_eligible"))),
            "ase_calibration_eligible": int(truthy(
                manifest_row.get("calibration_eligible"))),
            "cell_group_source": clean(
                manifest_row.get("cell_group_source")) or "LIBRARY_FALLBACK",
            "cell_group_target_chromosome_excluded": int(truthy(
                manifest_row.get("cell_group_target_chromosome_excluded"))),
            "biological_block": block, "biological_block_source": block_source,
            "crossfit_fold": stable_fold(block, args.crossfit_folds,
                                         b"tetra-cell-fold"),
            "ase_present": int(ase_row is not None),
            "expression_present": int(expression_row is not None),
            "expression_model_present": int(model_row is not None),
            "ase_status": ase_status, "expression_status": expression_status,
            "ase_depth_bin": int(math.floor(math.log2(ase_weight + 1.0))),
            "ase_ambient_bin": int(math.floor(ase_ambient / 0.10)),
            "expression_depth_bin": int(math.floor(math.log2(depth + 1.0))),
            "expression_breadth_bin": int(math.floor(math.log2(breadth + 1.0))),
        }
        for field in ASE_COPY_FIELDS:
            row[f"ase_{field}"] = (
                clean(ase_row.get(field)) or "NA" if ase_row is not None else "NA")
        for field in EXPRESSION_COPY_FIELDS:
            row[f"legacy_expression_{field}"] = (
                clean(expression_row.get(field)) or "NA"
                if expression_row is not None else "NA")
        for field in MODEL_COPY_FIELDS:
            row[f"expression_model_{field}"] = (
                clean(model_row.get(field)) or "NA"
                if model_row is not None else "NA")
        row.update(loco_summary(ase_by_cell, barcode, chromosome))
        for orientation in ("ref", "alt", "mixed"):
            field = f"ase_calibration_observations_{orientation}"
            # Genome-wide payloads live once in the compact sidecar below.
            # Keeping them out of a wide chromosome row prevents every MODEL
            # task from rereading hundreds of unrelated columns per cell.
            row[field] = "NA"
        confounding = []
        if not truthy(manifest_row.get("cell_group_target_chromosome_excluded")):
            confounding.append("UNADJUSTED_CELL_STATE")
        row["confounding_flags"] = join_flags(confounding)
        row["qc_flags"] = join_flags(block_flags)
        row["schema_version"] = HYBRID_SHARD_SCHEMA
        rows_by_chromosome[chromosome].append(row)
        status_counts[(ase_status, expression_status)] += 1

    output_dir.mkdir(parents=True, exist_ok=True)
    calibration_rows = []
    for barcode in sorted(calibration_payloads, key=natural_key):
        payload = calibration_payloads[barcode]
        if not any(payload.get(
                f"ase_calibration_observations_{orientation}", "[]") != "[]"
                for orientation in ("ref", "alt", "mixed")):
            continue
        calibration_rows.append({
            "library": library, "barcode": barcode,
            **{f"ase_calibration_observations_{orientation}": payload.get(
                f"ase_calibration_observations_{orientation}", "[]")
               for orientation in ("ref", "alt", "mixed")},
            "schema_version": ASE_CALIBRATION_SCHEMA,
        })
    calibration_count = write_tsv_atomic(
        str(ase_calibration_path), calibration_rows, ASE_CALIBRATION_FIELDS,
        deterministic_gzip=True)
    written_by_chromosome: dict[str, int] = {}
    for chromosome in chromosomes:
        rows = sorted(
            rows_by_chromosome[chromosome],
            key=lambda row: (
                natural_key(row["library"]), natural_key(row["barcode"]),
                natural_key(row["arm"])))
        written_by_chromosome[chromosome] = write_tsv_atomic(
            str(shard_paths[chromosome]), rows, SHARD_FIELDS,
            deterministic_gzip=True)

    total_rows = sum(written_by_chromosome.values())
    terminal = "NONE" if total_rows else "PASS_NO_ELIGIBLE_ROWS"
    qc_rows = [
        {"metric": "schema_version", "value": "tetra_arm_hybrid_shard_qc_v1"},
        {"metric": "library", "value": library},
        {"metric": "chromosomes", "value": len(chromosomes)},
        {"metric": "rows", "value": total_rows},
        {"metric": "ase_calibration_cells", "value": calibration_count},
        {"metric": "terminal_state", "value": terminal},
        {"metric": "status", "value": "PASS"},
    ]
    write_tsv_atomic(str(qc_path), qc_rows, ["metric", "value"])
    write_json_atomic(str(contract_path), {
        "schema_version": "tetra_arm_hybrid_shard_contract_v1",
        "release": PROGRAM_VERSION,
        "library": library,
        "inputs": {field: file_record(input_row[field]) for field in INPUT_HEADER[1:]},
        "outputs": {
            "shards": {chromosome: str(shard_paths[chromosome])
                       for chromosome in chromosomes},
            "qc": str(qc_path), "contract": str(contract_path),
            "ase_calibration": str(ase_calibration_path),
        },
        "output_schema": HYBRID_SHARD_SCHEMA,
        "fields": SHARD_FIELDS,
        "chromosome_rows": written_by_chromosome,
        "rows": total_rows,
        "crossfit_folds": args.crossfit_folds,
        "ase_loco_calibration_policy": NATIVE_ASE_CALIBRATION_POLICY,
        "ase_calibration_carrier_chromosome": calibration_carrier_chromosome,
        "ase_calibration_payload_schema": ASE_CALIBRATION_SCHEMA,
        "ase_calibration_payload_rows": calibration_count,
        "ase_calibration_observation_fields": [
            "ase_calibration_observations_ref",
            "ase_calibration_observations_alt",
            "ase_calibration_observations_mixed",
        ],
        "terminal_state": terminal,
        "status": "PASS",
    })
    return 0


def self_test() -> int:
    block, source, flags = biological_block(
        {"uid": "U1", "donor_pair": "A+B", "library": "1", "barcode": "x"})
    if block != "UID:U1" or source != "UID" or flags:
        raise AssertionError("biological block self-test failed")
    if stable_fold(block, 5, b"tetra-cell-fold") != stable_fold(
            block, 5, b"tetra-cell-fold"):
        raise AssertionError("fold determinism self-test failed")
    print("PASS tetra_arm_hybrid_shard self-test")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Join one library's hybrid inputs and shard by chromosome.")
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {PROGRAM_VERSION}")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--input-manifest")
    parser.add_argument("--library", type=int)
    parser.add_argument("--output-dir")
    parser.add_argument("--chromosomes", nargs="*", default=[])
    parser.add_argument("--crossfit-folds", type=int, default=5)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.self_test:
            return self_test()
        missing = [name for name in ("input_manifest", "library", "output_dir")
                   if getattr(args, name) in {None, ""}]
        if missing:
            raise ValueError("missing required option(s): " + ", ".join(missing))
        return main_impl(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
