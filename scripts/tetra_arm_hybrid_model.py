#!/usr/bin/env python3
"""Fit one chromosome's independent expression and ASE hybrid-v1 branches."""

from __future__ import annotations

import argparse
from array import array
import bisect
import csv
import json
import math
import os
import sys
import time
from collections import Counter, OrderedDict, defaultdict
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence

import numpy as np

from tetra_arm_common import (
    HYBRID_CALIBRATION_SCHEMA,
    HYBRID_COMPONENT_SCHEMA,
    HYBRID_SCORE_SCHEMA,
    HYBRID_SHARD_SCHEMA,
    clean,
    file_record,
    finite_float,
    natural_key,
    open_text,
    read_tsv,
    require_file,
    require_outputs_absent,
    truthy,
    write_json_atomic,
    write_tsv_atomic,
)
from tetra_arm_hybrid_common import (
    ASE_DIRECTIONS,
    AXES_TO_STATE,
    BASELINE_LEVELS,
    COPY_STATES,
    EXACT_STATES,
    EXPRESSION_PATTERN_SHIFTS,
    EXPRESSION_PATTERNS,
    PROGRAM_VERSION,
    ambient_adjusted_fraction,
    arm_side,
    ase_direction_statistics,
    ase_state_log_bfs,
    best_expression_component,
    bounded,
    canonical_chromosome,
    conjunction_pvalues,
    direction_from_log_bfs,
    empirical_p_floor,
    expression_direction_log_bfs,
    fit_ase_bias_component,
    fit_fixed_center_mixture,
    format_number,
    is_autosomal_chromosome,
    join_flags,
    present_identifier,
    resolve_evidence,
    robust_bivariate_fit,
    split_flags,
    validate_hybrid_table,
)


SHARD_MANIFEST_HEADER = ("library", "shard", "qc", "contract")

ASE_CALIBRATION_PAYLOAD_FIELDS = tuple(
    f"ase_calibration_observations_{orientation}"
    for orientation in ("ref", "alt", "mixed"))
ASE_CALIBRATION_PAYLOAD_SCHEMA = "tetra_arm_hybrid_ase_calibration_payload_v1"
ASE_ORIENTATIONS = ("ref", "alt", "mixed")


def model_timings(args) -> dict[str, float]:
    timings = getattr(args, "_hybrid_model_timings", None)
    if timings is None:
        timings = defaultdict(float)
        setattr(args, "_hybrid_model_timings", timings)
    return timings


def parse_ase_calibration_payload(
        row: Mapping[str, object], path: str,
        ) -> dict[str, list[tuple[str, str, float, float, float, float]]]:
    result: dict[str, list[tuple[str, str, float, float, float, float]]] = {
        orientation: [] for orientation in ("ref", "alt", "mixed")}
    for orientation in result:
        raw = clean(row.get(f"ase_calibration_observations_{orientation}"))
        if not raw or raw.upper() in {"NA", "N/A", "NONE", "NULL", "."}:
            continue
        try:
            values = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise ValueError(f"malformed ASE calibration payload in {path}") from exc
        if not isinstance(values, list):
            raise ValueError(f"ASE calibration payload is not a list in {path}")
        previous_key: tuple[object, ...] | None = None
        for item in values:
            if not isinstance(item, list) or len(item) != 6:
                raise ValueError(f"invalid ASE calibration observation in {path}")
            chromosome = canonical_chromosome(item[0])
            arm = clean(item[1])
            fraction = finite_float(item[2])
            weight = finite_float(item[3])
            contamination = finite_float(item[4])
            ambient = finite_float(item[5])
            if (not is_autosomal_chromosome(chromosome) or not arm
                    or not all(math.isfinite(value) for value in (
                        fraction, weight, contamination, ambient))
                    or not 0.0 <= fraction <= 1.0 or weight <= 0.0
                    or not 0.0 <= contamination < 1.0
                    or not 0.0 <= ambient <= 1.0):
                raise ValueError(f"invalid ASE calibration observation in {path}")
            key = (natural_key(chromosome), natural_key(arm))
            if previous_key is not None and key <= previous_key:
                raise ValueError(
                    f"ASE calibration observations are not uniquely sorted in {path}")
            previous_key = key
            result[orientation].append((
                chromosome, arm, fraction, weight, contamination, ambient))
    return result


class PackedAseCalibration:
    """Compact genome-wide ASE calibration observations for one model task."""

    def __init__(self) -> None:
        self.chromosome_ids = array("H")
        self.arm_ids = array("H")
        self.fractions = array("d")
        self.weights = array("d")
        self.contaminations = array("d")
        self.ambients = array("d")
        self.chromosomes: list[str] = []
        self.arms: list[str] = []
        self._chromosome_index: dict[str, int] = {}
        self._arm_index: dict[str, int] = {}
        self.index: dict[tuple[str, str], tuple[tuple[int, int], ...]] = {}

    @staticmethod
    def _intern(value: str, values: list[str], index: dict[str, int]) -> int:
        identifier = index.get(value)
        if identifier is None:
            identifier = len(values)
            if identifier >= 65535:
                raise ValueError("ASE calibration string dictionary overflow")
            index[value] = identifier
            values.append(value)
        return identifier

    def add(self, cell: tuple[str, str], payload: Mapping[
            str, Sequence[tuple[str, str, float, float, float, float]]]) -> None:
        if cell in self.index:
            raise ValueError(
                f"duplicate ASE calibration payload cell: {cell[0]}/{cell[1]}")
        ranges = []
        for orientation in ASE_ORIENTATIONS:
            start = len(self.fractions)
            for chromosome, arm, fraction, weight, contamination, ambient in (
                    payload.get(orientation, ())):
                self.chromosome_ids.append(self._intern(
                    chromosome, self.chromosomes, self._chromosome_index))
                self.arm_ids.append(self._intern(
                    arm, self.arms, self._arm_index))
                self.fractions.append(float(fraction))
                self.weights.append(float(weight))
                self.contaminations.append(float(contamination))
                self.ambients.append(float(ambient))
            ranges.append((start, len(self.fractions)))
        self.index[cell] = tuple(ranges)

    def observations(
            self, cell: tuple[str, str], orientation: str,
            ) -> Iterable[tuple[str, str, float, float, float, float]]:
        try:
            orientation_index = ASE_ORIENTATIONS.index(orientation)
            start, end = self.index[cell][orientation_index]
        except (KeyError, ValueError):
            return ()
        return (
            (self.chromosomes[self.chromosome_ids[index]],
             self.arms[self.arm_ids[index]], self.fractions[index],
             self.weights[index], self.contaminations[index],
             self.ambients[index])
            for index in range(start, end)
        )


def calibration_observations(
        calibration: object, cell: tuple[str, str], orientation: str,
        ) -> Iterable[tuple[str, str, float, float, float, float]]:
    if isinstance(calibration, PackedAseCalibration):
        return calibration.observations(cell, orientation)
    payload = calibration.get(cell, {})  # type: ignore[union-attr]
    return payload.get(orientation, ())

SCORE_FIELDS = [
    "library", "barcode", "calibration_group", "uid", "donor_a", "donor_b",
    "donor_pair", "chromosome", "arm", "hybrid_target_eligible",
    "ase_eligible", "expression_eligible", "ase_status", "expression_status",
    "ase_calibration_level", "expression_calibration_level",
    "ase_calibration_cells", "expression_calibration_cells",
    "ase_log_bf_DONOR_A_LOSS", "ase_log_bf_DONOR_B_LOSS",
    "ase_log_bf_DONOR_A_GAIN", "ase_log_bf_DONOR_B_GAIN",
    "ase_best_exact_state", "ase_direction", "ase_p_DONOR_A_DEPLETED",
    "ase_p_DONOR_A_ENRICHED", "ase_p_floor",
    "expression_component", "expression_best_copy_state",
    "expression_log_bf_LOSS", "expression_log_bf_GAIN",
    "expression_p_LOSS", "expression_p_GAIN", "expression_fold0_state",
    "expression_fold1_state", "expression_fold_replication_status",
    "expression_top1_fraction", "expression_top5_fraction",
    "expression_top10_fraction", "expression_nonzero_genes",
    "expression_p_floor", "expression_input_state",
    "expression_depth_bin", "expression_breadth_bin",
    "hybrid_conjunction_p_DONOR_A_LOSS",
    "hybrid_conjunction_p_DONOR_B_LOSS",
    "hybrid_conjunction_p_DONOR_A_GAIN",
    "hybrid_conjunction_p_DONOR_B_GAIN",
    "provisional_evidence_class", "provisional_resolved_state",
    "discordance_reason", "confounding_flags", "qc_flags", "schema_version",
]

COMPONENT_FIELDS = [
    "chromosome", "fold", "calibration_key", "calibration_level",
    "expression_input_state", "depth_bin", "breadth_bin", "component",
    "mean_shift_p", "mean_shift_q", "weight", "center_p", "center_q",
    "covariance_pp", "covariance_pq", "covariance_qq",
    "calibration_cells", "discovery_cells", "iterations", "status",
    "schema_version",
]

CALIBRATION_FIELDS = [
    "chromosome", "calibration_key", "branch", "fold", "calibration_level",
    "calibration_cells", "null_cells", "input_state", "depth_bin",
    "breadth_bin", "center_p", "center_q", "covariance_pp",
    "covariance_pq", "covariance_qq", "group_baseline_logit",
    "mapping_offset_ref", "mapping_offset_alt", "mapping_offset_mixed",
    "rho_ref", "rho_alt", "rho_mixed", "ref_alt_calibration_level",
    "mixed_calibration_level", "calibration_cells_ref",
    "calibration_cells_alt", "calibration_cells_mixed", "status",
    "schema_version",
]

MODEL_ROW_CORE_FIELDS = {
    "library", "barcode", "chromosome", "arm", "arm_side", "uid",
    "donor_a", "donor_b", "donor_pair", "calibration_group",
    "hybrid_target_eligible", "expression_reference_eligible",
    "ase_calibration_eligible", "cell_group_source",
    "cell_group_target_chromosome_excluded", "biological_block",
    "biological_block_source", "crossfit_fold", "ase_status",
    "expression_status", "ase_depth_bin", "ase_ambient_bin",
    "expression_depth_bin", "expression_breadth_bin", "confounding_flags",
    "qc_flags", "schema_version",
}
MODEL_ROW_ASE_FIELDS = {
    "ase_ambient_c", "ase_ambient_c_se", "ase_n_sites",
    "ase_effective_a_ref", "ase_effective_weight_ref", "ase_ambient_a_ref",
    "ase_effective_a_alt", "ase_effective_weight_alt", "ase_ambient_a_alt",
    "ase_effective_a_mixed", "ase_effective_weight_mixed",
    "ase_ambient_a_mixed", "ase_ambient_genotyped_mass",
    "ase_qname_fallback_fraction", "ase_evidence_basis",
    "ase_evidence_status", "ase_loco_arms",
}
MODEL_ROW_EXPRESSION_FIELDS = {
    "expression_model_p_arm", "expression_model_q_arm",
    "expression_model_p_score", "expression_model_q_score",
    "expression_model_p_fold0_score", "expression_model_p_fold1_score",
    "expression_model_q_fold0_score", "expression_model_q_fold1_score",
    "expression_model_p_nonzero_genes", "expression_model_q_nonzero_genes",
    "expression_model_p_top1_fraction", "expression_model_p_top5_fraction",
    "expression_model_p_top10_fraction", "expression_model_q_top1_fraction",
    "expression_model_q_top5_fraction", "expression_model_q_top10_fraction",
    "expression_model_expression_input_state",
    "expression_model_mapped_autosomal_library_size",
    "expression_model_ambient_c",
}
MODEL_ROW_FIELDS = (
    MODEL_ROW_CORE_FIELDS | MODEL_ROW_ASE_FIELDS | MODEL_ROW_EXPRESSION_FIELDS)


def require_absolute(path: str, label: str) -> str:
    if not os.path.isabs(path):
        raise ValueError(f"{label} must be absolute: {path}")
    return require_file(path, label)


def parse_target_libraries(value: object) -> set[str]:
    text = clean(value)
    if not text:
        return set()
    result: set[str] = set()
    for token in text.split(","):
        token = token.strip().lower().removeprefix("lib")
        if not token.isdigit() or int(token) < 1:
            raise ValueError(f"invalid --target-libraries token: {token}")
        normalized = str(int(token))
        if normalized in result:
            raise ValueError(f"duplicate --target-libraries token: {token}")
        result.add(normalized)
    return result


def load_shards(path: str, chromosome: str) -> tuple[
        list[dict[str, str]], list[str], PackedAseCalibration, set[str]]:
    manifest_path = require_file(path, "chromosome shard manifest")
    with open(manifest_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != SHARD_MANIFEST_HEADER:
            raise ValueError(
                f"shard manifest header must be {SHARD_MANIFEST_HEADER}: {path}")
        manifest_rows = list(reader)
    if not manifest_rows:
        raise ValueError(f"chromosome shard manifest has no libraries: {path}")
    rows: list[dict[str, str]] = []
    paths: list[str] = []
    seen_libraries: set[str] = set()
    seen_keys: set[tuple[str, str, str]] = set()
    calibration_by_cell = PackedAseCalibration()
    expected_crossfit_folds: int | None = None
    load_started = time.monotonic()
    for manifest_index, manifest_row in enumerate(manifest_rows, 1):
        library = clean(manifest_row.get("library")).lower().removeprefix("lib")
        if not library or library in seen_libraries:
            raise ValueError(f"empty/duplicate library in shard manifest: {path}")
        seen_libraries.add(library)
        shard = require_absolute(manifest_row["shard"], "hybrid shard")
        qc_path = require_absolute(manifest_row["qc"], "hybrid shard QC")
        contract_path = require_absolute(
            manifest_row["contract"], "hybrid shard contract")
        with open(contract_path, "r", encoding="utf-8") as handle:
            contract = json.load(handle)
        if (not isinstance(contract, dict)
                or contract.get("schema_version")
                   != "tetra_arm_hybrid_shard_contract_v1"
                or clean(contract.get("status")).upper() != "PASS"
                or clean(contract.get("library")).lower().removeprefix("lib")
                   != library
                or contract.get("output_schema") != HYBRID_SHARD_SCHEMA):
            raise ValueError(f"invalid hybrid shard contract: {contract_path}")
        outputs = contract.get("outputs", {})
        contract_shards = outputs.get("shards", {}) if isinstance(outputs, dict) else {}
        carrier_chromosome = canonical_chromosome(
            contract.get("ase_calibration_carrier_chromosome"))
        carrier_value = (str(contract_shards.get(carrier_chromosome, ""))
                         if isinstance(contract_shards, dict) else "")
        carrier = os.path.abspath(carrier_value) if carrier_value else ""
        if (not isinstance(contract_shards, dict)
                or os.path.abspath(str(contract_shards.get(chromosome, "")))
                   != os.path.abspath(shard)
                or not carrier_chromosome or not carrier
                or os.path.abspath(str(outputs.get("qc", "")))
                   != os.path.abspath(qc_path)
                or os.path.abspath(str(outputs.get("contract", "")))
                   != os.path.abspath(contract_path)):
            raise ValueError(
                f"hybrid shard contract output paths disagree: {contract_path}")
        folds = int(contract.get("crossfit_folds", 0))
        if folds < 2 or (expected_crossfit_folds is not None
                         and folds != expected_crossfit_folds):
            raise ValueError("hybrid shard crossfit-fold contracts disagree")
        expected_crossfit_folds = folds
        with open(qc_path, "r", encoding="utf-8", newline="") as handle:
            qc_reader = csv.DictReader(handle, delimiter="\t")
            if tuple(qc_reader.fieldnames or ()) != ("metric", "value"):
                raise ValueError(f"invalid hybrid shard QC header: {qc_path}")
            qc_rows = list(qc_reader)
        qc = {row.get("metric", ""): row.get("value", "") for row in qc_rows}
        if (len(qc) != len(qc_rows)
                or qc.get("schema_version") != "tetra_arm_hybrid_shard_qc_v1"
                or qc.get("status", "").upper() != "PASS"
                or clean(qc.get("library")).lower().removeprefix("lib")
                   != library):
            raise ValueError(f"invalid hybrid shard QC: {qc_path}")
        shard_header, _shard_count = validate_hybrid_table(
            shard, HYBRID_SHARD_SCHEMA)
        if not set(ASE_CALIBRATION_PAYLOAD_FIELDS).issubset(shard_header):
            raise ValueError(f"hybrid shard lacks ASE calibration payload: {shard}")
        if list(contract.get("fields", [])) != shard_header:
            raise ValueError(
                f"hybrid shard header disagrees with contract: {contract_path}")
        rows_before = len(rows)
        for row in read_tsv(shard):
            observed_chromosome = canonical_chromosome(row.get("chromosome"))
            if observed_chromosome != chromosome:
                raise ValueError(
                    f"wrong chromosome in shard {shard}: {observed_chromosome}")
            if clean(row.get("library")).lower().removeprefix("lib") != library:
                raise ValueError(f"wrong library in shard: {shard}")
            row = {field: clean(row.get(field)) for field in MODEL_ROW_FIELDS}
            row["library"], row["chromosome"] = library, chromosome
            key = library, clean(row.get("barcode")), clean(row.get("arm"))
            if not key[1] or not key[2] or key in seen_keys:
                raise ValueError(f"empty/duplicate score key {key}: {shard}")
            seen_keys.add(key)
            rows.append(row)
        chromosome_rows = contract.get("chromosome_rows", {})
        if (not isinstance(chromosome_rows, dict)
                or int(chromosome_rows.get(chromosome, -1))
                   != len(rows) - rows_before):
            raise ValueError(
                f"hybrid shard row count disagrees with contract: {shard}")
        paths.append(shard)

        sidecar_value = clean(outputs.get("ase_calibration"))
        if sidecar_value:
            sidecar = require_absolute(
                sidecar_value, "compact ASE calibration payload")
            with open_text(sidecar, "rt") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                expected_header = (
                    "library", "barcode", *ASE_CALIBRATION_PAYLOAD_FIELDS,
                    "schema_version")
                if tuple(reader.fieldnames or ()) != expected_header:
                    raise ValueError(
                        f"invalid compact ASE calibration header: {sidecar}")
                observed_payload_rows = 0
                for payload_row in reader:
                    observed_payload_rows += 1
                    barcode = clean(payload_row.get("barcode"))
                    if (clean(payload_row.get("library")).lower().removeprefix(
                            "lib") != library
                            or payload_row.get("schema_version")
                               != ASE_CALIBRATION_PAYLOAD_SCHEMA
                            or not barcode):
                        raise ValueError(
                            f"invalid compact ASE calibration row: {sidecar}")
                    payload = parse_ase_calibration_payload(payload_row, sidecar)
                    if not any(payload.values()):
                        raise ValueError(
                            f"empty compact ASE calibration row: {sidecar}")
                    calibration_by_cell.add((library, barcode), payload)
            if int(contract.get("ase_calibration_payload_rows", -1)) != (
                    observed_payload_rows):
                raise ValueError(
                    f"compact ASE calibration row count disagrees: {sidecar}")
            paths.append(sidecar)
            print(
                f"MODEL chromosome {chromosome}: loaded library "
                f"{manifest_index}/{len(manifest_rows)} ({library}), "
                f"{len(rows)} arm rows, "
                f"{len(calibration_by_cell.index)} ASE payload cells, "
                f"{time.monotonic() - load_started:.1f}s",
                flush=True)
            continue

        # Backward-compatible read path for v1.0.1 shards and hand-built
        # header-only fixtures. New shards always use the compact sidecar.
        carrier = require_absolute(carrier, "ASE calibration carrier shard")
        carrier_header, _carrier_count = validate_hybrid_table(
            carrier, HYBRID_SHARD_SCHEMA)
        if carrier_header != shard_header:
            raise ValueError(
                f"carrier/target hybrid shard headers disagree: {carrier}")
        carrier_rows = contract.get("chromosome_rows", {})
        observed_carrier_rows = 0
        for carrier_row in read_tsv(carrier):
            observed_carrier_rows += 1
            if (clean(carrier_row.get("library")).lower().removeprefix("lib")
                    != library
                    or canonical_chromosome(carrier_row.get("chromosome"))
                    != carrier_chromosome):
                raise ValueError(f"wrong library/chromosome in carrier: {carrier}")
            payload = parse_ase_calibration_payload(carrier_row, carrier)
            if not any(payload.values()):
                continue
            barcode = clean(carrier_row.get("barcode"))
            cell = (library, barcode)
            if not barcode or cell in calibration_by_cell.index:
                raise ValueError(
                    f"empty/duplicate ASE calibration carrier cell: {carrier}")
            calibration_by_cell.add(cell, payload)
        if (not isinstance(carrier_rows, dict)
                or int(carrier_rows.get(carrier_chromosome, -1))
                   != observed_carrier_rows):
            raise ValueError(
                f"ASE calibration carrier row count disagrees: {carrier}")
        if carrier not in paths:
            paths.append(carrier)
        print(
            f"MODEL chromosome {chromosome}: loaded library "
            f"{manifest_index}/{len(manifest_rows)} ({library}), "
            f"{len(rows)} arm rows, "
            f"{len(calibration_by_cell.index)} ASE payload cells, "
            f"{time.monotonic() - load_started:.1f}s",
            flush=True)
    return rows, paths, calibration_by_cell, seen_libraries


def cell_key(row: Mapping[str, object]) -> tuple[str, str]:
    return clean(row.get("library")), clean(row.get("barcode"))


def real_uid(row: Mapping[str, object]) -> str:
    return present_identifier(row.get("uid"))


def real_pair(row: Mapping[str, object]) -> str:
    return present_identifier(row.get("donor_pair"))


class CalibrationPoolIndex:
    """Index nuisance candidates by the four ordered fallback scopes."""

    def __init__(self, cells: Sequence[dict[str, str]]) -> None:
        self.source = cells
        self.all = tuple(sorted(cells, key=lambda row: (
            natural_key(row.get("library", "")),
            natural_key(row.get("barcode", "")))))
        by_library: dict[str, list[dict[str, str]]] = defaultdict(list)
        valid_by_group: dict[str, list[dict[str, str]]] = defaultdict(list)
        valid_by_library_group: dict[
            tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for row in self.all:
            library = clean(row.get("library"))
            group = clean(row.get("calibration_group"))
            by_library[library].append(row)
            if truthy(row.get("cell_group_target_chromosome_excluded")):
                valid_by_group[group].append(row)
                valid_by_library_group[(library, group)].append(row)
        self.by_library = {
            key: tuple(values) for key, values in by_library.items()}
        self.valid_by_group = {
            key: tuple(values) for key, values in valid_by_group.items()}
        self.valid_by_library_group = {
            key: tuple(values)
            for key, values in valid_by_library_group.items()}
        self._order = {id(row): index for index, row in enumerate(self.all)}
        self._expression_strata = None
        self._expression_candidate_cache: OrderedDict[
            tuple[object, ...], tuple[dict[str, str], ...]] = OrderedDict()
        self._expression_candidate_indices: dict[
            tuple[object, ...], np.ndarray] = {}
        self._expression_filter_columns = None
        self._expression_point_cache = None

    def expression_point_cache(self):
        """Parse each immutable cell's three fold-specific points only once."""
        if self._expression_point_cache is None:
            self._expression_point_cache = {
                id(row): (row_point(row), row_point(row, "0"),
                          row_point(row, "1"))
                for row in self.all
            }
        return self._expression_point_cache

    def _expression_filters(self):
        if self._expression_filter_columns is None:
            folds = [clean(row.get("crossfit_fold")) for row in self.all]
            pairs = [real_pair(row) for row in self.all]
            uids = [real_uid(row) for row in self.all]
            fold_codes = {value: index for index, value in enumerate(set(folds))}
            pair_codes = {value: index for index, value in enumerate(set(pairs))}
            uid_codes = {value: index for index, value in enumerate(set(uids))}
            self._expression_filter_columns = (
                np.asarray([fold_codes[value] for value in folds], dtype=np.int32),
                np.asarray([pair_codes[value] for value in pairs], dtype=np.int32),
                np.asarray([uid_codes[value] for value in uids], dtype=np.int32),
                np.asarray([
                    truthy(row.get("expression_reference_eligible"))
                    for row in self.all], dtype=np.bool_),
                fold_codes, pair_codes, uid_codes,
            )
        return self._expression_filter_columns

    def _expression_candidate_key(
            self, level: str, query: Mapping[str, object], fold: str):
        side = arm_side(query.get("arm")) or clean(query.get("arm_side"))
        stratum = expression_exact_stratum(query, side)
        if stratum is None:
            return None
        library = clean(query.get("library"))
        group = clean(query.get("calibration_group"))
        if level == "LIBRARY_VALID_GROUP_EXTERNAL_PAIR":
            scope = (library, group)
        elif level == "LIBRARY_EXTERNAL_PAIR":
            scope = (library,)
        elif level == "COHORT_VALID_GROUP_EXTERNAL_PAIR":
            scope = (group,)
        elif level == "COHORT_EXTERNAL_PAIR":
            scope = ()
        else:
            raise ValueError(f"unknown baseline level: {level}")
        return (level, scope, fold, side, *stratum)

    def filtered_expression_candidates(
            self, level: str, query: Mapping[str, object], fold: str,
            reference_only: bool) -> list[dict[str, str]]:
        """Apply the original external-block exclusions in stable pool order."""
        candidates = self.expression_candidates(level, query, fold)
        if not candidates:
            return []
        key = self._expression_candidate_key(level, query, fold)
        indices = self._expression_candidate_indices[key]
        (folds, pairs, uids, eligible, fold_codes, pair_codes,
         uid_codes) = self._expression_filters()
        query_fold = clean(query.get("crossfit_fold"))
        query_pair = real_pair(query)
        query_uid = real_uid(query)
        selected = folds[indices] != fold_codes.get(query_fold, -1)
        if query_pair:
            selected &= pairs[indices] != pair_codes.get(query_pair, -1)
        if query_uid:
            selected &= uids[indices] != uid_codes.get(query_uid, -1)
        if reference_only:
            selected &= eligible[indices]
        return [candidates[index] for index in np.flatnonzero(selected)]

    def expression_reference_pool(
            self, query: Mapping[str, object], fold: str,
            minimum: int) -> tuple[str, list[dict[str, str]]]:
        fallback_level = BASELINE_LEVELS[-1]
        fallback: list[dict[str, str]] = []
        for level in BASELINE_LEVELS:
            selected = self.filtered_expression_candidates(
                level, query, fold, reference_only=True)
            if len(selected) > len(fallback):
                fallback_level, fallback = level, selected
            if len(selected) >= minimum:
                return level, selected
        return fallback_level, fallback

    def candidates(self, level: str, query: Mapping[str, object]
                   ) -> Sequence[dict[str, str]]:
        library = clean(query.get("library"))
        group = clean(query.get("calibration_group"))
        group_valid = truthy(
            query.get("cell_group_target_chromosome_excluded"))
        if level == "LIBRARY_VALID_GROUP_EXTERNAL_PAIR":
            return (self.valid_by_library_group.get((library, group), ())
                    if group_valid else ())
        if level == "LIBRARY_EXTERNAL_PAIR":
            return self.by_library.get(library, ())
        if level == "COHORT_VALID_GROUP_EXTERNAL_PAIR":
            return self.valid_by_group.get(group, ()) if group_valid else ()
        if level == "COHORT_EXTERNAL_PAIR":
            return self.all
        raise ValueError(f"unknown baseline level: {level}")

    def expression_candidates(
            self, level: str, query: Mapping[str, object], fold: str
            ) -> Sequence[dict[str, str]]:
        """Return the original ordered pool restricted to matching bins.

        The nine neighboring depth/breadth bins are indexed once per task.
        Cells with only a sister-arm value for this fold remain in the pool,
        as in expression_stratum_match plus row_has_expression.
        """
        side = arm_side(query.get("arm")) or clean(query.get("arm_side"))
        stratum = expression_exact_stratum(query, side)
        if stratum is None:
            return []
        group_valid = truthy(
            query.get("cell_group_target_chromosome_excluded"))
        if level in (
                "LIBRARY_VALID_GROUP_EXTERNAL_PAIR",
                "COHORT_VALID_GROUP_EXTERNAL_PAIR") and not group_valid:
            return []
        if self._expression_strata is None:
            point_cache = self.expression_point_cache()
            indexes = {
                level_name: defaultdict(list)
                for level_name in BASELINE_LEVELS
            }
            for row in self.all:
                library = clean(row.get("library"))
                group = clean(row.get("calibration_group"))
                valid = truthy(
                    row.get("cell_group_target_chromosome_excluded"))
                for indexed_side in ("p", "q"):
                    exact = expression_exact_stratum(row, indexed_side)
                    if exact is None:
                        continue
                    for fold_index, indexed_fold in enumerate(
                            ("full", "0", "1")):
                        if not any(math.isfinite(value) for value in
                                   point_cache[id(row)][fold_index]):
                            continue
                        base = (indexed_fold, indexed_side, *exact)
                        indexes["COHORT_EXTERNAL_PAIR"][base].append(row)
                        indexes["LIBRARY_EXTERNAL_PAIR"][(library, *base)].append(row)
                        if valid:
                            indexes["COHORT_VALID_GROUP_EXTERNAL_PAIR"][(group, *base)].append(row)
                            indexes["LIBRARY_VALID_GROUP_EXTERNAL_PAIR"][(library, group, *base)].append(row)
            self._expression_strata = indexes
        library = clean(query.get("library"))
        group = clean(query.get("calibration_group"))
        scope: tuple[str, ...]
        if level == "LIBRARY_VALID_GROUP_EXTERNAL_PAIR":
            scope = (library, group)
        elif level == "LIBRARY_EXTERNAL_PAIR":
            scope = (library,)
        elif level == "COHORT_VALID_GROUP_EXTERNAL_PAIR":
            scope = (group,)
        elif level == "COHORT_EXTERNAL_PAIR":
            scope = ()
        else:
            raise ValueError(f"unknown baseline level: {level}")
        input_state, depth, breadth = stratum
        cache_key = (level, scope, fold, side, input_state, depth, breadth)
        cached = self._expression_candidate_cache.get(cache_key)
        if cached is not None:
            self._expression_candidate_cache.move_to_end(cache_key)
            return cached
        buckets = self._expression_strata[level]
        result = [row
                  for candidate_depth in range(depth - 1, depth + 2)
                  for candidate_breadth in range(breadth - 1, breadth + 2)
                  for row in buckets.get((
                      *scope, fold, side, input_state,
                      candidate_depth, candidate_breadth), ())]
        result.sort(key=lambda row: self._order[id(row)])
        ordered = tuple(result)
        self._expression_candidate_cache[cache_key] = ordered
        self._expression_candidate_indices[cache_key] = np.fromiter(
            (self._order[id(row)] for row in ordered),
            dtype=np.intp, count=len(ordered))
        if len(self._expression_candidate_cache) > 128:
            evicted, _ = self._expression_candidate_cache.popitem(last=False)
            self._expression_candidate_indices.pop(evicted)
        return ordered


def calibration_pool_index(
        cells: Sequence[dict[str, str]], args) -> CalibrationPoolIndex:
    index = getattr(args, "_hybrid_calibration_pool_index", None)
    if not isinstance(index, CalibrationPoolIndex) or index.source is not cells:
        index = CalibrationPoolIndex(cells)
        setattr(args, "_hybrid_calibration_pool_index", index)
    return index


class _IndexedNullPool:
    """One exact-stratum null pool with indexed biological exclusions."""

    def __init__(self, directions: Sequence[str]) -> None:
        self.directions = tuple(directions)
        self.values: dict[str, list[float] | array] = {
            direction: [] for direction in self.directions}
        self.by_pair: dict[object, dict[str, list[float] | array]] = {}
        self.by_uid: dict[object, dict[str, list[float] | array]] = {}
        self.by_pair_uid: dict[object, dict[str, list[float] | array]] = {}
        self.by_cell: dict[object, dict[str, list[float] | array]] = {}
        self.library_counts: Counter[str] = Counter()
        self.pair_library_counts: dict[object, Counter[str]] = {}
        self.uid_library_counts: dict[object, Counter[str]] = {}
        self.pair_uid_library_counts: dict[object, Counter[str]] = {}
        self.cell_library_counts: dict[object, Counter[str]] = {}
        self.finalized = False

    def _append_block(
            self, blocks: dict[object, dict[str, list[float] | array]],
            key: object, statistics: Mapping[str, float]) -> None:
        values = blocks.setdefault(
            key, {direction: [] for direction in self.directions})
        for direction in self.directions:
            values[direction].append(statistics[direction])  # type: ignore[union-attr]

    @staticmethod
    def _increment_library(
            blocks: dict[object, Counter[str]], key: object,
            library: str) -> None:
        blocks.setdefault(key, Counter())[library] += 1

    def add(self, row: Mapping[str, object],
            statistics: Mapping[str, float]) -> bool:
        if self.finalized:
            raise RuntimeError("cannot add to a finalized empirical-null pool")
        values = {
            direction: float(statistics.get(direction, math.nan))
            for direction in self.directions}
        if not all(math.isfinite(value) for value in values.values()):
            return False
        for direction in self.directions:
            self.values[direction].append(values[direction])  # type: ignore[union-attr]
        pair, uid, cell = real_pair(row), real_uid(row), cell_key(row)
        library = clean(row.get("library"))
        if pair:
            self._append_block(self.by_pair, pair, values)
            self._increment_library(self.pair_library_counts, pair, library)
        if uid:
            self._append_block(self.by_uid, uid, values)
            self._increment_library(self.uid_library_counts, uid, library)
        if pair and uid:
            pair_uid = (pair, uid)
            self._append_block(self.by_pair_uid, pair_uid, values)
            self._increment_library(
                self.pair_uid_library_counts, pair_uid, library)
        self._append_block(self.by_cell, cell, values)
        self._increment_library(self.cell_library_counts, cell, library)
        self.library_counts[library] += 1
        return True

    @staticmethod
    def _finalize_values(values: dict[str, list[float] | array]) -> None:
        for direction, observed in values.items():
            values[direction] = array("d", sorted(observed))

    @classmethod
    def _finalize_blocks(
            cls, blocks: dict[object, dict[str, list[float] | array]]) -> None:
        for values in blocks.values():
            cls._finalize_values(values)

    def finalize(self) -> None:
        if self.finalized:
            return
        self._finalize_values(self.values)
        for blocks in (
                self.by_pair, self.by_uid, self.by_pair_uid, self.by_cell):
            self._finalize_blocks(blocks)
        self.finalized = True

    def _exclusion_terms(self, query: Mapping[str, object]):
        pair, uid = real_pair(query), real_uid(query)
        if pair and uid:
            return (
                (-1, self.by_pair, self.pair_library_counts, pair),
                (-1, self.by_uid, self.uid_library_counts, uid),
                (1, self.by_pair_uid, self.pair_uid_library_counts,
                 (pair, uid)),
            )
        if pair:
            return ((-1, self.by_pair, self.pair_library_counts, pair),)
        if uid:
            return ((-1, self.by_uid, self.uid_library_counts, uid),)
        cell = cell_key(query)
        return ((-1, self.by_cell, self.cell_library_counts, cell),)

    def count(self, query: Mapping[str, object]) -> int:
        direction = self.directions[0]
        count = len(self.values[direction])
        for sign, blocks, _library_blocks, key in self._exclusion_terms(query):
            count += sign * len(blocks.get(key, {}).get(direction, ()))
        return max(0, count)

    @staticmethod
    def _tail_count(values: Sequence[float], statistic: float) -> int:
        return len(values) - bisect.bisect_left(values, statistic - 1e-12)

    def tail_count(self, statistic: float, direction: str,
                   query: Mapping[str, object]) -> int:
        count = self._tail_count(self.values[direction], statistic)
        for sign, blocks, _library_blocks, key in self._exclusion_terms(query):
            values = blocks.get(key, {}).get(direction, ())
            count += sign * self._tail_count(values, statistic)
        return max(0, count)

    def adjusted_library_counts(
            self, query: Mapping[str, object]) -> Counter[str]:
        counts = Counter(self.library_counts)
        for sign, _blocks, library_blocks, key in self._exclusion_terms(query):
            for library, count in library_blocks.get(key, {}).items():
                counts[library] += sign * count
        return counts


class IndexedEmpiricalNulls:
    """Calculate each pseudo-null once and query tails without materializing."""

    def __init__(self, directions: Sequence[str]) -> None:
        self.directions = tuple(directions)
        self.pools: dict[tuple[str, ...], _IndexedNullPool] = {}
        self.score_evaluations = 0
        self.retained_statistics = 0
        self.finalized = False

    def add(self, pool_key: Sequence[object], row: Mapping[str, object],
            statistics: Mapping[str, float]) -> bool:
        if self.finalized:
            raise RuntimeError("cannot add to a finalized empirical-null index")
        key = tuple(str(value) for value in pool_key)
        retained = self.pools.setdefault(
            key, _IndexedNullPool(self.directions)).add(row, statistics)
        if retained:
            self.retained_statistics += 1
        return retained

    def finalize(self) -> "IndexedEmpiricalNulls":
        if not self.finalized:
            for pool in self.pools.values():
                pool.finalize()
            self.finalized = True
        return self

    def _selected(self, pool_keys: Sequence[Sequence[object]]):
        if not self.finalized:
            raise RuntimeError("empirical-null index is not finalized")
        seen: set[tuple[str, ...]] = set()
        for pool_key in pool_keys:
            key = tuple(str(value) for value in pool_key)
            if key in seen:
                continue
            seen.add(key)
            pool = self.pools.get(key)
            if pool is not None:
                yield pool

    def count(self, pool_keys: Sequence[Sequence[object]],
              query: Mapping[str, object]) -> int:
        return sum(pool.count(query) for pool in self._selected(pool_keys))

    def pvalue(self, statistic: float, direction: str,
               pool_keys: Sequence[Sequence[object]],
               query: Mapping[str, object]) -> float:
        if not math.isfinite(statistic) or direction not in self.directions:
            return math.nan
        count = 0
        tail = 0
        for pool in self._selected(pool_keys):
            count += pool.count(query)
            tail += pool.tail_count(statistic, direction, query)
        return ((tail + 1.0) / (count + 1.0)) if count else math.nan

    def libraries(self, pool_keys: Sequence[Sequence[object]],
                  query: Mapping[str, object]) -> tuple[str, ...]:
        counts: Counter[str] = Counter()
        for pool in self._selected(pool_keys):
            counts.update(pool.adjusted_library_counts(query))
        return tuple(sorted(
            (library for library, count in counts.items() if count > 0),
            key=natural_key))


def row_point(row: Mapping[str, object], fold: str = "full") -> tuple[float, float]:
    prefix = "" if fold == "full" else f"fold{fold}_"
    return (
        finite_float(row.get(f"expression_model_p_{prefix}score")),
        finite_float(row.get(f"expression_model_q_{prefix}score")),
    )


def row_has_expression(row: Mapping[str, object], fold: str = "full") -> bool:
    return any(math.isfinite(value) for value in row_point(row, fold))


def row_has_side_expression(row: Mapping[str, object], side: str,
                            fold: str = "full") -> bool:
    if side not in {"p", "q"}:
        return False
    return math.isfinite(row_point(row, fold)[0 if side == "p" else 1])


def exclusion_match(query: Mapping[str, object],
                    candidate: Mapping[str, object]) -> bool:
    if clean(candidate.get("crossfit_fold")) == clean(query.get("crossfit_fold")):
        return False
    query_pair, candidate_pair = real_pair(query), real_pair(candidate)
    if query_pair and candidate_pair and query_pair == candidate_pair:
        return False
    query_uid, candidate_uid = real_uid(query), real_uid(candidate)
    if query_uid and candidate_uid and query_uid == candidate_uid:
        return False
    return True


def heldout_null_match(query: Mapping[str, object],
                       candidate: Mapping[str, object]) -> bool:
    """Select external biological blocks from the query's held-out fold."""
    if clean(candidate.get("crossfit_fold")) != clean(query.get("crossfit_fold")):
        return False
    query_pair, candidate_pair = real_pair(query), real_pair(candidate)
    if query_pair and candidate_pair and query_pair == candidate_pair:
        return False
    query_uid, candidate_uid = real_uid(query), real_uid(candidate)
    if query_uid and candidate_uid and query_uid == candidate_uid:
        return False
    return True


def self_crossfit_null_match(query: Mapping[str, object],
                             candidate: Mapping[str, object]) -> bool:
    """Select external blocks; the candidate supplies its own held-out fold."""
    query_cell, candidate_cell = cell_key(query), cell_key(candidate)
    if (all(query_cell) and all(candidate_cell)
            and query_cell == candidate_cell):
        return False
    query_pair, candidate_pair = real_pair(query), real_pair(candidate)
    if query_pair and candidate_pair and query_pair == candidate_pair:
        return False
    query_uid, candidate_uid = real_uid(query), real_uid(candidate)
    if query_uid and candidate_uid and query_uid == candidate_uid:
        return False
    return True


def level_match(level: str, query: Mapping[str, object],
                candidate: Mapping[str, object]) -> bool:
    same_library = clean(query.get("library")) == clean(candidate.get("library"))
    same_group = clean(query.get("calibration_group")) == clean(
        candidate.get("calibration_group"))
    group_valid = truthy(query.get("cell_group_target_chromosome_excluded"))
    candidate_group_valid = truthy(
        candidate.get("cell_group_target_chromosome_excluded"))
    if level == "LIBRARY_VALID_GROUP_EXTERNAL_PAIR":
        return same_library and same_group and group_valid and candidate_group_valid
    if level == "LIBRARY_EXTERNAL_PAIR":
        return same_library
    if level == "COHORT_VALID_GROUP_EXTERNAL_PAIR":
        return same_group and group_valid and candidate_group_valid
    if level == "COHORT_EXTERNAL_PAIR":
        return True
    raise ValueError(f"unknown baseline level: {level}")


def expression_exact_stratum(
        row: Mapping[str, object], side: str,
        ) -> tuple[str, int, int] | None:
    if side not in {"p", "q"}:
        return None
    try:
        depth = int(clean(row.get("expression_depth_bin")))
    except ValueError:
        return None
    breadth = int(math.floor(math.log2(max(0.0, finite_float(
        row.get(f"expression_model_{side}_nonzero_genes"), 0.0)) + 1.0)))
    return (
        clean(row.get("expression_model_expression_input_state")),
        depth, breadth)


def expression_stratum_match(query: Mapping[str, object],
                             candidate: Mapping[str, object]) -> bool:
    side = arm_side(query.get("arm")) or clean(query.get("arm_side"))
    query_stratum = expression_exact_stratum(query, side)
    candidate_stratum = expression_exact_stratum(candidate, side)
    return bool(
        query_stratum is not None and candidate_stratum is not None
        and query_stratum[0] == candidate_stratum[0]
        and abs(query_stratum[1] - candidate_stratum[1]) <= 1
        and abs(query_stratum[2] - candidate_stratum[2]) <= 1)


def expression_null_pool_keys(
        query: Mapping[str, object], fold: str, side: str,
        ) -> tuple[tuple[str, ...], ...]:
    stratum = expression_exact_stratum(query, side)
    if stratum is None:
        return ()
    input_state, depth, breadth = stratum
    return tuple(
        (fold, side, input_state, str(candidate_depth),
         str(candidate_breadth))
        for candidate_depth in range(depth - 1, depth + 2)
        for candidate_breadth in range(breadth - 1, breadth + 2))


def ase_null_pool_key(row: Mapping[str, object]) -> tuple[str, ...]:
    return (
        arm_side(row.get("arm")) or clean(row.get("arm_side")),
        clean(row.get("ase_depth_bin")), clean(row.get("ase_ambient_bin")))


def choose_pool(
        query: Mapping[str, object], candidates: Sequence[dict[str, str]],
        eligible: Callable[[Mapping[str, object]], bool], minimum: int,
        stratum: Callable[[Mapping[str, object], Mapping[str, object]], bool]
        | None = None,
        pool_index: CalibrationPoolIndex | None = None,
        candidate_provider: Callable[
            [str, Mapping[str, object]], Sequence[dict[str, str]]]
        | None = None,
        ordered_candidates: bool = False,
        ) -> tuple[str, list[dict[str, str]]]:
    fallback_level = BASELINE_LEVELS[-1]
    fallback: list[dict[str, str]] = []
    for level in BASELINE_LEVELS:
        level_candidates = (
            candidate_provider(level, query) if candidate_provider is not None
            else pool_index.candidates(level, query)
            if pool_index is not None else candidates)
        selected = [
            candidate for candidate in level_candidates
            if eligible(candidate) and exclusion_match(query, candidate)
            and level_match(level, query, candidate)
            and (stratum is None or stratum(query, candidate))
        ]
        if not ordered_candidates:
            selected.sort(key=lambda row: (
                natural_key(row.get("library", "")),
                natural_key(row.get("barcode", ""))))
        if len(selected) > len(fallback):
            fallback_level, fallback = level, selected
        if len(selected) >= minimum:
            return level, selected
    return fallback_level, fallback


def calibration_cache_key(query: Mapping[str, object], branch: str,
                          fold: str = "full") -> tuple[str, ...]:
    # With a donor pair, UID exclusion is a subset of the pair exclusion and
    # does not need a per-UID cache entry. If a malformed/exploratory query has
    # no pair, retain its UID so held-out-null membership remains exact.
    pair = real_pair(query)
    identity_exclusion = pair or (
        f"UID:{real_uid(query)}" if real_uid(query) else "NA")
    valid_group = truthy(query.get("cell_group_target_chromosome_excluded"))
    base = (
        branch, fold, clean(query.get("library")),
        clean(query.get("calibration_group")) if valid_group else "INVALID_GROUP",
        clean(query.get("crossfit_fold")), identity_exclusion,
    )
    if branch == "ASE":
        # Arm side is needed only because the matched empirical null is
        # side-specific. ASE depth and ambient bins preserve the native
        # matching principle without creating a per-cell fit. Active
        # orientation families keep an unused family from forcing a broader
        # nuisance fallback or being reused by a target that does need it. No
        # expression-derived field may partition ASE fits.
        active = ase_active_orientations(query)
        orientation_signature = "+".join((
            "REF_ALT" if "ref" in active or "alt" in active else "",
            "MIXED" if "mixed" in active else "",
        )).strip("+")
        # Direct calibration probes may omit target evidence fields. Preserve
        # their historical request for both orientation families; production
        # scoring rows always carry explicit effective weights.
        if not orientation_signature:
            orientation_signature = "REF_ALT+MIXED"
        return base + (
            clean(query.get("arm_side")) or "NA",
            clean(query.get("ase_depth_bin")) or "NA",
            clean(query.get("ase_ambient_bin")) or "NA",
            orientation_signature,
        )
    side = arm_side(query.get("arm")) or clean(query.get("arm_side"))
    breadth_bin = int(math.floor(math.log2(max(0.0, finite_float(
        query.get(f"expression_model_{side}_nonzero_genes"), 0.0)) + 1.0))) \
        if side in {"p", "q"} else -1
    return base + (
        clean(query.get("expression_model_expression_input_state")) or "NA",
        clean(query.get("expression_depth_bin")) or "NA",
        str(breadth_bin),
        clean(query.get("arm_side")) or "NA",
    )


def cache_key_text(key: Sequence[str]) -> str:
    return "|".join(key)


def best_side_row(rows: Sequence[dict[str, str]]) -> dict[str, str]:
    return sorted(rows, key=lambda row: (
        0 if arm_side(row.get("arm")) == "p" else 1,
        natural_key(row.get("arm", ""))))[0]


def query_for_side(row: Mapping[str, object], side: str) -> dict[str, str]:
    """Return a side-specific query without changing the shared cell record."""
    if side not in {"p", "q"}:
        raise ValueError(f"invalid arm side: {side}")
    result = {str(key): clean(value) for key, value in row.items()}
    model_arm = clean(row.get(f"expression_model_{side}_arm"))
    current_arm = clean(row.get("arm"))
    result["arm"] = model_arm or (
        current_arm if arm_side(current_arm) == side else
        f"{canonical_chromosome(row.get('chromosome'))}{side}")
    result["arm_side"] = side
    return result


def make_cell_records(rows: Sequence[dict[str, str]]) -> tuple[
        list[dict[str, str]], dict[tuple[str, str], list[dict[str, str]]]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[cell_key(row)].append(row)
    cells: list[dict[str, str]] = []
    for key in sorted(grouped, key=lambda value: tuple(natural_key(part) for part in value)):
        values = grouped[key]
        representative = best_side_row(values)
        fields = (
            "expression_model_p_score", "expression_model_q_score",
            "expression_model_p_fold0_score", "expression_model_p_fold1_score",
            "expression_model_q_fold0_score", "expression_model_q_fold1_score",
            "expression_model_expression_input_state",
            "expression_model_mapped_autosomal_library_size",
            "expression_depth_bin",
            "hybrid_target_eligible", "expression_reference_eligible",
            "ase_calibration_eligible", "crossfit_fold", "donor_pair", "uid",
            "biological_block", "biological_block_source",
            "calibration_group", "cell_group_target_chromosome_excluded",
        )
        for field in fields:
            observed = {clean(row.get(field)) for row in values}
            if len(observed) > 1:
                raise ValueError(f"cell-level field differs across arms for {key}: {field}")
        cells.append(representative)
    return cells, grouped


def expression_side_metrics(row: Mapping[str, object]) -> tuple[int, float, float, float]:
    side = arm_side(row.get("arm"))
    if not side:
        return 0, math.nan, math.nan, math.nan
    return (
        int(finite_float(row.get(f"expression_model_{side}_nonzero_genes"), 0.0)),
        finite_float(row.get(f"expression_model_{side}_top1_fraction")),
        finite_float(row.get(f"expression_model_{side}_top5_fraction")),
        finite_float(row.get(f"expression_model_{side}_top10_fraction")),
    )


def expression_breadth_pass(row: Mapping[str, object], args) -> bool:
    nonzero, top1, top5, top10 = expression_side_metrics(row)
    return (
        nonzero >= args.min_expression_nonzero_genes
        and math.isfinite(top1) and top1 <= args.max_expression_top1_fraction
        and math.isfinite(top5) and top5 <= args.max_expression_top5_fraction
        and math.isfinite(top10) and top10 <= args.max_expression_top10_fraction
    )


def ase_evidence(row: Mapping[str, object]) -> dict[str, object]:
    result = {
        "ambient_c": row.get("ase_ambient_c", "NA"),
        "ambient_c_se": row.get("ase_ambient_c_se", "NA"),
        "evidence_status": row.get("ase_evidence_status", "NO_DATA"),
    }
    for orientation in ("ref", "alt", "mixed"):
        result[f"effective_a_{orientation}"] = row.get(
            f"ase_effective_a_{orientation}", "0")
        result[f"effective_weight_{orientation}"] = row.get(
            f"ase_effective_weight_{orientation}", "0")
        result[f"ambient_a_{orientation}"] = row.get(
            f"ase_ambient_a_{orientation}", "0.5")
    return result


def total_ase_weight(row: Mapping[str, object]) -> float:
    return sum(max(0.0, finite_float(
        row.get(f"ase_effective_weight_{orientation}"), 0.0))
        for orientation in ("ref", "alt", "mixed"))


def ase_contamination(row: Mapping[str, object]) -> float:
    value = finite_float(row.get("ase_ambient_c"))
    if not math.isfinite(value):
        value = finite_float(row.get("expression_model_ambient_c"), 0.0)
    return bounded(value, 0.0, 0.999999)


def ase_scoring_safety_status(row: Mapping[str, object], args) -> str:
    """Defensively enforce native production ASE safety gates."""
    shard_status = clean(row.get("ase_status")) or "NO_DATA"
    evidence_status = clean(row.get("ase_evidence_status")) or "NO_DATA"
    evidence_basis = clean(row.get("ase_evidence_basis")).upper()
    if (shard_status == "UNSAFE_SITE_FALLBACK"
            or evidence_status == "PASS_SITE_FALLBACK"
            or "SITE" in evidence_basis):
        return "UNSAFE_SITE_FALLBACK"
    if shard_status != "PASS" or evidence_status != "PASS":
        return shard_status if shard_status != "PASS" else evidence_status
    contamination = finite_float(row.get("ase_ambient_c"))
    if math.isfinite(contamination) and contamination >= 1.0:
        return "UNSAFE_AMBIENT_C_AT_UPPER_BOUNDARY"
    if (finite_float(row.get("ase_qname_fallback_fraction"), 1.0)
            > args.max_ase_qname_fallback_fraction):
        return "UNSAFE_HIGH_QNAME_FALLBACK"
    if (finite_float(row.get("ase_ambient_genotyped_mass"), 0.0)
            < args.min_ase_ambient_genotyped_mass):
        return "UNSAFE_LOW_AMBIENT_GENOTYPED_MASS"
    standard_error = finite_float(row.get("ase_ambient_c_se"))
    if not math.isfinite(standard_error) or standard_error <= 0.0:
        return "UNSAFE_AMBIENT_UNCERTAINTY_UNAVAILABLE"
    return "PASS"


def fit_ase_nuisance(
        cells: Sequence[dict[str, str]], args,
        calibration_by_cell: Mapping[
            tuple[str, str],
            Mapping[str, Sequence[
                tuple[str, str, float, float, float, float]]]],
        excluded_chromosome: str) -> dict[str, object]:
    """Fit native-compatible cell-pooled bias and cell-arm rho."""
    observations: dict[str, list[tuple[float, float, float]]] = {
        orientation: [] for orientation in ("ref", "alt", "mixed")}
    rho_observations: dict[str, list[tuple[float, float, float]]] = {
        orientation: [] for orientation in ("ref", "alt", "mixed")}
    n_cells = {orientation: 0 for orientation in ("ref", "alt", "mixed")}
    for row in cells:
        for orientation in ("ref", "alt", "mixed"):
            retained = [
                (fraction, weight, bounded(
                    (1.0 - contamination) * 0.5 + contamination * ambient,
                    1e-8, 1.0 - 1e-8))
                for (source_chromosome, _arm, fraction, weight,
                     contamination, ambient)
                in calibration_observations(
                    calibration_by_cell, cell_key(row), orientation)
                if source_chromosome != excluded_chromosome
            ]
            if not retained:
                continue
            total_weight = sum(item[1] for item in retained)
            pooled_fraction = sum(
                fraction * weight for fraction, weight, _base in retained
            ) / total_weight
            pooled_base = sum(
                base * weight for _fraction, weight, base in retained
            ) / total_weight
            observations[orientation].append((
                bounded(pooled_fraction, 1e-8, 1.0 - 1e-8),
                total_weight, bounded(pooled_base, 1e-8, 1.0 - 1e-8)))
            rho_observations[orientation].extend(retained)
            n_cells[orientation] += 1

    fits = {
        orientation: fit_ase_bias_component(
            values, rho_observations=rho_observations[orientation],
            n_cells=n_cells[orientation],
            maximum_weight=args.ase_calibration_max_weight,
            huber_z=args.ase_calibration_huber_z,
            minimum_robust_weight=args.ase_calibration_min_robust_weight,
            default_rho=args.default_rho, minimum_rho=args.min_rho,
            maximum_rho=args.max_rho)
        for orientation, values in observations.items()
    }
    ref, alt, mixed = fits["ref"], fits["alt"], fits["mixed"]
    adequate = {
        orientation: (
            fit.observations > 0
            and fit.total_weight
            >= args.min_ase_orientation_calibration_effective_weight)
        for orientation, fit in fits.items()
    }
    candidates: list[tuple[float, float]] = []
    matched_ref_alt = adequate["ref"] and adequate["alt"]
    if matched_ref_alt:
        candidates.append(((ref.delta + alt.delta) / 2.0,
                           min(ref.total_weight, alt.total_weight)))
    if adequate["mixed"]:
        candidates.append((mixed.delta, mixed.total_weight))
    total = sum(weight for _value, weight in candidates)
    group_baseline = (
        sum(value * weight for value, weight in candidates) / total
        if total > 0.0 else 0.0)
    mapping = (ref.delta - alt.delta) / 2.0 if matched_ref_alt else 0.0
    offsets = {"ref": mapping, "alt": -mapping, "mixed": 0.0}
    rhos = {orientation: fit.rho for orientation, fit in fits.items()}
    orientation_cells = {
        orientation: fit.observations for orientation, fit in fits.items()}
    status = (
        "FULL_ORIENTATION_CALIBRATION"
        if matched_ref_alt and adequate["mixed"]
        else "PARTIAL_ORIENTATION_CALIBRATION"
        if matched_ref_alt or adequate["mixed"]
        else "NEUTRAL_PRIOR_FALLBACK")
    return {
        "group_baseline": group_baseline, "offsets": offsets, "rhos": rhos,
        "cells": max(orientation_cells.values(), default=0), "status": status,
        "fits": fits, "adequate": adequate,
        "orientation_cells": orientation_cells,
        "matched_ref_alt": matched_ref_alt,
    }


def ase_active_orientations(row: Mapping[str, object]) -> tuple[str, ...]:
    return tuple(
        orientation for orientation in ("ref", "alt", "mixed")
        if finite_float(row.get(
            f"ase_effective_weight_{orientation}"), 0.0) > 0.0)


def ase_orientation_calibration_status(
        row: Mapping[str, object], nuisance: Mapping[str, object], args,
        ) -> tuple[str, int]:
    """Require calibration only for orientations used by this observation."""
    active = ase_active_orientations(row)
    if not active:
        return "NO_ACTIVE_ORIENTATION", 0
    adequate = nuisance.get("adequate", {})
    counts = nuisance.get("orientation_cells", {})
    required_counts: list[int] = []
    missing = False
    if "ref" in active or "alt" in active:
        matched = bool(nuisance.get("matched_ref_alt"))
        missing = missing or not matched
        required_counts.append(min(
            int(counts.get("ref", 0)), int(counts.get("alt", 0))))
    if "mixed" in active:
        missing = missing or not bool(adequate.get("mixed"))
        required_counts.append(int(counts.get("mixed", 0)))
    cells = min(required_counts, default=0)
    if missing:
        return "PARTIAL_ORIENTATION_CALIBRATION", cells
    if cells < args.min_ase_calibration_cells:
        return "WEAK_ORIENTATION_CALIBRATION", cells
    return "PASS_EMPIRICAL_CALIBRATION", cells


def cell_ase_baseline(
        row: Mapping[str, object], nuisance: Mapping[str, object], args,
        calibration_by_cell: Mapping[
            tuple[str, str],
            Mapping[str, Sequence[
                tuple[str, str, float, float, float, float]]]],
        excluded_chromosome: str) -> tuple[float, str]:
    cache = getattr(args, "_hybrid_ase_baseline_cache", None)
    if cache is None:
        cache = {}
        setattr(args, "_hybrid_ase_baseline_cache", cache)
    cache_key = (cell_key(row), excluded_chromosome, id(nuisance))
    cached = cache.get(cache_key)
    if cached is not None and cached[0] is nuisance:
        model_timings(args)["ase_baseline_hits"] += 1
        return cached[1]
    baseline_started = time.monotonic()
    offsets = nuisance["offsets"]
    observations: list[tuple[str, str, float, float, float, float]] = []
    for orientation in ("ref", "alt", "mixed"):
        for (source_chromosome, arm, observed, weight, contamination,
             ambient) in calibration_observations(
                 calibration_by_cell, cell_key(row), orientation):
            if source_chromosome != excluded_chromosome:
                observations.append((
                    orientation, arm, observed, weight, contamination, ambient))
    group = float(nuisance["group_baseline"])
    if not observations:
        result = group, "GROUP_BASELINE_ONLY"
        cache[cache_key] = (nuisance, result)
        model_timings(args)["ase_baseline_seconds"] += (
            time.monotonic() - baseline_started)
        model_timings(args)["ase_baseline_fits"] += 1
        return result

    def fit_eta(robust_weights: Sequence[float]) -> float:
        fitted = [
            (orientation, observed,
             min(weight, args.ase_calibration_max_weight) * robust_weights[index],
             contamination, ambient)
            for index, (orientation, _arm, observed, weight, contamination,
                        ambient) in enumerate(observations)
        ]

        def score_and_slope(eta: float) -> tuple[float, float]:
            value = 0.0
            slope = 0.0
            cellular = 1.0 / (1.0 + math.exp(-eta)) if eta >= 0.0 else (
                math.exp(eta) / (1.0 + math.exp(eta)))
            for orientation, observed, weight, contamination, ambient in fitted:
                expected = ambient_adjusted_fraction(
                    "BALANCED", orientation, eta, contamination, ambient,
                    float(offsets[orientation]))
                value += weight * (observed - expected)
                mixed = ((1.0 - contamination) * cellular
                         + contamination * ambient)
                if 1e-12 < mixed < 1.0 - 1e-12:
                    slope -= (weight * expected * (1.0 - expected)
                              / (mixed * (1.0 - mixed))
                              * (1.0 - contamination) * cellular
                              * (1.0 - cellular))
            return value, slope

        lower, upper = -6.0, 6.0
        if score_and_slope(lower)[0] <= 0.0:
            return lower
        if score_and_slope(upper)[0] >= 0.0:
            return upper
        eta = bounded(group, lower, upper)
        for _ in range(64):
            value, slope = score_and_slope(eta)
            if value > 0.0:
                lower = eta
            else:
                upper = eta
            if slope < 0.0 and abs(value / slope) <= 1e-13:
                return eta
            if upper - lower <= 2e-14:
                break
            next_eta = eta - value / slope if slope < 0.0 else math.nan
            if not lower < next_eta < upper:
                next_eta = (lower + upper) / 2.0
            eta = next_eta
        return (lower + upper) / 2.0

    robust = [1.0] * len(observations)
    raw = group
    for _ in range(4):
        raw = fit_eta(robust)
        updated = []
        for (orientation, _arm, observed, weight, contamination,
             ambient) in observations:
            expected = ambient_adjusted_fraction(
                "BALANCED", orientation, raw, contamination, ambient,
                float(offsets[orientation]))
            rho = float(nuisance["rhos"][orientation])
            variance = max(expected * (1.0 - expected) * (
                rho + (1.0 - rho) / max(weight, 1e-8)), 1e-6)
            z_score = abs(observed - expected) / math.sqrt(variance)
            updated.append(min(
                1.0, args.ase_calibration_huber_z / max(z_score, 1e-12)))
        robust = updated
    retained = [weight >= args.ase_calibration_min_robust_weight
                for weight in robust]
    if sum(retained) >= 3:
        raw = fit_eta([1.0 if keep else 0.0 for keep in retained])
    retained_arms = len({
        arm for keep, (_orientation, arm, _observed, _weight,
                       _contamination, _ambient)
        in zip(retained, observations) if keep
    })
    shrinkage = retained_arms / (
        retained_arms + args.ase_cell_baseline_shrinkage_arms)
    baseline = shrinkage * raw + (1.0 - shrinkage) * group
    status = ("PASS_CELL_LOCO_BASELINE"
              if retained_arms >= args.min_ase_cell_baseline_arms
              else "SPARSE_CELL_LOCO_SHRUNK_TO_GROUP")
    result = baseline, status
    cache[cache_key] = (nuisance, result)
    model_timings(args)["ase_baseline_seconds"] += (
        time.monotonic() - baseline_started)
    model_timings(args)["ase_baseline_fits"] += 1
    return result


def expression_nuisance_for_query(
        query: dict[str, str], cells: Sequence[dict[str, str]], fold: str,
        args, cache: dict[tuple[str, ...], dict[str, object]],
        component_rows: dict[tuple[str, ...], dict[str, object]],
        chromosome: str) -> dict[str, object]:
    """Fit expression nuisance parameters with the query block excluded."""
    key = calibration_cache_key(query, "EXPRESSION", fold)
    if key in cache:
        return cache[key]
    pool_index = calibration_pool_index(cells, args)
    timings = model_timings(args)
    fit_started = time.monotonic()
    level, references = pool_index.expression_reference_pool(
        query, fold, args.min_expression_reference_cells)
    timings["expression_pool_seconds"] += time.monotonic() - fit_started
    fit_started = time.monotonic()
    point_cache = pool_index.expression_point_cache()
    fold_index = {"full": 0, "0": 1, "1": 2}[fold]
    fit = robust_bivariate_fit(
        [point_cache[id(row)][fold_index] for row in references],
        args.min_expression_sigma, args.max_expression_sigma,
        args.max_expression_abs_correlation, args.expression_df)
    timings["expression_robust_seconds"] += time.monotonic() - fit_started
    fit_started = time.monotonic()
    discoveries = pool_index.filtered_expression_candidates(
        level, query, fold, reference_only=False)
    timings["expression_discovery_seconds"] += time.monotonic() - fit_started
    fit_started = time.monotonic()
    mixture = fit_fixed_center_mixture(
        [point_cache[id(row)][fold_index] for row in discoveries], fit.center,
        fit.covariance, args.expression_df, args.expression_outlier_scale)
    timings["expression_mixture_seconds"] += time.monotonic() - fit_started
    timings["expression_nuisance_fits"] += 1
    timings["expression_reference_cells"] += len(references)
    timings["expression_discovery_cells"] += len(discoveries)
    calibration_key = cache_key_text(key)
    status = "PASS" if (
        len(references) >= args.min_expression_reference_cells
        and fit.status == "PASS" and mixture.status == "PASS") else (
            "WEAK_COMPONENT_LOW_SUPPORT" if references
            else "NO_EXTERNAL_REFERENCE")
    result = {
        "key": calibration_key, "level": level, "references": len(references),
        "fit": fit, "mixture": mixture, "status": status,
        "nuisance_status": status,
    }
    cache[key] = result
    for component in EXPRESSION_PATTERNS:
        shift = EXPRESSION_PATTERN_SHIFTS[component]
        component_rows[(calibration_key, fold, component)] = {
            "chromosome": chromosome, "fold": fold,
            "calibration_key": calibration_key, "calibration_level": level,
            "expression_input_state": clean(
                query.get("expression_model_expression_input_state")) or "NA",
            "depth_bin": clean(query.get("expression_depth_bin")) or "NA",
            "breadth_bin": key[-2],
            "component": component, "mean_shift_p": f"{shift[0]:.17g}",
            "mean_shift_q": f"{shift[1]:.17g}",
            "weight": f"{float(mixture.weights[component]):.17g}",
            "center_p": f"{fit.center[0]:.17g}",
            "center_q": f"{fit.center[1]:.17g}",
            "covariance_pp": f"{fit.covariance[0][0]:.17g}",
            "covariance_pq": f"{fit.covariance[0][1]:.17g}",
            "covariance_qq": f"{fit.covariance[1][1]:.17g}",
            "calibration_cells": len(references),
            "discovery_cells": len(discoveries), "iterations": mixture.iterations,
            "status": status, "schema_version": HYBRID_COMPONENT_SCHEMA,
        }
    return result


def build_expression_null_index(
        cells: Sequence[dict[str, str]], args,
        cache: dict[tuple[str, ...], dict[str, object]],
        component_rows: dict[tuple[str, ...], dict[str, object]],
        chromosome: str) -> IndexedEmpiricalNulls:
    """Score each expression pseudo-null once per fold, side, and stratum."""
    index = IndexedEmpiricalNulls(COPY_STATES)
    references = sorted(
        (row for row in cells
         if truthy(row.get("expression_reference_eligible"))),
        key=lambda row: (
            natural_key(row.get("library", "")),
            natural_key(row.get("barcode", ""))))
    started = next_progress = time.monotonic()
    for reference_index, reference in enumerate(references, 1):
        if time.monotonic() >= next_progress + 240.0:
            timings = model_timings(args)
            print(
                f"MODEL chromosome {chromosome}: expression null progress "
                f"{reference_index}/{len(references)} cells, "
                f"{index.score_evaluations} evaluations, "
                f"{len(cache)} fits, {time.monotonic() - started:.1f}s; "
                f"pool={timings['expression_pool_seconds']:.1f}s "
                f"discovery={timings['expression_discovery_seconds']:.1f}s "
                f"mixture={timings['expression_mixture_seconds']:.1f}s "
                f"score={timings['expression_null_score_seconds']:.1f}s",
                flush=True)
            next_progress = time.monotonic()
        for fold in ("full", "0", "1"):
            if not row_has_expression(reference, fold):
                continue
            point = row_point(reference, fold)
            for side in ("p", "q"):
                stratum = expression_exact_stratum(reference, side)
                if (stratum is None
                        or not row_has_side_expression(reference, side, fold)):
                    continue
                index.score_evaluations += 1
                reference_nuisance = expression_nuisance_for_query(
                    query_for_side(reference, side), cells, fold, args, cache,
                    component_rows, chromosome)
                if reference_nuisance["nuisance_status"] != "PASS":
                    continue
                score_started = time.monotonic()
                scores = expression_direction_log_bfs(
                    point, side, reference_nuisance["fit"].center,
                    reference_nuisance["fit"].covariance, args.expression_df,
                    reference_nuisance["mixture"].weights,
                    args.expression_outlier_scale)
                model_timings(args)["expression_null_score_seconds"] += (
                    time.monotonic() - score_started)
                index.add(
                    (fold, side, stratum[0], stratum[1], stratum[2]),
                    reference,
                    {direction: max(0.0, scores[direction])
                     for direction in COPY_STATES})
    print(
        f"MODEL chromosome {chromosome}: finalizing expression null index "
        f"({index.retained_statistics} statistics, {len(index.pools)} pools)",
        flush=True)
    finalize_started = time.monotonic()
    result = index.finalize()
    model_timings(args)["expression_finalize_seconds"] += (
        time.monotonic() - finalize_started)
    return result


def expression_model_for_query(
        query: dict[str, str], cells: Sequence[dict[str, str]], fold: str,
        args, null_index: IndexedEmpiricalNulls,
        cache: dict[tuple[str, ...], dict[str, object]],
        component_rows: dict[tuple[str, ...], dict[str, object]],
        calibration_rows: dict[tuple[str, ...], dict[str, object]],
        chromosome: str) -> dict[str, object]:
    """Combine a local nuisance fit with the calculate-once cohort null."""
    key = calibration_cache_key(query, "EXPRESSION", fold)
    model_key = ("EXPRESSION_MODEL",) + key
    if model_key in cache:
        return cache[model_key]
    nuisance = expression_nuisance_for_query(
        query, cells, fold, args, cache, component_rows, chromosome)
    query_side = arm_side(query.get("arm"))
    level = str(nuisance["level"])
    pool_keys = expression_null_pool_keys(query, fold, query_side)
    matched_null_count = null_index.count(pool_keys, query)
    null_libraries = null_index.libraries(pool_keys, query)
    calibration_key = str(nuisance["key"])
    status = "PASS" if (
        nuisance["nuisance_status"] == "PASS"
        and matched_null_count >= args.min_expression_reference_cells) else (
            "WEAK_COMPONENT_LOW_SUPPORT"
            if nuisance["references"] or matched_null_count
            else "NO_EXTERNAL_REFERENCE")
    result = {
        **nuisance, "null_references": matched_null_count,
        "null_count": matched_null_count, "null_pool_keys": pool_keys,
        "null_libraries": null_libraries, "status": status,
    }
    cache[model_key] = result
    calibration_rows[("EXPRESSION", calibration_key, fold)] = {
        "chromosome": chromosome, "calibration_key": calibration_key,
        "branch": "EXPRESSION", "fold": fold, "calibration_level": level,
        "calibration_cells": result["references"],
        "null_cells": matched_null_count,
        "input_state": clean(
            query.get("expression_model_expression_input_state")) or "NA",
        "depth_bin": clean(query.get("expression_depth_bin")) or "NA",
        "breadth_bin": key[-2],
        "center_p": f"{result['fit'].center[0]:.17g}",
        "center_q": f"{result['fit'].center[1]:.17g}",
        "covariance_pp": f"{result['fit'].covariance[0][0]:.17g}",
        "covariance_pq": f"{result['fit'].covariance[0][1]:.17g}",
        "covariance_qq": f"{result['fit'].covariance[1][1]:.17g}",
        "group_baseline_logit": "NA", "mapping_offset_ref": "NA",
        "mapping_offset_alt": "NA", "mapping_offset_mixed": "NA",
        "rho_ref": "NA", "rho_alt": "NA", "rho_mixed": "NA",
        "ref_alt_calibration_level": "NA", "mixed_calibration_level": "NA",
        "calibration_cells_ref": "NA", "calibration_cells_alt": "NA",
        "calibration_cells_mixed": "NA",
        "status": status, "schema_version": HYBRID_CALIBRATION_SCHEMA,
    }
    return result


def ase_nuisance_for_query(
        query: dict[str, str], cells: Sequence[dict[str, str]], args,
        calibration_by_cell: Mapping[
            tuple[str, str], Mapping[str, Sequence[
                tuple[str, str, float, float, float, float]]]],
        cache: dict[tuple[str, ...], dict[str, object]],
        chromosome: str) -> dict[str, object]:
    """Fit ASE nuisances with independent fallback for ref/alt and mixed."""
    full_key = calibration_cache_key(query, "ASE")
    # Side/depth/ambient define the matched null, not the nuisance fit.  Keeping
    # them out of this cache key avoids refitting the same exclusion set for
    # every arm-depth stratum in a chromosome.
    orientation_signature = full_key[-1]
    key = ("ASE_NUISANCE",) + full_key[1:6] + (orientation_signature,)
    if key in cache:
        return cache[key]
    pool_index = calibration_pool_index(cells, args)
    # The selected cells and all three orientation fits depend on the query
    # exclusion block and hierarchy level, not on the target's active
    # orientation signature.  Share those expensive fits across signatures.
    level_cache = getattr(args, "_hybrid_ase_level_fit_cache", None)
    source = getattr(args, "_hybrid_ase_level_fit_source", None)
    if (level_cache is None or source is None
            or source[0] is not cells
            or source[1] is not calibration_by_cell
            or source[2] != chromosome):
        level_cache = {}
        setattr(args, "_hybrid_ase_level_fit_cache", level_cache)
        setattr(args, "_hybrid_ase_level_fit_source", (
            cells, calibration_by_cell, chromosome))
    exclusion_key = full_key[1:6]
    require_ref_alt = "REF_ALT" in orientation_signature
    require_mixed = "MIXED" in orientation_signature
    level_fits: list[tuple[str, int, dict[str, object]]] = []
    ref_entry = None
    mixed_entry = None
    for level in BASELINE_LEVELS:
        level_key = (exclusion_key, level)
        entry = level_cache.get(level_key)
        if entry is None:
            fit_started = time.monotonic()
            selected = [
                row for row in pool_index.candidates(level, query)
                if truthy(row.get("ase_calibration_eligible"))
                and int(finite_float(row.get("ase_loco_arms"), 0.0)) > 0
                and exclusion_match(query, row) and level_match(level, query, row)
            ]
            selected.sort(key=lambda row: (
                natural_key(row.get("library", "")),
                natural_key(row.get("barcode", ""))))
            entry = (level, len(selected), fit_ase_nuisance(
                selected, args, calibration_by_cell, chromosome))
            level_cache[level_key] = entry
            timings = model_timings(args)
            timings["ase_level_fit_seconds"] += time.monotonic() - fit_started
            timings["ase_level_fits"] += 1
            timings["ase_selected_cells"] += len(selected)
        level_fits.append(entry)

        # The hierarchy is ordered from the narrowest local nuisance pool to
        # the broadest cohort fallback.  Once both independently calibrated
        # orientation families pass at a level, broader fits cannot be
        # selected: the code below always chooses the first passing entry.
        # Avoiding those mathematically irrelevant fits is critical because a
        # cohort fallback traverses every loaded cell and every retained
        # cell-arm observation.
        nuisance = entry[2]
        counts = nuisance["orientation_cells"]
        if (ref_entry is None and nuisance["matched_ref_alt"]
                and min(counts["ref"], counts["alt"])
                >= args.min_ase_calibration_cells):
            ref_entry = entry
        if (mixed_entry is None and nuisance["adequate"]["mixed"]
                and counts["mixed"] >= args.min_ase_calibration_cells):
            mixed_entry = entry
        if ((ref_entry is not None or not require_ref_alt)
                and (mixed_entry is not None or not require_mixed)):
            break

    def ref_alt_rank(entry) -> tuple[float, ...]:
        nuisance = entry[2]
        fits = nuisance["fits"]
        return (
            float(nuisance["matched_ref_alt"]),
            min(fits["ref"].total_weight, fits["alt"].total_weight),
            min(fits["ref"].observations, fits["alt"].observations),
            -float(BASELINE_LEVELS.index(entry[0])),
        )

    def mixed_rank(entry) -> tuple[float, ...]:
        nuisance = entry[2]
        fit = nuisance["fits"]["mixed"]
        return (float(nuisance["adequate"]["mixed"]), fit.total_weight,
                float(fit.observations),
                -float(BASELINE_LEVELS.index(entry[0])))

    if ref_entry is None:
        ref_entry = max(level_fits, key=ref_alt_rank)
    if mixed_entry is None:
        mixed_entry = max(level_fits, key=mixed_rank)
    ref_source, mixed_source = ref_entry[2], mixed_entry[2]
    ref_fit = ref_source["fits"]["ref"]
    alt_fit = ref_source["fits"]["alt"]
    mixed_fit = mixed_source["fits"]["mixed"]
    matched_ref_alt = bool(ref_source["matched_ref_alt"])
    adequate = {
        "ref": bool(ref_source["adequate"]["ref"]),
        "alt": bool(ref_source["adequate"]["alt"]),
        "mixed": bool(mixed_source["adequate"]["mixed"]),
    }
    orientation_cells = {
        "ref": ref_fit.observations, "alt": alt_fit.observations,
        "mixed": mixed_fit.observations,
    }
    candidates: list[tuple[float, float]] = []
    if matched_ref_alt:
        candidates.append(((ref_fit.delta + alt_fit.delta) / 2.0,
                           min(ref_fit.total_weight, alt_fit.total_weight)))
    if adequate["mixed"]:
        candidates.append((mixed_fit.delta, mixed_fit.total_weight))
    total = sum(weight for _value, weight in candidates)
    group_baseline = (sum(value * weight for value, weight in candidates) / total
                      if total > 0.0 else 0.0)
    mapping = ((ref_fit.delta - alt_fit.delta) / 2.0
               if matched_ref_alt else 0.0)
    selected_levels = {
        "ref_alt": ref_entry[0], "mixed": mixed_entry[0]}
    level = (ref_entry[0] if ref_entry[0] == mixed_entry[0]
             else f"REF_ALT={ref_entry[0]};MIXED={mixed_entry[0]}")
    nuisance = {
        "group_baseline": group_baseline,
        "offsets": {"ref": mapping, "alt": -mapping, "mixed": 0.0},
        "rhos": {"ref": ref_fit.rho, "alt": alt_fit.rho,
                 "mixed": mixed_fit.rho},
        "cells": max(orientation_cells.values(), default=0),
        "fits": {"ref": ref_fit, "alt": alt_fit, "mixed": mixed_fit},
        "adequate": adequate, "orientation_cells": orientation_cells,
        "matched_ref_alt": matched_ref_alt, "selected_levels": selected_levels,
        "level": level,
        "has_external_reference": any(entry[1] for entry in level_fits),
    }
    cache[key] = nuisance
    return nuisance


def build_ase_null_index(
        cells: Sequence[dict[str, str]],
        rows_by_cell: Mapping[tuple[str, str], list[dict[str, str]]], args,
        calibration_by_cell: Mapping[
            tuple[str, str], Mapping[str, Sequence[
                tuple[str, str, float, float, float, float]]]],
        cache: dict[tuple[str, ...], dict[str, object]],
        chromosome: str) -> IndexedEmpiricalNulls:
    """Score each ASE pseudo-null once per side/depth/ambient stratum."""
    index = IndexedEmpiricalNulls(ASE_DIRECTIONS)
    calibration_cells = sorted(
        (cell for cell in cells
         if truthy(cell.get("ase_calibration_eligible"))),
        key=lambda row: (
            natural_key(row.get("library", "")),
            natural_key(row.get("barcode", ""))))
    started = next_progress = time.monotonic()
    for cell_index, cell in enumerate(calibration_cells, 1):
        if time.monotonic() >= next_progress + 240.0:
            timings = model_timings(args)
            print(
                f"MODEL chromosome {chromosome}: ASE null progress "
                f"{cell_index}/{len(calibration_cells)} cells, "
                f"{index.score_evaluations} evaluations, "
                f"{len(cache)} fits, {time.monotonic() - started:.1f}s; "
                f"level_fit={timings['ase_level_fit_seconds']:.1f}s "
                f"baseline={timings['ase_baseline_seconds']:.1f}s "
                f"state={timings['ase_null_state_seconds']:.1f}s "
                f"baseline_fits={timings['ase_baseline_fits']:.0f}",
                flush=True)
            next_progress = time.monotonic()
        candidates = sorted(
            rows_by_cell[cell_key(cell)],
            key=lambda row: natural_key(row.get("arm", "")))
        for candidate in candidates:
            if (ase_scoring_safety_status(candidate, args) != "PASS"
                    or total_ase_weight(candidate) <= 0.0):
                continue
            index.score_evaluations += 1
            candidate_nuisance = ase_nuisance_for_query(
                candidate, cells, args, calibration_by_cell, cache, chromosome)
            candidate_status, _candidate_cells = (
                ase_orientation_calibration_status(
                    candidate, candidate_nuisance, args))
            if candidate_status != "PASS_EMPIRICAL_CALIBRATION":
                continue
            baseline, _status = cell_ase_baseline(
                candidate, candidate_nuisance, args, calibration_by_cell,
                chromosome)
            score_started = time.monotonic()
            log_bfs = ase_state_log_bfs(
                ase_evidence(candidate), baseline,
                candidate_nuisance["offsets"], candidate_nuisance["rhos"],
                args.site_fallback_weight)
            timings = model_timings(args)
            timings["ase_null_state_seconds"] += (
                time.monotonic() - score_started)
            timings["ase_null_state_calls"] += 1
            index.add(
                ase_null_pool_key(candidate), candidate,
                ase_direction_statistics(log_bfs))
    print(
        f"MODEL chromosome {chromosome}: finalizing ASE null index "
        f"({index.retained_statistics} statistics, {len(index.pools)} pools)",
        flush=True)
    finalize_started = time.monotonic()
    result = index.finalize()
    model_timings(args)["ase_finalize_seconds"] += (
        time.monotonic() - finalize_started)
    return result


def ase_model_for_query(
        query: dict[str, str], cells: Sequence[dict[str, str]], args,
        calibration_by_cell: Mapping[
            tuple[str, str], Mapping[str, Sequence[
                tuple[str, str, float, float, float, float]]]],
        null_index: IndexedEmpiricalNulls,
        cache: dict[tuple[str, ...], dict[str, object]],
        calibration_rows: dict[tuple[str, ...], dict[str, object]],
        chromosome: str) -> dict[str, object]:
    key = calibration_cache_key(query, "ASE")
    model_key = ("ASE_MODEL",) + key[1:]
    if model_key in cache:
        return cache[model_key]
    nuisance = ase_nuisance_for_query(
        query, cells, args, calibration_by_cell, cache, chromosome)
    orientation_status, target_calibration_cells = (
        ase_orientation_calibration_status(query, nuisance, args))
    pool_keys = (ase_null_pool_key(query),)
    null_count = null_index.count(pool_keys, query)
    calibration_key = cache_key_text(key)
    if (orientation_status == "PASS_EMPIRICAL_CALIBRATION"
            and null_count >= args.min_ase_null_cells):
        status = "PASS"
    elif not nuisance["has_external_reference"]:
        status = "NO_EXTERNAL_REFERENCE"
    elif orientation_status != "PASS_EMPIRICAL_CALIBRATION":
        status = orientation_status
    else:
        status = "WEAK_COMPONENT_LOW_SUPPORT"
    result = {
        "key": calibration_key, "level": nuisance["level"],
        "calibration_cells": target_calibration_cells, "nuisance": nuisance,
        "null_count": null_count, "null_rows": null_count,
        "null_pool_keys": pool_keys,
        "null_libraries": null_index.libraries(pool_keys, query),
        "status": status,
    }
    cache[model_key] = result
    calibration_rows[("ASE", calibration_key, "full")] = {
        "chromosome": chromosome, "calibration_key": calibration_key,
        "branch": "ASE", "fold": "full",
        "calibration_level": nuisance["level"],
        "calibration_cells": target_calibration_cells,
        "null_cells": null_count,
        "input_state": "MOLECULE_LEVEL_SOFT_ASE", "depth_bin": "NA",
        "breadth_bin": "NA", "center_p": "NA", "center_q": "NA",
        "covariance_pp": "NA", "covariance_pq": "NA", "covariance_qq": "NA",
        "group_baseline_logit": f"{float(nuisance['group_baseline']):.17g}",
        "mapping_offset_ref": f"{float(nuisance['offsets']['ref']):.17g}",
        "mapping_offset_alt": f"{float(nuisance['offsets']['alt']):.17g}",
        "mapping_offset_mixed": f"{float(nuisance['offsets']['mixed']):.17g}",
        "rho_ref": f"{float(nuisance['rhos']['ref']):.17g}",
        "rho_alt": f"{float(nuisance['rhos']['alt']):.17g}",
        "rho_mixed": f"{float(nuisance['rhos']['mixed']):.17g}",
        "ref_alt_calibration_level": nuisance["selected_levels"]["ref_alt"],
        "mixed_calibration_level": nuisance["selected_levels"]["mixed"],
        "calibration_cells_ref": nuisance["orientation_cells"]["ref"],
        "calibration_cells_alt": nuisance["orientation_cells"]["alt"],
        "calibration_cells_mixed": nuisance["orientation_cells"]["mixed"],
        "status": status, "schema_version": HYBRID_CALIBRATION_SCHEMA,
    }
    return result


def expression_scores_for_row(
        row: dict[str, str], cells: Sequence[dict[str, str]], args,
        null_index: IndexedEmpiricalNulls,
        cache: dict[tuple[str, ...], dict[str, object]],
        component_rows: dict[tuple[str, ...], dict[str, object]],
        calibration_rows: dict[tuple[str, ...], dict[str, object]],
    chromosome: str) -> dict[str, object]:
    side = arm_side(row.get("arm"))
    if not side or not row_has_side_expression(row, side):
        return {
            "eligible": False, "status": "NO_DATA", "level": "NA", "cells": 0,
            "component": "NO_DATA", "loss_lbf": math.nan, "gain_lbf": math.nan,
            "loss_p": math.nan, "gain_p": math.nan, "p_floor": math.nan,
            "fold0_state": "NO_DATA", "fold1_state": "NO_DATA",
            "fold_status": "INSUFFICIENT", "best_state": "NO_DATA",
        }
    models = {
        fold: expression_model_for_query(
            row, cells, fold, args, null_index, cache, component_rows,
            calibration_rows, chromosome)
        for fold in ("full", "0", "1")
    }
    full = models["full"]
    point = row_point(row)
    full_lbf = expression_direction_log_bfs(
        point, side, full["fit"].center, full["fit"].covariance,
        args.expression_df, full["mixture"].weights,
        args.expression_outlier_scale)
    component = best_expression_component(
        point, full["fit"].center, full["fit"].covariance,
        full["mixture"].weights, args.expression_df,
        args.expression_outlier_scale)
    fold_states: dict[str, str] = {}
    fold_p: dict[str, dict[str, float]] = {}
    for fold in ("0", "1"):
        model = models[fold]
        if not row_has_side_expression(row, side, fold):
            fold_states[fold] = "NO_DATA"
            fold_p[fold] = {direction: math.nan for direction in COPY_STATES}
            continue
        scores = expression_direction_log_bfs(
            row_point(row, fold), side, model["fit"].center,
            model["fit"].covariance, args.expression_df,
            model["mixture"].weights, args.expression_outlier_scale)
        fold_states[fold] = direction_from_log_bfs(
            scores["LOSS"], scores["GAIN"], args.min_expression_log_bf)
        fold_p[fold] = {
            direction: null_index.pvalue(
                max(0.0, scores[direction]), direction,
                model["null_pool_keys"], row)
            for direction in COPY_STATES
        }
    if fold_states["0"] == fold_states["1"] and fold_states["0"] in COPY_STATES:
        fold_status = "CONCORDANT"
        best_state = fold_states["0"]
    elif fold_states["0"] == fold_states["1"] == "BALANCED":
        fold_status = "CONCORDANT_BALANCED"
        best_state = "BALANCED"
    elif {fold_states["0"], fold_states["1"]} == {"LOSS", "GAIN"}:
        fold_status = "DISAGREE"
        best_state = "OUTLIER"
    else:
        fold_status = "INSUFFICIENT"
        best_state = direction_from_log_bfs(
            full_lbf["LOSS"], full_lbf["GAIN"], args.min_expression_log_bf)
    if component == "OUTLIER":
        # The broad component is an explicit fitted alternative to every
        # dosage state.  It must remain a veto even when both gene folds happen
        # to point in the same direction.
        fold_status = "OUTLIER_COMPONENT"
        best_state = "OUTLIER"
    p_values = {
        direction: max(fold_p["0"][direction], fold_p["1"][direction])
        if (math.isfinite(fold_p["0"][direction])
            and math.isfinite(fold_p["1"][direction])) else math.nan
        for direction in COPY_STATES
    }
    null_count = min(
        (int(models[fold]["null_count"]) for fold in ("0", "1")),
        default=0)
    breadth_ok = expression_breadth_pass(row, args)
    all_models_pass = all(models[fold]["status"] == "PASS"
                          for fold in ("full", "0", "1"))
    replicated = fold_status in {"CONCORDANT", "CONCORDANT_BALANCED"}
    status = "PASS" if all_models_pass and breadth_ok and replicated else (
        "EXPRESSION_OUTLIER"
        if fold_status in {"DISAGREE", "OUTLIER_COMPONENT"} else
        "INSUFFICIENT_GENE_BREADTH" if not breadth_ok else
        "INSUFFICIENT_FOLD_SUPPORT" if not replicated else
        "INSUFFICIENT_REFERENCE")
    return {
        "eligible": True, "status": status, "level": full["level"],
        "cells": full["references"], "component": component,
        "loss_lbf": full_lbf["LOSS"], "gain_lbf": full_lbf["GAIN"],
        "loss_p": p_values["LOSS"], "gain_p": p_values["GAIN"],
        "p_floor": empirical_p_floor(null_count) if null_count else math.nan,
        "fold0_state": fold_states["0"], "fold1_state": fold_states["1"],
        "fold_status": fold_status, "best_state": best_state,
        "breadth_ok": breadth_ok,
    }


def ase_scores_for_row(
        row: dict[str, str], cells: Sequence[dict[str, str]],
        args, null_index: IndexedEmpiricalNulls,
        calibration_by_cell: Mapping[
            tuple[str, str],
            Mapping[str, Sequence[
                tuple[str, str, float, float, float, float]]]],
        cache: dict[tuple[str, ...], dict[str, object]],
        calibration_rows: dict[tuple[str, ...], dict[str, object]],
        chromosome: str) -> dict[str, object]:
    if not truthy(row.get("hybrid_target_eligible")):
        return {"eligible": False, "status": "NOT_TARGET_ELIGIBLE"}
    safety_status = ase_scoring_safety_status(row, args)
    if safety_status != "PASS":
        return {"eligible": False, "status": safety_status}
    sites = int(finite_float(row.get("ase_n_sites"), 0.0))
    effective_weight = total_ase_weight(row)
    if sites < 1:
        return {"eligible": False, "status": "INSUFFICIENT_POWER_NO_SITES"}
    if effective_weight <= 0.0:
        return {
            "eligible": False,
            "status": "INSUFFICIENT_POWER_NO_EFFECTIVE_WEIGHT",
        }
    low_power = []
    if sites < args.min_ase_sites:
        low_power.append("LOW_SITES")
    if effective_weight < args.min_ase_effective_weight:
        low_power.append("LOW_EFFECTIVE_WEIGHT")
    model = ase_model_for_query(
        row, cells, args, calibration_by_cell, null_index, cache,
        calibration_rows, chromosome)
    nuisance = model["nuisance"]
    baseline, baseline_status = cell_ase_baseline(
        row, nuisance, args, calibration_by_cell, chromosome)
    log_bfs = ase_state_log_bfs(
        ase_evidence(row), baseline, nuisance["offsets"], nuisance["rhos"],
        args.site_fallback_weight)
    statistics = ase_direction_statistics(log_bfs)
    p_values = {
        direction: null_index.pvalue(
            statistics[direction], direction, model["null_pool_keys"], row)
        for direction in ASE_DIRECTIONS
    }
    null_count = int(model["null_count"])
    best_exact = max(EXACT_STATES, key=lambda state: (log_bfs[state], -EXACT_STATES.index(state)))
    best_direction = max(
        ASE_DIRECTIONS,
        key=lambda direction: (statistics[direction], -ASE_DIRECTIONS.index(direction)))
    if statistics[best_direction] <= args.min_ase_log_bf:
        best_direction = "BALANCED"
    status = (
        str(model["status"]) if model["status"] != "PASS"
        else "INSUFFICIENT_POWER_" + "_AND_".join(low_power)
        if low_power else "PASS")
    return {
        "eligible": True, "status": status, "level": model["level"],
        "cells": model["calibration_cells"], "baseline_status": baseline_status,
        "log_bfs": log_bfs, "statistics": statistics, "p_values": p_values,
        "p_floor": empirical_p_floor(null_count) if null_count else math.nan,
        "best_exact": best_exact, "best_direction": best_direction,
    }


def main_impl(args) -> int:
    model_started = time.monotonic()
    chromosome = canonical_chromosome(args.chromosome)
    if not chromosome:
        raise ValueError("--chromosome is empty")
    if args.min_expression_reference_cells < 1 or args.min_ase_calibration_cells < 1:
        raise ValueError("calibration cell minima must be positive")
    if (args.ase_calibration_max_weight <= 0.0
            or args.ase_calibration_huber_z <= 0.0
            or not 0.0 <= args.ase_calibration_min_robust_weight <= 1.0
            or args.min_ase_orientation_calibration_effective_weight <= 0.0):
        raise ValueError("invalid ASE robust-calibration setting")
    if not 0.0 < args.provisional_p_threshold <= 1.0:
        raise ValueError("--provisional-p-threshold must be in (0,1]")
    if (not 0.0 <= args.max_ase_qname_fallback_fraction <= 1.0
            or not 0.0 <= args.min_ase_ambient_genotyped_mass <= 1.0):
        raise ValueError("invalid ASE safety-gate threshold")
    rows, shard_paths, calibration_by_cell, observed_libraries = load_shards(
        args.shard_manifest, chromosome)
    target_libraries = parse_target_libraries(args.target_libraries)
    if target_libraries and not target_libraries <= observed_libraries:
        missing = sorted(target_libraries - observed_libraries, key=natural_key)
        raise ValueError(
            "--target-libraries are absent from the calibration cohort: "
            + ",".join(missing))
    if not target_libraries:
        target_libraries = observed_libraries
    cells, rows_by_cell = make_cell_records(rows)
    print(
        f"MODEL chromosome {chromosome}: loaded {len(rows)} arm rows, "
        f"{len(cells)} cells, {len(observed_libraries)} libraries",
        flush=True)

    output_dir = Path(os.path.abspath(args.output_dir))
    prefix = output_dir / f"tetra_arm_hybrid_{chromosome}"
    scores_path = Path(str(prefix) + ".scores.tsv.gz")
    components_path = Path(str(prefix) + ".expression_components.tsv.gz")
    calibration_path = Path(str(prefix) + ".calibration.tsv.gz")
    qc_path = Path(str(prefix) + ".qc.tsv")
    contract_path = Path(str(prefix) + ".contract.json")
    require_outputs_absent((scores_path, components_path, calibration_path,
                            qc_path, contract_path))

    expression_cache: dict[tuple[str, ...], dict[str, object]] = {}
    ase_cache: dict[tuple[str, ...], dict[str, object]] = {}
    component_rows: dict[tuple[str, ...], dict[str, object]] = {}
    calibration_rows: dict[tuple[str, ...], dict[str, object]] = {}
    phase_started = time.monotonic()
    expression_null_index = build_expression_null_index(
        cells, args, expression_cache, component_rows, chromosome)
    print(
        f"MODEL chromosome {chromosome}: expression null index complete in "
        f"{time.monotonic() - phase_started:.1f}s "
        f"({expression_null_index.score_evaluations} evaluations, "
        f"{len(expression_cache)} nuisance/model fits)",
        flush=True)
    timings = model_timings(args)
    print(
        f"MODEL chromosome {chromosome}: expression fit breakdown "
        f"keys={timings['expression_nuisance_fits']:.0f} "
        f"reference_cells={timings['expression_reference_cells']:.0f} "
        f"discovery_cells={timings['expression_discovery_cells']:.0f} "
        f"pool={timings['expression_pool_seconds']:.1f}s "
        f"robust={timings['expression_robust_seconds']:.1f}s "
        f"discovery={timings['expression_discovery_seconds']:.1f}s "
        f"mixture={timings['expression_mixture_seconds']:.1f}s "
        f"score={timings['expression_null_score_seconds']:.1f}s "
        f"finalize={timings['expression_finalize_seconds']:.1f}s",
        flush=True)
    phase_started = time.monotonic()
    ase_null_index = build_ase_null_index(
        cells, rows_by_cell, args, calibration_by_cell, ase_cache, chromosome)
    print(
        f"MODEL chromosome {chromosome}: ASE null index complete in "
        f"{time.monotonic() - phase_started:.1f}s "
        f"({ase_null_index.score_evaluations} evaluations, "
        f"{len(ase_cache)} nuisance/model fits)",
        flush=True)
    print(
        f"MODEL chromosome {chromosome}: ASE fit breakdown "
        f"level_fits={timings['ase_level_fits']:.0f} "
        f"selected_cells={timings['ase_selected_cells']:.0f} "
        f"level_fit={timings['ase_level_fit_seconds']:.1f}s "
        f"baseline_fits={timings['ase_baseline_fits']:.0f} "
        f"baseline_hits={timings['ase_baseline_hits']:.0f} "
        f"baseline={timings['ase_baseline_seconds']:.1f}s "
        f"state_calls={timings['ase_null_state_calls']:.0f} "
        f"state={timings['ase_null_state_seconds']:.1f}s "
        f"finalize={timings['ase_finalize_seconds']:.1f}s",
        flush=True)
    phase_started = time.monotonic()
    output_rows: list[dict[str, object]] = []
    status_counts = Counter()
    for row in sorted(rows, key=lambda value: (
            natural_key(value.get("library", "")),
            natural_key(value.get("barcode", "")),
            natural_key(value.get("arm", "")))):
        if (clean(row.get("library")) not in target_libraries
                or not truthy(row.get("hybrid_target_eligible"))):
            continue
        expression = expression_scores_for_row(
            row, cells, args, expression_null_index, expression_cache,
            component_rows,
            calibration_rows, chromosome)
        ase = ase_scores_for_row(
            row, cells, args, ase_null_index, calibration_by_cell,
            ase_cache, calibration_rows, chromosome)

        expression_p = {
            "LOSS": finite_float(expression.get("loss_p")),
            "GAIN": finite_float(expression.get("gain_p")),
        }
        ase_p = {
            "DONOR_A_DEPLETED": finite_float(
                (ase.get("p_values") or {}).get("DONOR_A_DEPLETED")),
            "DONOR_A_ENRICHED": finite_float(
                (ase.get("p_values") or {}).get("DONOR_A_ENRICHED")),
        }
        conjunction = conjunction_pvalues(expression_p, ase_p)
        expression_significant = (
            expression.get("status") == "PASS"
            and expression.get("best_state") in COPY_STATES
            and math.isfinite(expression_p[expression["best_state"]])
            and expression_p[expression["best_state"]]
            <= args.provisional_p_threshold)
        ase_direction = clean(ase.get("best_direction")) or "NO_DATA"
        ase_significant = (
            ase.get("status") == "PASS" and ase_direction in ASE_DIRECTIONS
            and math.isfinite(ase_p[ase_direction])
            and ase_p[ase_direction] <= args.provisional_p_threshold)
        provisional_state = AXES_TO_STATE.get((
            clean(expression.get("best_state")), ase_direction), "NO_CALL")
        conjunction_significant = (
            provisional_state in EXACT_STATES
            and math.isfinite(conjunction[provisional_state])
            and conjunction[provisional_state] <= args.provisional_p_threshold)

        confounding = split_flags(row.get("confounding_flags"))
        if not expression.get("breadth_ok", False) and expression.get("eligible"):
            confounding.add("INSUFFICIENT_GENE_BREADTH")
        ambient = finite_float(row.get("expression_model_ambient_c"))
        if (expression_significant and math.isfinite(ambient)
                and ambient >= args.high_ambient_threshold
                and clean(row.get("expression_model_expression_input_state"))
                == "OBSERVED_FILTERED_COUNTS"):
            confounding.add("HIGH_AMBIENT_EXPRESSION_SENSITIVITY")
        resolved = resolve_evidence(
            clean(expression.get("best_state")), expression_significant,
            ase_direction, ase_significant, provisional_state,
            conjunction_significant, clean(expression.get("fold_status")),
            bool(expression.get("status") == "PASS" and
                 ase.get("status") == "PASS"), confounding)
        if expression.get("fold_status") == "DISAGREE":
            discordance_reason = "GENE_FOLD_DIRECTION_DISAGREEMENT"
        elif expression.get("fold_status") == "OUTLIER_COMPONENT":
            discordance_reason = "FITTED_EXPRESSION_OUTLIER_COMPONENT"
        elif resolved["evidence_class"] == "DISCORDANT":
            discordance_reason = "BRANCH_SIGNIFICANCE_WITHOUT_EXACT_CONJUNCTION"
        else:
            discordance_reason = "NONE"

        nonzero, top1, top5, top10 = expression_side_metrics(row)
        log_bfs = ase.get("log_bfs") or {}
        output: dict[str, object] = {
            "library": row["library"], "barcode": row["barcode"],
            "calibration_group": row["calibration_group"], "uid": row["uid"],
            "donor_a": row["donor_a"], "donor_b": row["donor_b"],
            "donor_pair": row["donor_pair"], "chromosome": chromosome,
            "arm": row["arm"], "hybrid_target_eligible": 1,
            "ase_eligible": int(bool(ase.get("eligible"))),
            "expression_eligible": int(bool(expression.get("eligible"))),
            "ase_status": ase.get("status", "NO_DATA"),
            "expression_status": expression.get("status", "NO_DATA"),
            "ase_calibration_level": ase.get("level", "NA"),
            "expression_calibration_level": expression.get("level", "NA"),
            "ase_calibration_cells": ase.get("cells", 0),
            "expression_calibration_cells": expression.get("cells", 0),
            "ase_best_exact_state": ase.get("best_exact", "NO_DATA"),
            "ase_direction": ase_direction,
            "ase_p_DONOR_A_DEPLETED": format_number(
                ase_p["DONOR_A_DEPLETED"]),
            "ase_p_DONOR_A_ENRICHED": format_number(
                ase_p["DONOR_A_ENRICHED"]),
            "ase_p_floor": format_number(ase.get("p_floor")),
            "expression_component": expression.get("component", "NO_DATA"),
            "expression_best_copy_state": expression.get("best_state", "NO_DATA"),
            "expression_log_bf_LOSS": format_number(expression.get("loss_lbf")),
            "expression_log_bf_GAIN": format_number(expression.get("gain_lbf")),
            "expression_p_LOSS": format_number(expression_p["LOSS"]),
            "expression_p_GAIN": format_number(expression_p["GAIN"]),
            "expression_fold0_state": expression.get("fold0_state", "NO_DATA"),
            "expression_fold1_state": expression.get("fold1_state", "NO_DATA"),
            "expression_fold_replication_status": expression.get(
                "fold_status", "INSUFFICIENT"),
            "expression_top1_fraction": format_number(top1),
            "expression_top5_fraction": format_number(top5),
            "expression_top10_fraction": format_number(top10),
            "expression_nonzero_genes": nonzero,
            "expression_p_floor": format_number(expression.get("p_floor")),
            "expression_input_state": clean(
                row.get("expression_model_expression_input_state")) or "NA",
            "expression_depth_bin": row.get("expression_depth_bin", "NA"),
            "expression_breadth_bin": row.get("expression_breadth_bin", "NA"),
            "provisional_evidence_class": resolved["evidence_class"],
            "provisional_resolved_state": provisional_state,
            "discordance_reason": discordance_reason,
            "confounding_flags": resolved["confounding_flags"],
            "qc_flags": join_flags(
                [row.get("qc_flags"),
                 "ASE_AMBIENT_UNCERTAINTY_UNAVAILABLE"
                 if not math.isfinite(finite_float(row.get("ase_ambient_c_se")))
                 or finite_float(row.get("ase_ambient_c_se"), 0.0) <= 0.0
                 else "PASS"]),
            "schema_version": HYBRID_SCORE_SCHEMA,
        }
        for state in EXACT_STATES:
            output[f"ase_log_bf_{state}"] = format_number(log_bfs.get(state))
            output[f"hybrid_conjunction_p_{state}"] = format_number(
                conjunction[state])
        output_rows.append(output)
        status_counts[(output["ase_status"], output["expression_status"])] += 1

    output_rows.sort(key=lambda row: (
        natural_key(row["library"]), natural_key(row["barcode"]),
        natural_key(row["arm"])))
    output_dir.mkdir(parents=True, exist_ok=True)
    score_count = write_tsv_atomic(
        str(scores_path), output_rows, SCORE_FIELDS, deterministic_gzip=True)
    component_values = [component_rows[key] for key in sorted(
        component_rows, key=lambda key: tuple(natural_key(part) for part in key))]
    calibration_values = [calibration_rows[key] for key in sorted(
        calibration_rows, key=lambda key: tuple(natural_key(part) for part in key))]
    write_tsv_atomic(
        str(components_path), component_values, COMPONENT_FIELDS,
        deterministic_gzip=True)
    write_tsv_atomic(
        str(calibration_path), calibration_values, CALIBRATION_FIELDS,
        deterministic_gzip=True)
    terminal = "NONE" if score_count else "PASS_NO_TARGET_ROWS"
    qc_rows = [
        {"metric": "schema_version", "value": "tetra_arm_hybrid_model_qc_v1"},
        {"metric": "chromosome", "value": chromosome},
        {"metric": "input_rows", "value": len(rows)},
        {"metric": "calibration_libraries", "value": len(observed_libraries)},
        {"metric": "target_libraries", "value": len(target_libraries)},
        {"metric": "target_score_rows", "value": score_count},
        {"metric": "expression_calibration_models", "value": len(expression_cache)},
        {"metric": "ase_calibration_models", "value": len(ase_cache)},
        {"metric": "expression_null_score_evaluations",
         "value": expression_null_index.score_evaluations},
        {"metric": "expression_null_statistics",
         "value": expression_null_index.retained_statistics},
        {"metric": "expression_null_pools",
         "value": len(expression_null_index.pools)},
        {"metric": "ase_null_score_evaluations",
         "value": ase_null_index.score_evaluations},
        {"metric": "ase_null_statistics",
         "value": ase_null_index.retained_statistics},
        {"metric": "ase_null_pools", "value": len(ase_null_index.pools)},
        {"metric": "terminal_state", "value": terminal},
        {"metric": "status", "value": "PASS"},
    ]
    write_tsv_atomic(str(qc_path), qc_rows, ["metric", "value"])
    write_json_atomic(str(contract_path), {
        "schema_version": "tetra_arm_hybrid_model_contract_v1",
        "release": PROGRAM_VERSION, "chromosome": chromosome,
        "inputs": {"shard_manifest": file_record(args.shard_manifest),
                   "shards": [file_record(path) for path in shard_paths]},
        "outputs": {
            "scores": str(scores_path), "expression_components": str(components_path),
            "calibration": str(calibration_path), "qc": str(qc_path),
            "contract": str(contract_path),
        },
        "output_schemas": {
            "scores": HYBRID_SCORE_SCHEMA,
            "expression_components": HYBRID_COMPONENT_SCHEMA,
            "calibration": HYBRID_CALIBRATION_SCHEMA,
        },
        "score_fields": SCORE_FIELDS, "scores": score_count,
        "calibration_libraries": sorted(observed_libraries, key=natural_key),
        "target_libraries": sorted(target_libraries, key=natural_key),
        "component_rows": len(component_values),
        "calibration_rows": len(calibration_values),
        "multiplicity_correction": "NONE_CHROMOSOME_SHARD_P_VALUES_ONLY",
        "crossfit_contract": {
            "nuisance_and_component_fits": "EXCLUDE_QUERY_BIOLOGICAL_FOLD_PAIR_UID",
            "empirical_nulls": (
                "CALCULATE_ONCE_COHORT_INDEX_WITH_EACH_PSEUDO_NULL_SCORED_"
                "BY_ITS_OWN_PAIR_UID_FOLD_EXCLUDED_FIT_AND_QUERY_BLOCK_"
                "TAIL_SUBTRACTION"),
            "tail_queries": (
                "SORTED_EXACT_STRATUM_VECTORS_BINARY_SEARCH_PLUS_PAIR_UID_"
                "CELL_INCLUSION_EXCLUSION"),
            "ase_matching": "ARM_SIDE_LOG2_DEPTH_BIN_AMBIENT_0.10_BIN",
            "ase_bias_observations": "CELL_POOLED_LOCO",
            "ase_rho_and_cell_baseline_observations": "INDIVIDUAL_CELL_ARM_LOCO",
        },
        "expression_direction_statistic": (
            "COMMON_SISTER_ARM_NUISANCE_MARGINAL_COMPOSITE_LIKELIHOOD_"
            "RATIO_WITH_CROSSFITTED_POPULATION_MIXTURE_WEIGHTS"),
        "ase_orientation_calibration": (
            "TARGET_USED_ORIENTATIONS_ONLY_WITH_INDEPENDENT_REF_ALT_AND_"
            "MIXED_HIERARCHICAL_FALLBACK"),
        "empirical_p_value_resolution": (
            "PLUS_ONE_NULLS_FROM_RESOLVED_UPSTREAM_CALIBRATION_COHORT"),
        "empirical_null_index": {
            "expression_score_evaluations": (
                expression_null_index.score_evaluations),
            "expression_statistics": (
                expression_null_index.retained_statistics),
            "expression_pools": len(expression_null_index.pools),
            "ase_score_evaluations": ase_null_index.score_evaluations,
            "ase_statistics": ase_null_index.retained_statistics,
            "ase_pools": len(ase_null_index.pools),
        },
        "ase_target_safety_gates": (
            "REJECT_SITE_FALLBACK_HIGH_QNAME_LOW_AMBIENT_MASS_"
            "MISSING_OR_NONPOSITIVE_AMBIENT_SE"),
        "expression_ase_independence_contract": (
            "ASE_VALUES_DO_NOT_ENTER_EXPRESSION_FITTING_SELECTION_WEIGHTS_OR_PRIORS"),
        "terminal_state": terminal, "status": "PASS",
    })
    print(
        f"MODEL chromosome {chromosome}: wrote {score_count} score rows in "
        f"{time.monotonic() - phase_started:.1f}s; total "
        f"{time.monotonic() - model_started:.1f}s",
        flush=True)
    return 0


def self_test() -> int:
    expression_p = {"LOSS": 0.01, "GAIN": 0.8}
    ase_p = {"DONOR_A_DEPLETED": 0.02, "DONOR_A_ENRICHED": 0.9}
    values = conjunction_pvalues(expression_p, ase_p)
    if values["DONOR_A_LOSS"] != 0.02 or values["DONOR_B_GAIN"] != 0.8:
        raise AssertionError("partial-conjunction self-test failed")
    print("PASS tetra_arm_hybrid_model self-test")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fit one chromosome's independent hybrid expression/ASE branches.")
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {PROGRAM_VERSION}")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--chromosome")
    parser.add_argument("--shard-manifest")
    parser.add_argument("--output-dir")
    parser.add_argument(
        "--target-libraries", default="",
        help=("Comma-separated libraries to score; all shard libraries remain "
              "available as cross-fitted calibration/null evidence."))
    parser.add_argument("--min-expression-reference-cells", type=int, default=20)
    parser.add_argument("--min-expression-nonzero-genes", type=int, default=10)
    parser.add_argument("--max-expression-top1-fraction", type=float, default=0.50)
    parser.add_argument("--max-expression-top5-fraction", type=float, default=0.80)
    parser.add_argument("--max-expression-top10-fraction", type=float, default=0.95)
    parser.add_argument("--min-expression-log-bf", type=float, default=0.0)
    parser.add_argument("--min-expression-sigma", type=float, default=0.10)
    parser.add_argument("--max-expression-sigma", type=float, default=2.0)
    parser.add_argument("--max-expression-abs-correlation", type=float, default=0.95)
    parser.add_argument("--expression-df", type=float, default=4.0)
    parser.add_argument("--expression-outlier-scale", type=float, default=16.0)
    parser.add_argument("--min-ase-sites", type=int, default=3)
    parser.add_argument("--min-ase-effective-weight", type=float, default=8.0)
    parser.add_argument("--min-ase-calibration-cells", type=int, default=12)
    parser.add_argument("--min-ase-null-cells", type=int, default=20)
    parser.add_argument(
        "--min-ase-orientation-calibration-effective-weight",
        type=float, default=20.0)
    parser.add_argument("--min-ase-log-bf", type=float, default=0.0)
    parser.add_argument("--ase-calibration-max-weight", type=float, default=50.0)
    parser.add_argument("--ase-calibration-huber-z", type=float, default=2.5)
    parser.add_argument(
        "--ase-calibration-min-robust-weight", type=float, default=0.25)
    parser.add_argument("--ase-cell-baseline-shrinkage-arms", type=float, default=8.0)
    parser.add_argument("--min-ase-cell-baseline-arms", type=int, default=4)
    parser.add_argument("--default-rho", type=float, default=0.02)
    parser.add_argument("--min-rho", type=float, default=0.001)
    parser.add_argument("--max-rho", type=float, default=0.25)
    parser.add_argument("--site-fallback-weight", type=float, default=0.50)
    parser.add_argument(
        "--max-ase-qname-fallback-fraction", type=float, default=0.50)
    parser.add_argument(
        "--min-ase-ambient-genotyped-mass", type=float, default=0.50)
    parser.add_argument("--provisional-p-threshold", type=float, default=0.05)
    parser.add_argument("--high-ambient-threshold", type=float, default=0.25)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.self_test:
            return self_test()
        missing = [name for name in ("chromosome", "shard_manifest", "output_dir")
                   if getattr(args, name) in {None, ""}]
        if missing:
            raise ValueError("missing required option(s): " + ", ".join(missing))
        return main_impl(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
