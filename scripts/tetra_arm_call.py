#!/usr/bin/env python3
"""Empirically calibrated donor-directed chromosome-arm CNV caller.

The caller deliberately separates the two evidence types:

* molecule-aware soft ASE supplies donor direction and is modelled with an
  ambient-adjusted, orientation-specific quasi-likelihood;
* arm-level expression is an optional, supportive total-copy measurement used
  to distinguish a loss from the reciprocal donor gain.  In the absence of a
  matching ambient gene-expression profile it is attenuated, not claimed to be
  exactly decontaminated.

Calibration controls must be marked ``calibration_eligible`` in the prepared
cell manifest and must pass the caller's revalidated context gates.  Identity,
ploidy, occupancy, technical, species, mitochondrial, ATAC, review, and release
fields are propagated for audit and eligibility but are never multiplied into
the CNV likelihood.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import math
import os
import sqlite3
import statistics
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple

from tetra_arm_common import (
    CELL_MANIFEST_SCHEMA,
    EXPRESSION_SCHEMA,
    RELEASE,
    canonical_barcode,
    clean,
    finite_float,
    natural_key,
    open_text,
    require_file,
    require_outputs_absent,
    truthy,
    write_json_atomic,
    write_tsv_atomic,
)


PROGRAM_VERSION = RELEASE
ASE_INPUT_SCHEMA = "tetra_arm_ase_evidence_v2"
CALL_OUTPUT_SCHEMA = "tetra_arm_cnv_calls_v2"
CALIBRATION_OUTPUT_SCHEMA = "tetra_arm_calibration_v2"
UID_OUTPUT_SCHEMA = "tetra_arm_uid_chromosome_flags_v2"
PAIR_OUTPUT_SCHEMA = "tetra_arm_donor_pair_arm_summary_v2"
QC_OUTPUT_SCHEMA = "tetra_arm_call_qc_v2"
CONTRACT_OUTPUT_SCHEMA = "tetra_arm_call_contract_v2"

STATES: Tuple[str, ...] = (
    "BALANCED",
    "DONOR_A_LOSS",
    "DONOR_B_LOSS",
    "DONOR_A_GAIN",
    "DONOR_B_GAIN",
)

# Copy-number states act as odds multipliers relative to the cross-fitted
# cell-specific donor baseline.  The shift is applied to cellular donor odds
# before ambient mixing and orientation mapping bias.
STATE_LOG_ODDS_SHIFT = {
    "BALANCED": 0.0,
    "DONOR_A_LOSS": math.log(1.0 / 2.0),
    "DONOR_B_LOSS": math.log(2.0),
    "DONOR_A_GAIN": math.log(3.0 / 2.0),
    "DONOR_B_GAIN": math.log(2.0 / 3.0),
}

TOTAL_COPY_RATIO = {
    "BALANCED": 1.0,
    "DONOR_A_LOSS": 3.0 / 4.0,
    "DONOR_B_LOSS": 3.0 / 4.0,
    "DONOR_A_GAIN": 5.0 / 4.0,
    "DONOR_B_GAIN": 5.0 / 4.0,
}

ORIENTATIONS: Tuple[str, ...] = ("ref", "alt", "mixed")
ORIENTATION_FIELDS = {
    "ref": ("a_ref", "b_ref", "ambient_a_ref"),
    "alt": ("a_alt", "b_alt", "ambient_a_alt"),
    "mixed": ("a_mixed", "b_mixed", "ambient_a_mixed"),
}

# Stable producer contract from tetra_arm_ase.cpp.  These fields are copied
# byte-for-byte for observed rows.  Synthesized NO_DATA rows use the same keys
# and explicit NA/zero values rather than disappearing from the result.
ASE_RAW_FIELDS: Tuple[str, ...] = (
    "library", "barcode", "donor_a", "donor_b", "donor_pair", "arm",
    "chromosome", "arm_start", "arm_end", "ambient_c", "ambient_c_se",
    "n_sites", "n_molecules", "n_informative_molecules", "a_ref", "b_ref",
    "a_alt", "b_alt", "a_mixed", "b_mixed", "n_ambiguous", "soft_a",
    "soft_b", "n_molecules_ref", "n_molecules_alt", "n_molecules_mixed",
    "soft_a_ref", "soft_b_ref", "soft_a_alt", "soft_b_alt",
    "soft_a_mixed", "soft_b_mixed", "soft_a_sumsq_ref",
    "soft_a_sumsq_alt", "soft_a_sumsq_mixed", "effective_a_ref", "effective_b_ref",
    "effective_a_alt", "effective_b_alt", "effective_a_mixed",
    "effective_b_mixed", "effective_weight_ref", "effective_weight_alt",
    "effective_weight_mixed", "ambient_a_ref", "ambient_a_alt", "ambient_a_mixed",
    "ambient_genotyped_mass", "qname_fallback_fraction",
    "mean_sites_per_molecule", "evidence_basis", "model_eligible",
    "evidence_status", "schema_version",
)

EXPRESSION_VALUE_FIELDS: Tuple[str, ...] = (
    "arm_counts", "other_autosomal_counts", "total_autosomal_counts",
    "arm_fraction", "log2_arm_to_other", "mapped_genes_on_arm",
    "nonzero_genes_on_arm",
    "reference_autosomal_counts", "log2_arm_to_reference",
    "matrix_value_type", "expression_input_state", "schema_version",
)

MANIFEST_AUDIT_FIELDS: Tuple[str, ...] = (
    "production_assignment", "species_a", "species_b", "demux_score",
    "production_assignment_source", "application_state", "application_reason",
    "nuclear_warning_reasons", "mitochondrial_resolution_status",
    "known_line_relationship", "uid_resolution_status", "metadata_event_status",
    "library_exchange_status", "event_identity_evidence_disposition",
    "event_review_scope", "cell_exception_reasons", "ambient_evaluation_status",
    "ambient_production_minus_original_c", "ambient_exact_donor_burden_fields",
    "ambient_background_shift_fields",
    "current_ploidy_state", "ploidy_evidence_status", "nn_prob_tetraploid",
    "nn_qc_pass", "occupancy_state", "occupancy_evidence_status",
    "technical_state", "nuclear_reconciliation_status",
    "species_evidence_status", "mitochondrial_evidence_status",
    "atac_evidence_status", "review_required", "review_reasons",
    "downstream_release_status", "downstream_exclusion_reason",
    "identity_changed", "event_id", "ambient_c", "ambient_c_se",
    "ambient_arm", "ambient_profile_mode", "ambient_profile_status",
    "model_eligible", "calibration_eligible", "eligibility_reasons",
    "schema_version",
)

MODEL_FIELDS: Tuple[str, ...] = (
    "calibration_group", "uid", "calibration_level", "calibration_status",
    "calibration_cells_exact", "calibration_cells_effective",
    "total_effective_ase_weight",
    "calibration_crossfit_fold", "calibration_excluded_chromosome",
    "group_genomewide_baseline_logit", "cell_loo_baseline_logit",
    "cell_loo_baseline_donor_a_fraction", "cell_loo_baseline_arms",
    "cell_loo_retained_observations", "cell_loo_baseline_status",
    "orientation_mapping_logit_offset_ref",
    "orientation_mapping_logit_offset_alt",
    "orientation_mapping_logit_offset_mixed", "quasi_overdispersion_rho_ref",
    "quasi_overdispersion_rho_alt", "quasi_overdispersion_rho_mixed",
    "expected_a_balanced_ref", "expected_a_balanced_alt",
    "expected_a_balanced_mixed", "expected_a_best_ref",
    "expected_a_best_alt", "expected_a_best_mixed", "expression_used",
    "expression_role", "expression_metric", "expression_calibration_level",
    "expression_calibration_cells", "expression_baseline_log2_arm_metric",
    "expression_sigma", "best_state", "best_state_posterior",
    "posterior_interpretation", "event_posterior", "posterior_BALANCED", "posterior_DONOR_A_LOSS",
    "posterior_DONOR_B_LOSS", "posterior_DONOR_A_GAIN",
    "posterior_DONOR_B_GAIN", "ase_log_bf_best_vs_balanced",
    "expression_log_bf_best_vs_balanced", "joint_log_bf_best_vs_balanced",
    "ambient_sensitivity_status",
    "empirical_test_basis", "empirical_null_level", "empirical_null_cells",
    "empirical_p_floor", "empirical_q_resolution_floor",
    "empirical_p_DONOR_A_LOSS", "empirical_q_DONOR_A_LOSS",
    "empirical_p_DONOR_B_LOSS", "empirical_q_DONOR_B_LOSS",
    "empirical_p_DONOR_A_GAIN", "empirical_q_DONOR_A_GAIN",
    "empirical_p_DONOR_B_GAIN", "empirical_q_DONOR_B_GAIN",
    "empirical_p_value", "empirical_q_value", "bh_scope", "call_state",
    "call_status", "qc_flags", "call_schema_version",
)

CALL_FIELDS: Tuple[str, ...] = (
    ASE_RAW_FIELDS
    + tuple(f"manifest_{name}" for name in MANIFEST_AUDIT_FIELDS)
    + tuple(f"expression_{name}" for name in EXPRESSION_VALUE_FIELDS)
    + MODEL_FIELDS
)

CALIBRATION_FIELDS: Tuple[str, ...] = (
    "library", "barcode", "uid", "calibration_group", "donor_pair",
    "chromosome", "arm", "crossfit_fold", "excluded_chromosome",
    "target_effective_ase_weight",
    "group_genomewide_baseline_logit", "cell_loo_baseline_logit",
    "cell_loo_baseline_donor_a_fraction", "cell_loo_baseline_arms",
    "cell_loo_retained_observations", "cell_loo_baseline_status",
    "exact_calibration_cells", "ref_source_level", "ref_source_key",
    "ref_calibration_cells", "ref_orientation_mapping_logit_offset",
    "ref_calibration_effective_weight",
    "ref_retained_calibration_cells",
    "ref_observed_donor_a_fraction", "ref_quasi_overdispersion_rho",
    "alt_source_level", "alt_source_key",
    "alt_calibration_cells", "alt_orientation_mapping_logit_offset",
    "alt_calibration_effective_weight",
    "alt_retained_calibration_cells",
    "alt_observed_donor_a_fraction", "alt_quasi_overdispersion_rho",
    "mixed_source_level", "mixed_source_key",
    "mixed_calibration_cells", "mixed_orientation_mapping_logit_offset",
    "mixed_calibration_effective_weight",
    "mixed_retained_calibration_cells",
    "mixed_observed_donor_a_fraction", "mixed_quasi_overdispersion_rho",
    "expression_source_level",
    "expression_source_key", "expression_calibration_cells",
    "expression_retained_calibration_cells", "expression_input_state",
    "expression_metric", "ambient_handling",
    "expression_baseline_log2_arm_metric", "expression_sigma",
    "calibration_status", "schema_version",
)

UID_FIELDS: Tuple[str, ...] = (
    "uid", "donor_pair", "chromosome", "libraries", "calibration_groups",
    "n_cells", "p_arm", "p_rows", "p_total_cells", "p_evaluable_cells",
    "p_qc_blocked_cells", "p_supporting_cells", "p_concordant_cells", "p_concordance",
    "p_state", "p_best_state_posterior",
    "p_event_posterior", "p_empirical_p_value", "p_empirical_q_value",
    "p_empirical_p_floor", "p_empirical_q_resolution_floor",
    "p_state_p_DONOR_A_LOSS", "p_state_q_DONOR_A_LOSS",
    "p_state_p_DONOR_B_LOSS", "p_state_q_DONOR_B_LOSS",
    "p_state_p_DONOR_A_GAIN", "p_state_q_DONOR_A_GAIN",
    "p_state_p_DONOR_B_GAIN", "p_state_q_DONOR_B_GAIN",
    "q_arm", "q_rows", "q_total_cells", "q_evaluable_cells",
    "q_qc_blocked_cells", "q_supporting_cells", "q_concordant_cells", "q_concordance",
    "q_state", "q_best_state_posterior",
    "q_event_posterior", "q_empirical_p_value", "q_empirical_q_value",
    "q_empirical_p_floor", "q_empirical_q_resolution_floor",
    "q_state_p_DONOR_A_LOSS", "q_state_q_DONOR_A_LOSS",
    "q_state_p_DONOR_B_LOSS", "q_state_q_DONOR_B_LOSS",
    "q_state_p_DONOR_A_GAIN", "q_state_q_DONOR_A_GAIN",
    "q_state_p_DONOR_B_GAIN", "q_state_q_DONOR_B_GAIN",
    "paired_pq_concordant_cells",
    "whole_chromosome_flag", "whole_chromosome_state", "summary_status",
    "schema_version",
)

PAIR_SUMMARY_FIELDS: Tuple[str, ...] = (
    "library", "calibration_group", "donor_pair", "arm", "chromosome",
    "libraries",
    "total_cells", "individual_evaluable_cells", "aggregate_eligible_cells",
    "eligible_uid_blocks", "tested_cells", "tested_uid_blocks",
    "qc_blocked_cells", "supporting_cells",
    "concordant_cells", "concordance", "best_state", "best_state_posterior",
    "event_posterior", "aggregate_log_bf_best_vs_balanced",
    "pooled_sites", "pooled_soft_molecule_units", "pooled_effective_units",
    "pooled_hard_informative_molecules",
    "partial_conjunction_r", "partial_conjunction_method", "fisher_terms",
    "fisher_df", "dependence_assumption", "fdr_interpretation",
    "partial_conjunction_p_value",
    "partial_conjunction_q_value", "partial_conjunction_p_floor",
    "partial_conjunction_q_resolution_floor",
    "state_p_DONOR_A_LOSS", "state_q_DONOR_A_LOSS",
    "state_p_DONOR_B_LOSS", "state_q_DONOR_B_LOSS",
    "state_p_DONOR_A_GAIN", "state_q_DONOR_A_GAIN",
    "state_p_DONOR_B_GAIN", "state_q_DONOR_B_GAIN",
    "recurrence_flag", "summary_status",
    "schema_version",
)

ASE_CALIBRATION_LEVELS: Tuple[str, ...] = (
    "LIBRARY_GROUP_DONOR_PAIR_GENOMEWIDE_LOCO",
    "LIBRARY_DONOR_PAIR_GENOMEWIDE_LOCO",
    "LIBRARY_GROUP_GENOMEWIDE_LOCO",
    "LIBRARY_GENOMEWIDE_LOCO",
    "GLOBAL_GENOMEWIDE_LOCO",
)

EXTERNAL_EXPRESSION_LEVELS: Tuple[str, ...] = (
    "LIBRARY_GROUP_ARM_EXTERNAL_PAIR",
    "LIBRARY_ARM_EXTERNAL_PAIR",
    "GLOBAL_ARM_EXTERNAL_PAIR",
)

EMPIRICAL_NULL_LEVELS: Tuple[str, ...] = (
    "LIBRARY_GROUP_ARM_EXTERNAL_PAIR_SAME_FOLD_DEPTH_AMBIENT_MATCHED",
    "LIBRARY_ARM_EXTERNAL_PAIR_SAME_FOLD_DEPTH_AMBIENT_MATCHED",
    "GLOBAL_ARM_EXTERNAL_PAIR_SAME_FOLD_DEPTH_AMBIENT_MATCHED",
)


class ExpressionStore:
    """Disk-backed expression lookup; bounds RAM for global multi-library calls."""

    def __init__(self, directory: str):
        fd, self.path = tempfile.mkstemp(
            prefix=".tetra_arm_expression.", suffix=".sqlite3", dir=directory)
        os.close(fd)
        self.connection = sqlite3.connect(self.path)
        self.connection.execute("PRAGMA journal_mode=OFF")
        self.connection.execute("PRAGMA synchronous=OFF")
        self.connection.execute("PRAGMA temp_store=MEMORY")
        value_columns = ", ".join(f'"{name}" TEXT NOT NULL' for name in EXPRESSION_VALUE_FIELDS)
        self.connection.execute(
            "CREATE TABLE expression ("
            "library TEXT NOT NULL, barcode TEXT NOT NULL, arm TEXT NOT NULL, "
            "chromosome TEXT NOT NULL, " + value_columns + ", "
            "calibration_group TEXT NOT NULL, donor_pair TEXT NOT NULL, "
            "uid TEXT NOT NULL, fold INTEGER NOT NULL, "
            "calibration_eligible INTEGER NOT NULL, model_eligible INTEGER NOT NULL, "
            "PRIMARY KEY (library, barcode, arm)) WITHOUT ROWID")
        self._insert_sql = (
            "INSERT INTO expression VALUES (" +
            ",".join("?" for _ in range(10 + len(EXPRESSION_VALUE_FIELDS))) + ")")
        self._pending: List[Tuple[str, ...]] = []
        self.count = 0
        self._finished = False

    def add(self, library: str, barcode: str, arm: str, chromosome: str,
            row: Mapping[str, str], manifest: Mapping[str, str], fold: int) -> None:
        if self._finished:
            raise RuntimeError("cannot add expression rows after index finalization")
        context_eligible = not manifest_context_review_reasons(manifest)
        self._pending.append((
            library, barcode, arm, chromosome,
            *(str(row.get(name, "")) for name in EXPRESSION_VALUE_FIELDS),
            clean(manifest.get("calibration_group")) or f"lib{library}",
            pair_from_row(manifest), clean(manifest.get("uid")), str(fold),
            str(int(truthy(manifest.get("calibration_eligible")) and
                    context_eligible)),
            str(int(truthy(manifest.get("model_eligible")) and
                    context_eligible))))
        if len(self._pending) >= 10000:
            self.flush()

    def flush(self) -> None:
        if not self._pending:
            return
        try:
            self.connection.executemany(self._insert_sql, self._pending)
        except sqlite3.IntegrityError as exc:
            raise ValueError("duplicate expression cell-arm key") from exc
        self.count += len(self._pending)
        self._pending.clear()

    def finish(self) -> None:
        if self._finished:
            return
        self.flush()
        self.connection.commit()
        self.connection.execute(
            "CREATE INDEX IF NOT EXISTS expression_calibration_idx ON expression "
            "(arm, library, calibration_group, donor_pair, fold)")
        self._finished = True

    @staticmethod
    def _row_to_mapping(values: Sequence[str]) -> Dict[str, str]:
        result = {name: str(value) for name, value in zip(EXPRESSION_VALUE_FIELDS, values[4:])}
        result.update({
            "library": str(values[0]), "barcode": str(values[1]),
            "arm": str(values[2]), "chromosome": str(values[3]),
        })
        return result

    def get(self, key: Tuple[str, str, str]) -> Optional[Dict[str, str]]:
        self.finish()
        row = self.connection.execute(
            "SELECT * FROM expression WHERE library=? AND barcode=? AND arm=?", key
        ).fetchone()
        return self._row_to_mapping(row) if row is not None else None

    def get_cell(self, library: str, barcode: str) -> Dict[str, Dict[str, str]]:
        self.finish()
        rows = self.connection.execute(
            "SELECT * FROM expression WHERE library=? AND barcode=?",
            (library, barcode)).fetchall()
        return {str(row[2]): self._row_to_mapping(row) for row in rows}

    def calibration_values(
            self, pool_key: Tuple[str, Tuple[str, ...]], query_pair: str,
            query_fold: int, input_state: str, metric_field: str,
            args) -> List[float]:
        """Return an external-pair, held-out expression reference stratum."""
        self.finish()
        level, key = pool_key
        clauses = [
            "donor_pair<>?", "fold<>?", "calibration_eligible=1",
            "model_eligible=1",
        ]
        parameters: List[object] = [query_pair, query_fold]
        if level == EXTERNAL_EXPRESSION_LEVELS[0]:
            library, group, arm = key
            clauses.extend(["library=?", "calibration_group=?", "arm=?"])
            parameters.extend([library, group, arm])
        elif level == EXTERNAL_EXPRESSION_LEVELS[1]:
            library, arm = key
            clauses.extend(["library=?", "arm=?"])
            parameters.extend([library, arm])
        elif level == EXTERNAL_EXPRESSION_LEVELS[2]:
            (arm,) = key
            clauses.append("arm=?")
            parameters.append(arm)
        else:
            raise ValueError(f"unknown external expression level: {level}")
        clauses.append('"expression_input_state"=?')
        parameters.append(input_state)
        if metric_field not in {"log2_arm_to_reference", "log2_arm_to_other"}:
            return []
        query = (
            f'SELECT "{metric_field}", '
            '"total_autosomal_counts", "nonzero_genes_on_arm", '
            '"matrix_value_type" FROM expression WHERE ' +
            " AND ".join(clauses))
        values = []
        for value_text, total, genes, matrix_type in self.connection.execute(
                query, parameters):
            if not expression_matrix_supported(matrix_type, input_state):
                continue
            value = finite_float(value_text)
            if (math.isfinite(value) and finite_float(total, 0.0) >=
                    args.min_expression_counts and finite_float(genes, 0.0) >=
                    args.min_expression_genes):
                values.append(value)
        return values

    def close(self) -> None:
        try:
            self.connection.close()
        finally:
            try:
                os.unlink(self.path)
            except FileNotFoundError:
                pass


def canonical_library(value: object) -> str:
    text = clean(value)
    stripped = text.lower().removeprefix("lib")
    if stripped.isdigit():
        return str(int(stripped))
    if not text:
        raise ValueError("empty library value")
    return text


def format_number(value: float) -> str:
    return "NA" if not math.isfinite(value) else f"{value:.17g}"


def bounded(value: float, lower: float, upper: float) -> float:
    return min(upper, max(lower, value))


def logistic(value: float) -> float:
    if value >= 0:
        term = math.exp(-value)
        return 1.0 / (1.0 + term)
    term = math.exp(value)
    return term / (1.0 + term)


def logit(value: float) -> float:
    value = bounded(value, 1e-8, 1.0 - 1e-8)
    return math.log(value / (1.0 - value))


def logsumexp(values: Sequence[float]) -> float:
    maximum = max(values)
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def parse_count(value: object, field_name: str, context: str) -> int:
    number = finite_float(value)
    if not math.isfinite(number) or number < 0 or abs(number - round(number)) > 1e-7:
        raise ValueError(f"{context}: invalid integer count {field_name}={value!r}")
    return int(round(number))


def parse_nonnegative_float(value: object, field_name: str, context: str) -> float:
    number = finite_float(value)
    if not math.isfinite(number) or number < 0.0:
        raise ValueError(f"{context}: invalid nonnegative value {field_name}={value!r}")
    return number


def read_rows(path: str, label: str) -> Tuple[List[str], Iterable[Dict[str, str]]]:
    """Return a stable header and a generator that validates row widths."""
    target = require_file(path, label)
    handle = open_text(target)
    reader = csv.DictReader(handle, delimiter="\t")
    header = list(reader.fieldnames or [])
    if not header or len(header) != len(set(header)):
        handle.close()
        raise ValueError(f"{label} is headerless or has duplicate columns: {target}")

    def generate():
        try:
            for line_number, row in enumerate(reader, start=2):
                if None in row:
                    raise ValueError(f"malformed row {target}:{line_number}")
                yield {str(key): str(value) for key, value in row.items()}
        finally:
            handle.close()

    return header, generate()


def require_columns(header: Sequence[str], required: Iterable[str], path: str) -> None:
    missing = sorted(set(required) - set(header))
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")


def pair_from_row(row: Mapping[str, str]) -> str:
    pair = clean(row.get("donor_pair"))
    if pair:
        return pair
    donor_a, donor_b = clean(row.get("donor_a")), clean(row.get("donor_b"))
    return f"{donor_a}+{donor_b}" if donor_a and donor_b else "NA"


def arm_side(arm: str) -> str:
    lowered = clean(arm).lower()
    if lowered.endswith("p"):
        return "p"
    if lowered.endswith("q"):
        return "q"
    return ""


def ase_calibration_keys(library: str, group: str, donor_pair: str
                         ) -> List[Tuple[str, Tuple[str, ...]]]:
    return [
        (ASE_CALIBRATION_LEVELS[0], (library, group, donor_pair)),
        (ASE_CALIBRATION_LEVELS[1], (library, donor_pair)),
        (ASE_CALIBRATION_LEVELS[2], (library, group)),
        (ASE_CALIBRATION_LEVELS[3], (library,)),
        (ASE_CALIBRATION_LEVELS[4], ()),
    ]


def external_expression_keys(library: str, group: str, arm: str
                             ) -> List[Tuple[str, Tuple[str, ...]]]:
    return [
        (EXTERNAL_EXPRESSION_LEVELS[0], (library, group, arm)),
        (EXTERNAL_EXPRESSION_LEVELS[1], (library, arm)),
        (EXTERNAL_EXPRESSION_LEVELS[2], (arm,)),
    ]


def empirical_null_keys(library: str, group: str, arm: str
                        ) -> List[Tuple[str, Tuple[str, ...]]]:
    return [
        (EMPIRICAL_NULL_LEVELS[0], (library, group, arm)),
        (EMPIRICAL_NULL_LEVELS[1], (library, arm)),
        (EMPIRICAL_NULL_LEVELS[2], (arm,)),
    ]


def chromosome_key(value: object) -> str:
    text = clean(value).upper()
    return text[3:] if text.startswith("CHR") else text


def is_autosomal(value: object) -> bool:
    text = chromosome_key(value)
    return text.isdigit() and 1 <= int(text) <= 22


def calibration_entity(row: Mapping[str, str]) -> str:
    uid = clean(row.get("uid"))
    if uid:
        return "UID:" + uid
    return "CELL:" + canonical_library(row.get("library")) + ":" + canonical_barcode(
        row.get("barcode", ""))


def crossfit_fold(row: Mapping[str, str], folds: int) -> int:
    digest = hashlib.sha256(calibration_entity(row).encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % folds


def key_text(key: Tuple[str, ...]) -> str:
    return "|".join(key)


@dataclass(slots=True)
class CallRecord:
    library: str
    barcode: str
    arm: str
    chromosome: str
    group: str
    donor_pair: str
    manifest: Mapping[str, str]
    ase: Mapping[str, str]
    expression: Optional[Mapping[str, str]]
    counts: Dict[str, Tuple[int, int]]
    soft_evidence: Dict[str, Tuple[float, float, int]]
    n_sites: int
    n_informative: int
    n_soft_units: int
    n_effective_units: float
    row_flags: List[str] = field(default_factory=list)
    model: Optional["ResolvedCalibration"] = None
    posteriors: Dict[str, float] = field(default_factory=dict)
    ase_relative: Dict[str, float] = field(default_factory=dict)
    expression_relative: Dict[str, float] = field(default_factory=dict)
    joint_relative: Dict[str, float] = field(default_factory=dict)
    expected_balanced: Dict[str, float] = field(default_factory=dict)
    expected_best: Dict[str, float] = field(default_factory=dict)
    expression_used: bool = False
    expression_value: float = math.nan
    expression_sigma_effective: float = math.nan
    best_state: str = "NO_CALL"
    best_posterior: float = math.nan
    event_posterior: float = math.nan
    ase_test_stat: float = math.nan
    joint_test_stat: float = math.nan
    ase_state_test_stats: Dict[str, float] = field(default_factory=dict)
    empirical_basis: str = "NA"
    empirical_level: str = "NA"
    empirical_null_cells: int = 0
    empirical_p_floor: float = math.nan
    empirical_q_resolution_floor: float = math.nan
    empirical_p: float = math.nan
    empirical_q: float = math.nan
    empirical_state_p: Dict[str, float] = field(default_factory=dict)
    empirical_state_q: Dict[str, float] = field(default_factory=dict)
    ambient_sensitivity_status: str = "NOT_EVALUATED"
    call_state: str = "NO_CALL"
    call_status: str = "UNSCORED"

    @property
    def key(self) -> Tuple[str, str, str]:
        return self.library, self.barcode, self.arm

    @property
    def exact_key(self) -> Tuple[str, ...]:
        return self.library, self.group, self.donor_pair, self.arm

    @property
    def uid(self) -> str:
        return clean(self.manifest.get("uid"))

    @property
    def calibration_eligible(self) -> bool:
        return truthy(self.manifest.get("calibration_eligible"))

    @property
    def model_eligible(self) -> bool:
        return truthy(self.manifest.get("model_eligible")) and truthy(
            self.ase.get("model_eligible"))


@dataclass(slots=True)
class BiasFit:
    n_cells: int = 0
    total_count: float = 0.0
    delta: float = 0.0
    rho: float = 0.02
    observed_fraction: float = math.nan
    retained_cells: int = 0

@dataclass(slots=True)
class ResolvedBias:
    fit: BiasFit
    source_level: str
    source_key: Tuple[str, ...]
    exact_cells: int
    data_calibrated: bool


@dataclass(slots=True)
class ExpressionFit:
    n_cells: int = 0
    center: float = 0.0
    sigma: float = 0.50
    retained_cells: int = 0


@dataclass(slots=True)
class ResolvedExpression:
    fit: ExpressionFit = field(default_factory=ExpressionFit)
    source_level: str = "NO_EXPRESSION_CALIBRATION"
    source_key: Tuple[str, ...] = ()
    exact_cells: int = 0
    data_calibrated: bool = False
    metric_field: str = "NA"


@dataclass(slots=True)
class CellBaselineFit:
    logit_value: float = 0.0
    n_arms: int = 0
    retained_observations: int = 0
    status: str = "GROUP_BASELINE_ONLY"


@dataclass(slots=True)
class ResolvedCalibration:
    biases: Dict[str, ResolvedBias]
    expression: ResolvedExpression
    exact_cells: int
    level: str
    status: str
    crossfit_fold: int
    excluded_chromosome: str
    group_baseline_logit: float
    cell_baseline_logit: float
    cell_baseline_arms: int
    cell_baseline_retained_observations: int
    cell_baseline_status: str


class CalibrationResolver:
    def __init__(self, records: Sequence[CallRecord],
                 expression_store: ExpressionStore, args):
        self.args = args
        self.record_by_key = {record.key: record for record in records}
        self.records_by_cell: MutableMapping[
            Tuple[str, str], List[CallRecord]
        ] = defaultdict(list)
        self.bias_pools: MutableMapping[
            Tuple[str, Tuple[str, ...]], List[CallRecord]
        ] = defaultdict(list)
        self.expression_store = expression_store
        for record in records:
            self.records_by_cell[(record.library, record.barcode)].append(record)
            if not self._usable_calibration_record(record):
                continue
            for pool_key in ase_calibration_keys(
                    record.library, record.group, record.donor_pair):
                self.bias_pools[pool_key].append(record)
        self.bias_fit_cache: Dict[
            Tuple[str, Tuple[str, ...], str, str, int], BiasFit
        ] = {}
        self.expression_fit_cache: Dict[
            Tuple[str, Tuple[str, ...], str, int, str, str], ExpressionFit
        ] = {}
        self.cell_baseline_cache: Dict[
            Tuple[str, str, str], CellBaselineFit
        ] = {}
        self.resolved_cache: Dict[
            Tuple[str, str, str], ResolvedCalibration
        ] = {}

    def _usable_calibration_record(self, record: CallRecord) -> bool:
        if not record.calibration_eligible or not record.model_eligible:
            return False
        if manifest_event_review_reasons(record):
            return False
        if clean(record.ase.get("evidence_status")) != "PASS":
            return False
        if "INFORMATIVE_COUNT_MISMATCH" in record.row_flags:
            return False
        if not is_autosomal(record.chromosome):
            return False
        if (finite_float(record.ase.get("qname_fallback_fraction"), 1.0) >
                self.args.max_qname_fallback_fraction):
            return False
        if (finite_float(record.ase.get("ambient_genotyped_mass"), 0.0) <
                self.args.min_ambient_genotyped_mass):
            return False
        if not ambient_uncertainty_available(record):
            return False
        return (record.n_sites >= self.args.min_calibration_sites and
                record.n_effective_units >= self.args.min_calibration_effective_weight and
                record.n_effective_units > 0.0)

    def _usable_cell_baseline_record(self, record: CallRecord) -> bool:
        if (clean(record.ase.get("evidence_status")) != "PASS" or
                "INFORMATIVE_COUNT_MISMATCH" in record.row_flags or
                not is_autosomal(record.chromosome)):
            return False
        if (finite_float(record.ase.get("qname_fallback_fraction"), 1.0) >
                self.args.max_qname_fallback_fraction or
                finite_float(record.ase.get("ambient_genotyped_mass"), 0.0) <
                self.args.min_ambient_genotyped_mass):
            return False
        if not ambient_uncertainty_available(record):
            return False
        return (record.n_sites >= self.args.min_calibration_sites and
                record.n_effective_units >=
                self.args.min_calibration_effective_weight)

    def _orientation_observations(
            self, records: Sequence[CallRecord], orientation: str,
            excluded_chromosome: str, excluded_fold: int,
    ) -> Tuple[List[Tuple[float, float, float]],
               List[Tuple[float, float, float]], int]:
        """Return cell-capped baseline and target-scale arm observations.

        The donor baseline must not let a cell with many observed arms dominate,
        so its first-stage observations are pooled once per cell.  In contrast,
        quasi-overdispersion is a property of a cell-arm likelihood and is
        estimated from the individual held-out arms.  Collapsing those arms to
        one genome-wide fraction would erase neutral arm-to-arm heterogeneity
        and systematically underestimate rho.
        """
        by_cell: MutableMapping[Tuple[str, str], List[float]] = defaultdict(
            lambda: [0.0, 0.0, 0.0])
        arm_observations: List[Tuple[float, float, float]] = []
        for record in records:
            if (chromosome_key(record.chromosome) == excluded_chromosome or
                    crossfit_fold(record.manifest,
                                  self.args.calibration_crossfit_folds) == excluded_fold):
                continue
            successes, total, _raw_units = record.soft_evidence[orientation]
            if total <= 0:
                continue
            base = ambient_adjusted_fraction(record, "BALANCED", orientation, 0.0)
            arm_observations.append((successes / total, total, base))
            cell = by_cell[(record.library, record.barcode)]
            cell[0] += successes
            cell[1] += total - successes
            cell[2] += total * base
        cell_observations = [
            (a_count / (a_count + b_count), a_count + b_count,
             expected_sum / (a_count + b_count))
            for a_count, b_count, expected_sum in by_cell.values()
            if a_count + b_count > 0
        ]
        return cell_observations, arm_observations, len(by_cell)

    def fit_bias(self, pool_key: Tuple[str, Tuple[str, ...]],
                 orientation: str, excluded_chromosome: str,
                 excluded_fold: int) -> BiasFit:
        cache_key = pool_key + (orientation, excluded_chromosome, excluded_fold)
        cached = self.bias_fit_cache.get(cache_key)
        if cached is not None:
            return cached
        observations, arm_observations, n_cells = self._orientation_observations(
            self.bias_pools.get(pool_key, []), orientation,
            excluded_chromosome, excluded_fold)
        result = fit_bias_component(
            observations, self.args, rho_observations=arm_observations,
            n_cells=n_cells)
        self.bias_fit_cache[cache_key] = result
        return result

    def fit_expression(self, pool_key: Tuple[str, Tuple[str, ...]],
                       query_pair: str, query_fold: int,
                       input_state: str, metric_field: str) -> ExpressionFit:
        cache_key = pool_key + (query_pair, query_fold, input_state, metric_field)
        cached = self.expression_fit_cache.get(cache_key)
        if cached is not None:
            return cached
        values = self.expression_store.calibration_values(
            pool_key, query_pair, query_fold, input_state, metric_field, self.args)
        result = fit_expression_component(values, self.args)
        self.expression_fit_cache[cache_key] = result
        return result

    def resolve_bias(self, keys: Sequence[Tuple[str, Tuple[str, ...]]],
                     orientation: str, excluded_chromosome: str,
                     excluded_fold: int) -> ResolvedBias:
        exact_key = keys[0]
        exact = self.fit_bias(
            exact_key, orientation, excluded_chromosome, excluded_fold)
        candidates: List[Tuple[Tuple[str, Tuple[str, ...]], BiasFit]] = []
        for pool_key in keys[1:]:
            fit = self.fit_bias(
                pool_key, orientation, excluded_chromosome, excluded_fold)
            if fit.n_cells:
                candidates.append((pool_key, fit))
        adequate = [item for item in candidates if self.bias_fit_adequate(
            item[1], self.args.min_fallback_calibration_cells)]
        if adequate:
            parent_key, parent = adequate[0]
        else:
            parent_key, parent = ("NEUTRAL_PRIOR", ()), BiasFit()

        if (self.bias_fit_adequate(exact, 2) and self.bias_fit_adequate(
                parent, self.args.min_fallback_calibration_cells)):
            weight = exact.n_cells / (
                exact.n_cells + self.args.calibration_shrinkage_cells)
            delta = weight * exact.delta + (1.0 - weight) * parent.delta
            rho = weight * exact.rho + (1.0 - weight) * parent.rho
            total = (weight * exact.total_count +
                     (1.0 - weight) * parent.total_count)
            observed = exact.observed_fraction
            level = (exact_key[0] if parent.n_cells == 0 else
                     f"{exact_key[0]}_SHRUNK_TO_{parent_key[0]}")
            fit = BiasFit(exact.n_cells, total, delta, rho, observed,
                          exact.retained_cells)
            return ResolvedBias(fit, level,
                                exact_key[1] if parent.n_cells == 0 else parent_key[1],
                                exact.n_cells, True)
        if self.bias_fit_adequate(
                parent, self.args.min_fallback_calibration_cells):
            return ResolvedBias(parent, f"FALLBACK_{parent_key[0]}",
                                parent_key[1], 0, True)
        return ResolvedBias(BiasFit(), "NEUTRAL_PRIOR", (), 0, False)

    def bias_fit_adequate(self, fit: BiasFit, minimum_cells: int) -> bool:
        return (fit.n_cells >= minimum_cells and
                fit.total_count >=
                self.args.min_orientation_calibration_effective_weight)

    def enforce_matched_ref_alt(
            self, biases: Dict[str, ResolvedBias],
            keys: Sequence[Tuple[str, Tuple[str, ...]]],
            excluded_chromosome: str, excluded_fold: int,
            ) -> Dict[str, ResolvedBias]:
        """Make ref/alt nuisance estimates come from one identical stratum."""
        ref, alt = biases["ref"], biases["alt"]
        if (ref.data_calibrated and alt.data_calibrated and
                ref.source_level == alt.source_level and
                ref.source_key == alt.source_key):
            return biases
        for pool_key in keys[1:]:
            fits = {
                orientation: self.fit_bias(
                    pool_key, orientation, excluded_chromosome, excluded_fold)
                for orientation in ("ref", "alt")
            }
            if all(self.bias_fit_adequate(
                    fit, self.args.min_fallback_calibration_cells)
                    for fit in fits.values()):
                for orientation, fit in fits.items():
                    biases[orientation] = ResolvedBias(
                        fit, f"FALLBACK_{pool_key[0]}", pool_key[1], 0, True)
                return biases
        biases["ref"] = ResolvedBias(
            BiasFit(), "NO_MATCHED_REF_ALT_CALIBRATION", (), 0, False)
        biases["alt"] = ResolvedBias(
            BiasFit(), "NO_MATCHED_REF_ALT_CALIBRATION", (), 0, False)
        return biases

    def resolve_expression(
            self, keys: Sequence[Tuple[str, Tuple[str, ...]]],
            query_pair: str, query_fold: int, input_state: str,
            metric_field: str,
    ) -> ResolvedExpression:
        exact_key = keys[0]
        if input_state not in {
                "OBSERVED_FILTERED_COUNTS", "UPSTREAM_AMBIENT_CORRECTED"}:
            return ResolvedExpression()
        exact = self.fit_expression(
            exact_key, query_pair, query_fold, input_state, metric_field)
        candidates: List[Tuple[Tuple[str, Tuple[str, ...]], ExpressionFit]] = []
        for pool_key in keys[1:]:
            fit = self.fit_expression(
                pool_key, query_pair, query_fold, input_state, metric_field)
            if fit.n_cells:
                candidates.append((pool_key, fit))
        adequate = [item for item in candidates
                    if item[1].n_cells >= self.args.min_fallback_calibration_cells]
        if adequate:
            parent_key, parent = adequate[0]
        else:
            parent_key, parent = ("NO_EXPRESSION_CALIBRATION", ()), ExpressionFit()

        if (exact.n_cells >= 2 and
                parent.n_cells >= self.args.min_fallback_calibration_cells):
            weight = exact.n_cells / (
                exact.n_cells + self.args.calibration_shrinkage_cells)
            center = weight * exact.center + (1.0 - weight) * parent.center
            sigma = math.sqrt(
                weight * exact.sigma ** 2 + (1.0 - weight) * parent.sigma ** 2)
            level = (exact_key[0] if parent.n_cells == 0 else
                     f"{exact_key[0]}_SHRUNK_TO_{parent_key[0]}")
            return ResolvedExpression(
                ExpressionFit(exact.n_cells, center, sigma,
                              exact.retained_cells), level,
                exact_key[1] if parent.n_cells == 0 else parent_key[1],
                exact.n_cells, True, metric_field)
        if parent.n_cells >= self.args.min_fallback_calibration_cells:
            return ResolvedExpression(parent, f"FALLBACK_{parent_key[0]}",
                                      parent_key[1], 0, True, metric_field)
        return ResolvedExpression(ExpressionFit(), data_calibrated=False)

    def fit_cell_baseline(
            self, record: CallRecord, biases: Mapping[str, ResolvedBias],
            group_baseline: float) -> CellBaselineFit:
        cache_key = (record.library, record.barcode, chromosome_key(record.chromosome))
        cached = self.cell_baseline_cache.get(cache_key)
        if cached is not None:
            return cached
        observations = []
        for other in self.records_by_cell[(record.library, record.barcode)]:
            if (not self._usable_cell_baseline_record(other) or
                    chromosome_key(other.chromosome) == chromosome_key(record.chromosome)):
                continue
            for orientation in ORIENTATIONS:
                successes, total, _raw_units = other.soft_evidence[orientation]
                if total:
                    observations.append((
                        other, orientation, successes / total, total,
                        biases[orientation].fit.delta))
        if not observations:
            result = CellBaselineFit(group_baseline, 0, 0, "GROUP_BASELINE_ONLY")
            self.cell_baseline_cache[cache_key] = result
            return result

        robust = [1.0] * len(observations)

        def fit_eta(weights: Sequence[float]) -> float:
            def score(eta: float) -> float:
                value = 0.0
                for index, (source, orientation, fraction, count, delta) in enumerate(
                        observations):
                    expected = ambient_adjusted_fraction(
                        source, "BALANCED", orientation, delta,
                        cellular_baseline_logit=eta)
                    value += min(count, self.args.calibration_max_cell_weight) * (
                        weights[index]) * (fraction - expected)
                return value
            lower, upper = -6.0, 6.0
            for _ in range(80):
                midpoint = (lower + upper) / 2.0
                if score(midpoint) > 0.0:
                    lower = midpoint
                else:
                    upper = midpoint
            return (lower + upper) / 2.0

        raw_eta = group_baseline
        for _ in range(4):
            raw_eta = fit_eta(robust)
            updated = []
            for source, orientation, fraction, count, delta in observations:
                expected = ambient_adjusted_fraction(
                    source, "BALANCED", orientation, delta,
                    cellular_baseline_logit=raw_eta)
                rho = biases[orientation].fit.rho
                variance = max(expected * (1.0 - expected) * (
                    rho + (1.0 - rho) / max(count, 1e-8)), 1e-6)
                z_score = abs(fraction - expected) / math.sqrt(variance)
                updated.append(min(
                    1.0, self.args.calibration_huber_z / max(z_score, 1e-12)))
            robust = updated
        retained = [weight >= self.args.calibration_min_robust_weight
                    for weight in robust]
        if sum(retained) >= 3:
            raw_eta = fit_eta([1.0 if keep else 0.0 for keep in retained])
        retained_arms = {
            source.arm for keep, (source, _orientation, _fraction, _count, _delta)
            in zip(retained, observations) if keep
        }
        n_arms = len(retained_arms)
        weight = n_arms / (n_arms + self.args.cell_baseline_shrinkage_arms)
        effective_eta = weight * raw_eta + (1.0 - weight) * group_baseline
        status = ("PASS_CELL_LOCO_BASELINE" if
                  n_arms >= self.args.min_cell_baseline_arms else
                  "SPARSE_CELL_LOCO_SHRUNK_TO_GROUP")
        result = CellBaselineFit(
            effective_eta, n_arms, sum(retained), status)
        self.cell_baseline_cache[cache_key] = result
        return result

    @staticmethod
    def decompose_orientation_biases(
            biases: Mapping[str, ResolvedBias]
            ) -> Tuple[Dict[str, ResolvedBias], float]:
        ref, alt, mixed = (biases[name] for name in ORIENTATIONS)
        candidates: List[Tuple[float, float]] = []
        matched_ref_alt = (
            ref.data_calibrated and alt.data_calibrated and
            ref.source_level == alt.source_level and
            ref.source_key == alt.source_key)
        if matched_ref_alt:
            candidates.append(((ref.fit.delta + alt.fit.delta) / 2.0,
                               min(ref.fit.total_count, alt.fit.total_count)))
        mixed_matches_pair = (
            matched_ref_alt and mixed.data_calibrated and
            mixed.source_level == ref.source_level and
            mixed.source_key == ref.source_key)
        if mixed_matches_pair or (not matched_ref_alt and mixed.data_calibrated):
            candidates.append((mixed.fit.delta,
                               mixed.fit.total_count))
        total_weight = sum(weight for _value, weight in candidates)
        group_baseline = (sum(value * weight for value, weight in candidates) /
                          total_weight if total_weight else 0.0)
        if matched_ref_alt:
            mapping = (ref.fit.delta - alt.fit.delta) / 2.0
        else:
            mapping = 0.0
        offsets = {"ref": mapping, "alt": -mapping, "mixed": 0.0}
        decomposed = {}
        for orientation, resolved in biases.items():
            fit = resolved.fit
            usable = (resolved.data_calibrated and
                      (orientation == "mixed" or matched_ref_alt))
            decomposed[orientation] = ResolvedBias(
                BiasFit(fit.n_cells, fit.total_count, offsets[orientation],
                        fit.rho, fit.observed_fraction, fit.retained_cells),
                resolved.source_level, resolved.source_key,
                resolved.exact_cells, usable)
        return decomposed, group_baseline

    def resolve(self, record: CallRecord) -> ResolvedCalibration:
        cached = self.resolved_cache.get(record.key)
        if cached is not None:
            return cached
        excluded_chromosome = chromosome_key(record.chromosome)
        held_out_fold = crossfit_fold(
            record.manifest, self.args.calibration_crossfit_folds)
        bias_keys = ase_calibration_keys(
            record.library, record.group, record.donor_pair)
        raw_biases = {
            orientation: self.resolve_bias(
                bias_keys, orientation, excluded_chromosome, held_out_fold)
            for orientation in ORIENTATIONS
        }
        raw_biases = self.enforce_matched_ref_alt(
            raw_biases, bias_keys, excluded_chromosome, held_out_fold)
        biases, group_baseline = self.decompose_orientation_biases(raw_biases)
        cell_baseline = self.fit_cell_baseline(record, biases, group_baseline)
        expression = self.resolve_expression(
            external_expression_keys(record.library, record.group, record.arm),
            record.donor_pair, held_out_fold,
            clean(record.expression.get("expression_input_state")).upper()
            if record.expression is not None else "",
            expression_metric_field(record.expression))
        exact_cells = max(value.exact_cells for value in biases.values())
        levels = sorted(set(value.source_level for value in biases.values()))
        used = [orientation for orientation in ORIENTATIONS
                if biases[orientation].data_calibrated]
        if len(used) == len(ORIENTATIONS):
            status = "PASS_EMPIRICAL_CALIBRATION"
        elif used:
            status = "PARTIAL_ORIENTATION_CALIBRATION"
        else:
            status = "NO_EMPIRICAL_CALIBRATION"
        result = ResolvedCalibration(
            biases, expression, exact_cells, ";".join(levels), status,
            held_out_fold, excluded_chromosome, group_baseline,
            cell_baseline.logit_value, cell_baseline.n_arms,
            cell_baseline.retained_observations, cell_baseline.status)
        self.resolved_cache[record.key] = result
        return result


def fit_delta(observations: Sequence[Tuple[float, float, float]],
              robust_weights: Optional[Sequence[float]], max_weight: float) -> float:
    if not observations:
        return 0.0

    def score(delta: float) -> float:
        total = 0.0
        for index, (fraction, count, base) in enumerate(observations):
            weight = min(float(count), max_weight)
            if robust_weights is not None:
                weight *= robust_weights[index]
            expected = logistic(logit(base) + delta)
            total += weight * (fraction - expected)
        return total

    lower, upper = -6.0, 6.0
    if score(lower) <= 0.0:
        return lower
    if score(upper) >= 0.0:
        return upper
    for _ in range(80):
        midpoint = (lower + upper) / 2.0
        if score(midpoint) > 0.0:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0


def fit_bias_component(
        observations: Sequence[Tuple[float, float, float]], args,
        rho_observations: Optional[Sequence[Tuple[float, float, float]]] = None,
        n_cells: Optional[int] = None,
        ) -> BiasFit:
    if not observations:
        return BiasFit()

    rho_scale = list(rho_observations) if rho_observations is not None else list(
        observations)
    if n_cells is None:
        n_cells = len(observations)

    def estimate_rho(current_delta: float, weights: Sequence[float],
                     subset: Sequence[Tuple[float, float, float]]) -> float:
        numerator = 0.0
        denominator = 0.0
        for weight, (fraction, count, base) in zip(weights, subset):
            if count <= 1.0:
                continue
            expected = logistic(logit(base) + current_delta)
            binomial = expected * (1.0 - expected) / count
            numerator += weight * ((fraction - expected) ** 2 - binomial)
            denominator += weight * expected * (1.0 - expected) * (
                1.0 - 1.0 / count)
        if denominator <= 0.0:
            return args.default_rho
        return bounded(max(0.0, numerator) / denominator,
                       args.min_rho, args.max_rho)

    baseline_robust = [1.0] * len(observations)
    rho_robust = [1.0] * len(rho_scale)
    delta = 0.0
    rho = args.default_rho
    for _ in range(5):
        delta = fit_delta(
            observations, baseline_robust, args.calibration_max_cell_weight)
        rho = estimate_rho(delta, rho_robust, rho_scale)
        updated_baseline = []
        for fraction, count, base in observations:
            expected = logistic(logit(base) + delta)
            variance = max(expected * (1.0 - expected) * (
                rho + (1.0 - rho) / max(count, 1e-8)), 1e-6)
            z_score = abs(fraction - expected) / math.sqrt(variance)
            updated_baseline.append(min(
                1.0, args.calibration_huber_z / max(z_score, 1e-12)))
        updated_rho = []
        for fraction, count, base in rho_scale:
            expected = logistic(logit(base) + delta)
            variance = max(expected * (1.0 - expected) * (
                rho + (1.0 - rho) / max(count, 1e-8)), 1e-6)
            z_score = abs(fraction - expected) / math.sqrt(variance)
            updated_rho.append(min(
                1.0, args.calibration_huber_z / max(z_score, 1e-12)))
        baseline_robust = updated_baseline
        rho_robust = updated_rho

    retained_mask = [weight >= args.calibration_min_robust_weight
                     for weight in baseline_robust]
    if sum(retained_mask) >= 3:
        retained_observations = [observation for observation, keep in zip(
            observations, retained_mask) if keep]
        delta = fit_delta(
            retained_observations, None, args.calibration_max_cell_weight)
        working = [(1.0, observation) for observation in retained_observations]
    else:
        working = list(zip(baseline_robust, observations))
    retained_rho_mask = [weight >= args.calibration_min_robust_weight
                         for weight in rho_robust]
    if sum(retained_rho_mask) >= 3:
        retained_rho = [observation for observation, keep in zip(
            rho_scale, retained_rho_mask) if keep]
        rho = estimate_rho(delta, [1.0] * len(retained_rho), retained_rho)
    else:
        rho = estimate_rho(delta, rho_robust, rho_scale)
    if len(rho_scale) < 3:
        rho = args.default_rho
    rho = bounded(rho, args.min_rho, args.max_rho)
    total_count = sum(
        weight * min(count, args.calibration_max_cell_weight)
        for weight, (_fraction, count, _base) in working)
    observed_weight = sum(count for _fraction, count, _base in observations)
    observed = (sum(fraction * count for fraction, count, _ in observations) /
                observed_weight if observed_weight else math.nan)
    retained_cells = sum(retained_mask)
    return BiasFit(n_cells, total_count, delta, rho, observed,
                   retained_cells)


def fit_expression_component(values: Sequence[float], args) -> ExpressionFit:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return ExpressionFit()
    initial_count = len(finite)
    center = statistics.median(finite)
    deviations = [abs(value - center) for value in finite]
    mad_sigma = 1.4826 * statistics.median(deviations)
    if len(finite) >= 4:
        quartiles = statistics.quantiles(finite, n=4, method="inclusive")
        iqr_sigma = (quartiles[2] - quartiles[0]) / 1.349
    else:
        iqr_sigma = 0.0
    preliminary_sigma = max(args.min_expression_sigma, mad_sigma, iqr_sigma)
    retained = [value for value in finite
                if abs(value - center) <= args.calibration_expression_mad_cutoff *
                preliminary_sigma]
    if len(retained) >= 3:
        center = statistics.median(retained)
        deviations = [abs(value - center) for value in retained]
        mad_sigma = 1.4826 * statistics.median(deviations)
        if len(retained) >= 4:
            quartiles = statistics.quantiles(retained, n=4, method="inclusive")
            iqr_sigma = (quartiles[2] - quartiles[0]) / 1.349
    else:
        retained = finite
    sigma = max(args.min_expression_sigma, mad_sigma, iqr_sigma)
    sigma = min(args.max_expression_sigma, sigma)
    return ExpressionFit(initial_count, center, sigma, len(retained))


def ambient_a(record: CallRecord, orientation: str) -> float:
    field_name = ORIENTATION_FIELDS[orientation][2]
    conditional = finite_float(record.ase.get(field_name), 0.5)
    conditional = bounded(conditional, 0.0, 1.0)
    # Producer contract: ambient_a_* is already over the complete source
    # simplex, including 0.5 for untyped mass.  ambient_genotyped_mass is QC
    # context only; multiplying by it here would shrink the mixture twice.
    return conditional


def ambient_uncertainty_available(record: CallRecord) -> bool:
    """Production inference requires an estimable, positive ambient-c SE."""
    standard_error = finite_float(record.ase.get("ambient_c_se"))
    return math.isfinite(standard_error) and standard_error > 0.0


def ambient_adjusted_fraction(record: CallRecord, state: str,
                              orientation: str, mapping_delta: float,
                              contamination_override: Optional[float] = None,
                              cellular_baseline_logit: float = 0.0) -> float:
    contamination = bounded(
        finite_float(record.ase.get("ambient_c"), 0.0)
        if contamination_override is None else contamination_override,
        0.0, 0.999999)
    cellular = logistic(
        cellular_baseline_logit + STATE_LOG_ODDS_SHIFT[state])
    mixed = (1.0 - contamination) * cellular + contamination * ambient_a(
        record, orientation)
    return bounded(logistic(logit(mixed) + mapping_delta), 1e-8, 1.0 - 1e-8)


def quasi_ase_log_likelihood(
        effective_a: float, effective_weight: float,
        probability: float, rho: float) -> float:
    """Student-t working likelihood for a confidence-weighted ASE mean.

    This is deliberately a working likelihood for a confidence-weighted mean,
    not an exact fractional-count probability mass function. Held-out empirical
    p/q values, not this pseudo-likelihood alone, control production calls.
    """
    if effective_weight <= 0.0:
        return 0.0
    fraction = bounded(effective_a / effective_weight, 0.0, 1.0)
    probability = bounded(probability, 1e-8, 1.0 - 1e-8)
    variance = probability * (1.0 - probability) * (
        bounded(rho, 0.0, 1.0) +
        (1.0 - bounded(rho, 0.0, 1.0)) / effective_weight)
    return student_t_log_density(
        fraction, probability, math.sqrt(max(variance, 1e-8)))


def expression_matrix_supported(value: object, input_state: object = "") -> bool:
    normalized = clean(value).upper()
    state = clean(input_state).upper()
    raw_types = {
        "INTEGER_COUNTS", "RAW_INTEGER_COUNTS", "RAW_COUNTS", "UMI_COUNTS",
    }
    corrected_types = {
        "AMBIENT_CORRECTED_COUNTS", "DECONTAMINATED_COUNTS",
        "FRACTIONAL_CORRECTED_COUNTS",
    }
    if state == "OBSERVED_FILTERED_COUNTS":
        return normalized in raw_types
    if state == "UPSTREAM_AMBIENT_CORRECTED":
        return normalized in corrected_types
    return False


def expression_metric_field(expression: Optional[Mapping[str, str]]) -> str:
    if expression is None:
        return "NA"
    if math.isfinite(finite_float(expression.get("log2_arm_to_reference"))):
        return "log2_arm_to_reference"
    if math.isfinite(finite_float(expression.get("log2_arm_to_other"))):
        return "log2_arm_to_other"
    return "NA"


def expression_is_usable(record: CallRecord, model: ResolvedCalibration,
                         args) -> bool:
    if record.expression is None or not model.expression.data_calibrated:
        return False
    input_state = clean(record.expression.get("expression_input_state")).upper()
    if input_state not in {
            "OBSERVED_FILTERED_COUNTS", "UPSTREAM_AMBIENT_CORRECTED"}:
        return False
    if not expression_matrix_supported(
            record.expression.get("matrix_value_type"), input_state):
        return False
    value = finite_float(record.expression.get(model.expression.metric_field))
    total = finite_float(record.expression.get("total_autosomal_counts"), 0.0)
    genes = finite_float(record.expression.get("nonzero_genes_on_arm"), 0.0)
    return (math.isfinite(value) and total >= args.min_expression_counts and
            genes >= args.min_expression_genes)


def expression_shift(record: CallRecord, state: str,
                     contamination_override: Optional[float] = None) -> float:
    input_state = clean(
        record.expression.get("expression_input_state")
        if record.expression is not None else "").upper()
    if input_state == "UPSTREAM_AMBIENT_CORRECTED":
        return math.log2(TOTAL_COPY_RATIO[state])
    contamination = bounded(
        (finite_float(record.manifest.get("ambient_c"),
                      finite_float(record.ase.get("ambient_c"), 0.0))
         if contamination_override is None else contamination_override),
        0.0, 0.999999)
    # Ambient is assumed neutral only for attenuation.  Without a matching
    # gene-level ambient profile this is not an exact expression correction.
    observed_ratio = contamination + (1.0 - contamination) * TOTAL_COPY_RATIO[state]
    return math.log2(max(observed_ratio, 1e-8))


def student_t_log_density(value: float, mean: float, sigma: float,
                          degrees_freedom: float = 4.0) -> float:
    sigma = max(sigma, 1e-6)
    standardized = (value - mean) / sigma
    # State-invariant normalizing constants can be omitted.
    return -math.log(sigma) - 0.5 * (degrees_freedom + 1.0) * math.log1p(
        standardized * standardized / degrees_freedom)


def contamination_quadrature(record: CallRecord) -> List[Tuple[float, float]]:
    """Equal-probability quadrature over the actual bounded normal law."""
    center = bounded(finite_float(record.ase.get("ambient_c"), 0.0), 0.0, 0.999999)
    standard_error = finite_float(record.ase.get("ambient_c_se"))
    if not math.isfinite(standard_error) or standard_error <= 0.0:
        return [(center, 1.0)]
    normal = statistics.NormalDist()
    lower_probability = normal.cdf((0.0 - center) / standard_error)
    upper_probability = normal.cdf((1.0 - center) / standard_error)
    strata = 15
    interval = upper_probability - lower_probability
    if interval <= 1e-15:
        return [(center, 1.0)]
    result = []
    for index in range(strata):
        probability = lower_probability + (index + 0.5) * interval / strata
        probability = bounded(probability, 1e-15, 1.0 - 1e-15)
        value = center + standard_error * normal.inv_cdf(probability)
        result.append((bounded(value, 0.0, 0.999999), 1.0 / strata))
    return result


def expression_sigma_at_contamination(
        record: CallRecord, model: ResolvedCalibration,
        contamination: float) -> float:
    input_state = clean(
        record.expression.get("expression_input_state")
        if record.expression is not None else "").upper()
    if input_state == "UPSTREAM_AMBIENT_CORRECTED":
        return model.expression.fit.sigma
    sigma = model.expression.fit.sigma * (1.0 + contamination)
    standard_error = finite_float(record.ase.get("ambient_c_se"))
    if contamination > 0.0 and (not math.isfinite(standard_error) or standard_error <= 0.0):
        sigma *= 1.5
    return sigma


def posterior_at_contamination(record: CallRecord, model: ResolvedCalibration,
                               args, contamination: float) -> Tuple[str, float]:
    ase_values = {}
    for state in STATES:
        value = 0.0
        for orientation in ORIENTATIONS:
            a_count, total, _raw_units = record.soft_evidence[orientation]
            expected = ambient_adjusted_fraction(
                record, state, orientation, model.biases[orientation].fit.delta,
                contamination,
                cellular_baseline_logit=model.cell_baseline_logit)
            value += quasi_ase_log_likelihood(
                a_count, total, expected,
                model.biases[orientation].fit.rho)
        ase_values[state] = value
    evidence_weight = (args.site_fallback_likelihood_weight
                       if clean(record.ase.get("evidence_status")) == "PASS_SITE_FALLBACK"
                       else 1.0)
    relatives = {
        state: evidence_weight * (ase_values[state] - ase_values["BALANCED"])
        for state in STATES
    }
    if record.expression_used:
        expr_ll = {
            state: student_t_log_density(
                record.expression_value,
                model.expression.fit.center + expression_shift(
                    record, state, contamination),
                expression_sigma_at_contamination(record, model, contamination))
            for state in STATES
        }
        for state in STATES:
            expr_relative = bounded(
                expr_ll[state] - expr_ll["BALANCED"],
                -args.max_expression_log_bf, args.max_expression_log_bf)
            relatives[state] += args.expression_weight * expr_relative
    priors = {"BALANCED": 1.0 - args.event_prior}
    priors.update({state: args.event_prior / (len(STATES) - 1)
                   for state in STATES[1:]})
    logs = {state: math.log(priors[state]) + relatives[state] for state in STATES}
    normalizer = logsumexp(list(logs.values()))
    posterior = {state: math.exp(logs[state] - normalizer) for state in STATES}
    best = max(STATES, key=lambda state: posterior[state])
    return best, posterior[best]


def posterior_conclusion(state: str, posterior: float, args) -> str:
    threshold = (args.min_balanced_posterior if state == "BALANCED"
                 else args.min_event_posterior)
    return state if posterior >= threshold else "NO_CALL"


def score_record(record: CallRecord, resolver: CalibrationResolver, args) -> None:
    model = resolver.resolve(record)
    record.model = model
    contamination = bounded(
        finite_float(record.ase.get("ambient_c"), 0.0), 0.0, 0.999999)
    contamination_se = finite_float(record.ase.get("ambient_c_se"))
    if not math.isfinite(contamination_se) or contamination_se <= 0.0:
        record.ambient_sensitivity_status = "POINT_ESTIMATE_ONLY_MISSING_OR_ZERO_SE"
        record.row_flags.append("AMBIENT_C_SE_UNAVAILABLE")
    if not is_autosomal(record.chromosome):
        record.row_flags.append("NON_AUTOSOMAL_EXPLORATORY")
    if not record.model_eligible:
        record.call_status = "NOT_MODEL_ELIGIBLE"
        record.row_flags.append("NOT_MODEL_ELIGIBLE")
        return
    if "INFORMATIVE_COUNT_MISMATCH" in record.row_flags:
        record.call_status = "INVALID_INFORMATIVE_COUNT_SUM"
        return
    evidence_status = clean(record.ase.get("evidence_status"))
    if not evidence_status.startswith("PASS"):
        record.call_status = "ASE_EVIDENCE_NOT_PASS"
        record.row_flags.append(evidence_status or "ASE_EVIDENCE_NOT_PASS")
        return
    if record.n_sites < args.min_call_sites:
        record.row_flags.append("LOW_SITE_COUNT")
    if record.n_sites < 1:
        record.call_status = "NO_USABLE_SITES"
        return
    if record.n_effective_units < args.min_call_effective_weight:
        record.row_flags.append("LOW_EFFECTIVE_DIRECTIONAL_EXPOSURE")
    if record.n_effective_units <= 0.0:
        record.call_status = "NO_DIRECTIONAL_SOFT_EVIDENCE"
        record.row_flags.append("NO_DIRECTIONAL_SOFT_EVIDENCE")
        return

    for orientation in ORIENTATIONS:
        if (record.soft_evidence[orientation][1] > 0 and
                not model.biases[orientation].data_calibrated):
            record.row_flags.append(f"NO_{orientation.upper()}_EMPIRICAL_CALIBRATION")
    if any(flag.startswith("NO_") and flag.endswith("_EMPIRICAL_CALIBRATION")
           for flag in record.row_flags):
        record.call_status = "NO_EMPIRICAL_CALIBRATION"
        return

    evidence_weight = 1.0
    if evidence_status == "PASS_SITE_FALLBACK":
        evidence_weight = args.site_fallback_likelihood_weight
        record.row_flags.append("SITE_LEVEL_FALLBACK_DOWNWEIGHTED")

    record.expression_used = expression_is_usable(record, model, args)
    if record.expression_used:
        record.expression_value = finite_float(
            record.expression.get(
                model.expression.metric_field))  # type: ignore[union-attr]
        if model.expression.metric_field == "log2_arm_to_other":
            record.row_flags.append("EXPRESSION_LEGACY_SISTER_ARM_DENOMINATOR")
        point_contamination = bounded(
            finite_float(record.ase.get("ambient_c"), 0.0), 0.0, 0.999999)
        record.expression_sigma_effective = expression_sigma_at_contamination(
            record, model, point_contamination)
    elif record.expression is None:
        record.row_flags.append("EXPRESSION_UNAVAILABLE")
    elif not expression_matrix_supported(
            record.expression.get("matrix_value_type"),
            record.expression.get("expression_input_state")):
        record.row_flags.append("EXPRESSION_UNSUPPORTED_MATRIX_VALUE_TYPE")
    elif clean(record.expression.get("expression_input_state")).upper() not in {
            "OBSERVED_FILTERED_COUNTS", "UPSTREAM_AMBIENT_CORRECTED"}:
        record.row_flags.append("EXPRESSION_INPUT_STATE_UNSUPPORTED")
    elif not model.expression.data_calibrated:
        record.row_flags.append("NO_EXTERNAL_PAIR_EXPRESSION_CALIBRATION")
    else:
        record.row_flags.append("EXPRESSION_BELOW_QC")

    nodes = contamination_quadrature(record)
    ase_integrands: Dict[str, List[float]] = {state: [] for state in STATES}
    expression_integrands: Dict[str, List[float]] = {state: [] for state in STATES}
    joint_integrands: Dict[str, List[float]] = {state: [] for state in STATES}
    for contamination, node_weight in nodes:
        node_ase: Dict[str, float] = {}
        node_expression: Dict[str, float] = {}
        for state in STATES:
            value = 0.0
            for orientation in ORIENTATIONS:
                a_count, total, _raw_units = record.soft_evidence[orientation]
                expected = ambient_adjusted_fraction(
                    record, state, orientation,
                    model.biases[orientation].fit.delta, contamination,
                    cellular_baseline_logit=model.cell_baseline_logit)
                value += quasi_ase_log_likelihood(
                    a_count, total, expected,
                    model.biases[orientation].fit.rho)
            node_ase[state] = value
            if record.expression_used:
                node_expression[state] = student_t_log_density(
                    record.expression_value,
                    model.expression.fit.center + expression_shift(
                        record, state, contamination),
                    expression_sigma_at_contamination(
                        record, model, contamination))
            else:
                node_expression[state] = 0.0
        for state in STATES:
            log_weight = math.log(node_weight)
            expression_delta = bounded(
                node_expression[state] - node_expression["BALANCED"],
                -args.max_expression_log_bf, args.max_expression_log_bf)
            expression_value = node_expression["BALANCED"] + expression_delta
            ase_integrands[state].append(
                log_weight + evidence_weight * node_ase[state])
            expression_integrands[state].append(log_weight + expression_value)
            joint_integrands[state].append(
                log_weight + evidence_weight * node_ase[state] +
                args.expression_weight * expression_value)

    ase_ll = {state: logsumexp(values)
              for state, values in ase_integrands.items()}
    expression_ll = {state: logsumexp(values)
                     for state, values in expression_integrands.items()}
    joint_ll = {state: logsumexp(values)
                for state, values in joint_integrands.items()}
    record.ase_relative = {
        state: ase_ll[state] - ase_ll["BALANCED"] for state in STATES
    }
    record.expression_relative = {
        state: expression_ll[state] - expression_ll["BALANCED"]
        for state in STATES
    }
    record.joint_relative = {
        state: joint_ll[state] - joint_ll["BALANCED"] for state in STATES
    }

    expected_by_state: Dict[str, Dict[str, float]] = {}
    for state in STATES:
        expected_by_state[state] = {}
        for orientation in ORIENTATIONS:
            expected = ambient_adjusted_fraction(
                record, state, orientation, model.biases[orientation].fit.delta,
                cellular_baseline_logit=model.cell_baseline_logit)
            expected_by_state[state][orientation] = expected
    record.expected_balanced = expected_by_state["BALANCED"]
    priors = {"BALANCED": 1.0 - args.event_prior}
    for state in STATES[1:]:
        priors[state] = args.event_prior / (len(STATES) - 1)
    log_scores = {
        state: math.log(priors[state]) + record.joint_relative[state]
        for state in STATES
    }
    normalizer = logsumexp(list(log_scores.values()))
    record.posteriors = {
        state: math.exp(log_scores[state] - normalizer) for state in STATES
    }
    record.best_state = max(STATES, key=lambda state: record.posteriors[state])
    record.best_posterior = record.posteriors[record.best_state]
    record.event_posterior = 1.0 - record.posteriors["BALANCED"]
    record.expected_best = expected_by_state[record.best_state]
    record.ase_state_test_stats = {
        state: max(0.0, record.ase_relative[state])
        for state in STATES[1:]
    }
    record.ase_test_stat = max(
        record.ase_state_test_stats.values())
    record.joint_test_stat = max(
        0.0, max(record.joint_relative[state] for state in STATES[1:]))
    record.call_status = "SCORED_PENDING_EMPIRICAL_NULL"

    qname_fraction = finite_float(record.ase.get("qname_fallback_fraction"), 0.0)
    if qname_fraction > args.max_qname_fallback_fraction:
        record.row_flags.append("HIGH_QNAME_FALLBACK_FRACTION")
    genotyped_mass = finite_float(record.ase.get("ambient_genotyped_mass"), 0.0)
    if genotyped_mass < args.min_ambient_genotyped_mass:
        record.row_flags.append("LOW_AMBIENT_GENOTYPED_MASS")

    if math.isfinite(contamination_se) and contamination_se > 0.0:
        lower = bounded(contamination - 1.96 * contamination_se, 0.0, 0.999999)
        upper = bounded(contamination + 1.96 * contamination_se, 0.0, 0.999999)
        low_state, low_posterior = posterior_at_contamination(
            record, model, args, lower)
        high_state, high_posterior = posterior_at_contamination(
            record, model, args, upper)
        point_conclusion = posterior_conclusion(
            record.best_state, record.best_posterior, args)
        conclusions = {
            point_conclusion,
            posterior_conclusion(low_state, low_posterior, args),
            posterior_conclusion(high_state, high_posterior, args),
        }
        if len(conclusions) == 1:
            record.ambient_sensitivity_status = (
                "MARGINALIZED_TRUNCATED_NORMAL_STABLE_ACROSS_95_PERCENT_C_INTERVAL")
        else:
            record.ambient_sensitivity_status = (
                "MARGINALIZED_TRUNCATED_NORMAL_SENSITIVE_ACROSS_95_PERCENT_C_INTERVAL")
            record.row_flags.append("AMBIENT_C_SENSITIVE")


def empirical_p_value(statistic: float, null_values: Sequence[float]) -> float:
    if not math.isfinite(statistic) or not null_values:
        return math.nan
    more_extreme = sum(value >= statistic - 1e-12 for value in null_values)
    return (more_extreme + 1.0) / (len(null_values) + 1.0)


def benjamini_hochberg(values: Sequence[float]) -> List[float]:
    result = [math.nan] * len(values)
    finite_indices = [index for index, value in enumerate(values) if math.isfinite(value)]
    ordered = sorted(finite_indices, key=lambda index: values[index])
    running = 1.0
    total = len(ordered)
    for reverse_rank, index in enumerate(reversed(ordered), start=1):
        rank = total - reverse_rank + 1
        candidate = min(1.0, values[index] * total / rank)
        running = min(running, candidate)
        result[index] = running
    return result


def virtual_bh_q_bound(values: Sequence[float], replacement: float) -> float:
    """BH q for a virtual attainable-floor hypothesis, conservatively appended."""
    finite = sorted(value for value in values if math.isfinite(value))
    if not math.isfinite(replacement):
        return math.nan
    total = len(finite) + 1
    insertion = bisect.bisect_right(finite, replacement)
    best = total * replacement / (insertion + 1)
    for index in range(insertion, len(finite)):
        best = min(best, total * finite[index] / (index + 2))
    return min(1.0, best)


def partial_conjunction_p(values: Sequence[float], required: int) -> float:
    """Bonferroni partial conjunction valid without independence (UID use)."""
    finite = sorted(value for value in values if math.isfinite(value))
    if required < 1 or len(finite) < required:
        return math.nan
    return min(1.0, (len(finite) - required + 1) * finite[required - 1])


def fisher_partial_conjunction_p(
        values: Sequence[float], required: int) -> float:
    """Fisher test of at least ``required`` non-null inputs.

    The smallest ``required - 1`` p-values are discarded and standard Fisher
    combination is applied to p_(required), ..., p_(n).  Validity requires
    conditional independence; pair-summary output labels that assumption and
    its BH results as a working, approximate FDR rather than an exact guarantee.
    """
    finite = sorted(value for value in values
                    if math.isfinite(value) and 0.0 <= value <= 1.0)
    if required < 1 or len(finite) < required:
        return math.nan
    combined = [max(value, 1e-300) for value in finite[required - 1:]]
    terms = len(combined)
    gamma_argument = -sum(math.log(value) for value in combined)
    if gamma_argument <= 0.0:
        return 1.0
    log_argument = math.log(gamma_argument)
    log_terms = [
        index * log_argument - math.lgamma(index + 1.0)
        for index in range(terms)
    ]
    log_survival = -gamma_argument + logsumexp(log_terms)
    return bounded(math.exp(max(log_survival, -745.0)), 0.0, 1.0)


def pair_test_records(
        records: Sequence[CallRecord], args) -> Tuple[List[CallRecord], int]:
    """Apply deterministic score-blind per-UID and total pair-test caps."""
    blocks: MutableMapping[str, List[CallRecord]] = defaultdict(list)
    for record in records:
        block = record.uid or f"CELL:{record.library}:{record.barcode}"
        blocks[block].append(record)

    def digest(text: str) -> bytes:
        return hashlib.sha256(text.encode("utf-8")).digest()

    ordered_blocks: List[List[CallRecord]] = []
    for block, values in sorted(blocks.items(), key=lambda item: digest(item[0])):
        ordered = sorted(
            values,
            key=lambda record: digest(
                f"{record.library}\0{record.barcode}\0{record.arm}"))
        ordered_blocks.append(ordered[:args.max_pair_cells_per_uid])

    selected: List[CallRecord] = []
    depth = 0
    while len(selected) < args.max_pair_test_cells:
        added = False
        for block in ordered_blocks:
            if depth < len(block):
                selected.append(block[depth])
                added = True
                if len(selected) >= args.max_pair_test_cells:
                    break
        if not added:
            break
        depth += 1
    tested_blocks = len({
        record.uid or f"CELL:{record.library}:{record.barcode}"
        for record in selected
    })
    return selected, tested_blocks


def has_blocking_status(value: object) -> bool:
    """Recognize explicit unresolved/blocking audit states, not missing legacy fields."""
    text = clean(value).upper()
    if text in {
            "NOT_APPLICABLE", "NOT_REQUESTED", "NOT_DERIVED", "NO_EVENT",
            "NO_IMMEDIATE_REVIEW", "PASS", "RELEASED", "RESOLVED", "COMPLETE"}:
        return False
    return bool(text) and any(token in text for token in (
        "BLOCK", "FAIL", "REJECT", "EXCLUDE", "CONFLICT", "AMBIG",
        "UNRESOLVED", "REVIEW", "HOLD", "INVALID", "INELIGIBLE",
        "PENDING", "UNKNOWN", "MISSING", "NOT_RELEASE", "NOT_READY",
        "DEFER",
    ))


def has_meaningful_reason(value: object) -> bool:
    """Treat documented no-review sentinels as empty, not as reasons."""
    text = clean(value).upper()
    return bool(text) and text not in {
        "NOT_APPLICABLE", "NOT_REQUESTED", "NOT_DERIVED", "NO_EVENT",
        "NO_IMMEDIATE_REVIEW",
    }


def manifest_context_review_reasons(
        manifest: Mapping[str, str]) -> List[str]:
    """Return audit-only context gates without using them as likelihood terms."""
    reasons: List[str] = []
    identity_fields = (
        "application_state", "uid_resolution_status", "metadata_event_status",
        "library_exchange_status", "event_identity_evidence_disposition",
        "nuclear_reconciliation_status", "mitochondrial_resolution_status",
        "downstream_release_status",
    )
    if (truthy(manifest.get("review_required")) or
            has_meaningful_reason(manifest.get("downstream_exclusion_reason")) or
            has_meaningful_reason(manifest.get("cell_exception_reasons")) or
            any(has_blocking_status(manifest.get(name)) for name in identity_fields)):
        reasons.append("IDENTITY_CONTEXT")

    technical_fields = ("technical_state", "occupancy_evidence_status")
    if any(has_blocking_status(manifest.get(name)) for name in technical_fields):
        reasons.append("TECHNICAL_CONTEXT")

    ploidy = clean(manifest.get("current_ploidy_state")).upper()
    ploidy_supported = (
        ploidy in {"4", "4.0"} or
        (not ploidy.startswith("NOT_") and not has_blocking_status(ploidy) and any(
            token in ploidy for token in ("TETRA", "4N", "2:2", "2+2"))))
    ploidy_unsupported = bool(ploidy) and not ploidy_supported
    nn_qc = clean(manifest.get("nn_qc_pass"))
    if (ploidy_unsupported or has_blocking_status(
            manifest.get("ploidy_evidence_status")) or
            (bool(nn_qc) and not truthy(nn_qc))):
        reasons.append("PLOIDY_CONTEXT")
    return reasons


def manifest_event_review_reasons(record: CallRecord) -> List[str]:
    """Apply manifest context gates and persist their auditable row flags."""
    reasons = manifest_context_review_reasons(record.manifest)
    flags = {
        "IDENTITY_CONTEXT": "MANIFEST_IDENTITY_CONTEXT_REVIEW",
        "TECHNICAL_CONTEXT": "MANIFEST_TECHNICAL_CONTEXT_REVIEW",
        "PLOIDY_CONTEXT": "TETRAPLOID_BASELINE_UNSUPPORTED",
    }
    for reason in reasons:
        flag = flags[reason]
        if flag not in record.row_flags:
            record.row_flags.append(flag)
    return reasons


def assign_empirical_values(records: Sequence[CallRecord], args) -> None:
    pools: MutableMapping[
        Tuple[str, Tuple[str, ...]], List[CallRecord]
    ] = defaultdict(list)
    for record in records:
        if (not record.calibration_eligible or
                record.call_status != "SCORED_PENDING_EMPIRICAL_NULL"):
            continue
        if (manifest_event_review_reasons(record) or
                clean(record.ase.get("evidence_status")) != "PASS" or
                not ambient_uncertainty_available(record) or
                any(flag in record.row_flags for flag in (
                    "INFORMATIVE_COUNT_MISMATCH",
                    "HIGH_QNAME_FALLBACK_FRACTION",
                    "LOW_AMBIENT_GENOTYPED_MASS",
                    "SITE_LEVEL_FALLBACK_DOWNWEIGHTED",
                    "AMBIENT_C_SENSITIVE",
                    "NON_AUTOSOMAL_EXPLORATORY"))):
            continue
        for pool_key in empirical_null_keys(
                record.library, record.group, record.arm):
            pools[pool_key].append(record)

    for record in records:
        if record.call_status != "SCORED_PENDING_EMPIRICAL_NULL":
            continue
        query_fold = crossfit_fold(
            record.manifest, args.calibration_crossfit_folds)
        query_ambient = finite_float(record.ase.get("ambient_c"), 0.0)
        lower_depth = record.n_effective_units / args.empirical_depth_fold
        upper_depth = record.n_effective_units * args.empirical_depth_fold
        candidates = []
        for level, key in empirical_null_keys(
                record.library, record.group, record.arm):
            external = [
                candidate for candidate in pools.get((level, key), [])
                if candidate.donor_pair != record.donor_pair and
                not (record.uid and candidate.uid and
                     candidate.uid == record.uid) and
                crossfit_fold(candidate.manifest,
                              args.calibration_crossfit_folds) == query_fold and
                abs(finite_float(candidate.ase.get("ambient_c"), 0.0) -
                    query_ambient) <= args.empirical_ambient_window and
                lower_depth <= candidate.n_effective_units <= upper_depth
            ]
            if len(external) >= args.min_empirical_null_cells:
                candidates = [(level, external)]
                break
            if len(external) >= args.absolute_min_null_cells:
                candidates.append((level, external))
        chosen: Optional[Tuple[str, List[CallRecord]]] = None
        if candidates:
            chosen = max(candidates, key=lambda item: len(item[1]))
        if chosen is None:
            record.call_status = "NO_EXTERNAL_PAIR_EMPIRICAL_NULL"
            record.row_flags.append("NO_EXTERNAL_PAIR_EMPIRICAL_NULL")
            continue
        level, null_records = chosen
        null_values = [candidate.ase_test_stat for candidate in null_records]
        record.empirical_basis = (
            "EXTERNAL_PAIR_SAME_HELDOUT_FOLD_DEPTH_AND_AMBIENT_MATCHED_"
            "STATE_SPECIFIC_ASE_LOGBF")
        record.empirical_level = level
        record.empirical_null_cells = len(null_values)
        record.empirical_p_floor = 1.0 / (len(null_values) + 1.0)
        record.empirical_state_p = {
            state: empirical_p_value(
                record.ase_state_test_stats[state],
                [candidate.ase_state_test_stats[state]
                 for candidate in null_records])
            for state in STATES[1:]
        }
        record.empirical_p = min(
            1.0, (len(STATES) - 1) * min(record.empirical_state_p.values()))
        if len(null_values) < args.min_empirical_null_cells:
            record.row_flags.append("EMPIRICAL_NULL_BELOW_PREFERRED_SIZE")

    state_hypotheses = [
        (record, state, record.empirical_state_p[state])
        for record in records for state in STATES[1:]
        if state in record.empirical_state_p
    ]
    state_q_values = benjamini_hochberg(
        [p_value for _record, _state, p_value in state_hypotheses])
    for (record, state, _p_value), q_value in zip(
            state_hypotheses, state_q_values):
        record.empirical_state_q[state] = q_value

    q_values = benjamini_hochberg([record.empirical_p for record in records])
    all_state_p_values = [p_value for _record, _state, p_value in state_hypotheses]
    for record, q_value in zip(records, q_values):
        record.empirical_q = q_value
        if math.isfinite(record.empirical_p_floor):
            record.empirical_q_resolution_floor = virtual_bh_q_bound(
                all_state_p_values, record.empirical_p_floor)
        if record.call_status != "SCORED_PENDING_EMPIRICAL_NULL":
            continue
        review_reasons = []
        context_reasons = manifest_event_review_reasons(record)
        if "HIGH_QNAME_FALLBACK_FRACTION" in record.row_flags:
            review_reasons.append("HIGH_QNAME_FALLBACK")
        if "LOW_AMBIENT_GENOTYPED_MASS" in record.row_flags:
            review_reasons.append("LOW_AMBIENT_GENOTYPED_MASS")
        if "SITE_LEVEL_FALLBACK_DOWNWEIGHTED" in record.row_flags:
            review_reasons.append("SITE_LEVEL_FALLBACK")
        if "NON_AUTOSOMAL_EXPLORATORY" in record.row_flags:
            review_reasons.append("NON_AUTOSOMAL_EXPLORATORY")
        if "_SENSITIVE_" in record.ambient_sensitivity_status:
            review_reasons.append("AMBIENT_C_SENSITIVE")
        if record.ambient_sensitivity_status.startswith("POINT_ESTIMATE_ONLY_"):
            review_reasons.append("AMBIENT_C_UNCERTAINTY_UNAVAILABLE")
        low_coverage = ("LOW_SITE_COUNT" in record.row_flags or
                        "LOW_EFFECTIVE_DIRECTIONAL_EXPOSURE" in record.row_flags)
        if record.best_state == "BALANCED":
            if record.empirical_null_cells < args.min_empirical_null_cells:
                record.call_state = "NO_CALL"
                record.call_status = "INSUFFICIENT_EMPIRICAL_NULL_RESOLUTION"
            elif context_reasons or review_reasons:
                record.call_state = "NO_CALL"
                record.call_status = "REVIEW_BALANCED_" + "_AND_".join(
                    context_reasons + review_reasons)
            elif low_coverage:
                record.call_state = "NO_CALL"
                record.call_status = "AGGREGATE_ONLY_LOW_COVERAGE"
            elif record.best_posterior >= args.min_balanced_posterior:
                record.call_state = "BALANCED"
                record.call_status = "PASS_BALANCED"
            else:
                record.call_status = "LOW_BEST_STATE_POSTERIOR"
        elif (math.isfinite(record.empirical_q_resolution_floor) and
              record.empirical_q_resolution_floor > args.max_event_q):
            record.call_state = "NO_CALL"
            record.call_status = "INSUFFICIENT_EMPIRICAL_Q_RESOLUTION"
        elif record.empirical_null_cells < args.min_empirical_null_cells:
            record.call_state = "NO_CALL"
            record.call_status = "INSUFFICIENT_EMPIRICAL_NULL_RESOLUTION"
        elif context_reasons or review_reasons:
            record.call_state = "NO_CALL"
            record.call_status = "REVIEW_EVENT_" + "_AND_".join(
                context_reasons + review_reasons)
        elif low_coverage:
            record.call_state = "NO_CALL"
            record.call_status = "AGGREGATE_ONLY_LOW_COVERAGE"
        elif (record.best_posterior >= args.min_event_posterior and
              math.isfinite(record.empirical_state_q.get(record.best_state, math.nan)) and
              record.empirical_state_q[record.best_state] <= args.max_event_q):
            record.call_state = record.best_state
            record.call_status = "PASS_EVENT"
        elif record.best_posterior < args.min_event_posterior:
            record.call_status = "LOW_BEST_STATE_POSTERIOR"
        else:
            record.call_status = "EMPIRICAL_Q_NOT_SIGNIFICANT"


def load_manifests(paths: Sequence[str]) -> Tuple[Dict[Tuple[str, str], Dict[str, str]], List[str]]:
    manifests: Dict[Tuple[str, str], Dict[str, str]] = {}
    uid_pairs: Dict[str, str] = {}
    used_paths = []
    required = {
        "library", "barcode", "donor_a", "donor_b", "donor_pair",
        "calibration_group", "uid", "model_eligible", "calibration_eligible",
    }
    for raw_path in paths:
        path = require_file(raw_path, "cell manifest")
        header, rows = read_rows(path, "cell manifest")
        require_columns(header, required, path)
        for row in rows:
            if clean(row.get("schema_version")) not in {"", CELL_MANIFEST_SCHEMA}:
                raise ValueError(f"{path}: incompatible cell-manifest schema")
            library = canonical_library(row.get("library"))
            barcode = canonical_barcode(row.get("barcode", ""))
            if not barcode:
                raise ValueError(f"{path}: empty barcode")
            row["library"] = library
            row["barcode"] = barcode
            key = library, barcode
            if key in manifests:
                raise ValueError(f"duplicate manifest cell lib{library}/{barcode}")
            uid = clean(row.get("uid"))
            donor_pair = pair_from_row(row)
            if uid:
                previous_pair = uid_pairs.get(uid)
                if previous_pair is not None and previous_pair != donor_pair:
                    raise ValueError(
                        f"UID {uid!r} maps to multiple donor pairs: "
                        f"{previous_pair!r} and {donor_pair!r}")
                uid_pairs[uid] = donor_pair
            manifests[key] = row
        used_paths.append(path)
    if not manifests:
        raise ValueError("cell manifests contain no rows")
    return manifests, used_paths


def load_expression(
        paths: Sequence[str],
        manifests: Mapping[Tuple[str, str], Mapping[str, str]],
        target_cells: set,
        store: ExpressionStore,
        folds: int,
        ) -> Tuple[Dict[Tuple[str, str], str], List[str], int]:
    """Stream target-cell expression rows into a disk-backed lookup index."""
    chromosomes: Dict[Tuple[str, str], str] = {}
    used_paths = []
    skipped_nonheterotypic = 0
    required = {
        "library", "barcode", "arm", "chromosome", "log2_arm_to_other",
        "total_autosomal_counts", "matrix_value_type", "expression_input_state",
    }
    for raw_path in paths:
        path = require_file(raw_path, "expression evidence")
        header, rows = read_rows(path, "expression evidence")
        require_columns(header, required, path)
        for row in rows:
            if clean(row.get("schema_version")) not in {"", EXPRESSION_SCHEMA}:
                raise ValueError(f"{path}: incompatible expression schema")
            input_state = clean(row.get("expression_input_state")).upper()
            if input_state not in {
                    "OBSERVED_FILTERED_COUNTS", "UPSTREAM_AMBIENT_CORRECTED"}:
                raise ValueError(
                    f"{path}: invalid expression_input_state {input_state!r}")
            matrix_type = clean(row.get("matrix_value_type")).upper()
            if not expression_matrix_supported(matrix_type, input_state):
                raise ValueError(
                    f"{path}: matrix_value_type {matrix_type!r} is incompatible "
                    f"with expression_input_state {input_state!r}")
            library = canonical_library(row.get("library"))
            barcode = canonical_barcode(row.get("barcode", ""))
            arm = clean(row.get("arm"))
            if not barcode or not arm:
                raise ValueError(f"{path}: empty expression key")
            cell_key = library, barcode
            if cell_key not in manifests:
                raise ValueError(
                    f"expression cell absent from manifest: lib{library}/{barcode}")
            if cell_key not in target_cells:
                skipped_nonheterotypic += 1
                continue
            row["library"], row["barcode"] = library, barcode
            arm_key = library, arm
            chromosome = clean(row.get("chromosome"))
            previous = chromosomes.get(arm_key)
            if previous and chromosome and previous != chromosome:
                raise ValueError(f"conflicting chromosome for {arm_key}")
            if chromosome:
                chromosomes[arm_key] = chromosome
            manifest = manifests[cell_key]
            store.add(
                library, barcode, arm, chromosome, row, manifest,
                crossfit_fold(manifest, folds))
        used_paths.append(path)
    store.finish()
    return chromosomes, used_paths, skipped_nonheterotypic


def load_ase(paths: Sequence[str], manifests: Mapping[Tuple[str, str], Mapping[str, str]],
             target_cells: set, expression: ExpressionStore, args
             ) -> Tuple[List[CallRecord], Dict[Tuple[str, str], Dict[str, str]],
                        List[str], int]:
    records: List[CallRecord] = []
    seen = set()
    arm_metadata: Dict[Tuple[str, str], Dict[str, str]] = {}
    used_paths = []
    skipped_nonheterotypic = 0
    required = set(ASE_RAW_FIELDS)
    for raw_path in paths:
        path = require_file(raw_path, "ASE evidence")
        header, rows = read_rows(path, "ASE evidence")
        require_columns(header, required, path)
        for row in rows:
            if clean(row.get("schema_version")) != ASE_INPUT_SCHEMA:
                raise ValueError(
                    f"{path}: soft-primary caller requires {ASE_INPUT_SCHEMA}")
            library = canonical_library(row.get("library"))
            barcode = canonical_barcode(row.get("barcode", ""))
            arm = clean(row.get("arm"))
            key = library, barcode, arm
            if not barcode or not arm or key in seen:
                raise ValueError(f"empty/duplicate ASE key: {key}")
            seen.add(key)
            manifest = manifests.get((library, barcode))
            if manifest is None:
                raise ValueError(f"ASE cell absent from manifest: lib{library}/{barcode}")
            if (clean(row.get("donor_a")) != clean(manifest.get("donor_a")) or
                    clean(row.get("donor_b")) != clean(manifest.get("donor_b"))):
                raise ValueError(
                    f"ordered donor_a/donor_b mismatch for lib{library}/{barcode}")
            if pair_from_row(row) != pair_from_row(manifest):
                raise ValueError(f"donor-pair label mismatch for lib{library}/{barcode}")
            ambient_c_value = finite_float(row.get("ambient_c"))
            if (not math.isfinite(ambient_c_value) or
                    not 0.0 <= ambient_c_value < 1.0):
                raise ValueError(
                    f"invalid ambient_c for lib{library}/{barcode}: "
                    f"{row.get('ambient_c')!r}")
            ambient_se_value = finite_float(row.get("ambient_c_se"))
            if math.isfinite(ambient_se_value) and ambient_se_value < 0.0:
                raise ValueError(
                    f"negative ambient_c_se for lib{library}/{barcode}: "
                    f"{row.get('ambient_c_se')!r}")
            for ambient_field in ("ambient_c", "ambient_c_se"):
                ase_value = finite_float(row.get(ambient_field))
                manifest_value = finite_float(manifest.get(ambient_field))
                if (math.isfinite(ase_value) != math.isfinite(manifest_value) or
                        (math.isfinite(ase_value) and not math.isclose(
                            ase_value, manifest_value, rel_tol=1e-8, abs_tol=1e-10))):
                    raise ValueError(
                        f"{ambient_field} mismatch for lib{library}/{barcode}")
            if (library, barcode) not in target_cells:
                skipped_nonheterotypic += 1
                continue
            row["library"], row["barcode"] = library, barcode
            counts = {}
            soft_evidence: Dict[str, Tuple[float, float, int]] = {}
            n_soft_units = 0
            n_effective_units = 0.0
            context = f"lib{library}/{barcode}/{arm}"
            for orientation, (a_field, b_field, _) in ORIENTATION_FIELDS.items():
                counts[orientation] = (
                    parse_count(row.get(a_field), a_field, context),
                    parse_count(row.get(b_field), b_field, context),
                )
                unit_field = f"n_molecules_{orientation}"
                soft_a_field = f"soft_a_{orientation}"
                soft_b_field = f"soft_b_{orientation}"
                sumsq_field = f"soft_a_sumsq_{orientation}"
                effective_a_field = f"effective_a_{orientation}"
                effective_b_field = f"effective_b_{orientation}"
                effective_weight_field = f"effective_weight_{orientation}"
                units = parse_count(row.get(unit_field), unit_field, context)
                soft_a = parse_nonnegative_float(
                    row.get(soft_a_field), soft_a_field, context)
                soft_b = parse_nonnegative_float(
                    row.get(soft_b_field), soft_b_field, context)
                sumsq = parse_nonnegative_float(
                    row.get(sumsq_field), sumsq_field, context)
                tolerance = max(1e-7, 1e-6 * max(units, 1))
                if not math.isclose(
                        soft_a + soft_b, units, rel_tol=1e-6, abs_tol=tolerance):
                    raise ValueError(
                        f"{context}: {soft_a_field}+{soft_b_field} != {unit_field}")
                lower_sumsq = soft_a * soft_a / units if units else 0.0
                if sumsq < lower_sumsq - tolerance or sumsq > soft_a + tolerance:
                    raise ValueError(
                        f"{context}: {sumsq_field} violates fractional-support bounds")
                effective_a = parse_nonnegative_float(
                    row.get(effective_a_field), effective_a_field, context)
                effective_b = parse_nonnegative_float(
                    row.get(effective_b_field), effective_b_field, context)
                effective_units = parse_nonnegative_float(
                    row.get(effective_weight_field), effective_weight_field, context)
                if (not math.isclose(effective_a + effective_b, effective_units,
                                     rel_tol=1e-6, abs_tol=tolerance) or
                        effective_units > units + tolerance):
                    raise ValueError(
                        f"{context}: inconsistent effective-weight universe for {orientation}")
                identity_weight = 4.0 * sumsq - 4.0 * soft_a + units
                if not math.isclose(
                        effective_units, identity_weight,
                        rel_tol=1e-6, abs_tol=tolerance):
                    raise ValueError(
                        f"{context}: {effective_weight_field} violates the soft-sumsq identity")
                ambient_value = finite_float(
                    row.get(f"ambient_a_{orientation}"))
                if effective_units > 0.0 and (
                        not math.isfinite(ambient_value) or
                        not 0.0 <= ambient_value <= 1.0):
                    raise ValueError(
                        f"{context}: invalid effective-weighted ambient_a_{orientation}")
                soft_evidence[orientation] = (
                    effective_a, effective_units, units)
                n_soft_units += units
                n_effective_units += effective_units
            global_soft_a = parse_nonnegative_float(
                row.get("soft_a"), "soft_a", context)
            global_soft_b = parse_nonnegative_float(
                row.get("soft_b"), "soft_b", context)
            orientation_soft_a = sum(parse_nonnegative_float(
                row.get(f"soft_a_{orientation}"), f"soft_a_{orientation}", context)
                for orientation in ORIENTATIONS)
            orientation_soft_b = sum(parse_nonnegative_float(
                row.get(f"soft_b_{orientation}"), f"soft_b_{orientation}", context)
                for orientation in ORIENTATIONS)
            if (not math.isclose(global_soft_a, orientation_soft_a,
                                 rel_tol=1e-6, abs_tol=1e-6) or
                    not math.isclose(global_soft_b, orientation_soft_b,
                                     rel_tol=1e-6, abs_tol=1e-6)):
                raise ValueError(
                    f"{context}: global soft totals do not equal orientation totals")
            n_sites = parse_count(row.get("n_sites"), "n_sites", context)
            n_informative = parse_count(
                row.get("n_informative_molecules"), "n_informative_molecules", context)
            flags = []
            if sum(a + b for a, b in counts.values()) != n_informative:
                flags.append("INFORMATIVE_COUNT_MISMATCH")
            group = clean(manifest.get("calibration_group")) or f"lib{library}"
            chromosome = clean(row.get("chromosome"))
            expression_row = expression.get(key)
            if (expression_row is not None and chromosome_key(
                    expression_row.get("chromosome")) != chromosome_key(chromosome)):
                raise ValueError(
                    f"ASE/expression chromosome mismatch for {key}")
            record = CallRecord(
                library, barcode, arm, chromosome, group, pair_from_row(manifest),
                manifest, row, expression_row, counts, soft_evidence, n_sites,
                n_informative, n_soft_units, n_effective_units, flags)
            records.append(record)
            metadata_key = library, arm
            metadata = {
                "chromosome": chromosome,
                "arm_start": clean(row.get("arm_start")) or "NA",
                "arm_end": clean(row.get("arm_end")) or "NA",
            }
            previous = arm_metadata.get(metadata_key)
            if previous is not None and previous != metadata:
                raise ValueError(f"conflicting ASE arm metadata for {metadata_key}")
            arm_metadata[metadata_key] = metadata
        used_paths.append(path)
    return records, arm_metadata, used_paths, skipped_nonheterotypic


def calibration_rows(resolver: CalibrationResolver) -> Iterable[Dict[str, object]]:
    """Stream one query-level calibration row per resolved call model."""
    for record_key in sorted(resolver.resolved_cache, key=lambda value: tuple(
            natural_key(part) for part in value)):
        library, barcode, arm = record_key
        record = resolver.record_by_key[record_key]
        group, donor_pair = record.group, record.donor_pair
        model = resolver.resolved_cache[record_key]
        row: Dict[str, object] = {
            "library": library,
            "barcode": barcode,
            "uid": record.uid or "NA",
            "calibration_group": group,
            "donor_pair": donor_pair,
            "chromosome": record.chromosome,
            "arm": arm,
            "crossfit_fold": model.crossfit_fold,
            "excluded_chromosome": model.excluded_chromosome,
            "target_effective_ase_weight": format_number(
                record.n_effective_units),
            "group_genomewide_baseline_logit": format_number(
                model.group_baseline_logit),
            "cell_loo_baseline_logit": format_number(model.cell_baseline_logit),
            "cell_loo_baseline_donor_a_fraction": format_number(
                logistic(model.cell_baseline_logit)),
            "cell_loo_baseline_arms": model.cell_baseline_arms,
            "cell_loo_retained_observations": (
                model.cell_baseline_retained_observations),
            "cell_loo_baseline_status": model.cell_baseline_status,
            "exact_calibration_cells": model.exact_cells,
            "expression_source_level": model.expression.source_level,
            "expression_source_key": key_text(model.expression.source_key) or "NA",
            "expression_calibration_cells": model.expression.fit.n_cells,
            "expression_retained_calibration_cells": (
                model.expression.fit.retained_cells),
            "expression_input_state": clean(
                record.expression.get("expression_input_state")).upper()
            if record.expression is not None else "NA",
            "expression_metric": model.expression.metric_field,
            "ambient_handling": (
                "PRODUCER_FULL_SIMPLEX_AMBIENT_A;GENOTYPED_MASS_QC_ONLY"),
            "expression_baseline_log2_arm_metric": format_number(
                model.expression.fit.center) if model.expression.data_calibrated else "NA",
            "expression_sigma": format_number(
                model.expression.fit.sigma) if model.expression.data_calibrated else "NA",
            "calibration_status": model.status,
            "schema_version": CALIBRATION_OUTPUT_SCHEMA,
        }
        for orientation in ORIENTATIONS:
            resolved = model.biases[orientation]
            prefix = orientation
            row[f"{prefix}_source_level"] = resolved.source_level
            row[f"{prefix}_source_key"] = key_text(resolved.source_key) or "NA"
            row[f"{prefix}_calibration_cells"] = resolved.fit.n_cells
            row[f"{prefix}_calibration_effective_weight"] = format_number(
                resolved.fit.total_count)
            row[f"{prefix}_retained_calibration_cells"] = (
                resolved.fit.retained_cells)
            row[f"{prefix}_orientation_mapping_logit_offset"] = (
                format_number(resolved.fit.delta) if resolved.data_calibrated else "NA")
            row[f"{prefix}_observed_donor_a_fraction"] = format_number(
                resolved.fit.observed_fraction)
            row[f"{prefix}_quasi_overdispersion_rho"] = (
                format_number(resolved.fit.rho) if resolved.data_calibrated else "NA")
        yield row


def base_no_data_ase(library: str, barcode: str, arm: str,
                     manifest: Mapping[str, str],
                     metadata: Mapping[str, str]) -> Dict[str, str]:
    row = {field: "NA" for field in ASE_RAW_FIELDS}
    row.update({
        "library": library,
        "barcode": barcode,
        "donor_a": clean(manifest.get("donor_a")) or "NA",
        "donor_b": clean(manifest.get("donor_b")) or "NA",
        "donor_pair": pair_from_row(manifest),
        "arm": arm,
        "chromosome": clean(metadata.get("chromosome")) or "NA",
        "arm_start": clean(metadata.get("arm_start")) or "NA",
        "arm_end": clean(metadata.get("arm_end")) or "NA",
        "ambient_c": clean(manifest.get("ambient_c")) or "NA",
        "ambient_c_se": clean(manifest.get("ambient_c_se")) or "NA",
        "n_sites": "0", "n_molecules": "0", "n_informative_molecules": "0",
        "a_ref": "0", "b_ref": "0", "a_alt": "0", "b_alt": "0",
        "a_mixed": "0", "b_mixed": "0", "n_ambiguous": "0",
        "soft_a": "0", "soft_b": "0", "n_molecules_ref": "0",
        "n_molecules_alt": "0", "n_molecules_mixed": "0",
        "soft_a_ref": "0", "soft_b_ref": "0", "soft_a_alt": "0",
        "soft_b_alt": "0", "soft_a_mixed": "0", "soft_b_mixed": "0",
        "soft_a_sumsq_ref": "0", "soft_a_sumsq_alt": "0",
        "soft_a_sumsq_mixed": "0", "effective_a_ref": "0",
        "effective_b_ref": "0", "effective_a_alt": "0",
        "effective_b_alt": "0", "effective_a_mixed": "0",
        "effective_b_mixed": "0", "effective_weight_ref": "0",
        "effective_weight_alt": "0", "effective_weight_mixed": "0",
        "ambient_genotyped_mass": "NA",
        "qname_fallback_fraction": "NA", "mean_sites_per_molecule": "NA",
        "evidence_basis": "NO_DATA",
        "model_eligible": "1" if truthy(manifest.get("model_eligible")) else "0",
        "evidence_status": "NO_DATA", "schema_version": "NA",
    })
    return row


def output_row(record: Optional[CallRecord], library: str, barcode: str, arm: str,
               manifest: Mapping[str, str], expression: Optional[Mapping[str, str]],
               metadata: Mapping[str, str]) -> Dict[str, object]:
    row: Dict[str, object] = dict(record.ase) if record is not None else base_no_data_ase(
        library, barcode, arm, manifest, metadata)
    for field_name in ASE_RAW_FIELDS:
        row.setdefault(field_name, "NA")
    for field_name in MANIFEST_AUDIT_FIELDS:
        row[f"manifest_{field_name}"] = clean(manifest.get(field_name)) or "NA"
    for field_name in EXPRESSION_VALUE_FIELDS:
        row[f"expression_{field_name}"] = (
            clean(expression.get(field_name)) if expression is not None else "NA") or "NA"

    row["calibration_group"] = clean(manifest.get("calibration_group")) or f"lib{library}"
    row["uid"] = clean(manifest.get("uid")) or "NA"
    row["expression_role"] = "SUPPORTIVE_TOTAL_COPY_INPUT_STATE_AWARE"
    row["bh_scope"] = "INVOCATION"
    row["call_schema_version"] = CALL_OUTPUT_SCHEMA

    if record is None:
        for field_name in MODEL_FIELDS:
            row.setdefault(field_name, "NA")
        row.update({
            "calibration_group": clean(manifest.get("calibration_group")) or f"lib{library}",
            "uid": clean(manifest.get("uid")) or "NA",
            "expression_used": "0",
            "expression_role": "SUPPORTIVE_TOTAL_COPY_INPUT_STATE_AWARE",
            "best_state": "NO_CALL", "call_state": "NO_CALL",
            "call_status": "NO_DATA", "qc_flags": "NO_ASE_EVIDENCE",
            "bh_scope": "INVOCATION", "call_schema_version": CALL_OUTPUT_SCHEMA,
        })
        return row

    model = record.model
    row.update({
        "calibration_level": model.level if model else "NA",
        "calibration_status": model.status if model else "NA",
        "calibration_cells_exact": model.exact_cells if model else 0,
        "calibration_cells_effective": max(
            (bias.fit.n_cells for bias in model.biases.values()), default=0)
            if model else 0,
        "total_effective_ase_weight": format_number(record.n_effective_units),
        "calibration_crossfit_fold": model.crossfit_fold if model else "NA",
        "calibration_excluded_chromosome": (
            model.excluded_chromosome if model else "NA"),
        "group_genomewide_baseline_logit": format_number(
            model.group_baseline_logit) if model else "NA",
        "cell_loo_baseline_logit": format_number(
            model.cell_baseline_logit) if model else "NA",
        "cell_loo_baseline_donor_a_fraction": format_number(
            logistic(model.cell_baseline_logit)) if model else "NA",
        "cell_loo_baseline_arms": model.cell_baseline_arms if model else 0,
        "cell_loo_retained_observations": (
            model.cell_baseline_retained_observations if model else 0),
        "cell_loo_baseline_status": model.cell_baseline_status if model else "NA",
        "expression_used": int(record.expression_used),
        "expression_metric": model.expression.metric_field if model else "NA",
        "expression_calibration_level": (
            model.expression.source_level if model else "NA"),
        "expression_calibration_cells": (
            model.expression.fit.n_cells if model else 0),
        "expression_baseline_log2_arm_metric": (
            format_number(model.expression.fit.center)
            if model and model.expression.data_calibrated else "NA"),
        "expression_sigma": format_number(record.expression_sigma_effective),
        "best_state": record.best_state,
        "posterior_interpretation": "APPROXIMATE_WORKING_PSEUDO_POSTERIOR",
        "best_state_posterior": format_number(record.best_posterior),
        "event_posterior": format_number(record.event_posterior),
        "ase_log_bf_best_vs_balanced": format_number(
            record.ase_relative.get(record.best_state, math.nan)),
        "expression_log_bf_best_vs_balanced": format_number(
            record.expression_relative.get(record.best_state, math.nan)),
        "joint_log_bf_best_vs_balanced": format_number(
            record.joint_relative.get(record.best_state, math.nan)),
        "ambient_sensitivity_status": record.ambient_sensitivity_status,
        "empirical_test_basis": record.empirical_basis,
        "empirical_null_level": record.empirical_level,
        "empirical_null_cells": record.empirical_null_cells,
        "empirical_p_floor": format_number(record.empirical_p_floor),
        "empirical_q_resolution_floor": format_number(
            record.empirical_q_resolution_floor),
        "empirical_p_value": format_number(record.empirical_p),
        "empirical_q_value": format_number(record.empirical_q),
        "call_state": record.call_state,
        "call_status": record.call_status,
        "qc_flags": ";".join(sorted(set(record.row_flags), key=natural_key)) or "PASS",
    })
    for orientation in ORIENTATIONS:
        bias = model.biases[orientation] if model else None
        row[f"orientation_mapping_logit_offset_{orientation}"] = (
            format_number(bias.fit.delta) if bias and bias.data_calibrated else "NA")
        row[f"quasi_overdispersion_rho_{orientation}"] = (
            format_number(bias.fit.rho)
            if bias and bias.data_calibrated else "NA")
        row[f"expected_a_balanced_{orientation}"] = format_number(
            record.expected_balanced.get(orientation, math.nan))
        row[f"expected_a_best_{orientation}"] = format_number(
            record.expected_best.get(orientation, math.nan))
    for state in STATES:
        row[f"posterior_{state}"] = format_number(
            record.posteriors.get(state, math.nan))
        if state != "BALANCED":
            row[f"empirical_p_{state}"] = format_number(
                record.empirical_state_p.get(state, math.nan))
            row[f"empirical_q_{state}"] = format_number(
                record.empirical_state_q.get(state, math.nan))
    for field_name in MODEL_FIELDS:
        row.setdefault(field_name, "NA")
    return row


@dataclass(slots=True)
class UIDArmSummary:
    uid: str
    donor_pair: str
    chromosome: str
    side: str
    arm: str
    records: List[CallRecord]
    posterior: Dict[str, float]
    state: str
    best_posterior: float
    event_posterior: float
    p_value: float
    q_value: float = math.nan
    total_cells: int = 0
    evaluable_cells: int = 0
    supporting_cells: int = 0
    concordant_cells: int = 0
    concordance: float = 0.0
    qc_blocked_cells: int = 0
    state_p: Dict[str, float] = field(default_factory=dict)
    state_q: Dict[str, float] = field(default_factory=dict)
    p_floor: float = math.nan
    q_resolution_floor: float = math.nan


def aggregate_state_posterior(
        records: Sequence[CallRecord], args
        ) -> Tuple[Dict[str, float], str, float, float, float]:
    log_bayes_factors = {
        state: sum(bounded(record.joint_relative.get(state, 0.0),
                           -args.max_aggregate_cell_log_bf,
                           args.max_aggregate_cell_log_bf)
                   for record in records)
        for state in STATES
    }
    priors = {"BALANCED": 1.0 - args.event_prior}
    priors.update({state: args.event_prior / (len(STATES) - 1)
                   for state in STATES[1:]})
    log_scores = {
        state: math.log(priors[state]) + log_bayes_factors[state]
        for state in STATES
    }
    normalizer = logsumexp(list(log_scores.values()))
    posterior = {
        state: math.exp(log_scores[state] - normalizer) for state in STATES
    }
    state = max(STATES, key=lambda candidate: posterior[candidate])
    return (posterior, state, posterior[state], 1.0 - posterior["BALANCED"],
            log_bayes_factors[state])


def build_uid_rows(records: Sequence[CallRecord], args) -> List[Dict[str, object]]:
    grouped: MutableMapping[
        Tuple[str, str, str, str], List[CallRecord]
    ] = defaultdict(list)
    for record in records:
        side = arm_side(record.arm)
        if not record.uid or record.uid.upper() == "NA" or not side:
            continue
        grouped[(record.uid, record.donor_pair, record.chromosome, side)].append(record)

    summaries: List[UIDArmSummary] = []
    for (uid, donor_pair, chromosome, side), values in sorted(
            grouped.items(), key=lambda item: tuple(natural_key(part) for part in item[0])):
        evaluable = [
            record for record in values if aggregate_qc_eligible(record, args)]
        if evaluable:
            posterior, state, best_posterior, event_posterior, _aggregate_lbf = (
                aggregate_state_posterior(evaluable, args))
        else:
            posterior = {state_name: math.nan for state_name in STATES}
            state, best_posterior, event_posterior = "NO_CALL", math.nan, math.nan
        total_cells = len({(record.library, record.barcode) for record in values})
        evaluable_cells = len({(record.library, record.barcode)
                               for record in evaluable})
        supporting = [record for record in evaluable
                      if record.call_status == "PASS_EVENT"]
        concordant = [record for record in supporting
                      if record.call_state == state]
        supporting_cells = len({(record.library, record.barcode)
                                for record in supporting})
        concordant_cells = len({(record.library, record.barcode)
                                for record in concordant})
        concordance = concordant_cells / evaluable_cells if evaluable_cells else 0.0
        blocked_cells = total_cells - evaluable_cells
        summary = UIDArmSummary(
            uid, donor_pair, chromosome, side,
            ";".join(sorted({record.arm for record in values}, key=natural_key)),
            values, posterior, state, best_posterior, event_posterior,
            math.nan, math.nan, total_cells,
            evaluable_cells, supporting_cells, concordant_cells, concordance,
            blocked_cells)
        summary.state_p = {
            state_name: partial_conjunction_p(
                [record.empirical_state_p.get(state_name, math.nan)
                 for record in evaluable], args.min_uid_cells_per_arm)
            for state_name in STATES[1:]
        }
        summary.p_floor = partial_conjunction_p(
            [record.empirical_p_floor for record in evaluable],
            args.min_uid_cells_per_arm)
        summaries.append(summary)

    hypotheses = [
        (summary, state_name, summary.state_p[state_name])
        for summary in summaries for state_name in STATES[1:]
        if math.isfinite(summary.state_p[state_name])
    ]
    q_values = benjamini_hochberg(
        [p_value for _summary, _state, p_value in hypotheses])
    for (summary, state_name, _p_value), q_value in zip(hypotheses, q_values):
        summary.state_q[state_name] = q_value
    all_uid_state_p = [p_value for _summary, _state, p_value in hypotheses]
    for summary in summaries:
        if summary.state in STATES[1:]:
            summary.p_value = summary.state_p.get(summary.state, math.nan)
            summary.q_value = summary.state_q.get(summary.state, math.nan)
            summary.q_resolution_floor = virtual_bh_q_bound(
                all_uid_state_p, summary.p_floor)

    by_chromosome: MutableMapping[
        Tuple[str, str, str], Dict[str, UIDArmSummary]
    ] = defaultdict(dict)
    for summary in summaries:
        by_chromosome[(summary.uid, summary.donor_pair, summary.chromosome)][
            summary.side] = summary

    output = []
    for (uid, donor_pair, chromosome), sides in sorted(
            by_chromosome.items(),
            key=lambda item: tuple(natural_key(part) for part in item[0])):
        p_arm, q_arm = sides.get("p"), sides.get("q")
        all_records = ((p_arm.records if p_arm else []) +
                       (q_arm.records if q_arm else []))
        libraries = sorted({record.library for record in all_records}, key=natural_key)
        groups = sorted({record.group for record in all_records}, key=natural_key)
        cells = {(record.library, record.barcode) for record in all_records}
        evaluable = (
            p_arm is not None and q_arm is not None and
            p_arm.evaluable_cells >= args.min_uid_cells_per_arm and
            q_arm.evaluable_cells >= args.min_uid_cells_per_arm)
        flag = "NOT_EVALUABLE"
        state = "NA"
        status = "MISSING_P_OR_Q_ARM"
        paired_concordant = 0
        if evaluable:
            assert p_arm is not None and q_arm is not None
            same_event = p_arm.state == q_arm.state and p_arm.state != "BALANCED"
            p_concordant = {
                (record.library, record.barcode) for record in p_arm.records
                if record.call_status == "PASS_EVENT" and
                record.call_state == p_arm.state
            }
            q_concordant = {
                (record.library, record.barcode) for record in q_arm.records
                if record.call_status == "PASS_EVENT" and
                record.call_state == q_arm.state
            }
            paired_concordant = len(p_concordant & q_concordant)
            significant = (
                p_arm.best_posterior >= args.min_uid_arm_posterior and
                q_arm.best_posterior >= args.min_uid_arm_posterior and
                p_arm.q_value <= args.max_whole_chromosome_q and
                q_arm.q_value <= args.max_whole_chromosome_q and
                p_arm.q_resolution_floor <= args.max_whole_chromosome_q and
                q_arm.q_resolution_floor <= args.max_whole_chromosome_q and
                p_arm.concordant_cells >= args.min_uid_cells_per_arm and
                q_arm.concordant_cells >= args.min_uid_cells_per_arm and
                p_arm.concordance >= args.min_uid_state_concordance and
                q_arm.concordance >= args.min_uid_state_concordance)
            if (same_event and significant and
                    paired_concordant >= args.min_uid_cells_per_arm):
                flag, state, status = "TRUE", p_arm.state, "PASS_WHOLE_CHROMOSOME"
            elif same_event and significant:
                flag, state, status = (
                    "FALSE", p_arm.state,
                    "GROUP_LEVEL_PQ_SUPPORT_WITHOUT_PAIRED_CELL_OVERLAP")
            elif p_arm.state != q_arm.state:
                flag, status = "FALSE", "DISCORDANT_P_Q_STATES"
            else:
                flag, status = "FALSE", "P_Q_NOT_JOINTLY_SIGNIFICANT"

        def side_values(summary: Optional[UIDArmSummary], prefix: str) -> Dict[str, object]:
            if summary is None:
                return {
                    f"{prefix}_arm": "NA", f"{prefix}_rows": 0,
                    f"{prefix}_total_cells": 0,
                    f"{prefix}_evaluable_cells": 0,
                    f"{prefix}_qc_blocked_cells": 0,
                    f"{prefix}_supporting_cells": 0,
                    f"{prefix}_concordant_cells": 0,
                    f"{prefix}_concordance": "NA",
                    f"{prefix}_state": "NO_CALL",
                    f"{prefix}_best_state_posterior": "NA",
                    f"{prefix}_event_posterior": "NA",
                    f"{prefix}_empirical_p_value": "NA",
                    f"{prefix}_empirical_q_value": "NA",
                    f"{prefix}_empirical_p_floor": "NA",
                    f"{prefix}_empirical_q_resolution_floor": "NA",
                    **{
                        f"{prefix}_state_{kind}_{state_name}": "NA"
                        for state_name in STATES[1:] for kind in ("p", "q")
                    },
                }
            return {
                f"{prefix}_arm": summary.arm,
                f"{prefix}_rows": len(summary.records),
                f"{prefix}_total_cells": summary.total_cells,
                f"{prefix}_evaluable_cells": summary.evaluable_cells,
                f"{prefix}_qc_blocked_cells": summary.qc_blocked_cells,
                f"{prefix}_supporting_cells": summary.supporting_cells,
                f"{prefix}_concordant_cells": summary.concordant_cells,
                f"{prefix}_concordance": format_number(summary.concordance),
                f"{prefix}_state": summary.state,
                f"{prefix}_best_state_posterior": format_number(summary.best_posterior),
                f"{prefix}_event_posterior": format_number(summary.event_posterior),
                f"{prefix}_empirical_p_value": format_number(summary.p_value),
                f"{prefix}_empirical_q_value": format_number(summary.q_value),
                f"{prefix}_empirical_p_floor": format_number(summary.p_floor),
                f"{prefix}_empirical_q_resolution_floor": format_number(
                    summary.q_resolution_floor),
                **{
                    f"{prefix}_state_p_{state_name}": format_number(
                        summary.state_p.get(state_name, math.nan))
                    for state_name in STATES[1:]
                },
                **{
                    f"{prefix}_state_q_{state_name}": format_number(
                        summary.state_q.get(state_name, math.nan))
                    for state_name in STATES[1:]
                },
            }

        row: Dict[str, object] = {
            "uid": uid, "donor_pair": donor_pair, "chromosome": chromosome,
            "libraries": ";".join(libraries),
            "calibration_groups": ";".join(groups), "n_cells": len(cells),
            "whole_chromosome_flag": flag, "whole_chromosome_state": state,
            "paired_pq_concordant_cells": paired_concordant,
            "summary_status": status,
            "schema_version": UID_OUTPUT_SCHEMA,
        }
        row.update(side_values(p_arm, "p"))
        row.update(side_values(q_arm, "q"))
        output.append(row)
    return output


def aggregate_qc_eligible(record: CallRecord, args) -> bool:
    blocking_flags = {
        "HIGH_QNAME_FALLBACK_FRACTION", "LOW_AMBIENT_GENOTYPED_MASS",
        "SITE_LEVEL_FALLBACK_DOWNWEIGHTED", "AMBIENT_C_SENSITIVE",
        "MANIFEST_IDENTITY_CONTEXT_REVIEW",
        "MANIFEST_TECHNICAL_CONTEXT_REVIEW",
        "TETRAPLOID_BASELINE_UNSUPPORTED",
        "NON_AUTOSOMAL_EXPLORATORY",
    }
    ambient_uncertainty_block = not ambient_uncertainty_available(record)
    return (bool(record.posteriors) and math.isfinite(record.empirical_p) and
            record.empirical_null_cells >= args.min_empirical_null_cells and
            clean(record.ase.get("evidence_status")) == "PASS" and
            not ambient_uncertainty_block and
            not blocking_flags.intersection(record.row_flags))


def build_pair_summary_rows(
        records: Sequence[CallRecord], args) -> List[Dict[str, object]]:
    grouped: MutableMapping[
        Tuple[str, str, str, str, str], List[CallRecord]
    ] = defaultdict(list)
    for record in records:
        grouped[(record.library, record.group, record.donor_pair, record.arm,
                 record.chromosome)].append(record)

    provisional = []
    for (library, group, donor_pair, arm, chromosome), values in sorted(
            grouped.items(), key=lambda item: tuple(
                natural_key(part) for part in item[0])):
        aggregate_eligible = [
            record for record in values if aggregate_qc_eligible(record, args)]
        evaluable, tested_uid_blocks = pair_test_records(
            aggregate_eligible, args)
        if evaluable:
            posterior, state, best_posterior, event_posterior, aggregate_lbf = (
                aggregate_state_posterior(evaluable, args))
        else:
            posterior = {state_name: math.nan for state_name in STATES}
            state, best_posterior, event_posterior, aggregate_lbf = (
                "NO_CALL", math.nan, math.nan, math.nan)
        candidates = [
            record for record in evaluable
            if record.best_state != "BALANCED" and
            record.best_posterior >= args.min_pair_cell_posterior
        ]
        concordant = [record for record in candidates if record.best_state == state]
        total_cells = len({(record.library, record.barcode) for record in values})
        aggregate_eligible_cells = len({
            (record.library, record.barcode) for record in aggregate_eligible
        })
        eligible_uid_blocks = len({
            record.uid or f"CELL:{record.library}:{record.barcode}"
            for record in aggregate_eligible
        })
        tested_cells = len({(record.library, record.barcode)
                            for record in evaluable})
        individual_evaluable_cells = len({
            (record.library, record.barcode) for record in values
            if record.call_status in {"PASS_EVENT", "PASS_BALANCED"}
        })
        supporting_cells = len({(record.library, record.barcode)
                                for record in candidates})
        concordant_cells = len({(record.library, record.barcode)
                                for record in concordant})
        concordance = concordant_cells / tested_cells if tested_cells else 0.0
        qc_blocked = total_cells - len({
            (record.library, record.barcode) for record in aggregate_eligible
        })
        state_p_values = {
            state_name: fisher_partial_conjunction_p(
                [record.empirical_state_p.get(state_name, math.nan)
                 for record in evaluable], args.min_pair_recurrence_cells)
            for state_name in STATES[1:]
        }
        state_p_floors = {
            state_name: fisher_partial_conjunction_p(
                [record.empirical_p_floor for record in evaluable],
                args.min_pair_recurrence_cells)
            for state_name in STATES[1:]
        }
        fisher_terms = max(
            0, len(evaluable) - args.min_pair_recurrence_cells + 1)
        provisional.append({
            "library": library,
            "calibration_group": group,
            "donor_pair": donor_pair,
            "arm": arm,
            "chromosome": chromosome,
            "libraries": ";".join(sorted(
                {record.library for record in values}, key=natural_key)),
            "total_cells": total_cells,
            "individual_evaluable_cells": individual_evaluable_cells,
            "aggregate_eligible_cells": aggregate_eligible_cells,
            "eligible_uid_blocks": eligible_uid_blocks,
            "tested_cells": tested_cells,
            "tested_uid_blocks": tested_uid_blocks,
            "qc_blocked_cells": qc_blocked,
            "supporting_cells": supporting_cells,
            "concordant_cells": concordant_cells,
            "concordance": format_number(concordance),
            "best_state": state,
            "best_state_posterior": format_number(best_posterior),
            "event_posterior": format_number(event_posterior),
            "aggregate_log_bf_best_vs_balanced": format_number(aggregate_lbf),
            "pooled_sites": sum(record.n_sites for record in evaluable),
            "pooled_soft_molecule_units": sum(
                record.n_soft_units for record in evaluable),
            "pooled_effective_units": format_number(sum(
                record.n_effective_units for record in evaluable)),
            "pooled_hard_informative_molecules": sum(
                record.n_informative for record in evaluable),
            "partial_conjunction_r": args.min_pair_recurrence_cells,
            "partial_conjunction_method": (
                "FISHER_ON_P_ORDER_R_THROUGH_P_ORDER_N"),
            "fisher_terms": fisher_terms,
            "fisher_df": 2 * fisher_terms,
            "dependence_assumption": (
                "APPROXIMATE_CONDITIONAL_INDEPENDENCE_AFTER_SCORE_BLIND_UID_CAP"),
            "fdr_interpretation": "WORKING_FDR_APPROXIMATE",
            "_state_p_values": state_p_values,
            "_state_p_floors": state_p_floors,
            "_best_posterior": best_posterior,
            "_concordance": concordance,
            "_state": state,
            "_concordant": concordant_cells,
        })

    state_hypotheses = [
        (row, state_name, row["_state_p_values"][state_name])
        for row in provisional for state_name in STATES[1:]
        if math.isfinite(row["_state_p_values"][state_name])
    ]
    q_values = benjamini_hochberg(
        [p_value for _row, _state, p_value in state_hypotheses])
    for (row, state_name, _p_value), q_value in zip(state_hypotheses, q_values):
        row.setdefault("_state_q_values", {})[state_name] = q_value
    all_pair_state_p = [p_value for _row, _state, p_value in state_hypotheses]
    output = []
    for row in provisional:
        state = str(row.pop("_state"))
        best_posterior = float(row.pop("_best_posterior"))
        concordance = float(row.pop("_concordance"))
        concordant = int(row.pop("_concordant"))
        state_p_values = row.pop("_state_p_values")
        state_p_floors = row.pop("_state_p_floors")
        state_q_values = row.pop("_state_q_values", {})
        for state_name in STATES[1:]:
            row[f"state_p_{state_name}"] = format_number(
                state_p_values.get(state_name, math.nan))
            row[f"state_q_{state_name}"] = format_number(
                state_q_values.get(state_name, math.nan))
        if state in STATES[1:]:
            p_value = state_p_values.get(state, math.nan)
            q_value = state_q_values.get(state, math.nan)
            p_floor = state_p_floors.get(state, math.nan)
        else:
            finite_state_p = [value for value in state_p_values.values()
                              if math.isfinite(value)]
            p_value = min(1.0, 4.0 * min(finite_state_p)) if finite_state_p else math.nan
            q_value, p_floor = math.nan, math.nan
        q_resolution = virtual_bh_q_bound(all_pair_state_p, p_floor)
        row["partial_conjunction_p_value"] = format_number(p_value)
        row["partial_conjunction_q_value"] = format_number(q_value)
        row["partial_conjunction_p_floor"] = format_number(p_floor)
        row["partial_conjunction_q_resolution_floor"] = format_number(q_resolution)
        recurrent = (
            state in STATES[1:] and
            best_posterior >= args.min_pair_arm_posterior and
            math.isfinite(q_value) and q_value <= args.max_pair_arm_q and
            math.isfinite(q_resolution) and q_resolution <= args.max_pair_arm_q and
            concordant >= args.min_pair_recurrence_cells and
            concordance >= args.min_pair_state_concordance)
        row["recurrence_flag"] = "TRUE" if recurrent else "FALSE"
        if not math.isfinite(q_value):
            row["summary_status"] = "NO_EXTERNAL_PAIR_EMPIRICAL_EVIDENCE"
        elif recurrent:
            row["summary_status"] = (
                "PASS_WORKING_RECURRENT_PAIR_ARM_EVENT")
        elif state == "BALANCED":
            row["summary_status"] = "AGGREGATE_BALANCED_WORKING_FDR"
        else:
            row["summary_status"] = (
                "RECURRENCE_THRESHOLDS_NOT_MET_WORKING_FDR")
        row["schema_version"] = PAIR_OUTPUT_SCHEMA
        output.append(row)
    return output


def validate_args(args) -> None:
    probability_names = (
        "event_prior", "max_event_q", "max_whole_chromosome_q",
        "min_event_posterior", "min_balanced_posterior",
        "min_uid_arm_posterior", "min_uid_state_concordance",
        "min_pair_state_concordance", "min_pair_cell_posterior",
        "min_pair_arm_posterior",
        "max_pair_arm_q", "empirical_ambient_window",
        "max_qname_fallback_fraction",
        "min_ambient_genotyped_mass", "min_rho", "max_rho", "default_rho",
    )
    for name in probability_names:
        value = getattr(args, name)
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"--{name.replace('_', '-')} must be in [0,1]")
    if not 0.0 < args.event_prior < 1.0:
        raise ValueError("--event-prior must be strictly between zero and one")
    if args.min_rho <= 0.0 or args.max_rho <= args.min_rho:
        raise ValueError("invalid quasi-likelihood rho bounds")
    for name in (
        "min_call_sites", "min_calibration_sites",
        "min_fallback_calibration_cells",
        "min_empirical_null_cells", "absolute_min_null_cells",
        "min_expression_genes", "min_uid_cells_per_arm",
        "min_pair_recurrence_cells", "calibration_crossfit_folds",
        "min_cell_baseline_arms", "max_pair_test_cells",
        "max_pair_cells_per_uid",
    ):
        if getattr(args, name) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be positive")
    if args.absolute_min_null_cells > args.min_empirical_null_cells:
        raise ValueError("absolute minimum null cells cannot exceed preferred minimum")
    if args.calibration_crossfit_folds < 2:
        raise ValueError("--calibration-crossfit-folds must be at least two")
    if args.min_pair_recurrence_cells > args.max_pair_test_cells:
        raise ValueError(
            "--min-pair-recurrence-cells cannot exceed --max-pair-test-cells")
    if args.expression_weight < 0.0 or args.max_expression_log_bf < 0.0:
        raise ValueError("expression weights must be nonnegative")
    if args.empirical_depth_fold <= 1.0:
        raise ValueError("--empirical-depth-fold must exceed one")
    if (args.min_call_effective_weight <= 0.0 or
            args.min_calibration_effective_weight <= 0.0 or
            args.min_orientation_calibration_effective_weight <= 0.0):
        raise ValueError("effective-weight thresholds must be positive")
    if not 0.0 < args.site_fallback_likelihood_weight < 1.0:
        raise ValueError("--site-fallback-likelihood-weight must be in (0,1)")
    if not 0.0 < args.calibration_min_robust_weight <= 1.0:
        raise ValueError("--calibration-min-robust-weight must be in (0,1]")
    if (args.calibration_shrinkage_cells < 0.0 or
            args.calibration_max_cell_weight <= 0.0 or
            args.calibration_huber_z <= 0.0 or
            args.calibration_expression_mad_cutoff <= 0.0 or
            args.cell_baseline_shrinkage_arms <= 0.0 or
            args.max_aggregate_cell_log_bf <= 0.0):
        raise ValueError("calibration tuning values must be positive")


def main_impl(args) -> int:
    validate_args(args)
    output_prefix = os.path.abspath(args.output_prefix)
    require_outputs_absent((
        output_prefix + ".arm_calls.tsv.gz",
        output_prefix + ".calibration.tsv.gz",
        output_prefix + ".uid_chromosome_flags.tsv.gz",
        output_prefix + ".donor_pair_arm_summary.tsv.gz",
        output_prefix + ".qc.tsv",
        output_prefix + ".contract.json",
    ))
    Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)
    expression_store = ExpressionStore(str(Path(output_prefix).parent))
    try:
        return run_call(args, output_prefix, expression_store)
    finally:
        expression_store.close()


def run_call(args, output_prefix: str, expression_store: ExpressionStore) -> int:
    manifests, _manifest_paths = load_manifests(args.cell_manifest)
    target_manifests = {
        key: row for key, row in manifests.items()
        if clean(row.get("donor_a")) and clean(row.get("donor_b")) and
        clean(row.get("donor_a")) != clean(row.get("donor_b")) and
        pair_from_row(row) != "NA"
    }
    excluded_nonheterotypic_cells = len(manifests) - len(target_manifests)
    target_cells = set(target_manifests)
    expression_chromosomes, _expression_paths, skipped_nonheterotypic_expression = (
        load_expression(
            args.expression, manifests, target_cells, expression_store,
            args.calibration_crossfit_folds))
    records, arm_metadata, _ase_paths, skipped_nonheterotypic_ase = load_ase(
        args.ase, manifests, target_cells, expression_store, args)

    # Expression may provide arms with zero ASE evidence.  Preserve those in the
    # output universe and stream explicit NO_DATA calls later.
    for key, chromosome in expression_chromosomes.items():
        if key in arm_metadata and chromosome_key(
                arm_metadata[key].get("chromosome")) != chromosome_key(chromosome):
            raise ValueError(f"ASE/expression chromosome mismatch for {key}")
        arm_metadata.setdefault(key, {
            "chromosome": chromosome, "arm_start": "NA", "arm_end": "NA"})
    arms_by_library: MutableMapping[str, set] = defaultdict(set)
    for library, arm in arm_metadata:
        arms_by_library[library].add(arm)
    terminal_state = "NONE"
    if not target_manifests:
        terminal_state = "PASS_NO_HETEROTYPIC_TARGETS"
    elif not records:
        terminal_state = "PASS_NO_OBSERVED_ASE"
    elif not any(
            record.model_eligible and record.n_sites > 0 and
            record.n_effective_units > 0.0 and
            clean(record.ase.get("evidence_status")).startswith("PASS")
            for record in records):
        terminal_state = "PASS_NO_CALLABLE_ASE"

    resolver = CalibrationResolver(records, expression_store, args)
    for record in records:
        score_record(record, resolver, args)
    assign_empirical_values(records, args)

    record_by_key = {record.key: record for record in records}
    calls_path = output_prefix + ".arm_calls.tsv.gz"
    calibration_path = output_prefix + ".calibration.tsv.gz"
    uid_path = output_prefix + ".uid_chromosome_flags.tsv.gz"
    pair_summary_path = output_prefix + ".donor_pair_arm_summary.tsv.gz"
    qc_path = output_prefix + ".qc.tsv"
    contract_path = output_prefix + ".contract.json"

    no_data_rows = 0
    expression_only_rows = 0
    total_output_rows = 0

    def iter_calls():
        nonlocal no_data_rows, expression_only_rows, total_output_rows
        for (library, barcode), manifest in sorted(
                target_manifests.items(), key=lambda item: (
                    natural_key(item[0][0]), natural_key(item[0][1]))):
            cell_expression = expression_store.get_cell(library, barcode)
            for arm in sorted(arms_by_library[library], key=natural_key):
                key = library, barcode, arm
                record = record_by_key.get(key)
                expr = cell_expression.get(arm)
                if record is None:
                    no_data_rows += 1
                    if expr is not None:
                        expression_only_rows += 1
                total_output_rows += 1
                yield output_row(
                    record, library, barcode, arm, manifest, expr,
                    arm_metadata[(library, arm)])

    write_tsv_atomic(calls_path, iter_calls(), CALL_FIELDS)
    calibration_count = len(resolver.resolved_cache)
    write_tsv_atomic(
        calibration_path, calibration_rows(resolver), CALIBRATION_FIELDS)
    uid_output = build_uid_rows(records, args)
    write_tsv_atomic(uid_path, uid_output, UID_FIELDS)
    pair_summary_output = build_pair_summary_rows(records, args)
    write_tsv_atomic(
        pair_summary_path, pair_summary_output, PAIR_SUMMARY_FIELDS)

    status_counts = Counter(record.call_status for record in records)
    state_counts = Counter(record.call_state for record in records)
    calibration_counts = Counter(
        record.model.status if record.model else "UNRESOLVED" for record in records)
    qc_rows = [
        {"metric": "schema_version", "value": QC_OUTPUT_SCHEMA},
        {"metric": "tool_version", "value": PROGRAM_VERSION},
        {"metric": "input_libraries", "value": len(arms_by_library)},
        {"metric": "manifest_cells", "value": len(manifests)},
        {"metric": "target_heterotypic_cells", "value": len(target_manifests)},
        {"metric": "excluded_nonheterotypic_manifest_cells",
         "value": excluded_nonheterotypic_cells},
        {"metric": "ase_observed_rows", "value": len(records)},
        {"metric": "ase_total_effective_weight", "value": format_number(
            sum(record.n_effective_units for record in records))},
        {"metric": "aggregate_eligible_rows", "value": sum(
            aggregate_qc_eligible(record, args) for record in records)},
        {"metric": "expression_rows", "value": expression_store.count},
        {"metric": "skipped_nonheterotypic_ase_rows",
         "value": skipped_nonheterotypic_ase},
        {"metric": "skipped_nonheterotypic_expression_rows",
         "value": skipped_nonheterotypic_expression},
        {"metric": "expanded_output_rows", "value": total_output_rows},
        {"metric": "no_ase_data_rows", "value": no_data_rows},
        {"metric": "expression_only_rows", "value": expression_only_rows},
        {"metric": "calibration_rows", "value": calibration_count},
        {"metric": "uid_chromosome_rows", "value": len(uid_output)},
        {"metric": "donor_pair_arm_summary_rows", "value": len(pair_summary_output)},
        {"metric": "recurrent_donor_pair_arm_events", "value": sum(
            row.get("recurrence_flag") == "TRUE" for row in pair_summary_output)},
        {"metric": "call_status_counts", "value": ";".join(
            f"{key}={value}" for key, value in sorted(status_counts.items()))},
        {"metric": "call_state_counts", "value": ";".join(
            f"{key}={value}" for key, value in sorted(state_counts.items()))},
        {"metric": "calibration_status_counts", "value": ";".join(
            f"{key}={value}" for key, value in sorted(calibration_counts.items()))},
        {"metric": "bh_scope", "value": "INVOCATION"},
        {"metric": "recommended_execution_grain", "value": "MULTI_LIBRARY_GLOBAL_CALL"},
        {"metric": "arm_universe_source", "value": "UNION_OF_ASE_AND_EXPRESSION"},
        {"metric": "terminal_state", "value": terminal_state},
        {"metric": "status", "value": (
            terminal_state if terminal_state != "NONE" else "PASS")},
    ]
    write_tsv_atomic(qc_path, qc_rows, ("metric", "value"))

    write_json_atomic(contract_path, {
        "schema_version": CONTRACT_OUTPUT_SCHEMA,
        "release": PROGRAM_VERSION,
        "input_schemas": {
            "ase": ASE_INPUT_SCHEMA,
            "cell_manifest": CELL_MANIFEST_SCHEMA,
            "expression": EXPRESSION_SCHEMA,
        },
        "output_schemas": {
            "calls": CALL_OUTPUT_SCHEMA,
            "calibration": CALIBRATION_OUTPUT_SCHEMA,
            "uid_chromosome_flags": UID_OUTPUT_SCHEMA,
            "donor_pair_arm_summary": PAIR_OUTPUT_SCHEMA,
            "qc": QC_OUTPUT_SCHEMA,
            "contract": CONTRACT_OUTPUT_SCHEMA,
        },
        "output_suffixes": {
            "calls": ".arm_calls.tsv.gz",
            "calibration": ".calibration.tsv.gz",
            "uid_chromosome_flags": ".uid_chromosome_flags.tsv.gz",
            "donor_pair_arm_summary": ".donor_pair_arm_summary.tsv.gz",
            "qc": ".qc.tsv",
            "contract": ".contract.json",
        },
        "states": list(STATES),
        "model": {
            "effective_evidence": (
                "For each molecule-arm unit with normalized donor-A support p, "
                "w=(2p-1)^2, effective_a=sum(w*p), effective_b=sum(w*(1-p)), "
                "and W=effective_a+effective_b=sum(w); hard calls and raw soft "
                "totals are coverage/QC only"),
            "ase": (
                "The normalized support y=effective_a/W is scored with an "
                "ambient-adjusted Student-t working likelihood using variance "
                "theta*(1-theta)*(rho+(1-rho)/W); state outputs are approximate "
                "working pseudo-posteriors, not fractional-count PMF posteriors"),
            "state_odds_shifts": {
                state: STATE_LOG_ODDS_SHIFT[state] for state in STATES
            },
            "state_ordering": (
                "The state log-odds shift is applied to the cross-fitted cellular "
                "donor baseline before ambient mixing; ref mapping offset is +m, "
                "alt is -m, and mixed is zero after ambient mixing"),
            "ambient_ase": (
                "ambient_a_* is consumed as the producer's full-simplex expected "
                "donor-A fraction and is not shrunk again; ambient_genotyped_mass "
                "is QC only; likelihoods are marginalized over 15 equal-probability "
                "midpoint nodes of the [0,1]-truncated normal for ambient_c"),
            "expression": (
                "Independent supportive arm-total evidence with an external-pair, "
                "held-out reference. OBSERVED_FILTERED_COUNTS requires raw integer "
                "counts and uses conservative ambient_c attenuation/variance "
                "inflation; UPSTREAM_AMBIENT_CORRECTED requires an explicitly "
                "corrected count type and is not re-attenuated; minimum gene "
                "support uses per-cell nonzero_genes_on_arm, never the arm-wide "
                "mapped annotation-feature count"),
            "ase_calibration_hierarchy": list(ASE_CALIBRATION_LEVELS),
            "expression_calibration_hierarchy": list(EXTERNAL_EXPRESSION_LEVELS),
            "empirical_null_hierarchy": list(EMPIRICAL_NULL_LEVELS),
            "cross_fitting": (
                "Every query excludes its chromosome and deterministic UID/cell "
                "fold from nuisance calibration; empirical ranks compare external "
                "donor-pair cells in that same held-out fold"),
            "rho_calibration": (
                "Baseline centers use one capped genome-wide contribution per cell; "
                "rho uses robust residuals of individual held-out cell-arm observations"),
            "posterior_event_prior": args.event_prior,
            "empirical_p": (
                "state-specific plus-one upper-tail tests against held-out "
                "external-pair, effective-weight/depth- and ambient-matched nulls"),
            "multiple_testing": (
                "Benjamini-Hochberg across all cell-arm-event-state hypotheses "
                "in this invocation is retained as a deliberately conservative "
                "cell-level diagnostic; pair-primary BH is a separate family"),
            "uid_arm_p": (
                "state-specific Bonferroni partial-conjunction tests with "
                "configurable r, followed by BH over UID-arm-event-state hypotheses"),
            "primary_recurrence": (
                "state-specific Fisher partial-conjunction on p_(r)..p_(n), "
                "with deterministic score-blind cell/UID caps and BH over all "
                "library-group-pair-arm-event-state hypotheses; library is "
                "part of the inferential key because default expression "
                "cluster labels are not comparable across libraries. This is "
                "an approximate working FDR under conditional independence, "
                "not an exact/conservative dependence-robust guarantee"),
            "sex_chromosomes": (
                "Non-autosomal rows are exploratory only and cannot PASS or enter "
                "UID/donor-pair recurrence inference"),
        },
        "thresholds": {
            "min_call_sites": args.min_call_sites,
            "min_call_effective_weight": args.min_call_effective_weight,
            "min_calibration_sites": args.min_calibration_sites,
            "min_calibration_effective_weight": (
                args.min_calibration_effective_weight),
            "min_orientation_calibration_effective_weight": (
                args.min_orientation_calibration_effective_weight),
            "min_fallback_calibration_cells": (
                args.min_fallback_calibration_cells),
            "calibration_crossfit_folds": args.calibration_crossfit_folds,
            "min_event_posterior": args.min_event_posterior,
            "min_balanced_posterior": args.min_balanced_posterior,
            "max_event_q": args.max_event_q,
            "min_empirical_null_cells": args.min_empirical_null_cells,
            "absolute_min_null_cells": args.absolute_min_null_cells,
            "empirical_depth_fold": args.empirical_depth_fold,
            "empirical_ambient_window": args.empirical_ambient_window,
            "expression_weight": args.expression_weight,
            "max_expression_log_bf": args.max_expression_log_bf,
            "min_uid_cells_per_arm": args.min_uid_cells_per_arm,
            "min_uid_arm_posterior": args.min_uid_arm_posterior,
            "min_uid_state_concordance": args.min_uid_state_concordance,
            "max_whole_chromosome_q": args.max_whole_chromosome_q,
            "min_pair_recurrence_cells": args.min_pair_recurrence_cells,
            "max_pair_test_cells": args.max_pair_test_cells,
            "max_pair_cells_per_uid": args.max_pair_cells_per_uid,
            "min_pair_arm_posterior": args.min_pair_arm_posterior,
            "min_pair_state_concordance": args.min_pair_state_concordance,
            "max_pair_arm_q": args.max_pair_arm_q,
        },
        "limitations": [
            "Expression is supportive total-copy evidence, not donor-directional evidence.",
            "No gene-level ambient decontamination is claimed without a matching gene profile.",
            "Eligibility/status metadata are not independent likelihood terms.",
            "Empirical p-value resolution is limited by the matched calibration pool size.",
            "The calibrated donor baseline and antisymmetric orientation offset can still include biological allelic-expression and residual mapping effects that are not separately identifiable.",
            "An arm absent from both ASE and expression cannot be inferred; the emitted arm universe is their union.",
            "A multi-library invocation is preferred for study-wide BH control and cross-library UID aggregation.",
            "ASE records are held in memory once; optional expression rows are disk-backed in a temporary SQLite index.",
            "The approximate Fisher/BH pair result requires leave-one-donor-pair and UID-block negative-control QQ/type-I validation before its q-value is interpreted as calibrated FDR; the caller does not impute evidence when external p-values are scarce.",
        ],
        "terminal_state": terminal_state,
        "status": terminal_state if terminal_state != "NONE" else "PASS",
    })
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Call donor-directed chromosome-arm CNVs from calibrated ASE and "
            "optional arm-expression evidence."))
    parser.add_argument("--version", action="version", version=f"%(prog)s {PROGRAM_VERSION}")
    parser.add_argument("--ase", "--ase-evidence", nargs="+", required=True,
                        help="One or more tetra_arm_ase evidence TSV.GZ files")
    parser.add_argument("--cell-manifest", nargs="+", required=True,
                        help="Matching prepared cell-manifest TSV.GZ files")
    parser.add_argument("--expression", nargs="*", default=[],
                        help="Optional matching arm-expression TSV.GZ files")
    parser.add_argument("--output-prefix", required=True,
                        help="Prefix for atomic calls/calibration/UID/QC/contract outputs")

    parser.add_argument("--min-call-sites", type=int, default=3)
    parser.add_argument(
        "--min-call-effective-weight", "--min-call-molecules",
        dest="min_call_effective_weight", type=float, default=8.0,
        help="Minimum summed effective ASE weight for an individual call")
    parser.add_argument("--min-calibration-sites", type=int, default=2)
    parser.add_argument(
        "--min-calibration-effective-weight", "--min-calibration-molecules",
        dest="min_calibration_effective_weight", type=float, default=4.0,
        help="Minimum summed effective ASE weight for a calibration row")
    parser.add_argument(
        "--min-orientation-calibration-effective-weight", type=float,
        default=20.0,
        help=("Minimum capped robust effective weight required to calibrate "
              "one orientation"))
    parser.add_argument("--min-fallback-calibration-cells", type=int, default=12)
    parser.add_argument("--calibration-shrinkage-cells", type=float, default=20.0)
    parser.add_argument("--calibration-max-cell-weight", type=float, default=50.0)
    parser.add_argument("--calibration-huber-z", type=float, default=2.5)
    parser.add_argument("--calibration-min-robust-weight", type=float, default=0.25)
    parser.add_argument("--calibration-expression-mad-cutoff", type=float, default=4.0)
    parser.add_argument("--calibration-crossfit-folds", type=int, default=5)
    parser.add_argument("--min-cell-baseline-arms", type=int, default=4)
    parser.add_argument("--cell-baseline-shrinkage-arms", type=float, default=8.0)
    parser.add_argument("--default-rho", type=float, default=0.02)
    parser.add_argument("--min-rho", type=float, default=0.001)
    parser.add_argument("--max-rho", type=float, default=0.25)

    parser.add_argument("--event-prior", type=float, default=0.04,
                        help="Total prior probability assigned across four event states")
    parser.add_argument("--min-event-posterior", type=float, default=0.90)
    parser.add_argument("--min-balanced-posterior", type=float, default=0.80)
    parser.add_argument("--max-event-q", type=float, default=0.05)
    parser.add_argument("--min-empirical-null-cells", type=int, default=20)
    parser.add_argument("--absolute-min-null-cells", type=int, default=5)
    parser.add_argument("--empirical-depth-fold", type=float, default=2.0)
    parser.add_argument("--empirical-ambient-window", type=float, default=0.10)

    parser.add_argument("--min-expression-counts", type=float, default=100.0)
    parser.add_argument("--min-expression-genes", type=int, default=3)
    parser.add_argument("--min-expression-sigma", type=float, default=0.15)
    parser.add_argument("--max-expression-sigma", type=float, default=2.0)
    parser.add_argument("--expression-weight", type=float, default=0.50)
    parser.add_argument("--max-expression-log-bf", type=float, default=4.0)
    parser.add_argument("--site-fallback-likelihood-weight", type=float, default=0.50)

    parser.add_argument("--max-qname-fallback-fraction", type=float, default=0.50)
    parser.add_argument("--min-ambient-genotyped-mass", type=float, default=0.50)
    parser.add_argument("--min-uid-cells-per-arm", type=int, default=1)
    parser.add_argument("--min-uid-arm-posterior", type=float, default=0.80)
    parser.add_argument("--min-uid-state-concordance", type=float, default=0.67)
    parser.add_argument("--max-aggregate-cell-log-bf", type=float, default=8.0)
    parser.add_argument("--max-whole-chromosome-q", type=float, default=0.05)
    parser.add_argument("--min-pair-recurrence-cells", type=int, default=2)
    parser.add_argument(
        "--max-pair-test-cells", type=int, default=500,
        help="Score-blind cap on cells entering one Fisher pair-arm test")
    parser.add_argument(
        "--max-pair-cells-per-uid", type=int, default=100,
        help="Score-blind within-UID cap for one Fisher pair-arm test")
    parser.add_argument("--min-pair-state-concordance", type=float, default=0.60)
    parser.add_argument("--min-pair-cell-posterior", type=float, default=0.50)
    parser.add_argument("--min-pair-arm-posterior", type=float, default=0.90)
    parser.add_argument("--max-pair-arm-q", type=float, default=0.05)
    return parser


def main() -> int:
    try:
        return main_impl(build_parser().parse_args())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
