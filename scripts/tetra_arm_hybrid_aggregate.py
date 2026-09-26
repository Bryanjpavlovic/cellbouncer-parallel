#!/usr/bin/env python3
"""Globally gather hybrid chromosome scores and apply invocation-wide FDR."""

from __future__ import annotations

import argparse
import csv
import heapq
import json
import math
import os
import sqlite3
import struct
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np

from tetra_arm_common import (
    HYBRID_CALIBRATION_SCHEMA,
    HYBRID_CALL_SCHEMA,
    HYBRID_COMPONENT_SCHEMA,
    HYBRID_CONTRACT_SCHEMA,
    HYBRID_PAIR_SCHEMA,
    HYBRID_QC_SCHEMA,
    HYBRID_SCORE_SCHEMA,
    HYBRID_UID_SCHEMA,
    clean,
    file_record,
    finite_float,
    natural_key,
    read_tsv,
    require_file,
    require_outputs_absent,
    write_json_atomic,
    write_tsv_atomic,
)
from tetra_arm_hybrid_common import (
    ASE_DIRECTIONS,
    AXES_TO_STATE,
    COPY_STATES,
    EXACT_STATES,
    PROGRAM_VERSION,
    bh_adjust,
    by_adjust,
    canonical_chromosome,
    classify_pq_relationship,
    format_number,
    join_flags,
    present_identifier,
    resolve_evidence,
    split_flags,
    validate_hybrid_table,
)


MODEL_TASK_HEADER = (
    "task_index", "chromosome", "logical_arms", "shard_manifest",
    "output_directory", "qc", "contract")

FINAL_ADDITIONAL_FIELDS = [
    "ase_q_DONOR_A_DEPLETED", "ase_q_DONOR_A_ENRICHED",
    "expression_q_LOSS", "expression_q_GAIN",
    "hybrid_conjunction_q_DONOR_A_LOSS",
    "hybrid_conjunction_q_DONOR_B_LOSS",
    "hybrid_conjunction_q_DONOR_A_GAIN",
    "hybrid_conjunction_q_DONOR_B_GAIN",
    "hybrid_conjunction_by_DONOR_A_LOSS",
    "hybrid_conjunction_by_DONOR_B_LOSS",
    "hybrid_conjunction_by_DONOR_A_GAIN",
    "hybrid_conjunction_by_DONOR_B_GAIN",
    "evidence_class", "copy_state", "donor_origin", "resolved_state",
    "confidence_tier", "call_state", "call_status",
    "rna_dosage_interpretation", "orthogonal_validation_status",
    "bh_scope_ase", "bh_scope_expression", "bh_scope_hybrid",
    "score_schema_version", "schema_version",
]

UID_FIELDS = [
    "uid", "donor_pair", "chromosome", "libraries", "n_cells",
    "distinct_uid_blocks", "support_basis", "p_arm", "p_evidence_class",
    "p_resolved_state",
    "p_expression_support_cells", "p_ase_support_cells",
    "p_concordant_cells", "p_min_expression_q", "p_min_ase_q",
    "p_min_hybrid_q", "q_arm", "q_evidence_class", "q_resolved_state",
    "q_expression_support_cells", "q_ase_support_cells",
    "q_concordant_cells", "q_min_expression_q", "q_min_ase_q",
    "q_min_hybrid_q", "pq_relationship", "summary_status", "schema_version",
]

PAIR_FIELDS = [
    "donor_pair", "chromosome", "arm", "libraries", "calibration_groups",
    "total_cells", "distinct_uid_blocks", "support_basis",
    "expression_support_cells",
    "ase_support_cells", "concordant_cells", "expression_only_cells",
    "ase_only_cells", "discordant_cells", "balanced_cells",
    "expression_support_uid_blocks", "ase_support_uid_blocks",
    "concordant_uid_blocks", "expression_only_uid_blocks",
    "ase_only_uid_blocks", "discordant_uid_blocks", "balanced_uid_blocks",
    "best_resolved_state", "resolved_state_counts", "min_expression_q",
    "min_ase_q", "min_hybrid_q", "pairwide_identifiability_flag",
    "summary_status", "schema_version",
]


def read_model_tasks(path: str) -> list[dict[str, str]]:
    target = require_file(path, "hybrid model task manifest")
    with open(target, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != MODEL_TASK_HEADER:
            raise ValueError(
                f"model task manifest header must be {MODEL_TASK_HEADER}: {target}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"hybrid model task manifest has no chromosomes: {target}")
    seen_tasks: set[int] = set()
    seen_chromosomes: set[str] = set()
    for row in rows:
        task = int(clean(row.get("task_index")))
        chromosome = canonical_chromosome(row.get("chromosome"))
        if task in seen_tasks or chromosome in seen_chromosomes:
            raise ValueError("duplicate task index or chromosome in model manifest")
        for field in ("shard_manifest", "output_directory", "qc", "contract"):
            if not os.path.isabs(row.get(field, "")):
                raise ValueError(f"model task {field} paths must be absolute")
        if not clean(row.get("logical_arms")):
            raise ValueError("model task logical_arms must be nonempty")
        seen_tasks.add(task)
        seen_chromosomes.add(chromosome)
        row["chromosome"] = chromosome
    if sorted(seen_tasks) != list(range(len(rows))):
        raise ValueError("model task indices must be contiguous from zero")
    return sorted(rows, key=lambda row: int(row["task_index"]))


def model_paths(task: Mapping[str, str]) -> dict[str, str]:
    chromosome = task["chromosome"]
    prefix = os.path.join(
        task["output_directory"], f"tetra_arm_hybrid_{chromosome}")
    return {
        "scores": prefix + ".scores.tsv.gz",
        "components": prefix + ".expression_components.tsv.gz",
        "calibration": prefix + ".calibration.tsv.gz",
        "qc": task["qc"], "contract": task["contract"],
    }


def require_contract(path: str, chromosome: str) -> dict:
    target = require_file(path, "hybrid chromosome contract")
    with open(target, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if (not isinstance(payload, dict)
            or payload.get("schema_version") != "tetra_arm_hybrid_model_contract_v1"
            or canonical_chromosome(payload.get("chromosome")) != chromosome
            or clean(payload.get("status")).upper() != "PASS"):
        raise ValueError(f"invalid hybrid chromosome contract: {target}")
    schemas = payload.get("output_schemas", {})
    if (not isinstance(schemas, dict)
            or schemas.get("scores") != HYBRID_SCORE_SCHEMA
            or schemas.get("expression_components") != HYBRID_COMPONENT_SCHEMA
            or schemas.get("calibration") != HYBRID_CALIBRATION_SCHEMA):
        raise ValueError(f"hybrid chromosome schemas do not match: {target}")
    return payload


def selected_q(row: Mapping[str, object], namespace: str,
               state: str) -> float:
    if namespace == "expression" and state in COPY_STATES:
        return finite_float(row.get(f"expression_q_{state}"))
    if namespace == "ase" and state in ASE_DIRECTIONS:
        return finite_float(row.get(f"ase_q_{state}"))
    if namespace == "hybrid" and state in EXACT_STATES:
        return finite_float(row.get(f"hybrid_conjunction_q_{state}"))
    return math.nan


def model_selected_qs(row: Mapping[str, object]) -> tuple[float, float, float]:
    """Return q-values only for directions selected by the chromosome model."""
    expression_state = clean(row.get("expression_best_copy_state"))
    ase_direction = clean(row.get("ase_direction"))
    exact_state = AXES_TO_STATE.get((expression_state, ase_direction), "NO_CALL")
    return (
        selected_q(row, "expression", expression_state),
        selected_q(row, "ase", ase_direction),
        selected_q(row, "hybrid", exact_state),
    )


def minimum_finite(values: Iterable[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return min(finite) if finite else math.nan


def apply_family(rows: list[dict[str, object]], p_fields: Sequence[str],
                 q_fields: Sequence[str], by_fields: Sequence[str] | None = None
                 ) -> int:
    references: list[tuple[int, int]] = []
    values: list[float] = []
    for row_index, row in enumerate(rows):
        for field_index, p_field in enumerate(p_fields):
            value = finite_float(row.get(p_field))
            if math.isfinite(value):
                references.append((row_index, field_index))
                values.append(value)
    q_values = bh_adjust(values)
    by_values = by_adjust(values) if by_fields is not None else []
    for row in rows:
        for field in q_fields:
            row[field] = "NA"
        if by_fields is not None:
            for field in by_fields:
                row[field] = "NA"
    for index, (row_index, field_index) in enumerate(references):
        rows[row_index][q_fields[field_index]] = format_number(q_values[index])
        if by_fields is not None:
            rows[row_index][by_fields[field_index]] = format_number(by_values[index])
    return len(values)


FAMILY_P_FIELDS = {
    "ase": tuple(f"ase_p_{direction}" for direction in ASE_DIRECTIONS),
    "expression": tuple(f"expression_p_{state}" for state in COPY_STATES),
    "hybrid": tuple(
        f"hybrid_conjunction_p_{state}" for state in EXACT_STATES),
}
FAMILY_Q_FIELDS = {
    "ase": tuple(f"ase_q_{direction}" for direction in ASE_DIRECTIONS),
    "expression": tuple(f"expression_q_{state}" for state in COPY_STATES),
    "hybrid": tuple(
        f"hybrid_conjunction_q_{state}" for state in EXACT_STATES),
}
HYPOTHESIS_STRUCT = struct.Struct("<dQB")
HYPOTHESIS_DTYPE = np.dtype([
    ("p", "<f8"), ("location", "<u8"), ("field", "u1")])
ADJUSTMENT_CHUNK = 1_000_000


def family_row_eligible(row: Mapping[str, object], family: str) -> bool:
    expression_pass = (
        clean(row.get("expression_status")) == "PASS"
        and clean(row.get("expression_fold_replication_status"))
        in {"CONCORDANT", "CONCORDANT_BALANCED"})
    ase_pass = clean(row.get("ase_status")) == "PASS"
    if family == "ase":
        return ase_pass
    if family == "expression":
        return expression_pass
    if family == "hybrid":
        return ase_pass and expression_pass
    raise ValueError(f"unknown FDR family: {family}")


def append_hypotheses(
        row: Mapping[str, object], location: int,
        handles: Mapping[str, object], counts: dict[str, int]) -> None:
    """Append only fully calibrated hypotheses to narrow disk streams."""
    for family, p_fields in FAMILY_P_FIELDS.items():
        if not family_row_eligible(row, family):
            continue
        values = [finite_float(row.get(field)) for field in p_fields]
        if any(not math.isfinite(value) or not 0.0 <= value <= 1.0
               for value in values):
            raise ValueError(
                f"PASS {family} row has missing/invalid p-values at "
                f"{clean(row.get('library'))}/{clean(row.get('barcode'))}/"
                f"{clean(row.get('arm'))}")
        for field_index, value in enumerate(values):
            handles[family].write(HYPOTHESIS_STRUCT.pack(
                value, location, field_index))
            counts[family] += 1


def initialize_q_map(path: str, rows: int, fields: int) -> np.memmap | None:
    if rows == 0:
        return None
    result = np.memmap(
        path, dtype="<f8", mode="w+", shape=(rows, fields))
    for start in range(0, rows, ADJUSTMENT_CHUNK):
        result[start:min(rows, start + ADJUSTMENT_CHUNK), :] = np.nan
    result.flush()
    return result


def external_bh_map(
        hypothesis_path: str, q_path: str, hypothesis_count: int,
        total_rows: int, field_count: int) -> np.memmap | None:
    """Sort a narrow on-disk family and write a row-addressable q memmap."""
    q_map = initialize_q_map(q_path, total_rows, field_count)
    if hypothesis_count == 0:
        return q_map
    observed_size = os.path.getsize(hypothesis_path)
    expected_size = hypothesis_count * HYPOTHESIS_DTYPE.itemsize
    if observed_size != expected_size:
        raise ValueError(
            f"truncated external hypothesis stream: {hypothesis_path}")
    records = np.memmap(
        hypothesis_path, dtype=HYPOTHESIS_DTYPE, mode="r+",
        shape=(hypothesis_count,))
    records.sort(order=["p", "location", "field"], kind="quicksort")
    running = 1.0
    end = hypothesis_count
    while end:
        start = max(0, end - ADJUSTMENT_CHUNK)
        probabilities = np.asarray(records["p"][start:end], dtype=np.float64)
        ranks = np.arange(start + 1, end + 1, dtype=np.float64)
        adjusted = np.minimum(1.0, probabilities * hypothesis_count / ranks)
        adjusted = np.minimum(adjusted, running)
        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
        running = float(adjusted[0])
        locations = np.asarray(records["location"][start:end], dtype=np.uint64)
        fields = np.asarray(records["field"][start:end], dtype=np.uint8)
        if np.any(locations >= total_rows) or np.any(fields >= field_count):
            raise ValueError("external hypothesis address is outside q-map bounds")
        q_map[locations, fields] = adjusted
        end = start
    q_map.flush()
    records.flush()
    del records
    return q_map


def harmonic_number(count: int) -> float:
    total = 0.0
    for start in range(1, count + 1, ADJUSTMENT_CHUNK):
        stop = min(count + 1, start + ADJUSTMENT_CHUNK)
        total += float(np.sum(1.0 / np.arange(
            start, stop, dtype=np.float64), dtype=np.float64))
    return total


def attach_adjustments(
        row: dict[str, str], location: int,
        q_maps: Mapping[str, np.memmap | None], hybrid_harmonic: float,
        ) -> dict[str, str]:
    for family, q_fields in FAMILY_Q_FIELDS.items():
        q_map = q_maps.get(family)
        for field_index, field in enumerate(q_fields):
            value = (float(q_map[location, field_index])
                     if q_map is not None else math.nan)
            row[field] = format_number(value)
            if family == "hybrid":
                row[field.replace("_q_", "_by_")] = format_number(
                    min(1.0, value * hybrid_harmonic)
                    if math.isfinite(value) else math.nan)
    return row


def score_sort_key(row: Mapping[str, object]) -> tuple[object, ...]:
    return (
        natural_key(row.get("library", "")),
        natural_key(row.get("barcode", "")),
        natural_key(row.get("chromosome", "")),
        natural_key(row.get("arm", "")),
    )


def iter_adjusted_scores(
        tasks: Sequence[Mapping[str, str]], offsets: Mapping[int, int],
        q_maps: Mapping[str, np.memmap | None], hybrid_harmonic: float,
        ) -> Iterable[dict[str, str]]:
    def source(task_index: int, task: Mapping[str, str]):
        path = model_paths(task)["scores"]
        for row_index, row in enumerate(read_tsv(path)):
            yield (score_sort_key(row), task_index, row_index, row)

    sources = [source(index, task) for index, task in enumerate(tasks)]
    previous: tuple[object, ...] | None = None
    for key, task_index, row_index, row in heapq.merge(
            *sources, key=lambda value: (value[0], value[1], value[2])):
        if previous is not None and key <= previous:
            raise ValueError("global score rows are duplicate or unsorted")
        previous = key
        location = offsets[task_index] + row_index
        yield attach_adjustments(row, location, q_maps, hybrid_harmonic)


def independent_block(row: Mapping[str, object]) -> str:
    uid = present_identifier(row.get("uid"))
    if uid:
        return f"UID:{uid}"
    pair = present_identifier(row.get("donor_pair"))
    return f"PAIR:{pair}" if pair else (
        f"CELL:{clean(row.get('library'))}:{clean(row.get('barcode'))}")


def preliminary_axes(row: Mapping[str, object], threshold: float) -> dict[str, object]:
    # The chromosome model chooses direction using replicated gene folds and
    # ASE likelihoods.  Aggregation may threshold that selected hypothesis but
    # must never replace it with whichever competing q-value happens to be
    # smallest after invocation-wide adjustment.
    model_expression_state = clean(
        row.get("expression_best_copy_state"))
    expression_significant = bool(
        model_expression_state in COPY_STATES
        and clean(row.get("expression_status")) == "PASS"
        and clean(row.get("expression_fold_replication_status")) == "CONCORDANT"
        and selected_q(row, "expression", model_expression_state) <= threshold)
    expression_state = (
        model_expression_state if expression_significant else "BALANCED")
    model_ase_direction = clean(row.get("ase_direction"))
    ase_significant = bool(
        model_ase_direction in ASE_DIRECTIONS
        and clean(row.get("ase_status")) == "PASS"
        and selected_q(row, "ase", model_ase_direction) <= threshold)
    ase_direction = model_ase_direction if ase_significant else "BALANCED"
    expected_state = AXES_TO_STATE.get((expression_state, ase_direction))
    conjunction_significant = bool(
        expected_state in EXACT_STATES
        and selected_q(row, "hybrid", expected_state) <= threshold)
    conjunction_state = expected_state if conjunction_significant else "NO_CALL"
    return {
        "expression_state": expression_state,
        "expression_significant": expression_significant,
        "ase_direction": ase_direction,
        "ase_significant": ase_significant,
        "conjunction_state": conjunction_state,
        "conjunction_significant": conjunction_significant,
    }


def counter_choice(counter: Counter, allowed: set[str], default: str) -> str:
    candidates = [(value, count) for value, count in counter.items()
                  if value in allowed]
    if not candidates:
        return default
    maximum = max(count for _value, count in candidates)
    winners = [value for value, count in candidates if count == maximum]
    # A tied cell/UID vote is unresolved evidence, never an event selected by
    # lexical ordering of state names.
    if len(winners) != 1:
        return default
    return winners[0]


def pairwide_groups_from_stream(
        rows: Iterable[dict[str, str]], max_q: float, fraction: float,
        minimum_blocks: int) -> set[tuple[str, str, str, str]]:
    block_votes: dict[tuple[str, str, str, str], Counter] = defaultdict(Counter)
    for row in rows:
        pair = present_identifier(row.get("donor_pair"))
        if not pair:
            continue
        axes = preliminary_axes(row, threshold=max_q)
        state = clean(axes.get("conjunction_state"))
        if state not in EXACT_STATES:
            state = "NO_CALL"
        key = (pair, clean(row.get("chromosome")), clean(row.get("arm")),
               independent_block(row))
        block_votes[key][state] += 1
    grouped: dict[tuple[str, str, str], Counter] = defaultdict(Counter)
    for (pair, chromosome, arm, _block), votes in block_votes.items():
        selected = counter_choice(
            votes, set(EXACT_STATES) | {"NO_CALL"}, "NO_CALL")
        grouped[(pair, chromosome, arm)][selected] += 1
    flagged: set[tuple[str, str, str, str]] = set()
    for group, votes in grouped.items():
        total_blocks = sum(votes.values())
        if total_blocks < minimum_blocks:
            continue
        state = counter_choice(
            votes, set(EXACT_STATES) | {"NO_CALL", "BALANCED"}, "NO_CALL")
        if state in EXACT_STATES and votes[state] / total_blocks >= fraction:
            flagged.add((*group, state))
    return flagged


def update_block_summary(
        summaries: dict[tuple[str, str, str, str], dict[str, object]],
        row: Mapping[str, object]) -> None:
    pair = present_identifier(row.get("donor_pair"))
    if not pair:
        return
    chromosome, arm = clean(row.get("chromosome")), clean(row.get("arm"))
    block = independent_block(row)
    key = (pair, chromosome, arm, block)
    value = summaries.get(key)
    if value is None:
        value = {
            "pair": pair, "chromosome": chromosome, "arm": arm,
            "block": block, "uid": present_identifier(row.get("uid")),
            "libraries": set(), "groups": set(), "total_cells": 0,
            "classes": Counter(), "states": Counter(),
            "min_expression_q": math.nan, "min_ase_q": math.nan,
            "min_hybrid_q": math.nan, "pairwide": False,
        }
        summaries[key] = value
    value["libraries"].add(clean(row.get("library")))
    value["groups"].add(clean(row.get("calibration_group")))
    value["total_cells"] += 1
    value["classes"][clean(row.get("evidence_class")) or "NO_DATA"] += 1
    value["states"][clean(row.get("resolved_state")) or "NO_CALL"] += 1
    selected_values = model_selected_qs(row)
    for candidate, field in zip(selected_values, (
            "min_expression_q", "min_ase_q", "min_hybrid_q")):
        previous = float(value[field])
        if math.isfinite(candidate) and (
                not math.isfinite(previous) or candidate < previous):
            value[field] = candidate
    if "PAIRWIDE_EVENT_NOT_IDENTIFIABLE_FROM_RNA" in split_flags(
            row.get("confounding_flags")):
        value["pairwide"] = True


def block_effective_class(value: Mapping[str, object]) -> str:
    return counter_choice(
        value["classes"], {
            "CONCORDANT_BOTH", "EXPRESSION_ONLY", "ASE_ONLY", "DISCORDANT",
            "EXPRESSION_OUTLIER", "BALANCED", "INSUFFICIENT_EVIDENCE",
        }, "INSUFFICIENT_EVIDENCE")


def block_effective_state(value: Mapping[str, object]) -> str:
    return counter_choice(
        value["states"], set(EXACT_STATES) | {"NO_CALL", "BALANCED"},
        "NO_CALL")


def uid_summary_from_blocks(
        blocks: Mapping[tuple[str, str, str, str], Mapping[str, object]],
        uid_cell_counts: Mapping[tuple[str, str, str], int],
        ) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[Mapping[str, object]]] = defaultdict(list)
    for value in blocks.values():
        uid = clean(value.get("uid"))
        if uid:
            grouped[(uid, clean(value.get("pair")),
                     clean(value.get("chromosome")))].append(value)
    output: list[dict[str, object]] = []
    for key, values in sorted(
            grouped.items(), key=lambda item: tuple(
                natural_key(part) for part in item[0])):
        uid, pair, chromosome = key
        summary: dict[str, object] = {
            "uid": uid, "donor_pair": pair, "chromosome": chromosome,
            "libraries": ";".join(sorted({
                library for value in values for library in value["libraries"]
            }, key=natural_key)),
            "n_cells": uid_cell_counts.get(key, 0),
            "distinct_uid_blocks": 1,
            "support_basis": "WITHIN_UID_CLONE_REDUCED_ONCE",
        }
        selected: dict[str, str] = {}
        for side in ("p", "q"):
            side_values = [value for value in values
                           if clean(value.get("arm")).lower().endswith(side)]
            classes = Counter()
            states = Counter()
            for value in side_values:
                classes.update(value["classes"])
                states.update(value["states"])
            evidence_class = counter_choice(
                classes, {
                    "CONCORDANT_BOTH", "EXPRESSION_ONLY", "ASE_ONLY",
                    "DISCORDANT", "EXPRESSION_OUTLIER", "BALANCED",
                    "INSUFFICIENT_EVIDENCE",
                }, "NO_DATA")
            resolved_state = counter_choice(
                states, set(EXACT_STATES) | {"NO_CALL", "BALANCED"},
                "NO_CALL")
            selected[side] = resolved_state
            summary.update({
                f"{side}_arm": clean(side_values[0].get("arm"))
                    if side_values else "NA",
                f"{side}_evidence_class": evidence_class,
                f"{side}_resolved_state": resolved_state,
                f"{side}_expression_support_cells": (
                    classes["CONCORDANT_BOTH"] + classes["EXPRESSION_ONLY"]),
                f"{side}_ase_support_cells": (
                    classes["CONCORDANT_BOTH"] + classes["ASE_ONLY"]),
                f"{side}_concordant_cells": classes["CONCORDANT_BOTH"],
                f"{side}_min_expression_q": format_number(minimum_finite(
                    float(value["min_expression_q"]) for value in side_values)),
                f"{side}_min_ase_q": format_number(minimum_finite(
                    float(value["min_ase_q"]) for value in side_values)),
                f"{side}_min_hybrid_q": format_number(minimum_finite(
                    float(value["min_hybrid_q"]) for value in side_values)),
            })
        relationship = classify_pq_relationship(
            selected.get("p", "NO_CALL"), selected.get("q", "NO_CALL"))
        summary["pq_relationship"] = relationship
        summary["summary_status"] = (
            "RNA_CANDIDATE" if relationship != "NO_RESOLVED_EVENT"
            else "NO_CALL")
        summary["schema_version"] = HYBRID_UID_SCHEMA
        output.append(summary)
    return output


def pair_summary_from_blocks(
        blocks: Mapping[tuple[str, str, str, str], Mapping[str, object]],
        ) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[Mapping[str, object]]] = defaultdict(list)
    for value in blocks.values():
        grouped[(clean(value.get("pair")), clean(value.get("chromosome")),
                 clean(value.get("arm")))].append(value)
    output: list[dict[str, object]] = []
    for (pair, chromosome, arm), values in sorted(
            grouped.items(), key=lambda item: tuple(
                natural_key(part) for part in item[0])):
        cell_classes = Counter()
        block_classes = Counter()
        block_states = Counter()
        for value in values:
            cell_classes.update(value["classes"])
            block_classes[block_effective_class(value)] += 1
            block_states[block_effective_state(value)] += 1
        exact_block_states = Counter({
            state: count for state, count in block_states.items()
            if state in EXACT_STATES})
        best_state = counter_choice(
            block_states, set(EXACT_STATES) | {"NO_CALL", "BALANCED"},
            "NO_CALL")
        output.append({
            "donor_pair": pair, "chromosome": chromosome, "arm": arm,
            "libraries": ";".join(sorted({
                library for value in values for library in value["libraries"]
            }, key=natural_key)),
            "calibration_groups": ";".join(sorted({
                group for value in values for group in value["groups"]
            }, key=natural_key)),
            "total_cells": sum(int(value["total_cells"]) for value in values),
            "distinct_uid_blocks": len(values),
            "support_basis": "UID_BLOCK_FIRST",
            "expression_support_cells": (
                cell_classes["CONCORDANT_BOTH"]
                + cell_classes["EXPRESSION_ONLY"]),
            "ase_support_cells": (
                cell_classes["CONCORDANT_BOTH"] + cell_classes["ASE_ONLY"]),
            "concordant_cells": cell_classes["CONCORDANT_BOTH"],
            "expression_only_cells": cell_classes["EXPRESSION_ONLY"],
            "ase_only_cells": cell_classes["ASE_ONLY"],
            "discordant_cells": (
                cell_classes["DISCORDANT"] + cell_classes["EXPRESSION_OUTLIER"]),
            "balanced_cells": cell_classes["BALANCED"],
            "expression_support_uid_blocks": (
                block_classes["CONCORDANT_BOTH"]
                + block_classes["EXPRESSION_ONLY"]),
            "ase_support_uid_blocks": (
                block_classes["CONCORDANT_BOTH"] + block_classes["ASE_ONLY"]),
            "concordant_uid_blocks": block_classes["CONCORDANT_BOTH"],
            "expression_only_uid_blocks": block_classes["EXPRESSION_ONLY"],
            "ase_only_uid_blocks": block_classes["ASE_ONLY"],
            "discordant_uid_blocks": (
                block_classes["DISCORDANT"]
                + block_classes["EXPRESSION_OUTLIER"]),
            "balanced_uid_blocks": block_classes["BALANCED"],
            "best_resolved_state": best_state,
            "resolved_state_counts": ";".join(
                f"{state}={count}" for state, count in sorted(
                    exact_block_states.items(),
                    key=lambda item: natural_key(item[0]))) or "NONE",
            "min_expression_q": format_number(minimum_finite(
                float(value["min_expression_q"]) for value in values)),
            "min_ase_q": format_number(minimum_finite(
                float(value["min_ase_q"]) for value in values)),
            "min_hybrid_q": format_number(minimum_finite(
                float(value["min_hybrid_q"]) for value in values)),
            "pairwide_identifiability_flag": (
                "TRUE" if any(bool(value["pairwide"]) for value in values)
                else "FALSE"),
            "summary_status": "RNA_CANDIDATE"
                if best_state in EXACT_STATES else "NO_CALL",
            "schema_version": HYBRID_PAIR_SCHEMA,
        })
    return output


BLOCK_CLASSES = (
    "CONCORDANT_BOTH", "EXPRESSION_ONLY", "ASE_ONLY", "DISCORDANT",
    "EXPRESSION_OUTLIER", "BALANCED", "INSUFFICIENT_EVIDENCE",
)
BLOCK_STATES = (*EXACT_STATES, "NO_CALL", "BALANCED")
BLOCK_CLASS_COLUMNS = {
    value: f"class_{value.lower()}" for value in BLOCK_CLASSES}
BLOCK_STATE_COLUMNS = {
    value: f"state_{value.lower()}" for value in BLOCK_STATES}
SQLITE_BATCH_ROWS = 20_000


def natural_compare(left: str, right: str) -> int:
    left_key, right_key = natural_key(left), natural_key(right)
    return (left_key > right_key) - (left_key < right_key)


def configure_disk_database(connection: sqlite3.Connection) -> None:
    connection.create_collation("NATURALKEY", natural_compare)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("PRAGMA cache_size=-65536")


def pairwide_groups_from_stream_disk(
        rows: Iterable[dict[str, str]], max_q: float, fraction: float,
        minimum_blocks: int, database_path: str,
        ) -> set[tuple[str, str, str, str]]:
    """Reduce pair-wide votes by independent block using a disk table."""
    connection = sqlite3.connect(database_path)
    configure_disk_database(connection)
    connection.execute("""
        CREATE TABLE pairwide_votes (
            pair TEXT NOT NULL,
            chromosome TEXT NOT NULL,
            arm TEXT NOT NULL,
            block TEXT NOT NULL,
            state TEXT NOT NULL,
            votes INTEGER NOT NULL,
            PRIMARY KEY (pair, chromosome, arm, block, state)
        ) WITHOUT ROWID
    """)
    statement = """
        INSERT INTO pairwide_votes
            (pair, chromosome, arm, block, state, votes)
        VALUES (?, ?, ?, ?, ?, 1)
        ON CONFLICT (pair, chromosome, arm, block, state)
        DO UPDATE SET votes = votes + 1
    """
    pending: list[tuple[str, str, str, str, str]] = []
    connection.execute("BEGIN")
    for row in rows:
        pair = present_identifier(row.get("donor_pair"))
        if not pair:
            continue
        axes = preliminary_axes(row, threshold=max_q)
        state = clean(axes.get("conjunction_state"))
        if state not in EXACT_STATES:
            state = "NO_CALL"
        pending.append((
            pair, clean(row.get("chromosome")), clean(row.get("arm")),
            independent_block(row), state))
        if len(pending) >= SQLITE_BATCH_ROWS:
            connection.executemany(statement, pending)
            pending.clear()
    if pending:
        connection.executemany(statement, pending)
    connection.commit()

    flagged: set[tuple[str, str, str, str]] = set()
    current_block: tuple[str, str, str, str] | None = None
    block_votes: Counter = Counter()
    current_group: tuple[str, str, str] | None = None
    group_votes: Counter = Counter()

    def finish_group() -> None:
        if current_group is None:
            return
        total_blocks = sum(group_votes.values())
        state = counter_choice(
            group_votes, set(EXACT_STATES) | {"NO_CALL", "BALANCED"},
            "NO_CALL")
        if (total_blocks >= minimum_blocks and state in EXACT_STATES
                and group_votes[state] / total_blocks >= fraction):
            flagged.add((*current_group, state))

    def finish_block() -> None:
        nonlocal current_group, group_votes
        if current_block is None:
            return
        group = current_block[:3]
        if current_group is not None and group != current_group:
            finish_group()
            group_votes = Counter()
        current_group = group
        selected = counter_choice(
            block_votes, set(EXACT_STATES) | {"NO_CALL"}, "NO_CALL")
        group_votes[selected] += 1

    cursor = connection.execute("""
        SELECT pair, chromosome, arm, block, state, votes
        FROM pairwide_votes
        ORDER BY pair, chromosome, arm, block, state
    """)
    for pair, chromosome, arm, block, state, votes in cursor:
        key = (pair, chromosome, arm, block)
        if current_block is not None and key != current_block:
            finish_block()
            block_votes = Counter()
        current_block = key
        block_votes[state] += int(votes)
    finish_block()
    finish_group()
    connection.execute("DROP TABLE pairwide_votes")
    connection.commit()
    connection.close()
    return flagged


def nullable_minimum(previous: float | None, candidate: float | None
                     ) -> float | None:
    if candidate is None:
        return previous
    if previous is None:
        return candidate
    return min(previous, candidate)


class DiskBlockStore:
    """Bounded-memory block summaries for cohort-scale aggregation."""

    def __init__(self, database_path: str):
        self.connection = sqlite3.connect(database_path)
        configure_disk_database(self.connection)
        self.connection.row_factory = sqlite3.Row
        counter_columns = ",\n".join(
            f"{column} INTEGER NOT NULL"
            for column in (*BLOCK_CLASS_COLUMNS.values(),
                           *BLOCK_STATE_COLUMNS.values()))
        self.connection.execute(f"""
            CREATE TABLE block_summary (
                pair TEXT NOT NULL,
                chromosome TEXT NOT NULL,
                arm TEXT NOT NULL,
                block TEXT NOT NULL,
                uid TEXT NOT NULL,
                total_cells INTEGER NOT NULL,
                min_expression_q REAL,
                min_ase_q REAL,
                min_hybrid_q REAL,
                pairwide INTEGER NOT NULL,
                {counter_columns},
                PRIMARY KEY (pair, chromosome, arm, block)
            ) WITHOUT ROWID
        """)
        self.connection.execute("""
            CREATE TABLE block_libraries (
                pair TEXT NOT NULL, chromosome TEXT NOT NULL,
                arm TEXT NOT NULL, block TEXT NOT NULL, library TEXT NOT NULL,
                PRIMARY KEY (pair, chromosome, arm, block, library)
            ) WITHOUT ROWID
        """)
        self.connection.execute("""
            CREATE TABLE block_groups (
                pair TEXT NOT NULL, chromosome TEXT NOT NULL,
                arm TEXT NOT NULL, block TEXT NOT NULL,
                calibration_group TEXT NOT NULL,
                PRIMARY KEY (
                    pair, chromosome, arm, block, calibration_group)
            ) WITHOUT ROWID
        """)
        self.connection.execute("""
            CREATE TABLE uid_cell_counts (
                uid TEXT NOT NULL, pair TEXT NOT NULL,
                chromosome TEXT NOT NULL, n_cells INTEGER NOT NULL,
                PRIMARY KEY (uid, pair, chromosome)
            ) WITHOUT ROWID
        """)
        self.summary_fields = [
            "pair", "chromosome", "arm", "block", "uid", "total_cells",
            "min_expression_q", "min_ase_q", "min_hybrid_q", "pairwide",
            *BLOCK_CLASS_COLUMNS.values(), *BLOCK_STATE_COLUMNS.values(),
        ]
        additions = [
            "total_cells = block_summary.total_cells + excluded.total_cells",
            "min_expression_q = CASE "
            "WHEN excluded.min_expression_q IS NULL THEN block_summary.min_expression_q "
            "WHEN block_summary.min_expression_q IS NULL THEN excluded.min_expression_q "
            "ELSE MIN(block_summary.min_expression_q, excluded.min_expression_q) END",
            "min_ase_q = CASE "
            "WHEN excluded.min_ase_q IS NULL THEN block_summary.min_ase_q "
            "WHEN block_summary.min_ase_q IS NULL THEN excluded.min_ase_q "
            "ELSE MIN(block_summary.min_ase_q, excluded.min_ase_q) END",
            "min_hybrid_q = CASE "
            "WHEN excluded.min_hybrid_q IS NULL THEN block_summary.min_hybrid_q "
            "WHEN block_summary.min_hybrid_q IS NULL THEN excluded.min_hybrid_q "
            "ELSE MIN(block_summary.min_hybrid_q, excluded.min_hybrid_q) END",
            "pairwide = MAX(block_summary.pairwide, excluded.pairwide)",
            *(f"{column} = block_summary.{column} + excluded.{column}"
              for column in (*BLOCK_CLASS_COLUMNS.values(),
                             *BLOCK_STATE_COLUMNS.values())),
        ]
        placeholders = ",".join("?" for _field in self.summary_fields)
        self.summary_statement = (
            f"INSERT INTO block_summary ({','.join(self.summary_fields)}) "
            f"VALUES ({placeholders}) ON CONFLICT "
            "(pair, chromosome, arm, block) DO UPDATE SET "
            + ",".join(additions))
        self.summary_pending: list[tuple[object, ...]] = []
        self.library_pending: list[tuple[str, str, str, str, str]] = []
        self.group_pending: list[tuple[str, str, str, str, str]] = []
        self.uid_pending: list[tuple[str, str, str]] = []
        self.connection.execute("BEGIN")

    def add(self, row: Mapping[str, object], new_cell_chromosome: bool) -> None:
        pair = present_identifier(row.get("donor_pair"))
        if not pair:
            return
        chromosome = clean(row.get("chromosome"))
        arm = clean(row.get("arm"))
        block = independent_block(row)
        uid = present_identifier(row.get("uid"))
        evidence_class = clean(row.get("evidence_class"))
        if evidence_class not in BLOCK_CLASSES:
            evidence_class = "INSUFFICIENT_EVIDENCE"
        state = clean(row.get("resolved_state"))
        if state not in BLOCK_STATES:
            state = "NO_CALL"
        q_values = tuple(
            value if math.isfinite(value) else None
            for value in model_selected_qs(row))
        counters = [
            int(evidence_class == value) for value in BLOCK_CLASSES]
        counters.extend(int(state == value) for value in BLOCK_STATES)
        self.summary_pending.append((
            pair, chromosome, arm, block, uid, 1, *q_values,
            int("PAIRWIDE_EVENT_NOT_IDENTIFIABLE_FROM_RNA" in split_flags(
                row.get("confounding_flags"))), *counters))
        self.library_pending.append((
            pair, chromosome, arm, block,
            clean(row.get("library")) or "NA"))
        self.group_pending.append((
            pair, chromosome, arm, block,
            clean(row.get("calibration_group")) or "NA"))
        if uid and new_cell_chromosome:
            self.uid_pending.append((uid, pair, chromosome))
        if len(self.summary_pending) >= SQLITE_BATCH_ROWS:
            self.flush()

    def flush(self) -> None:
        if not self.summary_pending:
            return
        self.connection.executemany(
            self.summary_statement, self.summary_pending)
        self.connection.executemany("""
            INSERT OR IGNORE INTO block_libraries
                (pair, chromosome, arm, block, library)
            VALUES (?, ?, ?, ?, ?)
        """, self.library_pending)
        self.connection.executemany("""
            INSERT OR IGNORE INTO block_groups
                (pair, chromosome, arm, block, calibration_group)
            VALUES (?, ?, ?, ?, ?)
        """, self.group_pending)
        self.connection.executemany("""
            INSERT INTO uid_cell_counts (uid, pair, chromosome, n_cells)
            VALUES (?, ?, ?, 1)
            ON CONFLICT (uid, pair, chromosome)
            DO UPDATE SET n_cells = n_cells + 1
        """, self.uid_pending)
        self.summary_pending.clear()
        self.library_pending.clear()
        self.group_pending.clear()
        self.uid_pending.clear()

    def finish(self) -> None:
        self.flush()
        self.connection.commit()
        self.connection.execute("""
            CREATE INDEX block_uid_order
            ON block_summary (uid, pair, chromosome, arm, block)
        """)
        self.connection.commit()

    @staticmethod
    def _counter(row: Mapping[str, object],
                 columns: Mapping[str, str]) -> Counter:
        return Counter({
            value: int(row[column]) for value, column in columns.items()
            if int(row[column])
        })

    @staticmethod
    def _token_groups(cursor, key_length: int):
        current: tuple[str, ...] | None = None
        values: list[str] = []
        for row in cursor:
            key = tuple(str(row[index]) for index in range(key_length))
            if current is not None and key != current:
                yield current, values
                values = []
            current = key
            values.append(str(row[key_length]))
        if current is not None:
            yield current, values

    @staticmethod
    def _next_tokens(iterator, key: tuple[str, ...], label: str) -> list[str]:
        try:
            observed, values = next(iterator)
        except StopIteration as exc:
            raise ValueError(f"missing disk-backed {label} group for {key}") from exc
        if observed != key:
            raise ValueError(
                f"disk-backed {label} group mismatch: {observed} != {key}")
        return values

    def iter_pair_summaries(self) -> Iterable[dict[str, object]]:
        libraries = self._token_groups(self.connection.execute("""
            SELECT DISTINCT pair, chromosome, arm, library
            FROM block_libraries
            ORDER BY pair COLLATE NATURALKEY,
                     chromosome COLLATE NATURALKEY,
                     arm COLLATE NATURALKEY,
                     library COLLATE NATURALKEY
        """), 3)
        groups = self._token_groups(self.connection.execute("""
            SELECT DISTINCT pair, chromosome, arm, calibration_group
            FROM block_groups
            ORDER BY pair COLLATE NATURALKEY,
                     chromosome COLLATE NATURALKEY,
                     arm COLLATE NATURALKEY,
                     calibration_group COLLATE NATURALKEY
        """), 3)
        cursor = self.connection.execute("""
            SELECT * FROM block_summary
            ORDER BY pair COLLATE NATURALKEY,
                     chromosome COLLATE NATURALKEY,
                     arm COLLATE NATURALKEY,
                     block COLLATE NATURALKEY
        """)
        current: tuple[str, str, str] | None = None
        total_cells = distinct_blocks = 0
        cell_classes: Counter = Counter()
        block_classes: Counter = Counter()
        block_states: Counter = Counter()
        min_expression_q = min_ase_q = min_hybrid_q = None
        pairwide = False

        def build() -> dict[str, object]:
            if current is None:
                raise AssertionError("empty pair summary group")
            pair, chromosome, arm = current
            exact_states = Counter({
                state: count for state, count in block_states.items()
                if state in EXACT_STATES})
            best_state = counter_choice(
                block_states, set(EXACT_STATES) | {"NO_CALL", "BALANCED"},
                "NO_CALL")
            library_values = self._next_tokens(
                libraries, current, "pair library")
            group_values = self._next_tokens(groups, current, "pair calibration")
            return {
                "donor_pair": pair, "chromosome": chromosome, "arm": arm,
                "libraries": ";".join(library_values),
                "calibration_groups": ";".join(group_values),
                "total_cells": total_cells,
                "distinct_uid_blocks": distinct_blocks,
                "support_basis": "UID_BLOCK_FIRST",
                "expression_support_cells": (
                    cell_classes["CONCORDANT_BOTH"]
                    + cell_classes["EXPRESSION_ONLY"]),
                "ase_support_cells": (
                    cell_classes["CONCORDANT_BOTH"]
                    + cell_classes["ASE_ONLY"]),
                "concordant_cells": cell_classes["CONCORDANT_BOTH"],
                "expression_only_cells": cell_classes["EXPRESSION_ONLY"],
                "ase_only_cells": cell_classes["ASE_ONLY"],
                "discordant_cells": (
                    cell_classes["DISCORDANT"]
                    + cell_classes["EXPRESSION_OUTLIER"]),
                "balanced_cells": cell_classes["BALANCED"],
                "expression_support_uid_blocks": (
                    block_classes["CONCORDANT_BOTH"]
                    + block_classes["EXPRESSION_ONLY"]),
                "ase_support_uid_blocks": (
                    block_classes["CONCORDANT_BOTH"]
                    + block_classes["ASE_ONLY"]),
                "concordant_uid_blocks": block_classes["CONCORDANT_BOTH"],
                "expression_only_uid_blocks": block_classes["EXPRESSION_ONLY"],
                "ase_only_uid_blocks": block_classes["ASE_ONLY"],
                "discordant_uid_blocks": (
                    block_classes["DISCORDANT"]
                    + block_classes["EXPRESSION_OUTLIER"]),
                "balanced_uid_blocks": block_classes["BALANCED"],
                "best_resolved_state": best_state,
                "resolved_state_counts": ";".join(
                    f"{state}={count}" for state, count in sorted(
                        exact_states.items(),
                        key=lambda item: natural_key(item[0]))) or "NONE",
                "min_expression_q": format_number(min_expression_q),
                "min_ase_q": format_number(min_ase_q),
                "min_hybrid_q": format_number(min_hybrid_q),
                "pairwide_identifiability_flag": "TRUE" if pairwide else "FALSE",
                "summary_status": (
                    "RNA_CANDIDATE" if best_state in EXACT_STATES else "NO_CALL"),
                "schema_version": HYBRID_PAIR_SCHEMA,
            }

        for row in cursor:
            key = (str(row["pair"]), str(row["chromosome"]), str(row["arm"]))
            if current is not None and key != current:
                yield build()
                total_cells = distinct_blocks = 0
                cell_classes = Counter()
                block_classes = Counter()
                block_states = Counter()
                min_expression_q = min_ase_q = min_hybrid_q = None
                pairwide = False
            current = key
            classes = self._counter(row, BLOCK_CLASS_COLUMNS)
            states = self._counter(row, BLOCK_STATE_COLUMNS)
            total_cells += int(row["total_cells"])
            distinct_blocks += 1
            cell_classes.update(classes)
            block_classes[counter_choice(
                classes, set(BLOCK_CLASSES), "INSUFFICIENT_EVIDENCE")] += 1
            block_states[counter_choice(
                states, set(BLOCK_STATES), "NO_CALL")] += 1
            min_expression_q = nullable_minimum(
                min_expression_q, row["min_expression_q"])
            min_ase_q = nullable_minimum(min_ase_q, row["min_ase_q"])
            min_hybrid_q = nullable_minimum(min_hybrid_q, row["min_hybrid_q"])
            pairwide = pairwide or bool(row["pairwide"])
        if current is not None:
            yield build()

    def iter_uid_summaries(self) -> Iterable[dict[str, object]]:
        libraries = self._token_groups(self.connection.execute("""
            SELECT DISTINCT b.uid, b.pair, b.chromosome, l.library
            FROM block_summary AS b
            JOIN block_libraries AS l
              ON l.pair = b.pair AND l.chromosome = b.chromosome
             AND l.arm = b.arm AND l.block = b.block
            WHERE b.uid <> ''
            ORDER BY b.uid COLLATE NATURALKEY,
                     b.pair COLLATE NATURALKEY,
                     b.chromosome COLLATE NATURALKEY,
                     l.library COLLATE NATURALKEY
        """), 3)
        cell_counts = self.connection.execute("""
            SELECT uid, pair, chromosome, n_cells
            FROM uid_cell_counts
            ORDER BY uid COLLATE NATURALKEY,
                     pair COLLATE NATURALKEY,
                     chromosome COLLATE NATURALKEY
        """)
        count_iterator = iter(cell_counts)
        cursor = self.connection.execute("""
            SELECT * FROM block_summary
            WHERE uid <> ''
            ORDER BY uid COLLATE NATURALKEY,
                     pair COLLATE NATURALKEY,
                     chromosome COLLATE NATURALKEY,
                     arm COLLATE NATURALKEY,
                     block COLLATE NATURALKEY
        """)
        current: tuple[str, str, str] | None = None
        side_data: dict[str, dict[str, object]] = {}

        def empty_side() -> dict[str, object]:
            return {
                "arms": set(), "classes": Counter(), "states": Counter(),
                "min_expression_q": None, "min_ase_q": None,
                "min_hybrid_q": None,
            }

        def build() -> dict[str, object]:
            if current is None:
                raise AssertionError("empty UID summary group")
            uid, pair, chromosome = current
            try:
                count_row = next(count_iterator)
            except StopIteration as exc:
                raise ValueError(f"missing UID cell count for {current}") from exc
            count_key = tuple(str(count_row[index]) for index in range(3))
            if count_key != current:
                raise ValueError(
                    f"UID cell-count group mismatch: {count_key} != {current}")
            summary: dict[str, object] = {
                "uid": uid, "donor_pair": pair, "chromosome": chromosome,
                "libraries": ";".join(self._next_tokens(
                    libraries, current, "UID library")),
                "n_cells": int(count_row[3]), "distinct_uid_blocks": 1,
                "support_basis": "WITHIN_UID_CLONE_REDUCED_ONCE",
            }
            selected: dict[str, str] = {}
            for side in ("p", "q"):
                values = side_data.get(side, empty_side())
                classes = values["classes"]
                states = values["states"]
                evidence_class = counter_choice(
                    classes, set(BLOCK_CLASSES), "NO_DATA")
                resolved_state = counter_choice(
                    states, set(BLOCK_STATES), "NO_CALL")
                selected[side] = resolved_state
                arms = sorted(values["arms"], key=natural_key)
                summary.update({
                    f"{side}_arm": arms[0] if arms else "NA",
                    f"{side}_evidence_class": evidence_class,
                    f"{side}_resolved_state": resolved_state,
                    f"{side}_expression_support_cells": (
                        classes["CONCORDANT_BOTH"]
                        + classes["EXPRESSION_ONLY"]),
                    f"{side}_ase_support_cells": (
                        classes["CONCORDANT_BOTH"] + classes["ASE_ONLY"]),
                    f"{side}_concordant_cells": classes["CONCORDANT_BOTH"],
                    f"{side}_min_expression_q": format_number(
                        values["min_expression_q"]),
                    f"{side}_min_ase_q": format_number(values["min_ase_q"]),
                    f"{side}_min_hybrid_q": format_number(values["min_hybrid_q"]),
                })
            relationship = classify_pq_relationship(
                selected.get("p", "NO_CALL"), selected.get("q", "NO_CALL"))
            summary["pq_relationship"] = relationship
            summary["summary_status"] = (
                "RNA_CANDIDATE" if relationship != "NO_RESOLVED_EVENT"
                else "NO_CALL")
            summary["schema_version"] = HYBRID_UID_SCHEMA
            return summary

        for row in cursor:
            key = (str(row["uid"]), str(row["pair"]), str(row["chromosome"]))
            if current is not None and key != current:
                yield build()
                side_data = {}
            current = key
            arm = str(row["arm"])
            side = arm[-1:].lower()
            if side not in {"p", "q"}:
                continue
            values = side_data.setdefault(side, empty_side())
            values["arms"].add(arm)
            values["classes"].update(self._counter(
                row, BLOCK_CLASS_COLUMNS))
            values["states"].update(self._counter(
                row, BLOCK_STATE_COLUMNS))
            for field in (
                    "min_expression_q", "min_ase_q", "min_hybrid_q"):
                values[field] = nullable_minimum(values[field], row[field])
        if current is not None:
            yield build()

    def close(self) -> None:
        self.connection.close()


def iter_merged_diagnostics(
        tasks: Sequence[Mapping[str, str]], member: str,
        key_fields: Sequence[str]) -> Iterable[dict[str, str]]:
    def source(task_index: int, task: Mapping[str, str]):
        path = model_paths(task)[member]
        for row_index, row in enumerate(read_tsv(path)):
            key = tuple(natural_key(row.get(field, "")) for field in key_fields)
            yield key, task_index, row_index, row

    previous: tuple[object, ...] | None = None
    sources = [source(index, task) for index, task in enumerate(tasks)]
    for key, _task_index, _row_index, row in heapq.merge(
            *sources, key=lambda value: (value[0], value[1], value[2])):
        if previous is not None and key <= previous:
            raise ValueError(f"duplicate/unsorted merged {member} row")
        previous = key
        yield row


def main_impl(args) -> int:
    if not 0.0 < args.max_q <= 1.0:
        raise ValueError("--max-q must be in (0,1]")
    if not 0.0 < args.pairwide_fraction <= 1.0:
        raise ValueError("--pairwide-fraction must be in (0,1]")
    if args.min_pairwide_cells < 1:
        raise ValueError("--min-pairwide-cells must be positive")
    tasks = read_model_tasks(args.model_task_manifest)
    prefix = os.path.abspath(args.output_prefix)
    output_paths = {
        "calls": prefix + ".arm_calls.tsv.gz",
        "calibration": prefix + ".calibration.tsv.gz",
        "components": prefix + ".expression_components.tsv.gz",
        "uid": prefix + ".uid_chromosome_flags.tsv.gz",
        "pair": prefix + ".donor_pair_arm_summary.tsv.gz",
        "qc": prefix + ".qc.tsv",
        "contract": prefix + ".contract.json",
    }
    require_outputs_absent(output_paths.values())
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)

    input_records: list[dict[str, dict[str, object]]] = []
    score_header: list[str] | None = None
    component_header: list[str] | None = None
    calibration_header: list[str] | None = None
    offsets: dict[int, int] = {}
    component_count = calibration_count = total_rows = 0

    with tempfile.TemporaryDirectory(
            prefix=".tetra_arm_hybrid_aggregate.",
            dir=str(Path(prefix).parent)) as temporary:
        hypothesis_paths = {
            family: os.path.join(temporary, f"{family}.hypotheses.bin")
            for family in FAMILY_P_FIELDS}
        handles = {family: open(path, "wb")
                   for family, path in hypothesis_paths.items()}
        hypothesis_counts = {family: 0 for family in FAMILY_P_FIELDS}
        try:
            for task_index, task in enumerate(tasks):
                chromosome = task["chromosome"]
                paths = model_paths(task)
                contract = require_contract(paths["contract"], chromosome)
                contract_outputs = contract.get("outputs", {})
                expected_outputs = {
                    "scores": paths["scores"],
                    "expression_components": paths["components"],
                    "calibration": paths["calibration"], "qc": paths["qc"],
                    "contract": paths["contract"],
                }
                if (not isinstance(contract_outputs, dict) or any(
                        os.path.abspath(str(contract_outputs.get(name, "")))
                        != os.path.abspath(expected)
                        for name, expected in expected_outputs.items())):
                    raise ValueError(
                        "hybrid chromosome contract output paths disagree: "
                        f"{paths['contract']}")
                observed_score_header, validated_score_count = (
                    validate_hybrid_table(paths["scores"], HYBRID_SCORE_SCHEMA))
                observed_component_header, observed_component_count = (
                    validate_hybrid_table(
                        paths["components"], HYBRID_COMPONENT_SCHEMA))
                observed_calibration_header, observed_calibration_count = (
                    validate_hybrid_table(
                        paths["calibration"], HYBRID_CALIBRATION_SCHEMA))
                require_file(paths["qc"], "hybrid chromosome QC")
                with open(paths["qc"], "r", encoding="utf-8", newline="") as handle:
                    qc_reader = csv.DictReader(handle, delimiter="\t")
                    if tuple(qc_reader.fieldnames or ()) != ("metric", "value"):
                        raise ValueError(
                            f"invalid hybrid chromosome QC: {paths['qc']}")
                    qc_rows = list(qc_reader)
                qc = {row.get("metric", ""): row.get("value", "")
                      for row in qc_rows}
                if (len(qc) != len(qc_rows)
                        or qc.get("schema_version")
                           != "tetra_arm_hybrid_model_qc_v1"
                        or canonical_chromosome(qc.get("chromosome"))
                           != chromosome
                        or clean(qc.get("status")).upper() != "PASS"):
                    raise ValueError(
                        f"invalid hybrid chromosome QC: {paths['qc']}")
                if score_header is None:
                    score_header = observed_score_header
                    component_header = observed_component_header
                    calibration_header = observed_calibration_header
                elif (score_header != observed_score_header
                      or component_header != observed_component_header
                      or calibration_header != observed_calibration_header):
                    raise ValueError(
                        "chromosome output headers differ across model tasks")
                if list(contract.get("score_fields", [])) != observed_score_header:
                    raise ValueError(
                        f"score header disagrees with contract: {paths['contract']}")
                forbidden_q = [field for field in observed_score_header
                               if "_q_" in field or field.startswith(
                                   ("ase_q_", "expression_q_"))]
                if forbidden_q:
                    raise ValueError(
                        "chromosome score shards must not contain q-values: "
                        f"{forbidden_q}")
                if (int(contract.get("scores", -1)) != validated_score_count
                        or int(contract.get("component_rows", -1))
                           != observed_component_count
                        or int(contract.get("calibration_rows", -1))
                           != observed_calibration_count):
                    raise ValueError(
                        f"bundle row count disagrees with contract: "
                        f"{paths['contract']}")
                offsets[task_index] = total_rows
                local_count = 0
                previous_local_key: tuple[object, ...] | None = None
                for row in read_tsv(paths["scores"]):
                    if canonical_chromosome(row.get("chromosome")) != chromosome:
                        raise ValueError(
                            f"wrong chromosome in score table: {paths['scores']}")
                    local_key = (
                        natural_key(row.get("library", "")),
                        natural_key(row.get("barcode", "")),
                        natural_key(row.get("arm", "")),
                    )
                    if previous_local_key is not None and local_key <= previous_local_key:
                        raise ValueError(
                            f"score table is duplicate/unsorted: {paths['scores']}")
                    previous_local_key = local_key
                    append_hypotheses(
                        row, total_rows + local_count, handles,
                        hypothesis_counts)
                    local_count += 1
                if local_count != validated_score_count:
                    raise ValueError(
                        f"score table changed during aggregation: {paths['scores']}")
                total_rows += local_count
                component_count += observed_component_count
                calibration_count += observed_calibration_count
                input_records.append({
                    key: file_record(value) for key, value in paths.items()})
        finally:
            for handle in handles.values():
                handle.close()

        q_maps = {
            family: external_bh_map(
                hypothesis_paths[family],
                os.path.join(temporary, f"{family}.qmap.bin"),
                hypothesis_counts[family], total_rows,
                len(FAMILY_P_FIELDS[family]))
            for family in FAMILY_P_FIELDS
        }
        hybrid_harmonic = harmonic_number(hypothesis_counts["hybrid"])
        ase_count = hypothesis_counts["ase"]
        expression_count = hypothesis_counts["expression"]
        hybrid_count = hypothesis_counts["hybrid"]
        ase_scope = f"INVOCATION_WIDE_BH_{ase_count}_DIRECTIONAL_HYPOTHESES"
        expression_scope = (
            f"INVOCATION_WIDE_BH_{expression_count}_COPY_HYPOTHESES")
        hybrid_scope = (
            f"INVOCATION_WIDE_BH_{hybrid_count}_EXACT_STATE_HYPOTHESES")

        summary_database = os.path.join(temporary, "block_summaries.sqlite3")
        pairwide = pairwide_groups_from_stream_disk(
            iter_adjusted_scores(tasks, offsets, q_maps, hybrid_harmonic),
            args.max_q, args.pairwide_fraction, args.min_pairwide_cells,
            summary_database)
        block_store = DiskBlockStore(summary_database)
        counts: Counter = Counter()
        previous_cell_chromosome: tuple[str, str, str] | None = None

        def final_rows() -> Iterable[dict[str, object]]:
            nonlocal previous_cell_chromosome
            for row in iter_adjusted_scores(
                    tasks, offsets, q_maps, hybrid_harmonic):
                axes = preliminary_axes(row, args.max_q)
                flags = split_flags(row.get("confounding_flags"))
                pairwide_key = (
                    present_identifier(row.get("donor_pair")),
                    clean(row.get("chromosome")), clean(row.get("arm")),
                    clean(axes["conjunction_state"]),
                )
                if pairwide_key in pairwide:
                    flags.add("PAIRWIDE_EVENT_NOT_IDENTIFIABLE_FROM_RNA")
                flags.add("ORTHOGONAL_VALIDATION_NOT_ASSESSED")
                resolved = resolve_evidence(
                    clean(axes["expression_state"]),
                    bool(axes["expression_significant"]),
                    clean(axes["ase_direction"]),
                    bool(axes["ase_significant"]),
                    clean(axes["conjunction_state"]),
                    bool(axes["conjunction_significant"]),
                    clean(row.get("expression_fold_replication_status")),
                    clean(row.get("ase_status")) == "PASS" and
                    clean(row.get("expression_status")) == "PASS", flags)
                row.update(resolved)
                if (resolved["evidence_class"] == "INSUFFICIENT_EVIDENCE"
                        and clean(row.get("ase_status")) in {"NO_DATA", ""}
                        and clean(row.get("expression_status")) in {"NO_DATA", ""}):
                    row["call_status"] = "NO_DATA"
                row["orthogonal_validation_status"] = (
                    "ORTHOGONAL_VALIDATION_NOT_ASSESSED")
                row["bh_scope_ase"] = ase_scope
                row["bh_scope_expression"] = expression_scope
                row["bh_scope_hybrid"] = hybrid_scope
                row["score_schema_version"] = HYBRID_SCORE_SCHEMA
                row["schema_version"] = HYBRID_CALL_SCHEMA
                cell_chromosome = (
                    clean(row.get("library")), clean(row.get("barcode")),
                    clean(row.get("chromosome")))
                new_cell_chromosome = (
                    cell_chromosome != previous_cell_chromosome)
                if new_cell_chromosome:
                    previous_cell_chromosome = cell_chromosome
                block_store.add(row, new_cell_chromosome)
                counts[clean(row.get("evidence_class"))] += 1
                yield row

        final_fields = [field for field in (score_header or [])
                        if field != "schema_version"]
        final_fields.extend(FINAL_ADDITIONAL_FIELDS)
        call_count = write_tsv_atomic(
            output_paths["calls"], final_rows(), final_fields,
            deterministic_gzip=True)
        if call_count != total_rows:
            raise ValueError("streamed call row count changed during aggregation")
        block_store.finish()
        uid_count = write_tsv_atomic(
            output_paths["uid"], block_store.iter_uid_summaries(), UID_FIELDS,
            deterministic_gzip=True)
        pair_count = write_tsv_atomic(
            output_paths["pair"], block_store.iter_pair_summaries(), PAIR_FIELDS,
            deterministic_gzip=True)
        block_store.close()
        written_calibration = write_tsv_atomic(
            output_paths["calibration"], iter_merged_diagnostics(
                tasks, "calibration",
                ("chromosome", "branch", "calibration_key", "fold")),
            calibration_header or [], deterministic_gzip=True)
        written_components = write_tsv_atomic(
            output_paths["components"], iter_merged_diagnostics(
                tasks, "components",
                ("chromosome", "calibration_key", "fold", "component")),
            component_header or [], deterministic_gzip=True)
        if (written_calibration != calibration_count
                or written_components != component_count):
            raise ValueError("streamed diagnostic row count changed")
        for q_map in q_maps.values():
            if q_map is not None:
                q_map.flush()

    terminal = "NONE" if total_rows else "PASS_NO_TARGET_ROWS"
    qc_rows = [
        {"metric": "schema_version", "value": HYBRID_QC_SCHEMA},
        {"metric": "call_rows", "value": total_rows},
        {"metric": "uid_rows", "value": uid_count},
        {"metric": "donor_pair_rows", "value": pair_count},
        {"metric": "ase_bh_hypotheses", "value": ase_count},
        {"metric": "expression_bh_hypotheses", "value": expression_count},
        {"metric": "hybrid_bh_hypotheses", "value": hybrid_count},
        *({"metric": f"evidence_class_{key}", "value": value}
          for key, value in sorted(counts.items(), key=lambda item: natural_key(item[0]))),
        {"metric": "terminal_state", "value": terminal},
        {"metric": "status", "value": "PASS"},
    ]
    write_tsv_atomic(output_paths["qc"], qc_rows, ["metric", "value"])
    write_json_atomic(output_paths["contract"], {
        "schema_version": HYBRID_CONTRACT_SCHEMA,
        "release": PROGRAM_VERSION,
        "inputs": {"model_task_manifest": file_record(args.model_task_manifest),
                   "chromosome_bundles": input_records},
        "outputs": {key: value for key, value in output_paths.items()},
        "output_schemas": {
            "calls": HYBRID_CALL_SCHEMA,
            "calibration": HYBRID_CALIBRATION_SCHEMA,
            "expression_components": HYBRID_COMPONENT_SCHEMA,
            "uid_chromosome_flags": HYBRID_UID_SCHEMA,
            "donor_pair_arm_summary": HYBRID_PAIR_SCHEMA,
            "qc": HYBRID_QC_SCHEMA, "contract": HYBRID_CONTRACT_SCHEMA,
        },
        "call_fields": final_fields, "calls": total_rows,
        "calibration_rows": calibration_count,
        "expression_component_rows": component_count,
        "uid_rows": uid_count, "donor_pair_rows": pair_count,
        "bh_scope": {"ase": ase_scope, "expression": expression_scope,
                     "hybrid": hybrid_scope},
        "by_sensitivity": "HYBRID_EXACT_STATE_FAMILY_REPORTED_PER_HYPOTHESIS",
        "hybrid_statistic": "MAX_OF_BRANCH_EMPIRICAL_P_VALUES",
        "fdr_eligibility": "PASS_CALIBRATIONS_ONLY_COMPLETE_FAMILIES",
        "direction_selection": "MODEL_SELECTED_DIRECTION_ONLY",
        "summary_support_basis": "UID_BLOCK_FIRST",
        "aggregation_engine": (
            "DISK_BACKED_EXTERNAL_BH_AND_UID_BLOCK_STREAMING_OUTPUT_V1"),
        "posterior_claim": False,
        "orthogonal_validation_status": "ORTHOGONAL_VALIDATION_NOT_ASSESSED",
        "terminal_state": terminal, "status": "PASS",
    })
    return 0


def self_test() -> int:
    values = bh_adjust([0.01, 0.04, 0.03, math.nan])
    if values[:3] != [0.03, 0.04, 0.04] or not math.isnan(values[3]):
        raise AssertionError("BH self-test failed")
    if by_adjust([0.01])[0] != 0.01:
        raise AssertionError("BY self-test failed")
    print("PASS tetra_arm_hybrid_aggregate self-test")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Gather hybrid chromosome bundles and apply global BH.")
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {PROGRAM_VERSION}")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--model-task-manifest")
    parser.add_argument("--output-prefix")
    parser.add_argument("--max-q", type=float, default=0.05)
    parser.add_argument("--pairwide-fraction", type=float, default=0.80)
    parser.add_argument(
        "--min-pairwide-cells", type=int, default=2,
        help=("Minimum independent UID/fallback blocks for pairwide "
              "identifiability review; the legacy option name is retained"))
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.self_test:
            return self_test()
        missing = [name for name in ("model_task_manifest", "output_prefix")
                   if getattr(args, name) in {None, ""}]
        if missing:
            raise ValueError("missing required option(s): " + ", ".join(missing))
        return main_impl(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
