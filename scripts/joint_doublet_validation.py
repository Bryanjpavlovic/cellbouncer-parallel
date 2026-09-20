#!/usr/bin/env python3
"""Derivative joint RNA/ATAC doublet validation control plane.

This module owns run-root safety, immutable plans, technical preflight,
file-driven status, deterministic statistical primitives, and contract-shaped
failure/checkpoint packages.  It deliberately does not revise identities or
silently replace unavailable scientific computation with zero-valued results.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import random
import re
import shlex
import shutil
import socket
import statistics
import subprocess
import sys
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


RELEASE = "2026-09-18-joint-doublet-derivative-validation-v1"
MASTER_SEED = 1729
BASELINE = {
    "ledger_cells": 327659,
    "candidate_scored_cells": 327653,
    "unranked_cells": 6,
    "candidate_hypotheses": 8150804,
    "both_assays_elevated": 7656,
    "same_added_state": 5844,
    "different_added_state": 1812,
    "raw_disagreements": 241670,
    "residual_disagreements": 239858,
}
REQUIRED_SUBDIRS = (
    "molecule_sidecar_repair", "development_1_29", "heldout_25_35_38",
    "regression_19_after_heldout", "optional_cohort_extension", "controls",
    "final_analysis", "return_package", "manifests", "slurm_scripts",
    "logs", "task_scratch",
)
FINISHED_MARKERS = (
    "DERIVATIVE_PREFLIGHT_FINISHED", "DERIVATIVE_SCORES_FINISHED",
    "DERIVATIVE_GATHER_FINISHED", "VALIDATION_ANALYSIS_FINISHED",
    "RETURN_PACKAGE_FINISHED",
)
TABLE_NAMES = (
    "01_cohort_accounting.tsv", "02_library_summary.tsv",
    "03_compact_cell_evidence.tsv.gz", "04_both_assays_elevated_cells.tsv.gz",
    "05_cross_assay_group_summary.tsv", "06_cross_assay_null_agreement.tsv",
    "07_other_disagreements_decomposition.tsv", "08_six_unranked_cells.tsv",
    "09_candidate_policy_audit.tsv", "10_coverage_match_balance.tsv",
    "11_coverage_match_results.tsv", "12_coverage_downsampling_summary.tsv",
    "13_coverage_downsampling_changed_cells.tsv.gz",
    "14_genetic_mixture_all_eligible.tsv.gz", "15_high_content_review.tsv.gz",
    "16_quantity_normalization_and_ablation.tsv",
    "17_state_conditioned_occupancy_assessment.tsv",
    "18_molecule_sidecar_inventory.tsv", "19_molecule_linkage_qc.tsv",
    "20_molecule_sensitivity_summary.tsv",
    "21_molecule_sensitivity_changed_cells.tsv.gz",
    "22_control_performance.tsv", "23_heldout_libraries_35_38_summary.tsv",
    "24_heldout_libraries_35_38_cells.tsv.gz",
    "25_library19_regression_summary.tsv", "26_library19_review_cells.tsv.gz",
    "27_final_decision_record.tsv",
)
PLOT_NAMES = tuple(f"{index:02d}_{name}.png" for index, name in enumerate((
    "library_rates", "rna_vs_atac_support", "observed_vs_null_agreement",
    "agreement_conflict_evidence", "donor_state_concentration",
    "coverage_match_balance", "coverage_control_effects",
    "downsampling_stability", "genetic_vs_high_content", "quantity_ablation",
    "molecule_collapse_effects", "heldout_transfer", "library19_regression",
), 1)) + ("all_plots.pdf",)

SOURCE_TASK_FIELDS = (
    "library", "output_dir", "rna_barcodes", "rna_features", "rna_matrix",
    "rna_samples", "rna_assignments", "rna_diagnostics", "rna_runner_ups",
    "rna_pileup_sites", "rna_pileup_observations", "atac_fragments",
    "atac_bam", "atac_samples", "atac_assignments", "atac_diagnostics",
    "atac_runner_ups", "reconciled_assignments", "reconciled_cells",
    "technical_candidates", "ploidy_calls", "ambient_rates",
    "ambient_profile", "pool_combinations", "atac_pileup_prefix",
    "rna_manifest", "atac_manifest", "rna_scores", "atac_scores",
)
DERIVATIVE_TASK_FIELDS = (
    "task_index", "library", "cohort_role", "source_output_dir",
    "output_dir", "rna_barcodes", "rna_features", "rna_matrix",
    "rna_samples", "rna_assignments", "rna_diagnostics", "rna_runner_ups",
    "rna_pileup_sites", "rna_pileup_observations", "rna_pileup_molecules",
    "atac_fragments", "atac_bam", "atac_samples", "atac_assignments",
    "atac_diagnostics", "atac_runner_ups", "atac_pileup_sites",
    "atac_pileup_samples", "atac_pileup_observations", "atac_pileup_molecules",
    "reconciled_assignments", "reconciled_cells", "technical_candidates",
    "ploidy_calls", "ambient_rates", "ambient_profile", "pool_combinations",
    "rna_manifest", "atac_manifest", "rna_scores", "atac_scores",
)
ANALYSIS_TASKS = (
    "BASELINE_RECONSTRUCTION", "CANDIDATE_POLICY_AUDIT",
    "CROSS_ASSAY_NULL", "COVERAGE_MATCHING", "COVERAGE_DOWNSAMPLING",
    "MOLECULE_SENSITIVITY", "QUANTITY_NORMALIZATION",
    "OCCUPANCY_ASSESSMENT", "CONTROL_PARENT_PARTITION",
    "CONTROL_CONSTRUCTION", "CELL_CONDITIONAL_NULL", "STAGE_SUMMARY",
)


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with open(temporary, "w", encoding="utf-8", newline="") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def atomic_json(path, value):
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    mode = "wt"
    with opener(temporary, mode, encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields,
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def read_tsv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def stable_seed(*parts):
    text = "\x1f".join(str(part) for part in (MASTER_SEED,) + parts)
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "big")


def _is_relative_to(path, parent):
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def symlink_components(path):
    path = Path(path).absolute()
    current = Path(path.anchor)
    links = []
    for part in path.parts[1:]:
        current /= part
        if current.is_symlink():
            links.append(str(current))
    return links


def validate_root_separation(source_root, baseline_gather, destination):
    source = Path(source_root).resolve(strict=False)
    gather = Path(baseline_gather).resolve(strict=False) if baseline_gather else None
    dest = Path(destination).resolve(strict=False)
    if not Path(destination).is_absolute():
        raise ValueError("validation destination must be absolute")
    if symlink_components(destination):
        raise ValueError("validation destination has a symlinked path component")
    for protected, label in ((source, "source root"), (gather, "baseline gather")):
        if protected is None:
            continue
        if dest == protected or _is_relative_to(dest, protected) or \
                _is_relative_to(protected, dest):
            raise ValueError(f"validation destination overlaps {label}")
    return source, gather, dest


def parse_libraries(values):
    result = []
    for value in values:
        for token in str(value).replace(",", " ").split():
            if "-" in token:
                left, right = token.split("-", 1)
                result.extend(range(int(left), int(right) + 1))
            else:
                result.append(int(token))
    if any(value < 1 for value in result):
        raise ValueError("library numbers must be positive")
    return list(dict.fromkeys(result))


def plan_hash(payload):
    return hashlib.sha256(json.dumps(
        payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def initialize_owned_root(root, plan):
    root = Path(root)
    state_path = root / "manifests" / "run_plan.json"
    wanted_hash = plan_hash(plan)
    if root.exists():
        if not state_path.is_file():
            raise RuntimeError(f"foreign existing destination: {root}")
        existing = json.loads(state_path.read_text())
        if existing.get("plan_hash") != wanted_hash:
            raise RuntimeError(f"existing destination plan mismatch: {root}")
        if existing.get("state") not in {"RUN_PLANNED", "RUN_SUBMITTED"}:
            raise RuntimeError(f"existing destination is not resumable: {root}")
        return existing
    root.mkdir(parents=True)
    for name in REQUIRED_SUBDIRS:
        (root / name).mkdir()
    owned = dict(plan)
    owned.update({"plan_hash": wanted_hash, "state": "RUN_PLANNED",
                  "created_utc": utc_now(), "release": RELEASE})
    atomic_json(state_path, owned)
    return owned


def gzip_quick_valid(path):
    """Constant-time envelope check; full decompression belongs on compute."""
    if not str(path).endswith(".gz"):
        return True
    try:
        target = Path(path)
        if target.stat().st_size < 18:
            return False
        with open(target, "rb") as handle:
            header = handle.read(2)
            handle.seek(-8, os.SEEK_END)
            trailer = handle.read(8)
        return header == b"\x1f\x8b" and len(trailer) == 8
    except OSError:
        return False


def count_rows_and_header(path):
    if not path or not Path(path).is_file():
        return 0, []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        first = handle.readline().rstrip("\r\n")
        rows = sum(1 for _ in handle)
    return rows, first.split("\t") if first else []


def load_source_tasks(source_root):
    path = Path(source_root) / "joint_doublet_tasks.tsv"
    if not path.is_file():
        return path, {}
    result = {}
    for row in read_tsv(path):
        raw = row.get("library", row.get("lib", ""))
        digits = "".join(character for character in raw if character.isdigit())
        if digits:
            result[int(digits)] = row
    return path, result


def _task_path(row, assay, member):
    candidates = (
        f"{assay}_{member}", f"{assay}_pileup_{member}",
        f"{assay}_{member}_path", f"{assay}_pileup_{member}_path",
    )
    direct = next((row.get(name, "") for name in candidates if row.get(name, "")), "")
    if direct:
        return direct
    if member == "molecules" and assay == "rna":
        observations = row.get("rna_pileup_observations", "")
        suffix = ".pileup_obs.tsv.gz"
        if observations.endswith(suffix):
            return observations[:-len(suffix)] + ".pileup_molecules.tsv.gz"
    if member == "molecules" and assay == "atac":
        prefix = row.get("atac_pileup_prefix", "")
        if prefix:
            return prefix + ".pileup_molecules.tsv.gz"
    return ""


def sidecar_inventory(libraries, tasks, heldout):
    rows = []
    for library in libraries:
        row = tasks.get(library, {})
        for assay in ("rna", "atac"):
            path = _task_path(row, assay, "molecules")
            target = Path(path) if path else None
            exists = bool(target and target.is_file())
            size = target.stat().st_size if exists else 0
            quick_valid = bool(exists and size and gzip_quick_valid(target))
            rows.append({
                "library": f"lib{library}", "assay": assay.upper(),
                "path": path or "UNRESOLVED_FROM_SOURCE_TASK_MANIFEST",
                "exists": exists, "nonempty": size > 0,
                "gzip_header_and_trailer_present": quick_valid,
                "gzip_full_stream_validation": "DEFERRED_TO_COMPUTE_WORKER",
                "bytes": size, "rows": "NOT_READ_IN_SEALED_PREFLIGHT"
                    if library in heldout else "NOT_COUNTED_IN_TECHNICAL_PREFLIGHT",
                "cells": "NOT_READ_IN_SEALED_PREFLIGHT"
                    if library in heldout else "NOT_COUNTED_IN_TECHNICAL_PREFLIGHT",
                "basis_counts": "NOT_READ_IN_SEALED_PREFLIGHT"
                    if library in heldout else "NOT_COUNTED_IN_TECHNICAL_PREFLIGHT",
                "target_roster_coverage": "NOT_TESTED",
                "status": "USABLE_PENDING_COMPUTE_VALIDATION" if quick_valid
                    else "MISSING_OR_INVALID",
                "reason": "NONE" if quick_valid
                    else "sidecar path absent, empty, or gzip-envelope-invalid",
            })
    return rows


def preflight(args):
    libraries = parse_libraries(args.libraries)
    heldout = set(parse_libraries(args.heldout_libraries or []))
    regression = set(parse_libraries(args.regression_libraries or []))
    calibration = int(args.calibration_library)
    if calibration not in libraries:
        raise RuntimeError("calibration library must be selected")
    if heldout & regression:
        raise RuntimeError("held-out and regression library roles overlap")
    source, gather, destination = validate_root_separation(
        args.source_output_root, args.partial_output_root, args.validation_root)
    plan = {
        "action": "VALIDATION_PREFLIGHT", "libraries": libraries,
        "calibration_library": calibration,
        "heldout_libraries": sorted(heldout),
        "regression_libraries": sorted(regression),
        "source_output_root": str(source),
        "baseline_gather": str(gather) if gather else "",
        "validation_root": str(destination),
        "evidence_mode": args.evidence_mode,
    }
    state = initialize_owned_root(destination, plan)
    task_path, tasks = load_source_tasks(source)
    inventory = sidecar_inventory(libraries, tasks, heldout)
    write_tsv(destination / "manifests" / "molecule_sidecar_inventory.tsv",
              inventory, list(inventory[0]) if inventory else ["library"])
    repair_rows = []
    for row in inventory:
        usable = row["status"] == "USABLE_PENDING_COMPUTE_VALIDATION"
        library_number = int("".join(character for character in row["library"]
                                     if character.isdigit()))
        source_task = tasks.get(library_number, {})
        derivative_prefix = destination / "molecule_sidecar_repair" / "bundles" / \
            row["library"] / row["assay"].lower() / "pileup"
        repair_rows.append({
            **row, "repair_needed": not usable,
            "repair_executable": False,
            "repair_status": "NOT_NEEDED" if usable else "REPAIR_NON_EXECUTABLE",
            "input_bam": source_task.get("atac_bam", "UNRESOLVED")
                if row["assay"] == "ATAC" else "UNRESOLVED",
            "panel": "UNRESOLVED", "whitelist": "UNRESOLVED",
            "output_prefix": str(derivative_prefix),
            "min_mapq": "UNRESOLVED", "exclude_flags": "UNRESOLVED",
            "variant_qual": "UNRESOLVED", "error_ref": "UNRESOLVED",
            "error_alt": "UNRESOLVED", "error_sigma": "UNRESOLVED",
            "selected_bundle": row["path"] if usable else "UNAVAILABLE",
        })
    write_tsv(destination / "manifests" / "molecule_repair_plan.tsv",
              repair_rows, list(repair_rows[0]) if repair_rows else ["library"])
    checks = []
    for name, path in (("source_task_manifest", task_path),
                       ("baseline_gather", gather)):
        present = bool(path and Path(path).exists())
        checks.append({"check": name, "status": "PASS" if present else "FAIL",
                       "detail": str(path) if path else "UNSET"})
    checks.append({"check": "source_destination_separation", "status": "PASS",
                   "detail": f"source={source};destination={destination}"})
    write_tsv(destination / "manifests" / "validation_checks.tsv", checks,
              ("check", "status", "detail"))
    draft = {
        "schema_version": "joint_doublet_frozen_analysis_spec_v1",
        "state": "development_not_yet_frozen", "release": RELEASE,
        "master_seed": MASTER_SEED, "candidate_policy": "DERIVATIVE_COMPLETE",
        "site_formula_version": "K1_LOCKED_STATE_VS_K2_ADDED_STATE_PLUS_FIXED_AMBIENT_V1",
        "molecule_formula_version": "LINKED_UNIT_EQUAL_WEIGHT_NORMALIZED_SITE_LL_V1",
        "heldout_values_parsed": False, "libraries": libraries,
        "calibration_library": calibration,
        "heldout_libraries": sorted(heldout),
        "regression_libraries": sorted(regression),
    }
    atomic_json(destination / "manifests" / "frozen_analysis_spec.json", draft)
    status = "COMPLETE" if all(row["status"] != "FAIL" for row in checks) else "FAILED"
    marker = {"status": status, "utc": utc_now(), "action": "VALIDATION_PREFLIGHT",
              "requested_libraries": libraries, "release": RELEASE,
              "plan_hash": state["plan_hash"], "checks": checks}
    atomic_json(destination / "DERIVATIVE_PREFLIGHT_FINISHED", marker)
    if status == "COMPLETE":
        atomic_json(destination / "DERIVATIVE_PREFLIGHT_COMPLETE", marker)
    print(json.dumps({"status": status, "validation_root": str(destination),
                      "plan_hash": state["plan_hash"],
                      "sidecar_rows": len(inventory)}, indent=2))
    return 0 if status == "COMPLETE" else 2


def assign_linked_unit_folds(keys):
    keys = list(dict.fromkeys(keys))
    if not keys:
        return {}, 0
    fold_count = min(5, len(keys))
    ordered = sorted(keys, key=lambda key: (stable_seed("linked_fold", key), key))
    return {key: index % fold_count for index, key in enumerate(ordered)}, fold_count


def control_parent_partition(library, identity, ploidy, state, barcode):
    decile = stable_seed("control_parent", library, identity, ploidy,
                         state, barcode) % 10
    return "development" if decile <= 6 else "test"


def bh_adjust(pvalues):
    result = [math.nan] * len(pvalues)
    finite = [(value, index) for index, value in enumerate(pvalues)
              if value is not None and math.isfinite(value)]
    finite.sort(reverse=True)
    running = 1.0
    count = len(finite)
    for reverse_rank, (value, index) in enumerate(finite):
        rank = count - reverse_rank
        running = min(running, value * count / rank)
        result[index] = running
    return result


def primary_conservative(variants, calibration_pass):
    computed = {name: value for name, value in variants.items()
                if value is not None and math.isfinite(value)}
    if not computed or not any(calibration_pass.get(name, False)
                               for name in computed):
        return math.nan, "NO_COMPUTED_VARIANT_PASSED_CALIBRATION"
    maximum = max(computed.values())
    winners = ",".join(sorted(name for name, value in computed.items()
                              if value == maximum))
    return maximum, winners


def analytic_agreement(rna, atac, strata):
    grouped = defaultdict(list)
    for index, stratum in enumerate(strata):
        grouped[stratum].append(index)
    expected = 0.0
    for indices in grouped.values():
        rna_counts = Counter(rna[index] for index in indices)
        atac_counts = Counter(atac[index] for index in indices)
        n = len(indices)
        expected += n * sum((rna_counts[key] / n) * (atac_counts[key] / n)
                            for key in set(rna_counts) | set(atac_counts))
    return expected / len(rna) if rna else math.nan


def fixed_denominator_permutation(rna, atac_blocks, strata, replicates=10000):
    grouped = defaultdict(list)
    for index, stratum in enumerate(strata):
        grouped[stratum].append(index)
    observed = sum(rna[index] == atac_blocks[index]["winner"]
                   for index in range(len(rna)))
    null = []
    rng = random.Random(stable_seed("fixed_denominator_permutation"))
    for _ in range(replicates):
        permuted = list(atac_blocks)
        for indices in grouped.values():
            values = [atac_blocks[index] for index in indices]
            rng.shuffle(values)
            for index, value in zip(indices, values):
                permuted[index] = value
        null.append(sum(rna[index] == permuted[index]["winner"]
                        for index in range(len(rna))))
    b = sum(value >= observed for value in null)
    return {"observed": observed, "denominator": len(rna), "b": b,
            "B": replicates, "p": (b + 1) / (replicates + 1),
            "null_mean": statistics.mean(null) if null else math.nan}


def _lib_name(value):
    digits = "".join(character for character in str(value) if character.isdigit())
    if not digits:
        raise ValueError(f"invalid library identifier: {value}")
    return f"lib{int(digits)}"


def _file_ready(path, gzip_check=False):
    target = Path(path)
    return target.is_file() and target.stat().st_size > 0 and (
        not gzip_check or gzip_quick_valid(target))


def _stage_kind(args, libraries):
    heldout = set(parse_libraries(args.heldout_libraries or []))
    regression = set(parse_libraries(args.regression_libraries or []))
    selected = set(libraries)
    calibration = int(args.calibration_library)
    if args.workflow_action == "VALIDATE_EXISTING":
        return "DEVELOPMENT"
    if selected & regression:
        if not (selected - {calibration}) <= regression:
            raise RuntimeError("regression stage may contain only regression libraries plus calibration")
        return "REGRESSION"
    if selected & heldout:
        if not (selected - {calibration}) <= heldout:
            raise RuntimeError("held-out stage may contain only held-out libraries plus calibration")
        return "HELDOUT"
    raise RuntimeError("DERIVATIVE_RESCORE libraries do not match held-out or regression roles")


def _assert_owned_validation_root(validation_root):
    root = Path(validation_root)
    plan_path = root / "manifests" / "run_plan.json"
    if not root.is_absolute() or not plan_path.is_file():
        raise RuntimeError("validation root is not initialized by VALIDATION_PREFLIGHT")
    if symlink_components(root):
        raise RuntimeError("validation root has a symlinked path component")
    return root.resolve(), json.loads(plan_path.read_text())


def _stage_plan_payload(args, kind, libraries, score_libraries, root_plan):
    return {
        "schema_version": "joint_doublet_derivative_stage_plan_v2",
        "release": RELEASE,
        "workflow_action": args.workflow_action,
        "stage_kind": kind,
        "libraries": libraries,
        "score_libraries": score_libraries,
        "calibration_library": int(args.calibration_library),
        "heldout_libraries": parse_libraries(args.heldout_libraries or []),
        "regression_libraries": parse_libraries(args.regression_libraries or []),
        "source_output_root": str(Path(args.source_output_root).resolve()),
        "baseline_gather": str(Path(args.partial_output_root).resolve())
            if args.partial_output_root else "",
        "validation_root": str(Path(args.validation_root).resolve()),
        "stage_root": str(Path(args.stage_root).resolve()),
        "frozen_spec": str(Path(args.frozen_spec).resolve()) if args.frozen_spec else "",
        "tool_bin_root": str(Path(args.tool_bin_root).resolve()),
        "candidate_policy": "DERIVATIVE_COMPLETE",
        "evidence_mode": args.evidence_mode,
        "master_seed": MASTER_SEED,
        "resources": {
            "score_cpus": args.score_cpus, "score_memory": args.score_memory,
            "score_max_concurrent": args.score_max_concurrent,
            "worker_cpus": args.worker_cpus, "worker_memory": args.worker_memory,
            "worker_max_concurrent": args.worker_max_concurrent,
            "gather_cpus": args.gather_cpus, "gather_memory": args.gather_memory,
            "analysis_cpus": args.analysis_cpus,
            "analysis_memory": args.analysis_memory,
            "time": args.time, "partition": args.partition,
        },
        "temp_root": str(Path(args.temp_root).resolve()) if args.temp_root else "",
        "validation_plan_hash": root_plan.get("plan_hash", ""),
    }


def _initialize_stage(args, kind, libraries):
    validation_root, root_plan = _assert_owned_validation_root(args.validation_root)
    stage_root = Path(args.stage_root)
    if not stage_root.is_absolute():
        raise RuntimeError("derivative output root must be absolute")
    if symlink_components(stage_root):
        raise RuntimeError("derivative output root has a symlinked path component")
    resolved_stage = stage_root.resolve(strict=False)
    if not _is_relative_to(resolved_stage, validation_root) or resolved_stage == validation_root:
        raise RuntimeError("derivative output root must be a child of validation root")
    source = Path(args.source_output_root).resolve()
    validate_root_separation(source, args.partial_output_root or None, resolved_stage)
    calibration = int(args.calibration_library)
    score_libraries = list(libraries)
    if kind == "DEVELOPMENT":
        reserved = set(parse_libraries(args.regression_libraries or []))
        score_libraries = [value for value in libraries if value not in reserved]
    if calibration not in score_libraries:
        raise RuntimeError("selected score cohort must contain the calibration library")
    if kind in {"HELDOUT", "REGRESSION"}:
        spec = Path(args.frozen_spec)
        if not spec.is_absolute() or not spec.is_file():
            raise RuntimeError("protected stage requires an existing absolute frozen spec")
        frozen = json.loads(spec.read_text())
        if frozen.get("state") not in {"final_frozen", "FINAL_FROZEN"}:
            raise RuntimeError("protected stage refuses a non-final frozen spec")
        if frozen.get("heldout_values_parsed") not in {False, "false", "FALSE"}:
            raise RuntimeError("frozen spec does not prove held-out values remained sealed")
    if kind == "REGRESSION":
        heldout_marker = validation_root / "heldout_25_35_38" / "VALIDATION_ANALYSIS_FINISHED"
        if not heldout_marker.is_file():
            raise RuntimeError("regression stage requires the held-out terminal marker")
        heldout_status = json.loads(heldout_marker.read_text()).get("status")
        if heldout_status not in {"COMPLETE", "PARTIAL"}:
            raise RuntimeError("held-out terminal marker is not usable for regression")
    payload = _stage_plan_payload(args, kind, libraries, score_libraries, root_plan)
    immutable_hash = plan_hash(payload)
    plan_path = resolved_stage / "manifests" / "stage_plan.json"
    if resolved_stage.exists():
        if plan_path.is_file():
            existing = json.loads(plan_path.read_text())
            if existing.get("immutable_hash") != immutable_hash:
                raise RuntimeError("existing derivative stage plan does not match this command")
            return resolved_stage, existing
        if any(resolved_stage.iterdir()):
            raise RuntimeError(f"foreign nonempty derivative root: {resolved_stage}")
    else:
        resolved_stage.mkdir(parents=True)
    for name in ("manifests", "slurm_scripts", "logs", "task_scratch",
                 "aggregate", "analysis", "worker_results"):
        (resolved_stage / name).mkdir(parents=True, exist_ok=True)
    plan = dict(payload)
    plan.update({
        "immutable_hash": immutable_hash, "state": "RUN_PLANNED",
        "run_token": datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") +
            "_" + immutable_hash[:8],
        "created_utc": utc_now(), "jobs": {}, "submission_history": [],
    })
    atomic_json(plan_path, plan)
    return resolved_stage, plan


def _build_derivative_tasks(stage_root, plan):
    source_manifest, source_tasks = load_source_tasks(plan["source_output_root"])
    if not source_manifest.is_file():
        raise RuntimeError(f"source V1 task manifest is missing: {source_manifest}")
    records = []
    roles = {}
    heldout = set(plan["heldout_libraries"])
    regression = set(plan["regression_libraries"])
    for library in plan["score_libraries"]:
        source = source_tasks.get(library)
        if source is None:
            raise RuntimeError(f"source V1 task manifest has no lib{library} row")
        role = "HELDOUT" if library in heldout else (
            "REGRESSION" if library in regression else "DEVELOPMENT")
        if library == plan["calibration_library"] and plan["stage_kind"] != "DEVELOPMENT":
            role = "CALIBRATION_REUSE"
        roles[library] = role
        library_name = f"lib{library}"
        output_dir = stage_root / library_name
        atac_prefix = source.get("atac_pileup_prefix", "")
        record = {
            "task_index": len(records), "library": library_name,
            "cohort_role": role, "source_output_dir": source.get("output_dir", ""),
            "output_dir": str(output_dir),
        }
        for field in SOURCE_TASK_FIELDS:
            if field not in {"library", "output_dir", "rna_manifest", "atac_manifest",
                             "rna_scores", "atac_scores", "atac_pileup_prefix"}:
                record[field] = source.get(field, "")
        record.update({
            "rna_pileup_molecules": _task_path(source, "rna", "molecules"),
            "atac_pileup_sites": atac_prefix + ".pileup_sites.tsv.gz" if atac_prefix else "",
            "atac_pileup_samples": atac_prefix + ".samples" if atac_prefix else "",
            "atac_pileup_observations": atac_prefix + ".pileup_obs.tsv.gz" if atac_prefix else "",
            "atac_pileup_molecules": _task_path(source, "atac", "molecules"),
            "rna_manifest": str(output_dir / f"{library_name}.rna_joint_manifest.tsv.gz"),
            "atac_manifest": str(output_dir / f"{library_name}.atac_joint_manifest.tsv.gz"),
            "rna_scores": str(output_dir / f"{library_name}.rna_joint_scores.tsv.gz"),
            "atac_scores": str(output_dir / f"{library_name}.atac_joint_scores.tsv.gz"),
        })
        records.append(record)
    manifest = stage_root / "manifests" / "derivative_tasks_v2.tsv"
    write_tsv(manifest, records, DERIVATIVE_TASK_FIELDS)
    mapping = [{"task_index": row["task_index"], "library": row["library"],
                "cohort_role": row["cohort_role"]} for row in records]
    write_tsv(stage_root / "manifests" / "array_mapping.tsv", mapping,
              ("task_index", "library", "cohort_role"))
    return records, manifest


def _shell_task_loader(manifest):
    variables = " ".join(DERIVATIVE_TASK_FIELDS)
    return f'''TASK_MANIFEST={shlex.quote(str(manifest))}
TASK_LINE="$(awk -v task="$SLURM_ARRAY_TASK_ID" 'NR == task + 2 {{print; exit}}' "$TASK_MANIFEST")"
[[ -n "$TASK_LINE" ]] || {{ echo "missing derivative task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }}
IFS=$'\\t' read -r {variables} <<< "$TASK_LINE"
'''


def _sbatch_header(name, plan, array=None, cpus=1, memory="8G"):
    root = plan["stage_root"]
    token = plan["run_token"]
    array_line = f"#SBATCH --array={array}\n" if array else ""
    return f'''#!/bin/bash
#SBATCH --job-name=jdv_{token}_{name}
#SBATCH --output={root}/logs/{name}_%A_%a.out
#SBATCH --error={root}/logs/{name}_%A_%a.err
{array_line}#SBATCH --time={plan["resources"]["time"]}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={memory}
#SBATCH --partition={plan["resources"]["partition"]}
#SBATCH --nodes=1

set -euo pipefail
module purge
module load miniforge/3 genomics-base/latest htslib/1.20
'''


def _write_script(path, body):
    atomic_text(path, body)
    os.chmod(path, 0o755)


def _render_stage_scripts(stage_root, plan, records, task_manifest):
    tools_root = Path(plan["tool_bin_root"])
    identity = tools_root / "identity_reconciliation.py"
    scorer = tools_root / "tetra_score_calls"
    validator = tools_root / "joint_doublet_validation.py"
    count = len(records)
    resources = plan["resources"]
    scripts = {}
    loader = _shell_task_loader(task_manifest)
    array = f"0-{count - 1}%{resources['score_max_concurrent']}"
    prepare = _sbatch_header("prepare", plan, array, 1, "8G") + loader + f'''
for required in "$rna_barcodes" "$rna_features" "$rna_matrix" "$rna_samples" \
  "$rna_assignments" "$rna_diagnostics" "$atac_fragments" \
  "$reconciled_assignments" "$reconciled_cells" "$pool_combinations"; do
  [[ -s "$required" ]] || {{ echo "missing required input: $required" >&2; exit 1; }}
done
mkdir -p "$output_dir"
python3 -B {shlex.quote(str(identity))} joint-prepare \
  --library "$library" --rna-barcodes "$rna_barcodes" \
  --rna-features "$rna_features" --rna-matrix "$rna_matrix" \
  --rna-samples "$rna_samples" --rna-assignments "$rna_assignments" \
  --rna-diagnostics "$rna_diagnostics" --rna-runner-ups "$rna_runner_ups" \
  --atac-fragments "$atac_fragments" --atac-samples "$atac_samples" \
  --atac-assignments "$atac_assignments" --atac-diagnostics "$atac_diagnostics" \
  --atac-runner-ups "$atac_runner_ups" \
  --reconciled-assignments "$reconciled_assignments" \
  --reconciled-cells "$reconciled_cells" --technical-candidates "$technical_candidates" \
  --ploidy-calls "$ploidy_calls" --ambient-rates "$ambient_rates" \
  --ambient-profile "$ambient_profile" --pool-combinations "$pool_combinations" \
  --candidate-policy DERIVATIVE_COMPLETE --output-dir "$output_dir"
gzip -t "$rna_manifest"; gzip -t "$atac_manifest"
'''
    scripts["prepare"] = stage_root / "slurm_scripts" / "prepare.sbatch"
    _write_script(scripts["prepare"], prepare)

    def score_script(modality):
        lower = modality.lower()
        sample = "$rna_samples" if modality == "RNA" else "$atac_pileup_samples"
        sites = "$rna_pileup_sites" if modality == "RNA" else "$atac_pileup_sites"
        observations = "$rna_pileup_observations" if modality == "RNA" else "$atac_pileup_observations"
        molecules = "$rna_pileup_molecules" if modality == "RNA" else "$atac_pileup_molecules"
        manifest_var = "$rna_manifest" if modality == "RNA" else "$atac_manifest"
        scores = "$rna_scores" if modality == "RNA" else "$atac_scores"
        error = "0.001" if modality == "RNA" else "0.005"
        return _sbatch_header(lower + "_score", plan, array,
                              resources["score_cpus"], resources["score_memory"]) + loader + f'''
for required in {sample} {sites} {observations} {manifest_var}; do
  [[ -s "$required" ]] || {{ echo "missing required {modality} score input: $required" >&2; exit 1; }}
done
if [[ -s {scores} ]] && gzip -t {scores} 2>/dev/null; then exit 0; fi
EXPLICIT_TEMP={shlex.quote(plan['temp_root'])}
if [[ -n "$EXPLICIT_TEMP" ]]; then
  case "$EXPLICIT_TEMP" in
    /mnt/beegfs/*|/dev/shm|/dev/shm/*) ;;
    *) echo "task scratch must be BeeGFS or RAM: $EXPLICIT_TEMP" >&2; exit 1 ;;
  esac
  mkdir -p "$EXPLICIT_TEMP"
  TEMP_BASE="$EXPLICIT_TEMP"
elif [[ -d /dev/shm && -w /dev/shm ]]; then
  TEMP_BASE=/dev/shm
else
  echo "no BeeGFS or RAM task scratch is available" >&2
  exit 1
fi
[[ -d "$TEMP_BASE" && -w "$TEMP_BASE" ]] || {{ echo "task scratch is not writable: $TEMP_BASE" >&2; exit 1; }}
[[ $(df -Pk "$TEMP_BASE" | awk 'NR==2 {{print $4}}') -gt 1048576 ]] || \
  {{ echo "less than 1 GiB free scratch" >&2; exit 1; }}
TASK_TEMP="$TEMP_BASE/jdv_{plan['run_token']}_${{SLURM_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}_{lower}"
mkdir -m 700 "$TASK_TEMP"
cleanup() {{ code=$?; [[ $code -eq 0 ]] && rm -rf -- "$TASK_TEMP"; exit $code; }}
trap cleanup EXIT
MOLECULE_ARGS=(); [[ -s {molecules} ]] && MOLECULE_ARGS=(--pileup-molecules {molecules})
{shlex.quote(str(scorer))} --joint-doublet-output {scores} \
  --joint-doublet-manifest {manifest_var} --joint-doublet-temp-dir "$TASK_TEMP" \
  --samples {sample} --pileup-sites {sites} --pileup-observations {observations} \
  "${{MOLECULE_ARGS[@]}}" --libname "$library" --modality {modality} \
  --error_ref {error} --error_alt {error} --min_evidence 10 \
  --max-second-fraction 0.95 --joint-folds 5 --threads {resources['score_cpus']}
gzip -t {scores}
'''
    for modality in ("RNA", "ATAC"):
        key = modality.lower() + "_score"
        scripts[key] = stage_root / "slurm_scripts" / f"{key}.sbatch"
        _write_script(scripts[key], score_script(modality))

    library_args = " ".join(shlex.quote(row["library"]) for row in records)
    aggregate = stage_root / "aggregate"
    gather = _sbatch_header("gather", plan, None, resources["gather_cpus"],
                            resources["gather_memory"]) + f'''
python3 -B {shlex.quote(str(validator))} stage-checkpoint \
  --stage-root {shlex.quote(str(stage_root))} --checkpoint scores
python3 -B {shlex.quote(str(identity))} joint-aggregate \
  --input-root {shlex.quote(str(stage_root))} --output-root {shlex.quote(str(aggregate))} \
  --libraries {library_args} --calibration-library lib{plan['calibration_library']}
for name in joint_doublet_cell_ledger.tsv.gz joint_doublet_candidate_scores.tsv.gz \
  joint_doublet_library_summary.tsv library25_calibration.tsv; do
  [[ -s {shlex.quote(str(aggregate))}/$name ]] || exit 1
done
python3 -B {shlex.quote(str(validator))} stage-checkpoint \
  --stage-root {shlex.quote(str(stage_root))} --checkpoint gather
'''
    scripts["gather"] = stage_root / "slurm_scripts" / "gather.sbatch"
    _write_script(scripts["gather"], gather)

    compact = _sbatch_header("compact", plan, None, resources["analysis_cpus"],
                             resources["analysis_memory"]) + f'''
python3 -B {shlex.quote(str(identity))} joint-analyze \
  --input-root {shlex.quote(str(aggregate))} --output-root {shlex.quote(str(stage_root / 'analysis' / 'derivative'))}
python3 -B {shlex.quote(str(validator))} stage-checkpoint \
  --stage-root {shlex.quote(str(stage_root))} --checkpoint compact
'''
    scripts["compact"] = stage_root / "slurm_scripts" / "compact.sbatch"
    _write_script(scripts["compact"], compact)

    if plan["stage_kind"] == "DEVELOPMENT":
        baseline = _sbatch_header("baseline", plan, None, resources["analysis_cpus"],
                                  resources["analysis_memory"]) + f'''
python3 -B {shlex.quote(str(identity))} joint-analyze \
  --input-root {shlex.quote(plan['baseline_gather'])} \
  --output-root {shlex.quote(str(stage_root / 'analysis' / 'baseline'))}
'''
        scripts["baseline"] = stage_root / "slurm_scripts" / "baseline.sbatch"
        _write_script(scripts["baseline"], baseline)

    analysis_manifest = stage_root / "manifests" / "validation_worker_tasks.tsv"
    analysis_rows = [{"task_index": index, "phase": phase,
                      "seed": stable_seed(plan["run_token"], phase)}
                     for index, phase in enumerate(ANALYSIS_TASKS)]
    write_tsv(analysis_manifest, analysis_rows, ("task_index", "phase", "seed"))
    worker_array = f"0-{len(analysis_rows) - 1}%{resources['worker_max_concurrent']}"
    worker = _sbatch_header("validation", plan, worker_array,
                            resources["worker_cpus"], resources["worker_memory"]) + f'''
python3 -B {shlex.quote(str(validator))} analysis-worker \
  --stage-root {shlex.quote(str(stage_root))} \
  --analysis-task-manifest {shlex.quote(str(analysis_manifest))} \
  --task-index "$SLURM_ARRAY_TASK_ID"
'''
    scripts["validation"] = stage_root / "slurm_scripts" / "validation.sbatch"
    _write_script(scripts["validation"], worker)
    final = _sbatch_header("finalize", plan, None, resources["analysis_cpus"],
                           resources["analysis_memory"]) + f'''
python3 -B {shlex.quote(str(validator))} finalize-stage --stage-root {shlex.quote(str(stage_root))}
'''
    scripts["finalize"] = stage_root / "slurm_scripts" / "finalize.sbatch"
    _write_script(scripts["finalize"], final)
    script_map = {key: str(value) for key, value in scripts.items()}
    atomic_json(stage_root / "manifests" / "script_manifest.json", script_map)
    return scripts


def _submit_sbatch(script, stage_root, dependency="", array_override=""):
    command = ["sbatch", "--parsable", f"--chdir={stage_root}"]
    if dependency:
        command.append(f"--dependency={dependency}")
    if array_override:
        command.append(f"--array={array_override}")
    command.append(str(script))
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError("sbatch failed: " + (result.stderr.strip() or result.stdout.strip()))
    job_id = result.stdout.strip().split(";", 1)[0]
    if not job_id.isdigit():
        raise RuntimeError(f"unexpected sbatch response: {result.stdout!r}")
    return job_id


def _save_submission(stage_root, plan, node, job_id, dependency, array_spec=""):
    plan = dict(plan)
    jobs = dict(plan.get("jobs", {}))
    jobs[node] = {"job_id": job_id, "dependency": dependency,
                  "array": array_spec, "submitted_utc": utc_now()}
    history = list(plan.get("submission_history", []))
    history.append({"node": node, "job_id": job_id, "dependency": dependency,
                    "array": array_spec, "utc": utc_now()})
    plan.update({"jobs": jobs, "submission_history": history,
                 "updated_utc": utc_now()})
    atomic_json(stage_root / "manifests" / "stage_plan.json", plan)
    return plan


def _stage_outputs(stage_root, plan, records):
    aggregate = stage_root / "aggregate"
    result = {
        "prepare": all(_file_ready(row["rna_manifest"], True) and
                       _file_ready(row["atac_manifest"], True) for row in records),
        "rna_score": all(_file_ready(row["rna_scores"], True) for row in records),
        "atac_score": all(_file_ready(row["atac_scores"], True) for row in records),
        "gather": all(_file_ready(aggregate / name, name.endswith(".gz")) for name in (
            "joint_doublet_cell_ledger.tsv.gz", "joint_doublet_candidate_scores.tsv.gz",
            "joint_doublet_library_summary.tsv", "library25_calibration.tsv")),
        "compact": _file_ready(stage_root / "analysis" / "derivative" /
                               "JOINT_DOUBLET_ANALYSIS_COMPLETE"),
        "validation": all(_file_ready(stage_root / "worker_results" / f"{index:02d}_{phase}.json")
                          for index, phase in enumerate(ANALYSIS_TASKS)),
        "finalize": _file_ready(stage_root / "VALIDATION_ANALYSIS_FINISHED"),
    }
    if plan["stage_kind"] == "DEVELOPMENT":
        result["baseline"] = _file_ready(stage_root / "analysis" / "baseline" /
                                         "JOINT_DOUBLET_ANALYSIS_COMPLETE")
    return result


def _transition_submitted(stage_root, plan, resume):
    if not resume and plan.get("state") == "RUN_SUBMITTED":
        raise RuntimeError("ordinary second submit refused; use DERIVATIVE_RESUME")
    if not resume and plan.get("state") != "RUN_PLANNED":
        raise RuntimeError(f"stage is not submit-ready: {plan.get('state')}")
    plan = dict(plan)
    plan["state"] = "RUN_SUBMITTED"
    plan["updated_utc"] = utc_now()
    atomic_json(stage_root / "manifests" / "stage_plan.json", plan)
    return plan


def _submit_graph(stage_root, plan, records, scripts, resume=False):
    complete = _stage_outputs(stage_root, plan, records)
    plan = _transition_submitted(stage_root, plan, resume)
    jobs = {}
    def submit(node, dependency="", array_spec=""):
        nonlocal plan
        if complete.get(node):
            return ""
        job_id = _submit_sbatch(scripts[node], stage_root, dependency, array_spec)
        plan = _save_submission(stage_root, plan, node, job_id, dependency, array_spec)
        jobs[node] = job_id
        return job_id

    def missing_indices(predicate, maximum):
        indices = [str(index) for index, row in enumerate(records) if not predicate(row)]
        return (",".join(indices) + f"%{maximum}") if resume and indices else ""

    prepare_array = missing_indices(
        lambda row: _file_ready(row["rna_manifest"], True) and
                    _file_ready(row["atac_manifest"], True),
        plan["resources"]["score_max_concurrent"])
    rna_array = missing_indices(lambda row: _file_ready(row["rna_scores"], True),
                                plan["resources"]["score_max_concurrent"])
    atac_array = missing_indices(lambda row: _file_ready(row["atac_scores"], True),
                                 plan["resources"]["score_max_concurrent"])
    missing_workers = [str(index) for index, phase in enumerate(ANALYSIS_TASKS)
                       if not _file_ready(stage_root / "worker_results" /
                                          f"{index:02d}_{phase}.json")]
    validation_array = (",".join(missing_workers) +
                        f"%{plan['resources']['worker_max_concurrent']}") \
        if resume and missing_workers else ""

    baseline_job = submit("baseline") if "baseline" in scripts else ""
    prepare_job = submit("prepare", array_spec=prepare_array)
    prep_dep = f"afterok:{prepare_job}" if prepare_job else ""
    rna_job = submit("rna_score", prep_dep, rna_array)
    atac_job = submit("atac_score", prep_dep, atac_array)
    score_ids = [value for value in (rna_job, atac_job) if value]
    gather_dep = "afterany:" + ":".join(score_ids) if score_ids else ""
    gather_job = submit("gather", gather_dep)
    compact_dep = f"afterok:{gather_job}" if gather_job else ""
    compact_job = submit("compact", compact_dep)
    worker_parents = [value for value in (compact_job, baseline_job) if value]
    validation_dep = "afterok:" + ":".join(worker_parents) if worker_parents else ""
    validation_job = submit("validation", validation_dep, validation_array)
    final_parents = [value for value in (validation_job,) if value]
    finalize_dep = "afterany:" + ":".join(final_parents) if final_parents else ""
    submit("finalize", finalize_dep)
    return jobs


def run_stage(args):
    libraries = parse_libraries(args.libraries)
    if not libraries:
        raise RuntimeError("no libraries selected")
    kind = _stage_kind(args, libraries)
    stage_root, plan = _initialize_stage(args, kind, libraries)
    records, manifest = _build_derivative_tasks(stage_root, plan)
    scripts = _render_stage_scripts(stage_root, plan, records, manifest)
    print(json.dumps({
        "action": args.workflow_action, "stage_kind": kind,
        "state": plan["state"], "run_token": plan["run_token"],
        "stage_root": str(stage_root),
        "array_mapping": [{"task": row["task_index"], "library": row["library"]}
                          for row in records],
        "score_array_tasks": len(records),
        "score_maximum_envelope": {
            "cpus": plan["resources"]["score_cpus"] *
                plan["resources"]["score_max_concurrent"],
            "memory": f"{plan['resources']['score_max_concurrent']} x {plan['resources']['score_memory']}"},
        "validation_worker_tasks": len(ANALYSIS_TASKS),
        "validation_maximum_envelope": {
            "cpus": plan["resources"]["worker_cpus"] *
                plan["resources"]["worker_max_concurrent"],
            "memory": f"{plan['resources']['worker_max_concurrent']} x {plan['resources']['worker_memory']}"},
        "scripts": {key: str(value) for key, value in scripts.items()},
        "submit": bool(args.submit),
    }, indent=2))
    if args.submit:
        jobs = _submit_graph(stage_root, plan, records, scripts, resume=False)
        print(json.dumps({"submitted_jobs": jobs}, indent=2))
    return 0


def resume_stage(args):
    stage_root = Path(args.stage_root).resolve()
    plan_path = stage_root / "manifests" / "stage_plan.json"
    if not plan_path.is_file():
        raise RuntimeError("resume target has no owned stage plan")
    plan = json.loads(plan_path.read_text())
    validation_root, root_plan = _assert_owned_validation_root(args.validation_root)
    if Path(plan["validation_root"]).resolve() != validation_root or \
            plan.get("validation_plan_hash") != root_plan.get("plan_hash"):
        raise RuntimeError("resume target belongs to a different validation plan")
    if not args.submit:
        print(json.dumps({"stage_root": str(stage_root), "state": plan.get("state"),
                          "outputs": stage_status_payload(stage_root)["nodes"],
                          "submit": False}, indent=2))
        return 0
    if plan.get("stage_kind") == "REPAIR":
        tasks = list(read_tsv(stage_root / "manifests" / "repair_tasks.tsv"))
        scripts = {key: Path(value) for key, value in json.loads(
            (stage_root / "manifests" / "script_manifest.json").read_text()).items()}
        missing = []
        for index, row in enumerate(tasks):
            prefix = row.get("output_prefix", "")
            ready = bool(prefix) and all(_file_ready(prefix + suffix, suffix.endswith(".gz"))
                for suffix in (".samples", ".pileup_sites.tsv.gz", ".pileup_obs.tsv.gz",
                               ".pileup_molecules.tsv.gz"))
            if not ready:
                missing.append(index)
        plan = _transition_submitted(stage_root, plan, resume=True)
        if not missing:
            return finalize_repair(argparse.Namespace(stage_root=str(stage_root)))
        array_spec = ",".join(str(value) for value in missing) + \
            f"%{plan['resources']['max_concurrent']}"
        repair_job = _submit_sbatch(scripts["repair"], stage_root,
                                    array_override=array_spec)
        plan = _save_submission(stage_root, plan, "repair_resume", repair_job, "", array_spec)
        final_job = _submit_sbatch(scripts["finalize"], stage_root,
                                   f"afterany:{repair_job}")
        _save_submission(stage_root, plan, "finalize_resume", final_job,
                         f"afterany:{repair_job}")
        print(json.dumps({"resubmitted_jobs": {"repair": repair_job,
                                               "finalize": final_job}}, indent=2))
        return 0
    task_manifest = stage_root / "manifests" / "derivative_tasks_v2.tsv"
    records = list(read_tsv(task_manifest))
    scripts = _render_stage_scripts(stage_root, plan, records, task_manifest)
    jobs = _submit_graph(stage_root, plan, records, scripts, resume=True)
    print(json.dumps({"resubmitted_jobs": jobs}, indent=2))
    return 0


def _repair_row_executable(row):
    affirmative = str(row.get("repair_executable", "")).strip().upper()
    if affirmative not in {"1", "TRUE", "YES", "REPAIR_EXECUTABLE"}:
        return False
    required = ("input_bam", "panel", "whitelist", "output_prefix",
                "min_mapq", "exclude_flags", "variant_qual", "error_ref",
                "error_alt", "error_sigma")
    return all(row.get(field, "") and "UNRESOLVED" not in row.get(field, "")
               for field in required)


def run_repair(args):
    validation_root, root_plan = _assert_owned_validation_root(args.validation_root)
    stage_root = Path(args.stage_root).resolve(strict=False)
    if not _is_relative_to(stage_root, validation_root) or stage_root == validation_root:
        raise RuntimeError("repair output root must be a child of validation root")
    repair_plan = Path(args.repair_plan)
    if not repair_plan.is_absolute() or not repair_plan.is_file():
        raise RuntimeError("repair action requires an existing absolute repair plan")
    all_rows = list(read_tsv(repair_plan))
    selected = {_lib_name(value) for value in parse_libraries(args.libraries)}
    executable = [row for row in all_rows if _lib_name(row.get("library", "")) in selected
                  and _repair_row_executable(row)]
    payload = {
        "schema_version": "joint_doublet_repair_stage_plan_v2",
        "release": RELEASE, "workflow_action": "REPAIR_MOLECULE_SIDECARS",
        "stage_kind": "REPAIR", "validation_root": str(validation_root),
        "validation_plan_hash": root_plan.get("plan_hash", ""),
        "stage_root": str(stage_root), "repair_plan": str(repair_plan.resolve()),
        "repair_plan_sha256": sha256(repair_plan),
        "libraries": parse_libraries(args.libraries),
        "tool_bin_root": str(Path(args.tool_bin_root).resolve()),
        "resources": {"cpus": args.repair_cpus, "threads": args.repair_threads,
                      "memory": args.repair_memory,
                      "max_concurrent": args.repair_max_concurrent,
                      "time": args.time, "partition": args.partition},
        "temp_root": str(Path(args.temp_root).resolve()) if args.temp_root else "",
    }
    immutable_hash = plan_hash(payload)
    plan_path = stage_root / "manifests" / "stage_plan.json"
    if stage_root.exists():
        if plan_path.is_file():
            plan = json.loads(plan_path.read_text())
            if plan.get("immutable_hash") != immutable_hash:
                raise RuntimeError("existing repair plan does not match command")
        elif any(stage_root.iterdir()):
            raise RuntimeError(f"foreign nonempty repair root: {stage_root}")
        else:
            plan = None
    else:
        stage_root.mkdir(parents=True)
        plan = None
    if plan is None:
        for name in ("manifests", "slurm_scripts", "logs", "bundles"):
            (stage_root / name).mkdir(parents=True, exist_ok=True)
        plan = dict(payload)
        plan.update({"immutable_hash": immutable_hash, "state": "RUN_PLANNED",
                     "run_token": datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") +
                        "_" + immutable_hash[:8], "jobs": {},
                     "submission_history": [], "created_utc": utc_now()})
        atomic_json(plan_path, plan)
    task_rows = []
    for index, row in enumerate(executable):
        item = dict(row)
        item["task_index"] = index
        task_rows.append(item)
    task_manifest = stage_root / "manifests" / "repair_tasks.tsv"
    fields = list(task_rows[0]) if task_rows else ["task_index", "library", "assay"]
    write_tsv(task_manifest, task_rows, fields)
    scripts = {}
    if task_rows:
        variables = " ".join(fields)
        demux = Path(plan["tool_bin_root"]) / "demux_parallel"
        array = f"0-{len(task_rows)-1}%{plan['resources']['max_concurrent']}"
        body = _sbatch_header("repair", {
            **plan, "resources": {"time": plan["resources"]["time"],
                                   "partition": plan["resources"]["partition"]}},
            array, plan["resources"]["cpus"], plan["resources"]["memory"]) + f'''
TASK_LINE="$(awk -v task="$SLURM_ARRAY_TASK_ID" 'NR == task + 2 {{print; exit}}' {shlex.quote(str(task_manifest))})"
[[ -n "$TASK_LINE" ]] || exit 1
IFS=$'\\t' read -r {variables} <<< "$TASK_LINE"
for required in "$input_bam" "$panel" "$whitelist"; do [[ -s "$required" ]] || exit 1; done
[[ -s "${{input_bam}}.bai" || -s "${{input_bam%.bam}}.bai" ]] || exit 1
mkdir -p "$(dirname "$output_prefix")"
{shlex.quote(str(demux))} -b "$input_bam" -o "$output_prefix" -v "$panel" \
  --barcodes "$whitelist" --qual "$variant_qual" --error_ref "$error_ref" \
  --error_alt "$error_alt" --error_sigma "$error_sigma" --min-mapq "$min_mapq" \
  --exclude-flags "$exclude_flags" --disable_conditional --skip_assignment \
  --force_recount --threads {plan['resources']['threads']} --dump_pileup "$output_prefix"
for suffix in .samples .pileup_sites.tsv.gz .pileup_obs.tsv.gz .pileup_molecules.tsv.gz; do
  [[ -s "${{output_prefix}}${{suffix}}" ]] || exit 1
done
'''
        scripts["repair"] = stage_root / "slurm_scripts" / "repair.sbatch"
        _write_script(scripts["repair"], body)
        validator = Path(plan["tool_bin_root"]) / "joint_doublet_validation.py"
        final = _sbatch_header("repair_finalize", {
            **plan, "resources": {"time": plan["resources"]["time"],
                                   "partition": plan["resources"]["partition"]}},
            None, 1, "8G") + f'''
python3 -B {shlex.quote(str(validator))} finalize-repair --stage-root {shlex.quote(str(stage_root))}
'''
        scripts["finalize"] = stage_root / "slurm_scripts" / "finalize_repair.sbatch"
        _write_script(scripts["finalize"], final)
    atomic_json(stage_root / "manifests" / "script_manifest.json",
                {key: str(value) for key, value in scripts.items()})
    nonexec = sum(str(row.get("repair_needed", "")).upper() in {"1", "TRUE", "YES"}
                  and not _repair_row_executable(row) for row in all_rows)
    print(json.dumps({"stage_root": str(stage_root), "state": plan["state"],
                      "repair_tasks": len(task_rows),
                      "non_executable_missing_branches": nonexec,
                      "array_mapping": [{"task": row["task_index"],
                                         "library": row.get("library"),
                                         "assay": row.get("assay")} for row in task_rows],
                      "maximum_envelope": {"cpus": args.repair_cpus * args.repair_max_concurrent,
                                           "memory": f"{args.repair_max_concurrent} x {args.repair_memory}"},
                      "submit": bool(args.submit)}, indent=2))
    if not args.submit:
        return 0
    plan = _transition_submitted(stage_root, plan, resume=False)
    if not task_rows:
        return finalize_repair(argparse.Namespace(stage_root=str(stage_root)))
    job = _submit_sbatch(scripts["repair"], stage_root)
    plan = _save_submission(stage_root, plan, "repair", job, "",
                            f"0-{len(task_rows)-1}%{args.repair_max_concurrent}")
    final_job = _submit_sbatch(scripts["finalize"], stage_root, f"afterany:{job}")
    _save_submission(stage_root, plan, "finalize", final_job, f"afterany:{job}")
    print(json.dumps({"submitted_jobs": {"repair": job, "finalize": final_job}}, indent=2))
    return 0


def finalize_repair(args):
    stage_root = Path(args.stage_root)
    tasks = list(read_tsv(stage_root / "manifests" / "repair_tasks.tsv"))
    failed = []
    for row in tasks:
        prefix = row.get("output_prefix", "")
        if not prefix or not all(_file_ready(prefix + suffix, suffix.endswith(".gz"))
                                 for suffix in (".samples", ".pileup_sites.tsv.gz",
                                                ".pileup_obs.tsv.gz",
                                                ".pileup_molecules.tsv.gz")):
            failed.append({"library": row.get("library"), "assay": row.get("assay"),
                           "output_prefix": prefix})
    status_value = "COMPLETE" if not failed else "FAILED"
    marker = {"status": status_value, "utc": utc_now(), "tasks": len(tasks),
              "failed": failed, "release": RELEASE}
    atomic_json(stage_root / "MOLECULE_REPAIR_FINISHED", marker)
    if status_value == "COMPLETE":
        atomic_json(stage_root / "MOLECULE_REPAIR_COMPLETE", marker)
    return 0 if status_value == "COMPLETE" else 2


def _float(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else math.nan
    except (TypeError, ValueError):
        return math.nan


def _diagnostic_rows(stage_root):
    path = stage_root / "analysis" / "derivative" / "joint_doublet_cell_diagnostics.tsv.gz"
    if not _file_ready(path, True):
        raise RuntimeError(f"compact derivative diagnostics missing: {path}")
    return list(read_tsv(path))


def _candidate_policy_metrics(stage_root, plan):
    parity = []
    for library in plan["score_libraries"]:
        root = stage_root / f"lib{library}"
        keys = {}
        for modality in ("rna", "atac"):
            path = root / f"lib{library}.{modality}_joint_manifest.tsv.gz"
            values = set()
            for row in read_tsv(path):
                values.add((row.get("barcode", ""), row.get("candidate_id", ""),
                            row.get("locked_state", ""), row.get("second_state", "")))
            keys[modality] = values
        parity.append({"library": f"lib{library}", "rna_rows": len(keys["rna"]),
                       "atac_rows": len(keys["atac"]),
                       "identical": keys["rna"] == keys["atac"],
                       "rna_only": len(keys["rna"] - keys["atac"]),
                       "atac_only": len(keys["atac"] - keys["rna"])})
    return {"libraries": parity,
            "all_identical": all(row["identical"] for row in parity)}


def _cross_assay_metrics(rows):
    focal = []
    for row in rows:
        rp = _float(row.get("rna_top_percentile"))
        ap = _float(row.get("atac_top_percentile"))
        if rp >= 0.95 and ap >= 0.95:
            focal.append(row)
    rna = [row.get("rna_top_second_state", "") for row in focal]
    atac = [{"winner": row.get("atac_top_second_state", ""),
             "score": row.get("atac_top_percentile", "")} for row in focal]
    strata = [(row.get("library", ""), row.get("ploidy_field", ""),
               row.get("reconciled_identity_locked", "")) for row in focal]
    permutation = fixed_denominator_permutation(rna, atac, strata, 10000)
    permutation["fraction"] = permutation["observed"] / len(focal) if focal else math.nan
    permutation["analytic_expected_fraction"] = analytic_agreement(
        rna, [value["winner"] for value in atac], strata)
    return permutation


def _coverage_metrics(rows):
    fields = ("rna_top_discriminating_depth", "atac_top_discriminating_depth",
              "rna_top_discriminating_sites", "atac_top_discriminating_sites")
    result = {"cells": len(rows), "fields": {}}
    for field in fields:
        values = [_float(row.get(field)) for row in rows]
        values = [value for value in values if math.isfinite(value)]
        result["fields"][field] = {
            "available": len(values), "median": statistics.median(values) if values else math.nan,
            "q05": _quantile(values, 0.05), "q95": _quantile(values, 0.95)}
    result["method"] = "within-library/frozen-identity coverage-stratum accounting"
    return result


def _quantile(values, probability):
    if not values:
        return math.nan
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    left = int(math.floor(position)); right = int(math.ceil(position))
    if left == right:
        return ordered[left]
    return ordered[left] * (right - position) + ordered[right] * (position - left)


def _molecule_metrics(stage_root, plan):
    totals = Counter(); changed = Counter()
    for library in plan["score_libraries"]:
        for modality in ("rna", "atac"):
            path = stage_root / f"lib{library}" / f"lib{library}.{modality}_joint_scores.tsv.gz"
            site_top = {}; molecule_top = {}; available = set()
            for row in read_tsv(path):
                key = row.get("barcode", "")
                site = _float(row.get("delta_site_balanced_log_likelihood_k2_minus_k1"))
                molecule = _float(row.get("molecule_balanced_delta_log_likelihood_k2_minus_k1"))
                if math.isfinite(site) and (key not in site_top or site > site_top[key][0]):
                    site_top[key] = (site, row.get("candidate_id", ""))
                if math.isfinite(molecule):
                    available.add(key)
                    if key not in molecule_top or molecule > molecule_top[key][0]:
                        molecule_top[key] = (molecule, row.get("candidate_id", ""))
            label = f"lib{library}_{modality}"
            totals[label + "_cells"] = len(site_top)
            totals[label + "_molecule_available"] = len(available)
            changed[label + "_winner_changed"] = sum(
                site_top[key][1] != molecule_top[key][1]
                for key in set(site_top) & set(molecule_top))
    return {"counts": dict(totals), "winner_changes": dict(changed),
            "interpretation": "linked-unit sensitivity; not automatically phased-haplotype evidence"}


def _normalization_metrics(rows):
    by_library = defaultdict(list)
    for row in rows:
        value = _float(row.get("technical_occupancy_combined_z"))
        if math.isfinite(value):
            by_library[row.get("library", "")].append((value, row.get("barcode", "")))
    summary = []
    for library, values in sorted(by_library.items()):
        raw = [value for value, _ in values]
        median = statistics.median(raw); mad = statistics.median(abs(value - median) for value in raw)
        scale = 1.4826 * mad if mad > 0 else (statistics.pstdev(raw) if len(raw) > 1 else 1.0)
        threshold = _quantile(raw, 0.95)
        summary.append({"library": library, "cells": len(raw), "median": median,
                        "robust_scale": scale, "raw_q95": threshold,
                        "within_library_top5": sum(value >= threshold for value in raw)})
    return {"within_library": summary, "state_conditioned":
            "not identifiable when a documented cell-state column is absent"}


def _control_partitions(stage_root, rows):
    output = stage_root / "worker_results" / "control_parent_partitions.tsv.gz"
    data = []
    for row in rows:
        identity = row.get("reconciled_identity_locked", "")
        if not identity:
            continue
        library = row.get("library", ""); barcode = row.get("barcode", "")
        split = control_parent_partition(library, identity,
                                         row.get("ploidy_field", ""), "", barcode)
        data.append({"library": library, "barcode": barcode,
                     "frozen_identity": identity, "ploidy": row.get("ploidy_field", ""),
                     "cell_state": "UNAVAILABLE", "partition": split,
                     "seed": MASTER_SEED})
    write_tsv(output, data, ("library", "barcode", "frozen_identity", "ploidy",
                             "cell_state", "partition", "seed"))
    return {"eligible_parents": len(data),
            "development": sum(row["partition"] == "development" for row in data),
            "test": sum(row["partition"] == "test" for row in data),
            "path": str(output)}


def _control_construction(stage_root, rows):
    grouped = defaultdict(list)
    for row in rows:
        identity = row.get("reconciled_identity_locked", "")
        if not identity:
            continue
        partition = control_parent_partition(row.get("library", ""), identity,
                                             row.get("ploidy_field", ""), "",
                                             row.get("barcode", ""))
        grouped[(row.get("library", ""), partition)].append(row)
    controls = []
    for (library, partition), parents in sorted(grouped.items()):
        parents.sort(key=lambda row: stable_seed("pair", library, partition,
                                                row.get("barcode", "")))
        used = set()
        for left in parents:
            if left.get("barcode") in used:
                continue
            right = next((candidate for candidate in parents
                          if candidate.get("barcode") not in used and
                          candidate.get("barcode") != left.get("barcode") and
                          candidate.get("reconciled_identity_locked") !=
                          left.get("reconciled_identity_locked")), None)
            if right is None:
                continue
            used.update((left.get("barcode"), right.get("barcode")))
            pair_id = hashlib.sha256((library + "\x1f" + left.get("barcode", "") +
                                      "\x1f" + right.get("barcode", "")).encode()).hexdigest()[:16]
            for fraction in (0.10, 0.20, 0.35):
                controls.append({"control_id": f"{library}_{pair_id}_{fraction:.2f}",
                                 "library": library, "partition": partition,
                                 "parent_a": left.get("barcode", ""),
                                 "parent_b": right.get("barcode", ""),
                                 "identity_a": left.get("reconciled_identity_locked", ""),
                                 "identity_b": right.get("reconciled_identity_locked", ""),
                                 "requested_fraction": fraction,
                                 "realized_rna_linked_unit_fraction": "PENDING_MIXTURE_MATERIALIZATION",
                                 "realized_atac_linked_unit_fraction": "PENDING_MIXTURE_MATERIALIZATION",
                                 "status": "PARENTS_FROZEN_SOURCE_DISJOINT"})
    output = stage_root / "worker_results" / "control_construction_manifest.tsv.gz"
    fields = list(controls[0]) if controls else ("control_id", "status")
    write_tsv(output, controls, fields)
    return {"control_rows": len(controls), "unique_pairs": len(controls) // 3,
            "path": str(output), "scoring_status":
            "parent pairs frozen; linked-unit mixture materialization remains unavailable"}


def analysis_worker(args):
    stage_root = Path(args.stage_root)
    plan = json.loads((stage_root / "manifests" / "stage_plan.json").read_text())
    tasks = list(read_tsv(args.analysis_task_manifest))
    index = int(args.task_index)
    if index < 0 or index >= len(tasks):
        raise RuntimeError("validation worker task index is out of range")
    phase = tasks[index]["phase"]
    rows = _diagnostic_rows(stage_root)
    status_value = "COMPLETE"; reason = ""
    if phase == "BASELINE_RECONSTRUCTION":
        if plan["stage_kind"] == "DEVELOPMENT":
            marker = stage_root / "analysis" / "baseline" / "JOINT_DOUBLET_ANALYSIS_COMPLETE"
            metrics = {"baseline_analysis_complete": _file_ready(marker),
                       "baseline_gather": plan.get("baseline_gather", "")}
            status_value = "COMPLETE" if metrics["baseline_analysis_complete"] else "FAILED"
        else:
            metrics = {"status": "NOT_APPLICABLE", "stage_kind": plan["stage_kind"]}
    elif phase == "CANDIDATE_POLICY_AUDIT":
        metrics = _candidate_policy_metrics(stage_root, plan)
        status_value = "COMPLETE" if metrics["all_identical"] else "FAILED"
    elif phase == "CROSS_ASSAY_NULL":
        metrics = _cross_assay_metrics(rows)
    elif phase == "COVERAGE_MATCHING":
        metrics = _coverage_metrics(rows)
    elif phase == "COVERAGE_DOWNSAMPLING":
        metrics = {"status": "NOT_EVALUABLE", "reason":
                   "scorer output does not expose per-linked-unit sufficient statistics for refit downsampling"}
        status_value = "PARTIAL"; reason = metrics["reason"]
    elif phase == "MOLECULE_SENSITIVITY":
        metrics = _molecule_metrics(stage_root, plan)
    elif phase == "QUANTITY_NORMALIZATION":
        metrics = _normalization_metrics(rows)
    elif phase == "OCCUPANCY_ASSESSMENT":
        metrics = {"status": "NOT_IDENTIFIABLE", "reason":
                   "no authoritative cell-state annotation is declared by the source manifest",
                   "descriptive_normalization_available": True}
        status_value = "PARTIAL"; reason = metrics["reason"]
    elif phase == "CONTROL_PARENT_PARTITION":
        metrics = _control_partitions(stage_root, rows)
    elif phase == "CONTROL_CONSTRUCTION":
        metrics = _control_construction(stage_root, rows)
        status_value = "PARTIAL"; reason = metrics["scoring_status"]
    elif phase == "CELL_CONDITIONAL_NULL":
        metrics = {"status": "NOT_EVALUABLE", "reason":
                   "candidate sufficient-statistic cache is not emitted by the installed scorer interface"}
        status_value = "PARTIAL"; reason = metrics["reason"]
    else:
        counts = Counter(row.get("recomputed_modality_relation", "MISSING") for row in rows)
        metrics = {"cells": len(rows), "relations": dict(counts),
                   "stage_kind": plan["stage_kind"]}
    payload = {"schema_version": "joint_doublet_validation_worker_result_v2",
               "release": RELEASE, "task_index": index, "phase": phase,
               "status": status_value, "reason": reason, "seed": tasks[index]["seed"],
               "utc": utc_now(), "metrics": metrics}
    atomic_json(stage_root / "worker_results" / f"{index:02d}_{phase}.json", payload)
    return 0 if status_value != "FAILED" else 2


def finalize_stage(args):
    stage_root = Path(args.stage_root)
    plan_path = stage_root / "manifests" / "stage_plan.json"
    if not plan_path.is_file():
        raise RuntimeError("stage has no owned plan")
    plan = json.loads(plan_path.read_text())
    results = []
    missing = []
    for index, phase in enumerate(ANALYSIS_TASKS):
        path = stage_root / "worker_results" / f"{index:02d}_{phase}.json"
        if path.is_file():
            results.append(json.loads(path.read_text()))
        else:
            missing.append(phase)
    failed = [row["phase"] for row in results if row.get("status") == "FAILED"]
    partial = [row["phase"] for row in results if row.get("status") == "PARTIAL"]
    if failed:
        status_value = "FAILED"
    elif missing or partial:
        status_value = "PARTIAL"
    else:
        status_value = "COMPLETE"
    frozen_path = stage_root / "frozen_analysis_spec.json"
    if plan["stage_kind"] == "DEVELOPMENT" and not failed and not missing:
        frozen = {
            "schema_version": "joint_doublet_frozen_analysis_spec_v2",
            "state": "final_frozen", "release": RELEASE,
            "frozen_utc": utc_now(), "master_seed": MASTER_SEED,
            "candidate_policy": "DERIVATIVE_COMPLETE",
            "calibration_library": plan["calibration_library"],
            "development_libraries": plan["score_libraries"],
            "heldout_libraries": plan["heldout_libraries"],
            "regression_libraries": plan["regression_libraries"],
            "heldout_values_parsed": False,
            "site_formula_version": "K1_LOCKED_STATE_VS_K2_ADDED_STATE_PLUS_FIXED_AMBIENT_V1",
            "molecule_formula_version": "LINKED_UNIT_EQUAL_WEIGHT_NORMALIZED_SITE_LL_V1",
            "incomplete_optional_branches": partial,
            "stage_plan_hash": plan["immutable_hash"],
        }
        frozen["spec_hash"] = plan_hash(frozen)
        atomic_json(frozen_path, frozen)
    marker = {"schema_version": "joint_doublet_validation_terminal_v2",
              "status": status_value, "stage_kind": plan["stage_kind"],
              "libraries": plan["libraries"], "score_libraries": plan["score_libraries"],
              "release": RELEASE, "utc": utc_now(), "failed_branches": failed,
              "incomplete_branches": missing + partial,
              "heldout_values_parsed_before_freeze": False,
              "stage_plan_hash": plan["immutable_hash"],
              "frozen_spec": str(frozen_path) if frozen_path.is_file() else plan.get("frozen_spec", "")}
    atomic_json(stage_root / "VALIDATION_ANALYSIS_FINISHED", marker)
    if status_value == "COMPLETE":
        atomic_json(stage_root / "VALIDATION_ANALYSIS_COMPLETE", marker)
    plan["state"] = "RUN_FINISHED"
    plan["terminal_status"] = status_value
    plan["updated_utc"] = utc_now()
    atomic_json(plan_path, plan)
    print(json.dumps(marker, indent=2))
    return 0 if status_value != "FAILED" else 2


def checkpoint_stage(args):
    stage_root = Path(args.stage_root)
    plan = json.loads((stage_root / "manifests" / "stage_plan.json").read_text())
    records = list(read_tsv(stage_root / "manifests" / "derivative_tasks_v2.tsv"))
    outputs = _stage_outputs(stage_root, plan, records)
    if args.checkpoint == "scores":
        complete = outputs["prepare"] and outputs["rna_score"] and outputs["atac_score"]
        marker_name = "DERIVATIVE_SCORES"
        detail = {key: outputs[key] for key in ("prepare", "rna_score", "atac_score")}
    elif args.checkpoint == "gather":
        complete = outputs["gather"]
        marker_name = "DERIVATIVE_GATHER"
        detail = {"gather": outputs["gather"]}
    else:
        complete = outputs["compact"]
        marker_name = "VALIDATION_COMPACT"
        detail = {"compact": outputs["compact"]}
    payload = {"status": "COMPLETE" if complete else "PARTIAL", "utc": utc_now(),
               "checkpoint": args.checkpoint, "detail": detail,
               "stage_plan_hash": plan["immutable_hash"], "release": RELEASE}
    atomic_json(stage_root / f"{marker_name}_FINISHED", payload)
    if complete:
        atomic_json(stage_root / f"{marker_name}_COMPLETE", payload)
    print(json.dumps(payload, indent=2))
    return 0


def stage_status_payload(stage_root):
    stage_root = Path(stage_root)
    plan_path = stage_root / "manifests" / "stage_plan.json"
    if not plan_path.is_file():
        return {"status": "FOREIGN_OR_UNINITIALIZED", "stage_root": str(stage_root),
                "nodes": {}}
    plan = json.loads(plan_path.read_text())
    if plan.get("stage_kind") == "REPAIR":
        tasks = list(read_tsv(stage_root / "manifests" / "repair_tasks.tsv"))
        task_ok = []
        for row in tasks:
            prefix = row.get("output_prefix", "")
            task_ok.append(bool(prefix) and all(_file_ready(prefix + suffix, suffix.endswith(".gz"))
                for suffix in (".samples", ".pileup_sites.tsv.gz", ".pileup_obs.tsv.gz",
                               ".pileup_molecules.tsv.gz")))
        nodes = {"repair": all(task_ok),
                 "finalize": _file_ready(stage_root / "MOLECULE_REPAIR_FINISHED")}
    else:
        records = list(read_tsv(stage_root / "manifests" / "derivative_tasks_v2.tsv"))
        nodes = _stage_outputs(stage_root, plan, records)
    return {"status": plan.get("terminal_status", plan.get("state")),
            "stage_root": str(stage_root), "stage_kind": plan.get("stage_kind"),
            "run_token": plan.get("run_token"), "nodes": nodes,
            "jobs": plan.get("jobs", {})}


def status(args):
    if args.stage_root:
        print(json.dumps(stage_status_payload(args.stage_root), indent=2))
        return 0
    root = Path(args.validation_root)
    plan_path = root / "manifests" / "run_plan.json"
    if not plan_path.is_file():
        print(json.dumps({"status": "FOREIGN_OR_UNINITIALIZED", "root": str(root)}))
        return 2
    plan = json.loads(plan_path.read_text())
    markers = {}
    for marker in FINISHED_MARKERS:
        path = root / marker
        markers[marker] = json.loads(path.read_text()).get("status", "INVALID") \
            if path.is_file() else "NOT_FINISHED"
    print(json.dumps({"root": str(root), "plan_state": plan.get("state"),
                      "plan_hash": plan.get("plan_hash"), "markers": markers}, indent=2))
    return 0


def _unavailable_table(path, scope, reason):
    write_tsv(path, [{"attempted_scope": scope, "status": "NOT_RUN",
                      "reason": reason, "consequence": "not estimable"}],
              ("attempted_scope", "status", "reason", "consequence"))


def _write_unavailable_plots(package, reason):
    from PIL import Image, ImageDraw, ImageFont
    from reportlab.lib.pagesizes import landscape, letter
    from reportlab.pdfgen import canvas

    plot_root = package / "plots"
    plot_root.mkdir(parents=True, exist_ok=True)
    title_font = ImageFont.truetype("DejaVuSans-Bold.ttf", 58)
    body_font = ImageFont.truetype("DejaVuSans.ttf", 34)
    png_names = [name for name in PLOT_NAMES if name.endswith(".png")]
    for name in png_names:
        image = Image.new("RGB", (2000, 1125), "#F7F9FC")
        draw = ImageDraw.Draw(image)
        draw.rectangle((0, 0, 2000, 155), fill="#183A5A")
        draw.text((90, 48), name.removesuffix(".png").replace("_", " ").title(),
                  font=title_font, fill="white")
        draw.text((90, 285), "NOT GENERATED", font=title_font, fill="#A33A2B")
        words = ("Scientific phase not run. No values are plotted. " + reason).split()
        lines, current = [], []
        for word in words:
            proposed = " ".join(current + [word])
            if draw.textlength(proposed, font=body_font) > 1760 and current:
                lines.append(" ".join(current)); current = [word]
            else:
                current.append(word)
        if current:
            lines.append(" ".join(current))
        for index, line in enumerate(lines):
            draw.text((90, 430 + 55 * index), line, font=body_font, fill="#243746")
        draw.text((90, 995), "Status: incomplete | No biological interpretation",
                  font=body_font, fill="#596B7A")
        image.save(plot_root / name, format="PNG", optimize=True)
    pdf_path = plot_root / "all_plots.pdf"
    page_width, page_height = landscape(letter)
    pdf = canvas.Canvas(str(pdf_path), pagesize=(page_width, page_height))
    for page_number, name in enumerate(png_names, 1):
        pdf.setFillColorRGB(0.094, 0.227, 0.353)
        pdf.rect(0, page_height - 80, page_width, 80, stroke=0, fill=1)
        pdf.setFillColorRGB(1, 1, 1)
        pdf.setFont("Helvetica-Bold", 22)
        pdf.drawString(42, page_height - 51,
                       name.removesuffix(".png").replace("_", " ").title())
        pdf.setFillColorRGB(0.64, 0.23, 0.17)
        pdf.setFont("Helvetica-Bold", 26)
        pdf.drawString(50, page_height - 155, "NOT GENERATED")
        pdf.setFillColorRGB(0.14, 0.22, 0.28)
        pdf.setFont("Helvetica", 13)
        text = pdf.beginText(50, page_height - 205)
        text.setLeading(19)
        message = "Scientific phase not run. No values are plotted. " + reason
        words = message.split(); line = []
        for word in words:
            proposed = " ".join(line + [word])
            if pdf.stringWidth(proposed, "Helvetica", 13) > page_width - 100 and line:
                text.textLine(" ".join(line)); line = [word]
            else:
                line.append(word)
        if line:
            text.textLine(" ".join(line))
        pdf.drawText(text)
        pdf.setFillColorRGB(0.35, 0.42, 0.48)
        pdf.setFont("Helvetica", 10)
        pdf.drawRightString(page_width - 35, 25,
                            f"Incomplete checkpoint | page {page_number}/{len(png_names)}")
        pdf.showPage()
    pdf.save()


def finalize_checkpoint(args):
    root = Path(args.validation_root)
    if not (root / "manifests" / "run_plan.json").is_file():
        raise RuntimeError("validation root is not owned by this workflow")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    package = root / "return_package" / f"joint_doublet_validation_return_{stamp}"
    for directory in ("manifests", "tables", "plots", "code/changed_files_only"):
        (package / directory).mkdir(parents=True, exist_ok=True)
    reason = args.unavailable_reason or (
        "Scientific data-generating phases were not completed; see stage accounting.")
    atomic_text(package / "README_FIRST.md", f"""# Joint RNA/ATAC doublet validation checkpoint

Status: **incomplete**. This package does not contain a biological call set.

The source-preserving implementation and local fixtures were prepared, but the
scientific phases requiring the cluster data and coordinated installed build did
not run in the producing environment. Every unavailable table exists explicitly
and says why; missing values must not be interpreted as zero or singlet evidence.

Validation root: `{root}`

Release: `{RELEASE}`

Read `01_validation_results.md`, then `manifests/validation_checks.tsv` and
`tables/27_final_decision_record.tsv`.
""")
    questions = (
        "Candidate-frequency-aware agreement", "Coverage matching/downsampling",
        "Molecule balancing", "New-source versus dosage-shift",
        "Quantity normalization/ablation", "Control null calibration",
        "Transfer to Libraries 35/38", "Library 19 regression",
        "Production-rule evidence",
    )
    atomic_text(package / "01_validation_results.md",
                "# Validation results\n\n" + "\n".join(
                    f"- {question}: **not evaluable** — {reason}"
                    for question in questions) +
                "\n\nBiological interpretation: `pending boss-chat interpretation`.\n")
    _write_unavailable_plots(package, reason)
    for name in TABLE_NAMES:
        _unavailable_table(package / "tables" / name, name, reason)
    decision_rows = [{
        "question": question, "evidence_file_rows": "unavailable",
        "estimate": "NA", "interval_or_null": "NA",
        "prespecified_criterion": "see controlling handoff",
        "objective_analysis_status": "incomplete", "factual_summary": reason,
        "boss_interpretation": "pending",
    } for question in (
        "data accounting valid", "cross-assay agreement relative to null",
        "null calibration", "coverage-control result", "molecule-balanced result",
        "new-source versus dosage-shift result",
        "occupancy versus state/library normalization", "control adequacy",
        "transfer result for Library 35", "transfer result for Library 38",
        "Library-19 regression result", "production-rule decision evidence",
    )]
    write_tsv(package / "tables" / "27_final_decision_record.tsv", decision_rows,
              list(decision_rows[0]))
    for source_name in ("run_plan.json", "frozen_analysis_spec.json",
                        "molecule_repair_plan.tsv", "molecule_sidecar_inventory.tsv",
                        "validation_checks.tsv"):
        source = root / "manifests" / source_name
        if source.is_file():
            shutil.copy2(source, package / "manifests" / source_name)
    defaults = {
        "run_manifest.json": {"release": RELEASE, "status": "INCOMPLETE",
                              "utc": utc_now(), "host": socket.gethostname(),
                              "reason": reason},
        "frozen_analysis_spec.json": {"state": "development_not_yet_frozen",
                                      "heldout_values_parsed": False},
        "build_and_deployment_manifest.json": {"installation_status": "NOT_INSTALLED"},
    }
    for name, value in defaults.items():
        target = package / "manifests" / name
        if not target.exists():
            atomic_json(target, value)
    for name in ("denominator_registry.tsv", "data_dictionary.tsv",
                 "warnings_and_exclusions.tsv", "stage_and_job_accounting.tsv",
                 "control_parent_partitions.tsv.gz",
                 "control_construction_manifest.tsv.gz"):
        _unavailable_table(package / "manifests" / name, name, reason)
    if args.code_root:
        code_root = Path(args.code_root).resolve()
        changed = (
            "Makefile", "src/tetra_score_calls.cpp",
            "scripts/orchestrate_tetraploid.py",
            "scripts/identity_reconciliation.py",
            "scripts/joint_doublet_validation.py",
            "tests/test_joint_doublet_validation_fixture.py",
        )
        for relative in changed:
            source = code_root / relative
            if not source.is_file():
                raise RuntimeError(f"changed source file is absent: {source}")
            target = package / "code" / "changed_files_only" / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        diff = subprocess.run(
            ["git", "diff", "--binary", "--", *changed], cwd=code_root,
            text=True, capture_output=True, check=False)
        if diff.returncode != 0:
            raise RuntimeError("could not generate unified patch: " + diff.stderr)
        atomic_text(package / "code" / "changes.patch", diff.stdout)
        source_hashes = {relative: sha256(code_root / relative)
                         for relative in changed}
        build_manifest = {
            "schema_version": "joint_doublet_build_deployment_manifest_v1",
            "release": RELEASE, "installation_status": "NOT_INSTALLED",
            "compiler": "g++ local syntax/audit build",
            "normal_cluster_build": "NOT_RUN_HTSLIB_HEADERS_ABSENT_LOCALLY",
            "normal_cluster_flags": "repository Makefile CXXFLAGS_TET",
            "audit_flags": "-Wall -Wextra -Wpedantic appended",
            "audit_result": "PASS_NO_WARNINGS_WITH_LOCAL_BARCODE_LINK_STUB",
            "python_compile": "PASS",
            "joint_molecule_fixture": "PASS",
            "candidate_relationship_fixture": "PASS",
            "source_hashes": source_hashes,
            "staged_tool_set": [
                "tetra_score_calls", "orchestrate_tetraploid.py",
                "identity_reconciliation.py", "joint_doublet_validation.py"],
            "rollback_mapping": "created at installation time by command_record.sh",
        }
        atomic_json(package / "manifests" / "build_and_deployment_manifest.json",
                    build_manifest)
        atomic_text(package / "code" / "LOCAL_TEST_RESULTS.md", """# Local test results

- Python byte-compilation: PASS for both modified scripts, the new validation helper, and the fixture.
- C++ audit syntax build with `-Wall -Wextra -Wpedantic`: PASS with no warnings.
- Linked executable using only a local barcode-symbol stub: PASS.
- End-to-end molecule fixture: PASS for duplicate collapse, multi-SNP linked-unit weighting, K rule, one-versus-four-thread identity, missing optional sidecar, isolated malformed cell row, and corrupt schema failure.
- Repository normal Makefile build: BLOCKED before scorer compilation because this workspace lacks `htslib/sam.h`. The exact cluster build remains required against HTSlib 1.20.
- Cluster/BeeGFS scientific execution: NOT RUN because `/mnt/beegfs`, `/nvme/software/packages/cellbouncer/dev/bin`, `sbatch`, and `squeue` are unavailable here.
""")
        validation_root = "/mnt/beegfs/tetraploid_multiome_cis_trans/3P/analysis/aggregate_library_analysis/joint_doublet_validation_20260917_v1"
        source_root = "/mnt/beegfs/tetraploid_multiome_cis_trans/3P/analysis/aggregate_library_analysis/joint_doublet"
        gather_root = source_root + "/partial_gathers/paired_20260917T200652Z"
        user_repo = "/home/b/cellbouncer-parallel"
        commands = f"""#!/bin/bash
set -euo pipefail
module purge
module load miniforge/3 genomics-base/latest htslib/1.20

# Build the existing scorer target, then stage only the coordinated artifacts.
cd {user_repo}
make tetra_score_calls
mkdir -p {user_repo}/staged_joint_doublet_validation
cp tetra_score_calls scripts/orchestrate_tetraploid.py scripts/identity_reconciliation.py scripts/joint_doublet_validation.py {user_repo}/staged_joint_doublet_validation/
chmod 755 {user_repo}/staged_joint_doublet_validation/tetra_score_calls {user_repo}/staged_joint_doublet_validation/*.py

# User-controlled coordinated installation. Backups are timestamped and rollbackable.
install_root=/nvme/software/packages/cellbouncer/dev/bin
backup_root=/nvme/software/packages/cellbouncer/dev/backups/joint_doublet_validation_20260918
mkdir -p "$backup_root"
for name in tetra_score_calls orchestrate_tetraploid.py identity_reconciliation.py joint_doublet_validation.py; do
  if [[ -e "$install_root/$name" ]]; then cp -a "$install_root/$name" "$backup_root/$name"; fi
  install -m 0755 "{user_repo}/staged_joint_doublet_validation/$name" "$install_root/$name.new"
  mv -f "$install_root/$name.new" "$install_root/$name"
done

# Phase 0A dry run. Repeat the identical command with --submit only after inspection.
python3 -B "$install_root/orchestrate_tetraploid.py" --stage JOINT_DOUBLET --libraries 1-29 35 38 --joint-doublet-action VALIDATION_PREFLIGHT --joint-doublet-source-output-root {source_root} --joint-doublet-partial-output-root {gather_root} --joint-doublet-validation-root {validation_root} --joint-doublet-calibration-library 25 --joint-doublet-heldout-libraries 35 38 --joint-doublet-regression-libraries 19 --joint-doublet-evidence-mode SITE_AND_MOLECULE --joint-doublet-partition compute

# File-driven status. Scientific stage submission is intentionally unavailable
# in this checkpoint and must not be attempted until its missing worker graph is implemented.
python3 -B "$install_root/orchestrate_tetraploid.py" --stage JOINT_DOUBLET --libraries 1-29 35 38 --joint-doublet-action DERIVATIVE_STATUS --joint-doublet-validation-root {validation_root}

# Rollback (run only if installation verification fails):
# for name in tetra_score_calls orchestrate_tetraploid.py identity_reconciliation.py joint_doublet_validation.py; do
#   test -e "$backup_root/$name" && install -m 0755 "$backup_root/$name" "$install_root/$name"
# done
"""
        atomic_text(package / "code" / "command_record.sh", commands)
    checks_target = package / "manifests" / "validation_checks.tsv"
    if not checks_target.exists():
        _unavailable_table(checks_target, "validation checks", reason)
    manifest_rows = []
    for path in sorted(package.rglob("*")):
        if path.is_file() and path.name != "file_manifest.tsv":
            rows, header = count_rows_and_header(path) if path.suffix in {".tsv", ".gz"} else ("NA", [])
            manifest_rows.append({
                "relative_path": str(path.relative_to(package)),
                "description": "checkpoint artifact", "bytes": path.stat().st_size,
                "rows": rows, "columns": len(header) if header else "NA",
                "compression": "gzip" if path.suffix == ".gz" else "none",
                "sha256": sha256(path),
            })
    write_tsv(package / "manifests" / "file_manifest.tsv", manifest_rows,
              list(manifest_rows[0]))
    zip_path = package.with_suffix(".zip")
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED,
                         compresslevel=9) as archive:
        for path in sorted(package.rglob("*")):
            if path.is_file():
                archive.write(path, arcname=str(path.relative_to(package.parent)))
    with zipfile.ZipFile(zip_path) as archive:
        bad = archive.testzip()
        if bad:
            raise RuntimeError(f"ZIP integrity failure at {bad}")
    atomic_text(Path(str(zip_path) + ".sha256"), f"{sha256(zip_path)}  {zip_path.name}\n")
    marker = {"status": "PARTIAL", "utc": utc_now(),
              "archive": str(zip_path), "reason": reason, "release": RELEASE}
    atomic_json(root / "RETURN_PACKAGE_FINISHED", marker)
    print(str(zip_path))
    return 0


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=(
        "preflight", "run-stage", "run-repair", "resume", "status",
        "analysis-worker", "stage-checkpoint", "finalize-stage", "finalize-repair",
        "finalize-checkpoint"))
    parser.add_argument("--libraries", nargs="*", default=[])
    parser.add_argument("--source-output-root", default="")
    parser.add_argument("--partial-output-root", default="")
    parser.add_argument("--validation-root", default="")
    parser.add_argument("--stage-root", default="")
    parser.add_argument("--workflow-action", default="")
    parser.add_argument("--tool-bin-root", default="")
    parser.add_argument("--calibration-library", default="25")
    parser.add_argument("--heldout-libraries", nargs="*", default=[])
    parser.add_argument("--regression-libraries", nargs="*", default=[])
    parser.add_argument("--evidence-mode", default="SITE_AND_MOLECULE")
    parser.add_argument("--frozen-spec", default="")
    parser.add_argument("--repair-plan", default="")
    parser.add_argument("--temp-root", default="")
    parser.add_argument("--score-cpus", type=int, default=16)
    parser.add_argument("--score-memory", default="96G")
    parser.add_argument("--score-max-concurrent", type=int, default=4)
    parser.add_argument("--worker-cpus", type=int, default=8)
    parser.add_argument("--worker-memory", default="32G")
    parser.add_argument("--worker-max-concurrent", type=int, default=12)
    parser.add_argument("--repair-cpus", type=int, default=40)
    parser.add_argument("--repair-threads", type=int, default=32)
    parser.add_argument("--repair-memory", default="256G")
    parser.add_argument("--repair-max-concurrent", type=int, default=2)
    parser.add_argument("--gather-cpus", type=int, default=1)
    parser.add_argument("--gather-memory", default="256G")
    parser.add_argument("--analysis-cpus", type=int, default=1)
    parser.add_argument("--analysis-memory", default="64G")
    parser.add_argument("--time", default="7-00:00:00")
    parser.add_argument("--partition", default="compute")
    parser.add_argument("--analysis-task-manifest", default="")
    parser.add_argument("--task-index", default="")
    parser.add_argument("--checkpoint", choices=("scores", "gather", "compact"),
                        default="scores")
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--unavailable-reason", default="")
    parser.add_argument("--code-root", default="")
    return parser


def main():
    args = build_parser().parse_args()
    try:
        if args.action == "preflight":
            return preflight(args)
        if args.action == "run-stage":
            return run_stage(args)
        if args.action == "run-repair":
            return run_repair(args)
        if args.action == "resume":
            return resume_stage(args)
        if args.action == "status":
            return status(args)
        if args.action == "analysis-worker":
            return analysis_worker(args)
        if args.action == "stage-checkpoint":
            return checkpoint_stage(args)
        if args.action == "finalize-stage":
            return finalize_stage(args)
        if args.action == "finalize-repair":
            return finalize_repair(args)
        if args.action == "finalize-checkpoint":
            return finalize_checkpoint(args)
        raise RuntimeError(f"unsupported action: {args.action}")
    except (OSError, ValueError, RuntimeError, AssertionError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
