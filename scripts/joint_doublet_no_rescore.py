#!/usr/bin/env python3
"""No-rescore completion analysis and targeted-workload renderer.

This program consumes only completed joint-doublet ledgers, manifests, and
score tables.  It never prepares cells, reads BAMs/fragments during the
no-rescore action, or invokes tetra_score_calls.  It also renders (but does not
submit) the later cache-backed targeted validation workload.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import heapq
import json
import math
import os
import random
import re
import resource
import shlex
import shutil
import statistics
import subprocess
import sys
import time
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


SCIENTIFIC_METHOD_VERSION = "JOINT_DOUBLET_TARGETED_V4_20260920"
SCHEMA = "joint_doublet_no_rescore_completion_v4"
CACHE_SCHEMA = "joint_doublet_normalized_indexed_cache_v5"
WORKLOAD_SCHEMA = "joint_doublet_targeted_workload_v5"
TASK_MARKER_SCHEMA = "joint_doublet_artifact_marker_v5"
IMPLEMENTATION_VERSION = "targeted_bounded_20260926_v5"
SEED = 1729
ALLOWED_LIBRARIES = (7, 9, 12, 17, 20, 25, 29)
PRIMARY_LIBRARIES = (7, 9, 12, 17, 20, 29)
CALIBRATION_LIBRARY = 25
PROTECTED_LIBRARIES = (19, 35, 38)
TARGET_LIBRARIES = (12, 20, 29)
EARLY_EXCESS_LIBRARIES = (7, 12, 20)
LOW_EXCESS_LIBRARIES = (9, 17, 29)
INPUT_ROLES = (
    "cell_ledger", "rna_manifest", "atac_manifest", "rna_scores",
    "atac_scores",
)
INPUT_BASENAMES = {
    "cell_ledger": "{library}.cell_ledger.tsv.gz",
    "rna_manifest": "{library}.rna_joint_manifest.tsv.gz",
    "atac_manifest": "{library}.atac_joint_manifest.tsv.gz",
    "rna_scores": "{library}.rna_joint_scores.tsv.gz",
    "atac_scores": "{library}.atac_joint_scores.tsv.gz",
}
MANIFEST_MEANING_FIELDS = (
    "candidate_policy", "physical_pool_state", "component_only_state",
    "structural_added_state_relationship", "candidate_origin",
    "exhaustive_fallback", "nomination_modalities",
)
SITE_METRICS = {
    "status": "score_status",
    "percentile": "library25_empirical_percentile",
    "locked": "k1_site_balanced_log_likelihood",
    "interior": "k2_site_balanced_log_likelihood",
    "contributor": "contributor_only_fraction1_site_balanced_log_likelihood",
    "preferred_model": "preferred_site_model_interpretation",
    "delta": "delta_site_balanced_log_likelihood_k2_minus_k1",
    "fraction": "fitted_second_fraction",
    "profile_low": "fitted_second_fraction_profile_low",
    "profile_high": "fitted_second_fraction_profile_high",
    "fold_support": "leave_one_fold_out_support_fraction",
    "fold_minimum": "minimum_leave_one_fold_out_balanced_delta",
    "influence": "maximum_single_site_absolute_balanced_delta_fraction",
    "units": "n_discriminating_sites",
    "effective_units": "n_discriminating_sites",
    "warnings": "warnings",
}
MOLECULE_METRICS = {
    "status": "molecule_score_status",
    "percentile": "molecule_library25_empirical_percentile",
    "locked": "molecule_balanced_k1_log_likelihood",
    "interior": "molecule_balanced_k2_log_likelihood",
    "contributor":
        "molecule_balanced_contributor_only_fraction1_log_likelihood",
    "preferred_model": "preferred_molecule_model_interpretation",
    "delta": "molecule_balanced_delta_log_likelihood_k2_minus_k1",
    "fraction": "molecule_balanced_fitted_second_fraction",
    "profile_low": "molecule_balanced_fitted_second_fraction_profile_low",
    "profile_high": "molecule_balanced_fitted_second_fraction_profile_high",
    "fold_support": "molecule_heldout_fold_support_fraction",
    "fold_minimum": "molecule_heldout_minimum_fold_mean_delta",
    "influence": "maximum_single_linked_unit_absolute_contribution_fraction",
    "units": "n_discriminating_linked_units",
    "effective_units": "effective_linked_unit_count",
    "warnings": "molecule_warnings",
}
CONTROL_FRACTIONS = (0.10, 0.20, 0.35)
DOWNSAMPLE_FRACTIONS = (0.25, 0.50, 0.75)
DOWNSAMPLE_REPLICATES = 100
CELL_NULL_REPLICATES = 1000
BOOTSTRAP_REPLICATES = 2000
MIN_CALIBRATION_REFERENCE = 20
TARGET_MATCH_MAX_ABS_Z = 0.75
TARGET_MATCH_MAX_DISTANCE = 1.50
DECOY_MIN_SHARED_SITES = 20
DECOY_RELATIVE_TOLERANCE = 0.10
CONTROL_SOURCE_CELL_REUSE_CAP = 3
CONTROL_SOURCE_DONOR_REUSE_CAP = 20
CONTROL_RECIPIENT_REUSE_CAP = 3
FROZEN_TARGET_BASENAME = "joint_doublet_frozen_targets_20260920.tsv"
STRICT_SITE_TARGETS = {
    "lib12:GCCTAATAGCATTTCT", "lib12:TGTTACTTCAAGGACA",
    "lib12:TCCATTGTCTTTGACT", "lib20:CCTGAGTCATGTTGCA",
    "lib29:TCAATCGCAATAACCT",
}
STRICT_SITE_TARGET_CONTRIBUTORS = {
    "lib12:GCCTAATAGCATTTCT": ("C40280", "JOS3C1"),
    "lib12:TGTTACTTCAAGGACA": ("H25576", "JOS3C1"),
    "lib12:TCCATTGTCTTTGACT": ("CongoA4B", "C3624"),
    "lib20:CCTGAGTCATGTTGCA": ("H20961", "JOS3C1"),
    "lib29:TCAATCGCAATAACCT": ("H1", "KOLF"),
}
LEGACY_EXACT_MOLECULE_TARGETS = {
    "lib12:ATTCATGAGGACAATG", "lib12:GCCTAATAGCATTTCT",
    "lib12:TGTTACTTCAAGGACA", "lib20:ACACTAGGTGTTGTAG",
    "lib20:AGCGTGCTCTGTGCAG", "lib20:CCTGAGTCATGTTGCA",
    "lib20:TTACGTTTCGCGCTAA",
}
LEGACY_EXACT_MOLECULE_CONTRIBUTORS = {
    "lib12:ATTCATGAGGACAATG": ("KOLF", "C3624"),
    "lib12:GCCTAATAGCATTTCT": ("C40280", "JOS3C1"),
    "lib12:TGTTACTTCAAGGACA": ("H25576", "JOS3C1"),
    "lib20:ACACTAGGTGTTGTAG": ("H21792", "C3624"),
    "lib20:AGCGTGCTCTGTGCAG": ("H25576", "C8861"),
    "lib20:CCTGAGTCATGTTGCA": ("H20961", "JOS3C1"),
    "lib20:TTACGTTTCGCGCTAA": ("H25576", "JOS3C1"),
}
HISTORICAL_SITE_DUAL_P95_CONTRIBUTORS = {
    "lib12:ATTCATGAGGACAATG": "C3624",
    "lib12:CTTCAATTCCTTGTTG": "C3624",
    "lib12:GCCCAAATCCCTGGAA": "C3651",
    "lib12:GCCTAATAGCATTTCT": "JOS3C1",
    "lib20:AGCGTGCTCTGTGCAG": "C8861",
    "lib20:GCTTATCGTCCGCTGT": "C3624",
    "lib20:TGTGGCGGTGGACCTG": "C3651",
    "lib20:TGTTGCACAAGATTCT": "C40280",
    "lib29:TCAATCGCAATAACCT": "KOLF",
}
HISTORICAL_MOLECULE_RAW_DUAL_P95 = set(
    HISTORICAL_SITE_DUAL_P95_CONTRIBUTORS) | {"lib12:TGCACACCAATGAGGT"}
ENDPOINT_LABELS = {
    "SITE_UNRESTRICTED": "all legal proposals",
    "SITE_GENETIC_SOURCE": "genetically distinguishable new contributor",
    "SITE_STRICT_MIXTURE": "strict two-source mixture",
    "SITE_ADDITION_UNCERTAIN": "addition-compatible but weakly separated",
    "SITE_REPLACEMENT_BOUNDARY": "replacement or boundary fit",
    "MOLECULE_GENETIC_SOURCE": "genetically distinguishable new contributor",
}
PILOT_SCORER_HOURS = {
    "lib7": {"RNA": 1.72, "ATAC": 1.97},
    "lib9": {"RNA": 1.44, "ATAC": 2.08},
    "lib12": {"RNA": 0.60, "ATAC": 1.36},
    "lib17": {"RNA": 1.95, "ATAC": 2.75},
    "lib20": {"RNA": 0.65, "ATAC": 1.38},
    "lib25": {"RNA": 3.79, "ATAC": 2.34},
    "lib29": {"RNA": 1.89, "ATAC": 4.22},
}


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def clean(value):
    text = "" if value is None else str(value).strip()
    return "" if text.lower() in {"na", "nan", "none", "null"} else text


def finite(value, default=math.nan):
    try:
        result = float(clean(value))
        return result if math.isfinite(result) else default
    except (TypeError, ValueError):
        return default


def truthy(value):
    return clean(value).lower() in {"1", "true", "yes", "y"}


def library_name(value):
    digits = "".join(character for character in str(value) if character.isdigit())
    if not digits:
        raise ValueError(f"invalid library identifier: {value}")
    return f"lib{int(digits)}"


def stable_seed(*parts):
    payload = "\x1f".join(str(part) for part in (SEED,) + parts)
    return int.from_bytes(hashlib.sha256(payload.encode()).digest()[:8], "big")


class PhaseAccounting:
    def __init__(self):
        self.started = time.monotonic()
        self.phase_started = self.started
        self.rows = []

    def finish(self, phase, **counts):
        now = time.monotonic()
        row = {
            "phase": phase,
            "phase_seconds": now - self.phase_started,
            "wall_seconds": now - self.phase_started,
            "elapsed_seconds": now - self.started,
            **counts,
        }
        self.rows.append(row)
        self.phase_started = now
        return row


def load_frozen_targets(path=""):
    if not path:
        raise RuntimeError(
            "--frozen-targets must name the authoritative frozen-target manifest")
    authoritative = Path(path).resolve()
    if not authoritative.is_file():
        raise RuntimeError(
            f"authoritative frozen-target manifest is missing: {authoritative}")

    rows = list(read_tsv(authoritative))
    if not rows or not {"target_id", "library"}.issubset(rows[0]) or not (
            "barcode" in rows[0] or "target_barcode" in rows[0]):
        raise RuntimeError("frozen manifest lacks target identity columns")
    parsed = []
    seen = set()
    for row in rows:
        library = library_name(row["library"])
        barcode = clean(row.get("target_barcode", row.get("barcode", "")))
        key = (library, barcode)
        if int(library[3:]) not in ALLOWED_LIBRARIES or not re.fullmatch(
                r"[ACGT]{16}", barcode) or row["target_id"] != f"{library}:{barcode}":
            raise RuntimeError(f"invalid frozen target identity: {row['target_id']}")
        if key in seen:
            raise RuntimeError(f"duplicate frozen target: {row['target_id']}")
        seen.add(key)
        parsed.append({"target_id": row["target_id"],
                       "library": int(library[3:]), "barcode": barcode})
    return authoritative, parsed


def load_frozen_comparisons(path, frozen_keys):
    """Read the authoritative selection, including its original comparison map."""
    if not path or not Path(path).is_absolute():
        raise RuntimeError("--frozen-comparisons must be the absolute path to "
                           "target_cells_and_matched_comparisons.tsv")
    rows = list(read_tsv(path))
    required = {"target_id", "library", "target_barcode", "locked_identity",
                "proposed_contributor", "matched_comparison_barcode",
                "selection_rule"}
    if not rows or not required.issubset(rows[0]):
        raise RuntimeError("frozen comparison manifest lacks selection/provenance columns")
    targets, comparisons = set(), set()
    for row in rows:
        library = library_name(row["library"])
        row["library"] = library
        key = (library, clean(row["target_barcode"]))
        comparison = (library, clean(row["matched_comparison_barcode"]))
        if row["target_id"] != ":".join(key) or key in targets or not re.fullmatch(
                r"[ACGT]{16}", comparison[1]) or comparison in comparisons:
            raise RuntimeError("duplicate/invalid frozen target or comparison identity")
        if not clean(row["locked_identity"]) or not clean(row["proposed_contributor"]):
            raise RuntimeError("frozen selection lacks locked identity/contributor provenance")
        targets.add(key)
        comparisons.add(comparison)
    if targets != frozen_keys or targets & comparisons:
        raise RuntimeError("frozen comparison identities differ from frozen target selection")
    return rows


def reuse_frozen_comparisons(cells, frozen_rows):
    by_key = {(row["library"], row["barcode"]): row for row in cells}
    targets, comparisons, manifest, balance = [], [], [], []
    for frozen in frozen_rows:
        row = dict(frozen)  # historical fields are immutable provenance
        for role, field, output in (("target", "target_barcode", targets),
                                    ("comparison", "matched_comparison_barcode", comparisons)):
            key = (row["library"], row[field])
            current = by_key.get(key)
            if current is None:
                # Retain the identity. A missing complete menu is recorded by the
                # planner as unavailable, never replaced by another cell.
                current = {"library": key[0], "barcode": key[1], "menu_size": 0,
                           "_candidates": [], "evidence_status": "UNAVAILABLE_RETAINED_CELL"}
                cells.append(current)
                by_key[key] = current
            if role == "target" and clean(current.get("reconciled_identity_locked")) and \
                    current["reconciled_identity_locked"] != row["locked_identity"]:
                raise RuntimeError(f"locked identity changed for frozen target {row['target_id']}")
            row[f"current_{role}_evidence_status"] = current.get(
                "evidence_status", "RETAINED_SCORES_AVAILABLE")
            row[f"current_{role}_rna_winner"] = current.get("rna_site_genetic_second_state", "")
            row[f"current_{role}_atac_winner"] = current.get("atac_site_genetic_second_state", "")
            output.append(current)
        row["comparison_mapping_source"] = "AUTHORITATIVE_FROZEN_MANIFEST"
        manifest.append(row)
    for library in ["ALL_TARGET"] + sorted({row["library"] for row in manifest}):
        pairs = [(left, right) for left, right in zip(targets, comparisons)
                 if library == "ALL_TARGET" or left["library"] == library]
        for index, field in enumerate(COVERAGE_COVARIATES):
            values = [(coverage_vector(left)[index], coverage_vector(right)[index]) for left, right in pairs]
            complete = [(a,b) for a,b in values if math.isfinite(a) and math.isfinite(b)]
            difference = standardized_difference([a for a,b in complete], [b for a,b in complete])
            balance.append({"balance_stage": "AFTER_MATCHING", "library": library, "covariate": field,
                            "requested_pairs": len(pairs), "complete_pairs": len(complete),
                            "unavailable_pairs": len(pairs)-len(complete),
                            "absolute_standardized_difference": abs(difference),
                            "balance_status": "PASS" if math.isfinite(difference) and abs(difference)<=0.1 else "FAIL",
                            "interpretation": "DESCRIPTIVE_ORIGINAL_FROZEN_PAIRS_NO_RESELECTION"})
    return targets, comparisons, manifest, balance


def endpoint_label(endpoint):
    return ENDPOINT_LABELS.get(endpoint, endpoint.replace("_", " ").lower())


def target_id(row):
    return f"{row['library']}:{row['barcode']}"


def open_text(path, mode="rt"):
    return gzip.open(path, mode, encoding="utf-8", newline="") \
        if str(path).endswith(".gz") else \
        open(path, mode, encoding="utf-8", newline="")


def read_tsv(path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


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


def write_tsv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if fields is None:
        fields = []
        seen = set()
        for row in rows:
            for field in row:
                if field not in seen:
                    fields.append(field)
                    seen.add(field)
    fields = list(fields)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(temporary, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore",
            lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def quantile(values, probability):
    ordered = sorted(value for value in values if math.isfinite(value))
    if not ordered:
        return math.nan
    position = probability * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def percentile(value, reference):
    if not math.isfinite(value) or not reference:
        return math.nan
    return bisect.bisect_right(reference, value) / len(reference)


def median(values):
    values = [value for value in values if math.isfinite(value)]
    return statistics.median(values) if values else math.nan


def mean(values):
    values = [value for value in values if math.isfinite(value)]
    return statistics.mean(values) if values else math.nan


def donor_components(state):
    state = clean(state)
    if not state or state.startswith("M{"):
        return []
    return [token for token in state.split("+") if token]


def component_shape(state):
    components = donor_components(state)
    return len(components), len(set(components))


def gzip_envelope_valid(path):
    try:
        path = Path(path)
        if path.stat().st_size < 18:
            return False
        with open(path, "rb") as handle:
            header = handle.read(2)
            handle.seek(-8, os.SEEK_END)
            trailer = handle.read(8)
        return header == b"\x1f\x8b" and len(trailer) == 8
    except OSError:
        return False


def read_header(path):
    try:
        with open_text(path) as handle:
            return handle.readline().rstrip("\r\n").split("\t")
    except (OSError, EOFError, gzip.BadGzipFile):
        return []


def _path_mentions_protected(path):
    lowered = str(Path(path)).lower()
    return any(
        re.search(rf"(?<![a-z0-9])lib0*{library}(?![0-9])", lowered) or
        re.search(rf"(?<![a-z0-9])library0*{library}(?![0-9])", lowered) or
        re.search(
            rf"tet_2025_[^/]*[_-]0*{library}(?![0-9])", lowered)
        for library in PROTECTED_LIBRARIES)


def _is_within(path, root):
    try:
        Path(path).resolve(strict=False).relative_to(Path(root).resolve())
        return True
    except ValueError:
        return False


def safe_find_named(root, basename, allowed_library, access_audit=None):
    root = Path(root)
    if not root.is_dir():
        return []
    matches = []
    protected_names = {f"lib{library}" for library in PROTECTED_LIBRARIES}
    for current, directories, files in os.walk(root):
        current_path = Path(current)
        directories[:] = [
            name for name in directories
            if name.lower() not in protected_names and
            not re_library_is_protected(name) and
            not _path_mentions_protected(name) and
            not (name.lower().startswith("lib") and
                 "".join(character for character in name if character.isdigit()) and
                 int("".join(character for character in name
                             if character.isdigit())) not in
                 {int(allowed_library.removeprefix("lib"))})]
        if access_audit is not None:
            access_audit.append({
                "event": "DIRECTORY_ENUMERATION",
                "library": allowed_library,
                "path": str(current_path.resolve()),
                "protected_content_opened": False,
                "detail": "library-scoped traversal; protected and other-library directories pruned",
            })
        if basename in files:
            candidate = current_path / basename
            if not _path_mentions_protected(candidate) and \
                    (allowed_library in candidate.name.lower() or
                     allowed_library in {part.lower() for part in candidate.parts}):
                matches.append(candidate.resolve())
    return sorted(set(matches), key=str)


def re_library_is_protected(value):
    lowered = str(value).lower()
    digits = "".join(character for character in lowered if character.isdigit())
    return lowered.startswith("lib") and digits and int(digits) in PROTECTED_LIBRARIES


def load_task_rows(existing_stage_root, expected_libraries=None):
    """Load and validate the complete upstream task manifest.

    This is deliberately performed before discovery, stat, gzip probing, or
    content reads of any biological source.  The task manifest is the allow
    list; discovery is never allowed to invent a source contract first and
    validate it afterwards.
    """
    candidates = (
        Path(existing_stage_root) / "manifests" / "derivative_tasks_v2.tsv",
        Path(existing_stage_root) / "joint_doublet_tasks.tsv",
    )
    for path in candidates:
        if path.is_file():
            raw_rows = list(read_tsv(path))
            if not raw_rows:
                raise RuntimeError(f"upstream task manifest is empty: {path}")
            required = {
                "library", "output_dir", "rna_manifest", "atac_manifest",
                "rna_scores", "atac_scores",
            }
            if not required.issubset(raw_rows[0]):
                raise RuntimeError(
                    "upstream task manifest schema is incomplete: " +
                    ",".join(sorted(required - set(raw_rows[0]))))
            expected = set(expected_libraries or ALLOWED_LIBRARIES)
            rows = {}
            seen_paths = set()
            indices = []
            for row in raw_rows:
                library = library_name(row.get("library", ""))
                number = int(library.removeprefix("lib"))
                if number not in expected:
                    # Rows outside the explicitly authorized no-rescore set are
                    # rejected, rather than silently filtered from an allowlist.
                    raise RuntimeError(
                        f"upstream task manifest contains unauthorized {library}")
                if library in rows:
                    raise RuntimeError(
                        f"upstream task manifest duplicates {library}")
                if clean(row.get("task_index", "")):
                    try:
                        indices.append(int(row["task_index"]))
                    except ValueError as error:
                        raise RuntimeError(
                            "upstream task manifest has a nonnumeric task index") \
                            from error
                for field in ("output_dir", "rna_manifest", "atac_manifest",
                              "rna_scores", "atac_scores"):
                    value = clean(row.get(field, ""))
                    if not value or not Path(value).is_absolute() or \
                            _path_mentions_protected(value):
                        raise RuntimeError(
                            f"invalid upstream task path: {library}:{field}")
                    normalized = os.path.abspath(value)
                    key = (field, normalized)
                    if key in seen_paths:
                        raise RuntimeError(
                            f"upstream task manifest duplicates {field} path")
                    seen_paths.add(key)
                rows[library] = dict(row)
            if set(rows) != {f"lib{value}" for value in expected}:
                raise RuntimeError(
                    "upstream task manifest does not contain exactly the "
                    "authorized no-rescore libraries")
            if indices and (len(indices) != len(raw_rows) or
                            sorted(indices) != list(range(len(raw_rows))) or
                            len(indices) != len(set(indices))):
                raise RuntimeError(
                    "upstream task indices must be present, unique, and contiguous")
            return path.resolve(), rows
    raise RuntimeError(
        "a complete upstream joint-doublet task manifest is required before "
        "source discovery")


def validate_no_rescore_candidate_manifests(task_rows, libraries):
    """Validate both assay candidate menus in full before other source opens."""
    required = {
        "library", "barcode", "candidate_id", "locked_state",
        "locked_copy_vector", "second_state", "second_copy_vector",
        *MANIFEST_MEANING_FIELDS,
    }
    validated = {}
    for number in libraries:
        library = f"lib{number}"
        row = task_rows[library]
        assay_keys = {}
        assay_semantics = {}
        for modality, field in (("RNA", "rna_manifest"),
                                ("ATAC", "atac_manifest")):
            path = Path(row[field])
            # Opening these two files is itself the manifest-validation step;
            # no ledger, score table, pileup, or cache has yet been discovered
            # or statted.
            try:
                manifest_rows = list(read_tsv(path))
            except (OSError, EOFError, gzip.BadGzipFile,
                    UnicodeDecodeError) as error:
                raise RuntimeError(
                    f"candidate manifest cannot be fully read: {path}: {error}") \
                    from error
            if not manifest_rows or not required.issubset(manifest_rows[0]):
                missing = required - set(manifest_rows[0] if manifest_rows else {})
                raise RuntimeError(
                    f"candidate manifest schema invalid: {path}:" +
                    ",".join(sorted(missing)))
            keys = []
            semantics = {}
            cells = defaultdict(list)
            for candidate in manifest_rows:
                if library_name(candidate.get("library", "")) != library:
                    raise RuntimeError(
                        f"candidate manifest mixes libraries: {path}")
                barcode = clean(candidate.get("barcode", ""))
                candidate_id = clean(candidate.get("candidate_id", ""))
                if not barcode or not candidate_id:
                    raise RuntimeError(
                        f"candidate manifest has an empty key: {path}")
                if clean(candidate.get("candidate_policy", "")) != \
                        "DERIVATIVE_COMPLETE":
                    raise RuntimeError(
                        f"candidate manifest is not derivative-complete: {path}")
                key = (library, barcode, candidate_id)
                keys.append(key)
                semantic = tuple(clean(candidate.get(field_name, ""))
                                 for field_name in (
                                     "locked_state", "locked_copy_vector",
                                     "second_state", "second_copy_vector",
                                     *MANIFEST_MEANING_FIELDS))
                semantics[key] = semantic
                cells[barcode].append(candidate_id)
            if len(keys) != len(set(keys)):
                raise RuntimeError(
                    f"candidate manifest contains duplicate keys: {path}")
            if any(not menu for menu in cells.values()):
                raise RuntimeError(
                    f"candidate manifest contains an empty menu: {path}")
            assay_keys[modality] = set(keys)
            assay_semantics[modality] = semantics
            validated[(library, modality)] = {
                "path": str(path), "rows": len(keys), "cells": len(cells),
                "content_sha256": _small_file_sha256(path),
            }
        if assay_keys["RNA"] != assay_keys["ATAC"] or \
                assay_semantics["RNA"] != assay_semantics["ATAC"]:
            raise RuntimeError(
                f"RNA/ATAC candidate menu mismatch for {library}")
    return validated


def discover_inputs(existing_stage_root, task_search_root, libraries,
                    task_manifest, task_rows):
    existing_stage_root = Path(existing_stage_root).resolve()
    task_search_root = Path(task_search_root).resolve()
    sources = []
    inventory = []
    warnings = []
    access_audit = []
    required_headers = {
        "cell_ledger": {
            "library", "barcode", "reconciled_identity_locked",
            "biological_ploidy", "biological_state",
            "reconciled_droplet_state", "identity_disposition",
            "demux_modality_relation", "rna_total_counts",
            "rna_detected_features", "atac_fragments",
            "atac_fragment_records", "atac_cut_sites",
        },
        "rna_manifest": {
            "library", "barcode", "candidate_id", "locked_state",
            "locked_copy_vector", "second_state", "second_copy_vector",
            *MANIFEST_MEANING_FIELDS,
        },
        "atac_manifest": {
            "library", "barcode", "candidate_id", "locked_state",
            "locked_copy_vector", "second_state", "second_copy_vector",
            *MANIFEST_MEANING_FIELDS,
        },
        "rna_scores": {
            "library", "barcode", "candidate_id", "score_status",
            "delta_site_balanced_log_likelihood_k2_minus_k1",
            "fitted_second_fraction", "fitted_second_fraction_profile_low",
            "fitted_second_fraction_profile_high",
            "leave_one_fold_out_support_fraction",
            "minimum_leave_one_fold_out_balanced_delta",
            "maximum_single_site_absolute_balanced_delta_fraction",
            "n_discriminating_sites", "discriminating_depth",
            "replacement_like", "molecule_score_status",
            "molecule_balanced_delta_log_likelihood_k2_minus_k1",
            "molecule_balanced_fitted_second_fraction",
            "molecule_balanced_fitted_second_fraction_profile_low",
            "molecule_balanced_fitted_second_fraction_profile_high",
            "molecule_heldout_fold_support_fraction",
            "maximum_single_linked_unit_absolute_contribution_fraction",
            "n_discriminating_linked_units", "effective_linked_unit_count",
        },
        "atac_scores": {
            "library", "barcode", "candidate_id", "score_status",
            "delta_site_balanced_log_likelihood_k2_minus_k1",
            "fitted_second_fraction", "fitted_second_fraction_profile_low",
            "fitted_second_fraction_profile_high",
            "leave_one_fold_out_support_fraction",
            "minimum_leave_one_fold_out_balanced_delta",
            "maximum_single_site_absolute_balanced_delta_fraction",
            "n_discriminating_sites", "discriminating_depth",
            "replacement_like", "molecule_score_status",
            "molecule_balanced_delta_log_likelihood_k2_minus_k1",
            "molecule_balanced_fitted_second_fraction",
            "molecule_balanced_fitted_second_fraction_profile_low",
            "molecule_balanced_fitted_second_fraction_profile_high",
            "molecule_heldout_fold_support_fraction",
            "maximum_single_linked_unit_absolute_contribution_fraction",
            "n_discriminating_linked_units", "effective_linked_unit_count",
        },
    }
    for number in libraries:
        library = f"lib{number}"
        source = {"library": library}
        task_row = task_rows.get(library, {})
        for role in INPUT_ROLES:
            basename = INPUT_BASENAMES[role].format(library=library)
            direct = []
            task_value = clean(task_row.get(role, ""))
            manifest_declared = bool(task_value)
            if task_value:
                declared = Path(task_value)
                if not declared.is_absolute():
                    raise RuntimeError(
                        f"input manifest path is not absolute: {library}:{role}")
                resolved_declared = declared.resolve(strict=False)
                if _path_mentions_protected(declared) or \
                        _path_mentions_protected(resolved_declared):
                    raise RuntimeError(
                        f"input manifest names protected content: {library}:{role}")
                if not any(_is_within(resolved_declared, allowed_root)
                           for allowed_root in (existing_stage_root,
                                                task_search_root)):
                    raise RuntimeError(
                        f"input manifest path escapes authorized roots: "
                        f"{library}:{role}:{declared}")
                direct.append(declared)
            else:
                direct.extend((
                    existing_stage_root / library / basename,
                    existing_stage_root / basename,
                    task_search_root / library / basename,
                    task_search_root / basename,
                ))
            usable = []
            seen_paths = set()
            for path in direct:
                resolved = path.resolve(strict=False)
                if resolved in seen_paths or _path_mentions_protected(resolved):
                    continue
                seen_paths.add(resolved)
                if path.is_file() and path.stat().st_size > 0:
                    usable.append(path.resolve())
            if not usable and not manifest_declared:
                usable = safe_find_named(
                    existing_stage_root, basename, library, access_audit)
            if not usable and not manifest_declared:
                usable = safe_find_named(
                    task_search_root, basename, library, access_audit)
            usable = sorted(set(usable), key=str)
            gzip_usable = [path for path in usable if gzip_envelope_valid(path)]
            chosen = (gzip_usable[0] if gzip_usable else usable[0]) \
                if usable else None
            if len(usable) > 1:
                warnings.append({
                    "library": library, "warning": "AMBIGUOUS_INPUT_MATCH",
                    "detail": f"{role}:" + "|".join(str(path) for path in usable),
                })
            gzip_ok = bool(chosen and gzip_envelope_valid(chosen))
            source[role] = str(chosen) if gzip_ok else ""
            header = read_header(chosen) if chosen else []
            if chosen:
                access_audit.append({
                    "event": "INPUT_OPEN",
                    "library": library,
                    "path": str(chosen),
                    "protected_content_opened": False,
                    "detail": role,
                })
            missing_fields = sorted(required_headers[role] - set(header))
            inventory.append({
                "record_type": "SOURCE_FILE", "library": library,
                "assay": "RNA" if role.startswith("rna_") else
                    "ATAC" if role.startswith("atac_") else "BOTH",
                "input_role": role,
                "absolute_path": str(chosen) if chosen else "UNRESOLVED",
                "exists": bool(chosen),
                "bytes": chosen.stat().st_size if chosen else 0,
                "gzip_envelope_status": (
                    "PASS" if gzip_ok else
                    "FAIL" if chosen else "MISSING"),
                "schema_status": (
                    "PASS" if chosen and not missing_fields else
                    "MISSING_FIELDS:" + ",".join(missing_fields)
                    if chosen else "MISSING"),
                "unique_key_status": "CHECKED_DURING_AGGREGATION",
                "count": "",
                "exclusion_reason": "" if chosen else "MISSING_EXPECTED_INPUT",
            })
            if not chosen:
                warnings.append({
                    "library": library, "warning": "MISSING_EXPECTED_INPUT",
                    "detail": f"{role}:{basename}",
                })
            elif not gzip_ok:
                warnings.append({
                    "library": library, "warning": "INVALID_GZIP_ENVELOPE",
                    "detail": f"{role}:{chosen}",
                })
            if chosen and missing_fields:
                warnings.append({
                    "library": library, "warning": "MISSING_REQUIRED_FIELDS",
                    "detail": f"{role}:" + ",".join(missing_fields),
                })
        sources.append(source)
    accessed = sorted({
        int(row["library"].removeprefix("lib"))
        for row in access_audit if row["event"] == "INPUT_OPEN"
    })
    if set(accessed) - set(libraries):
        raise RuntimeError(
            f"input-open audit contains unexpected libraries: {accessed}")
    return (sources, inventory, warnings, task_manifest, task_rows,
            access_audit, accessed)


def canonical_menu_signature(candidates):
    pairs = sorted({
        (candidate["second_state"], candidate["second_copy_vector"])
        for candidate in candidates
    })
    return ";".join(f"{state}|{copy}" for state, copy in pairs)


def evaluable_mask(candidates, modality, evidence="site"):
    return ";".join(sorted(
        f"{candidate['second_state']}|{candidate['second_copy_vector']}"
        for candidate in candidates
        if genetic_candidate(candidate, candidate[f"{modality}_{evidence}"])))


def _anonymous_copy_multiplicities(copy_vector):
    values = []
    for token in clean(copy_vector).split(","):
        if not token:
            continue
        try:
            _donor, copies = token.rsplit(":", 1)
            values.append(float(copies))
        except (ValueError, TypeError):
            return ("UNAVAILABLE",)
    return tuple(sorted(values))


def _anonymous_candidate_shape(candidate):
    locked_components, locked_distinct = component_shape(
        candidate.get("locked_state", ""))
    candidate_components, candidate_distinct = component_shape(
        candidate.get("second_state", ""))
    return (
        locked_components, locked_distinct,
        _anonymous_copy_multiplicities(candidate.get("locked_copy_vector", "")),
        candidate_components, candidate_distinct,
        _anonymous_copy_multiplicities(candidate.get("second_copy_vector", "")),
        candidate.get("relationship", "UNAVAILABLE"),
        candidate.get("candidate_origin", "UNAVAILABLE"),
        int(candidate.get("physical", False)),
        int(candidate.get("component_only", False)),
        int(candidate.get("new_donor", False)),
        int(candidate.get("reweighting", False)),
        int(candidate.get("state_genotype_equivalent", False)),
    )


def transferable_geometry_signature(row):
    """Donor-name-invariant candidate geometry for external transfer checks."""
    locked_components, locked_distinct = component_shape(
        row.get("reconciled_identity_locked", ""))
    shapes = [_anonymous_candidate_shape(candidate)
              for candidate in row.get("_candidates", [])]
    return json.dumps({
        "locked_components": locked_components,
        "locked_distinct_donors": locked_distinct,
        "locked_copy_multiplicities": sorted({
            _anonymous_copy_multiplicities(candidate.get(
                "locked_copy_vector", "UNAVAILABLE"))
            for candidate in row.get("_candidates", [])}),
        "canonical_menu_size": int(row.get("menu_size", 0) or 0),
        "candidate_shapes": sorted(shapes),
    }, sort_keys=True, separators=(",", ":"))


def candidate_origin_class(candidate):
    if candidate["fallback"]:
        return "EXHAUSTIVE_FALLBACK"
    if candidate["physical"]:
        return "PHYSICAL_POOL_STATE"
    if candidate["component_only"]:
        return "COMPONENT_ONLY_NONPHYSICAL"
    return "OTHER_LEGAL_STATE"


def extract_metric(row, modality, evidence):
    fields = SITE_METRICS if evidence == "site" else MOLECULE_METRICS
    status = clean(row.get(f"{modality}_{fields['status']}", "")).upper()
    value = {
        "status": status or "MISSING",
        "percentile": finite(row.get(f"{modality}_{fields['percentile']}", "")),
        "locked": finite(row.get(f"{modality}_{fields['locked']}", "")),
        "interior": finite(row.get(f"{modality}_{fields['interior']}", "")),
        "contributor": finite(row.get(
            f"{modality}_{fields['contributor']}", "")),
        "preferred_model": clean(row.get(
            f"{modality}_{fields['preferred_model']}", "")),
        "delta": finite(row.get(f"{modality}_{fields['delta']}", "")),
        "fraction": finite(row.get(f"{modality}_{fields['fraction']}", "")),
        "profile_low": finite(row.get(f"{modality}_{fields['profile_low']}", "")),
        "profile_high": finite(row.get(f"{modality}_{fields['profile_high']}", "")),
        "fold_support": finite(row.get(f"{modality}_{fields['fold_support']}", "")),
        "fold_minimum": finite(row.get(f"{modality}_{fields['fold_minimum']}", "")),
        "influence": finite(row.get(f"{modality}_{fields['influence']}", "")),
        "units": finite(row.get(f"{modality}_{fields['units']}", "")),
        "effective_units": finite(
            row.get(f"{modality}_{fields['effective_units']}", "")),
        "warnings": clean(row.get(f"{modality}_{fields['warnings']}", "")),
        "genotype_equivalent": (
            status == "GENOTYPE_EQUIVALENT" or
            truthy(row.get(f"{modality}_genotype_equivalent", ""))),
        "discriminating_sites": finite(
            row.get(f"{modality}_n_discriminating_sites", "")),
        "discriminating_depth": finite(
            row.get(f"{modality}_discriminating_depth", "")),
        "replacement_like": truthy(
            row.get(f"{modality}_replacement_like", "")),
        "umi_gene_fraction": finite(
            row.get(f"{modality}_rna_umi_gene_basis_fraction", "")),
        "qname_fraction": finite(
            row.get(f"{modality}_query_name_fallback_basis_fraction", "")),
    }
    if not value["preferred_model"] and all(math.isfinite(value[field])
            for field in ("locked", "interior", "contributor")):
        value["preferred_model"] = "LOCKED_IDENTITY"
        preferred = value["locked"]
        if value["fraction"] > 1e-12 and value["interior"] > preferred:
            value["preferred_model"] = "INTERIOR_LOCKED_PLUS_CONTRIBUTOR"
            preferred = value["interior"]
        if value["contributor"] > preferred:
            value["preferred_model"] = "CONTRIBUTOR_ONLY_FRACTION_1"
    value["available"] = status == "AVAILABLE" and math.isfinite(value["delta"])
    return value


def compact_candidate(row):
    relationship = clean(row.get("structural_added_state_relationship", ""))
    candidate = {
        "candidate_id": clean(row.get("candidate_id", "")),
        "locked_state": clean(row.get("locked_state", "")),
        "locked_copy_vector": clean(row.get("locked_copy_vector", "")),
        "second_state": clean(row.get("second_state", "")),
        "second_copy_vector": clean(row.get("second_copy_vector", "")),
        "candidate_origin": clean(row.get("candidate_origin", "")),
        "candidate_policy": clean(row.get("candidate_policy", "")),
        "physical": truthy(row.get("physical_pool_state", "")),
        "component_only": truthy(row.get("component_only_state", "")),
        "relationship": relationship,
        "fallback": truthy(row.get("exhaustive_fallback", "")),
        "nomination_modalities": clean(row.get("nomination_modalities", "")),
        "genotype_equivalent": truthy(row.get("genotype_equivalent_candidate", "")),
    }
    candidate["state_genotype_equivalent"] = bool(
        candidate["locked_copy_vector"] and
        candidate["locked_copy_vector"] == candidate["second_copy_vector"])
    candidate["origin_class"] = candidate_origin_class(candidate)
    candidate["new_donor"] = "CONTAINS_NEW_DONOR" in relationship
    candidate["reweighting"] = "EXISTING_COMPONENT_REWEIGHTING" in relationship
    candidate["genetically_identical"] = \
        "GENETICALLY_IDENTICAL_NO_CHANGE" in relationship
    for modality in ("rna", "atac"):
        candidate[f"{modality}_site"] = extract_metric(row, modality, "site")
        candidate[f"{modality}_molecule"] = extract_metric(
            row, modality, "molecule")
    candidate["_cache_record_counts"] = tuple(int(max(finite(row.get(
        f"{modality}_{field}", ""), 0.0), 0.0))
        for modality in ("rna", "atac")
        for field in ("n_common_nuclear_sites", "total_snps_in_linked_units"))
    return candidate


def top_two(candidates, modality, evidence, predicate=None):
    values = []
    for candidate in candidates:
        metric = candidate[f"{modality}_{evidence}"]
        if not metric["available"] or not math.isfinite(metric["percentile"]):
            continue
        if predicate is not None and not predicate(candidate, metric):
            continue
        values.append((candidate, metric))
    values.sort(key=lambda item: (
        item[1]["percentile"], item[1]["delta"], item[0]["candidate_id"]),
        reverse=True)
    return values[:2]


def flatten_top(prefix, ranked):
    if not ranked:
        return {
            f"{prefix}_candidate_id": "", f"{prefix}_second_state": "",
            f"{prefix}_runner_up_candidate_id": "",
            f"{prefix}_runner_up_second_state": "",
            f"{prefix}_second_copy_vector": "", f"{prefix}_percentile": math.nan,
            f"{prefix}_runner_up_percentile": math.nan,
            f"{prefix}_margin": math.nan, f"{prefix}_delta": math.nan,
            f"{prefix}_fraction": math.nan, f"{prefix}_profile_low": math.nan,
            f"{prefix}_profile_high": math.nan, f"{prefix}_fold_support": math.nan,
            f"{prefix}_fold_minimum": math.nan, f"{prefix}_influence": math.nan,
            f"{prefix}_units": math.nan, f"{prefix}_effective_units": math.nan,
            f"{prefix}_warnings": "", f"{prefix}_physical_pool_state": "",
            f"{prefix}_component_only_state": "", f"{prefix}_relationship": "",
            f"{prefix}_origin_class": "", f"{prefix}_umi_gene_fraction": math.nan,
            f"{prefix}_qname_fraction": math.nan,
            f"{prefix}_preferred_model": "UNAVAILABLE",
            f"{prefix}_replacement_like": False,
            f"{prefix}_genotype_equivalent": False,
        }
    candidate, metric = ranked[0]
    runner = ranked[1][1]["percentile"] if len(ranked) > 1 else math.nan
    return {
        f"{prefix}_candidate_id": candidate["candidate_id"],
        f"{prefix}_second_state": candidate["second_state"],
        f"{prefix}_runner_up_candidate_id": (
            ranked[1][0]["candidate_id"] if len(ranked) > 1 else ""),
        f"{prefix}_runner_up_second_state": (
            ranked[1][0]["second_state"] if len(ranked) > 1 else ""),
        f"{prefix}_second_copy_vector": candidate["second_copy_vector"],
        f"{prefix}_percentile": metric["percentile"],
        f"{prefix}_runner_up_percentile": runner,
        f"{prefix}_margin": metric["percentile"] - runner
            if math.isfinite(runner) else math.nan,
        f"{prefix}_delta": metric["delta"],
        f"{prefix}_fraction": metric["fraction"],
        f"{prefix}_profile_low": metric["profile_low"],
        f"{prefix}_profile_high": metric["profile_high"],
        f"{prefix}_fold_support": metric["fold_support"],
        f"{prefix}_fold_minimum": metric["fold_minimum"],
        f"{prefix}_influence": metric["influence"],
        f"{prefix}_units": metric["units"],
        f"{prefix}_effective_units": metric["effective_units"],
        f"{prefix}_warnings": metric["warnings"],
        f"{prefix}_physical_pool_state": candidate["physical"],
        f"{prefix}_component_only_state": candidate["component_only"],
        f"{prefix}_relationship": candidate["relationship"],
        f"{prefix}_origin_class": candidate["origin_class"],
        f"{prefix}_umi_gene_fraction": metric["umi_gene_fraction"],
        f"{prefix}_qname_fraction": metric["qname_fraction"],
        f"{prefix}_preferred_model": metric["preferred_model"],
        f"{prefix}_replacement_like": metric["replacement_like"],
        f"{prefix}_genotype_equivalent": (
            candidate["state_genotype_equivalent"] or
            metric["genotype_equivalent"]),
    }


def genetic_candidate(candidate, metric):
    return (candidate["physical"] and candidate["new_donor"] and
            not candidate["state_genotype_equivalent"] and
            not metric["genotype_equivalent"] and
            math.isfinite(metric["discriminating_sites"]) and
            metric["discriminating_sites"] > 0 and metric["available"])


def component_new_donor_candidate(candidate, metric):
    return (candidate["component_only"] and candidate["new_donor"] and
            not candidate["state_genotype_equivalent"] and
            not metric["genotype_equivalent"] and
            math.isfinite(metric["discriminating_sites"]) and
            metric["discriminating_sites"] > 0 and metric["available"])


def evidence_category(row, prefix):
    rna_state = clean(row.get(f"rna_{prefix}_second_state", ""))
    atac_state = clean(row.get(f"atac_{prefix}_second_state", ""))
    if not rna_state or not atac_state:
        evidence = "molecule" if prefix.startswith("molecule") else "site"
        statuses = {
            candidate[f"{modality}_{evidence}"].get("status", "")
            for candidate in row.get("_candidates", [])
            for modality in ("rna", "atac")
        }
        return "low/limited evidence" if statuses & {
            "LOW_EVIDENCE", "LIMITED_EVIDENCE"} else "unavailable"
    if rna_state != atac_state:
        return "conflicting assay/evidence result"
    if not truthy(row.get(f"rna_{prefix}_physical_pool_state", "")) or \
            not truthy(row.get(f"atac_{prefix}_physical_pool_state", "")) or \
            "CONTAINS_NEW_DONOR" not in clean(
                row.get(f"rna_{prefix}_relationship", "")):
        return "unavailable"
    preferred_models = {
        modality: clean(row.get(f"{modality}_{prefix}_preferred_model", ""))
        for modality in ("rna", "atac")
    }
    if "CONTRIBUTOR_ONLY_FRACTION_1" in preferred_models.values():
        return "replacement-like or upper-boundary fit"
    if preferred_models["rna"] != preferred_models["atac"]:
        return "conflicting assay/evidence result"
    if preferred_models["rna"] == "LOCKED_IDENTITY":
        return "locked-source preferred"
    values = {}
    for modality in ("rna", "atac"):
        values[modality] = {
            field: finite(row.get(f"{modality}_{prefix}_{field}", ""))
            for field in (
                "fraction", "profile_low", "profile_high", "margin",
                "fold_support", "influence")
        }
    boundary = any(
        truthy(row.get(f"{modality}_{prefix}_replacement_like", "")) or
        values[modality]["fraction"] >= 0.90 or
        values[modality]["profile_high"] >= 0.94
        for modality in ("rna", "atac")
        if truthy(row.get(f"{modality}_{prefix}_replacement_like", "")) or
        math.isfinite(values[modality]["fraction"]) or
        math.isfinite(values[modality]["profile_high"]))
    if boundary:
        return "replacement-like or upper-boundary fit"
    stable = []
    for modality in ("rna", "atac"):
        item = values[modality]
        stable.append(
            math.isfinite(item["fraction"]) and 0.01 < item["fraction"] <= 0.50 and
            math.isfinite(item["profile_low"]) and item["profile_low"] > 0.0 and
            math.isfinite(item["profile_high"]) and item["profile_high"] < 0.90 and
            math.isfinite(item["margin"]) and item["margin"] > 0.01 and
            math.isfinite(item["fold_support"]) and item["fold_support"] >= 0.80 and
            math.isfinite(item["influence"]) and item["influence"] <= 0.50)
    if preferred_models["rna"] != "INTERIOR_LOCKED_PLUS_CONTRIBUTOR":
        return "unavailable"
    return "mixture-compatible in both assays" if all(stable) else \
        "addition-compatible but uncertain"


def molecule_basis(row, modality):
    umi = finite(row.get(f"{modality}_molecule_genetic_umi_gene_fraction", ""))
    qname = finite(row.get(f"{modality}_molecule_genetic_qname_fraction", ""))
    units = finite(row.get(f"{modality}_molecule_genetic_units", ""))
    if (math.isfinite(units) and units <= 0) or not all(
            math.isfinite(value) and 0.0 <= value <= 1.0
            for value in (umi, qname)) or not math.isclose(
                umi + qname, 1.0, rel_tol=0.0, abs_tol=1e-9):
        return "UNAVAILABLE"
    if umi >= 0.999:
        return "UMI_GENE_LINKED"
    if qname >= 0.999:
        return "QNAME_FALLBACK"
    return "MIXED_LINKAGE_BASIS"


def midranks(values):
    finite_pairs = sorted(
        ((value, index) for index, value in enumerate(values)
         if math.isfinite(value)), key=lambda item: (item[0], item[1]))
    result = [math.nan] * len(values)
    total = len(finite_pairs)
    cursor = 0
    while cursor < total:
        end = cursor + 1
        while end < total and finite_pairs[end][0] == finite_pairs[cursor][0]:
            end += 1
        rank = ((cursor + 1) + end) / 2.0 / total
        for _, index in finite_pairs[cursor:end]:
            result[index] = rank
        cursor = end
    return result


def assign_coverage_quintiles(cells):
    by_library = defaultdict(list)
    for index, row in enumerate(cells):
        by_library[row["library"]].append(index)
    edges = []
    for library, indices in sorted(by_library.items()):
        for modality in ("rna", "atac"):
            specifications = (
                ("site", (
                    f"{modality}_candidate_independent_log1p_sites",
                    f"{modality}_candidate_independent_log1p_depth")),
                ("molecule", (
                    f"{modality}_candidate_independent_log1p_molecule_units",
                    f"{modality}_candidate_independent_log1p_effective_units")),
            )
            for channel, fields in specifications:
                complete = [index for index in indices if all(
                    math.isfinite(finite(cells[index].get(field, "")))
                    for field in fields)]
                values_by_field = [[finite(cells[index].get(field, ""))
                                    for index in complete] for field in fields]
                ranks_by_field = [midranks(values) for values in values_by_field]
                rank_by_cell = {}
                for position, cell_index in enumerate(complete):
                    rank_by_cell[cell_index] = (
                        ranks_by_field[0][position],
                        ranks_by_field[1][position])
                field_cutpoints = {
                    field: [quantile(values, probability)
                            for probability in (0.2, 0.4, 0.6, 0.8)]
                    for field, values in zip(fields, values_by_field)
                }
                edge_json = json.dumps(field_cutpoints, sort_keys=True)
                for cell_index in indices:
                    prefix = f"{modality}_candidate_independent_coverage" \
                        if channel == "site" else \
                        f"{modality}_candidate_independent_molecule_opportunity"
                    pair = rank_by_cell.get(cell_index, (math.nan, math.nan))
                    value = (pair[0] + pair[1]) / 2.0 \
                        if all(math.isfinite(item) for item in pair) else math.nan
                    cells[cell_index][f"{prefix}_field1_rank"] = pair[0]
                    cells[cell_index][f"{prefix}_field2_rank"] = pair[1]
                    cells[cell_index][f"{prefix}_rank"] = value
                    cells[cell_index][f"{prefix}_eligible_population_n"] = \
                        len(complete)
                    cells[cell_index][f"{prefix}_missing"] = not math.isfinite(value)
                    cells[cell_index][f"{prefix}_quintile"] = (
                        min(5, max(1, int(math.ceil(value * 5.0))))
                        if math.isfinite(value) else "UNAVAILABLE")
                    cells[cell_index][f"{prefix}_quintile_edges"] = edge_json
                edges.append({
                    "library": library, "modality": modality.upper(),
                    "coverage_channel": channel.upper(),
                    "candidate_independent_fields": ",".join(fields),
                    "field1_q20": field_cutpoints[fields[0]][0],
                    "field1_q40": field_cutpoints[fields[0]][1],
                    "field1_q60": field_cutpoints[fields[0]][2],
                    "field1_q80": field_cutpoints[fields[0]][3],
                    "field2_q20": field_cutpoints[fields[1]][0],
                    "field2_q40": field_cutpoints[fields[1]][1],
                    "field2_q60": field_cutpoints[fields[1]][2],
                    "field2_q80": field_cutpoints[fields[1]][3],
                    "eligible_cells": len(complete),
                    "formation_rule": (
                        "within-library mean midrank; no outcome or winner identity"),
                })
    return edges


def load_ledgers(path):
    rows = []
    by_key = {}
    for raw in read_tsv(path):
        row = dict(raw)
        key = (clean(row.get("library", "")), clean(row.get("barcode", "")))
        if key in by_key:
            raise RuntimeError(f"duplicate aggregate ledger cell: {key[0]}:{key[1]}")
        by_key[key] = row
        rows.append(row)
    return rows, by_key


def _candidate_menu_counts(candidates):
    return {
        "menu_size": len(candidates),
        "menu_physical_states": sum(candidate["physical"] for candidate in candidates),
        "menu_component_only_states": sum(
            candidate["component_only"] for candidate in candidates),
        "menu_new_donor_states": sum(candidate["new_donor"] for candidate in candidates),
        "menu_reweighting_states": sum(candidate["reweighting"] for candidate in candidates),
        "menu_genetically_identical_states": sum(
            candidate["genetically_identical"] for candidate in candidates),
        "menu_fallback_states": sum(candidate["fallback"] for candidate in candidates),
    }


def analyze_candidate_group(key, candidates, ledger, opportunity, candidate_counts):
    library, barcode = key
    locked = clean(ledger.get("reconciled_identity_locked", ""))
    menu_signature = canonical_menu_signature(candidates)
    menu_counts = _candidate_menu_counts(candidates)
    row = dict(ledger)
    row.update({
        "schema_version": SCHEMA,
        "library": library,
        "barcode": barcode,
        "canonical_menu_signature": menu_signature,
        "canonical_menu_signature_sha256": hashlib.sha256(
            menu_signature.encode()).hexdigest(),
        **menu_counts,
    })
    rankings = {}
    for modality in ("rna", "atac"):
        for evidence in ("site", "molecule"):
            rankings[(modality, evidence, "unrestricted")] = top_two(
                candidates, modality, evidence)
            rankings[(modality, evidence, "genetic")] = top_two(
                candidates, modality, evidence, genetic_candidate)
            rankings[(modality, evidence, "component_new_donor")] = top_two(
                candidates, modality, evidence, component_new_donor_candidate)
            for endpoint in ("unrestricted", "genetic", "component_new_donor"):
                prefix = f"{modality}_{evidence}_{endpoint}"
                row.update(flatten_top(prefix, rankings[(modality, evidence, endpoint)]))
            row[f"{modality}_{evidence}_genetic_evaluable_mask"] = \
                evaluable_mask(candidates, modality, evidence)
            row[f"{modality}_{evidence}_genetic_evaluable_count"] = len([
                candidate for candidate in candidates
                if genetic_candidate(
                    candidate, candidate[f"{modality}_{evidence}"])
            ])

    for modality in ("rna", "atac"):
        site_values = []
        depth_values = []
        for candidate in candidates:
            metric = candidate[f"{modality}_site"]
            if not candidate["state_genotype_equivalent"] and \
                    not metric["genotype_equivalent"]:
                if math.isfinite(metric["discriminating_sites"]):
                    site_values.append(math.log1p(max(
                        metric["discriminating_sites"], 0.0)))
                if math.isfinite(metric["discriminating_depth"]):
                    depth_values.append(math.log1p(max(
                        metric["discriminating_depth"], 0.0)))
        row[f"{modality}_candidate_independent_log1p_sites"] = median(site_values)
        row[f"{modality}_candidate_independent_log1p_depth"] = median(depth_values)
        row[f"{modality}_genotype_distinguishable_candidate_rows"] = max(
            len(site_values), len(depth_values))
        molecule_units = []
        molecule_effective = []
        for candidate in candidates:
            metric = candidate[f"{modality}_molecule"]
            if genetic_candidate(candidate, metric):
                if math.isfinite(metric["units"]):
                    molecule_units.append(math.log1p(max(metric["units"], 0.0)))
                if math.isfinite(metric["effective_units"]):
                    molecule_effective.append(math.log1p(max(
                        metric["effective_units"], 0.0)))
        row[f"{modality}_candidate_independent_log1p_molecule_units"] = \
            median(molecule_units)
        row[f"{modality}_candidate_independent_log1p_effective_units"] = \
            median(molecule_effective)

    rna_site = rankings[("rna", "site", "unrestricted")]
    atac_site = rankings[("atac", "site", "unrestricted")]
    rna_genetic = rankings[("rna", "site", "genetic")]
    atac_genetic = rankings[("atac", "site", "genetic")]
    rna_molecule = rankings[("rna", "molecule", "genetic")]
    atac_molecule = rankings[("atac", "molecule", "genetic")]
    row["original_site_both_assays_p95"] = bool(
        rna_site and atac_site and
        rna_site[0][1]["percentile"] >= 0.95 and
        atac_site[0][1]["percentile"] >= 0.95)
    row["original_site_exact_winner_agreement"] = bool(
        row["original_site_both_assays_p95"] and
        rna_site[0][0]["second_state"] == atac_site[0][0]["second_state"])
    row["genetic_source_both_assays_p95"] = bool(
        rna_genetic and atac_genetic and
        rna_genetic[0][1]["percentile"] >= 0.95 and
        atac_genetic[0][1]["percentile"] >= 0.95)
    row["genetic_source_exact_winner_agreement"] = bool(
        row["genetic_source_both_assays_p95"] and
        rna_genetic[0][0]["second_state"] == atac_genetic[0][0]["second_state"])
    row["molecule_genetic_both_assays_p95"] = bool(
        rna_molecule and atac_molecule and
        rna_molecule[0][1]["percentile"] >= 0.95 and
        atac_molecule[0][1]["percentile"] >= 0.95)
    row["molecule_genetic_exact_winner_agreement"] = bool(
        row["molecule_genetic_both_assays_p95"] and
        rna_molecule[0][0]["second_state"] ==
        atac_molecule[0][0]["second_state"])
    # Category assignment needs the complete status distribution so it can
    # distinguish a low/limited menu from a wholly unavailable one.
    row["_candidates"] = candidates
    row["site_existing_evidence_category"] = evidence_category(row, "site_genetic")
    row["molecule_existing_evidence_category"] = evidence_category(
        row, "molecule_genetic")
    row["site_to_molecule_category_changed"] = (
        row["site_existing_evidence_category"] !=
        row["molecule_existing_evidence_category"])
    row["rna_molecule_retains_unrestricted_site_winner"] = bool(
        rna_site and rankings[("rna", "molecule", "unrestricted")] and
        rna_site[0][0]["second_state"] ==
        rankings[("rna", "molecule", "unrestricted")][0][0]["second_state"])
    row["atac_molecule_retains_unrestricted_site_winner"] = bool(
        atac_site and rankings[("atac", "molecule", "unrestricted")] and
        atac_site[0][0]["second_state"] ==
        rankings[("atac", "molecule", "unrestricted")][0][0]["second_state"])
    row["molecule_retains_original_concordant_winner"] = bool(
        row["original_site_exact_winner_agreement"] and
        rankings[("rna", "molecule", "unrestricted")] and
        rankings[("atac", "molecule", "unrestricted")] and
        rna_site[0][0]["second_state"] ==
        rankings[("rna", "molecule", "unrestricted")][0][0]["second_state"] ==
        rankings[("atac", "molecule", "unrestricted")][0][0]["second_state"])
    row["rna_molecule_evidence_basis_class"] = molecule_basis(row, "rna")
    row["atac_molecule_evidence_basis_class"] = molecule_basis(row, "atac")

    winner_ids = {}
    for modality in ("rna", "atac"):
        ranked = rankings[(modality, "site", "unrestricted")]
        winner_ids[modality] = ranked[0][0]["candidate_id"] if ranked else ""
    concordant_id = winner_ids["rna"] \
        if winner_ids["rna"] and winner_ids["rna"] == winner_ids["atac"] else ""
    for candidate in candidates:
        candidate_counts[(library, candidate["second_state"])] += 1
        group_key = (
            library, locked, candidate["second_state"],
            candidate["origin_class"], candidate["candidate_origin"],
            candidate["candidate_policy"], candidate["physical"],
            candidate["component_only"], candidate["relationship"],
            candidate["fallback"], candidate["nomination_modalities"],
            menu_signature,
        )
        values = opportunity[group_key]
        values["offered_cells"] += 1
        values["rna_wins"] += candidate["candidate_id"] == winner_ids["rna"]
        values["atac_wins"] += candidate["candidate_id"] == winner_ids["atac"]
        values["concordant_wins"] += candidate["candidate_id"] == concordant_id
        values["menu_size_total"] += menu_counts["menu_size"]
        values["menu_physical_total"] += menu_counts["menu_physical_states"]
        values["menu_component_total"] += menu_counts["menu_component_only_states"]

    row["transferable_geometry_signature"] = transferable_geometry_signature(row)
    return row


def stream_candidate_analysis(candidates_path, ledger_by_key):
    cells = []
    opportunity = defaultdict(Counter)
    opportunity_denominators = Counter()
    candidate_counts = Counter()
    score_rows_by_library = Counter()
    finished_cell_keys = set()
    current_key = None
    current_candidate_ids = set()
    current_candidates = []

    def finish(key, candidates):
        if key is None:
            return
        ledger = ledger_by_key.get(key)
        if ledger is None:
            raise RuntimeError(
                f"candidate cell absent from aggregate ledger: {key[0]}:{key[1]}")
        analyzed = analyze_candidate_group(
            key, candidates, ledger, opportunity, candidate_counts)
        cells.append(analyzed)
        if candidates:
            opportunity_denominators[(
                key[0], clean(analyzed.get("reconciled_identity_locked", ""))
            )] += 1

    for raw in read_tsv(candidates_path):
        key = (clean(raw.get("library", "")), clean(raw.get("barcode", "")))
        candidate_id = clean(raw.get("candidate_id", ""))
        if current_key is not None and key != current_key:
            finish(current_key, current_candidates)
            finished_cell_keys.add(current_key)
            current_candidates = []
            current_candidate_ids = set()
        if key in finished_cell_keys:
            raise RuntimeError(
                "aggregate candidate table has a noncontiguous cell: " +
                ":".join(key))
        if candidate_id in current_candidate_ids:
            raise RuntimeError(
                "duplicate aggregate candidate key: " +
                ":".join(key + (candidate_id,)))
        current_key = key
        current_candidate_ids.add(candidate_id)
        current_candidates.append(compact_candidate(raw))
        score_rows_by_library[key[0]] += 1
    finish(current_key, current_candidates)

    cells_by_key = {(row["library"], row["barcode"]): row for row in cells}
    for key, ledger in ledger_by_key.items():
        if key not in cells_by_key:
            row = dict(ledger)
            row.update({
                "schema_version": SCHEMA,
                "library": key[0], "barcode": key[1],
                "canonical_menu_signature": "",
                "canonical_menu_signature_sha256": "",
                "menu_size": 0,
                "original_site_both_assays_p95": False,
                "original_site_exact_winner_agreement": False,
                "genetic_source_both_assays_p95": False,
                "genetic_source_exact_winner_agreement": False,
                "molecule_genetic_both_assays_p95": False,
                "molecule_genetic_exact_winner_agreement": False,
                "site_existing_evidence_category": "insufficient evidence",
                "molecule_existing_evidence_category": "insufficient evidence",
                "site_to_molecule_category_changed": False,
                "rna_candidate_independent_log1p_sites": math.nan,
                "rna_candidate_independent_log1p_depth": math.nan,
                "atac_candidate_independent_log1p_sites": math.nan,
                "atac_candidate_independent_log1p_depth": math.nan,
                "rna_candidate_independent_log1p_molecule_units": math.nan,
                "rna_candidate_independent_log1p_effective_units": math.nan,
                "atac_candidate_independent_log1p_molecule_units": math.nan,
                "atac_candidate_independent_log1p_effective_units": math.nan,
                "_candidates": [],
            })
            cells.append(row)
    cells.sort(key=lambda row: (row["library"], row["barcode"]))
    coverage_edges = assign_coverage_quintiles(cells)
    return (cells, opportunity, opportunity_denominators, candidate_counts,
            score_rows_by_library, coverage_edges)


def opportunity_rows(opportunity, denominators):
    rows = []
    for key, values in sorted(opportunity.items()):
        (library, locked, candidate, origin_class, candidate_origin,
         candidate_policy, physical, component_only, relationship, fallback,
         nomination_modalities, signature) = key
        offered = values["offered_cells"]
        eligible = denominators[(library, locked)]
        expected = (values["rna_wins"] * values["atac_wins"] / offered) \
            if offered else math.nan
        rows.append({
            "analysis_type": "CANDIDATE_OPPORTUNITY_AND_CONCENTRATION",
            "endpoint": "SITE_UNRESTRICTED_WINNER",
            "scope": "LIBRARY_LOCKED_CANDIDATE_ORIGIN_MENU",
            "library": library,
            "locked_identity": locked,
            "proposed_contributor": candidate,
            "candidate_origin_class": origin_class,
            "candidate_origin": candidate_origin,
            "candidate_policy": candidate_policy,
            "physical_pool_state": physical,
            "component_only_state": component_only,
            "structural_added_state_relationship": relationship,
            "exhaustive_fallback": fallback,
            "nomination_modalities": nomination_modalities,
            "canonical_menu_signature": signature,
            "offered_cells": offered,
            "scoreable_cells_with_locked_identity": eligible,
            "candidate_opportunity_rate": offered / eligible
                if eligible else math.nan,
            "rna_wins": values["rna_wins"],
            "atac_wins": values["atac_wins"],
            "observed_concordant_wins": values["concordant_wins"],
            "expected_concordant_wins": expected,
            "excess_concordant_wins": values["concordant_wins"] - expected,
            "mean_menu_size": values["menu_size_total"] / offered if offered else math.nan,
            "mean_physical_menu_states": values["menu_physical_total"] / offered
                if offered else math.nan,
            "mean_component_only_menu_states": values["menu_component_total"] / offered
                if offered else math.nan,
        })
    return rows


def endpoint_top(row, endpoint, modality):
    prefix = {
        "SITE_UNRESTRICTED": "site_unrestricted",
        "SITE_GENETIC_SOURCE": "site_genetic",
        "MOLECULE_GENETIC_SOURCE": "molecule_genetic",
    }[endpoint]
    state = clean(row.get(f"{modality}_{prefix}_second_state", ""))
    return {
        "winner": state,
        "candidate_id": clean(row.get(f"{modality}_{prefix}_candidate_id", "")),
        "percentile": finite(row.get(f"{modality}_{prefix}_percentile", "")),
        "fraction": finite(row.get(f"{modality}_{prefix}_fraction", "")),
        "profile_low": finite(row.get(f"{modality}_{prefix}_profile_low", "")),
        "profile_high": finite(row.get(f"{modality}_{prefix}_profile_high", "")),
        "margin": finite(row.get(f"{modality}_{prefix}_margin", "")),
        "warnings": clean(row.get(f"{modality}_{prefix}_warnings", "")),
        "available": bool(state),
    }


def documented_technical_state(row):
    fields = (
        "reconciled_droplet_state", "current_droplet_flag",
        "occupancy_resolution_status", "competing_technical_state",
        "explicit_multiplet_evidence",
    )
    values = []
    for field in fields:
        value = clean(row.get(field, ""))
        if value and value.upper() not in {"NONE", "UNAVAILABLE", "UNKNOWN"}:
            values.append(f"{field}={value}")
    return "|".join(values) if values else "UNAVAILABLE"


def _base_block(row, endpoint, include_technical, coarse_identity):
    locked = clean(row.get("reconciled_identity_locked", ""))
    identity = component_shape(locked) if coarse_identity else locked
    key = [
        row.get("library", ""), identity,
        clean(row.get("ploidy_field", "") or row.get("ploidy_call", "")) or
            "UNAVAILABLE",
        clean(row.get("biological_state", "") or
              row.get("biological_tetraploidy_field", "")) or "UNAVAILABLE",
        clean(row.get("sensitivity_menu_signature", "") or
              row.get("canonical_menu_signature", "")) or "UNAVAILABLE",
        row.get("rna_candidate_independent_coverage_quintile", "UNAVAILABLE"),
        row.get("atac_candidate_independent_coverage_quintile", "UNAVAILABLE"),
    ]
    if include_technical:
        key.append(documented_technical_state(row))
    if endpoint.startswith("MOLECULE"):
        key.extend((
            row.get("rna_molecule_evidence_basis_class", "UNAVAILABLE"),
            row.get("atac_molecule_evidence_basis_class", "UNAVAILABLE"),
        ))
    return tuple(key)


def sparse_permutation_blocks(rows, endpoint):
    full_keys = [_base_block(row, endpoint, True, False) for row in rows]
    full_counts = Counter(full_keys)
    levels = ["FULL" if full_counts[key] >= 20 else "DROP_TECHNICAL"
              for key in full_keys]
    drop_keys = [_base_block(row, endpoint, False, False) for row in rows]
    drop_counts = Counter(
        key for key, level in zip(drop_keys, levels)
        if level == "DROP_TECHNICAL")
    for index, level in enumerate(levels):
        if level == "DROP_TECHNICAL" and drop_counts[drop_keys[index]] < 20:
            levels[index] = "COARSE_IDENTITY"
    coarse_keys = [_base_block(row, endpoint, False, True) for row in rows]
    final_keys = []
    for index, level in enumerate(levels):
        key = (level,) + (
            full_keys[index] if level == "FULL" else
            drop_keys[index] if level == "DROP_TECHNICAL" else
            coarse_keys[index])
        final_keys.append(key)
    final_counts = Counter(final_keys)
    accounting = {
        "cells": len(rows),
        "initial_blocks": len(full_counts),
        "initial_cells_in_blocks_ge20": sum(
            full_counts[key] for key in full_counts if full_counts[key] >= 20),
        "drop_technical_blocks": len({
            drop_keys[index] for index, level in enumerate(levels)
            if level == "DROP_TECHNICAL"}),
        "drop_technical_cells": sum(level == "DROP_TECHNICAL" for level in levels),
        "coarse_identity_blocks": len({
            coarse_keys[index] for index, level in enumerate(levels)
            if level == "COARSE_IDENTITY"}),
        "coarse_identity_cells": sum(level == "COARSE_IDENTITY" for level in levels),
        "final_blocks": len(final_counts),
        "exchangeable_multirow_blocks": sum(count > 1 for count in final_counts.values()),
        "exchangeable_multirow_cells": sum(
            count for count in final_counts.values() if count > 1),
        "singleton_blocks": sum(count == 1 for count in final_counts.values()),
        "singleton_cells": sum(count for count in final_counts.values() if count == 1),
    }
    return final_keys, accounting


def prescribed_block_keys(rows, endpoint, definition):
    if definition == "ADAPTIVE_PREDECLARED":
        return sparse_permutation_blocks(rows, endpoint)
    if definition == "FULL_EXACT":
        keys = [("FULL_EXACT",) + _base_block(row, endpoint, True, False)
                for row in rows]
    elif definition == "WITHOUT_TECHNICAL_STATE":
        keys = [("WITHOUT_TECHNICAL_STATE",) +
                _base_block(row, endpoint, False, False) for row in rows]
    elif definition == "COARSE_IDENTITY":
        keys = [("COARSE_IDENTITY",) +
                _base_block(row, endpoint, False, True) for row in rows]
    else:
        raise ValueError(f"unsupported block definition: {definition}")
    counts = Counter(keys)
    return keys, {
        "cells": len(rows), "initial_blocks": len(counts),
        "initial_cells_in_blocks_ge20": sum(
            count for count in counts.values() if count >= 20),
        "drop_technical_blocks": 0, "drop_technical_cells": 0,
        "coarse_identity_blocks": 0, "coarse_identity_cells": 0,
        "final_blocks": len(counts),
        "exchangeable_multirow_blocks": sum(
            count > 1 for count in counts.values()),
        "exchangeable_multirow_cells": sum(
            count for count in counts.values() if count > 1),
        "singleton_blocks": sum(count == 1 for count in counts.values()),
        "singleton_cells": sum(count == 1 for count in counts.values()),
    }


def block_bootstrap_interval(rows, rna, atac_winners, strata, replicates,
                             seed_label, exchangeable_only=False):
    blocks = defaultdict(list)
    for index, key in enumerate(strata):
        blocks[(rows[index]["library"], key)].append(index)
    by_library = defaultdict(list)
    for (library, _), indices in blocks.items():
        if exchangeable_only and len(indices) == 1:
            continue
        size = len(indices)
        observed = sum(rna[index] == atac_winners[index] for index in indices)
        rna_counts = Counter(rna[index] for index in indices)
        atac_counts = Counter(atac_winners[index] for index in indices)
        expected = size * sum(
            (rna_counts[value] / size) * (atac_counts[value] / size)
            for value in set(rna_counts) | set(atac_counts))
        by_library[library].append((size, observed, expected))
    rng = random.Random(stable_seed("block_bootstrap", seed_label,
                                    exchangeable_only))
    effects = []
    for _ in range(replicates):
        total_n = total_observed = 0.0
        total_expected = 0.0
        for library in sorted(by_library):
            source = by_library[library]
            for _ in range(len(source)):
                size, observed, expected = source[rng.randrange(len(source))]
                total_n += size
                total_observed += observed
                total_expected += expected
        if total_n:
            effects.append((total_observed - total_expected) / total_n)
    return {
        "bootstrap_replicates": replicates,
        "bootstrap_excess_95pct_low": quantile(effects, 0.025),
        "bootstrap_excess_95pct_high": quantile(effects, 0.975),
        "bootstrap_unit": "final blocks resampled with replacement within library",
    }


def analytic_agreement(rna, atac, strata):
    groups = defaultdict(list)
    for index, key in enumerate(strata):
        groups[key].append(index)
    expected = 0.0
    for indices in groups.values():
        rna_counts = Counter(rna[index] for index in indices)
        atac_counts = Counter(atac[index] for index in indices)
        size = len(indices)
        expected += size * sum(
            (rna_counts[value] / size) * (atac_counts[value] / size)
            for value in set(rna_counts) | set(atac_counts))
    return expected


def block_permutation_agreement_counts(rna, atac, groups, permutations, seed):
    """Sample complete-block permutations through their exact contingency law."""
    rng = np.random.default_rng(seed)
    totals = np.zeros(permutations, dtype=np.int64)
    for indices in groups.values():
        if len(indices) == 1:
            totals += int(rna[indices[0]] == atac[indices[0]]["winner"])
            continue
        rna_values = [rna[index] for index in indices]
        atac_values = [atac[index]["winner"] for index in indices]
        labels = sorted(set(rna_values) | set(atac_values))
        label_index = {label: index for index, label in enumerate(labels)}
        rna_counts = np.zeros(len(labels), dtype=np.int64)
        atac_counts = np.zeros(len(labels), dtype=np.int64)
        for value in rna_values:
            rna_counts[label_index[value]] += 1
        for value in atac_values:
            atac_counts[label_index[value]] += 1
        remaining = np.broadcast_to(
            atac_counts, (permutations, len(labels))).copy()
        block_matches = np.zeros(permutations, dtype=np.int64)
        for rna_index, row_size in enumerate(rna_counts):
            if row_size <= 0:
                continue
            allocations = np.zeros_like(remaining)
            draws_left = np.full(permutations, row_size, dtype=np.int64)
            colors_left = remaining.sum(axis=1)
            for atac_index in range(len(labels) - 1):
                good = remaining[:, atac_index]
                bad = colors_left - good
                draw = rng.hypergeometric(good, bad, draws_left)
                allocations[:, atac_index] = draw
                draws_left -= draw
                colors_left -= good
            allocations[:, -1] = draws_left
            block_matches += allocations[:, rna_index]
            remaining -= allocations
        totals += block_matches
    return totals.tolist()


def permutation_agreement(rows, endpoint, permutations, seed_label,
                          block_definition="ADAPTIVE_PREDECLARED"):
    input_cells = len(rows)
    rows = [row for row in rows
            if endpoint_top(row, endpoint, "rna")["available"] and
            endpoint_top(row, endpoint, "atac")["available"]]
    rna = [endpoint_top(row, endpoint, "rna")["winner"] for row in rows]
    atac = [endpoint_top(row, endpoint, "atac") for row in rows]
    strata, accounting = prescribed_block_keys(
        rows, endpoint, block_definition)
    groups = defaultdict(list)
    for index, key in enumerate(strata):
        groups[key].append(index)
    observed = sum(rna[index] == atac[index]["winner"]
                   for index in range(len(rows)))
    expected = analytic_agreement(
        rna, [item["winner"] for item in atac], strata)
    singleton_indices = {
        index for indices in groups.values() if len(indices) == 1
        for index in indices
    }
    singleton_agreement = sum(
        rna[index] == atac[index]["winner"] for index in singleton_indices)
    null = block_permutation_agreement_counts(
        rna, atac, groups, permutations,
        stable_seed("block_permutation", seed_label))
    exceedance = sum(value >= observed for value in null)
    denominator = len(rows)
    exchangeable_index_set = set(range(len(rows))) - singleton_indices
    exchangeable_indices = sorted(exchangeable_index_set)
    exchangeable_groups = {}
    for key, indices in groups.items():
        kept = [index for index in indices if index in exchangeable_index_set]
        if kept:
            exchangeable_groups[key] = kept
    exchangeable_observed = sum(
        rna[index] == atac[index]["winner"] for index in exchangeable_indices)
    exchangeable_expected = analytic_agreement(
        [rna[index] for index in exchangeable_indices],
        [atac[index]["winner"] for index in exchangeable_indices],
        [strata[index] for index in exchangeable_indices]) \
        if exchangeable_indices else 0.0
    exchangeable_null = block_permutation_agreement_counts(
        rna, atac, exchangeable_groups, permutations,
        stable_seed("block_permutation", seed_label, "exchangeable_only")) \
        if exchangeable_groups else []
    exchangeable_exceedance = sum(
        value >= exchangeable_observed for value in exchangeable_null)
    result = {
        "observed_agreements": observed,
        "input_cells": input_cells,
        "excluded_missing_endpoint": input_cells - len(rows),
        "denominator": denominator,
        "observed_agreement_fraction": observed / denominator
            if denominator else math.nan,
        "analytic_expected_agreements": expected,
        "analytic_expected_fraction": expected / denominator
            if denominator else math.nan,
        "observed_minus_expected_fraction": (observed - expected) / denominator
            if denominator else math.nan,
        "observed_expected_ratio": observed / expected
            if expected > 0 else math.nan,
        "permutation_mean_agreements": statistics.mean(null) if null else math.nan,
        "permutation_mean_fraction": statistics.mean(null) / denominator
            if null and denominator else math.nan,
        "permutation_2.5pct_agreements": quantile(null, 0.025),
        "permutation_97.5pct_agreements": quantile(null, 0.975),
        "permutation_exceedance_count": exceedance,
        "permutations": permutations,
        "empirical_upper_tail_probability": (exceedance + 1) / (permutations + 1),
        "empirical_probability_reporting": (
            f"Monte Carlo upper-tail probability; resolution 1/{permutations + 1}"),
        "singleton_fixed_agreements": singleton_agreement,
        "singleton_fixed_cells": len(singleton_indices),
        "exchangeable_only_denominator": len(exchangeable_indices),
        "exchangeable_only_observed_agreements": exchangeable_observed,
        "exchangeable_only_expected_agreements": exchangeable_expected,
        "exchangeable_only_excess_fraction": (
            (exchangeable_observed - exchangeable_expected) /
            len(exchangeable_indices) if exchangeable_indices else math.nan),
        "exchangeable_only_empirical_upper_tail_probability": (
            (exchangeable_exceedance + 1) / (permutations + 1)
            if exchangeable_null else math.nan),
        "block_definition": block_definition,
        "fixed_singleton_materiality_warning": (
            "YES" if denominator and singleton_indices and
            abs((observed - expected) / denominator -
                ((exchangeable_observed - exchangeable_expected) /
                 len(exchangeable_indices)
                 if exchangeable_indices else 0.0)) >= 0.02 else "NO"),
        "permutation_engine": (
            "COMPLETE_ATAC_BLOCK_PERMUTATION_EQUIVALENT_"
            "MULTIVARIATE_HYPERGEOMETRIC_AGREEMENT_COUNTS"),
        **accounting,
    }
    result.update(block_bootstrap_interval(
        rows, rna, [item["winner"] for item in atac], strata,
        BOOTSTRAP_REPLICATES, seed_label, False))
    exchangeable_bootstrap = block_bootstrap_interval(
        rows, rna, [item["winner"] for item in atac], strata,
        BOOTSTRAP_REPLICATES, seed_label, True)
    result.update({
        f"exchangeable_only_{key}": value
        for key, value in exchangeable_bootstrap.items()
    })
    return result


def null_result_row(result, endpoint, scope, library="", exclusion="NONE"):
    return {
        "analysis_type": "CANDIDATE_AND_COVERAGE_AWARE_CROSS_ASSAY_NULL",
        "endpoint": endpoint,
        "endpoint_plain_language": endpoint_label(endpoint),
        "scope": scope,
        "library": library,
        "exclusion": exclusion,
        **result,
        "seed": SEED,
        "coverage_definition": (
            "within-library average midrank of median log1p discriminating "
            "sites and depth across all legal genotype-distinguishable candidates"),
        "sparse_block_rule": (
            "drop technical state below 20; then replace exact identity with "
            "component and distinct-donor counts; singleton blocks fixed"),
        "interval_definition": (
            "95% block-bootstrap interval for observed-minus-expected excess; "
            "permutation quantiles are null-distribution diagnostics only"),
    }


def endpoint_focal(cells, endpoint):
    if endpoint == "SITE_UNRESTRICTED":
        return [row for row in cells if truthy(
            row.get("original_site_both_assays_p95", ""))]
    if endpoint == "SITE_GENETIC_SOURCE":
        return [row for row in cells if truthy(
            row.get("genetic_source_both_assays_p95", ""))]
    if endpoint == "SITE_MIXTURE_OR_UNCERTAIN":
        return [row for row in cells
                if truthy(row.get("genetic_source_both_assays_p95", "")) and
                row.get("site_existing_evidence_category") in {
                    "mixture-compatible in both assays",
                    "addition-compatible but uncertain"}]
    if endpoint == "MOLECULE_GENETIC_SOURCE":
        return [row for row in cells if truthy(
            row.get("molecule_genetic_both_assays_p95", ""))]
    raise ValueError(f"unsupported endpoint: {endpoint}")


def run_null_suite(cells, permutations):
    rows = []
    endpoint_map = {
        "SITE_UNRESTRICTED": "SITE_UNRESTRICTED",
        "SITE_GENETIC_SOURCE": "SITE_GENETIC_SOURCE",
        "SITE_MIXTURE_OR_UNCERTAIN": "SITE_GENETIC_SOURCE",
        "MOLECULE_GENETIC_SOURCE": "MOLECULE_GENETIC_SOURCE",
    }
    for endpoint, winner_endpoint in endpoint_map.items():
        focal = endpoint_focal(cells, endpoint)
        primary = [row for row in focal if int(row["library"].removeprefix("lib"))
                   in PRIMARY_LIBRARIES]
        if endpoint == "SITE_MIXTURE_OR_UNCERTAIN":
            rows.append({
                "analysis_type": "SELECTION_CONDITIONED_AUDIT",
                "endpoint": endpoint,
                "endpoint_plain_language": (
                    "strict two-source mixture or addition-compatible but weakly separated"),
                "scope": "PRIMARY_SIX_LIBRARIES",
                "denominator": len(primary),
                "observed_agreements": len(primary),
                "inferential_status": (
                    "agreement fixed by selection; no valid agreement test"),
                "empirical_upper_tail_probability": math.nan,
            })
            continue
        calibration_name = f"lib{CALIBRATION_LIBRARY}"
        calibration = [row for row in focal
                       if row["library"] == calibration_name]
        result = permutation_agreement(
            primary, winner_endpoint, permutations, f"{endpoint}:PRIMARY")
        rows.append(null_result_row(result, endpoint, "PRIMARY_SIX_LIBRARIES"))
        if calibration:
            result = permutation_agreement(
                calibration, winner_endpoint, permutations,
                f"{endpoint}:CALIBRATION_{calibration_name.upper()}")
            rows.append(null_result_row(
                result, endpoint, "CALIBRATION_SENSITIVITY", calibration_name))
        per_library = []
        for library in (f"lib{value}" for value in PRIMARY_LIBRARIES):
            subset = [row for row in primary if row["library"] == library]
            if subset:
                result = permutation_agreement(
                    subset, winner_endpoint, permutations,
                    f"{endpoint}:{library}")
                output = null_result_row(
                    result, endpoint, "PER_LIBRARY", library)
                per_library.append(output)
                rows.append(output)
        tested = len(per_library)
        for output in per_library:
            raw = finite(output.get("empirical_upper_tail_probability", ""))
            output["multiplicity_adjustment_rule"] = \
                f"Bonferroni across {tested} primary-library tests"
            output["multiplicity_adjusted_probability"] = min(
                1.0, raw * tested) if math.isfinite(raw) else math.nan
        for omitted in (f"lib{value}" for value in PRIMARY_LIBRARIES):
            subset = [row for row in primary if row["library"] != omitted]
            if subset:
                result = permutation_agreement(
                    subset, winner_endpoint, permutations,
                    f"{endpoint}:LEAVE_OUT:{omitted}")
                rows.append(null_result_row(
                    result, endpoint, "LEAVE_ONE_LIBRARY_OUT", omitted,
                    f"OMIT_{omitted}"))
        if endpoint == "SITE_GENETIC_SOURCE" and primary:
            for definition in (
                    "FULL_EXACT", "WITHOUT_TECHNICAL_STATE",
                    "COARSE_IDENTITY"):
                result = permutation_agreement(
                    primary, winner_endpoint, permutations,
                    f"{endpoint}:BLOCK_SENSITIVITY:{definition}", definition)
                rows.append(null_result_row(
                    result, endpoint, "PREDECLARED_BLOCK_SENSITIVITY",
                    exclusion=definition))
    return rows


def _delta_rank(candidates, modality, evidence, allowed):
    ranked = []
    for candidate in candidates:
        metric = candidate[f"{modality}_{evidence}"]
        if genetic_candidate(candidate, metric) and allowed(candidate):
            copied = dict(metric)
            copied["percentile"] = math.nan
            ranked.append((candidate, copied))
    ranked.sort(key=lambda item: (
        item[1]["delta"], item[0]["candidate_id"]), reverse=True)
    return ranked[:2]


def rerank_endpoint_population(cells, exclusion_kind, exclusion_value,
                               frozen_keys):
    frozen_ids = {f"{library}:{barcode}" for library, barcode in frozen_keys}
    population = []
    for source in cells:
        if int(source["library"].removeprefix("lib")) not in PRIMARY_LIBRARIES:
            continue
        locked = clean(source.get("reconciled_identity_locked", ""))

        def allowed(candidate):
            if exclusion_kind == "candidate" and \
                    candidate["second_state"] == exclusion_value:
                return False
            if exclusion_kind == "pattern" and \
                    (locked, candidate["second_state"]) == exclusion_value:
                return False
            return True

        candidates = [candidate for candidate in source.get("_candidates", [])
                      if allowed(candidate)]
        row = dict(source)
        row["_sensitivity_candidates"] = candidates
        row["_candidates"] = candidates
        opportunity = [candidate for candidate in candidates if any(
            genetic_candidate(candidate, candidate[f"{assay}_{mode}"])
            for assay in ("rna", "atac") for mode in ("site", "molecule"))]
        row["sensitivity_menu_signature"] = canonical_menu_signature(opportunity)
        row["canonical_menu_signature"] = row["sensitivity_menu_signature"]
        row["menu_size"] = len(opportunity)
        for evidence in ("site", "molecule"):
            for modality in ("rna", "atac"):
                ranked = _delta_rank(candidates, modality, evidence, allowed)
                row[f"_sensitivity_{modality}_{evidence}_ranked"] = ranked
                row.update(flatten_top(
                    f"{modality}_{evidence}_genetic", ranked))
                row[f"{modality}_{evidence}_genetic_evaluable_mask"] = \
                    evaluable_mask(candidates, modality, evidence)
                row[f"{modality}_{evidence}_genetic_evaluable_count"] = len([
                    candidate for candidate in candidates
                    if genetic_candidate(
                        candidate, candidate[f"{modality}_{evidence}"])
                ])
        population.append(row)

    for evidence in ("site", "molecule"):
        exact = defaultdict(list)
        coarsened = defaultdict(list)
        for row in population:
            coverage = (
                row.get("rna_candidate_independent_coverage_quintile"),
                row.get("atac_candidate_independent_coverage_quintile"),
            ) if evidence == "site" else (
                row.get("rna_candidate_independent_coverage_quintile"),
                row.get("atac_candidate_independent_coverage_quintile"),
                row.get("rna_candidate_independent_molecule_opportunity_quintile"),
                row.get("atac_candidate_independent_molecule_opportunity_quintile"),
                row.get("rna_molecule_evidence_basis_class"),
                row.get("atac_molecule_evidence_basis_class"),
            )
            exact_key = (
                row["library"], row["sensitivity_menu_signature"],
                row.get(f"rna_{evidence}_genetic_evaluable_mask", ""),
                row.get(f"atac_{evidence}_genetic_evaluable_mask", ""),
                *coverage,
            )
            count_key = (
                row["library"], row["sensitivity_menu_signature"],
                row.get(f"rna_{evidence}_genetic_evaluable_count", 0),
                row.get(f"atac_{evidence}_genetic_evaluable_count", 0),
                *coverage,
            )
            row[f"_sensitivity_{evidence}_exact_key"] = exact_key
            row[f"_sensitivity_{evidence}_count_key"] = count_key
            for modality in ("rna", "atac"):
                ranked = row.get(
                    f"_sensitivity_{modality}_{evidence}_ranked", [])
                if ranked:
                    value = ranked[0][1]["delta"]
                    exact[(exact_key, modality)].append(
                        (row["barcode"], value))
                    coarsened[(count_key, modality)].append(
                        (row["barcode"], value))
        for row in population:
            for modality in ("rna", "atac"):
                ranked = row.get(
                    f"_sensitivity_{modality}_{evidence}_ranked", [])
                if not ranked:
                    continue
                reference = [item for item in exact.get((
                    row[f"_sensitivity_{evidence}_exact_key"], modality), [])
                    if item[0] != row["barcode"] and f"{row['library']}:{item[0]}" not in frozen_ids]
                level = "EXACT_MENU_EVALUABLE_MASK_COVERAGE"
                values = sorted(item[1] for item in reference)
                target_value = ranked[0][1]["delta"]
                runner_value = ranked[1][1]["delta"] if len(ranked) > 1 \
                    else math.nan
                prefix = f"{modality}_{evidence}_genetic"
                row[f"{prefix}_percentile"] = percentile(target_value, values) \
                    if len(values) >= MIN_CALIBRATION_REFERENCE else math.nan
                row[f"{prefix}_runner_up_percentile"] = percentile(
                    runner_value, values) if len(values) >= MIN_CALIBRATION_REFERENCE \
                    and math.isfinite(runner_value) else math.nan
                row[f"{prefix}_margin"] = (
                    row[f"{prefix}_percentile"] -
                    row[f"{prefix}_runner_up_percentile"]
                    if math.isfinite(row[f"{prefix}_percentile"]) and
                    math.isfinite(row[f"{prefix}_runner_up_percentile"])
                    else math.nan)
                row[f"{prefix}_calibration_level"] = level \
                    if len(values) >= MIN_CALIBRATION_REFERENCE else "UNAVAILABLE"
                row[f"{prefix}_reference_cells"] = len(values)
    for row in population:
        row["genetic_source_both_assays_p95"] = all(
            finite(row.get(f"{modality}_site_genetic_percentile", "")) >= 0.95
            for modality in ("rna", "atac"))
        row["genetic_source_exact_winner_agreement"] = bool(
            row["genetic_source_both_assays_p95"] and
            clean(row.get("rna_site_genetic_second_state", "")) ==
            clean(row.get("atac_site_genetic_second_state", "")))
        row["molecule_genetic_both_assays_p95"] = all(
            finite(row.get(f"{modality}_molecule_genetic_percentile", "")) >= 0.95
            for modality in ("rna", "atac"))
        row["molecule_genetic_exact_winner_agreement"] = bool(
            row["molecule_genetic_both_assays_p95"] and
            clean(row.get("rna_molecule_genetic_second_state", "")) ==
            clean(row.get("atac_molecule_genetic_second_state", "")))
        row["site_existing_evidence_category"] = evidence_category(
            row, "site_genetic")
        row["molecule_existing_evidence_category"] = evidence_category(
            row, "molecule_genetic")
    return population


def dominant_sensitivity_rows(cells, permutations, frozen_keys):
    focal = [row for row in endpoint_focal(cells, "SITE_GENETIC_SOURCE")
             if int(row["library"].removeprefix("lib")) in PRIMARY_LIBRARIES]
    candidate_wins = Counter()
    pattern_wins = Counter()
    for row in focal:
        rna = endpoint_top(row, "SITE_GENETIC_SOURCE", "rna")["winner"]
        atac = endpoint_top(row, "SITE_GENETIC_SOURCE", "atac")["winner"]
        locked = clean(row.get("reconciled_identity_locked", ""))
        for winner in {rna, atac} - {""}:
            candidate_wins[winner] += 1
            pattern_wins[(locked, winner)] += 1
    candidates = {value for value, count in candidate_wins.items() if count >= 20}
    candidates.update({"H28126", "JOS3C1"})
    patterns = {value for value, count in pattern_wins.items() if count >= 20}
    patterns.add(("H21792+H29089", "H28126"))
    rows = []
    for kind, values in (("candidate", sorted(candidates)),
                         ("pattern", sorted(patterns))):
        for value in values:
            all_reranked = rerank_endpoint_population(
                cells, kind, value, frozen_keys)
            reranked = [row for row in all_reranked
                        if truthy(row.get("genetic_source_both_assays_p95", ""))]
            label = value if isinstance(value, str) else "->".join(value)
            historical_keys = {(row["library"], row["barcode"]) for row in focal}
            current_keys = {(row["library"], row["barcode"]) for row in reranked}
            fixed = [row for row in reranked if (row["library"], row["barcode"]) in historical_keys]
            for cohort, members in (("HISTORICAL_FIXED_FOCAL", fixed),
                                    ("NEWLY_ELIGIBLE_RERANKED", reranked)):
                for scope_library in [""] + [f"lib{x}" for x in PRIMARY_LIBRARIES]:
                    subset = [row for row in members if not scope_library or row["library"] == scope_library]
                    scope_historical = {key for key in historical_keys if not scope_library or key[0] == scope_library}
                    scope_current = {key for key in current_keys if not scope_library or key[0] == scope_library}
                    scope_population = [row for row in all_reranked if not scope_library or row["library"] == scope_library]
                    result = permutation_agreement(subset, "SITE_GENETIC_SOURCE", permutations,
                        f"LEAVE_ONE_{kind.upper()}:{label}:{cohort}:{scope_library}")
                    output = null_result_row(result, "SITE_GENETIC_SOURCE",
                        f"LEAVE_ONE_{kind.upper()}_RERANK", library=scope_library, exclusion=label)
                    output.update({
                        "cohort": cohort, "original_focal_cells": len(scope_historical),
                        "original_focal_wins": candidate_wins[value] if kind == "candidate" else pattern_wins[value],
                        "historical_focal_retained": len(scope_historical & scope_current),
                        "historical_focal_attrition": len(scope_historical-scope_current),
                        "newly_eligible_cells": len(scope_current-scope_historical),
                        "reranked_both_assays_p95_cells": len(subset),
                        "remaining_menu_eligible_cells": sum(bool(row.get("menu_size")) for row in scope_population),
                        "compatible_reference_unavailable_cells": sum(any(
                            row.get(f"{assay}_site_genetic_calibration_level") == "UNAVAILABLE"
                            for assay in ("rna", "atac")) for row in scope_population),
                        "percentile_reference_definition": "EXACT_REMAINING_MENU_MASK_COVERAGE_WITHIN_LIBRARY_LEAVE_ONE_CELL_OUT_DESCRIPTIVE",
                        "remaining_complete_menu_signatures": len({row["sensitivity_menu_signature"] for row in scope_population}),
                    })
                    rows.append(output)
    for label, libraries in (
            ("LIBRARIES_7_12_20", EARLY_EXCESS_LIBRARIES),
            ("LIBRARIES_9_17_29", LOW_EXCESS_LIBRARIES)):
        subset = [row for row in focal
                  if int(row["library"].removeprefix("lib")) in libraries]
        if subset:
            result = permutation_agreement(
                subset, "SITE_GENETIC_SOURCE", permutations, label)
            rows.append(null_result_row(
                result, "SITE_GENETIC_SOURCE", label))
    return rows, candidate_wins, pattern_wins


def _genetic_maximum(row, modality, evidence):
    ranked = _delta_rank(
        row.get("_candidates", []), modality, evidence, lambda _candidate: True)
    if not ranked:
        return math.nan, "", ""
    return ranked[0][1]["delta"], ranked[0][0]["second_state"], \
        ranked[0][0]["candidate_id"]


def _unrestricted_maximum(row, modality, evidence):
    ranked = []
    for candidate in row.get("_candidates", []):
        metric = candidate[f"{modality}_{evidence}"]
        if metric.get("available") and math.isfinite(metric.get("delta", math.nan)):
            ranked.append((candidate, metric))
    ranked.sort(key=lambda item: (
        item[1]["delta"], item[0]["candidate_id"]), reverse=True)
    if not ranked:
        return math.nan, "", ""
    return ranked[0][1]["delta"], ranked[0][0]["second_state"], \
        ranked[0][0]["candidate_id"]


def _calibration_key(row, evidence, level, external=False):
    site_coverage = (
        row.get("rna_candidate_independent_coverage_quintile", "UNAVAILABLE"),
        row.get("atac_candidate_independent_coverage_quintile", "UNAVAILABLE"),
    )
    molecule_fields = () if evidence == "site" else (
        row.get("rna_molecule_evidence_basis_class", "UNAVAILABLE"),
        row.get("atac_molecule_evidence_basis_class", "UNAVAILABLE"),
        row.get("rna_candidate_independent_molecule_opportunity_quintile",
                "UNAVAILABLE"),
        row.get("atac_candidate_independent_molecule_opportunity_quintile",
                "UNAVAILABLE"),
    )
    library = () if external else (row["library"],)
    # A geometry-only match changes the candidate opportunity experiment.
    # Exact canonical donor/copy menus are required even across libraries.
    menu = row.get("canonical_menu_signature", "")
    if level == "WHOLE_LIBRARY":
        return library
    if level == "EXACT_MENU":
        return library + (menu,)
    if level == "MENU_SITE_COVERAGE":
        return library + (menu,) + site_coverage
    if level == "PRIMARY_EXACT":
        return library + (
            menu,
            row.get(f"rna_{evidence}_genetic_evaluable_mask", ""),
            row.get(f"atac_{evidence}_genetic_evaluable_mask", ""),
        ) + site_coverage + molecule_fields
    if level == "PRIMARY_COUNTS":
        return library + (
            menu,
            int(row.get(f"rna_{evidence}_genetic_evaluable_count", 0) or 0),
            int(row.get(f"atac_{evidence}_genetic_evaluable_count", 0) or 0),
        ) + site_coverage + molecule_fields
    if level == "AUDIT_PROTOTYPE":
        basis = () if evidence == "site" else (
            row.get("rna_molecule_evidence_basis_class", "UNAVAILABLE"),
            row.get("atac_molecule_evidence_basis_class", "UNAVAILABLE"),
        )
        return library + (menu,) + site_coverage + basis
    raise ValueError(level)


def apply_calibration_sensitivities(cells, frozen_keys, excluded_reference_keys=None):
    reference_exclusions = excluded_reference_keys or frozen_keys
    external_scheme = f"EXTERNAL_LIBRARY_{CALIBRATION_LIBRARY}_TRANSFER"
    global_scheme = \
        f"GLOBAL_LIBRARY_{CALIBRATION_LIBRARY}_UNSTRATIFIED"
    schemes = (
        "WITHIN_LIBRARY_WHOLE", "WITHIN_LIBRARY_EXACT_MENU",
        "WITHIN_LIBRARY_EXACT_MENU_SITE_COVERAGE",
        "PRIMARY_EXACT_MASK_AND_COVERAGE",
        "LOO_AUDIT_PROTOTYPE", external_scheme, global_scheme,
    )
    maps = {scheme: defaultdict(list) for scheme in schemes}
    legacy_within = defaultdict(list)
    legacy_external = defaultdict(list)
    row_by_id = {target_id(row): row for row in cells}
    for row in cells:
        cell_id = target_id(row)
        for evidence in ("site", "molecule"):
            for modality in ("rna", "atac"):
                value, winner, candidate_id = _genetic_maximum(
                    row, modality, evidence)
                row[f"{modality}_{evidence}_genetic_maximum_delta"] = value
                row[f"{modality}_{evidence}_genetic_maximum_winner"] = winner
                row[f"{modality}_{evidence}_genetic_maximum_candidate_id"] = \
                    candidate_id
                if not math.isfinite(value):
                    continue
                entries = (
                    ("WITHIN_LIBRARY_WHOLE", _calibration_key(
                        row, evidence, "WHOLE_LIBRARY")),
                    ("WITHIN_LIBRARY_EXACT_MENU", _calibration_key(
                        row, evidence, "EXACT_MENU")),
                    ("WITHIN_LIBRARY_EXACT_MENU_SITE_COVERAGE", _calibration_key(
                        row, evidence, "MENU_SITE_COVERAGE")),
                    ("PRIMARY_EXACT_MASK_AND_COVERAGE", _calibration_key(
                        row, evidence, "PRIMARY_EXACT")),
                    ("PRIMARY_COUNTS", _calibration_key(
                        row, evidence, "PRIMARY_COUNTS")),
                    ("LOO_AUDIT_PROTOTYPE", _calibration_key(
                        row, evidence, "AUDIT_PROTOTYPE")),
                    (external_scheme, _calibration_key(
                        row, evidence, "PRIMARY_EXACT", external=True)),
                    (global_scheme, ()),
                )
                for scheme, key in entries:
                    if scheme in {external_scheme, global_scheme} and \
                            row["library"] != f"lib{CALIBRATION_LIBRARY}":
                        continue
                    maps.setdefault(scheme, defaultdict(list))[
                        (evidence, modality, key)].append((cell_id, value))
        for modality in ("rna", "atac"):
            legacy_value, legacy_winner, legacy_candidate = \
                _unrestricted_maximum(row, modality, "site")
            row[f"{modality}_legacy_unrestricted_maximum_delta"] = legacy_value
            row[f"{modality}_legacy_unrestricted_maximum_winner"] = legacy_winner
            row[f"{modality}_legacy_unrestricted_maximum_candidate_id"] = \
                legacy_candidate
            if math.isfinite(legacy_value):
                legacy_within[(row["library"], modality)].append(
                    (cell_id, legacy_value))
                if row["library"] == f"lib{CALIBRATION_LIBRARY}":
                    legacy_external[modality].append((cell_id, legacy_value))

    for row in cells:
        indexed = target_id(row)
        for modality in ("rna", "atac"):
            target_value = finite(row.get(
                f"{modality}_legacy_unrestricted_maximum_delta", ""))
            within_reference = sorted(
                value for cell_id, value in legacy_within.get(
                    (row["library"], modality), []) if cell_id != indexed)
            external_reference = sorted(
                value for _cell_id, value in legacy_external.get(modality, []))
            row[f"{modality}_simple_within_library_percentile"] = percentile(
                target_value, within_reference)
            row[f"{modality}_global_calibration_library_unstratified_percentile"] = \
                percentile(target_value, external_reference)
        row["simple_within_library_dual_p95"] = all(
            finite(row.get(
                f"{modality}_simple_within_library_percentile", "")) >= 0.95
            for modality in ("rna", "atac"))
        row["global_calibration_library_unstratified_dual_p95"] = all(
            finite(row.get(
                f"{modality}_global_calibration_library_unstratified_percentile",
                "")) >= 0.95
            for modality in ("rna", "atac"))
        row["within_library_dual_p95"] = row[
            "simple_within_library_dual_p95"]

        # Materialize the prespecified primary calibration for every cell so
        # subsequently selected negative comparisons carry the same
        # exact-stratum/coverage evidence as the frozen targets.  Reference
        # membership is fixed before comparison matching and excludes the
        # entire frozen discovery set plus the indexed non-target cell.
        excluded_ids = {
            f"{library_name}:{frozen_barcode}"
            for library_name, frozen_barcode in reference_exclusions
        } | {indexed}
        for evidence in ("site", "molecule"):
            opportunity_available = all(
                str(row.get(
                    f"{assay}_candidate_independent_molecule_opportunity_quintile",
                    "UNAVAILABLE")) != "UNAVAILABLE"
                for assay in ("rna", "atac"))
            primary_pass = {}
            primary_contributors = {}
            for modality in ("rna", "atac"):
                target_ranked = _delta_rank(
                    row.get("_candidates", []), modality, evidence,
                    lambda _candidate: True)
                target_value = target_ranked[0][1]["delta"] \
                    if target_ranked else math.nan
                runner_value = target_ranked[1][1]["delta"] \
                    if len(target_ranked) > 1 else math.nan
                exact_key = _calibration_key(row, evidence, "PRIMARY_EXACT")
                reference = [item for item in maps[
                    "PRIMARY_EXACT_MASK_AND_COVERAGE"].get(
                        (evidence, modality, exact_key), [])
                    if item[0] not in excluded_ids]
                coarsening = "NONE"
                structural_unavailable = evidence == "molecule" and \
                    not opportunity_available
                available = math.isfinite(target_value) and \
                    len(reference) >= MIN_CALIBRATION_REFERENCE and \
                    not structural_unavailable
                reference_values = sorted(item[1] for item in reference)
                winner_count = bisect.bisect_right(reference_values, target_value) \
                    if available else None
                runner_count = bisect.bisect_right(reference_values, runner_value) \
                    if available and math.isfinite(runner_value) else None
                calibrated = winner_count / len(reference_values) \
                    if winner_count is not None else math.nan
                runner_cdf = runner_count / len(reference_values) \
                    if runner_count is not None else math.nan
                p95_pass = bool(available and
                    20 * winner_count >= 19 * len(reference_values))
                margin_pass = bool(available and runner_count is not None and
                    100 * (winner_count - runner_count) > len(reference_values))
                prefix = f"{modality}_{evidence}_primary_exact_mask_and_coverage"
                row[f"{prefix}_percentile"] = calibrated
                row[f"{prefix}_runner_percentile"] = runner_cdf
                row[f"{prefix}_winner_count"] = winner_count
                row[f"{prefix}_runner_count"] = runner_count
                row[f"{prefix}_winner_delta"] = target_value
                row[f"{prefix}_runner_delta"] = runner_value
                row[f"{prefix}_raw_delta_margin"] = target_value - runner_value \
                    if all(math.isfinite(value) for value in
                           (target_value, runner_value)) else math.nan
                row[f"{prefix}_primary_winner_p95"] = p95_pass
                row[f"{prefix}_reference_cdf_margin_pass"] = margin_pass
                row[f"{prefix}_reference_count"] = len(reference)
                row[f"{prefix}_coarsening"] = coarsening
                row[f"{prefix}_status"] = "AVAILABLE" if available else \
                    "UNAVAILABLE_MOLECULE_OPPORTUNITY" \
                    if structural_unavailable else \
                    "UNAVAILABLE_REFERENCE_LT20" \
                    if len(reference) < MIN_CALIBRATION_REFERENCE else \
                    "UNAVAILABLE_TARGET_SCORE"
                primary_pass[modality] = p95_pass
                primary_contributors[modality] = clean(row.get(
                    f"{modality}_{evidence}_genetic_maximum_winner", ""))
            contributor_concordant = bool(primary_contributors.get("rna")) and \
                primary_contributors.get("rna") == primary_contributors.get("atac")
            row[f"primary_{evidence}_coverage_calibrated_contributor"] = \
                primary_contributors.get("rna", "") if contributor_concordant else ""
            row[f"primary_{evidence}_coverage_calibrated_contributor_concordant"] = \
                contributor_concordant
            row[f"primary_{evidence}_coverage_calibrated_dual_p95"] = \
                contributor_concordant and all(
                    primary_pass.get(modality, False)
                    for modality in ("rna", "atac"))

    frozen_rows = []
    for library, barcode in sorted(frozen_keys):
        row = row_by_id.get(f"{library}:{barcode}")
        if row is None:
            raise RuntimeError(f"frozen target absent from completed ledger: {library}:{barcode}")
        frozen_rows.append(row)

    audit_rows = []
    roster_rows = []
    primary_status = defaultdict(dict)
    for row in frozen_rows:
        indexed = target_id(row)
        for evidence in ("site", "molecule"):
            molecule_opportunity_available = all(
                str(row.get(
                    f"{modality}_candidate_independent_molecule_opportunity_quintile",
                    "UNAVAILABLE")) != "UNAVAILABLE"
                for modality in ("rna", "atac"))
            for modality in ("rna", "atac"):
                target_ranked = _delta_rank(
                    row.get("_candidates", []), modality, evidence,
                    lambda _candidate: True)
                target_value = target_ranked[0][1]["delta"] \
                    if target_ranked else math.nan
                runner_value = target_ranked[1][1]["delta"] \
                    if len(target_ranked) > 1 else math.nan
                eligible_runner_set = [item[0]["candidate_id"]
                                       for item in target_ranked[1:]]
                requested = [
                    ("WITHIN_LIBRARY_WHOLE", "WHOLE_LIBRARY", False),
                    ("WITHIN_LIBRARY_EXACT_MENU", "EXACT_MENU", False),
                    ("WITHIN_LIBRARY_EXACT_MENU_SITE_COVERAGE",
                     "MENU_SITE_COVERAGE", False),
                    ("PRIMARY_EXACT_MASK_AND_COVERAGE", "PRIMARY_EXACT", False),
                    ("LOO_AUDIT_PROTOTYPE", "AUDIT_PROTOTYPE", False),
                    (external_scheme, "PRIMARY_EXACT", True),
                    (global_scheme, "WHOLE_LIBRARY", True),
                ]
                for scheme, level, external in requested:
                    if external:
                        key = () if scheme == global_scheme \
                            else _calibration_key(
                                row, evidence, level, external=True)
                    else:
                        key = _calibration_key(row, evidence, level)
                    reference = list(maps[scheme].get(
                        (evidence, modality, key), []))
                    exclusion = "INDEXED_TARGET_ONLY" \
                        if scheme == "LOO_AUDIT_PROTOTYPE" else \
                        "FROZEN_TARGETS_COMPARISONS_AND_SELECTED_CONTROL_PARENTS"
                    excluded_ids = {indexed} if exclusion == \
                        "INDEXED_TARGET_ONLY" else {
                            f"{library_name}:{frozen_barcode}"
                            for library_name, frozen_barcode in reference_exclusions}
                    reference = [item for item in reference
                                 if item[0] not in excluded_ids]
                    coarsening = "NONE"
                    structural_unavailable = (
                        evidence == "molecule" and
                        scheme == "PRIMARY_EXACT_MASK_AND_COVERAGE" and
                        not molecule_opportunity_available)
                    available = math.isfinite(target_value) and \
                        len(reference) >= MIN_CALIBRATION_REFERENCE and \
                        not structural_unavailable
                    values = sorted(item[1] for item in reference)
                    value_percentile = percentile(target_value, values) \
                        if available else math.nan
                    winner_count = bisect.bisect_right(values, target_value) \
                        if available else None
                    runner_count = bisect.bisect_right(values, runner_value) \
                        if available and math.isfinite(runner_value) else None
                    runner_percentile = runner_count / len(values) \
                        if runner_count is not None else math.nan
                    primary_winner_p95 = bool(available and
                        20 * winner_count >= 19 * len(values))
                    cdf_margin_pass = bool(available and runner_count is not None and
                        100 * (winner_count - runner_count) > len(values))
                    status = "AVAILABLE" if available else \
                        "UNAVAILABLE_MOLECULE_OPPORTUNITY" \
                        if structural_unavailable else \
                        "UNAVAILABLE_REFERENCE_LT20" \
                        if len(reference) < MIN_CALIBRATION_REFERENCE else \
                        "UNAVAILABLE_TARGET_SCORE"
                    audit_rows.append({
                        "target_id": indexed, "library": row["library"],
                        "barcode": row["barcode"], "evidence_channel": evidence.upper(),
                        "modality": modality.upper(), "scheme": scheme,
                        "target_statistic": (
                            "maximum score improvement among physical, new-donor, "
                            "genotype-distinguishable usable candidates"),
                        "target_value": target_value,
                        "winner_delta": target_value,
                        "runner_delta": runner_value,
                        "raw_delta_margin": target_value - runner_value
                            if all(math.isfinite(value) for value in
                                   (target_value, runner_value)) else math.nan,
                        "winner_count": winner_count,
                        "runner_count": runner_count,
                        "runner_empirical_percentile": runner_percentile,
                        "reference_cdf_margin": value_percentile - runner_percentile
                            if all(math.isfinite(value) for value in
                                   (value_percentile, runner_percentile)) else math.nan,
                        "primary_winner_p95": primary_winner_p95,
                        "reference_cdf_margin_pass": cdf_margin_pass,
                        "eligible_runner_set": json.dumps(eligible_runner_set),
                        "tie_rule": "DESCENDING_CANDIDATE_ID",
                        "evaluated_contributor": row.get(
                            f"{modality}_{evidence}_genetic_maximum_winner", ""),
                        "reference_exclusion": exclusion,
                        "reference_count": len(reference),
                        "reference_stratum": json.dumps(key, default=str),
                        "coarsening": coarsening,
                        "status": status,
                        "empirical_percentile": value_percentile,
                        "passes_p95": primary_winner_p95,
                        "rna_evaluable_mask_size": row.get(
                            f"rna_{evidence}_genetic_evaluable_count", 0),
                        "atac_evaluable_mask_size": row.get(
                            f"atac_{evidence}_genetic_evaluable_count", 0),
                        "molecule_coverage_conditioning": (
                            "SITE_COVERAGE_AND_MOLECULE_BASIS_ONLY"
                            if scheme == "LOO_AUDIT_PROTOTYPE" and
                            evidence == "molecule" else
                            "MOLECULE_OPPORTUNITY_CONDITIONED"
                            if evidence == "molecule" and available else
                            "NOT_APPLICABLE"),
                    })
                    if scheme in {"PRIMARY_EXACT_MASK_AND_COVERAGE",
                                  "LOO_AUDIT_PROTOTYPE",
                                  external_scheme}:
                        for reference_id, reference_value in reference:
                            roster_rows.append({
                                "target_id": indexed, "library": row["library"],
                                "barcode": row["barcode"],
                                "evidence_channel": evidence.upper(),
                                "modality": modality.upper(), "scheme": scheme,
                                "reference_cell_id": reference_id,
                                "reference_value": reference_value,
                                "reference_stratum": json.dumps(key, default=str),
                                "coarsening": coarsening,
                            })
                    field = scheme.lower()
                    row[f"{modality}_{evidence}_{field}_percentile"] = \
                        value_percentile
                    row[f"{modality}_{evidence}_{field}_reference_count"] = \
                        len(reference)
                    row[f"{modality}_{evidence}_{field}_status"] = status
                    if scheme == "PRIMARY_EXACT_MASK_AND_COVERAGE":
                        primary_status[(indexed, evidence)][modality] = \
                            (primary_winner_p95,
                             clean(row.get(
                                 f"{modality}_{evidence}_genetic_maximum_winner", "")))

        row["within_library_dual_p95"] = all(
            finite(row.get(
                f"{modality}_site_within_library_whole_percentile", "")) >= 0.95
            for modality in ("rna", "atac"))
        for evidence, output_field in (
                ("site", "menu_coverage_calibrated_dual_p95"),
                ("molecule", "primary_molecule_calibrated_dual_p95")):
            states = primary_status[(indexed, evidence)]
            rna_pass, rna_contributor = states.get("rna", (False, ""))
            atac_pass, atac_contributor = states.get("atac", (False, ""))
            row[output_field] = bool(
                rna_pass and atac_pass and rna_contributor and
                rna_contributor == atac_contributor)
            row[f"{output_field}_contributor"] = \
                rna_contributor if row[output_field] else ""
            row[f"{output_field}_contributor_concordant"] = bool(
                rna_contributor and rna_contributor == atac_contributor)

    audit_checks = []
    grouped_audit = {evidence: defaultdict(dict)
                     for evidence in ("site", "molecule")}
    for audit in audit_rows:
        evidence = audit.get("evidence_channel", "").lower()
        if audit.get("scheme") == "LOO_AUDIT_PROTOTYPE" and \
                evidence in grouped_audit:
            grouped_audit[evidence][audit["target_id"]][audit["modality"]] = \
                audit

    # These are immutable historical oracles.  They are deliberately not
    # derived from the rows whose reproduction is being audited.
    expected_members = {
        "site": set(HISTORICAL_SITE_DUAL_P95_CONTRIBUTORS),
        "molecule": set(HISTORICAL_MOLECULE_RAW_DUAL_P95),
    }
    expected_concordant = {
        "site": set(HISTORICAL_SITE_DUAL_P95_CONTRIBUTORS),
        "molecule": set(HISTORICAL_SITE_DUAL_P95_CONTRIBUTORS),
    }
    reproduced_members = {"site": set(), "molecule": set()}
    reproduced_concordant = {"site": set(), "molecule": set()}
    for evidence in ("site", "molecule"):
        contributor_mismatches = []
        for cell_id, assays in grouped_audit[evidence].items():
            rna, atac = assays.get("RNA", {}), assays.get("ATAC", {})
            passes = bool(rna.get("passes_p95") and atac.get("passes_p95"))
            rna_contributor = clean(rna.get("evaluated_contributor", ""))
            atac_contributor = clean(atac.get("evaluated_contributor", ""))
            if passes:
                reproduced_members[evidence].add(cell_id)
            if passes and rna_contributor and rna_contributor == atac_contributor:
                reproduced_concordant[evidence].add(cell_id)
            expected_contributor = HISTORICAL_SITE_DUAL_P95_CONTRIBUTORS.get(
                cell_id)
            if expected_contributor and (rna_contributor != expected_contributor or
                                         atac_contributor != expected_contributor):
                contributor_mismatches.append(
                    f"{cell_id}:expected={expected_contributor}:"
                    f"rna={rna_contributor or 'MISSING'}:"
                    f"atac={atac_contributor or 'MISSING'}")
        expected = expected_members[evidence]
        observed = reproduced_members[evidence]
        fixed_count = 9 if evidence == "site" else 10
        audit_checks.append({
            "scheme": "LOO_AUDIT_PROTOTYPE",
            "evidence_channel": evidence.upper(),
            "dual_p95_targets": len(observed),
            "contributor_concordant_targets": len(
                reproduced_concordant[evidence]),
            "expected_audit_count": fixed_count,
            "expected_target_ids": ",".join(sorted(expected)),
            "target_ids": ",".join(sorted(observed)),
            "contributor_concordant_target_ids": ",".join(sorted(
                reproduced_concordant[evidence])),
            "missing_expected_target_ids": ",".join(sorted(expected - observed)),
            "unexpected_target_ids": ",".join(sorted(observed - expected)),
            "contributor_mapping_mismatches": ";".join(
                sorted(contributor_mismatches)),
            "per_library_composition": json.dumps(Counter(
                value.split(":", 1)[0] for value in observed), sort_keys=True),
            "reproduction_status": "REPRODUCED_EXACT_IDENTITIES"
                if observed == expected and
                reproduced_concordant[evidence] ==
                expected_concordant[evidence] and not contributor_mismatches and
                len(observed) == fixed_count
                else "DIFFERENCE_REQUIRES_CELL_LEVEL_REVIEW",
        })

    def joint_member_from_assays(site_assays, molecule_assays):
        site_rna, site_atac = site_assays.get("RNA", {}), \
            site_assays.get("ATAC", {})
        molecule_rna, molecule_atac = molecule_assays.get("RNA", {}), \
            molecule_assays.get("ATAC", {})
        contributors = [clean(item.get("evaluated_contributor", "")) for item in
                        (site_rna, site_atac, molecule_rna, molecule_atac)]
        return bool(site_rna.get("passes_p95") and
                    site_atac.get("passes_p95") and
                    molecule_rna.get("passes_p95") and
                    molecule_atac.get("passes_p95") and contributors[0] and
                    len(set(contributors)) == 1)

    observed_overlap = {
        cell_id for cell_id in set(grouped_audit["site"]) |
        set(grouped_audit["molecule"])
        if joint_member_from_assays(grouped_audit["site"].get(cell_id, {}),
                                   grouped_audit["molecule"].get(cell_id, {}))}
    expected_overlap = set(HISTORICAL_SITE_DUAL_P95_CONTRIBUTORS)
    overlap_composition = Counter(
        value.split(":", 1)[0] for value in observed_overlap)
    audit_checks.append({
        "scheme": "LOO_AUDIT_PROTOTYPE",
        "evidence_channel": "SITE_AND_MOLECULE",
        "dual_p95_targets": len(observed_overlap),
        "expected_audit_count": 9,
        "expected_target_ids": ",".join(sorted(expected_overlap)),
        "target_ids": ",".join(sorted(observed_overlap)),
        "missing_expected_target_ids": ",".join(
            sorted(expected_overlap - observed_overlap)),
        "unexpected_target_ids": ",".join(
            sorted(observed_overlap - expected_overlap)),
        "per_library_composition": json.dumps(
            overlap_composition, sort_keys=True),
        "reproduction_status": "REPRODUCED_EXACT_IDENTITIES"
            if observed_overlap == expected_overlap and
            len(observed_overlap) == 9 and overlap_composition == Counter(
                {"lib12": 4, "lib20": 4, "lib29": 1}) else
            "DIFFERENCE_REQUIRES_CELL_LEVEL_REVIEW",
    })
    for row in frozen_rows:
        cell_id = target_id(row)
        site = grouped_audit["site"].get(cell_id, {})
        molecule = grouped_audit["molecule"].get(cell_id, {})
        site_rna, site_atac = site.get("RNA", {}), site.get("ATAC", {})
        molecule_rna, molecule_atac = molecule.get("RNA", {}), \
            molecule.get("ATAC", {})
        site_rna_winner = clean(site_rna.get("evaluated_contributor", ""))
        site_atac_winner = clean(site_atac.get("evaluated_contributor", ""))
        molecule_rna_winner = clean(molecule_rna.get(
            "evaluated_contributor", ""))
        molecule_atac_winner = clean(molecule_atac.get(
            "evaluated_contributor", ""))
        audit_checks.append({
            "scheme": "LOO_AUDIT_IDENTITY_RECONCILIATION",
            "evidence_channel": "CELL_LEVEL",
            "target_id": cell_id, "library": row["library"],
            "barcode": row["barcode"],
            "site_dual_p95": bool(site_rna.get("passes_p95") and
                                  site_atac.get("passes_p95")),
            "molecule_dual_p95": bool(molecule_rna.get("passes_p95") and
                                      molecule_atac.get("passes_p95")),
            "site_rna_winner": site_rna_winner,
            "site_atac_winner": site_atac_winner,
            "molecule_rna_winner": molecule_rna_winner,
            "molecule_atac_winner": molecule_atac_winner,
            "site_contributor_agreement": bool(
                site_rna_winner and site_rna_winner == site_atac_winner),
            "molecule_contributor_agreement": bool(
                molecule_rna_winner and
                molecule_rna_winner == molecule_atac_winner),
            "molecule_equals_site_contributor": bool(
                site_rna_winner and len({site_rna_winner, site_atac_winner,
                    molecule_rna_winner, molecule_atac_winner}) == 1),
            "joint_preservation_member": cell_id in observed_overlap,
            "strict_five_member": cell_id in STRICT_SITE_TARGETS,
            "legacy_exact_seven_member": cell_id in
                LEGACY_EXACT_MOLECULE_TARGETS,
            "overlap_three_member": cell_id in (
                STRICT_SITE_TARGETS & LEGACY_EXACT_MOLECULE_TARGETS),
            "reproduction_status": "CELL_LEVEL_RECONCILED",
        })
    def literal_roster_check(label, expected_mapping, evidence):
        missing = []
        contributor_mismatch = []
        locked_mismatch = []
        for cell_id, (expected_locked, expected_contributor) in \
                sorted(expected_mapping.items()):
            current = row_by_id.get(cell_id)
            if current is None:
                missing.append(cell_id)
                continue
            actual_locked = clean(current.get("reconciled_identity_locked", ""))
            if actual_locked != expected_locked:
                locked_mismatch.append(
                    f"{cell_id}:expected={expected_locked}:actual={actual_locked}")
            actual = [clean(current.get(
                f"{modality}_{evidence}_genetic_maximum_winner", ""))
                for modality in ("rna", "atac")]
            if actual != [expected_contributor, expected_contributor]:
                contributor_mismatch.append(
                    f"{cell_id}:expected={expected_contributor}:"
                    f"rna={actual[0] or 'MISSING'}:atac={actual[1] or 'MISSING'}")
        passed = not (missing or locked_mismatch or contributor_mismatch)
        return {
            "scheme": "FROZEN_MEMBERSHIP_AUDIT",
            "evidence_channel": label,
            "expected_audit_count": len(expected_mapping),
            "dual_p95_targets": len(expected_mapping) - len(missing),
            "expected_target_ids": ",".join(sorted(expected_mapping)),
            "target_ids": ",".join(sorted(
                set(expected_mapping) - set(missing))),
            "missing_expected_target_ids": ",".join(missing),
            "unexpected_target_ids": "",
            "locked_mapping_mismatches": ";".join(locked_mismatch),
            "contributor_mapping_mismatches": ";".join(contributor_mismatch),
            "reproduction_status": "REPRODUCED_EXACT_IDENTITIES_AND_CONTRIBUTORS"
                if passed else "DIFFERENCE_REQUIRES_CELL_LEVEL_REVIEW",
        }

    audit_checks.append(literal_roster_check(
        "STRICT_FIVE", STRICT_SITE_TARGET_CONTRIBUTORS, "site"))
    audit_checks.append(literal_roster_check(
        "EXACT_SEVEN", LEGACY_EXACT_MOLECULE_CONTRIBUTORS, "molecule"))
    strict_exact_overlap = set(STRICT_SITE_TARGET_CONTRIBUTORS) & \
        set(LEGACY_EXACT_MOLECULE_CONTRIBUTORS)
    expected_strict_exact_overlap = {
        "lib12:GCCTAATAGCATTTCT", "lib12:TGTTACTTCAAGGACA",
        "lib20:CCTGAGTCATGTTGCA",
    }
    audit_checks.append({
        "scheme": "FROZEN_MEMBERSHIP_AUDIT",
        "evidence_channel": "OVERLAP_THREE",
        "expected_audit_count": 3,
        "dual_p95_targets": len(strict_exact_overlap),
        "expected_target_ids": ",".join(sorted(expected_strict_exact_overlap)),
        "target_ids": ",".join(sorted(strict_exact_overlap)),
        "missing_expected_target_ids": ",".join(sorted(
            expected_strict_exact_overlap - strict_exact_overlap)),
        "unexpected_target_ids": ",".join(sorted(
            strict_exact_overlap - expected_strict_exact_overlap)),
        "reproduction_status": "REPRODUCED_EXACT_IDENTITIES"
            if strict_exact_overlap == expected_strict_exact_overlap else
            "DIFFERENCE_REQUIRES_REVIEW",
    })

    changed = []
    for row in cells:
        changes = []
        old_dual = truthy(row.get("original_site_both_assays_p95", ""))
        genetic_dual = truthy(row.get("genetic_source_both_assays_p95", ""))
        molecule_dual = truthy(row.get("molecule_genetic_both_assays_p95", ""))
        if old_dual != genetic_dual:
            changes.append("UNRESTRICTED_TO_GENETIC_SOURCE_P95_STATUS")
        if genetic_dual != molecule_dual:
            changes.append("SITE_TO_MOLECULE_GENETIC_P95_STATUS")
        if (row["library"], row["barcode"]) in frozen_keys:
            if old_dual != truthy(row.get("menu_coverage_calibrated_dual_p95", "")):
                changes.append("LEGACY_TO_PRIMARY_EXACT_CALIBRATION_STATUS")
            if truthy(row.get("site_to_molecule_category_changed", "")):
                changes.append("SITE_TO_MOLECULE_CATEGORY")
        if changes:
            changed.append({
                "schema_version": SCHEMA, "library": row["library"],
                "barcode": row["barcode"], "change_reasons": ",".join(changes),
                "legacy_library25_dual_p95": old_dual,
                "genetic_source_dual_p95": genetic_dual,
                "molecule_genetic_dual_p95": molecule_dual,
                "primary_site_calibrated_dual_p95": row.get(
                    "menu_coverage_calibrated_dual_p95", ""),
                "primary_molecule_calibrated_dual_p95": row.get(
                    "primary_molecule_calibrated_dual_p95", ""),
                "site_category": row.get("site_existing_evidence_category", ""),
                "molecule_category": row.get(
                    "molecule_existing_evidence_category", ""),
            })
    return changed, audit_rows, roster_rows, audit_checks


def outcome_value(row, name):
    if name == "DUAL_P95_ENTRY":
        return truthy(row.get("original_site_both_assays_p95", ""))
    if name == "SITE_DISCORDANCE":
        return truthy(row.get("original_site_both_assays_p95", "")) and not \
            truthy(row.get("original_site_exact_winner_agreement", ""))
    if name == "GENETIC_SOURCE_AGREEMENT":
        return truthy(row.get("genetic_source_exact_winner_agreement", ""))
    if name == "MIXTURE_COMPATIBLE":
        return row.get("site_existing_evidence_category") == \
            "mixture-compatible in both assays"
    if name == "REPLACEMENT_LIKE":
        return row.get("site_existing_evidence_category") == \
            "replacement-like or upper-boundary fit"
    if name in {"WITHIN_LIBRARY_ONLY_P95_INFLATION", "LIBRARY25_POSITIVE_WITHIN_LIBRARY_NEGATIVE"}:
        return truthy(row.get("global_calibration_library_unstratified_dual_p95",
                              row.get("original_site_both_assays_p95", ""))) and not truthy(
            row.get("simple_within_library_dual_p95", row.get("within_library_dual_p95", "")))
    if name == "WITHIN_LIBRARY_POSITIVE_LIBRARY25_NEGATIVE":
        return truthy(row.get("simple_within_library_dual_p95", "")) and not truthy(
            row.get("global_calibration_library_unstratified_dual_p95", ""))
    raise ValueError(name)


COVERAGE_COVARIATES = (
    "rna_total_counts", "rna_detected_features", "atac_fragments",
    "atac_fragment_records", "atac_cut_sites",
    "rna_candidate_independent_log1p_sites",
    "rna_candidate_independent_log1p_depth",
    "atac_candidate_independent_log1p_sites",
    "atac_candidate_independent_log1p_depth",
    "rna_candidate_independent_log1p_molecule_units",
    "atac_candidate_independent_log1p_molecule_units",
)


def coverage_vector(row):
    values = []
    for field in COVERAGE_COVARIATES:
        value = finite(row.get(field, ""))
        if field in {"rna_total_counts", "rna_detected_features", "atac_fragments",
                     "atac_fragment_records", "atac_cut_sites"} and math.isfinite(value):
            value = math.log1p(max(value, 0.0))
        values.append(value)
    return values


def standardized_difference(left, right, weights_left=None, weights_right=None):
    def moments(values, weights):
        pairs = [(value, 1.0 if weights is None else weights[index])
                 for index, value in enumerate(values)
                 if math.isfinite(value) and
                 (weights is None or (math.isfinite(weights[index]) and
                                      weights[index] > 0))]
        total = sum(weight for _, weight in pairs)
        if total <= 0:
            return math.nan, math.nan
        average = sum(value * weight for value, weight in pairs) / total
        variance = sum(weight * (value - average) ** 2
                       for value, weight in pairs) / total
        return average, variance
    left_mean, left_var = moments(left, weights_left)
    right_mean, right_var = moments(right, weights_right)
    if not all(math.isfinite(value) for value in
               (left_mean, right_mean, left_var, right_var)):
        return math.nan
    pooled_variance = (left_var + right_var) / 2.0
    if pooled_variance <= 1e-15:
        return 0.0 if abs(left_mean - right_mean) <= 1e-15 else math.nan
    return (left_mean - right_mean) / math.sqrt(pooled_variance)


def coverage_exposure(row):
    quintiles = [row.get(f"{modality}_candidate_independent_coverage_quintile")
                 for modality in ("rna", "atac")]
    try:
        return min(int(value) for value in quintiles) <= 2
    except (TypeError, ValueError):
        return None


def coverage_stratum(row):
    return (
        row.get("library", ""),
        clean(row.get("reconciled_identity_locked", "")) or "UNAVAILABLE",
        clean(row.get("ploidy_field", "") or row.get("ploidy_call", "")) or
            "UNAVAILABLE",
        clean(row.get("biological_state", "") or
              row.get("biological_tetraploidy_field", "")) or "UNAVAILABLE",
    )


def effective_sample_size(weights):
    if len(weights) == 0 or any(not math.isfinite(value) or value < 0
                                for value in weights):
        return math.nan
    total = sum(weights)
    squares = sum(value * value for value in weights)
    return total * total / squares if squares > 0 else math.nan


def coverage_adjustment_rows(cells, frozen_keys=None, calibration_audit=None):
    """Coverage-conditional target support and descriptive associations.

    Coverage is not dichotomized into an exposure that is then forced to
    balance on itself.  All rows are explicitly descriptive except the
    prespecified target calibration inherited from the calibration audit.
    """
    frozen_keys = frozen_keys or set()
    calibration_audit = calibration_audit or []
    result_rows = [{"analysis_type": "HISTORICAL_LOW_HIGH_COVERAGE_ADJUSTMENT",
                    "status": "UNRESOLVED_DESCRIPTIVE_ONLY",
                    "reason": "NO_VALIDATED_POST_WEIGHT_BALANCE_OR_COMMON_SUPPORT",
                    "interpretation": "Historical broad matching is not evidence of coverage adjustment; frozen target/comparison matching is separate"}]
    for row in calibration_audit:
        if row.get("scheme") != "PRIMARY_EXACT_MASK_AND_COVERAGE":
            continue
        result_rows.append({
            "analysis_type": "TARGET_COVERAGE_CONDITIONAL_SUPPORT",
            "endpoint": endpoint_label(
                "SITE_GENETIC_SOURCE" if row["evidence_channel"] == "SITE"
                else "MOLECULE_GENETIC_SOURCE"),
            "scope": "FROZEN_TARGET",
            **row,
            "interpretation": (
                "support relative to same-library, pre-outcome coverage and "
                "candidate-opportunity stratum; not a causal coverage adjustment"),
        })

    primary = [row for row in cells
               if int(row["library"].removeprefix("lib")) in PRIMARY_LIBRARIES]
    complete = [(row, coverage_exposure(row), coverage_vector(row)) for row in primary]
    complete = [(row,label,values) for row,label,values in complete
                if label is not None and all(math.isfinite(x) for x in values)]
    if complete and len({label for _,label,_ in complete}) == 2:
        values = np.asarray([x for _,_,x in complete], dtype=float)
        labels = np.asarray([label for _,label,_ in complete], dtype=float)
        scale = values.std(axis=0);scale[scale<1e-12] = 1.0
        design = np.column_stack([np.ones(len(values)), (values-values.mean(axis=0))/scale])
        beta = np.zeros(design.shape[1]);converged = False
        for iteration in range(100):
            probability = 1/(1+np.exp(-np.clip(design@beta, -35,35)))
            weights = probability*(1-probability)
            penalty = np.eye(design.shape[1])*1e-6;penalty[0,0]=0
            gradient = design.T@(labels-probability)-penalty@beta
            hessian = design.T@(design*weights[:,None])+penalty
            step = np.linalg.lstsq(hessian,gradient,rcond=None)[0]
            beta += step
            if np.max(np.abs(step))<1e-8:
                converged = True;break
        probability = 1/(1+np.exp(-np.clip(design@beta,-35,35)))
        low = labels.astype(bool);high = ~low
        lower = max(float(probability[low].min()),float(probability[high].min()))
        upper = min(float(probability[low].max()),float(probability[high].max()))
        support = (probability>=lower)&(probability<=upper) if lower<=upper else np.zeros(len(values),dtype=bool)
        overlap = np.where(low,1-probability,probability)*support
        for index, field in enumerate(COVERAGE_COVARIATES):
            difference = standardized_difference(values[low,index],values[high,index],overlap[low],overlap[high])
            result_rows.append({"analysis_type": "BROAD_LOW_HIGH_OVERLAP_WEIGHTED_BALANCE", "covariate": field,
                "total_primary_cells": len(primary), "complete_case_cells": len(complete),
                "missing_covariate_attrition": len(primary)-len(complete),
                "common_support_lower": lower, "common_support_upper": upper,
                "common_support_low_cells": int(np.sum(support&low)), "common_support_high_cells": int(np.sum(support&high)),
                "common_support_attrition": int(np.sum(~support)),
                "low_effective_sample_size": effective_sample_size(overlap[low]),
                "high_effective_sample_size": effective_sample_size(overlap[high]),
                "weighted_absolute_standardized_difference": abs(difference),
                "balance_threshold": 0.1, "balance_status": "PASS" if math.isfinite(difference) and abs(difference)<=0.1 else "FAIL",
                "propensity_converged": converged,
                "interpretation": "UNRESOLVED_DESCRIPTIVE_COVERAGE_DEFINES_EXPOSURE_NOT_INDEPENDENT_VALIDATION"})
    else:
        result_rows.append({"analysis_type": "BROAD_LOW_HIGH_OVERLAP_WEIGHTED_BALANCE", "status": "UNAVAILABLE_COMPLETE_CASE_OVERLAP",
                            "complete_case_cells": len(complete), "total_primary_cells": len(primary)})
    strata = defaultdict(list)
    for row in primary:
        key = coverage_stratum(row) + (
            row.get("canonical_menu_signature", ""),
            row.get("rna_molecule_evidence_basis_class", "UNAVAILABLE"),
            row.get("atac_molecule_evidence_basis_class", "UNAVAILABLE"),
        )
        strata[key].append(row)

    outcome_fields = {
        "RNA_GENETIC_MAX_SCORE": lambda row: finite(
            row.get("rna_site_genetic_maximum_delta", "")),
        "ATAC_GENETIC_MAX_SCORE": lambda row: finite(
            row.get("atac_site_genetic_maximum_delta", "")),
        "RNA_WINNER_MARGIN": lambda row: finite(
            row.get("rna_site_genetic_margin", "")),
        "ATAC_WINNER_MARGIN": lambda row: finite(
            row.get("atac_site_genetic_margin", "")),
        "REPLACEMENT_OR_BOUNDARY": lambda row: float(
            row.get("site_existing_evidence_category") ==
            "replacement-like or upper-boundary fit"),
        "GENETIC_WINNER_AGREEMENT": lambda row: float(
            truthy(row.get("genetic_source_exact_winner_agreement", ""))),
    }

    def within_stratum_association(covariate_index, outcome_function):
        x_residual = []
        y_residual = []
        contributing = 0
        for members in strata.values():
            pairs = []
            for row in members:
                x = coverage_vector(row)[covariate_index]
                y = outcome_function(row)
                if math.isfinite(x) and math.isfinite(y):
                    pairs.append((x, y))
            if len(pairs) < 3:
                continue
            contributing += 1
            x_mean = mean([pair[0] for pair in pairs])
            y_mean = mean([pair[1] for pair in pairs])
            x_residual.extend(pair[0] - x_mean for pair in pairs)
            y_residual.extend(pair[1] - y_mean for pair in pairs)
        if len(x_residual) < 3:
            return math.nan, math.nan, len(x_residual), contributing
        cross = sum(x * y for x, y in zip(x_residual, y_residual))
        x_square = sum(x * x for x in x_residual)
        y_square = sum(y * y for y in y_residual)
        correlation = cross / math.sqrt(x_square * y_square) \
            if x_square > 0 and y_square > 0 else math.nan
        slope = cross / x_square if x_square > 0 else math.nan
        return correlation, slope, len(x_residual), contributing

    for covariate_index, covariate in enumerate(COVERAGE_COVARIATES):
        values = [coverage_vector(row)[covariate_index] for row in primary]
        finite_values = [value for value in values if math.isfinite(value)]
        for outcome, function in outcome_fields.items():
            correlation, slope, n_cells, n_strata = \
                within_stratum_association(covariate_index, function)
            result_rows.append({
                "analysis_type": "DESCRIPTIVE_CONTINUOUS_COVERAGE_ASSOCIATION",
                "endpoint": outcome, "scope": "PRIMARY_SIX_LIBRARIES",
                "coverage_covariate": covariate,
                "transformation": (
                    "log1p" if covariate in {
                        "rna_total_counts", "rna_detected_features",
                        "atac_fragments", "atac_fragment_records",
                        "atac_cut_sites"} else "precomputed log1p median"),
                "within_compatible_stratum_correlation": correlation,
                "within_compatible_stratum_slope": slope,
                "contributing_cells": n_cells,
                "contributing_strata": n_strata,
                "coverage_min": min(finite_values) if finite_values else math.nan,
                "coverage_max": max(finite_values) if finite_values else math.nan,
                "interpretation": "DESCRIPTIVE_ONLY_NOT_CAUSAL_ADJUSTMENT",
            })

    target_ids = {f"{library}:{barcode}" for library, barcode in frozen_keys}
    for label, membership in (
            ("ALL_FROZEN_59", target_ids),
            ("STRICT_SITE_FIVE", STRICT_SITE_TARGETS),
            ("LEGACY_EXACT_MOLECULE_SEVEN", LEGACY_EXACT_MOLECULE_TARGETS),
            ("STRICT_EXACT_OVERLAP_THREE",
             STRICT_SITE_TARGETS & LEGACY_EXACT_MOLECULE_TARGETS)):
        relevant = [row for row in result_rows
                    if row.get("analysis_type") ==
                    "TARGET_COVERAGE_CONDITIONAL_SUPPORT" and
                    row.get("target_id") in membership]
        combinations = defaultdict(set)
        contributors = defaultdict(lambda: defaultdict(dict))
        for row in relevant:
            if truthy(row.get("passes_p95", "")):
                combinations[row["target_id"]].add((
                    row["evidence_channel"], row["modality"]))
                contributors[row["target_id"]][row["evidence_channel"]][
                    row["modality"]] = clean(row.get(
                        "evaluated_contributor", ""))
        def preserved(target_id_value, channel):
            assays = contributors[target_id_value][channel]
            return bool(assays.get("RNA")) and \
                assays.get("RNA") == assays.get("ATAC") and \
                {(channel, "RNA"), (channel, "ATAC")}.issubset(
                    combinations[target_id_value])
        site_numerator = sum(preserved(value, "SITE")
                             for value in membership)
        molecule_numerator = sum(preserved(value, "MOLECULE")
                                 for value in membership)
        denominator = len(membership)
        def wilson_interval(successes, total):
            if not total:
                return math.nan, math.nan
            z = 1.959963984540054
            p = successes / total
            scale = 1 + z * z / total
            center = (p + z * z / (2 * total)) / scale
            radius = z * math.sqrt(
                p * (1 - p) / total + z * z / (4 * total * total)) / scale
            return center - radius, center + radius
        site_low, site_high = wilson_interval(site_numerator, denominator)
        molecule_low, molecule_high = wilson_interval(
            molecule_numerator, denominator)
        result_rows.append({
            "analysis_type": "TARGET_COVERAGE_CONDITIONAL_SUBSET_SUMMARY",
            "endpoint": "genetically distinguishable new contributor",
            "scope": label, "targets": len(membership),
            "site_contributor_concordant_numerator": site_numerator,
            "site_denominator": denominator,
            "site_fraction_retained": site_numerator / denominator
                if denominator else math.nan,
            "site_fraction_wilson_low": site_low,
            "site_fraction_wilson_high": site_high,
            "molecule_contributor_concordant_numerator": molecule_numerator,
            "molecule_denominator": denominator,
            "molecule_fraction_retained": molecule_numerator / denominator
                if denominator else math.nan,
            "molecule_fraction_wilson_low": molecule_low,
            "molecule_fraction_wilson_high": molecule_high,
            "preserved_site_library_distribution": json.dumps(Counter(
                value.split(":", 1)[0] for value in membership
                if preserved(value, "SITE")), sort_keys=True),
            "preserved_molecule_library_distribution": json.dumps(Counter(
                value.split(":", 1)[0] for value in membership
                if preserved(value, "MOLECULE")), sort_keys=True),
            "prespecified_preservation_criterion":
                "DUAL_P95_SAME_EVALUATED_CONTRIBUTOR",
            "preservation_status": "UNRESOLVED_PENDING_ALL_COVARIATE_BALANCE_PASS",
            "interpretation": "COVERAGE_CONDITIONAL_SENSITIVITY",
        })
    return result_rows, []


def grouped_evidence_rows(cells):
    rows = []
    group_fields = (
        "library", "reconciled_identity_locked",
        "rna_site_genetic_second_state",
    )
    groups = defaultdict(list)
    for row in cells:
        groups[tuple(clean(row.get(field, "")) for field in group_fields)].append(row)
    for (library, locked, candidate), values in sorted(groups.items()):
        if not candidate:
            continue
        categories = Counter(row["site_existing_evidence_category"] for row in values)
        molecule_categories = Counter(
            row["molecule_existing_evidence_category"] for row in values)
        rows.append({
            "analysis_type": "EXISTING_EVIDENCE_CATEGORY",
            "endpoint": "SITE_AND_MOLECULE",
            "scope": "LIBRARY_LOCKED_IDENTITY_PROPOSED_CONTRIBUTOR",
            "library": library, "locked_identity": locked,
            "proposed_contributor": candidate, "cells": len(values),
            "site_mixture_compatible": categories[
                "mixture-compatible in both assays"],
            "site_addition_uncertain": categories[
                "addition-compatible but uncertain"],
            "site_replacement_or_boundary": categories[
                "replacement-like or upper-boundary fit"],
            "site_rna_atac_conflict": categories["RNA/ATAC conflict"],
            "site_insufficient": categories["insufficient evidence"],
            "molecule_mixture_compatible": molecule_categories[
                "mixture-compatible in both assays"],
            "molecule_addition_uncertain": molecule_categories[
                "addition-compatible but uncertain"],
            "molecule_replacement_or_boundary": molecule_categories[
                "replacement-like or upper-boundary fit"],
            "molecule_rna_atac_conflict": molecule_categories["RNA/ATAC conflict"],
            "molecule_insufficient": molecule_categories["insufficient evidence"],
            "rna_molecule_winner_retained": sum(truthy(row.get(
                "rna_molecule_retains_unrestricted_site_winner", ""))
                for row in values),
            "atac_molecule_winner_retained": sum(truthy(row.get(
                "atac_molecule_retains_unrestricted_site_winner", ""))
                for row in values),
        })
    return rows


def occupancy_descriptive_rows(cells):
    groups = defaultdict(list)
    for row in cells:
        genotype_equivalent_top = any(truthy(row.get(
            f"{modality}_site_unrestricted_genotype_equivalent", ""))
            for modality in ("rna", "atac"))
        occupancy_ranked = clean(row.get("ranking_evidence_basis", "")) == \
            "OCCUPANCY_GENOTYPE_EQUIVALENT"
        if not genotype_equivalent_top and not occupancy_ranked:
            continue
        candidate = clean(row.get("best_second_state", "") or
                          row.get("rna_site_unrestricted_second_state", "") or
                          row.get("atac_site_unrestricted_second_state", ""))
        groups[(row["library"],
                clean(row.get("reconciled_identity_locked", "")),
                candidate)].append(row)
    output = []
    for (library, locked, candidate), rows in sorted(groups.items()):
        occupancy = [finite(row.get(
            "technical_occupancy_library25_empirical_percentile", ""))
            for row in rows]
        output.append({
            "analysis_type": "OCCUPANCY_HIGH_CONTENT_DESCRIPTIVE",
            "endpoint": "GENOTYPE_EQUIVALENT_NOT_GENETIC_DOUBLET_EVIDENCE",
            "scope": "LIBRARY_LOCKED_IDENTITY_CANDIDATE",
            "library": library, "locked_identity": locked,
            "proposed_contributor": candidate, "cells": len(rows),
            "priority_p99_cells": sum(
                clean(row.get("discovery_priority_tier", "")) == "PRIORITY_P99"
                for row in rows),
            "review_p95_cells": sum(
                clean(row.get("discovery_priority_tier", "")) == "REVIEW_P95"
                for row in rows),
            "median_technical_occupancy_library25_percentile":
                median(occupancy),
            "interpretation": (
                "DESCRIPTIVE_HIGH_CONTENT_OR_OCCUPANCY_ONLY;EXCLUDED_FROM_"
                "GENETIC_SOURCE_ENDPOINT"),
        })
    return output


def denominator_rows(cells, score_rows_by_library):
    rows = []
    scopes = (
        ("LEDGER_CELLS", lambda row: True),
        ("SCOREABLE_CELLS", lambda row: int(row.get("menu_size", 0) or 0) > 0),
        ("SITE_EVALUABLE_BOTH_ASSAYS", lambda row:
            bool(clean(row.get("rna_site_unrestricted_second_state", ""))) and
            bool(clean(row.get("atac_site_unrestricted_second_state", "")))),
        ("MOLECULE_EVALUABLE_BOTH_ASSAYS", lambda row:
            bool(clean(row.get("rna_molecule_unrestricted_second_state", ""))) and
            bool(clean(row.get("atac_molecule_unrestricted_second_state", "")))),
        ("ORIGINAL_BOTH_ASSAYS_P95", lambda row:
            truthy(row.get("original_site_both_assays_p95", ""))),
        ("ORIGINAL_EXACT_WINNER_AGREEMENT", lambda row:
            truthy(row.get("original_site_exact_winner_agreement", ""))),
        ("GENETIC_SOURCE_BOTH_ASSAYS_P95", lambda row:
            truthy(row.get("genetic_source_both_assays_p95", ""))),
        ("GENETIC_SOURCE_EXACT_AGREEMENT", lambda row:
            truthy(row.get("genetic_source_exact_winner_agreement", ""))),
        ("SITE_MIXTURE_COMPATIBLE", lambda row:
            row.get("site_existing_evidence_category") ==
            "mixture-compatible in both assays"),
        ("SITE_ADDITION_UNCERTAIN", lambda row:
            row.get("site_existing_evidence_category") ==
            "addition-compatible but uncertain"),
        ("SITE_REPLACEMENT_OR_BOUNDARY", lambda row:
            row.get("site_existing_evidence_category") ==
            "replacement-like or upper-boundary fit"),
    )
    libraries = [f"lib{value}" for value in ALLOWED_LIBRARIES]
    for scope, predicate in scopes:
        for library in libraries + ["PRIMARY_SIX", "ALL_SEVEN"]:
            subset = cells if library == "ALL_SEVEN" else \
                [row for row in cells
                 if (int(row["library"].removeprefix("lib")) in PRIMARY_LIBRARIES
                     if library == "PRIMARY_SIX" else row["library"] == library)]
            rows.append({
                "record_type": "DENOMINATOR", "library": library,
                "assay": "BOTH", "input_role": scope,
                "absolute_path": "DERIVED_FROM_AGGREGATE",
                "exists": True, "bytes": "", "gzip_envelope_status": "",
                "schema_status": "", "unique_key_status": "PASS",
                "count": sum(predicate(row) for row in subset),
                "exclusion_reason": "",
            })
    for library, count in sorted(score_rows_by_library.items()):
        rows.append({
            "record_type": "DENOMINATOR", "library": library,
            "assay": "BOTH", "input_role": "UNIQUE_CANDIDATE_ROWS",
            "absolute_path": "DERIVED_FROM_AGGREGATE", "exists": True,
            "bytes": "", "gzip_envelope_status": "", "schema_status": "",
            "unique_key_status": "PASS", "count": count,
            "exclusion_reason": "",
        })
    return rows


def target_and_match_cells(cells, frozen_keys):
    by_key = {(row["library"], row["barcode"]): row for row in cells}
    missing = sorted(frozen_keys - set(by_key))
    if missing:
        raise RuntimeError(
            "frozen targets absent from completed cells: " +
            ",".join(f"{library}:{barcode}" for library, barcode in missing))
    targets = [by_key[key] for key in sorted(frozen_keys)]
    target_keys = {(row["library"], row["barcode"]) for row in targets}
    library_stats = {}
    for library in (f"lib{value}" for value in TARGET_LIBRARIES):
        rows = [row for row in cells if row["library"] == library and
                all(math.isfinite(value) for value in coverage_vector(row))]
        means = []
        scales = []
        for index in range(len(COVERAGE_COVARIATES)):
            values = [coverage_vector(row)[index] for row in rows]
            average = mean(values)
            variance = mean([(value - average) ** 2 for value in values]) \
                if values else math.nan
            means.append(average)
            scales.append(math.sqrt(variance) if math.isfinite(variance) and
                          variance > 1e-15 else 0.0)
        library_stats[library] = (means, scales)

    def standardized(row):
        means, scales = library_stats[row["library"]]
        return [
            (value - average) / scale if math.isfinite(value) and
            math.isfinite(average) and scale > 0 else
            0.0 if math.isfinite(value) and math.isfinite(average) and
            abs(value - average) <= 1e-15 else math.nan
            for value, average, scale in zip(coverage_vector(row), means, scales)
        ]

    def compatibility(row):
        return (
            row["library"], clean(row.get("reconciled_identity_locked", "")),
            clean(row.get("ploidy_field", "") or row.get("ploidy_call", "")) or
                "UNAVAILABLE",
            clean(row.get("biological_state", "") or
                  row.get("biological_tetraploidy_field", "")) or "UNAVAILABLE",
            row.get("canonical_menu_signature", ""),
            row.get("rna_molecule_evidence_basis_class", "UNAVAILABLE"),
            row.get("atac_molecule_evidence_basis_class", "UNAVAILABLE"),
        )

    by_stratum = defaultdict(list)
    for row in cells:
        key = (row["library"], row["barcode"])
        if int(row["library"].removeprefix("lib")) not in TARGET_LIBRARIES or \
                key in target_keys:
            continue
        exact_agreement = truthy(row.get(
            "genetic_source_exact_winner_agreement", ""))
        addition_call = row.get("site_existing_evidence_category") in {
            "mixture-compatible in both assays",
            "addition-compatible but uncertain",
        }
        if exact_agreement or addition_call:
            continue
        row["_comparison_class"] = (
            "HARD_NEGATIVE_CONFLICTING_DUAL_P95"
            if truthy(row.get("genetic_source_both_assays_p95", "")) else
            "ORDINARY_NEGATIVE")
        by_stratum[compatibility(row)].append(row)
    def hungarian(cost):
        """Minimum-cost rectangular assignment; rows <= columns."""
        n = len(cost)
        m = len(cost[0]) if n else 0
        if n > m:
            raise RuntimeError("matching assignment has fewer columns than rows")
        u = [0.0] * (n + 1)
        v = [0.0] * (m + 1)
        p = [0] * (m + 1)
        way = [0] * (m + 1)
        for i in range(1, n + 1):
            p[0] = i
            minv = [math.inf] * (m + 1)
            used = [False] * (m + 1)
            j0 = 0
            while True:
                used[j0] = True
                i0 = p[j0]
                delta = math.inf
                j1 = 0
                for j in range(1, m + 1):
                    if used[j]:
                        continue
                    current = cost[i0 - 1][j - 1] - u[i0] - v[j]
                    if current < minv[j] - 1e-15:
                        minv[j] = current
                        way[j] = j0
                    if minv[j] < delta - 1e-15 or (
                            abs(minv[j] - delta) <= 1e-15 and j < j1):
                        delta, j1 = minv[j], j
                if not math.isfinite(delta):
                    raise RuntimeError("matching assignment is infeasible")
                for j in range(m + 1):
                    if used[j]:
                        u[p[j]] += delta
                        v[j] -= delta
                    else:
                        minv[j] -= delta
                j0 = j1
                if p[j0] == 0:
                    break
            while True:
                j1 = way[j0]
                p[j0] = p[j1]
                j0 = j1
                if j0 == 0:
                    break
        assignment = [-1] * n
        for j in range(1, m + 1):
            if p[j]:
                assignment[p[j] - 1] = j - 1
        return assignment

    sorted_targets = sorted(targets, key=lambda row: (
        row["library"], row["barcode"]))
    eligible_edges = {}
    target_diagnostics = {}
    comparison_by_key = {}
    for target in sorted_targets:
        t_key = (target["library"], target["barcode"])
        vector = standardized(target)
        choices = sorted(by_stratum.get(compatibility(target), []),
                         key=lambda row: (row["library"], row["barcode"]))
        target_missing = [COVERAGE_COVARIATES[index]
                          for index, value in enumerate(vector)
                          if not math.isfinite(value)]
        complete_choices = [candidate for candidate in choices
                            if all(math.isfinite(value)
                                   for value in standardized(candidate))]
        bounds = []
        if complete_choices:
            raw_vectors = [coverage_vector(candidate)
                           for candidate in complete_choices]
            bounds = [(min(values), max(values)) for values in zip(*raw_vectors)]
        target_raw = coverage_vector(target)
        support_pass = not target_missing and bool(bounds) and all(
            low <= value <= high
            for value, (low, high) in zip(target_raw, bounds))
        edges = []
        if support_pass:
            for candidate in complete_choices:
                c_vector = standardized(candidate)
                absolute = [abs(left - right)
                            for left, right in zip(vector, c_vector)]
                distance = math.sqrt(sum(value * value for value in absolute))
                if max(absolute) <= TARGET_MATCH_MAX_ABS_Z and \
                        distance <= TARGET_MATCH_MAX_DISTANCE:
                    key = (candidate["library"], candidate["barcode"])
                    comparison_by_key[key] = candidate
                    edges.append((key, distance, max(absolute)))
        eligible_edges[t_key] = edges
        target_diagnostics[t_key] = {
            "choices": len(choices), "complete": len(complete_choices),
            "target_missing": target_missing,
            "common_support_pass": support_pass,
            "common_support_bounds": bounds,
        }

    comparison_keys = sorted(comparison_by_key)
    comparison_index = {key: index for index, key in enumerate(comparison_keys)}
    penalty = 1_000_000.0
    invalid = 2_000_000.0
    costs = []
    edge_lookup = {}
    for row_index, target in enumerate(sorted_targets):
        t_key = (target["library"], target["barcode"])
        row = [invalid] * len(comparison_keys) + [penalty] * len(sorted_targets)
        for key, distance, maximum in eligible_edges[t_key]:
            row[comparison_index[key]] = distance
            edge_lookup[(t_key, key)] = (distance, maximum)
        # Each target has its own deterministic unmatched column.
        for dummy in range(len(sorted_targets)):
            row[len(comparison_keys) + dummy] = penalty if dummy == row_index \
                else invalid
        costs.append(row)
    assigned_columns = hungarian(costs)
    assigned = {}
    matches = []
    for row_index, (target, column) in enumerate(zip(
            sorted_targets, assigned_columns)):
        t_key = (target["library"], target["barcode"])
        if 0 <= column < len(comparison_keys) and costs[
                row_index][column] < penalty:
            key = comparison_keys[column]
            assigned[t_key] = key
            matches.append(comparison_by_key[key])

    manifest = []
    for target in sorted_targets:
        t_key = (target["library"], target["barcode"])
        match_key = assigned.get(t_key)
        match = comparison_by_key.get(match_key) if match_key else None
        edge = edge_lookup.get((t_key, match_key), (math.nan, math.nan))
        diagnostic = target_diagnostics[t_key]
        manifest.append({
            "target_id": f"{target['library']}:{target['barcode']}",
            "library": target["library"], "target_barcode": target["barcode"],
            "locked_identity": target.get("reconciled_identity_locked", ""),
            "proposed_contributor": target.get("rna_site_genetic_second_state", ""),
            "site_existing_evidence_category": target[
                "site_existing_evidence_category"],
            "molecule_existing_evidence_category": target[
                "molecule_existing_evidence_category"],
            "matched_comparison_barcode": match["barcode"] if match else "",
            "match_distance": edge[0] if match else math.nan,
            "maximum_absolute_standardized_coverage_difference":
                edge[1] if match else math.nan,
            "matching_stratum": json.dumps(compatibility(target)),
            "eligible_comparisons_in_stratum": diagnostic["choices"],
            "common_support_comparisons": diagnostic["complete"]
                if diagnostic["common_support_pass"] else 0,
            "common_support_pass": diagnostic["common_support_pass"],
            "common_support_bounds": json.dumps(
                diagnostic["common_support_bounds"]),
            "missing_covariate_policy": (
                "COMPLETE_CASE_ALL_PREDECLARED_COVARIATES;NO_PARTIAL_DISTANCE"),
            "target_missing_covariates": ",".join(
                diagnostic["target_missing"]),
            "comparison_missing_covariate_exclusions": diagnostic["choices"] -
                diagnostic["complete"],
            "comparison_class": match.get("_comparison_class", "") if match else "",
            "match_status": "MATCHED" if match else
                "UNMATCHED_TARGET_MISSING_COVARIATE"
                if diagnostic["target_missing"] else
                "UNMATCHED_COMMON_SUPPORT" if not diagnostic[
                    "common_support_pass"] else
                "UNMATCHED_CALIPER_OR_GLOBAL_CAPACITY",
            "coverage_caliper": (
                f"max_abs_z<={TARGET_MATCH_MAX_ABS_Z};"
                f"euclidean_z<={TARGET_MATCH_MAX_DISTANCE}"),
            "selection_rule": (
                "MAX_CARDINALITY_THEN_MINIMUM_TOTAL_DISTANCE_THEN_"
                "LEXICOGRAPHIC_ASSIGNMENT;"
                "NO_EXACT_GENETIC_AGREEMENT;NO_LEGACY_ADDITION_CALL"),
            "strict_site_five": target_id(target) in STRICT_SITE_TARGETS,
            "legacy_dual_molecule_p95_nine": truthy(target.get(
                "molecule_genetic_both_assays_p95", "")),
            "legacy_exact_molecule_supported_seven":
                target_id(target) in LEGACY_EXACT_MOLECULE_TARGETS,
            "strict_exact_overlap_three": target_id(target) in (
                STRICT_SITE_TARGETS & LEGACY_EXACT_MOLECULE_TARGETS),
            "target_rna_molecule_evidence_basis": target.get(
                "rna_molecule_evidence_basis_class", "UNAVAILABLE"),
            "target_atac_molecule_evidence_basis": target.get(
                "atac_molecule_evidence_basis_class", "UNAVAILABLE"),
        })
        proposed_id = clean(target.get("rna_site_genetic_candidate_id", ""))
        proposed = next((candidate for candidate in target.get("_candidates", [])
                         if candidate.get("candidate_id") == proposed_id), {})
        manifest[-1].update({
            "proposed_candidate_id": proposed_id,
            "proposed_locked_copy_vector": proposed.get(
                "locked_copy_vector", ""),
            "proposed_second_copy_vector": proposed.get(
                "second_copy_vector", ""),
            "proposed_relationship": proposed.get("relationship", ""),
            "proposed_physical_pool_state": proposed.get("physical", False),
            "proposed_component_only_state": proposed.get(
                "component_only", False),
            "proposed_new_donor_state": proposed.get("new_donor", False),
            "proposed_genotype_distinguishable": not (
                proposed.get("state_genotype_equivalent", False) or
                proposed.get("genotype_equivalent", False)),
            "proposed_candidate_origin_class": proposed.get(
                "origin_class", ""),
        })
        manifest[-1].update({
            "complete_candidate_menu_signature": target.get(
                "canonical_menu_signature", ""),
            "complete_candidate_menu_size": target.get("menu_size", 0),
        })
        if match is not None:
            comparison_id = clean(match.get(
                "rna_site_genetic_candidate_id", ""))
            comparison_candidate = next((candidate for candidate in
                match.get("_candidates", [])
                if candidate.get("candidate_id") == comparison_id), {})
            manifest[-1].update({
                "comparison_locked_identity": match.get(
                    "reconciled_identity_locked", ""),
                "comparison_proposed_contributor": match.get(
                    "rna_site_genetic_second_state", ""),
                "comparison_proposed_candidate_id": comparison_id,
                "comparison_locked_copy_vector": comparison_candidate.get(
                    "locked_copy_vector", ""),
                "comparison_second_copy_vector": comparison_candidate.get(
                    "second_copy_vector", ""),
                "comparison_relationship": comparison_candidate.get(
                    "relationship", ""),
                "comparison_physical_pool_state": comparison_candidate.get(
                    "physical", False),
                "comparison_component_only_state": comparison_candidate.get(
                    "component_only", False),
                "comparison_new_donor_state": comparison_candidate.get(
                    "new_donor", False),
                "comparison_genotype_distinguishable": not (
                    comparison_candidate.get("state_genotype_equivalent", False)
                    or comparison_candidate.get("genotype_equivalent", False)),
                "comparison_complete_candidate_menu_signature": match.get(
                    "canonical_menu_signature", ""),
                "comparison_complete_candidate_menu_size": match.get(
                    "menu_size", 0),
                "comparison_rna_molecule_evidence_basis": match.get(
                    "rna_molecule_evidence_basis_class", "UNAVAILABLE"),
                "comparison_atac_molecule_evidence_basis": match.get(
                    "atac_molecule_evidence_basis_class", "UNAVAILABLE"),
            })
        for member_name, member in (("target", target), ("comparison", match)):
            if member is None:
                continue
            for modality in ("rna", "atac"):
                for evidence in ("site", "molecule"):
                    source_prefix = \
                        f"{modality}_{evidence}_primary_exact_mask_and_coverage"
                    output_prefix = f"{member_name}_{source_prefix}"
                    for field in ("percentile", "reference_count", "coarsening",
                                  "status"):
                        manifest[-1][f"{output_prefix}_{field}"] = member.get(
                            f"{source_prefix}_{field}", "")
    balance = []
    paired_targets = [by_key[(row["library"], row["target_barcode"])]
                      for row in manifest if row["match_status"] == "MATCHED"]
    paired_matches = [by_key[(row["library"], row["matched_comparison_barcode"])]
                      for row in manifest if row["match_status"] == "MATCHED"]
    comparison_pool = [row for values in by_stratum.values() for row in values]

    def complete(rows):
        return [row for row in rows
                if all(math.isfinite(value) for value in coverage_vector(row))]

    def weighted_moments(values, weights):
        total = sum(weights)
        if total <= 0 or not values:
            return math.nan, math.nan
        average = sum(value * weight for value, weight in zip(values, weights)) / total
        variance = sum(weight * (value - average) ** 2
                       for value, weight in zip(values, weights)) / total
        return average, variance

    def overlap_weights(target_rows, comparison_rows):
        target_complete = complete(target_rows)
        comparison_complete = complete(comparison_rows)
        combined = target_complete + comparison_complete
        labels = np.array([1.0] * len(target_complete) +
                          [0.0] * len(comparison_complete))
        if not combined or len(set(labels.tolist())) < 2:
            return {}, {}, {"status": "UNAVAILABLE", "converged": False}
        matrix = np.array([coverage_vector(row) for row in combined], dtype=float)
        means = matrix.mean(axis=0)
        scales = matrix.std(axis=0)
        scales[scales == 0] = 1.0
        design = np.column_stack((np.ones(len(matrix)), (matrix - means) / scales))
        beta = np.zeros(design.shape[1])
        converged = False
        iterations = 0
        try:
            for iteration in range(1, 51):
                probability = np.clip(1.0 / (1.0 + np.exp(-design @ beta)),
                                      1e-6, 1 - 1e-6)
                weights = probability * (1 - probability)
                penalty_matrix = np.diag(
                    [0.0] + [1e-6] * (design.shape[1] - 1))
                hessian = design.T @ (design * weights[:, None]) + \
                    penalty_matrix
                gradient = design.T @ (labels - probability) - \
                    np.r_[0.0, 1e-6 * beta[1:]]
                step = np.linalg.solve(hessian, gradient)
                beta += step
                iterations = iteration
                if np.max(np.abs(step)) < 1e-9:
                    converged = True
                    break
        except (np.linalg.LinAlgError, FloatingPointError):
            return {}, {}, {"status": "UNAVAILABLE_NUMERIC_FAILURE",
                            "converged": False, "iterations": iterations}
        probability = np.clip(1.0 / (1.0 + np.exp(-design @ beta)),
                              1e-6, 1 - 1e-6)
        target_probability = probability[:len(target_complete)]
        comparison_probability = probability[len(target_complete):]
        lower = max(float(np.min(target_probability)),
                    float(np.min(comparison_probability)))
        upper = min(float(np.max(target_probability)),
                    float(np.max(comparison_probability)))
        if lower > upper:
            return {}, {}, {
                "status": "UNAVAILABLE_NO_PROPENSITY_COMMON_SUPPORT",
                "converged": converged, "iterations": iterations,
                "coefficients": beta.tolist(), "support_lower": lower,
                "support_upper": upper,
            }
        target_weights, comparison_weights = {}, {}
        for row, label, probability_value in zip(combined, labels, probability):
            key = (row["library"], row["barcode"])
            if probability_value < lower or probability_value > upper:
                continue
            if label:
                target_weights[key] = 1.0 - probability_value
            else:
                comparison_weights[key] = probability_value
        return target_weights, comparison_weights, {
            "status": "DESCRIPTIVE_ONLY" if converged else
                "UNAVAILABLE_NOT_CONVERGED",
            "converged": converged, "iterations": iterations,
            "coefficients": beta.tolist(), "support_lower": lower,
            "support_upper": upper,
            "target_probability_range": [float(np.min(target_probability)),
                                         float(np.max(target_probability))],
            "comparison_probability_range": [
                float(np.min(comparison_probability)),
                float(np.max(comparison_probability))],
        }

    for library in [f"lib{value}" for value in TARGET_LIBRARIES] + ["ALL_TARGET"]:
        choose = (lambda row: True) if library == "ALL_TARGET" else \
            (lambda row, selected=library: row["library"] == selected)
        groups = (
            ("BEFORE_MATCHING", [row for row in targets if choose(row)],
             [row for row in comparison_pool if choose(row)]),
            ("AFTER_MATCHING", [row for row in paired_targets if choose(row)],
             [row for row in paired_matches if choose(row)]),
        )
        for stage, target_group, comparison_group in groups:
            target_complete = complete(target_group)
            comparison_complete = complete(comparison_group)
            for index, field in enumerate(COVERAGE_COVARIATES):
                difference = standardized_difference(
                    [coverage_vector(row)[index] for row in target_complete],
                    [coverage_vector(row)[index] for row in comparison_complete])
                balance.append({
                    "balance_stage": stage, "library": library,
                    "covariate": field,
                    "target_sample_size": len(target_group),
                    "comparison_sample_size": len(comparison_group),
                    "target_complete_case_size": len(target_complete),
                    "comparison_complete_case_size": len(comparison_complete),
                    "missing_exclusions": len(target_group) + len(comparison_group) -
                        len(target_complete) - len(comparison_complete),
                    "absolute_standardized_difference": abs(difference)
                        if math.isfinite(difference) else math.nan,
                    "balance_threshold": 0.10,
                    "balance_status": "PASS" if math.isfinite(difference) and
                        abs(difference) <= 0.10 else "FAIL",
                    "missing_covariate_policy":
                        "COMPLETE_CASE_ALL_PREDECLARED_COVARIATES",
                })
        target_group = [row for row in targets if choose(row)]
        comparison_group = [row for row in comparison_pool if choose(row)]
        target_weights, comparison_weights, overlap_info = overlap_weights(
            target_group, comparison_group)
        target_complete = [row for row in complete(target_group)
                           if (row["library"], row["barcode"])
                           in target_weights]
        comparison_complete = [row for row in complete(comparison_group)
                               if (row["library"], row["barcode"])
                               in comparison_weights]
        target_weight_values = [target_weights.get(
            (row["library"], row["barcode"]), 0.0) for row in target_complete]
        comparison_weight_values = [comparison_weights.get(
            (row["library"], row["barcode"]), 0.0)
            for row in comparison_complete]
        target_ess = effective_sample_size(target_weight_values)
        comparison_ess = effective_sample_size(comparison_weight_values)
        for index, field in enumerate(COVERAGE_COVARIATES):
            t_mean, t_var = weighted_moments(
                [coverage_vector(row)[index] for row in target_complete],
                target_weight_values)
            c_mean, c_var = weighted_moments(
                [coverage_vector(row)[index] for row in comparison_complete],
                comparison_weight_values)
            pooled = math.sqrt((t_var + c_var) / 2) \
                if math.isfinite(t_var) and math.isfinite(c_var) and \
                t_var + c_var > 0 else math.nan
            difference = (t_mean - c_mean) / pooled \
                if math.isfinite(pooled) and pooled > 0 else math.nan
            balance.append({
                "balance_stage": "OVERLAP_WEIGHTED", "library": library,
                "covariate": field,
                "target_sample_size": len(target_group),
                "comparison_sample_size": len(comparison_group),
                "target_complete_case_size": len(target_complete),
                "comparison_complete_case_size": len(comparison_complete),
                "target_effective_sample_size": target_ess,
                "comparison_effective_sample_size": comparison_ess,
                "missing_exclusions": len(target_group) + len(comparison_group) -
                    len(target_complete) - len(comparison_complete),
                "absolute_standardized_difference": abs(difference)
                    if math.isfinite(difference) else math.nan,
                "balance_threshold": 0.10,
                "balance_status": "PASS" if math.isfinite(difference) and
                    abs(difference) <= 0.10 else "FAIL",
                "method_status": overlap_info.get("status", "UNAVAILABLE"),
                "interpretation": "DESCRIPTIVE_ONLY",
                "propensity_converged": overlap_info.get("converged", False),
                "propensity_iterations": overlap_info.get("iterations", 0),
                "propensity_coefficients": json.dumps(
                    overlap_info.get("coefficients", [])),
                "propensity_common_support_lower": overlap_info.get(
                    "support_lower", math.nan),
                "propensity_common_support_upper": overlap_info.get(
                    "support_upper", math.nan),
                "missing_covariate_policy":
                    "COMPLETE_CASE_ALL_PREDECLARED_COVARIATES",
            })
        for group_name, group_rows, group_weights in (
                ("TARGET", target_complete, target_weights),
                ("COMPARISON", comparison_complete, comparison_weights)):
            for row in group_rows:
                balance.append({
                    "balance_stage": "OVERLAP_WEIGHT",
                    "library": library,
                    "covariate": "ROW_WEIGHT",
                    "group": group_name,
                    "cell_id": target_id(row),
                    "weight": group_weights[(row["library"], row["barcode"])],
                    "method_status": overlap_info.get(
                        "status", "UNAVAILABLE"),
                    "interpretation": "DESCRIPTIVE_ONLY",
                })
    return targets, matches, manifest, balance


def frozen_target_evidence_rows(targets, changed_rows):
    changed = {(row["library"], row["barcode"]): row for row in changed_rows}
    rows = []
    for row in sorted(targets, key=lambda item: (item["library"], item["barcode"])):
        indexed = target_id(row)
        output = {
            "target_id": indexed, "library": row["library"],
            "barcode": row["barcode"],
            "locked_source": row.get("reconciled_identity_locked", ""),
            "current_proposed_contributor": row.get(
                "rna_site_genetic_second_state", ""),
            "complete_candidate_menu_signature": row.get(
                "canonical_menu_signature", ""),
            "complete_candidate_menu_size": row.get("menu_size", 0),
            "site_level_rna_atac_agreement": (
                clean(row.get("rna_site_genetic_second_state", "")) ==
                clean(row.get("atac_site_genetic_second_state", "")) and
                bool(clean(row.get("rna_site_genetic_second_state", "")))),
            "molecule_level_rna_atac_agreement": (
                clean(row.get("rna_molecule_genetic_second_state", "")) ==
                clean(row.get("atac_molecule_genetic_second_state", "")) and
                bool(clean(row.get("rna_molecule_genetic_second_state", "")))),
            "strict_site_five": indexed in STRICT_SITE_TARGETS,
            "legacy_dual_molecule_p95_nine": truthy(row.get(
                "molecule_genetic_both_assays_p95", "")),
            "legacy_exact_molecule_supported_seven":
                indexed in LEGACY_EXACT_MOLECULE_TARGETS,
            "strict_exact_overlap_three": indexed in (
                STRICT_SITE_TARGETS & LEGACY_EXACT_MOLECULE_TARGETS),
            "site_plain_language_category": row.get(
                "site_existing_evidence_category", "insufficient evidence"),
            "molecule_plain_language_category": row.get(
                "molecule_existing_evidence_category", "insufficient evidence"),
            "primary_site_coverage_calibrated_dual_p95": row.get(
                "menu_coverage_calibrated_dual_p95", ""),
            "primary_molecule_coverage_calibrated_dual_p95": row.get(
                "primary_molecule_calibrated_dual_p95", ""),
            "changed_conclusions": changed.get(
                (row["library"], row["barcode"]), {}).get(
                    "change_reasons", "NONE"),
        }
        for modality in ("rna", "atac"):
            for evidence in ("site", "molecule"):
                prefix = f"{modality}_{evidence}_genetic"
                output.update({
                    f"{prefix}_winner": row.get(f"{prefix}_second_state", ""),
                    f"{prefix}_runner_up": row.get(
                        f"{prefix}_runner_up_second_state", ""),
                    f"{prefix}_margin": row.get(f"{prefix}_margin", ""),
                    f"{prefix}_score_improvement": row.get(
                        f"{prefix}_delta", ""),
                    f"{prefix}_fitted_fraction": row.get(
                        f"{prefix}_fraction", ""),
                    f"{prefix}_interval_low": row.get(
                        f"{prefix}_profile_low", ""),
                    f"{prefix}_interval_high": row.get(
                        f"{prefix}_profile_high", ""),
                    f"{prefix}_percentile": row.get(
                        f"{prefix}_percentile", ""),
                    f"{prefix}_availability": bool(clean(row.get(
                        f"{prefix}_second_state", ""))),
                    f"{prefix}_winner_equals_site_contributor": (
                        clean(row.get(f"{prefix}_second_state", "")) ==
                        clean(row.get(f"{modality}_site_genetic_second_state", ""))
                        if evidence == "molecule" else True),
                })
            output[f"{modality}_molecule_evidence_basis"] = row.get(
                f"{modality}_molecule_evidence_basis_class", "UNAVAILABLE")
            output[f"{modality}_molecule_basis_explanation"] = (
                "RNA molecules grouped using UMI/gene linkage"
                if modality == "rna" else
                "ATAC read-name-based units used when a stronger molecule identifier is unavailable")
        rows.append(output)
    return rows


def control_parent_rows(cells, excluded_keys):
    parents = []
    for row in cells:
        if int(row["library"].removeprefix("lib")) not in ALLOWED_LIBRARIES:
            continue
        components, distinct = component_shape(
            row.get("reconciled_identity_locked", ""))
        both_site = bool(clean(row.get("rna_site_unrestricted_second_state", ""))) and \
            bool(clean(row.get("atac_site_unrestricted_second_state", "")))
        if not components:
            continue
        exclusions = []
        if not both_site:
            exclusions.append("COMPLETE_SITE_MENU_NOT_EVALUABLE")
        if clean(row.get("demux_modality_relation", "")).upper() != "AGREE":
            exclusions.append("RNA_ATAC_IDENTITY_NOT_CONCORDANT")
        droplet_state = clean(row.get("reconciled_droplet_state", "")).upper()
        if droplet_state != "SINGLE_CELL":
            exclusions.append("NOT_RECONCILED_SINGLE_CELL")
        if clean(row.get("identity_disposition", "")).upper() == "REVIEW_NEEDED":
            exclusions.append("IDENTITY_REVIEW_NEEDED")
        if (row["library"], row["barcode"]) in excluded_keys:
            exclusions.append("FROZEN_TARGET_OR_MATCHED_COMPARISON")
        legal_evidence = []
        for candidate in row.get("_candidates", []):
            if all(candidate[f"{modality}_{evidence}"]["available"]
                   for modality in ("rna", "atac")
                   for evidence in ("site", "molecule")):
                legal_evidence.append(candidate["second_state"])
        if not legal_evidence:
            exclusions.append("NO_LEGAL_CONTRIBUTOR_WITH_SITE_AND_MOLECULE_BOTH_ASSAYS")
        is_eligible = not exclusions
        parents.append({
            "library": row["library"], "barcode": row["barcode"],
            "locked_identity": row.get("reconciled_identity_locked", ""),
            "component_count": components, "distinct_donor_count": distinct,
            "simple_donor": distinct == 1,
            "parent_role": "SOURCE_OR_RECIPIENT" if is_eligible else "NOT_ELIGIBLE",
            "control_parent_eligible": is_eligible,
            "exclusion_reason": ";".join(exclusions),
            "rna_coverage_quintile": row.get(
                "rna_candidate_independent_coverage_quintile", ""),
            "atac_coverage_quintile": row.get(
                "atac_candidate_independent_coverage_quintile", ""),
            "eligibility": (
                "CONFIDENT_SIMPLE_DONOR" if distinct == 1 else
                "COMPOSITE_RECIPIENT_FOR_DOSAGE_SHIFT")
                if is_eligible else "INELIGIBLE",
            "legal_evidence_contributors": ",".join(sorted(set(legal_evidence))),
            "_source_row": row,
        })
    return parents


def select_control_pairs(parent_rows, maximum_per_class=100):
    """Deterministic role split and balanced parent reuse, without a Cartesian graph."""
    pairs, accounting = [], []
    eligible = [row for row in parent_rows if truthy(row.get("control_parent_eligible"))]
    source_usage, recipient_usage, donor_usage = Counter(), Counter(), Counter()
    for library in (f"lib{value}" for value in ALLOWED_LIBRARIES):
        rows = [row for row in eligible if row["library"] == library]
        by_donor = defaultdict(list)
        for row in rows:
            if row["simple_donor"]:
                by_donor[row["locked_identity"]].append(row)
        source_keys = set()
        for donor, members in sorted(by_donor.items()):
            members.sort(key=lambda row: (stable_seed("CONTROL_ROLE", library, donor, row["barcode"]), row["barcode"]))
            # The split is fixed before recovery is measured and across all classes.
            source_keys.update(row["barcode"] for row in members[:len(members)//2])
        sources = [row for row in rows if row["barcode"] in source_keys]
        recipients = [row for row in rows if row["barcode"] not in source_keys]
        for control_class in ("PLANTED_NEW_SOURCE", "ALREADY_PRESENT_DONOR_DOSAGE_SHIFT", "SAME_SOURCE_NULL"):
            def compatible(recipient, source):
                contributor = source["locked_identity"]
                locked = set(donor_components(recipient["locked_identity"]))
                legal = set(clean(recipient.get("legal_evidence_contributors", "")).split(","))
                valid = (contributor in legal and contributor not in locked
                         if control_class == "PLANTED_NEW_SOURCE" else
                         contributor in legal and len(locked) >= 2 and contributor in locked
                         if control_class == "ALREADY_PRESENT_DONOR_DOSAGE_SHIFT" else
                         recipient["simple_donor"] and contributor == recipient["locked_identity"])
                return valid and all(source["_source_row"].get(f"{assay}_molecule_evidence_basis_class", "UNAVAILABLE") ==
                                     recipient["_source_row"].get(f"{assay}_molecule_evidence_basis_class", "UNAVAILABLE")
                                     for assay in ("rna", "atac"))
            def source_group(source):
                return (source["locked_identity"], *(source["_source_row"].get(
                    f"{assay}_molecule_evidence_basis_class", "UNAVAILABLE") for assay in ("rna", "atac")))
            group_examples = {source_group(source): source for source in sources}
            recipient_groups = {key: [row for row in recipients if compatible(row, source)]
                                for key, source in group_examples.items()}
            selected, seen = [], set()
            for _ in range(maximum_per_class):
                best = None
                ordered_recipients = {key: sorted((row for row in members if recipient_usage[(library,row["barcode"])] < CONTROL_RECIPIENT_REUSE_CAP),
                    key=lambda row: (recipient_usage[(library,row["barcode"])],row["barcode"]))
                    for key,members in recipient_groups.items()}
                for source in sources:
                    skey = (library, source["barcode"])
                    dkey = (library, source["locked_identity"])
                    if source_usage[skey] >= CONTROL_SOURCE_CELL_REUSE_CAP or donor_usage[dkey] >= CONTROL_SOURCE_DONOR_REUSE_CAP:
                        continue
                    recipient = next((row for row in ordered_recipients[source_group(source)]
                                      if (row["barcode"],source["barcode"]) not in seen), None)
                    if recipient is None:
                        continue
                    rkey = (library, recipient["barcode"])
                    rank = (source_usage[skey], recipient_usage[rkey], donor_usage[dkey], source["barcode"],recipient["barcode"])
                    if best is None or rank < best[0]:
                        best = (rank, recipient, source)
                if best is None:
                    break
                _, recipient, source = best
                contributor = source["locked_identity"]
                key = (control_class, library, recipient["barcode"], source["barcode"], contributor)
                seen.add((recipient["barcode"], source["barcode"]))
                source_usage[(library, source["barcode"])] += 1
                recipient_usage[(library, recipient["barcode"])] += 1
                donor_usage[(library, contributor)] += 1
                row = {
                    "schema_version": "joint_doublet_control_manifest_v4",
                    "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
                    "calibration_library": CALIBRATION_LIBRARY,
                    "control_id": f"{library}:{control_class}:{stable_seed(*key):016x}",
                    "control_class": control_class, "library": library,
                    "recipient_barcode": recipient["barcode"], "source_barcode": source["barcode"],
                    "recipient_identity": recipient["locked_identity"], "source_identity": contributor,
                    "expected_contributor": contributor,
                    "legal_evidence_contributors": recipient["legal_evidence_contributors"],
                    "requested_fractions": ",".join(f"{value:.2f}" for value in CONTROL_FRACTIONS),
                    "genotype_distance_matched_decoy": "PENDING_CACHE_FINALIZER" if control_class == "PLANTED_NEW_SOURCE" else "",
                    "decoy_status": "PROVISIONAL_PENDING_NORMALIZED_CACHE" if control_class == "PLANTED_NEW_SOURCE" else "NOT_APPLICABLE",
                    "parent_role_separation": "GLOBAL_DISJOINT_ROLE_SPLIT_BEFORE_OUTCOMES",
                    "seed": stable_seed("control", *key), "source_cluster": contributor,
                    "source_cell_cluster": f"{library}:{source['barcode']}", "source_donor_cluster": contributor,
                }
                for assay in ("rna", "atac"):
                    row[f"{assay}_molecule_evidence_basis"] = recipient["_source_row"].get(f"{assay}_molecule_evidence_basis_class", "UNAVAILABLE")
                    row[f"source_{assay}_molecule_evidence_basis"] = source["_source_row"].get(f"{assay}_molecule_evidence_basis_class", "UNAVAILABLE")
                selected.append(row)
            pairs.extend(selected)
            for contributor in ["ALL"] + sorted({row["expected_contributor"] for row in selected}):
                members = [row for row in selected if contributor == "ALL" or row["expected_contributor"] == contributor]
                usage = Counter(row["source_barcode"] for row in members)
                rusage = Counter(row["recipient_barcode"] for row in members)
                accounting.append({
                    "library": library, "control_class": control_class, "contributor": contributor,
                    "requested_pairs": maximum_per_class if contributor == "ALL" else "CLASS_SHARED_QUOTA",
                    "selected_pairs": len(members), "unique_sources": len(usage), "unique_recipients": len(rusage),
                    "capacity_limited": len(selected) < maximum_per_class,
                    "capacity_reason": "ELIGIBILITY_DISJOINT_ROLES_AND_GLOBAL_REUSE_CAPS" if len(selected) < maximum_per_class else "NONE",
                    "maximum_source_reuse": max(usage.values(), default=0),
                    "source_reuse_distribution": json.dumps(dict(sorted(usage.items()))),
                    "recipient_reuse_distribution": json.dumps(dict(sorted(rusage.items()))),
                    "effective_source_cluster_count": effective_sample_size(list(usage.values())),
                    "source_cell_reuse_cap": CONTROL_SOURCE_CELL_REUSE_CAP,
                    "source_donor_reuse_cap": CONTROL_SOURCE_DONOR_REUSE_CAP,
                    "recipient_reuse_cap": CONTROL_RECIPIENT_REUSE_CAP,
                    "independent_parent_information": "REPORT_UNIQUE_PARENTS_NOT_FRACTIONS_AS_REPLICATES",
                    "selection_rule": "FIXED_ROLE_SPLIT_THEN_MINIMUM_GLOBAL_SOURCE_RECIPIENT_REUSE",
                })
    return pairs, accounting


def assign_genotype_distance_decoys(by_parent, control_pairs):
    """Choose an actual legal-menu decoy matched on panel discrimination."""
    planted = {
        (row["library"], row["recipient_barcode"], row["expected_contributor"])
        for row in control_pairs
        if row["control_class"] == "PLANTED_NEW_SOURCE"
    }
    assignments = {}
    for library, barcode, expected in planted:
        choices = by_parent.get((library, barcode), [])
        expected_rows = [row for row in choices
                         if row["second_state"] == expected]
        if not expected_rows:
            continue
        target_distance = max(row["distance"] for row in expected_rows)
        recipient_identity = next((
            row["recipient_identity"] for row in control_pairs
            if row["library"] == library and
            row["recipient_barcode"] == barcode and
            row["expected_contributor"] == expected), "")
        recipient_donors = set(donor_components(recipient_identity))
        alternatives = [
            row for row in choices
            if row["second_state"] != expected and row["physical"] and
            not row["component_only"] and
            set(donor_components(row["second_state"])) - recipient_donors
        ]
        if not alternatives:
            continue
        decoy = min(alternatives, key=lambda row: (
            abs(row["distance"] - target_distance), row["second_state"]))
        assignments[(library, barcode, expected)] = (
            decoy["second_state"], target_distance, decoy["distance"],
            abs(decoy["distance"] - target_distance))
    for row in control_pairs:
        if row["control_class"] != "PLANTED_NEW_SOURCE":
            row["genotype_distance_metric"] = "NOT_APPLICABLE"
            row["genotype_distance_matched_decoy"] = ""
            row["expected_contributor_genotype_distance"] = math.nan
            row["decoy_genotype_distance"] = math.nan
            row["decoy_distance_absolute_difference"] = math.nan
            continue
        assignment = assignments.get((
            row["library"], row["recipient_barcode"],
            row["expected_contributor"]))
        row["genotype_distance_metric"] = \
            "MEAN_RNA_ATAC_DISCRIMINATING_SITES"
        if assignment is None:
            row["genotype_distance_matched_decoy"] = ""
            row["expected_contributor_genotype_distance"] = math.nan
            row["decoy_genotype_distance"] = math.nan
            row["decoy_distance_absolute_difference"] = math.nan
        else:
            (row["genotype_distance_matched_decoy"],
             row["expected_contributor_genotype_distance"],
             row["decoy_genotype_distance"],
             row["decoy_distance_absolute_difference"]) = assignment


def reconstruct_target_manifests(candidates_path, selected_keys, targeted_root,
                                 control_pairs):
    handles = {}
    writers = {}
    counts = Counter()
    fields = [
        "schema_version", "scientific_method_version", "calibration_library",
        "library", "barcode", "candidate_id",
        "locked_state", "locked_copy_vector", "second_state",
        "second_copy_vector", "candidate_origin", "exhaustive_fallback",
        "nomination_modalities", "rho", "ambient_copy_vector", "ambient_status",
        "candidate_policy", "physical_pool_state", "component_only_state",
        "structural_added_state_relationship",
    ]
    try:
        for number in ALLOWED_LIBRARIES:
            library = f"lib{number}"
            for modality in ("rna", "atac"):
                path = targeted_root / "manifests" / \
                    f"{library}.{modality}_targeted_manifest.tsv.gz"
                path.parent.mkdir(parents=True, exist_ok=True)
                temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
                handle = gzip.open(temporary, "wt", encoding="utf-8", newline="")
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                        lineterminator="\n")
                writer.writeheader()
                handles[(library, modality)] = (handle, temporary, path)
                writers[(library, modality)] = writer
        for row in read_tsv(candidates_path):
            library = clean(row.get("library", ""))
            barcode = clean(row.get("barcode", ""))
            if (library, barcode) not in selected_keys:
                continue
            for modality in ("rna", "atac"):
                output = {
                    "schema_version": "joint_doublet_candidate_manifest_v4",
                    "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
                    "calibration_library": CALIBRATION_LIBRARY,
                    "library": library, "barcode": barcode,
                    "candidate_id": row.get("candidate_id", ""),
                    "locked_state": row.get("locked_state", ""),
                    "locked_copy_vector": row.get("locked_copy_vector", ""),
                    "second_state": row.get("second_state", ""),
                    "second_copy_vector": row.get("second_copy_vector", ""),
                    "candidate_origin": row.get("candidate_origin", ""),
                    "exhaustive_fallback": row.get("exhaustive_fallback", ""),
                    "nomination_modalities": row.get("nomination_modalities", ""),
                    "rho": row.get(f"{modality}_rho_requested", "0"),
                    "ambient_copy_vector": row.get(
                        f"{modality}_ambient_copy_vector", ""),
                    "ambient_status": row.get(f"{modality}_ambient_status", ""),
                    "candidate_policy": row.get("candidate_policy", ""),
                    "physical_pool_state": row.get("physical_pool_state", ""),
                    "component_only_state": row.get("component_only_state", ""),
                    "structural_added_state_relationship": row.get(
                        "structural_added_state_relationship", ""),
                }
                writers[(library, modality)].writerow(output)
                counts[(library, modality, "candidate_rows")] += 1
                counts[(library, modality, "site_cache_rows")] += int(max(
                    finite(row.get(
                        f"{modality}_n_common_nuclear_sites", ""), 0.0), 0.0))
                counts[(library, modality, "linked_units")] += int(max(
                    finite(row.get(
                        f"{modality}_n_independent_linked_units", ""), 0.0), 0.0))
                counts[(library, modality, "molecule_cache_rows")] += int(max(
                    finite(row.get(
                        f"{modality}_total_snps_in_linked_units", ""), 0.0), 0.0))
    finally:
        for handle, temporary, path in handles.values():
            handle.close()
            os.replace(temporary, path)
    return counts


def _task_input(task_row, modality, name):
    return clean(task_row.get(f"{modality}_pileup_{name}", "") or
                 task_row.get(f"{modality}_{name}", ""))




def _small_file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_fnv1a64(path):
    """Native bounded checksum; never loop over gigabytes one Python byte at a time."""
    scorer = Path(__file__).resolve().parent / "tetra_score_calls"
    result = subprocess.run([str(scorer), "--joint-doublet-cache-digest", str(path)],
                            capture_output=True, text=True, check=False)
    value = result.stdout.strip()
    if result.returncode or not re.fullmatch(r"fnv1a64:[0-9a-f]{16}", value):
        raise RuntimeError("native cache integrity check failed: " + result.stderr.strip())
    return value


def _manifest_content_digest(path):
    """FNV-1a over decompressed manifest bytes, matching the C++ reader."""
    value = 1469598103934665603
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            for byte in block:
                value ^= byte
                value = (value * 1099511628211) & ((1 << 64) - 1)
    return f"fnv1a64:{value:016x}"


def _candidate_record_estimate(row, modality):
    offset = 0 if modality == "rna" else 2
    counts = [candidate.get("_cache_record_counts", (0, 0, 0, 0))
              for candidate in row.get("_candidates", [])]
    # Candidate evidence sets need not be nested.  Summing candidate counts is
    # a conservative upper bound on the complete-menu union; using the largest
    # single candidate would be an unsafe lower bound for memory planning.
    # Linked units and linked SNP observations are different quantities. Only
    # the latter project fixed-width molecule records. Evidence absent from all
    # candidate supports is measured by extraction and is not estimated here.
    return sum(value[offset] for value in counts), sum(value[offset+1] for value in counts)


def _validated_targeted_resources(policy):
    resolved = {
        "extraction_wall_time": policy.get("extraction_wall_time", "12:00:00"),
        "extraction_memory": policy.get("extraction_memory", "64G"),
        # Normalized extraction is a one-pass streaming stage; its scientific
        # work is serial.  Request one CPU until a measured parallel worker
        # implementation satisfies the E06 utilization contract.
        "extraction_cpus": 1,
        "analysis_wall_time": policy.get("analysis_wall_time", "08:00:00"),
        "analysis_memory": policy.get("analysis_memory", "48G"),
        "analysis_cpus": int(policy.get("analysis_cpus", 8)),
        "finalizer_wall_time": policy.get("finalizer_wall_time", "01:00:00"),
        "finalizer_memory": policy.get("finalizer_memory", "16G"),
        "finalizer_cpus": 1,
        "gather_wall_time": policy.get("gather_wall_time", "04:00:00"),
        "gather_memory": policy.get("gather_memory", "64G"),
        "gather_cpus": 1,
        "partition": clean(policy.get("partition", "compute")),
    }
    time_pattern = re.compile(r"^(?:\d+-)?\d{2}:\d{2}:\d{2}$")
    memory_pattern = re.compile(r"^[1-9]\d*[KMGTP]$")
    for key in ("extraction_wall_time", "analysis_wall_time",
                "finalizer_wall_time", "gather_wall_time"):
        if not time_pattern.fullmatch(str(resolved[key])):
            raise RuntimeError(f"invalid targeted Slurm time: {key}={resolved[key]}")
    for key in ("extraction_memory", "analysis_memory", "finalizer_memory",
                "gather_memory"):
        value = str(resolved[key]).upper()
        if not memory_pattern.fullmatch(value):
            raise RuntimeError(f"invalid targeted Slurm memory: {key}={value}")
        amount = int(value[:-1])
        unit = value[-1]
        powers = {"K": 10, "M": 20, "G": 30, "T": 40, "P": 50}
        raw_bytes = amount * (1 << powers[unit])
        mib = (raw_bytes + (1 << 20) - 1) // (1 << 20)
        scheduler_bytes = mib * (1 << 20)
        reserve = max(16 * (1 << 20), (scheduler_bytes + 9) // 10)
        if reserve >= scheduler_bytes:
            raise RuntimeError(f"targeted Slurm memory is too small: {key}")
        resolved[key] = f"{mib}M"
        prefix = key.removesuffix("_memory")
        resolved[f"{prefix}_scheduler_memory_bytes"] = scheduler_bytes
        resolved[f"{prefix}_launcher_runtime_reserve_bytes"] = reserve
        resolved[f"{prefix}_worker_memory_budget_bytes"] = \
            scheduler_bytes - reserve
    for key in ("extraction_cpus", "analysis_cpus", "finalizer_cpus",
                "gather_cpus"):
        if resolved[key] < 1 or resolved[key] > 256:
            raise RuntimeError(f"invalid targeted Slurm CPU count: {key}")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", resolved["partition"]):
        raise RuntimeError("invalid targeted Slurm partition")
    return resolved


def _write_targeted_scripts(root, tool_bin_root, generation_id,
                            extraction_tasks, analysis_tasks, benchmark=False,
                            production_threads=8, resource_policy=None):
    root = Path(root)
    control_root = root
    resource_policy = _validated_targeted_resources(resource_policy or {})
    extraction_time = resource_policy["extraction_wall_time"]
    extraction_memory = resource_policy["extraction_memory"]
    extraction_cpus = resource_policy["extraction_cpus"]
    analysis_time = resource_policy["analysis_wall_time"]
    analysis_memory = resource_policy["analysis_memory"]
    analysis_cpus = resource_policy["analysis_cpus"]
    finalizer_time = resource_policy["finalizer_wall_time"]
    finalizer_memory = resource_policy["finalizer_memory"]
    finalizer_cpus = resource_policy["finalizer_cpus"]
    gather_time = resource_policy["gather_wall_time"]
    gather_memory = resource_policy["gather_memory"]
    gather_cpus = resource_policy["gather_cpus"]
    partition = resource_policy["partition"]
    scorer = Path(tool_bin_root) / "tetra_score_calls"
    helper = Path(tool_bin_root) / "joint_doublet_no_rescore.py"
    extraction_manifest = root / "manifests" / "extraction_tasks.tsv"
    analysis_manifest = root / "manifests" / "analysis_tasks.tsv"
    extraction_loader = f'''TASK_MANIFEST={shlex.quote(str(extraction_manifest))}
TASK_LINE="$(awk -v task="$SLURM_ARRAY_TASK_ID" 'NR == task + 2 {{print; exit}}' "$TASK_MANIFEST" | tr '\\t' '\\034')"
[[ -n "$TASK_LINE" ]] || {{ echo "missing extraction task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }}
IFS=$'\\034' read -r task_index scientific_method_version calibration_library workload_generation_id cache_generation_id library modality candidate_manifest manifest_digest samples pileup_sites pileup_observations pileup_molecules cache_prefix marker input_status input_issues selected_cells projected_observation_records projected_molecule_records error_ref error_alt min_evidence max_second_fraction projected_peak_bytes bounded_memory_limit_bytes scheduler_memory_bytes launcher_runtime_reserve_bytes worker_memory_budget_bytes cache_policy promoted_cache_source_generation promoted_cache_proof_sha256 <<< "$TASK_LINE"
'''
    extract_script = f'''#!/bin/bash
#SBATCH --job-name=jd_norm_extract
#SBATCH --output={root}/logs/extract_%A_%a.out
#SBATCH --error={root}/logs/extract_%A_%a.err
#SBATCH --time={extraction_time}
#SBATCH --cpus-per-task={extraction_cpus}
#SBATCH --mem={extraction_memory}
#SBATCH --partition={partition}
#SBATCH --nodes=1

set -uo pipefail
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
module purge
module load miniforge/3 genomics-base/latest htslib/1.20 || exit $?
{extraction_loader}
python3 -B {shlex.quote(str(helper))} validate-task \
    --task-kind extraction --task-manifest "$TASK_MANIFEST" \
    --task-index "$SLURM_ARRAY_TASK_ID" --input-only --quiet || exit $?
if python3 -B {shlex.quote(str(helper))} validate-task \
    --task-kind extraction --task-manifest "$TASK_MANIFEST" \
    --task-index "$SLURM_ARRAY_TASK_ID" --quiet; then
  exit 0
fi
[[ "$input_status" == "READY" ]] || {{ echo "$input_issues" >&2; exit 2; }}
task_temp={root}/task_scratch/extract_${{SLURM_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}
mkdir -p "$task_temp" || exit $?
time_file={root}/analysis/extract_${{task_index}}.time.txt
interval_file={root}/analysis/extract_${{task_index}}.interval.tsv
start_ns="$(date +%s%N)" || exit $?
/usr/bin/time -v -o "$time_file" {shlex.quote(str(scorer))} \
  --joint-doublet-normalized-cache-prefix "$cache_prefix" \
  --joint-doublet-manifest "$candidate_manifest" \
  --joint-doublet-manifest-digest "$manifest_digest" \
  --joint-doublet-workload-generation "$workload_generation_id" \
  --joint-doublet-temp-dir "$task_temp" --samples "$samples" \
  --pileup-sites "$pileup_sites" --pileup-observations "$pileup_observations" \
  --pileup-molecules "$pileup_molecules" --libname "$library" \
  --modality "$modality" --error_ref "$error_ref" --error_alt "$error_alt" \
  --min_evidence "$min_evidence" --max-second-fraction "$max_second_fraction" \
  --joint-doublet-scheduler-memory-bytes "$scheduler_memory_bytes" \
  --joint-doublet-launcher-runtime-reserve-bytes "$launcher_runtime_reserve_bytes" \
  --joint-doublet-worker-memory-budget-bytes "$worker_memory_budget_bytes" \
  --threads {extraction_cpus}
status=$?
end_ns="$(date +%s%N)" || exit $?
printf 'start_ns\tend_ns\tstatus\n%s\t%s\t%s\n' \
  "$start_ns" "$end_ns" "$status" > "$interval_file.tmp.$$" && \
  mv -f "$interval_file.tmp.$$" "$interval_file" || exit $?
if [[ $status -eq 0 ]]; then
  python3 -B {shlex.quote(str(helper))} validate-task \
    --task-kind extraction --task-manifest "$TASK_MANIFEST" \
    --task-index "$SLURM_ARRAY_TASK_ID" --write-marker
  status=$?
fi
rm -rf -- "$task_temp"
exit $status
'''
    analysis_loader = f'''if [[ -n "${{JOINT_DOUBLET_TASK_INDEX_MAP:-}}" ]]; then
  actual_index="$(awk -v idx="$SLURM_ARRAY_TASK_ID" 'NR == idx + 1 {{print; exit}}' "$JOINT_DOUBLET_TASK_INDEX_MAP")"
  [[ "$actual_index" =~ ^[0-9]+$ ]] || {{ echo "invalid stable task index mapping" >&2; exit 1; }}
  export SLURM_ARRAY_TASK_ID="$actual_index"
fi
TASK_MANIFEST={shlex.quote(str(analysis_manifest))}
TASK_LINE="$(awk -v task="$SLURM_ARRAY_TASK_ID" 'NR == task + 2 {{print; exit}}' "$TASK_MANIFEST" | tr '\\t' '\\034')"
[[ -n "$TASK_LINE" ]] || {{ echo "missing analysis task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }}
IFS=$'\\034' read -r task_index scientific_method_version calibration_library workload_generation_id cache_generation_id action library modality role barcode control_ids cache_prefix candidate_manifest control_manifest analysis_output marker predicted_evidence_records predicted_likelihood_work error_ref error_alt min_evidence max_second_fraction scheduler_memory_bytes launcher_runtime_reserve_bytes worker_memory_budget_bytes <<< "$TASK_LINE"
'''
    if benchmark:
        analysis_body = f'''#!/bin/bash
#SBATCH --job-name=jd_benchmark
#SBATCH --output={root}/logs/analysis_%A_%a.out
#SBATCH --error={root}/logs/analysis_%A_%a.err
#SBATCH --time={analysis_time}
#SBATCH --cpus-per-task={analysis_cpus}
#SBATCH --mem={analysis_memory}
#SBATCH --partition={partition}
#SBATCH --nodes=1

set -uo pipefail
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
module purge
module load miniforge/3 genomics-base/latest htslib/1.20 || exit $?
{analysis_loader}
python3 -B {shlex.quote(str(helper))} validate-task \
  --task-kind analysis --task-manifest "$TASK_MANIFEST" \
  --task-index "$SLURM_ARRAY_TASK_ID" --input-only --quiet || exit $?
extraction_index="$(awk -F '\t' -v lib="$library" -v mod="$modality" \
  'NR > 1 && $6 == lib && $7 == mod {{print $1; exit}}' \
  {shlex.quote(str(root / 'manifests' / 'extraction_tasks.tsv'))})"
[[ -n "$extraction_index" ]] || {{ echo "missing matching extraction task" >&2; exit 1; }}
python3 -B {shlex.quote(str(helper))} validate-task \
  --task-kind extraction \
  --task-manifest {shlex.quote(str(root / 'manifests' / 'extraction_tasks.tsv'))} \
  --task-index "$extraction_index" --quiet || exit $?
single_output="${{analysis_output%.tsv.gz}}.threads1.tsv.gz"
production_output="${{analysis_output%.tsv.gz}}.threads{production_threads}.tsv.gz"
common=(--joint-doublet-targeted-analysis-manifest "$TASK_MANIFEST" \
  --joint-doublet-targeted-analysis-index "$SLURM_ARRAY_TASK_ID" \
  --joint-doublet-normalized-cache-prefix "$cache_prefix" \
  --joint-doublet-manifest "$candidate_manifest" \
  --joint-doublet-control-manifest "$control_manifest" \
  --joint-doublet-seed {SEED} \
  --joint-doublet-downsample-replicates {DOWNSAMPLE_REPLICATES} \
  --joint-doublet-null-replicates {CELL_NULL_REPLICATES} \
  --error_ref "$error_ref" --error_alt "$error_alt" \
  --min_evidence "$min_evidence" --max-second-fraction "$max_second_fraction" \
  --joint-doublet-reference-mode)
single_interval="${{single_output%.tsv.gz}}.interval.tsv"
start_ns="$(date +%s%N)" || exit $?
/usr/bin/time -v -o "${{single_output%.tsv.gz}}.time.txt" \
  {shlex.quote(str(scorer))} "${{common[@]}}" \
  --joint-doublet-targeted-analysis-output "$single_output" --threads 1
status=$?
end_ns="$(date +%s%N)" || exit $?
printf 'start_ns\tend_ns\tstatus\n%s\t%s\t%s\n' \
  "$start_ns" "$end_ns" "$status" > "$single_interval.tmp.$$" && \
  mv -f "$single_interval.tmp.$$" "$single_interval" || exit $?
if [[ $status -eq 0 ]]; then
  production_interval="${{production_output%.tsv.gz}}.interval.tsv"
  start_ns="$(date +%s%N)" || exit $?
  /usr/bin/time -v -o "${{production_output%.tsv.gz}}.time.txt" \
    {shlex.quote(str(scorer))} "${{common[@]}}" \
    --joint-doublet-targeted-analysis-output "$production_output" \
    --threads {production_threads}
  status=$?
  end_ns="$(date +%s%N)" || exit $?
  printf 'start_ns\tend_ns\tstatus\n%s\t%s\t%s\n' \
    "$start_ns" "$end_ns" "$status" > "$production_interval.tmp.$$" && \
    mv -f "$production_interval.tmp.$$" "$production_interval" || exit $?
fi
if [[ $status -eq 0 ]]; then
  python3 -B {shlex.quote(str(helper))} compare-equivalence \
    --left "$single_output" --right "$production_output" \
    --output "$analysis_output" --marker "$marker" \
    --task-manifest "$TASK_MANIFEST" \
    --legacy-reference {shlex.quote(str(root / 'manifests' / 'benchmark_completed_score_reference.tsv.gz'))} \
    --generation "$workload_generation_id" --task-index "$task_index"
  status=$?
fi
exit $status
'''
    else:
        analysis_body = f'''#!/bin/bash
#SBATCH --job-name=jd_target_analysis
#SBATCH --output={root}/logs/analysis_%A_%a.out
#SBATCH --error={root}/logs/analysis_%A_%a.err
#SBATCH --time={analysis_time}
#SBATCH --cpus-per-task={analysis_cpus}
#SBATCH --mem={analysis_memory}
#SBATCH --partition={partition}
#SBATCH --nodes=1

set -uo pipefail
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
module purge
module load miniforge/3 genomics-base/latest htslib/1.20 || exit $?
{analysis_loader}
python3 -B {shlex.quote(str(helper))} validate-task \
  --task-kind analysis --task-manifest "$TASK_MANIFEST" \
  --task-index "$SLURM_ARRAY_TASK_ID" --input-only --quiet || exit $?
extraction_index="$(awk -F '\t' -v lib="$library" -v mod="$modality" \
  'NR > 1 && $6 == lib && $7 == mod {{print $1; exit}}' \
  {shlex.quote(str(root / 'manifests' / 'extraction_tasks.tsv'))})"
[[ -n "$extraction_index" ]] || {{ echo "missing matching extraction task" >&2; exit 1; }}
python3 -B {shlex.quote(str(helper))} validate-task \
  --task-kind extraction \
  --task-manifest {shlex.quote(str(root / 'manifests' / 'extraction_tasks.tsv'))} \
  --task-index "$extraction_index" --quiet || exit $?
if python3 -B {shlex.quote(str(helper))} validate-task \
    --task-kind analysis --task-manifest "$TASK_MANIFEST" \
    --task-index "$SLURM_ARRAY_TASK_ID" --quiet; then
  exit 0
fi
interval_file="${{analysis_output%.tsv.gz}}.interval.tsv"
start_ns="$(date +%s%N)" || exit $?
/usr/bin/time -v -o "${{analysis_output%.tsv.gz}}.time.txt" \
  {shlex.quote(str(scorer))} \
  --joint-doublet-targeted-analysis-output "$analysis_output" \
  --joint-doublet-targeted-analysis-manifest "$TASK_MANIFEST" \
  --joint-doublet-targeted-analysis-index "$SLURM_ARRAY_TASK_ID" \
  --joint-doublet-normalized-cache-prefix "$cache_prefix" \
  --joint-doublet-manifest "$candidate_manifest" \
  --joint-doublet-control-manifest "$control_manifest" \
  --joint-doublet-seed {SEED} \
  --joint-doublet-downsample-replicates {DOWNSAMPLE_REPLICATES} \
  --joint-doublet-null-replicates {CELL_NULL_REPLICATES} \
  --error_ref "$error_ref" --error_alt "$error_alt" \
  --min_evidence "$min_evidence" --max-second-fraction "$max_second_fraction" \
  --threads {production_threads}
status=$?
end_ns="$(date +%s%N)" || exit $?
printf 'start_ns\tend_ns\tstatus\n%s\t%s\t%s\n' \
  "$start_ns" "$end_ns" "$status" > "$interval_file.tmp.$$" && \
  mv -f "$interval_file.tmp.$$" "$interval_file" || exit $?
if [[ $status -eq 0 ]]; then
  python3 -B {shlex.quote(str(helper))} validate-task \
    --task-kind analysis --task-manifest "$TASK_MANIFEST" \
    --task-index "$SLURM_ARRAY_TASK_ID" --write-marker
  status=$?
fi
exit $status
'''
    finalizer = f'''#!/bin/bash
#SBATCH --job-name=jd_control_finalize
#SBATCH --output={root}/logs/control_finalize_%j.out
#SBATCH --error={root}/logs/control_finalize_%j.err
#SBATCH --time={finalizer_time}
#SBATCH --cpus-per-task={finalizer_cpus}
#SBATCH --mem={finalizer_memory}
#SBATCH --partition={partition}
#SBATCH --nodes=1

set -uo pipefail
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
module purge
module load miniforge/3 || exit $?
python3 -B {shlex.quote(str(helper))} finalize-controls \
  --targeted-root {shlex.quote(str(control_root))} --generation {shlex.quote(generation_id)}
'''
    gather = f'''#!/bin/bash
#SBATCH --job-name=jd_target_gather
#SBATCH --output={root}/logs/gather_%j.out
#SBATCH --error={root}/logs/gather_%j.err
#SBATCH --time={gather_time}
#SBATCH --cpus-per-task={gather_cpus}
#SBATCH --mem={gather_memory}
#SBATCH --partition={partition}
#SBATCH --nodes=1

set -uo pipefail
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
module purge
module load miniforge/3 || exit $?
python3 -B {shlex.quote(str(helper))} targeted-gather \
  --targeted-root {shlex.quote(str(root))} --generation {shlex.quote(generation_id)}
'''
    scripts = {
        "extract": root / "slurm_scripts" / "01_normalized_extract.sbatch",
        "finalize_controls": root / "slurm_scripts" / "02_finalize_controls.sbatch",
        "analysis": root / "slurm_scripts" / "03_batched_analysis.sbatch",
        "gather": root / "slurm_scripts" / "04_gather.sbatch",
    }
    for name, body in (("extract", extract_script),
                       ("finalize_controls", finalizer),
                       ("analysis", analysis_body), ("gather", gather)):
        atomic_text(scripts[name], body)
        os.chmod(scripts[name], 0o755)
    return {key: str(value) for key, value in scripts.items()}


def render_targeted_workload(targeted_root, candidates_path, cells, task_rows,
                             tool_bin_root, frozen_keys, frozen_manifest,
                             resource_policy=None, model_parameters=None,
                             calibration_audit=None,
                             calibration_roster=None, frozen_comparisons=None, control_design=None):
    targeted_root = Path(targeted_root)
    resolved_resources = _validated_targeted_resources(resource_policy or {})
    model_parameters = dict(model_parameters or {})
    required_parameters = {
        "RNA": ("rna_error_ref", "rna_error_alt"),
        "ATAC": ("atac_error_ref", "atac_error_alt"),
    }
    min_evidence = int(model_parameters.get("min_evidence", -1))
    max_second_fraction = finite(
        model_parameters.get("max_second_fraction", ""))
    if min_evidence < 0 or not 0 < max_second_fraction <= 1:
        raise RuntimeError("targeted model parameter contract is incomplete")
    for modality, fields in required_parameters.items():
        ref = finite(model_parameters.get(fields[0], ""))
        alt = finite(model_parameters.get(fields[1], ""))
        if not (0 <= ref <= 1 and 0 <= alt <= 1 and ref + alt < 1):
            raise RuntimeError(
                f"invalid targeted {modality} error parameter contract")
    for name in ("manifests", "slurm_scripts", "logs", "task_scratch",
                 "scores", "cache", "analysis", "markers", "generations"):
        (targeted_root / name).mkdir(parents=True, exist_ok=True)
    if frozen_comparisons is None:
        raise RuntimeError("authoritative frozen comparison mappings are required")
    targets, matches, target_manifest, match_balance = reuse_frozen_comparisons(
        cells, frozen_comparisons)
    excluded_parent_keys = frozen_keys | {
        (row["library"], row["barcode"]) for row in matches}
    if control_design is None:
        parents = control_parent_rows(cells, excluded_parent_keys)
        control_pairs, control_accounting = select_control_pairs(parents)
    else:
        parents, control_pairs, control_accounting = control_design
    selected_keys = {(row["library"], row["barcode"])
                     for row in targets + matches}
    selected_keys.update({
        (pair["library"], barcode) for pair in control_pairs
        for barcode in (pair["recipient_barcode"], pair["source_barcode"])
    })
    unavailable_cells = [{"library": row["library"], "barcode": row["barcode"],
                          "status": "UNAVAILABLE_COMPLETE_RETAINED_MENU",
                          "reason": "Original cell retained; no replacement selected"}
                         for row in targets + matches if not row.get("_candidates")]
    unavailable_keys = {(row["library"], row["barcode"]) for row in unavailable_cells}
    selected_keys -= unavailable_keys
    write_tsv(targeted_root / "unavailable_cells.tsv", unavailable_cells,
              ("library", "barcode", "status", "reason"))

    manifest_counts = reconstruct_target_manifests(
        candidates_path, selected_keys, targeted_root, control_pairs)
    write_tsv(targeted_root / "frozen_targets_20260920.tsv",
              load_frozen_targets(frozen_manifest)[1])
    write_tsv(targeted_root / "frozen_target_comparisons.tsv", frozen_comparisons)
    # The targeted final decision must use the frozen primary reference CDF,
    # not a full-data score percentile copied from a selected target.  Keep the
    # exact audit rows and reference-cell maxima with the rendered generation.
    write_tsv(targeted_root / "primary_calibration_audit.tsv.gz",
              calibration_audit or [])
    write_tsv(targeted_root / "primary_calibration_reference_roster.tsv.gz",
              calibration_roster or [])
    write_tsv(targeted_root / "target_cells_and_matched_comparisons.tsv",
              target_manifest)
    write_tsv(targeted_root / "target_comparison_balance.tsv", match_balance)
    write_tsv(targeted_root / "control_parent_eligibility.tsv.gz",
              (public_row(row) for row in parents))
    write_tsv(targeted_root / "control_capacity_and_reuse.tsv",
              control_accounting)
    provisional_controls = targeted_root / "manifests" / \
        "control_pairs_provisional.tsv"
    final_controls = targeted_root / "manifests" / "control_pairs_final.tsv"
    write_tsv(provisional_controls, control_pairs)

    generation_payload = {
        "schema": WORKLOAD_SCHEMA,
        "implementation_version": IMPLEMENTATION_VERSION,
        "frozen_comparisons": frozen_comparisons,
        "unavailable_cells": unavailable_cells,
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "frozen_targets": sorted(f"{library}:{barcode}"
                                 for library, barcode in frozen_keys),
        "selected_cells": sorted(f"{library}:{barcode}"
                                 for library, barcode in selected_keys),
        "control_ids": sorted(row["control_id"] for row in control_pairs),
        "fractions": CONTROL_FRACTIONS,
        "downsample": [DOWNSAMPLE_FRACTIONS, DOWNSAMPLE_REPLICATES],
        "null_replicates": CELL_NULL_REPLICATES,
        "model_parameters": model_parameters,
        "resource_policy": resolved_resources,
    }
    generation_id = "prebenchmark_" + hashlib.sha256(json.dumps(
        generation_payload, sort_keys=True).encode()).hexdigest()[:16]
    cell_by_key = {(row["library"], row["barcode"]): row for row in cells}
    extraction_tasks = []
    for number in sorted({int(key[0][3:]) for key in selected_keys}):
        library = f"lib{number}"
        task_row = task_rows.get(library, {})
        for modality_lower in ("rna", "atac"):
            modality = modality_lower.upper()
            manifest = targeted_root / "manifests" / \
                f"{library}.{modality_lower}_targeted_manifest.tsv.gz"
            inputs = {
                "samples": _task_input(task_row, modality_lower, "samples"),
                "pileup_sites": _task_input(task_row, modality_lower, "sites"),
                "pileup_observations": _task_input(
                    task_row, modality_lower, "observations"),
                "pileup_molecules": _task_input(
                    task_row, modality_lower, "molecules"),
            }
            issues = [
                f"{key}={value or 'UNRESOLVED'}" for key, value in inputs.items()
                if not value or not Path(value).is_absolute() or
                _path_mentions_protected(value)
            ]
            cell_keys = [key for key in selected_keys if key[0] == library]
            estimates = [_candidate_record_estimate(
                cell_by_key[key], modality_lower) for key in cell_keys]
            ref_field, alt_field = required_parameters[modality]
            error_ref = finite(model_parameters[ref_field])
            error_alt = finite(model_parameters[alt_field])
            projected_peak_bytes = 256 * 1024 * 1024
            bounded_limit = resolved_resources["extraction_worker_memory_budget_bytes"]
            prefix = targeted_root / "cache" / \
                f"{library}.{modality_lower}.{generation_id}"
            extraction_tasks.append({
                "task_index": len(extraction_tasks),
                "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
                "calibration_library": CALIBRATION_LIBRARY,
                "workload_generation_id": generation_id,
                "cache_generation_id": generation_id,
                "library": library, "modality": modality,
                "candidate_manifest": str(manifest),
                "manifest_digest": _manifest_content_digest(manifest),
                **inputs, "cache_prefix": str(prefix),
                "marker": str(targeted_root / "markers" /
                              f"extract_{library}_{modality_lower}.complete.json"),
                "input_status": "READY" if not issues else "INCOMPLETE",
                "input_issues": ";".join(issues) or "NONE",
                "selected_cells": len(cell_keys),
                "projected_observation_records": sum(value[0] for value in estimates),
                "projected_molecule_records": sum(value[1] for value in estimates),
                "error_ref": error_ref, "error_alt": error_alt,
                "min_evidence": min_evidence,
                "max_second_fraction": max_second_fraction,
                "projected_peak_bytes": projected_peak_bytes,
                "bounded_memory_limit_bytes": bounded_limit,
                "scheduler_memory_bytes": resolved_resources[
                    "extraction_scheduler_memory_bytes"],
                "launcher_runtime_reserve_bytes": resolved_resources[
                    "extraction_launcher_runtime_reserve_bytes"],
                "worker_memory_budget_bytes": resolved_resources[
                    "extraction_worker_memory_budget_bytes"],
            })
    extraction_fields = (
        "task_index", "scientific_method_version", "calibration_library",
        "workload_generation_id", "cache_generation_id", "library", "modality",
        "candidate_manifest", "manifest_digest", "samples", "pileup_sites",
        "pileup_observations", "pileup_molecules", "cache_prefix", "marker",
        "input_status", "input_issues", "selected_cells",
        "projected_observation_records", "projected_molecule_records",
        "error_ref", "error_alt", "min_evidence", "max_second_fraction",
        "projected_peak_bytes", "bounded_memory_limit_bytes",
        "scheduler_memory_bytes", "launcher_runtime_reserve_bytes",
        "worker_memory_budget_bytes")
    write_tsv(targeted_root / "manifests" / "extraction_tasks.tsv",
              extraction_tasks, extraction_fields)

    extraction_by_key = {(row["library"], row["modality"]): row
                         for row in extraction_tasks}
    analysis_tasks = []
    target_roles = {(row["library"], row["barcode"]): "FROZEN_TARGET"
                    for row in targets}
    target_roles.update({(row["library"], row["barcode"]): "MATCHED_COMPARISON"
                         for row in matches})
    for (library, barcode), role in sorted(target_roles.items()):
        if (library, barcode) in unavailable_keys:
            continue
        cell = cell_by_key[(library, barcode)]
        for modality in ("RNA", "ATAC"):
            extraction = extraction_by_key[(library, modality)]
            site_records, molecule_records = _candidate_record_estimate(
                cell, modality.lower())
            menu_size = int(cell.get("menu_size", 0) or 0)
            analysis_tasks.append({
                "task_index": len(analysis_tasks),
                "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
                "calibration_library": CALIBRATION_LIBRARY,
                "workload_generation_id": generation_id,
                "cache_generation_id": generation_id,
                "action": "CELL", "library": library,
                "modality": modality, "role": role, "barcode": barcode,
                "control_ids": "", "cache_prefix": extraction["cache_prefix"],
                "candidate_manifest": extraction["candidate_manifest"],
                "control_manifest": str(final_controls),
                "analysis_output": str(targeted_root / "analysis" /
                    f"cell_{len(analysis_tasks):04d}.tsv.gz"),
                "marker": str(targeted_root / "markers" /
                    f"analysis_{len(analysis_tasks):04d}.complete.json"),
                "predicted_evidence_records": site_records + molecule_records,
                "predicted_likelihood_work": menu_size * (
                    (site_records + molecule_records) *
                    (1 + sum(DOWNSAMPLE_FRACTIONS) * DOWNSAMPLE_REPLICATES +
                     2 * CELL_NULL_REPLICATES)),
                "error_ref": extraction["error_ref"],
                "error_alt": extraction["error_alt"],
                "min_evidence": extraction["min_evidence"],
                "max_second_fraction": extraction["max_second_fraction"],
                "scheduler_memory_bytes": resolved_resources[
                    "analysis_scheduler_memory_bytes"],
                "launcher_runtime_reserve_bytes": resolved_resources[
                    "analysis_launcher_runtime_reserve_bytes"],
                "worker_memory_budget_bytes": resolved_resources[
                    "analysis_worker_memory_budget_bytes"],
            })

    pair_by_id = {row["control_id"]: row for row in control_pairs}
    for library in (f"lib{value}" for value in ALLOWED_LIBRARIES):
        library_pairs = [row for row in control_pairs if row["library"] == library]
        for modality in ("RNA", "ATAC"):
            weighted = []
            for pair in library_pairs:
                recipient = cell_by_key[(library, pair["recipient_barcode"])]
                source = cell_by_key[(library, pair["source_barcode"])]
                r_site, r_molecule = _candidate_record_estimate(
                    recipient, modality.lower())
                s_site, s_molecule = _candidate_record_estimate(
                    source, modality.lower())
                work = (int(recipient.get("menu_size", 0) or 0) *
                        (r_site + r_molecule + s_site + s_molecule) *
                        len(CONTROL_FRACTIONS))
                weighted.append((work, pair["control_id"]))
            weighted.sort(key=lambda item: (-item[0], item[1]))
            bins = []
            for work, control_id in weighted:
                if not bins or all(len(item["ids"]) >= 8 for item in bins):
                    bins.append({"work": 0, "ids": []})
                target_bin = min(
                    (item for item in bins if len(item["ids"]) < 8),
                    key=lambda item: (item["work"], len(item["ids"])))
                target_bin["ids"].append(control_id)
                target_bin["work"] += work
            if not bins:
                continue
            extraction = extraction_by_key[(library, modality)]
            for item in bins:
                analysis_tasks.append({
                    "task_index": len(analysis_tasks),
                    "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
                    "calibration_library": CALIBRATION_LIBRARY,
                    "workload_generation_id": generation_id,
                    "cache_generation_id": generation_id,
                    "action": "CONTROL_BIN", "library": library,
                    "modality": modality, "role": "SYNTHETIC_CONTROL",
                    "barcode": "", "control_ids": ",".join(item["ids"]),
                    "cache_prefix": extraction["cache_prefix"],
                    "candidate_manifest": extraction["candidate_manifest"],
                    "control_manifest": str(final_controls),
                    "analysis_output": str(targeted_root / "analysis" /
                        f"control_{len(analysis_tasks):04d}.tsv.gz"),
                    "marker": str(targeted_root / "markers" /
                        f"analysis_{len(analysis_tasks):04d}.complete.json"),
                    "predicted_evidence_records": "",
                    "predicted_likelihood_work": item["work"],
                    "error_ref": extraction["error_ref"],
                    "error_alt": extraction["error_alt"],
                    "min_evidence": extraction["min_evidence"],
                    "max_second_fraction": extraction["max_second_fraction"],
                    "scheduler_memory_bytes": resolved_resources[
                        "analysis_scheduler_memory_bytes"],
                    "launcher_runtime_reserve_bytes": resolved_resources[
                        "analysis_launcher_runtime_reserve_bytes"],
                    "worker_memory_budget_bytes": resolved_resources[
                        "analysis_worker_memory_budget_bytes"],
                })
    analysis_fields = (
        "task_index", "scientific_method_version", "calibration_library",
        "workload_generation_id", "cache_generation_id", "action", "library",
        "modality", "role", "barcode", "control_ids", "cache_prefix",
        "candidate_manifest", "control_manifest", "analysis_output", "marker",
        "predicted_evidence_records", "predicted_likelihood_work",
        "error_ref", "error_alt", "min_evidence", "max_second_fraction",
        "scheduler_memory_bytes", "launcher_runtime_reserve_bytes",
        "worker_memory_budget_bytes")
    write_tsv(targeted_root / "manifests" / "analysis_tasks.tsv",
              analysis_tasks, analysis_fields)
    scripts = _write_targeted_scripts(
        targeted_root, tool_bin_root, generation_id, extraction_tasks,
        analysis_tasks, benchmark=False,
        production_threads=resolved_resources["analysis_cpus"],
        resource_policy=resolved_resources)

    benchmark_root = targeted_root / "benchmark_unsubmitted"
    for name in ("manifests", "slurm_scripts", "logs", "task_scratch",
                 "analysis", "markers"):
        (benchmark_root / name).mkdir(parents=True, exist_ok=True)
    benchmark_extract = [dict(row) for row in extraction_tasks
                         if row["library"] == "lib12"]
    for index, row in enumerate(benchmark_extract):
        row["task_index"] = index
    write_tsv(benchmark_root / "manifests" / "extraction_tasks.tsv",
              benchmark_extract, extraction_fields)
    benchmark_controls = [row for row in control_pairs
                          if row["library"] == "lib12"]
    write_tsv(benchmark_root / "manifests" /
              "control_pairs_provisional.tsv", benchmark_controls)
    lib12_cell_tasks = [row for row in analysis_tasks
                        if row["library"] == "lib12" and row["action"] == "CELL"]
    branch_cells = []
    lib12_rows = [cell_by_key[(row["library"], row["barcode"])]
                  for row in lib12_cell_tasks]
    benchmark_branch_coverage = []

    def add_benchmark_cell_branch(label, candidates):
        if not candidates:
            benchmark_branch_coverage.append({"branch": label, "coverage_status": "UNAVAILABLE"})
            return
        selected = min(candidates, key=lambda row: row["barcode"])
        branch_cells.append(selected)
        benchmark_branch_coverage.append({
            "schema_version": "joint_doublet_benchmark_branch_v4",
            "workload_generation_id": generation_id,
            "branch": label, "eligible_cells": len(candidates),
            "selected_library": selected["library"],
            "selected_barcode": selected["barcode"],
            "coverage_status": "SELECTED",
        })

    add_benchmark_cell_branch(
        "STRICT_SITE_MIXTURE",
        [row for row in lib12_rows if target_id(row) in STRICT_SITE_TARGETS])
    add_benchmark_cell_branch(
        "UNCERTAIN_ADDITION",
        [row for row in lib12_rows if row.get(
            "site_existing_evidence_category") ==
            "addition-compatible but uncertain"])
    add_benchmark_cell_branch(
        "REPLACEMENT_OR_BOUNDARY",
        [row for row in lib12_rows if row.get(
            "site_existing_evidence_category") ==
            "replacement-like or upper-boundary fit"])
    for modality in ("rna", "atac"):
        field = f"{modality}_molecule_evidence_basis_class"
        for basis in sorted({clean(row.get(field, "")) for row in lib12_rows} -
                            {"", "UNAVAILABLE"}):
            add_benchmark_cell_branch(
                f"{modality.upper()}_MOLECULE_BASIS_{basis}",
                [row for row in lib12_rows if clean(row.get(field, "")) == basis])
    if not lib12_rows:
        raise RuntimeError("Library-12 benchmark has no target/comparison cells")
    smallest = min(lib12_rows, key=lambda row: (
        int(row.get("menu_size", 0) or 0), row["barcode"]))
    largest = max(lib12_rows, key=lambda row: (
        int(row.get("menu_size", 0) or 0), row["barcode"]))
    for label, selected in (("SMALLEST_COMPLETE_MENU", smallest),
                            ("LARGEST_COMPLETE_MENU", largest)):
        branch_cells.append(selected)
        benchmark_branch_coverage.append({
            "schema_version": "joint_doublet_benchmark_branch_v4",
            "workload_generation_id": generation_id,
            "branch": label, "eligible_cells": len(lib12_rows),
            "selected_library": selected["library"],
            "selected_barcode": selected["barcode"],
            "selected_menu_size": selected.get("menu_size", 0),
            "coverage_status": "SELECTED",
        })
    branch_keys = {(row["library"], row["barcode"]) for row in branch_cells}
    for role in ("FROZEN_TARGET", "MATCHED_COMPARISON"):
        candidates = [row for row in lib12_cell_tasks if row["role"] == role]
        if not candidates:
            raise RuntimeError(
                f"Library-12 benchmark lacks required {role} shard")
        highest = max(candidates, key=lambda row: (
            float(row["predicted_likelihood_work"]), row["barcode"]))
        branch_keys.add((highest["library"], highest["barcode"]))
        benchmark_branch_coverage.append({
            "schema_version": "joint_doublet_benchmark_branch_v4",
            "workload_generation_id": generation_id,
            "branch": role, "eligible_cells": len(candidates),
            "selected_library": highest["library"],
            "selected_barcode": highest["barcode"],
            "coverage_status": "SELECTED",
        })
    benchmark_tasks = [dict(row) for row in lib12_cell_tasks
                       if (row["library"], row["barcode"]) in branch_keys]
    control_candidates = [row for row in analysis_tasks
                          if row["library"] == "lib12" and
                          row["action"] == "CONTROL_BIN"]
    selected_control = max(control_candidates, key=lambda row: (
        float(row["predicted_likelihood_work"] or 0), row["task_index"])) if control_candidates else {}
    if selected_control:
        benchmark_tasks.append(dict(selected_control))
    benchmark_branch_coverage.append({
        "schema_version": "joint_doublet_benchmark_branch_v4",
        "workload_generation_id": generation_id,
        "branch": "HIGHEST_WORK_CONTROL_SHARD",
        "eligible_cells": len(control_candidates),
        "selected_library": selected_control.get("library", ""),
        "selected_barcode": "CONTROL_BIN",
        "selected_control_ids": selected_control.get("control_ids", ""),
        "coverage_status": "SELECTED",
    })
    for index, row in enumerate(benchmark_tasks):
        row["task_index"] = index
        row["control_manifest"] = str(
            benchmark_root / "manifests" / "control_pairs_final.tsv")
        row["analysis_output"] = str(benchmark_root / "analysis" /
            f"benchmark_{index:03d}.tsv.gz")
        row["marker"] = str(benchmark_root / "markers" /
            f"analysis_{index:03d}.complete.json")
    write_tsv(benchmark_root / "manifests" / "analysis_tasks.tsv",
              benchmark_tasks, analysis_fields)
    benchmark_scripts = _write_targeted_scripts(
        benchmark_root, tool_bin_root, generation_id, benchmark_extract,
        benchmark_tasks, benchmark=True,
        production_threads=resolved_resources["analysis_cpus"],
        resource_policy=resolved_resources)
    write_tsv(benchmark_root / "manifests" / "benchmark_panel.tsv",
              benchmark_tasks)
    write_tsv(benchmark_root / "manifests" /
              "benchmark_branch_coverage.tsv", benchmark_branch_coverage)
    write_tsv(benchmark_root / "frozen_targets_20260920.tsv", load_frozen_targets(frozen_manifest)[1])
    legacy_reference_rows = []
    for task in benchmark_tasks:
        if task["action"] != "CELL":
            continue
        cell = cell_by_key[(task["library"], task["barcode"])]
        modality = task["modality"].lower()
        for candidate in cell.get("_candidates", []):
            for evidence in ("site", "molecule"):
                metric = candidate[f"{modality}_{evidence}"]
                legacy_reference_rows.append({
                    "workload_generation_id": generation_id,
                    "library": task["library"],
                    "modality": task["modality"],
                    "barcode": task["barcode"],
                    "candidate_id": candidate["candidate_id"],
                    "evidence_channel": evidence.upper(),
                    "legacy_status": metric["status"],
                    "legacy_locked_log_likelihood": metric["locked"],
                    "legacy_interior_log_likelihood": metric["interior"],
                    "legacy_contributor_only_log_likelihood":
                        metric["contributor"],
                    "legacy_delta_log_likelihood": metric["delta"],
                    "legacy_fitted_fraction": metric["fraction"],
                    "legacy_fitted_fraction_profile_low":
                        metric["profile_low"],
                    "legacy_fitted_fraction_profile_high":
                        metric["profile_high"],
                    "legacy_evidence_units": metric["units"],
                    "source": "COMPLETED_JOINT_DOUBLET_CANDIDATE_SCORE_ROW",
                })
    write_tsv(benchmark_root / "manifests" /
              "benchmark_completed_score_reference.tsv.gz",
              legacy_reference_rows)

    expanded_rows = sum(manifest_counts[(f"lib{number}", modality, key)]
                        for number in ALLOWED_LIBRARIES
                        for modality in ("rna", "atac")
                        for key in ("site_cache_rows", "molecule_cache_rows"))
    unique_records = sum(row["projected_observation_records"] +
                         row["projected_molecule_records"]
                         for row in extraction_tasks)
    accounting = [{
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "workload_generation_id": generation_id,
        "selected_target_cells": len(targets),
        "requested_comparisons": len(targets),
        "matched_comparison_cells": len(matches),
        "unmatched_target_cells": len(targets) - len(matches),
        "selected_unique_cells_including_control_parents": len(selected_keys),
        "extraction_libraries": ",".join(sorted({row["library"] for row in extraction_tasks})),
        "full_source_scans": len(extraction_tasks),
        "full_source_scan_definition": "one per target-library/modality pair",
        "candidate_support_record_upper_bound": unique_records,
        "unique_evidence_records": "MEASURED_DURING_EXTRACTION",
        "unique_units": "MEASURED_DURING_EXTRACTION",
        "legacy_candidate_expanded_rows_avoided": expanded_rows,
        "projected_reduction_factor": expanded_rows / unique_records
            if unique_records else math.nan,
        "normalized_cache_support_bytes_projection": unique_records * 40,
        "normalized_cache_total_bytes": "UNKNOWN_UNTIL_EXTRACTION_INCLUDES_NONDISCRIMINATING_EVIDENCE_AND_DICTIONARY",
        "compressed_cache_bytes": "NOT_COMPRESSED_NATIVE_RECORDS;GENOTYPE_DICTIONARY_MEASURED",
        "temporary_storage_bound": "extraction numeric records <=3*selected pre-deduplication fixed-width bytes, plus dictionary/output files; measured during extraction. Analysis spills compact candidate coefficients and whole-unit metadata; observed fields are interned once",
        "analysis_metadata_memory_bound": "(observations+molecule_rows)*(2048+32*menu+128*workers)+8*menu*spill_limit+128MiB; dictionaries additional",
        "extraction_sort_memory_bound_bytes": 256*1024*1024,
        "replicate_block_size": resolved_resources["analysis_cpus"],
        "gather_task_hours": "UNKNOWN_UNMEASURED_THROUGHPUT",
        "extraction_task_hours": "UNKNOWN_UNMEASURED_THROUGHPUT",
        "analysis_task_hours": "UNKNOWN_UNMEASURED_THROUGHPUT",
        "normalized_cache_projection_formula":
            "candidate support sum * 40 bytes; not an estimate of the unique raw observation universe",
        "target_comparison_tasks": sum(row["action"] == "CELL"
                                       for row in analysis_tasks),
        "control_bin_tasks": sum(row["action"] == "CONTROL_BIN"
                                 for row in analysis_tasks),
        "control_unique_pairs": len(control_pairs),
        "control_fractions": ",".join(map(str, CONTROL_FRACTIONS)),
        "downsample_replicates_per_fraction": DOWNSAMPLE_REPLICATES,
        "site_null_replicates_per_cell": CELL_NULL_REPLICATES,
        "molecule_null_replicates_per_cell": CELL_NULL_REPLICATES,
        "extraction_cpus_each": resolved_resources["extraction_cpus"],
        "extraction_memory_each": resolved_resources["extraction_memory"],
        "analysis_cpus_each": resolved_resources["analysis_cpus"],
        "analysis_memory_each": resolved_resources["analysis_memory"],
        "extraction_wall_time": resolved_resources["extraction_wall_time"],
        "analysis_wall_time": resolved_resources["analysis_wall_time"],
        "partition": resolved_resources["partition"],
        "prebenchmark_concurrency": 4,
        "concurrency_basis": "CONSERVATIVE_PREBENCHMARK_ASSUMPTION",
        "projected_task_hours": "UNKNOWN_UNMEASURED_THROUGHPUT",
        "projected_core_hours": "UNKNOWN_UNMEASURED_THROUGHPUT",
        "projected_throttled_critical_path":
            "UNKNOWN_UNMEASURED_THROUGHPUT",
        "queue_time": "EXCLUDED",
        "memory_safety_margin": 1.50,
        "submission_status": "RENDERED_UNSUBMITTED",
        "invalid_legacy_0.312_hour_estimate":
            "AUDIT_ONLY;REMOVED_FROM_WORKLOAD_PROJECTION",
    }]
    write_tsv(targeted_root / "workload_accounting.tsv", accounting)
    benchmark_projection = [{
        **accounting[0],
        "scope": "BOUNDED_LIBRARY12_PRODUCTION_CODE_BENCHMARK",
        "selected_target_cells": len({row.get("barcode", "") for row in
                                      benchmark_tasks
                                      if row.get("action") == "CELL"}),
        "extraction_libraries": "12",
        "full_source_scans": len(benchmark_extract),
        "target_comparison_tasks": sum(
            row.get("action") == "CELL" for row in benchmark_tasks),
        "control_bin_tasks": sum(
            row.get("action") == "CONTROL_BIN" for row in benchmark_tasks),
        "submission_status": "RENDERED_UNSUBMITTED",
    }]
    write_tsv(benchmark_root / "workload_accounting.tsv",
              benchmark_projection)
    provisional_projection = {
        "schema_version": "joint_doublet_resource_projection_v4",
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "workload_generation_id": generation_id,
        "projection_status": "PREBENCHMARK_PROVISIONAL",
        "submission_status": "RENDERED_UNSUBMITTED",
        "resource_policy": resolved_resources,
        "array_concurrency": 4,
        "memory_safety_margin": 1.50,
        "maximum_extraction_projected_peak_bytes": max((int(row.get(
            "projected_peak_bytes", 0) or 0) for row in extraction_tasks),
            default=0),
        "projected_task_hours": "UNKNOWN_UNMEASURED_THROUGHPUT",
        "projected_core_hours": "UNKNOWN_UNMEASURED_THROUGHPUT",
        "projected_throttled_critical_path":
            "UNKNOWN_UNMEASURED_THROUGHPUT",
        "queue_time": "EXCLUDED",
    }
    atomic_json(targeted_root / "resource_projection.json",
                provisional_projection)
    atomic_json(benchmark_root / "resource_projection.json", {
        **provisional_projection,
        "scope": "BOUNDED_LIBRARY12_PRODUCTION_CODE_BENCHMARK",
        "maximum_extraction_projected_peak_bytes": max((int(row.get(
            "projected_peak_bytes", 0) or 0) for row in benchmark_extract),
            default=0),
    })
    blueprint = {
        "schema_version": WORKLOAD_SCHEMA,
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "workload_generation_id": generation_id,
        "state": "RENDERED_UNSUBMITTED",
        "created_utc": utc_now(), "seed": SEED,
        "target_libraries": sorted({int(row["library"][3:]) for row in target_manifest}),
        "unavailable_cells": unavailable_cells,
        "protected_libraries": list(PROTECTED_LIBRARIES),
        "extraction_tasks": extraction_tasks,
        "analysis_tasks": analysis_tasks,
        "scripts": scripts,
        "benchmark_root": str(benchmark_root),
        "benchmark_scripts": benchmark_scripts,
        "resource_policy": {**accounting[0], **resolved_resources},
        "control_finalization": (
            "provisional assignments become executable only after both modality "
            "genotype dictionaries validate and decoy tolerances pass"),
    }
    atomic_json(targeted_root / "workload_blueprint.json", blueprint)
    atomic_json(benchmark_root / "workload_blueprint.json", {
        "schema_version": WORKLOAD_SCHEMA,
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "workload_generation_id": generation_id,
        "state": "RENDERED_UNSUBMITTED",
        "created_utc": utc_now(), "seed": SEED,
        "scope": "BOUNDED_LIBRARY12_BENCHMARK",
        "parent_workload_root": str(targeted_root.resolve()),
        "extraction_tasks": benchmark_extract,
        "analysis_tasks": benchmark_tasks,
        "scripts": benchmark_scripts,
        "resource_policy": resolved_resources,
        "required_measurements": {
            "extraction_rows": len(benchmark_extract),
            "analysis_rows_per_thread_configuration": len(benchmark_tasks),
            "thread_configurations": [1, resolved_resources["analysis_cpus"]],
            "requires_target_and_comparison": True,
            "requires_control": True,
            "requires_compiled_scalar_equivalence": True,
        },
    })
    atomic_json(targeted_root / "CURRENT_GENERATION.json", {
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "workload_generation_id": generation_id,
        "root": str(targeted_root), "state": "RENDERED_UNSUBMITTED",
    })
    atomic_json(targeted_root / "TARGETED_WORKLOAD_RENDERED_UNSUBMITTED", {
        "schema_version": WORKLOAD_SCHEMA,
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "status": "RENDERED_UNSUBMITTED", "utc": utc_now(),
        "workload_generation_id": generation_id,
        "input_status": "READY" if all(row["input_status"] == "READY"
                                       for row in extraction_tasks) else "INCOMPLETE",
        "extraction_tasks": len(extraction_tasks),
        "analysis_tasks": len(analysis_tasks),
        "manifest_libraries": list(ALLOWED_LIBRARIES),
        "protected_libraries_accessed": sorted(
            set(ALLOWED_LIBRARIES) & set(PROTECTED_LIBRARIES)),
    })
    atomic_json(benchmark_root / "BENCHMARK_RENDERED_UNSUBMITTED", {
        "schema_version": WORKLOAD_SCHEMA,
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "status": "RENDERED_UNSUBMITTED", "utc": utc_now(),
        "workload_generation_id": generation_id,
        "complete_library12_cache": True,
        "immutable_cache_reuse_required": True,
        "analysis_tasks": len(benchmark_tasks),
    })
    orchestrator = Path(tool_bin_root) / "orchestrate_tetraploid.py"
    common = ("module purge\nmodule load miniforge/3\n" +
              f"python3 -B {shlex.quote(str(orchestrator))} "
              f"--stage JOINT_DOUBLET --libraries 7 9 12 17 20 25 29 "
              f"--joint-doublet-targeted-output-root {shlex.quote(str(targeted_root))} "
              f"--joint-doublet-tool-bin-root {shlex.quote(str(tool_bin_root))}")
    commands = {
        "benchmark_launch": common +
            " --joint-doublet-action TARGETED_BENCHMARK_LAUNCH --submit",
        "benchmark_status": common +
            " --joint-doublet-action TARGETED_BENCHMARK_STATUS",
        "benchmark_resume": common +
            " --joint-doublet-action TARGETED_BENCHMARK_RESUME --submit",
        "reproject": common + " --joint-doublet-action TARGETED_REPROJECT",
        "full_launch": common +
            " --joint-doublet-action TARGETED_LAUNCH --submit",
        "full_status": common + " --joint-doublet-action TARGETED_STATUS",
        "full_resume": common +
            " --joint-doublet-action TARGETED_RESUME --submit",
    }
    atomic_json(targeted_root / "exact_commands.json", commands)
    atomic_json(benchmark_root / "exact_commands.json", commands)
    return accounting[0], blueprint




def submit_sbatch(script, root, dependency="", array="", index_map=""):
    command = ["sbatch", "--parsable", f"--chdir={root}"]
    if dependency:
        command.append(f"--dependency={dependency}")
    if array:
        command.append(f"--array={array}")
    if index_map:
        command.append(f"--export=ALL,JOINT_DOUBLET_TASK_INDEX_MAP={index_map}")
    command.append(str(script))
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or result.stdout.strip())
    job_id = result.stdout.strip().split(";", 1)[0]
    if not job_id.isdigit():
        raise RuntimeError(f"unexpected sbatch response: {result.stdout!r}")
    return job_id




def _full_gzip_audit(path):
    digest = hashlib.sha256()
    rows = 0
    header = []
    with gzip.open(path, "rb") as handle:
        for line in handle:
            digest.update(line)
            if rows == 0:
                header = line.decode("utf-8").rstrip("\r\n").split("\t")
            rows += 1
    return header, max(rows - 1, 0), digest.hexdigest()


def validate_targeted_manifest_boundary(rows, task_kind):
    allowed = {f"lib{value}" for value in ALLOWED_LIBRARIES}
    for row in rows:
        library = clean(row.get("library", ""))
        modality = clean(row.get("modality", ""))
        if library not in allowed:
            raise RuntimeError(
                f"{task_kind} manifest rejected unexpected library {library!r}")
        if modality not in {"RNA", "ATAC"}:
            raise RuntimeError(
                f"{task_kind} manifest rejected modality {modality!r}")
        if task_kind == "analysis" and row.get("action") not in {
                "CELL", "CONTROL_BIN"}:
            raise RuntimeError(
                f"analysis manifest rejected action {row.get('action')!r}")
        for field, value in row.items():
            if value and ("path" in field or field in {
                    "candidate_manifest", "control_manifest", "cache_prefix",
                    "analysis_output", "marker", "samples", "pileup_sites",
                    "pileup_observations", "pileup_molecules"}) and \
                    _path_mentions_protected(value):
                raise RuntimeError(
                    f"{task_kind} manifest path names a protected library: "
                    f"{field}={value}")


def validate_targeted_input_contract(task_manifest, task_kind, task_index):
    """Validate every generated manifest before any biological source open."""
    manifest_path = Path(task_manifest).resolve(strict=True)
    if manifest_path.name != f"{task_kind}_tasks.tsv" or \
            manifest_path.parent.name != "manifests":
        raise RuntimeError("targeted task manifest has an unexpected location")
    root = manifest_path.parent.parent.resolve(strict=True)
    try:
        manifest_path.relative_to(root)
    except ValueError as error:
        raise RuntimeError("targeted task manifest escapes workload root") from error
    tasks = list(read_tsv(manifest_path))
    if not tasks:
        raise RuntimeError("targeted task manifest is empty")
    required_task_fields = {
        "task_index", "scientific_method_version", "calibration_library",
        "workload_generation_id", "cache_generation_id",
        "library", "modality", "candidate_manifest", "cache_prefix",
        "marker", "error_ref", "error_alt", "min_evidence",
        "max_second_fraction", "scheduler_memory_bytes",
        "launcher_runtime_reserve_bytes", "worker_memory_budget_bytes",
    }
    if task_kind == "extraction":
        required_task_fields.update({
            "manifest_digest", "samples", "pileup_sites",
            "pileup_observations", "pileup_molecules", "input_status",
            "selected_cells", "projected_observation_records",
            "projected_molecule_records", "projected_peak_bytes",
            "bounded_memory_limit_bytes",
        })
    else:
        required_task_fields.update({
            "action", "role", "barcode", "control_ids",
            "control_manifest", "analysis_output",
            "predicted_likelihood_work",
        })
    if not required_task_fields.issubset(tasks[0]):
        raise RuntimeError(
            "targeted task manifest schema missing: " + ",".join(sorted(
                required_task_fields - set(tasks[0]))))
    validate_targeted_manifest_boundary(tasks, task_kind)
    indices = [int(row.get("task_index", -1)) for row in tasks]
    if indices != list(range(len(tasks))) or len(indices) != len(set(indices)):
        raise RuntimeError("targeted task indices must be unique and contiguous")
    generations = {clean(row.get("workload_generation_id", "")) for row in tasks}
    if len(generations) != 1 or "" in generations:
        raise RuntimeError("targeted manifest mixes workload generations")
    if task_index < 0 or task_index >= len(tasks):
        raise RuntimeError("task index outside targeted manifest")
    generation = next(iter(generations))
    if not re.fullmatch(r"(?:prebenchmark|reprojected)_[0-9a-f]{16}",
                        generation):
        raise RuntimeError("targeted manifest generation identifier is invalid")
    blueprint_path = root / "workload_blueprint.json"
    if not blueprint_path.is_file():
        raise RuntimeError("targeted workload blueprint is missing")
    blueprint = json.loads(blueprint_path.read_text())
    if blueprint.get("workload_generation_id") != generation:
        raise RuntimeError("workload blueprint generation mismatch")
    if blueprint.get("scientific_method_version") != \
            SCIENTIFIC_METHOD_VERSION or int(blueprint.get(
                "calibration_library", -1)) != CALIBRATION_LIBRARY:
        raise RuntimeError("workload blueprint method/calibration mismatch")
    for row in tasks:
        if row.get("scientific_method_version") != \
                SCIENTIFIC_METHOD_VERSION or row.get(
                    "calibration_library") != str(CALIBRATION_LIBRARY):
            raise RuntimeError("targeted task method/calibration mismatch")
        scheduler = int(row.get("scheduler_memory_bytes", 0) or 0)
        reserve = int(row.get("launcher_runtime_reserve_bytes", 0) or 0)
        worker = int(row.get("worker_memory_budget_bytes", 0) or 0)
        if scheduler <= 0 or reserve != max(
                16 * (1 << 20), (scheduler + 9) // 10) or \
                worker != scheduler-reserve or not 0 < worker <= scheduler:
            raise RuntimeError("targeted task memory contract mismatch")
    blueprint_key = "extraction_tasks" if task_kind == "extraction" else \
        "analysis_tasks"
    declared_tasks = blueprint.get(blueprint_key, [])
    if len(declared_tasks) != len(tasks):
        raise RuntimeError("task manifest count differs from workload blueprint")
    for observed, declared in zip(tasks, declared_tasks):
        normalized_declared = {key: str(value) for key, value in declared.items()}
        if not set(normalized_declared).issubset(observed) or any(
                observed.get(key, "") != value
                for key, value in normalized_declared.items()):
            raise RuntimeError(
                f"task manifest content differs from workload blueprint at "
                f"index {observed.get('task_index')}")
    expected_libraries = {"lib12"} if root.name == "benchmark_unsubmitted" \
        else {row["library"] for row in blueprint["extraction_tasks"]}
    observed_pairs = {(row["library"], row["modality"]) for row in tasks}
    if {library for library, _modality in observed_pairs} != expected_libraries or \
            {modality for _library, modality in observed_pairs} != {"RNA", "ATAC"}:
        raise RuntimeError(
            "targeted task manifest does not cover every required library/modality")
    if task_kind == "extraction" and observed_pairs != {
            (library, modality) for library in expected_libraries
            for modality in ("RNA", "ATAC")}:
        raise RuntimeError(
            "extraction task manifest must contain exactly one row per "
            "required library/modality")

    # Validate the companion task manifest as part of the same input boundary.
    # Extraction workers therefore cannot begin source discovery after checking
    # only their own row while a stale/mixed analysis manifest is present.
    companion_kind = "analysis" if task_kind == "extraction" else "extraction"
    companion_path = root / "manifests" / f"{companion_kind}_tasks.tsv"
    companion_tasks = list(read_tsv(companion_path))
    if not companion_tasks:
        raise RuntimeError("targeted companion task manifest is empty")
    companion_required = {
        "task_index", "scientific_method_version", "calibration_library",
        "workload_generation_id", "cache_generation_id",
        "library", "modality", "candidate_manifest", "cache_prefix",
        "marker", "error_ref", "error_alt", "min_evidence",
        "max_second_fraction", "scheduler_memory_bytes",
        "launcher_runtime_reserve_bytes", "worker_memory_budget_bytes",
    }
    if companion_kind == "extraction":
        companion_required.update({
            "manifest_digest", "samples", "pileup_sites",
            "pileup_observations", "pileup_molecules", "input_status",
            "selected_cells", "projected_observation_records",
            "projected_molecule_records", "projected_peak_bytes",
            "bounded_memory_limit_bytes",
        })
    else:
        companion_required.update({
            "action", "role", "barcode", "control_ids",
            "control_manifest", "analysis_output",
            "predicted_likelihood_work",
        })
    if not companion_required.issubset(companion_tasks[0]):
        raise RuntimeError(
            "targeted companion task manifest schema missing: " +
            ",".join(sorted(companion_required - set(companion_tasks[0]))))
    validate_targeted_manifest_boundary(companion_tasks, companion_kind)
    companion_indices = [int(row.get("task_index", -1))
                         for row in companion_tasks]
    if companion_indices != list(range(len(companion_tasks))) or \
            len(companion_indices) != len(set(companion_indices)):
        raise RuntimeError(
            "targeted companion task indices must be unique and contiguous")
    if {clean(row.get("workload_generation_id", ""))
            for row in companion_tasks} != {generation}:
        raise RuntimeError("targeted companion manifest generation mismatch")
    for row in companion_tasks:
        if row.get("scientific_method_version") != \
                SCIENTIFIC_METHOD_VERSION or row.get(
                    "calibration_library") != str(CALIBRATION_LIBRARY):
            raise RuntimeError(
                "targeted companion method/calibration mismatch")
        scheduler = int(row.get("scheduler_memory_bytes", 0) or 0)
        reserve = int(row.get("launcher_runtime_reserve_bytes", 0) or 0)
        worker = int(row.get("worker_memory_budget_bytes", 0) or 0)
        if scheduler <= 0 or reserve != max(
                16 * (1 << 20), (scheduler + 9) // 10) or \
                worker != scheduler-reserve or not 0 < worker <= scheduler:
            raise RuntimeError("targeted companion memory contract mismatch")
    companion_declared = blueprint.get(
        "analysis_tasks" if companion_kind == "analysis" else
        "extraction_tasks", [])
    if len(companion_declared) != len(companion_tasks):
        raise RuntimeError(
            "companion task manifest count differs from workload blueprint")
    for observed, declared in zip(companion_tasks, companion_declared):
        if not set(declared).issubset(observed) or any(
                observed.get(key, "") != str(value)
                for key, value in declared.items()):
            raise RuntimeError(
                "companion task manifest content differs from workload blueprint")
    all_extraction_tasks = tasks if task_kind == "extraction" else \
        companion_tasks
    all_analysis_tasks = tasks if task_kind == "analysis" else companion_tasks
    if {(row["library"], row["modality"])
            for row in all_extraction_tasks} != {
            (library, modality) for library in expected_libraries
            for modality in ("RNA", "ATAC")}:
        raise RuntimeError(
            "companion extraction manifest lacks exact library/modality coverage")
    if {row["library"] for row in all_analysis_tasks} != expected_libraries or \
            {row["modality"] for row in all_analysis_tasks} != {"RNA", "ATAC"}:
        raise RuntimeError(
            "analysis task manifest lacks required library/modality coverage")
    allowed_root = root.parent.resolve() if root.name == \
        "benchmark_unsubmitted" else root.resolve()
    workload_container = root.parent.parent.resolve() if \
        root.parent.name == "generations" else allowed_root
    extraction_contract_by_pair = {
        (row["library"], row["modality"]): row
        for row in all_extraction_tasks}
    generated_path_fields = {"cache_prefix", "marker", "analysis_output",
                             "candidate_manifest", "control_manifest"}
    for manifest_kind, manifest_rows in (
            (task_kind, tasks), (companion_kind, companion_tasks)):
        for row in manifest_rows:
            for field in generated_path_fields & set(row):
                value = clean(row.get(field, ""))
                if not value:
                    continue
                path = Path(value)
                if not path.is_absolute():
                    raise RuntimeError(
                        f"targeted generated path is not absolute: {field}")
                field_root = allowed_root
                extraction_contract = extraction_contract_by_pair.get(
                    (row.get("library", ""), row.get("modality", "")), {})
                if field == "cache_prefix" and extraction_contract.get(
                        "cache_policy") == "IMMUTABLE_BENCHMARK_REUSE_REQUIRED":
                    field_root = workload_container
                try:
                    path.resolve(strict=False).relative_to(field_root)
                except ValueError as error:
                    raise RuntimeError(
                        f"targeted generated path escapes workload root: {field}") \
                        from error
            if manifest_kind == "extraction":
                cache_generation = clean(row.get("cache_generation_id", ""))
                immutable_reuse = row.get("cache_policy") == \
                    "IMMUTABLE_BENCHMARK_REUSE_REQUIRED"
                if cache_generation != generation and not (
                        immutable_reuse and re.fullmatch(
                            r"prebenchmark_[0-9a-f]{16}", cache_generation)):
                    raise RuntimeError("extraction cache generation is not bound")
                for field in ("samples", "pileup_sites", "pileup_observations",
                              "pileup_molecules"):
                    value = clean(row.get(field, ""))
                    if not value or not Path(value).is_absolute() or \
                            _path_mentions_protected(value):
                        raise RuntimeError(
                            f"invalid extraction source path contract: {field}")
                if int(row.get("min_evidence", -1)) < 0 or not 0 < finite(
                        row.get("max_second_fraction", "")) <= 1 or not (
                        0 <= finite(row.get("error_ref", "")) <= 1 and
                        0 <= finite(row.get("error_alt", "")) <= 1 and
                        finite(row.get("error_ref", "")) + finite(
                            row.get("error_alt", "")) < 1):
                    raise RuntimeError(
                        "invalid extraction model parameter contract")
                if int(row.get("projected_peak_bytes", 0) or 0) > int(
                        row.get("bounded_memory_limit_bytes", 0) or 0):
                    raise RuntimeError(
                        "extraction projected peak exceeds memory bound")
    frozen_path = root / "frozen_targets_20260920.tsv"
    if not frozen_path.is_file() and root.name == "benchmark_unsubmitted":
        frozen_path = root.parent / "frozen_targets_20260920.tsv"
    if frozen_path.is_file():
        _frozen, frozen_rows = load_frozen_targets(frozen_path)
        frozen_counts = Counter(int(row["library"]) for row in frozen_rows)
        if frozen_counts != Counter({12: 32, 20: 26, 29: 1}):
            raise RuntimeError("frozen target distribution mismatch")
    else:
        raise RuntimeError("targeted workload lacks frozen target manifest")
    candidate_paths = sorted({Path(row["candidate_manifest"]).absolute()
                              for row in all_extraction_tasks
                              if row.get("candidate_manifest")})
    candidate_rows_by_path = {}
    candidate_semantics_by_pair = {}
    for candidate_path in candidate_paths:
        rows = list(read_tsv(candidate_path))
        required = {"schema_version", "scientific_method_version",
                    "calibration_library", "library", "barcode", "candidate_id",
                    "locked_state", "locked_copy_vector", "second_state",
                    "second_copy_vector", "candidate_policy",
                    "physical_pool_state", "component_only_state",
                    "structural_added_state_relationship",
                    "candidate_origin", "exhaustive_fallback",
                    "nomination_modalities"}
        if not rows or not required.issubset(rows[0]):
            raise RuntimeError(f"candidate manifest schema invalid: {candidate_path}")
        keys = []
        for row in rows:
            if row.get("schema_version") != \
                    "joint_doublet_candidate_manifest_v4" or row.get(
                        "scientific_method_version") != \
                    SCIENTIFIC_METHOD_VERSION or row.get(
                        "calibration_library") != str(CALIBRATION_LIBRARY):
                raise RuntimeError(
                    "candidate manifest method/calibration provenance mismatch")
            if row.get("library") not in {f"lib{x}" for x in ALLOWED_LIBRARIES}:
                raise RuntimeError("candidate manifest contains invalid library")
            if not re.fullmatch(r"[ACGT]{16}", row.get("barcode", "")):
                raise RuntimeError("candidate manifest contains invalid barcode")
            if row.get("candidate_policy") != "DERIVATIVE_COMPLETE":
                raise RuntimeError("candidate manifest policy is not derivative-complete")
            keys.append((row["library"], row["barcode"], row["candidate_id"]))
        if len(keys) != len(set(keys)):
            raise RuntimeError("candidate manifest has duplicate candidate keys")
        by_cell = defaultdict(list)
        for row in rows:
            by_cell[(row["library"], row["barcode"])].append(row)
        if any(not values for values in by_cell.values()):
            raise RuntimeError("candidate manifest contains an empty legal menu")
        candidate_rows_by_path[candidate_path.resolve()] = rows
    for row in all_extraction_tasks:
        candidate_path = Path(row["candidate_manifest"]).resolve()
        semantic_fields = (
            "schema_version", "library", "barcode", "candidate_id",
            "locked_state", "locked_copy_vector", "second_state",
            "second_copy_vector", "candidate_origin", "exhaustive_fallback",
            "nomination_modalities", "candidate_policy",
            "physical_pool_state", "component_only_state",
            "structural_added_state_relationship",
        )
        candidate_semantics_by_pair[(row["library"], row["modality"])] = {
            tuple(clean(candidate.get(field, "")) for field in semantic_fields)
            for candidate in candidate_rows_by_path[candidate_path]
        }
    for library in expected_libraries:
        if candidate_semantics_by_pair[(library, "RNA")] != \
                candidate_semantics_by_pair[(library, "ATAC")]:
            raise RuntimeError(
                f"RNA/ATAC complete candidate menus differ for {library}")
    target_manifest = root / "target_cells_and_matched_comparisons.tsv"
    if not target_manifest.is_file() and root.name == "benchmark_unsubmitted":
        target_manifest = root.parent / "target_cells_and_matched_comparisons.tsv"
    target_rows = list(read_tsv(target_manifest))
    original_path = root / "frozen_target_comparisons.tsv"
    if not original_path.is_file() and root.name == "benchmark_unsubmitted":
        original_path = root.parent / "frozen_target_comparisons.tsv"
    original = load_frozen_comparisons(original_path, {(f"lib{row['library']}", row['barcode']) for row in frozen_rows})
    by_original = {row["target_id"]: row for row in original}
    for row in target_rows:
        source = by_original.get(row.get("target_id"), {})
        if any(row.get(key, "") != value for key, value in source.items()):
            raise RuntimeError("original target/comparison mapping or historical provenance changed")
    target_required = {"target_id", "library", "target_barcode", "locked_identity",
                       "proposed_contributor", "matched_comparison_barcode", "selection_rule"}
    if not target_rows or not target_required.issubset(target_rows[0]) or \
            len({row.get("target_id") for row in target_rows}) != \
            len(target_rows):
        raise RuntimeError("target/comparison manifest is invalid")
    frozen_ids = {row["target_id"] for row in frozen_rows}
    target_ids = {row.get("target_id", "") for row in target_rows}
    if target_ids != frozen_ids:
        raise RuntimeError("target/comparison manifest target identities mismatch")
    selected_by_library = defaultdict(set)
    comparison_keys = set()
    for row in target_rows:
        library = row.get("library", "")
        target_barcode = row.get("target_barcode", "")
        if library not in {f"lib{value}" for value in ALLOWED_LIBRARIES} or \
                not re.fullmatch(r"[ACGT]{16}", target_barcode) or \
                row.get("target_id") != f"{library}:{target_barcode}":
            raise RuntimeError("target/comparison manifest has an invalid target key")
        selected_by_library[library].add(target_barcode)
        comparison = clean(row.get("matched_comparison_barcode", ""))
        if comparison:
            if not re.fullmatch(r"[ACGT]{16}", comparison):
                raise RuntimeError(
                    "target/comparison manifest has an invalid matched comparison")
            selected_by_library[library].add(comparison)
            if (library, comparison) in comparison_keys or \
                    f"{library}:{comparison}" in frozen_ids:
                raise RuntimeError(
                    "target/comparison manifest reuses a comparison or target")
            comparison_keys.add((library, comparison))
        elif row.get("match_status") == "MATCHED":
            raise RuntimeError(
                "target/comparison manifest marks an empty comparison as matched")
    control_paths = {Path(row["control_manifest"]).absolute() for row in tasks
                     if row.get("action") == "CONTROL_BIN" and
                     row.get("control_manifest")}
    provisional = root / "manifests" / "control_pairs_provisional.tsv"
    if provisional.is_file():
        control_paths.add(provisional)
    for control_path in control_paths:
        if not control_path.is_file():
            # The finalized control manifest is expected only after extraction;
            # provisional content still validates the complete planned controls.
            if control_path.name == "control_pairs_final.tsv" and provisional.is_file():
                continue
            raise RuntimeError(f"control manifest missing: {control_path}")
        rows = list(read_tsv(control_path))
        required = {"schema_version", "scientific_method_version",
                    "calibration_library", "control_id", "control_class", "library",
                    "recipient_barcode", "source_barcode", "source_identity",
                    "expected_contributor", "requested_fractions", "seed"}
        if rows and not required.issubset(rows[0]):
            raise RuntimeError("control manifest schema invalid")
        identifiers = [row.get("control_id", "") for row in rows]
        if any(not value for value in identifiers) or \
                len(identifiers) != len(set(identifiers)):
            raise RuntimeError("control manifest identifiers are invalid")
        if any(row.get("library") not in {f"lib{x}" for x in ALLOWED_LIBRARIES}
               for row in rows):
            raise RuntimeError("control manifest contains invalid library")
        source_usage = Counter()
        donor_usage = Counter()
        recipient_usage = Counter()
        for row in rows:
            if row.get("schema_version") != \
                    "joint_doublet_control_manifest_v4" or row.get(
                        "scientific_method_version") != \
                    SCIENTIFIC_METHOD_VERSION or row.get(
                        "calibration_library") != str(CALIBRATION_LIBRARY):
                raise RuntimeError(
                    "control manifest method/calibration provenance mismatch")
            if row.get("control_class") not in {
                    "PLANTED_NEW_SOURCE", "SAME_SOURCE_NULL",
                    "ALREADY_PRESENT_DONOR_DOSAGE_SHIFT"}:
                raise RuntimeError("control manifest contains invalid class")
            if row.get("requested_fractions") != ",".join(
                    f"{value:.2f}" for value in CONTROL_FRACTIONS):
                raise RuntimeError("control fraction contract mismatch")
            library = row["library"]
            recipient = row.get("recipient_barcode", "")
            source = row.get("source_barcode", "")
            if not re.fullmatch(r"[ACGT]{16}", recipient) or not re.fullmatch(
                    r"[ACGT]{16}", source) or recipient == source:
                raise RuntimeError("control parent separation is invalid")
            if (library, recipient) in comparison_keys or \
                    (library, source) in comparison_keys or \
                    f"{library}:{recipient}" in frozen_ids or \
                    f"{library}:{source}" in frozen_ids:
                raise RuntimeError("control parents overlap targets/comparisons")
            selected_by_library[library].update((recipient, source))
            source_usage[(library, source)] += 1
            donor_usage[(library, row.get("source_identity", ""))] += 1
            recipient_usage[(library, recipient)] += 1
            for modality in ("rna", "atac"):
                if row.get(f"{modality}_molecule_evidence_basis", "") != \
                        row.get(f"source_{modality}_molecule_evidence_basis", ""):
                    raise RuntimeError("control source/recipient molecule basis mismatch")
            if control_path.name == "control_pairs_final.tsv" and \
                    row.get("control_class") == "PLANTED_NEW_SOURCE":
                eligible = truthy(row.get("decoy_comparison_eligible", ""))
                legal_states = {candidate.get("second_state", "")
                                for candidate_rows in candidate_rows_by_path.values()
                                for candidate in candidate_rows
                                if candidate.get("library") == library and
                                candidate.get("barcode") == recipient}
                if eligible and row.get(
                        "genotype_distance_matched_decoy", "") not in legal_states:
                    raise RuntimeError("final control decoy is not legal")
                if eligible:
                    for modality in ("rna", "atac"):
                        if int(row.get(
                                f"{modality}_shared_callable_sites", 0) or 0) < \
                                DECOY_MIN_SHARED_SITES or finite(row.get(
                                f"{modality}_distance_relative_mismatch", ""),
                                math.inf) > DECOY_RELATIVE_TOLERANCE or finite(
                                row.get(
                                    f"{modality}_opportunity_relative_mismatch",
                                    ""), math.inf) > DECOY_RELATIVE_TOLERANCE:
                            raise RuntimeError("final control decoy violates tolerance")
                elif not row.get("decoy_exclusion_reason"):
                    raise RuntimeError("excluded control decoy lacks reason")
        if set(source_usage) & set(recipient_usage):
            raise RuntimeError("control source and recipient roles overlap across classes")
        if max(source_usage.values(), default=0) > CONTROL_SOURCE_CELL_REUSE_CAP or \
                max(donor_usage.values(), default=0) > \
                CONTROL_SOURCE_DONOR_REUSE_CAP or max(
                    recipient_usage.values(), default=0) > \
                CONTROL_RECIPIENT_REUSE_CAP:
            raise RuntimeError("control reuse cap exceeded")
    unavailable_path = root / "unavailable_cells.tsv"
    if root.name == "benchmark_unsubmitted":
        unavailable_path = root.parent / "unavailable_cells.tsv"
    unavailable = list(read_tsv(unavailable_path)) if unavailable_path.is_file() else []
    if root.name != "benchmark_unsubmitted" and unavailable != blueprint.get("unavailable_cells", []):
        raise RuntimeError("unavailable-cell roster differs from declared workload")
    for row in unavailable:
        key = (row.get("library"), row.get("barcode"))
        if key not in comparison_keys and f"{key[0]}:{key[1]}" not in frozen_ids:
            raise RuntimeError("unavailable roster contains a foreign cell")
        selected_by_library[key[0]].discard(key[1])
    # Candidate manifests must cover the exact selected cell set for their
    # library; this proves the menu used by cache extraction is neither stale
    # nor a subset of the requested panel.
    for candidate_path, rows in candidate_rows_by_path.items():
        libraries_in_path = {row["library"] for row in rows}
        if len(libraries_in_path) != 1:
            raise RuntimeError("candidate manifest mixes libraries")
        library = next(iter(libraries_in_path))
        cells = {row["barcode"] for row in rows}
        if cells != selected_by_library[library]:
            raise RuntimeError(
                f"candidate manifest selected-cell set mismatch for {library}")
    ledger = _read_submission_ledger(root)
    current_ledger = [row for row in ledger
                      if row.get("workload_generation_id") == generation]
    valid_indices = set(indices)
    for record in current_ledger:
        if record.get("node") not in {
                "extraction", "control_finalizer", "analysis", "gather"} or \
                not str(record.get("job_id", "")).isdigit():
            raise RuntimeError("submission ledger record is malformed")
        if record.get("node") == task_kind and any(
                int(value) not in valid_indices
                for value in record.get("task_indices", [])):
            raise RuntimeError("submission ledger references foreign task index")
    return tasks[task_index], generation


def _task_output_paths(task, task_kind):
    if task_kind == "extraction":
        return [Path(str(task["cache_prefix"]) + suffix) for suffix in (
            ".metadata.json", ".cells.tsv", ".observations.bin", ".molecules.bin",
            ".sites.bin", ".genotypes.tsv.gz", ".samples.tsv")]
    return [Path(task["analysis_output"])]


def _file_states(paths):
    return {str(path): {"bytes": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns}
            for path in paths}


def _task_dependency_states(task, task_kind):
    if task_kind != "analysis":
        return {}
    paths = [Path(task["cache_prefix"] + ".metadata.json"), Path(task["cache_prefix"] + ".cells.tsv")]
    if task.get("action") == "CONTROL_BIN":
        paths.append(Path(task["control_manifest"]))
    return _file_states(paths)


def _validated_marker_unchanged(task, task_kind):
    try:
        marker = json.loads(Path(task["marker"]).read_text())
        valid = (marker.get("schema_version") == TASK_MARKER_SCHEMA and
                 marker.get("operational_status") == "COMPLETE" and
                 marker.get("implementation_version") == IMPLEMENTATION_VERSION and
                 marker.get("cache_schema") == CACHE_SCHEMA and
                 marker.get("task_kind") == task_kind and
                 marker.get("task_configuration") == dict(task) and
                 marker.get("audited_output_files") == _file_states(_task_output_paths(task, task_kind)) and
                 marker.get("audited_dependencies") == _task_dependency_states(task, task_kind))
        return valid, [] if valid else ["TERMINAL_MARKER_OR_AUDITED_OUTPUT_CHANGED"], marker
    except (OSError, ValueError, TypeError) as error:
        return False, [f"TERMINAL_MARKER_UNAVAILABLE:{error}"], {}


def validate_task_artifact(task, task_kind, require_marker=True, deep=True):
    if require_marker and not deep:
        return _validated_marker_unchanged(task, task_kind)
    validate_targeted_manifest_boundary([task], task_kind)
    problems = []
    generation = clean(task.get("workload_generation_id", ""))
    cache_generation = clean(task.get("cache_generation_id", generation))
    marker_path = Path(task["marker"])
    details = {
        "task_kind": task_kind, "task_index": int(task["task_index"]),
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "workload_generation_id": cache_generation
            if task_kind == "extraction" else generation,
        "library": task["library"], "modality": task["modality"],
        "action": task.get("action", "NORMALIZED_EXTRACTION")
            if task_kind == "analysis" else "NORMALIZED_EXTRACTION",
        "implementation_version": IMPLEMENTATION_VERSION,
        "task_configuration": dict(task),
        "cache_schema": CACHE_SCHEMA,
    }
    if task_kind == "extraction":
        prefix = Path(task["cache_prefix"])
        metadata_path = Path(str(prefix) + ".metadata.json")
        expected = [
            metadata_path, Path(str(prefix) + ".cells.tsv"),
            Path(str(prefix) + ".observations.bin"),
            Path(str(prefix) + ".molecules.bin"),
            Path(str(prefix) + ".sites.bin"),
            Path(str(prefix) + ".genotypes.tsv.gz"),
            Path(str(prefix) + ".samples.tsv"),
        ]
        for path in expected:
            if not path.is_file() or (path.stat().st_size == 0 and path.suffix != ".bin"):
                problems.append(f"MISSING_OR_EMPTY:{path}")
        metadata = {}
        if metadata_path.is_file():
            try:
                metadata = json.loads(metadata_path.read_text())
            except (OSError, json.JSONDecodeError) as error:
                problems.append(f"MALFORMED_METADATA:{error}")
        if metadata:
            checks = {
                "schema_version": CACHE_SCHEMA,
                "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
                "calibration_library": CALIBRATION_LIBRARY,
                "workload_generation_id": cache_generation,
                "library": task["library"], "modality": task["modality"],
                "manifest_digest": task["manifest_digest"],
                "error_ref": float(task["error_ref"]),
                "error_alt": float(task["error_alt"]),
                "min_evidence": int(task["min_evidence"]),
                "max_second_fraction": float(task["max_second_fraction"]),
                "scheduler_memory_bytes": int(task["scheduler_memory_bytes"]),
                "launcher_runtime_reserve_bytes": int(
                    task["launcher_runtime_reserve_bytes"]),
                "worker_memory_budget_bytes": int(
                    task["worker_memory_budget_bytes"]),
            }
            for field, expected_value in checks.items():
                actual = metadata.get(field, "")
                equal = math.isclose(float(actual), float(expected_value),
                                     rel_tol=0, abs_tol=1e-15) \
                    if field in {"error_ref", "error_alt",
                                 "max_second_fraction"} else \
                    str(actual) == str(expected_value)
                if not equal:
                    problems.append(
                        f"METADATA_MISMATCH:{field}:{metadata.get(field)}!={expected_value}")
            sizes = (
                ("observation_records", "observation_record_bytes",
                 Path(str(prefix) + ".observations.bin")),
                ("molecule_records", "molecule_record_bytes",
                 Path(str(prefix) + ".molecules.bin")),
            )
            for count_field, size_field, path in sizes:
                if path.is_file() and metadata.get(count_field) is not None:
                    expected_size = int(metadata[count_field]) * int(metadata[size_field])
                    if path.stat().st_size != expected_size:
                        problems.append(
                            f"BINARY_SIZE_MISMATCH:{path}:{path.stat().st_size}!={expected_size}")
            for source_field in ("samples_path", "sites_path",
                                 "observations_path", "molecules_path"):
                source = metadata.get(source_field, "")
                if not source or not Path(source).is_absolute():
                    problems.append(f"METADATA_SOURCE_NOT_ABSOLUTE:{source_field}")
            source_contract = {
                "samples_path": "samples", "sites_path": "pileup_sites",
                "observations_path": "pileup_observations",
                "molecules_path": "pileup_molecules",
            }
            for metadata_field, task_field in source_contract.items():
                if os.path.abspath(str(metadata.get(metadata_field, ""))) != \
                        os.path.abspath(str(task.get(task_field, ""))):
                    problems.append(
                        f"METADATA_SOURCE_PATH_MISMATCH:{metadata_field}")
            if metadata.get("source_scans") != {
                    "sites": 1, "observations": 1, "molecules": 1}:
                problems.append("METADATA_SOURCE_SCAN_COUNT_INVALID")
            if metadata.get("atomic_publication") != "metadata published last":
                problems.append("CACHE_ATOMIC_PUBLICATION_UNPROVEN")
            for field in ("samples_source_bytes", "sites_source_bytes",
                          "observations_source_bytes", "molecules_source_bytes"):
                if int(metadata.get(field, 0) or 0) <= 0:
                    problems.append(f"METADATA_SOURCE_SIZE_INVALID:{field}")
            for field in ("samples_content_digest", "sites_content_digest",
                          "observations_content_digest", "molecules_content_digest"):
                if not re.fullmatch(r"fnv1a64:[0-9a-f]{16}",
                                    str(metadata.get(field, ""))):
                    problems.append(f"METADATA_DIGEST_INVALID:{field}")
            projected = metadata.get("bucket_projected_peak_bytes", [])
            measured = metadata.get("bucket_measured_bytes", [])
            limit = int(metadata.get("bucket_memory_limit_bytes", 0) or 0)
            if not projected or len(projected) != len(measured) or limit <= 0:
                problems.append("BUCKET_MEMORY_ACCOUNTING_INVALID")
            elif any(int(value) > limit for value in projected):
                problems.append("BUCKET_PROJECTED_PEAK_EXCEEDS_LIMIT")
            if int(task.get("projected_peak_bytes", 0) or 0) > int(
                    task.get("bounded_memory_limit_bytes", 0) or 0):
                problems.append("WORKLOAD_PROJECTED_PEAK_EXCEEDS_DECLARED_LIMIT")
            component_sizes = {
                "observations_binary_bytes": Path(str(prefix) + ".observations.bin"),
                "molecules_binary_bytes": Path(str(prefix) + ".molecules.bin"),
                "sites_binary_bytes": Path(str(prefix) + ".sites.bin"),
                "cells_index_bytes": Path(str(prefix) + ".cells.tsv"),
                "samples_dictionary_bytes": Path(str(prefix) + ".samples.tsv"),
                "genotypes_dictionary_bytes": Path(
                    str(prefix) + ".genotypes.tsv.gz"),
            }
            for field, path in component_sizes.items():
                if not path.is_file() or int(metadata.get(field, -1)) != \
                        path.stat().st_size:
                    problems.append(f"CACHE_COMPONENT_SIZE_MISMATCH:{field}")
                digest_field = field.removesuffix("_bytes") + \
                    "_content_digest"
                if path.is_file() and metadata.get(digest_field, "") != \
                        _file_fnv1a64(path):
                    problems.append(
                        f"CACHE_COMPONENT_DIGEST_MISMATCH:{digest_field}")
            if metadata.get("cache_component_digest_definition") != (
                    "LEGACY_PROJECT_FNV1A64 over exact published component bytes"):
                problems.append("CACHE_COMPONENT_DIGEST_DEFINITION_INVALID")
        cells_path = Path(str(prefix) + ".cells.tsv")
        if cells_path.is_file():
            try:
                cell_rows = list(read_tsv(cells_path))
                required_cell_fields = {
                    "schema_version", "workload_generation_id", "library",
                    "modality", "barcode", "encoded_barcode",
                    "observation_offset", "observation_count",
                    "molecule_offset", "molecule_count",
                }
                if cell_rows and not required_cell_fields.issubset(cell_rows[0]):
                    problems.append("CELL_INDEX_SCHEMA_MISSING_FIELDS")
                keys = [(row.get("library"), row.get("modality"),
                         row.get("barcode")) for row in cell_rows]
                if len(keys) != len(set(keys)):
                    problems.append("CELL_INDEX_DUPLICATE_KEYS")
                if metadata and len(cell_rows) != int(metadata.get(
                        "selected_cells", -1)):
                    problems.append(
                        f"CELL_INDEX_COUNT:{len(cell_rows)}!={metadata.get('selected_cells')}")
                candidate_cells = {row.get("barcode", "") for row in
                                   read_tsv(task["candidate_manifest"])}
                indexed_cells = {row.get("barcode", "") for row in cell_rows}
                if candidate_cells != indexed_cells:
                    problems.append("CELL_INDEX_SELECTED_MANIFEST_CELLS_MISMATCH")
                if int(task.get("selected_cells", -1) or -1) != len(indexed_cells):
                    problems.append("CELL_INDEX_TASK_SELECTED_COUNT_MISMATCH")
                for row in cell_rows:
                    if row.get("schema_version") != CACHE_SCHEMA or \
                            row.get("workload_generation_id") != cache_generation or \
                            row.get("library") != task["library"] or \
                            row.get("modality") != task["modality"]:
                        problems.append("CELL_INDEX_PROVENANCE_MISMATCH")
                        break
                for channel, total_field in (("observation", "observation_records"),
                                             ("molecule", "molecule_records")):
                    intervals = []
                    try:
                        for row in cell_rows:
                            offset = int(row[f"{channel}_offset"])
                            count = int(row[f"{channel}_count"])
                            if offset < 0 or count < 0:
                                raise ValueError("negative offset/count")
                            intervals.append((offset, offset + count))
                        nonempty = sorted(item for item in intervals
                                          if item[1] > item[0])
                        if any(left[1] > right[0]
                               for left, right in zip(nonempty, nonempty[1:])):
                            problems.append(
                                f"CELL_INDEX_OVERLAPPING_{channel.upper()}_SLICES")
                        total = int(metadata.get(total_field, -1)) \
                            if metadata else -1
                        if total >= 0 and (sum(end - start for start, end in intervals)
                                           != total or
                                           any(end > total for _, end in intervals)):
                            problems.append(
                                f"CELL_INDEX_{channel.upper()}_COUNT_MISMATCH")
                    except (KeyError, TypeError, ValueError):
                        problems.append(
                            f"CELL_INDEX_{channel.upper()}_OFFSET_INVALID")
                details["selected_cells"] = len(cell_rows)
            except (OSError, UnicodeDecodeError, csv.Error) as error:
                problems.append(f"CELL_INDEX_INVALID:{error}")
        genotype_path = Path(str(prefix) + ".genotypes.tsv.gz")
        if genotype_path.is_file():
            try:
                header, rows, digest = _full_gzip_audit(genotype_path)
                required = {"schema_version", "library", "modality", "tid", "pos",
                            "sample", "genotype"}
                if not required.issubset(header):
                    problems.append("GENOTYPE_SCHEMA_MISSING_FIELDS")
                details.update({"genotype_rows": rows,
                                "genotype_content_sha256": digest})
            except (OSError, EOFError, gzip.BadGzipFile, UnicodeDecodeError) as error:
                problems.append(f"GENOTYPE_GZIP_INVALID:{error}")
        samples_path = Path(str(prefix) + ".samples.tsv")
        if samples_path.is_file():
            try:
                sample_rows = list(read_tsv(samples_path))
                required_samples = {"cache_donor_slot", "sample_index", "sample"}
                if sample_rows and not required_samples.issubset(sample_rows[0]):
                    problems.append("SAMPLE_DICTIONARY_SCHEMA_MISSING_FIELDS")
                for field in required_samples:
                    values = [row.get(field, "") for row in sample_rows]
                    if any(not value for value in values) or \
                            len(values) != len(set(values)):
                        problems.append(
                            f"SAMPLE_DICTIONARY_KEYS_INVALID:{field}")
                if metadata and len(sample_rows) != int(metadata.get("donors", -1)):
                    problems.append("SAMPLE_DICTIONARY_DONOR_COUNT_MISMATCH")
            except (OSError, UnicodeDecodeError, csv.Error, ValueError) as error:
                problems.append(f"SAMPLE_DICTIONARY_INVALID:{error}")
        details.update({
            "expected_outputs": [str(path) for path in expected],
            "cache_generation_id": cache_generation,
        })
    elif task_kind == "analysis":
        output = Path(task["analysis_output"])
        if not output.is_file() or output.stat().st_size == 0:
            problems.append(f"MISSING_OR_EMPTY:{output}")
        else:
            try:
                header, rows, digest = _full_gzip_audit(output)
                required = {"schema_version", "workload_generation_id",
                            "scientific_method_version", "calibration_library",
                            "task_index", "result_key", "result_class",
                            "action", "library", "modality"}
                if not required.issubset(header):
                    problems.append("ANALYSIS_SCHEMA_MISSING_FIELDS")
                if rows == 0:
                    problems.append("ANALYSIS_ZERO_ROWS")
                result_rows = list(read_tsv(output))
                keys = [row.get("result_key", "") for row in result_rows]
                if any(not key for key in keys) or len(keys) != len(set(keys)):
                    problems.append("ANALYSIS_RESULT_KEYS_INVALID")
                for row in result_rows:
                    expected_values = {
                        "schema_version": "joint_doublet_benchmark_equivalence_v4" if row.get("result_class") == "THREAD_AND_REFERENCE_EQUIVALENCE" else "joint_doublet_compiled_targeted_analysis_v5",
                        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
                        "calibration_library": str(CALIBRATION_LIBRARY),
                        "workload_generation_id": generation,
                        "task_index": str(task["task_index"]),
                        "library": task["library"],
                        "modality": task["modality"],
                        "action": task["action"],
                    }
                    for field, expected_value in expected_values.items():
                        if field in row and row.get(field, "") != expected_value:
                            problems.append(
                                f"ANALYSIS_ROW_MISMATCH:{field}:{row.get(field)}!={expected_value}")
                            break
                    for field in ("error_ref", "error_alt", "min_evidence",
                                  "max_second_fraction"):
                        if field in row and str(row.get(field, "")) != str(
                                task.get(field, "")):
                            try:
                                if not math.isclose(float(row.get(field, "")),
                                                    float(task.get(field, "")),
                                                    rel_tol=0, abs_tol=1e-15):
                                    problems.append(
                                        f"ANALYSIS_MODEL_PARAMETER_MISMATCH:{field}")
                            except ValueError:
                                problems.append(
                                    f"ANALYSIS_MODEL_PARAMETER_INVALID:{field}")
                classes = Counter(row.get("result_class", "")
                                  for row in result_rows)
                benchmark_equivalence = set(classes) == {
                    "THREAD_AND_REFERENCE_EQUIVALENCE"}
                if not benchmark_equivalence:
                    if classes.get("COMPILED_SHARD_ACCOUNTING", 0) != 1:
                        problems.append("ANALYSIS_ACCOUNTING_ROW_COUNT_INVALID")
                    if task["action"] == "CELL":
                        if classes.get("FULL_MENU_DOWNSAMPLING", 0) != \
                                2 * len(DOWNSAMPLE_FRACTIONS):
                            problems.append("ANALYSIS_DOWNSAMPLE_ROW_COUNT_INVALID")
                        if classes.get("FULL_MENU_DOWNSAMPLING_REPLICATE", 0) != \
                                2 * len(DOWNSAMPLE_FRACTIONS) * \
                                DOWNSAMPLE_REPLICATES:
                            problems.append(
                                "ANALYSIS_DOWNSAMPLE_REPLICATE_COUNT_INVALID")
                        if classes.get("CELL_CONDITIONAL_LOCKED_MODEL_NULL", 0) != 2:
                            problems.append("ANALYSIS_NULL_ROW_COUNT_INVALID")
                        if classes.get("CELL_CONDITIONAL_NULL_REPLICATE", 0) != \
                                2 * CELL_NULL_REPLICATES:
                            problems.append("ANALYSIS_NULL_REPLICATE_COUNT_INVALID")
                        if classes.get("OBSERVED_FULL_MENU_CANDIDATE", 0) == 0:
                            problems.append("ANALYSIS_OBSERVED_ROWS_MISSING")
                        observed = [row for row in result_rows
                                    if row.get("result_class") ==
                                    "OBSERVED_FULL_MENU_CANDIDATE" and
                                    row.get("engine") == "OPTIMIZED_BATCHED"]
                        menu_sizes = {int(row.get("legal_menu_candidates", 0) or 0)
                                      for row in observed}
                        if len(menu_sizes) != 1 or not menu_sizes or \
                                len(observed) != 2 * next(iter(menu_sizes), 0):
                            problems.append(
                                "ANALYSIS_OBSERVED_FULL_MENU_COUNT_INVALID")
                        expected_candidates = {row["candidate_id"] for row in read_tsv(task["candidate_manifest"])
                                               if row["barcode"] == task["barcode"] and row["library"] == task["library"]}
                        if {(row.get("evidence_channel"),row.get("candidate_id")) for row in observed} != {
                                (channel,candidate) for channel in ("SITE","MOLECULE") for candidate in expected_candidates}:
                            problems.append("ANALYSIS_FULL_MENU_IDENTITIES_MISMATCH")
                        for result_class, fractions, count in (("FULL_MENU_DOWNSAMPLING_REPLICATE",DOWNSAMPLE_FRACTIONS,DOWNSAMPLE_REPLICATES),
                                                                 ("CELL_CONDITIONAL_NULL_REPLICATE",(None,),CELL_NULL_REPLICATES)):
                            actual = {(row.get("evidence_channel"),float(row["fraction"]) if fractions != (None,) else None,
                                       int(row.get("replicate_index",-1))) for row in result_rows if row.get("result_class") == result_class}
                            expected = {(channel,fraction,replicate) for channel in ("SITE","MOLECULE") for fraction in fractions for replicate in range(count)}
                            if actual != expected:
                                problems.append("ANALYSIS_REPLICATE_IDENTITIES_MISMATCH:"+result_class)
                        if any(row.get("barcode") != task["barcode"] for row in result_rows if row.get("result_class") != "COMPILED_SHARD_ACCOUNTING"):
                            problems.append("ANALYSIS_CELL_IDENTITY_MISMATCH")
                        downsample_summaries = [row for row in result_rows
                                                if row.get("result_class") ==
                                                "FULL_MENU_DOWNSAMPLING"]
                        combinations = {(row.get("evidence_channel"),
                                         float(row.get("fraction", 0)))
                                        for row in downsample_summaries}
                        expected_combinations = {(channel, fraction)
                            for channel in ("SITE", "MOLECULE")
                            for fraction in DOWNSAMPLE_FRACTIONS}
                        if combinations != expected_combinations:
                            problems.append("ANALYSIS_DOWNSAMPLE_CONTRACT_INCOMPLETE")
                        for summary in downsample_summaries:
                            requested = int(summary.get(
                                "requested_replicates", 0) or 0)
                            attempted = int(summary.get(
                                "attempted_replicates", 0) or 0)
                            successful = int(summary.get(
                                "successful_replicates", 0) or 0)
                            unavailable = int(summary.get(
                                "scientifically_unavailable_replicates", 0) or 0)
                            failures = int(summary.get(
                                "technical_failures", 0) or 0)
                            if requested != DOWNSAMPLE_REPLICATES or \
                                    requested != attempted + failures or \
                                    attempted != successful + unavailable or \
                                    failures != 0:
                                problems.append("DOWNSAMPLE_SUCCESS_ACCOUNTING_INVALID")
                            expected_eligible = successful == \
                                DOWNSAMPLE_REPLICATES
                            if truthy(summary.get("decision_eligible", "")) != \
                                    expected_eligible or summary.get(
                                        "calibrated_support_status") != \
                                    "CALIBRATED_SUPPORT_NOT_EVALUATED" or \
                                    (successful == 0 and summary.get(
                                        "status") == "AVAILABLE"):
                                problems.append(
                                    "DOWNSAMPLE_DECISION_OR_CALIBRATION_INVALID")
                            try:
                                model_counts = json.loads(summary.get(
                                    "model_counts", "{}") or "{}")
                                if "" in model_counts or sum(
                                        int(value) for value in
                                        model_counts.values()) != successful:
                                    problems.append(
                                        "DOWNSAMPLE_MODEL_ACCOUNTING_INVALID")
                            except (ValueError, TypeError, json.JSONDecodeError):
                                problems.append("DOWNSAMPLE_MODEL_COUNTS_INVALID")
                        null_summaries = [row for row in result_rows
                                          if row.get("result_class") ==
                                          "CELL_CONDITIONAL_LOCKED_MODEL_NULL"]
                        if {row.get("evidence_channel") for row in null_summaries} != {"SITE", "MOLECULE"}:
                            problems.append("ANALYSIS_NULL_CHANNELS_INCOMPLETE")
                        for summary in null_summaries:
                            requested = int(summary.get(
                                "requested_replicates", 0) or 0)
                            attempted = int(summary.get(
                                "attempted_replicates", 0) or 0)
                            successful = int(summary.get(
                                "successful_replicates", 0) or 0)
                            unavailable = int(summary.get(
                                "scientifically_unavailable_replicates", 0) or 0)
                            failures = int(summary.get(
                                "technical_failures", 0) or 0)
                            if requested != CELL_NULL_REPLICATES or \
                                    requested != attempted + failures or \
                                    attempted != successful + unavailable or \
                                    failures != 0:
                                problems.append("NULL_SUCCESS_ACCOUNTING_INVALID")
                            expected_eligible = successful == \
                                CELL_NULL_REPLICATES
                            if truthy(summary.get("decision_eligible", "")) != \
                                    expected_eligible or summary.get(
                                        "evidence_category") != \
                                    "SUPPORT_CATEGORY_NOT_EVALUATED" or \
                                    summary.get("calibrated_support_status") != \
                                    "CALIBRATED_SUPPORT_NOT_EVALUATED" or \
                                    (successful == 0 and any(clean(summary.get(
                                        field, "")) not in {"", "NA"} for field in (
                                            "empirical_upper_tail_probability",
                                            "empirical_p_numerator",
                                            "empirical_p_denominator"))):
                                problems.append(
                                    "NULL_DECISION_OR_CALIBRATION_INVALID")
                        for channel in ("SITE", "MOLECULE"):
                            channel_rows = [row for row in observed
                                            if row.get("evidence_channel") == channel]
                            candidates = [row.get("candidate_id", "")
                                          for row in channel_rows]
                            if any(not value for value in candidates) or \
                                    len(candidates) != len(set(candidates)):
                                problems.append(
                                    f"ANALYSIS_{channel}_CANDIDATE_KEYS_INVALID")
                    else:
                        expected_controls = len([
                            value for value in task.get("control_ids", "").split(",")
                            if value])
                        summaries = [row for row in result_rows
                                     if row.get("result_class") ==
                                     "SOURCE_DISJOINT_CONTROL_SUMMARY" and
                                     row.get("engine") == "OPTIMIZED_BATCHED"]
                        if len(summaries) != expected_controls * 2 * \
                                len(CONTROL_FRACTIONS):
                            problems.append(
                                "ANALYSIS_CONTROL_SUMMARY_COUNT_INVALID")
                        control_ids = set(filter(None,task.get("control_ids", "").split(",")))
                        expected_keys = {(cid,channel,fraction) for cid in control_ids for channel in ("SITE","MOLECULE") for fraction in CONTROL_FRACTIONS}
                        if {(row.get("control_id"),row.get("evidence_channel"),float(row.get("fraction",0))) for row in summaries} != expected_keys:
                            problems.append("ANALYSIS_CONTROL_IDENTITIES_MISMATCH")
                        recipients = {row["control_id"]:row["recipient_barcode"] for row in read_tsv(task["control_manifest"]) if row["control_id"] in control_ids}
                        candidate_by_barcode = defaultdict(set)
                        for row in read_tsv(task["candidate_manifest"]):
                            candidate_by_barcode[row["barcode"]].add(row["candidate_id"])
                        expected_candidate_keys = {(cid,channel,fraction,candidate) for cid,channel,fraction in expected_keys
                                                   for candidate in candidate_by_barcode[recipients[cid]]}
                        actual_candidate_keys = {(row.get("control_id"),row.get("evidence_channel"),float(row.get("fraction",0)),row.get("candidate_id"))
                                                 for row in result_rows if row.get("result_class") == "SOURCE_DISJOINT_CONTROL_CANDIDATE" and row.get("engine") == "OPTIMIZED_BATCHED"}
                        if actual_candidate_keys != expected_candidate_keys:
                            problems.append("ANALYSIS_CONTROL_COMPLETE_MENU_MISMATCH")
                        per_control = Counter(row.get("control_id", "")
                                              for row in summaries)
                        if set(per_control.values()) != ({6} if per_control else set()):
                            problems.append(
                                "ANALYSIS_CONTROL_FRACTION_CHANNEL_COUNT_INVALID")
                        if any(row.get("calibrated_support_status") !=
                               "CALIBRATED_SUPPORT_NOT_EVALUATED"
                               for row in summaries):
                            problems.append(
                                "ANALYSIS_CONTROL_CALIBRATED_SUPPORT_FORBIDDEN")
                elif not Path(str(output) +
                              ".field_differences.tsv").is_file():
                    problems.append("EQUIVALENCE_FIELD_DIFFERENCE_TABLE_MISSING")
                details.update({"output": str(output), "rows": rows,
                                "output_content_sha256": digest,
                                "result_classes": dict(classes)})
            except (OSError, EOFError, gzip.BadGzipFile, UnicodeDecodeError) as error:
                problems.append(f"ANALYSIS_GZIP_INVALID:{error}")
    else:
        raise ValueError(f"unsupported task kind: {task_kind}")

    marker = None
    if require_marker and marker_path.is_file():
        try:
            marker = json.loads(marker_path.read_text())
        except (OSError, json.JSONDecodeError) as error:
            problems.append(f"MALFORMED_MARKER:{error}")
        if marker:
            expected_marker = {
                "schema_version": TASK_MARKER_SCHEMA,
                "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
                "calibration_library": CALIBRATION_LIBRARY,
                "operational_status": "COMPLETE",
                "task_kind": task_kind,
                "workload_generation_id": cache_generation
                    if task_kind == "extraction" else generation,
                "task_index": int(task["task_index"]),
                "library": task["library"], "modality": task["modality"],
                "action": task.get("action", "NORMALIZED_EXTRACTION")
                    if task_kind == "analysis" else "NORMALIZED_EXTRACTION",
            }
            for field, expected_value in expected_marker.items():
                if marker.get(field) != expected_value:
                    problems.append(
                        f"MARKER_MISMATCH:{field}:{marker.get(field)}!={expected_value}")
            if marker.get("audited_dependencies") != _task_dependency_states(task, task_kind):
                problems.append("MARKER_UPSTREAM_DEPENDENCY_CHANGED")
            if marker.get("implementation_version") != IMPLEMENTATION_VERSION or marker.get("task_configuration") != dict(task):
                problems.append("MARKER_TASK_CONFIGURATION_MISMATCH")
            if task_kind == "analysis" and details.get(
                    "output_content_sha256") and marker.get(
                        "output_content_sha256") != details[
                            "output_content_sha256"]:
                problems.append("MARKER_OUTPUT_DIGEST_MISMATCH")
            if task_kind == "extraction":
                for field in ("cache_generation_id", "selected_cells",
                              "genotype_content_sha256"):
                    if field in details and marker.get(field) != details[field]:
                        problems.append(
                            f"MARKER_OUTPUT_DETAIL_MISMATCH:{field}")
            field_difference = Path(str(task.get("analysis_output", "")) +
                                    ".field_differences.tsv")
            if task_kind == "analysis" and field_difference.is_file() and \
                    marker.get("field_difference_content_sha256") != \
                    _small_file_sha256(field_difference):
                problems.append("MARKER_FIELD_DIFFERENCE_DIGEST_MISMATCH")
            if task_kind == "analysis" and marker.get("single_thread_output"):
                for path_field, digest_field in (
                        ("single_thread_output", "single_thread_content_sha256"),
                        ("production_thread_output",
                         "production_thread_content_sha256")):
                    benchmark_output = Path(marker.get(path_field, ""))
                    if not benchmark_output.is_file() or \
                            benchmark_output.stat().st_size == 0:
                        problems.append(
                            f"BENCHMARK_RAW_OUTPUT_MISSING:{path_field}")
                        continue
                    try:
                        _header, benchmark_rows, _digest = _full_gzip_audit(
                            benchmark_output)
                        if benchmark_rows == 0:
                            problems.append(
                                f"BENCHMARK_RAW_OUTPUT_EMPTY:{path_field}")
                    except (OSError, EOFError, gzip.BadGzipFile,
                            UnicodeDecodeError) as error:
                        problems.append(
                            f"BENCHMARK_RAW_OUTPUT_INVALID:{path_field}:{error}")
                    if marker.get(digest_field) != _small_file_sha256(
                            benchmark_output):
                        problems.append(
                            f"BENCHMARK_RAW_OUTPUT_DIGEST_MISMATCH:{path_field}")
                if marker.get("thread_material_differences") != 0 or \
                        marker.get("optimized_reference_failures"):
                    problems.append("BENCHMARK_EQUIVALENCE_NOT_CLEAN")
                if marker.get("legacy_completed_score_failures"):
                    problems.append("BENCHMARK_LEGACY_SCORE_EQUIVALENCE_NOT_CLEAN")
                if task.get("action") == "CELL" and int(marker.get(
                        "legacy_completed_score_fields_compared", 0) or 0) <= 0:
                    problems.append("BENCHMARK_LEGACY_SCORE_COMPARISON_MISSING")
                legacy_reference = Path(marker.get("legacy_reference", ""))
                if not legacy_reference.is_file() or \
                        marker.get("legacy_reference_content_sha256") != \
                        _small_file_sha256(legacy_reference):
                    problems.append("BENCHMARK_LEGACY_REFERENCE_DIGEST_MISMATCH")
    elif require_marker:
        problems.append("MISSING_MARKER")
    if not problems:
        details["audited_output_files"] = _file_states(_task_output_paths(task, task_kind))
        details["audited_dependencies"] = _task_dependency_states(task, task_kind)
    return not problems, problems, details


def validate_task_action(args):
    index = int(args.task_index)
    task, generation = validate_targeted_input_contract(
        args.task_manifest, args.task_kind, index)
    if args.input_only:
        payload = {"valid": True, "input_contract": "COMPLETE",
                   "workload_generation_id": generation,
                   "task_index": index, "task_kind": args.task_kind}
        if not args.quiet:
            print(json.dumps(payload, indent=2))
        return 0
    valid, problems, details = validate_task_artifact(
        task, args.task_kind, require_marker=not args.write_marker,
        deep=args.write_marker or not args.quiet)
    if args.write_marker:
        if not valid and problems == ["MISSING_MARKER"]:
            valid = True
            problems = []
        if valid:
            atomic_json(task["marker"], {
                "schema_version": TASK_MARKER_SCHEMA,
                "operational_status": "COMPLETE", "utc": utc_now(),
                **details,
            })
            valid, problems, details = validate_task_artifact(
                task, args.task_kind, require_marker=True, deep=False)
    payload = {"valid": valid, "problems": problems, **details}
    if not args.quiet:
        print(json.dumps(payload, indent=2))
    return 0 if valid else 2


def _parse_composition(text):
    result = {}
    for token in clean(text).split(","):
        if not token:
            continue
        donor, copies = token.rsplit(":", 1)
        result[donor] = result.get(donor, 0.0) + float(copies)
    return result


def _expected_alt_fraction(composition, genotypes):
    if not composition:
        return math.nan
    numerator = denominator = 0.0
    for donor, copies in composition.items():
        genotype = genotypes.get(donor)
        if genotype is None or genotype < 0:
            return math.nan
        numerator += copies * genotype / 2.0
        denominator += copies
    return numerator / denominator if denominator else math.nan


def _control_distance_cache(tasks, controls, candidate_by_library):
    """One bounded numeric pass per genotype dictionary; reuse composition vectors."""
    requested = defaultdict(set)
    for control in controls:
        if control["control_class"] != "PLANTED_NEW_SOURCE":
            continue
        candidates = candidate_by_library[control["library"]][control["recipient_barcode"]]
        planted = next((row for row in candidates if row["second_state"] == control["expected_contributor"]), None)
        if not candidates or planted is None:
            continue
        for decoy in candidates:
            if decoy["second_state"] != control["expected_contributor"] and truthy(decoy.get("physical_pool_state")) and "CONTAINS_NEW_DONOR" in clean(decoy.get("structural_added_state_relationship")):
                requested[control["library"]].add((candidates[0]["locked_copy_vector"], planted["second_copy_vector"], decoy["second_copy_vector"]))
    output = {}
    for task in tasks:
        library, modality = task["library"], task["modality"]
        triples = sorted(requested[library])
        if not triples:
            continue
        samples = [row["sample"] for row in read_tsv(task["cache_prefix"] + ".samples.tsv")]
        sample_index = {value:index for index,value in enumerate(samples)}
        definitions = {key:_parse_composition(key) for triple in triples for key in triple}
        sums = {triple:[0, np.longdouble(0), np.longdouble(0), 0, 0, np.longdouble(0)] for triple in triples}
        block = []
        def reduce_block():
            if not block:
                return
            matrix = np.asarray(block, dtype=np.int8)
            probabilities = {}
            for key, composition in definitions.items():
                if not composition or not set(composition).issubset(sample_index):
                    probabilities[key] = np.full(len(matrix), np.nan)
                    continue
                columns = [sample_index[donor] for donor in composition]
                weights = np.asarray(list(composition.values()), dtype=float)
                values = matrix[:,columns]
                probability = (values@weights)/(2*weights.sum())
                probability[np.any(values<0,axis=1)] = np.nan
                probabilities[key] = probability
            for triple in triples:
                locked, planted, decoy = (probabilities[key] for key in triple)
                valid = np.isfinite(locked)&np.isfinite(planted)&np.isfinite(decoy)
                first = np.abs(locked[valid]-planted[valid]);second = np.abs(locked[valid]-decoy[valid])
                acc = sums[triple]
                acc[0] += int(valid.sum());acc[1] += first.sum(dtype=np.longdouble);acc[2] += second.sum(dtype=np.longdouble)
                acc[3] += int(np.count_nonzero(first));acc[4] += int(np.count_nonzero(second))
                acc[5] += np.abs(planted[valid]-decoy[valid]).sum(dtype=np.longdouble)
            block.clear()
        current = None
        genotypes = None
        for row in read_tsv(task["cache_prefix"] + ".genotypes.tsv.gz"):
            key = (row["tid"],row["pos"])
            if key != current:
                if genotypes is not None:
                    block.append(genotypes)
                    if len(block)>=4096:
                        reduce_block()
                current = key;genotypes = np.full(len(samples), -1, dtype=np.int8)
            if row["sample"] not in sample_index:
                raise RuntimeError("genotype donor absent from normalized sample dictionary")
            genotypes[sample_index[row["sample"]]] = int(row["genotype"])
        if genotypes is not None:
            block.append(genotypes)
        reduce_block()
        for triple, (shared,first,second,first_n,second_n,separation) in sums.items():
            first = float(first/shared) if shared else math.nan
            second = float(second/shared) if shared else math.nan
            output[(library,modality,*triple)] = {
                "shared_callable_sites": shared, "genotype_site_count": shared,
                "planted_distance": first, "decoy_distance": second,
                "expected_decoy_genotype_distance": float(separation/shared) if shared else math.nan,
                "distance_absolute_mismatch": abs(first-second),
                "distance_relative_mismatch": abs(first-second)/first if first>0 else math.inf,
                "planted_discriminating_sites": first_n, "decoy_discriminating_sites": second_n,
                "opportunity_absolute_mismatch": abs(first_n-second_n),
                "opportunity_relative_mismatch": abs(first_n-second_n)/first_n if first_n else math.inf,
            }
    return output


def finalize_controls(args):
    root = Path(args.targeted_root).resolve()
    extraction_path = root / "manifests" / "extraction_tasks.tsv"
    analysis_path = root / "manifests" / "analysis_tasks.tsv"
    tasks = list(read_tsv(extraction_path))
    analysis_tasks = list(read_tsv(analysis_path))
    if not tasks:
        raise RuntimeError("no normalized-cache extraction tasks")
    # This shared boundary reads and validates targets, comparisons, complete
    # candidate menus, provisional controls, generation metadata, and both
    # task manifests before the first cache/genotype artifact is inspected.
    validate_targeted_input_contract(extraction_path, "extraction", 0)
    if analysis_tasks:
        validate_targeted_input_contract(analysis_path, "analysis", 0)
    validate_targeted_manifest_boundary(tasks, "extraction")
    validate_targeted_manifest_boundary(analysis_tasks, "analysis")
    for task in tasks:
        valid, problems, _ = validate_task_artifact(
            task, "extraction", require_marker=True)
        if not valid:
            raise RuntimeError(
                f"extraction task {task['task_index']} invalid: {';'.join(problems)}")
    if any(task["workload_generation_id"] != args.generation for task in tasks):
        raise RuntimeError("control finalizer generation does not match workload tasks")
    provisional = root / "manifests" / "control_pairs_provisional.tsv"
    provisional_rows = list(read_tsv(provisional))
    if not provisional_rows:
        raise RuntimeError("provisional control manifest is empty")
    for row in provisional_rows:
        library = int(library_name(row.get("library", "")).removeprefix("lib"))
        if library not in ALLOWED_LIBRARIES or library in PROTECTED_LIBRARIES:
            raise RuntimeError(
                f"control finalizer rejected unexpected library {library}")
    for task in tasks:
        for row in read_tsv(task["candidate_manifest"]):
            if row.get("library") != task["library"] or \
                    int(library_name(row.get("library", "")).removeprefix(
                        "lib")) not in ALLOWED_LIBRARIES:
                raise RuntimeError(
                    "candidate manifest library boundary failed before genotype open")
    candidate_by_library = defaultdict(lambda: defaultdict(list))
    seen_libraries = set()
    for task in tasks:
        if task["library"] in seen_libraries:
            continue
        seen_libraries.add(task["library"])
        manifest = task["candidate_manifest"]
        for row in read_tsv(manifest):
            candidate_by_library[row["library"]][row["barcode"]].append(row)

    distance_cache = _control_distance_cache(tasks, provisional_rows, candidate_by_library)
    def distances(library, recipient, planted, decoy, modality):
        return distance_cache[(library,modality,recipient["locked_copy_vector"],
                               planted["second_copy_vector"],decoy["second_copy_vector"])]

    output_rows = []
    decoy_diagnostics = []
    for control in provisional_rows:
        result = dict(control)
        result["schema_version"] = "joint_doublet_control_manifest_v4"
        result["scientific_method_version"] = SCIENTIFIC_METHOD_VERSION
        result["calibration_library"] = CALIBRATION_LIBRARY
        result["workload_generation_id"] = args.generation
        if control["control_class"] != "PLANTED_NEW_SOURCE":
            result["decoy_status"] = "NOT_APPLICABLE"
            result["decoy_comparison_eligible"] = False
            result["decoy_exclusion_reason"] = "NOT_APPLICABLE"
            output_rows.append(result)
            continue
        candidates = candidate_by_library[control["library"]][
            control["recipient_barcode"]]
        recipient = candidates[0] if candidates else None
        planted = next((row for row in candidates
                        if row["second_state"] == control["expected_contributor"]), None)
        alternatives = [row for row in candidates
                        if row["second_state"] != control["expected_contributor"] and
                        row["second_state"] in clean(control.get("legal_evidence_contributors")).split(",") and
                        truthy(row.get("physical_pool_state", "")) and
                        "CONTAINS_NEW_DONOR" in clean(row.get(
                            "structural_added_state_relationship", ""))]
        eligible = []
        if recipient and planted:
            for decoy in alternatives:
                per_modality = {
                    modality: distances(control["library"], recipient,
                                        planted, decoy, modality)
                    for modality in ("RNA", "ATAC")}
                valid = all(
                    item["shared_callable_sites"] >= DECOY_MIN_SHARED_SITES and
                    item["expected_decoy_genotype_distance"] > 0 and
                    item["planted_distance"] > 0 and
                    item["distance_relative_mismatch"] <= DECOY_RELATIVE_TOLERANCE and
                    item["opportunity_relative_mismatch"] <= DECOY_RELATIVE_TOLERANCE
                    for item in per_modality.values())
                score = max(
                    max(item["distance_relative_mismatch"],
                        item["opportunity_relative_mismatch"])
                    for item in per_modality.values())
                decoy_diagnostics.append({
                    "control_id": control["control_id"],
                    "decoy_contributor": decoy["second_state"],
                    "eligible": valid, "maximum_relative_mismatch": score,
                    "coverage_relative_mismatch": 0.0,
                    "coverage_tolerance": DECOY_RELATIVE_TOLERANCE,
                    "opportunity_tolerance": DECOY_RELATIVE_TOLERANCE,
                    **{f"{modality.lower()}_{key}": value
                       for modality, item in per_modality.items()
                       for key, value in item.items()},
                    "distance_formula": (
                        "mean absolute difference in expected alternate-allele "
                        "fraction from recipient locked source over the same shared "
                        "callable selected sites"),
                })
                if valid:
                    eligible.append((score, decoy["second_state"], per_modality))
        if eligible:
            _, decoy_name, metrics = min(eligible, key=lambda item: (item[0], item[1]))
            result["genotype_distance_matched_decoy"] = decoy_name
            result["decoy_status"] = "MATCHED_WITHIN_PREDECLARED_TOLERANCES"
            result["decoy_comparison_eligible"] = True
            for modality, item in metrics.items():
                for key, value in item.items():
                    result[f"{modality.lower()}_{key}"] = value
                result[f"{modality.lower()}_continuous_coverage_mismatch"] = 0.0
            result["continuous_coverage_mismatch_explanation"] = (
                "planted and decoy candidates are evaluated in the same recipient "
                "cell and therefore share identical observed coverage")
            result["decoy_exclusion_reason"] = "NONE"
        else:
            result["genotype_distance_matched_decoy"] = ""
            result["decoy_status"] = "UNMATCHED_EXCLUDED_FROM_DECOY_COMPARISON"
            result["decoy_comparison_eligible"] = False
            result["decoy_exclusion_reason"] = (
                "NO_LEGAL_DECOY_MEETS_SHARED_CALLABLE_SITE_GENOTYPE_DISTANCE_"
                "OPPORTUNITY_AND_COVERAGE_TOLERANCES")
        output_rows.append(result)
    output_rows.sort(key=lambda row: row["control_id"])
    final_path = root / "manifests" / "control_pairs_final.tsv"
    control_ids = [row["control_id"] for row in output_rows]
    if len(control_ids) != len(set(control_ids)) or \
            len(output_rows) != sum(1 for _ in read_tsv(provisional)):
        raise RuntimeError("final control manifest key/count validation failed")
    write_tsv(final_path, output_rows)
    write_tsv(root / "analysis" / "decoy_match_diagnostics.tsv",
              decoy_diagnostics)
    contract = _control_finalizer_contract(root, args.generation,
                                           analysis_tasks)
    marker = {
        "schema_version": TASK_MARKER_SCHEMA,
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "operational_status": "COMPLETE", "utc": utc_now(),
        "workload_generation_id": args.generation,
        "controls": len(output_rows),
        "planted_controls": sum(row["control_class"] == "PLANTED_NEW_SOURCE"
                                for row in output_rows),
        "decoy_matched": sum(truthy(row.get("decoy_comparison_eligible", ""))
                             for row in output_rows),
        "output": str(final_path),
        "output_content_sha256": _small_file_sha256(final_path),
        "implementation_version": IMPLEMENTATION_VERSION,
        "cache_schema": CACHE_SCHEMA,
        "finalizer_configuration": contract,
    }
    atomic_json(root / "markers" / "control_finalizer.complete.json", marker)
    print(json.dumps(marker, indent=2))
    return 0


def compare_equivalence(args):
    def load(path):
        rows = list(read_tsv(path))
        by_key = {}
        for row in rows:
            key = row.get("result_key", "")
            if not key or key in by_key:
                raise RuntimeError(f"{path}: missing or duplicate result_key {key!r}")
            by_key[key] = row
        return rows, by_key
    left_rows, left = load(args.left)
    right_rows, right = load(args.right)
    differences = []
    tolerance = 1e-9
    keys = sorted(set(left) | set(right))
    excluded = {"threads", "wall_seconds", "cpu_seconds", "result_key"}
    for key in keys:
        if key not in left or key not in right:
            differences.append({"result_key": key, "field": "ROW_PRESENCE",
                                "left": key in left, "right": key in right,
                                "within_tolerance": False})
            continue
        for field in sorted(set(left[key]) | set(right[key])):
            if field in excluded:
                continue
            a = clean(left[key].get(field, ""))
            b = clean(right[key].get(field, ""))
            if a == b:
                continue
            a_number = finite(a)
            b_number = finite(b)
            within = math.isfinite(a_number) and math.isfinite(b_number) and \
                abs(a_number - b_number) <= tolerance * max(
                    1.0, abs(a_number), abs(b_number))
            differences.append({
                "result_key": key, "field": field, "left": a, "right": b,
                "absolute_difference": abs(a_number - b_number)
                    if math.isfinite(a_number) and math.isfinite(b_number)
                    else math.nan,
                "within_tolerance": within,
            })
    for difference in differences:
        difference["comparison_source"] = "THREAD_INVARIANCE"
    thread_material = [row for row in differences
                       if not truthy(row["within_tolerance"])]
    optimized_reference_failures = []
    optimized_reference_differences = []
    for name, rows in (("threads1", left_rows), ("production_threads", right_rows)):
        engines = defaultdict(dict)
        for row in rows:
            if row.get("result_class", "") not in {
                    "OBSERVED_FULL_MENU_CANDIDATE",
                    "SOURCE_DISJOINT_CONTROL_CANDIDATE"}:
                continue
            engine = row.get("engine", "")
            if engine in {"OPTIMIZED_BATCHED", "REFERENCE_SCALAR"}:
                comparison_key = row.get("equivalence_key", "")
                engines[comparison_key][engine] = row
        for comparison_key, pair in engines.items():
            if set(pair) != {"OPTIMIZED_BATCHED", "REFERENCE_SCALAR"}:
                optimized_reference_failures.append(
                    f"{name}:{comparison_key}:MISSING_ENGINE")
                optimized_reference_differences.append({
                    "result_key": comparison_key,
                    "field": "ENGINE_PRESENCE", "left": ",".join(sorted(pair)),
                    "right": "OPTIMIZED_BATCHED,REFERENCE_SCALAR",
                    "within_tolerance": False,
                    "comparison_source": f"OPTIMIZED_VS_SCALAR_{name}",
                })
                continue
            optimized_definition = pair["OPTIMIZED_BATCHED"].get(
                "reference_definition", "")
            scalar_definition = pair["REFERENCE_SCALAR"].get(
                "reference_definition", "")
            if optimized_definition != (
                    "normalized-cache compiled batch with shared deterministic "
                    "unit masks") or scalar_definition != (
                    "completed joint-doublet scalar scorer using identical "
                    "cached observations") or optimized_definition == \
                    scalar_definition:
                optimized_reference_failures.append(
                    f"{name}:{comparison_key}:EVALUATOR_CODE_PATH_UNPROVEN")
                optimized_reference_differences.append({
                    "result_key": comparison_key,
                    "field": "reference_definition",
                    "left": optimized_definition, "right": scalar_definition,
                    "within_tolerance": False,
                    "comparison_source": f"OPTIMIZED_VS_SCALAR_{name}",
                })
            for field in ("locked_log_likelihood", "interior_log_likelihood",
                          "contributor_only_log_likelihood", "fitted_fraction",
                          "fitted_fraction_profile_low",
                          "fitted_fraction_profile_high",
                          "delta_log_likelihood", "preferred_model", "status",
                          "evidence_category", "rank", "winner", "runner_up",
                          "winner_margin"):
                a = clean(pair["OPTIMIZED_BATCHED"].get(field, ""))
                b = clean(pair["REFERENCE_SCALAR"].get(field, ""))
                if a == b:
                    continue
                av, bv = finite(a), finite(b)
                if not (math.isfinite(av) and math.isfinite(bv) and
                        abs(av - bv) <= tolerance * max(1.0, abs(av), abs(bv))):
                    optimized_reference_failures.append(
                        f"{name}:{comparison_key}:{field}:{a}!={b}")
                    optimized_reference_differences.append({
                        "result_key": comparison_key, "field": field,
                        "left": a, "right": b,
                        "absolute_difference": abs(av - bv)
                            if math.isfinite(av) and math.isfinite(bv)
                            else math.nan,
                        "within_tolerance": False,
                        "comparison_source": f"OPTIMIZED_VS_SCALAR_{name}",
                    })

    legacy_failures = []
    legacy_differences = []
    legacy_fields_compared = 0
    action = left_rows[0].get("action", "") if left_rows else ""
    if action == "CELL":
        legacy_path = Path(args.legacy_reference)
        if not legacy_path.is_file() or legacy_path.stat().st_size == 0:
            raise RuntimeError(
                "cell benchmark requires the completed-score legacy reference")
        legacy_rows = list(read_tsv(legacy_path))
        legacy = {}
        for row in legacy_rows:
            if row.get("workload_generation_id") != args.generation:
                raise RuntimeError("legacy reference workload generation mismatch")
            key = (row.get("library", ""), row.get("modality", ""),
                   row.get("barcode", ""), row.get("candidate_id", ""),
                   row.get("evidence_channel", ""))
            if key in legacy:
                raise RuntimeError(f"duplicate completed-score legacy key: {key}")
            legacy[key] = row
        mappings = (
            ("status", "legacy_status", False),
            ("locked_log_likelihood", "legacy_locked_log_likelihood", True),
            ("interior_log_likelihood", "legacy_interior_log_likelihood", True),
            ("contributor_only_log_likelihood",
             "legacy_contributor_only_log_likelihood", True),
            ("delta_log_likelihood", "legacy_delta_log_likelihood", True),
            ("fitted_fraction", "legacy_fitted_fraction", True),
            ("fitted_fraction_profile_low",
             "legacy_fitted_fraction_profile_low", True),
            ("fitted_fraction_profile_high",
             "legacy_fitted_fraction_profile_high", True),
            ("evidence_units", "legacy_evidence_units", True),
        )
        optimized_rows = [row for row in left_rows
                          if row.get("result_class") ==
                          "OBSERVED_FULL_MENU_CANDIDATE" and
                          row.get("engine") == "OPTIMIZED_BATCHED"]
        for row in optimized_rows:
            key = (row.get("library", ""), row.get("modality", ""),
                   row.get("barcode", ""), row.get("candidate_id", ""),
                   row.get("evidence_channel", ""))
            reference = legacy.get(key)
            result_key = row.get("equivalence_key", "")
            if reference is None:
                detail = f"{result_key}:MISSING_COMPLETED_SCORE_ROW"
                legacy_failures.append(detail)
                legacy_differences.append({
                    "result_key": result_key, "field": "ROW_PRESENCE",
                    "left": "OPTIMIZED_PRESENT", "right": "LEGACY_MISSING",
                    "within_tolerance": False,
                    "comparison_source": "OPTIMIZED_VS_COMPLETED_SCORE",
                })
                continue
            for optimized_field, legacy_field, numeric in mappings:
                observed_value = clean(row.get(optimized_field, ""))
                legacy_value = clean(reference.get(legacy_field, ""))
                if numeric and not math.isfinite(finite(legacy_value)):
                    continue
                if not numeric and not legacy_value:
                    continue
                legacy_fields_compared += 1
                if observed_value == legacy_value:
                    continue
                observed_number = finite(observed_value)
                legacy_number = finite(legacy_value)
                within = numeric and math.isfinite(observed_number) and \
                    math.isfinite(legacy_number) and abs(
                        observed_number - legacy_number) <= tolerance * max(
                            1.0, abs(observed_number), abs(legacy_number))
                if within:
                    continue
                legacy_failures.append(
                    f"{result_key}:{optimized_field}:"
                    f"{observed_value}!={legacy_value}")
                legacy_differences.append({
                    "result_key": result_key, "field": optimized_field,
                    "left": observed_value, "right": legacy_value,
                    "absolute_difference": abs(
                        observed_number - legacy_number)
                        if math.isfinite(observed_number) and
                        math.isfinite(legacy_number) else math.nan,
                    "within_tolerance": False,
                    "comparison_source": "OPTIMIZED_VS_COMPLETED_SCORE",
                })
    differences.extend(optimized_reference_differences)
    differences.extend(legacy_differences)
    material = [row for row in differences if not truthy(row["within_tolerance"])]
    status = "COMPLETE" if not material and not optimized_reference_failures and \
        not legacy_failures \
        else "FAILED"
    output_rows = []
    for index, difference in enumerate(differences):
        detail = dict(difference)
        compared_result_key = detail.pop("result_key", "")
        output_rows.append({
            "schema_version": "joint_doublet_benchmark_equivalence_v4",
            "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
            "calibration_library": CALIBRATION_LIBRARY,
            "workload_generation_id": args.generation,
            "task_index": args.task_index,
            "result_key": f"equivalence:{args.task_index}:{index}",
            "result_class": "THREAD_AND_REFERENCE_EQUIVALENCE",
            "action": left_rows[0].get("action", "") if left_rows else "",
            "role": left_rows[0].get("role", "") if left_rows else "",
            "library": left_rows[0].get("library", "") if left_rows else "",
            "modality": left_rows[0].get("modality", "") if left_rows else "",
            "compared_result_key": compared_result_key,
            **detail,
        })
    if not output_rows:
        output_rows = [{
            "schema_version": "joint_doublet_benchmark_equivalence_v4",
            "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
            "calibration_library": CALIBRATION_LIBRARY,
            "workload_generation_id": args.generation,
            "task_index": args.task_index,
            "result_key": f"equivalence:{args.task_index}:PASS",
            "result_class": "THREAD_AND_REFERENCE_EQUIVALENCE",
            "action": left_rows[0].get("action", "") if left_rows else "",
            "role": left_rows[0].get("role", "") if left_rows else "",
            "library": left_rows[0].get("library", "") if left_rows else "",
            "modality": left_rows[0].get("modality", "") if left_rows else "",
            "field": "ALL_SCIENTIFIC_FIELDS", "left": "MATCH",
            "right": "MATCH", "within_tolerance": True,
        }]
    write_tsv(args.output, output_rows)
    difference_rows = differences or [{
        "result_key": "ALL_COMPARED_RESULTS",
        "field": "ALL_SCIENTIFIC_FIELDS",
        "left": "MATCH", "right": "MATCH",
        "absolute_difference": 0.0,
        "within_tolerance": True,
    }]
    write_tsv(str(args.output) + ".field_differences.tsv", difference_rows)
    manifest_rows = list(read_tsv(args.task_manifest))
    task_index = int(args.task_index)
    if task_index < 0 or task_index >= len(manifest_rows) or \
            int(manifest_rows[task_index].get("task_index", -1)) != task_index:
        raise RuntimeError("equivalence task manifest/index mismatch")
    marker = {
        "schema_version": TASK_MARKER_SCHEMA,
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "operational_status": status, "utc": utc_now(),
        "task_kind": "analysis", "workload_generation_id": args.generation,
        "task_index": int(args.task_index),
        "library": output_rows[0]["library"],
        "modality": output_rows[0]["modality"],
        "action": output_rows[0]["action"],
        "implementation_version": IMPLEMENTATION_VERSION,
        "cache_schema": CACHE_SCHEMA,
        "task_configuration": manifest_rows[task_index],
        "audited_output_files": _file_states([Path(args.output)]),
        "audited_dependencies": _task_dependency_states(manifest_rows[task_index], "analysis"),
        "numerical_tolerance": tolerance,
        "thread_material_differences": len(thread_material),
        "optimized_reference_failures": optimized_reference_failures,
        "legacy_completed_score_fields_compared": legacy_fields_compared,
        "legacy_completed_score_failures": legacy_failures,
        "legacy_reference": args.legacy_reference,
        "legacy_reference_content_sha256": _small_file_sha256(
            args.legacy_reference) if args.legacy_reference else "",
        "field_difference_output": str(args.output) +
            ".field_differences.tsv",
        "single_thread_output": str(args.left),
        "production_thread_output": str(args.right),
        "output_content_sha256": _small_file_sha256(args.output),
        "field_difference_content_sha256": _small_file_sha256(
            str(args.output) + ".field_differences.tsv"),
        "single_thread_content_sha256": _small_file_sha256(args.left),
        "production_thread_content_sha256": _small_file_sha256(args.right),
    }
    atomic_json(args.marker, marker)
    return 0 if status == "COMPLETE" else 2


def _read_submission_ledger(root):
    path = Path(root) / "submission_ledger.jsonl"
    if not path.is_file():
        return []
    rows = []
    for line in path.read_text().splitlines():
        if line.strip():
            rows.append(json.loads(line))
    return rows


def _append_submission(root, record):
    path = Path(root) / "submission_ledger.jsonl"
    rows = _read_submission_ledger(root)
    rows.append(record)
    atomic_text(path, "".join(json.dumps(row, sort_keys=True) + "\n"
                              for row in rows))


def _active_scheduler_jobs(job_ids):
    job_ids = sorted({str(job_id) for job_id in job_ids if str(job_id).isdigit()})
    if not job_ids:
        return {}, "NONE"
    try:
        result = subprocess.run(
            ["squeue", "-h", "-j", ",".join(job_ids), "-o", "%A|%T"],
            capture_output=True, text=True, check=False)
    except OSError as error:
        return None, f"SCHEDULER_QUERY_FAILED:{error}"
    if result.returncode:
        return None, "SCHEDULER_QUERY_FAILED:" + (
            result.stderr.strip() or result.stdout.strip() or
            f"returncode={result.returncode}")
    return ({line.split("|", 1)[0]: line.split("|", 1)[1]
             for line in result.stdout.splitlines() if "|" in line}, "OK")


def _resolve_targeted_root(root, scope):
    root = Path(root).resolve()
    if scope == "benchmark":
        candidate = (root / "benchmark_unsubmitted").resolve()
        candidate.relative_to(root)
        return candidate, root
    current = root / "CURRENT_GENERATION.json"
    if current.is_file():
        payload = json.loads(current.read_text())
        if payload.get("scientific_method_version") != \
                SCIENTIFIC_METHOD_VERSION or int(payload.get(
                    "calibration_library", -1)) != CALIBRATION_LIBRARY:
            raise RuntimeError(
                "CURRENT_GENERATION method/calibration mismatch")
        candidate = Path(payload.get("root", ""))
        if not candidate.is_absolute():
            candidate = (root / candidate).resolve()
        else:
            candidate = candidate.resolve()
        try:
            candidate.relative_to(root)
        except ValueError as error:
            raise RuntimeError("CURRENT_GENERATION root escapes workload root") from error
        generation = clean(payload.get("workload_generation_id", ""))
        if not generation or not re.fullmatch(r"(?:prebenchmark|reprojected)_[0-9a-f]{16}",
                                               generation):
            raise RuntimeError("CURRENT_GENERATION identifier is invalid")
        if generation.startswith("prebenchmark_"):
            if candidate != root:
                raise RuntimeError(
                    "prebenchmark CURRENT_GENERATION must resolve to workload root")
        elif candidate == root or candidate.parent != \
                (root / "generations").resolve() or candidate.name != generation:
            raise RuntimeError(
                "reprojected CURRENT_GENERATION root/generation mismatch")
        blueprint_path = candidate / "workload_blueprint.json"
        if not candidate.is_dir() or not blueprint_path.is_file():
            raise RuntimeError("CURRENT_GENERATION root lacks workload blueprint")
        blueprint = json.loads(blueprint_path.read_text())
        if blueprint.get("workload_generation_id") != generation:
            raise RuntimeError("CURRENT_GENERATION blueprint mismatch")
        if blueprint.get("scientific_method_version") != \
                SCIENTIFIC_METHOD_VERSION or int(blueprint.get(
                    "calibration_library", -1)) != CALIBRATION_LIBRARY:
            raise RuntimeError(
                "CURRENT_GENERATION blueprint method/calibration mismatch")
        return candidate, root
    return root, root


def _control_finalizer_contract(root, generation, analysis_tasks):
    root = Path(root)
    provisional = root / "manifests" / "control_pairs_provisional.tsv"
    extraction = root / "manifests" / "extraction_tasks.tsv"
    # The finalizer publishes the complete provisional control population.
    # A benchmark may deliberately analyze only one bounded control bin, so
    # using analysis-shard membership as the expected final-manifest set would
    # incorrectly reject a complete finalizer output.
    provisional_rows = list(read_tsv(provisional)) if provisional.is_file() else []
    provisional_ids = sorted(row.get("control_id", "") for row in provisional_rows)
    extraction_rows = list(read_tsv(extraction)) if extraction.is_file() else []
    analysis_control_ids = sorted({
        identifier for row in analysis_tasks
        if row.get("action") == "CONTROL_BIN"
        for identifier in row.get("control_ids", "").split(",") if identifier
    })
    payload = {
        "implementation_version": IMPLEMENTATION_VERSION,
        "cache_schema": CACHE_SCHEMA,
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "workload_generation_id": generation,
        "control_ids": provisional_ids,
        "analysis_control_ids": analysis_control_ids,
        "provisional_controls": provisional_rows,
        "extraction_tasks": extraction_rows,
        "audited_cache_dependencies": _file_states([
            Path(row["cache_prefix"] + suffix) for row in extraction_rows
            for suffix in (".metadata.json", ".genotypes.tsv.gz", ".samples.tsv")]),
    }
    return payload


def _validate_control_finalizer(root, generation, analysis_tasks):
    if not any(row.get("action") == "CONTROL_BIN" for row in analysis_tasks):
        return True, []
    marker_path = Path(root) / "markers" / "control_finalizer.complete.json"
    output = Path(root) / "manifests" / "control_pairs_final.tsv"
    problems = []
    try:
        marker = json.loads(marker_path.read_text())
        if marker.get("schema_version") != TASK_MARKER_SCHEMA or \
                marker.get("implementation_version") != IMPLEMENTATION_VERSION or \
                marker.get("cache_schema") != CACHE_SCHEMA:
            problems.append("CONTROL_FINALIZER_SCHEMA_OR_IMPLEMENTATION_MISMATCH")
        if marker.get("operational_status") != "COMPLETE":
            problems.append("CONTROL_FINALIZER_NOT_COMPLETE")
        if marker.get("scientific_method_version") != \
                SCIENTIFIC_METHOD_VERSION or marker.get(
                    "calibration_library") != CALIBRATION_LIBRARY:
            problems.append("CONTROL_FINALIZER_METHOD_OR_CALIBRATION_MISMATCH")
        if marker.get("workload_generation_id") != generation:
            problems.append("CONTROL_FINALIZER_GENERATION_MISMATCH")
        contract = _control_finalizer_contract(root, generation, analysis_tasks)
        if marker.get("finalizer_configuration") != contract:
            problems.append("CONTROL_FINALIZER_CONTRACT_MISMATCH")
        if not output.is_file() or output.stat().st_size == 0:
            problems.append("FINAL_CONTROL_MANIFEST_MISSING")
        else:
            if marker.get("output_content_sha256") != _small_file_sha256(output):
                problems.append("FINAL_CONTROL_MANIFEST_DIGEST_MISMATCH")
            rows = list(read_tsv(output))
            required = {
                "schema_version", "scientific_method_version",
                "calibration_library", "workload_generation_id", "control_id",
                "control_class", "library", "recipient_barcode",
                "source_barcode", "source_identity", "expected_contributor",
                "requested_fractions", "seed", "decoy_status",
                "decoy_comparison_eligible", "decoy_exclusion_reason",
                "rna_molecule_evidence_basis",
                "atac_molecule_evidence_basis",
                "source_rna_molecule_evidence_basis",
                "source_atac_molecule_evidence_basis",
            }
            if not rows or not required.issubset(rows[0]):
                problems.append("FINAL_CONTROL_MANIFEST_SCHEMA_INVALID")
            identifiers = [row.get("control_id", "") for row in rows]
            if identifiers != sorted(identifiers) or any(not value for value in
                                                         identifiers) or \
                    len(identifiers) != len(set(identifiers)):
                problems.append("FINAL_CONTROL_MANIFEST_KEYS_INVALID")
            if set(identifiers) != set(contract["control_ids"]):
                problems.append("FINAL_CONTROL_MANIFEST_CONTROL_SET_MISMATCH")
            for row in rows:
                if row.get("schema_version") != \
                        "joint_doublet_control_manifest_v4" or \
                        row.get("scientific_method_version") != \
                        SCIENTIFIC_METHOD_VERSION or row.get(
                            "calibration_library") != str(CALIBRATION_LIBRARY) or \
                        row.get("workload_generation_id") != generation:
                    problems.append("FINAL_CONTROL_MANIFEST_PROVENANCE_MISMATCH")
                    break
                if row.get("library") not in {
                        f"lib{value}" for value in ALLOWED_LIBRARIES}:
                    problems.append("FINAL_CONTROL_MANIFEST_LIBRARY_INVALID")
                    break
                if row.get("control_class") == "PLANTED_NEW_SOURCE":
                    eligible = truthy(row.get("decoy_comparison_eligible", ""))
                    if eligible and (not row.get(
                            "genotype_distance_matched_decoy") or
                            row.get("decoy_exclusion_reason") != "NONE"):
                        problems.append("FINAL_CONTROL_DECOY_CONTRACT_INVALID")
                        break
                    if not eligible and not row.get("decoy_exclusion_reason"):
                        problems.append("FINAL_CONTROL_DECOY_EXCLUSION_MISSING")
                        break
            if int(marker.get("controls", -1)) != len(rows):
                problems.append("CONTROL_FINALIZER_COUNT_MISMATCH")
    except (OSError, ValueError, TypeError) as error:
        problems.append(f"CONTROL_FINALIZER_MARKER_INVALID:{error}")
    return not problems, problems


def targeted_control(args):
    run_root, _ = _resolve_targeted_root(args.targeted_root, args.scope)
    top_root = Path(args.targeted_root).resolve()
    extraction = list(read_tsv(run_root / "manifests" / "extraction_tasks.tsv"))
    analysis = list(read_tsv(run_root / "manifests" / "analysis_tasks.tsv"))
    if extraction:
        validate_targeted_input_contract(
            run_root / "manifests" / "extraction_tasks.tsv", "extraction", 0)
    if analysis:
        validate_targeted_input_contract(
            run_root / "manifests" / "analysis_tasks.tsv", "analysis", 0)
    validate_targeted_manifest_boundary(extraction, "extraction")
    validate_targeted_manifest_boundary(analysis, "analysis")
    generation = analysis[0]["workload_generation_id"] if analysis else \
        extraction[0]["workload_generation_id"] if extraction else ""
    ledger = [row for row in _read_submission_ledger(run_root)
              if row.get("workload_generation_id") == generation]
    active, scheduler_query = _active_scheduler_jobs(
        row.get("job_id", "") for row in ledger)
    if active is None:
        payload = {
            "operational_status": "TECHNICAL_UNRESOLVED_SCHEDULER_QUERY",
            "scientific_status": "UNASSESSED", "scope": args.scope,
            "run_root": str(run_root), "workload_generation_id": generation,
            "scheduler_query": scheduler_query,
            "submission_authorized": False,
        }
        if args.operation == "status":
            print(json.dumps(payload, indent=2))
            return 2
        raise RuntimeError(
            "scheduler query failed; resubmission is prohibited until a successful check")
    extraction_status = {}
    analysis_status = {}
    for row in extraction:
        valid, problems, _ = validate_task_artifact(row, "extraction", True, deep=False)
        extraction_status[int(row["task_index"])] = (valid, problems)
    for row in analysis:
        valid, problems, _ = validate_task_artifact(row, "analysis", True, deep=False)
        analysis_status[int(row["task_index"])] = (valid, problems)
    finalizer_valid, finalizer_problems = _validate_control_finalizer(
        run_root, generation, analysis)
    invalid_sources = {(row["library"], row["modality"]) for row in extraction
                       if not extraction_status[int(row["task_index"])][0]}
    if invalid_sources:
        finalizer_valid = False
        finalizer_problems.append("UPSTREAM_CACHE_INCOMPLETE")
    for row in analysis:
        if (row["library"], row["modality"]) in invalid_sources or (row["action"] == "CONTROL_BIN" and not finalizer_valid):
            analysis_status[int(row["task_index"])] = (False, ["UPSTREAM_DEPENDENCY_INCOMPLETE"])

    gather_marker = run_root / "TARGETED_VALIDATION_FINISHED"
    try:
        gather = json.loads(gather_marker.read_text()) if gather_marker.is_file() else {}
    except (OSError, ValueError):
        gather = {}
    combined = Path(gather.get("combined_output", "")) \
        if gather.get("combined_output") else None
    completion_problems = _targeted_completion_problems(
        run_root, args.scope == "benchmark", deep=False) if gather else [
            "COMPLETION_MARKER_MISSING"]
    gather_valid = not completion_problems
    active_states = Counter(active.values())
    if gather_valid:
        operational = "EXECUTION_COMPLETE"
    elif active_states.get("RUNNING") or active_states.get("COMPLETING"):
        operational = "RUNNING"
    elif active:
        operational = "QUEUED"
    elif not ledger:
        operational = "RENDERED_UNSUBMITTED"
    elif all(valid for valid, _ in extraction_status.values()) and \
            all(valid for valid, _ in analysis_status.values()) and \
            finalizer_valid:
        operational = "EXECUTION_COMPLETE_PENDING_GATHER"
    else:
        operational = "TECHNICAL_FAILURE_OR_INCOMPLETE"
    payload = {
        "operational_status": operational,
        "scientific_status": gather.get("scientific_status", "UNASSESSED"),
        "scope": args.scope, "run_root": str(run_root),
        "workload_generation_id": generation,
        "extraction": {str(key): {"valid": value[0], "problems": value[1]}
                       for key, value in extraction_status.items()},
        "analysis": {str(key): {"valid": value[0], "problems": value[1]}
                     for key, value in analysis_status.items()},
        "control_finalizer": {"valid": finalizer_valid,
                              "problems": finalizer_problems},
        "active_scheduler_states": dict(active_states),
        "scheduler_query": scheduler_query,
        "submission_records": len(ledger), "gather_valid": gather_valid,
        "completion_problems": completion_problems,
        "gather": gather,
    }
    if args.operation == "status":
        job_ids = sorted({row.get("job_id", "") for row in ledger
                          if str(row.get("job_id", "")).isdigit()})
        if job_ids:
            try:
                sacct = subprocess.run([
                    "sacct", "-P", "-n", "-j", ",".join(job_ids),
                    "--format=JobID,JobName,State,ExitCode,Elapsed,TotalCPU,AllocCPUS,MaxRSS,ReqMem,NodeList"],
                    capture_output=True, text=True, check=False)
                payload["sacct_snapshot"] = sacct.stdout.strip() \
                    if sacct.returncode == 0 else sacct.stderr.strip()
            except OSError as error:
                payload["sacct_snapshot"] = f"UNAVAILABLE:{error}"
        print(json.dumps(payload, indent=2))
        return 0
    if not args.submit:
        raise RuntimeError(f"targeted {args.operation} requires explicit --submit")
    for name in ("logs", "task_scratch", "scores", "cache", "analysis", "markers"):
        (run_root / name).mkdir(parents=True, exist_ok=True)
    scripts = {
        "extract": run_root / "slurm_scripts" / "01_normalized_extract.sbatch",
        "finalize": run_root / "slurm_scripts" / "02_finalize_controls.sbatch",
        "analysis": run_root / "slurm_scripts" / "03_batched_analysis.sbatch",
        "gather": run_root / "slurm_scripts" / "04_gather.sbatch",
    }
    resource_projection_path = run_root / "resource_projection.json"
    resource_projection = json.loads(resource_projection_path.read_text()) \
        if resource_projection_path.is_file() else {}
    array_concurrency = max(1, int(resource_projection.get(
        "array_concurrency", 4)))
    active_records = [row for row in ledger
                      if str(row.get("job_id", "")) in active]
    extraction_jobs = {}
    submitted = {}
    pending_extract = []
    for task in extraction:
        index = int(task["task_index"])
        if extraction_status[index][0]:
            continue
        if clean(task.get("cache_policy", "")) == "IMMUTABLE_BENCHMARK_REUSE_REQUIRED":
            raise RuntimeError(f"immutable promoted cache invalid for extraction task {index}")
        existing = next((row for row in active_records if row.get("node") == "extraction" and index in row.get("task_indices", [])), None)
        if existing:
            extraction_jobs[(task["library"],task["modality"])] = existing["job_id"]
        else:
            pending_extract.append(task)
    if pending_extract:
        indices = [int(row["task_index"]) for row in pending_extract]
        active_dependencies = sorted(set(extraction_jobs.values()))
        dependency = "afterany:"+":".join(active_dependencies) if active_dependencies else ""
        job_id = submit_sbatch(scripts["extract"],run_root,dependency=dependency,
                               array=",".join(map(str,indices))+f"%{array_concurrency}")
        _append_submission(run_root,{"utc":utc_now(),"workload_generation_id":generation,
            "node":"extraction","job_id":job_id,"task_indices":indices,"dependency":dependency})
        submitted["extraction"] = job_id
        for task in pending_extract:
            extraction_jobs[(task["library"],task["modality"])] = job_id
    finalizer_job = ""
    # The finalizer is an independently repairable graph node.  A corrupt or
    # missing final control manifest must be retried even when every upstream
    # extraction/analysis artifact happens to be valid (for example after an
    # interrupted atomic publication or manual corruption).
    if any(row["action"] == "CONTROL_BIN" for row in analysis) and \
            not finalizer_valid:
        existing = next((row for row in active_records
                         if row.get("node") == "control_finalizer"), None)
        if existing:
            finalizer_job = existing["job_id"]
        else:
            dependencies = sorted(set(extraction_jobs.values()))
            dependency = "afterany:" + ":".join(dependencies) \
                if dependencies else ""
            finalizer_job = submit_sbatch(
                scripts["finalize"], run_root, dependency=dependency)
            _append_submission(run_root, {
                "utc": utc_now(), "workload_generation_id": generation,
                "node": "control_finalizer", "job_id": finalizer_job,
                "task_indices": [], "dependency": dependency})
            submitted["control_finalizer"] = finalizer_job
    analysis_jobs = []
    grouped = defaultdict(list)
    for task in analysis:
        index = int(task["task_index"])
        if analysis_status[index][0]:
            continue
        if any(row.get("node") == "analysis" and
               index in row.get("task_indices", []) for row in active_records):
            continue
        grouped[("ALL", "ALL", "ALL")].append(index)
    for (library, modality, action), indices in sorted(grouped.items()):
        dependencies = list(extraction_jobs.values())
        dependencies.extend(row["job_id"] for row in active_records if row.get("node") == "analysis")
        if finalizer_job:
            dependencies.append(finalizer_job)
        # Run after every prerequisite reaches a terminal state, then let the
        # shared artifact validator decide whether this shard is executable.
        # This prevents a failed extraction/finalizer from leaving dependent
        # arrays (and the afterany gather) permanently dependency-blocked.
        dependency = "afterany:" + ":".join(sorted(set(dependencies))) \
            if dependencies else ""
        for offset in range(0, len(indices), 999):
            shard_indices = indices[offset:offset+999]
            mapping = run_root / "manifests" / f"submit_{library}_{modality}_{action}_{offset}_{len(ledger)+len(submitted)}.indices.txt"
            atomic_text(mapping, "\n".join(map(str, shard_indices)) + "\n")
            array = f"0-{len(shard_indices)-1}%{array_concurrency}"
            job_id = submit_sbatch(
                scripts["analysis"], run_root, dependency=dependency, array=array, index_map=str(mapping))
            analysis_jobs.append(job_id)
            # At most one new analysis array is active at a time, including when
            # the workload is split across the scheduler's 999-element limit.
            next_dependency = "afterany:" + job_id
            _append_submission(run_root, {
                "utc": utc_now(), "workload_generation_id": generation,
                "node": "analysis", "job_id": job_id,
                "task_indices": shard_indices, "dependency": dependency})
            submitted[f"analysis_{library}_{modality}_{action}_{offset}"] = job_id
            dependency = next_dependency
    gather_active = next((row for row in active_records
                          if row.get("node") == "gather"), None)
    # Likewise, repair the terminal gather itself without requiring a newly
    # submitted analysis shard.  validate_task_artifact() above determines the
    # precise upstream set; valid shards are never replayed merely to repair a
    # stale/malformed gather.
    if not gather_active and not gather_valid:
        parents = analysis_jobs + [
            row["job_id"] for row in active_records if row.get("node") == "analysis"]
        if finalizer_job:
            parents.append(finalizer_job)
        dependency = "afterany:" + ":".join(sorted(set(parents))) \
            if parents else ""
        gather_job = submit_sbatch(
            scripts["gather"], run_root, dependency=dependency)
        _append_submission(run_root, {
            "utc": utc_now(), "workload_generation_id": generation,
            "node": "gather", "job_id": gather_job,
            "task_indices": [], "dependency": dependency})
        submitted["gather"] = gather_job
    print(json.dumps({"submitted_jobs": submitted,
                      "workload_generation_id": generation}, indent=2))
    return 0


def _parse_time_v(path):
    metrics = {}
    for line in Path(path).read_text(errors="replace").splitlines():
        stripped = line.strip()
        if stripped.startswith("Elapsed (wall clock) time"):
            match = re.search(r"\):\s*([0-9:.]+)$", stripped)
            if match:
                parts = [float(item) for item in match.group(1).split(":")]
                metrics["elapsed_seconds"] = sum(
                    item * 60 ** index
                    for index, item in enumerate(reversed(parts)))
            continue
        if ":" not in line:
            continue
        name, value = [item.strip() for item in line.split(":", 1)]
        if name == "Maximum resident set size (kbytes)":
            metrics["max_rss_kb"] = int(value)
        elif name == "User time (seconds)":
            metrics["user_seconds"] = float(value)
        elif name == "System time (seconds)":
            metrics["system_seconds"] = float(value)
        elif name == "Percent of CPU this job got":
            metrics["cpu_utilization_percent"] = finite(
                value.rstrip("%"))
        elif name == "File system inputs":
            metrics["filesystem_input_blocks"] = int(value)
        elif name == "File system outputs":
            metrics["filesystem_output_blocks"] = int(value)
    return metrics


def _parse_process_interval(path):
    rows = list(read_tsv(path))
    if len(rows) != 1:
        raise RuntimeError(f"invalid process interval record: {path}")
    try:
        start = int(rows[0]["start_ns"])
        end = int(rows[0]["end_ns"])
        status = int(rows[0]["status"])
    except (KeyError, TypeError, ValueError) as error:
        raise RuntimeError(f"invalid process interval values: {path}: {error}")
    if start <= 0 or end < start:
        raise RuntimeError(f"invalid process interval bounds: {path}")
    return {"process_start_ns": start, "process_end_ns": end,
            "process_exit_status": status}


def _measured_io_capacity(rows):
    """Return peak overlapping logical throughput and maximum task demand."""
    intervals = []
    maximum_task = 0.0
    for row in rows:
        demand = finite(row.get("logical_io_mib_per_second", ""))
        start = int(row.get("process_start_ns", 0) or 0)
        end = int(row.get("process_end_ns", 0) or 0)
        if not math.isfinite(demand) or demand <= 0 or start <= 0 or end < start:
            continue
        maximum_task = max(maximum_task, demand)
        intervals.append((start, 1, demand))
        intervals.append((end, -1, demand))
    # End events precede starts at identical timestamps, so touching but
    # non-overlapping processes are not falsely counted as concurrent.
    active_bandwidth = 0.0
    active_count = 0
    maximum_bandwidth = 0.0
    maximum_count = 0
    for _, kind, demand in sorted(intervals, key=lambda item: (item[0], item[1])):
        if kind < 0:
            active_bandwidth -= demand
            active_count -= 1
        else:
            active_bandwidth += demand
            active_count += 1
            if (active_bandwidth, active_count) > (
                    maximum_bandwidth, maximum_count):
                maximum_bandwidth = active_bandwidth
                maximum_count = active_count
    return maximum_bandwidth, maximum_task, maximum_count


def _resource_concurrency_caps(memory_budget_gib, task_memory_gib,
                               usable_cpus, task_cpus,
                               safe_io_mib_per_second,
                               task_io_mib_per_second, hard_cap=8):
    values = (memory_budget_gib, task_memory_gib, usable_cpus, task_cpus,
              safe_io_mib_per_second, task_io_mib_per_second, hard_cap)
    if any(not math.isfinite(float(value)) or float(value) <= 0
           for value in values):
        raise RuntimeError(
            "memory, CPU, I/O, and hard concurrency limits must be positive")
    caps = {
        "MEMORY": int(float(memory_budget_gib) // float(task_memory_gib)),
        "CPU": int(float(usable_cpus) // float(task_cpus)),
        "IO": int(float(safe_io_mib_per_second) //
                  float(task_io_mib_per_second)),
        "HARD_CAP": int(hard_cap),
    }
    if any(value < 1 for value in caps.values()):
        raise RuntimeError("measured resource capacity cannot run one task")
    selected = min(caps.values())
    limiting = sorted(key for key, value in caps.items() if value == selected)
    return selected, caps, limiting


def _benchmark_measurement_problems(root, generation, extraction_tasks,
                                    analysis_tasks, measurements):
    """Require the entire prescribed benchmark panel, not a convenient subset."""
    problems = []
    try:
        blueprint = json.loads((Path(root) / "workload_blueprint.json").read_text())
        resources = _validated_targeted_resources(
            blueprint.get("resource_policy", {}))
        extraction_threads = resources["extraction_cpus"]
        production_threads = resources["analysis_cpus"]
    except (OSError, json.JSONDecodeError, RuntimeError, ValueError) as error:
        return [f"BENCHMARK_RESOURCE_CONTRACT_INVALID:{error}"]
    expected = {
        ("EXTRACTION", int(task["task_index"]),
         task["modality"], extraction_threads)
        for task in extraction_tasks
    }
    expected.update({
        ("ANALYSIS", int(task["task_index"]), task["modality"], threads)
        for task in analysis_tasks for threads in (1, production_threads)
    })
    observed = []
    for row in measurements:
        try:
            key = (row.get("stage", ""), int(row.get("task_index", -1)),
                   row.get("modality", ""), int(row.get("threads", 0)))
        except (TypeError, ValueError):
            problems.append("BENCHMARK_MEASUREMENT_KEY_INVALID")
            continue
        observed.append(key)
        if row.get("workload_generation_id") != generation or row.get(
                "scientific_method_version") != SCIENTIFIC_METHOD_VERSION or \
                row.get("calibration_library") != str(CALIBRATION_LIBRARY):
            problems.append(f"BENCHMARK_MEASUREMENT_GENERATION_MISMATCH:{key}")
        for field in ("elapsed_seconds", "max_rss_kb", "user_seconds",
                      "system_seconds", "cpu_utilization_percent",
                      "filesystem_input_blocks", "filesystem_output_blocks",
                      "process_start_ns", "process_end_ns",
                      "process_exit_status"):
            if not math.isfinite(finite(row.get(field, ""))):
                problems.append(f"BENCHMARK_MEASUREMENT_MISSING:{key}:{field}")
        if int(row.get("process_exit_status", -1) or -1) != 0 or \
                int(row.get("process_end_ns", 0) or 0) < int(
                    row.get("process_start_ns", 0) or 0):
            problems.append(f"BENCHMARK_PROCESS_INTERVAL_INVALID:{key}")
        if key[0] == "EXTRACTION":
            for field in ("measured_source_bytes", "measured_selected_cells",
                          "measured_candidate_rows",
                          "measured_unique_evidence_records",
                          "measured_uncompressed_evidence_bytes",
                          "written_cache_bytes", "measured_source_scans",
                          "logical_io_bytes",
                          "logical_io_mib_per_second"):
                if finite(row.get(field, ""), -1) < 0:
                    problems.append(
                        f"BENCHMARK_EXTRACTION_MEASUREMENT_MISSING:{key}:{field}")
            if int(row.get("measured_source_scans", 0) or 0) != 3:
                problems.append(f"BENCHMARK_SOURCE_SCAN_COUNT_INVALID:{key}")
        else:
            for field in ("output_bytes", "optimized_fit_calls",
                          "reference_fit_calls", "optimizer_calls",
                          "likelihood_evaluations", "derivative_evaluations",
                          "row_scans", "candidate_owned_raw_evidence_copies",
                          "per_fraction_string_rehashes",
                          "per_fraction_string_sorts",
                          "unchanged_80_pass_refits",
                          "predicted_likelihood_work"):
                if finite(row.get(field, ""), -1) < 0:
                    problems.append(
                        f"BENCHMARK_ANALYSIS_MEASUREMENT_MISSING:{key}:{field}")
            if int(row.get("candidate_owned_raw_evidence_copies", -1) or -1) != 0 or \
                    int(row.get("per_fraction_string_rehashes", -1) or -1) != 0 or \
                    int(row.get("per_fraction_string_sorts", -1) or -1) != 0 or \
                    int(row.get("unchanged_80_pass_refits", -1) or -1) != 0:
                problems.append(
                    f"BENCHMARK_COMPILED_EFFICIENCY_INVARIANT_FAILED:{key}")
            if key[3] == production_threads and production_threads > 1 and \
                    finite(row.get("cpu_utilization_percent", ""), 0) <= 125:
                problems.append(
                    f"BENCHMARK_PARALLEL_CPU_UTILIZATION_INSUFFICIENT:{key}")
        matching_tasks = extraction_tasks if key[0] == "EXTRACTION" else \
            analysis_tasks
        task = next((item for item in matching_tasks
                     if int(item.get("task_index", -1)) == key[1] and
                     item.get("modality") == key[2]), None)
        if task is not None:
            budget = int(task.get("worker_memory_budget_bytes", 0) or 0)
            measured_rss = int(finite(row.get("max_rss_kb", ""), 0) * 1024)
            if budget <= 0 or measured_rss > budget:
                problems.append(
                    f"BENCHMARK_WORKER_MEMORY_BUDGET_EXCEEDED:{key}:"
                    f"{measured_rss}>{budget}")
    if len(observed) != len(set(observed)):
        problems.append("BENCHMARK_MEASUREMENT_KEYS_DUPLICATED")
    missing = sorted(expected - set(observed))
    unexpected = sorted(set(observed) - expected)
    if missing:
        problems.append("BENCHMARK_MEASUREMENTS_MISSING:" + repr(missing))
    if unexpected:
        problems.append("BENCHMARK_MEASUREMENTS_UNEXPECTED:" + repr(unexpected))
    if {task.get("modality") for task in extraction_tasks} != {"RNA", "ATAC"}:
        problems.append("BENCHMARK_EXTRACTION_MODALITIES_INCOMPLETE")
    roles = {task.get("role") for task in analysis_tasks
             if task.get("action") == "CELL"}
    if not {"FROZEN_TARGET", "MATCHED_COMPARISON"}.issubset(roles):
        problems.append("BENCHMARK_TARGET_COMPARISON_PANEL_INCOMPLETE")
    if not any(task.get("action") == "CONTROL_BIN" for task in analysis_tasks):
        problems.append("BENCHMARK_CONTROL_PANEL_MISSING")
    branch_path = Path(root) / "manifests" / "benchmark_branch_coverage.tsv"
    try:
        branch_rows = list(read_tsv(branch_path))
        required_branches = {
            "STRICT_SITE_MIXTURE", "UNCERTAIN_ADDITION",
            "REPLACEMENT_OR_BOUNDARY", "SMALLEST_COMPLETE_MENU",
            "LARGEST_COMPLETE_MENU", "FROZEN_TARGET",
            "MATCHED_COMPARISON", "HIGHEST_WORK_CONTROL_SHARD",
        }
        observed_branches = {row.get("branch", "") for row in branch_rows}
        if not required_branches.issubset(observed_branches) or any(
                row.get("workload_generation_id") != generation or
                row.get("coverage_status") != "SELECTED"
                for row in branch_rows):
            problems.append("BENCHMARK_BRANCH_COVERAGE_INCOMPLETE")
        for modality in ("RNA", "ATAC"):
            if not any(branch.startswith(f"{modality}_MOLECULE_BASIS_")
                       for branch in observed_branches):
                problems.append(
                    f"BENCHMARK_{modality}_MOLECULE_BASIS_COVERAGE_MISSING")
    except (OSError, csv.Error) as error:
        problems.append(f"BENCHMARK_BRANCH_COVERAGE_INVALID:{error}")
    for task in analysis_tasks:
        marker_path = Path(task.get("marker", ""))
        try:
            marker = json.loads(marker_path.read_text())
        except (OSError, json.JSONDecodeError) as error:
            problems.append(
                f"BENCHMARK_EQUIVALENCE_MARKER_INVALID:{task.get('task_index')}:{error}")
            continue
        if marker.get("workload_generation_id") != generation or \
                marker.get("scientific_method_version") != \
                SCIENTIFIC_METHOD_VERSION or marker.get(
                    "calibration_library") != CALIBRATION_LIBRARY or \
                marker.get("operational_status") != "COMPLETE" or \
                marker.get("thread_material_differences") != 0 or \
                marker.get("optimized_reference_failures") or \
                marker.get("legacy_completed_score_failures"):
            problems.append(
                f"BENCHMARK_EQUIVALENCE_OR_INVARIANCE_FAILED:{task.get('task_index')}")
    for task in extraction_tasks:
        metadata_path = Path(task.get("cache_prefix", "") + ".metadata.json")
        try:
            metadata = json.loads(metadata_path.read_text())
        except (OSError, json.JSONDecodeError) as error:
            problems.append(
                f"BENCHMARK_CACHE_METADATA_INVALID:{task.get('task_index')}:{error}")
            continue
        if metadata.get("workload_generation_id") != task.get(
                "cache_generation_id") or metadata.get(
                    "scientific_method_version") != \
                SCIENTIFIC_METHOD_VERSION or metadata.get(
                    "calibration_library") != CALIBRATION_LIBRARY or \
                metadata.get("source_scans") != {
                    "sites": 1, "observations": 1, "molecules": 1}:
            problems.append(
                f"BENCHMARK_CACHE_PROVENANCE_INVALID:{task.get('task_index')}")
    accounting_path = Path(root) / "benchmark_accounting.tsv"
    if accounting_path.is_file():
        accounting = list(read_tsv(accounting_path))
        if len(accounting) != 1 or accounting[0].get(
                "workload_generation_id") != generation or accounting[0].get(
                    "scientific_method_version") != \
                SCIENTIFIC_METHOD_VERSION or accounting[0].get(
                    "calibration_library") != str(CALIBRATION_LIBRARY) or int(
                accounting[0].get("measurement_rows", -1)) != len(measurements):
            problems.append("BENCHMARK_ACCOUNTING_BINDING_INVALID")
        elif any(not math.isfinite(finite(accounting[0].get(field, ""))) or
                 finite(accounting[0].get(field, "")) <= 0 for field in (
                     "configured_usable_cluster_memory_gib",
                     "configured_usable_cluster_cpus",
                     "measured_safe_aggregate_io_mib_per_second",
                     "measured_maximum_task_io_mib_per_second")) or int(
                         accounting[0].get(
                             "measured_io_overlap_processes", 0) or 0) < 1:
            problems.append("BENCHMARK_CONCURRENCY_MEASUREMENTS_INVALID")
    else:
        problems.append("BENCHMARK_ACCOUNTING_MISSING")
    return sorted(set(problems))


def _targeted_archive_problems(archive_path, root, generation, benchmark):
    """Validate the complete compact-return contract and current-file binding."""
    root = Path(root)
    archive_path = Path(archive_path)
    required = {
        "analysis/task_artifact_inventory.tsv",
        "analysis/decoy_match_diagnostics.tsv",
        "manifests/control_pairs_final.tsv",
        "manifests/extraction_tasks.tsv", "manifests/analysis_tasks.tsv",
        "targeted_validation_compact_summary.tsv",
        "targeted_validation_compact_results.tsv.gz",
        "warnings_and_exclusions.tsv", "targeted_provenance.tsv",
        "operational_status.json", "workload_blueprint.json",
        "exact_commands.json", "workload_accounting.tsv",
        "resource_projection.json", "frozen_targets_20260920.tsv", "frozen_target_comparisons.tsv",
    }
    if benchmark:
        required.update({
            "manifests/benchmark_panel.tsv", "benchmark_measurements.tsv",
            "benchmark_accounting.tsv",
            "manifests/benchmark_branch_coverage.tsv",
            "manifests/benchmark_completed_score_reference.tsv.gz",
            "manifests/lib12.rna_targeted_manifest.tsv.gz",
            "manifests/lib12.atac_targeted_manifest.tsv.gz",
            "BENCHMARK_RENDERED_UNSUBMITTED",
        })
    else:
        required.update({
            "target_decision_evidence.tsv",
            "control_performance_cluster_summary.tsv",
            "control_capacity_and_reuse.tsv", "target_comparison_balance.tsv",
            "target_cells_and_matched_comparisons.tsv",
            "primary_calibration_audit.tsv.gz",
            "primary_calibration_reference_roster.tsv.gz",
            "frozen_targets_20260920.tsv", "workload_accounting.tsv",
            "exact_commands.json",
            "TARGETED_WORKLOAD_RENDERED_UNSUBMITTED",
        })
        required.update({
            f"manifests/lib{library}.{modality}_targeted_manifest.tsv.gz"
            for library in ALLOWED_LIBRARIES for modality in ("rna", "atac")
        })
    problems = []
    try:
        with zipfile.ZipFile(archive_path) as archive:
            corrupt = archive.testzip()
            names = archive.namelist()
            if len(names) != len(set(names)):
                problems.append("RETURN_ARCHIVE_DUPLICATE_MEMBER")
            for info in archive.infolist():
                if info.filename.startswith("/") or ".." in Path(
                        info.filename).parts:
                    problems.append(
                        f"RETURN_ARCHIVE_UNSAFE_MEMBER:{info.filename}")
                if ((info.external_attr >> 16) & 0o170000) == 0o120000:
                    problems.append(
                        f"RETURN_ARCHIVE_SYMLINK_MEMBER:{info.filename}")
            if corrupt:
                problems.append(f"RETURN_ARCHIVE_CORRUPT:{corrupt}")
            missing = sorted(required - set(names))
            if missing:
                problems.append("RETURN_ARCHIVE_MISSING:" + ",".join(missing))
            allowed_exact = set(required) | {
                "control_parent_eligibility.tsv.gz",
            }
            allowed_prefixes = ("manifests/", "logs/", "markers/")
            if benchmark:
                allowed_prefixes += ("analysis/",)
            for name in names:
                if name not in allowed_exact and not name.startswith(
                        allowed_prefixes):
                    problems.append(f"RETURN_ARCHIVE_UNEXPECTED_ALIAS:{name}")
            # Every packaged regular file must be byte-identical to the current
            # artifact, preventing stale mixed-generation ZIP acceptance.
            for name in names:
                live = root / name
                try:
                    resolved = live.resolve(strict=True)
                    resolved.relative_to(root.resolve())
                except (OSError, ValueError):
                    problems.append(f"RETURN_ARCHIVE_UNBOUND_MEMBER:{name}")
                    continue
                if not resolved.is_file() or live.is_symlink():
                    problems.append(f"RETURN_ARCHIVE_UNBOUND_MEMBER:{name}")
                    continue
                archived_digest = hashlib.sha256(
                    archive.read(name)).hexdigest()
                if archived_digest != _small_file_sha256(resolved):
                    problems.append(f"RETURN_ARCHIVE_STALE_ENTRY:{name}")
    except (OSError, zipfile.BadZipFile) as error:
        return [f"RETURN_ARCHIVE_INVALID:{error}"]
    frozen = root / "frozen_targets_20260920.tsv"
    if not frozen.is_file() and benchmark:
        frozen = root.parent / "frozen_targets_20260920.tsv"
    try:
        _path, frozen_rows = load_frozen_targets(frozen)
        if Counter(int(row["library"]) for row in frozen_rows) != Counter(
                {12: 32, 20: 26, 29: 1}):
            problems.append("RETURN_ARCHIVE_FROZEN_TARGETS_INVALID")
    except (OSError, RuntimeError, ValueError) as error:
        problems.append(f"RETURN_ARCHIVE_FROZEN_TARGETS_INVALID:{error}")
    status_path = root / "operational_status.json"
    try:
        status = json.loads(status_path.read_text())
        if status.get("workload_generation_id") != generation or \
                status.get("operational_status") != "EXECUTION_COMPLETE" or \
                status.get("scientific_method_version") != \
                SCIENTIFIC_METHOD_VERSION or status.get(
                    "calibration_library") != CALIBRATION_LIBRARY:
            problems.append("RETURN_ARCHIVE_OPERATIONAL_STATUS_INVALID")
    except (OSError, json.JSONDecodeError) as error:
        problems.append(f"RETURN_ARCHIVE_OPERATIONAL_STATUS_INVALID:{error}")
    if not benchmark:
        decision_rows = list(read_tsv(root / "target_decision_evidence.tsv"))
        mappings = list(read_tsv(
            root / "target_cells_and_matched_comparisons.tsv"))
        expected_decisions = len(mappings) + sum(bool(clean(row.get(
            "matched_comparison_barcode", ""))) for row in mappings)
        required_decision = {
            "library", "barcode", "locked_source",
            "legacy_proposed_contributor", "candidate_relationship",
            "locked_copy_vector", "candidate_copy_vector",
            "complete_menu_signature", "complete_menu_size",
            "rna_site_winner", "atac_site_winner", "rna_molecule_winner",
            "atac_molecule_winner", "rna_site_null_observed_statistic",
            "rna_site_null_median", "rna_site_null_maximum",
            "rna_site_null_p95", "rna_site_null_p99",
            "rna_site_fit_interior_eligible",
            "rna_site_primary_calibration_available",
            "rna_site_primary_winner_p95",
            "rna_site_reference_cdf_margin_pass",
            "rna_site_supported_interior_mixture",
            "joint_four_channel_supported", "joint_four_channel_state",
            "plain_language_evidence_category", "unresolved_conflicts",
        }
        if len(decision_rows) != expected_decisions or not decision_rows or \
                not required_decision.issubset(decision_rows[0]):
            problems.append("RETURN_ARCHIVE_DECISION_TABLE_INVALID")
    return sorted(set(problems))


def _targeted_completion_problems(root, benchmark, deep=True):
    """One completion validator shared by status, resume, and launch gates."""
    root = Path(root).resolve()
    problems = []
    blueprint_path = root / "workload_blueprint.json"
    try:
        blueprint = json.loads(blueprint_path.read_text())
        generation = clean(blueprint.get("workload_generation_id", ""))
        if blueprint.get("scientific_method_version") != \
                SCIENTIFIC_METHOD_VERSION or int(blueprint.get(
                    "calibration_library", -1)) != CALIBRATION_LIBRARY:
            problems.append("COMPLETION_BLUEPRINT_METHOD_OR_CALIBRATION_MISMATCH")
    except (OSError, json.JSONDecodeError, ValueError) as error:
        return [f"COMPLETION_BLUEPRINT_INVALID:{error}"]
    extraction_path = root / "manifests" / "extraction_tasks.tsv"
    analysis_path = root / "manifests" / "analysis_tasks.tsv"
    try:
        extraction = list(read_tsv(extraction_path))
        analysis = list(read_tsv(analysis_path))
        if not extraction or not analysis:
            problems.append("COMPLETION_TASK_MANIFEST_EMPTY")
        else:
            validate_targeted_input_contract(extraction_path, "extraction", 0)
            validate_targeted_input_contract(analysis_path, "analysis", 0)
        for task in extraction:
            valid, task_problems, _ = validate_task_artifact(
                task, "extraction", True, deep=deep)
            if not valid:
                problems.extend(
                    f"COMPLETION_EXTRACTION_{task.get('task_index')}:" + item
                    for item in task_problems)
        for task in analysis:
            valid, task_problems, _ = validate_task_artifact(
                task, "analysis", True, deep=deep)
            if not valid:
                problems.extend(
                    f"COMPLETION_ANALYSIS_{task.get('task_index')}:" + item
                    for item in task_problems)
        finalizer_valid, finalizer_problems = _validate_control_finalizer(
            root, generation, analysis)
        if not finalizer_valid:
            problems.extend(finalizer_problems)
    except (OSError, RuntimeError, ValueError, csv.Error) as error:
        problems.append(f"COMPLETION_TASK_GRAPH_INVALID:{error}")
    terminal_path = root / "TARGETED_VALIDATION_FINISHED"
    complete_path = root / "TARGETED_VALIDATION_COMPLETE"
    try:
        terminal_bytes = terminal_path.read_bytes()
        complete_bytes = complete_path.read_bytes()
        if terminal_bytes != complete_bytes:
            problems.append("COMPLETION_TERMINAL_MARKERS_DIFFER")
        terminal = json.loads(terminal_bytes)
        if terminal.get("operational_status") != "EXECUTION_COMPLETE" or \
                terminal.get("workload_generation_id") != generation or \
                terminal.get("scientific_method_version") != \
                SCIENTIFIC_METHOD_VERSION or terminal.get(
                    "calibration_library") != CALIBRATION_LIBRARY:
            problems.append("COMPLETION_TERMINAL_PROVENANCE_INVALID")
        if not deep:
            outputs = [Path(terminal.get(field, "")) for field in ("combined_output", "return_archive")]
            if terminal.get("audited_terminal_files") != _file_states(outputs):
                problems.append("TERMINAL_OUTPUT_CHANGED_OR_UNAUDITED")
        else:
            combined = Path(terminal.get("combined_output", ""))
            if not combined.is_file() or terminal.get(
                    "combined_content_sha256") != _small_file_sha256(combined):
                problems.append("COMPLETION_COMBINED_OUTPUT_INVALID")
            archive = Path(terminal.get("return_archive", ""))
            if not archive.is_file() or terminal.get(
                    "return_archive_content_sha256") != _small_file_sha256(archive):
                problems.append("COMPLETION_RETURN_ARCHIVE_DIGEST_INVALID")
            else:
                problems.extend(_targeted_archive_problems(
                    archive, root, generation, benchmark))
    except (OSError, json.JSONDecodeError, TypeError) as error:
        problems.append(f"COMPLETION_TERMINAL_INVALID:{error}")
    return sorted(set(problems))


def targeted_compact_summary(rows):
    summaries = []
    by_group = defaultdict(list)
    for row in rows:
        by_group[(row.get("result_class", ""), row.get("role", ""),
                  row.get("modality", ""), row.get("evidence_channel", ""),
                  row.get("fraction", ""),
                  row.get("control_class", ""))].append(row)

    def wilson(successes, denominator):
        if not denominator:
            return math.nan, math.nan
        z = 1.959963984540054
        proportion = successes / denominator
        scale = 1 + z * z / denominator
        center = (proportion + z * z / (2 * denominator)) / scale
        radius = z * math.sqrt(
            proportion * (1 - proportion) / denominator +
            z * z / (4 * denominator * denominator)) / scale
        return center - radius, center + radius

    for key, values in sorted(by_group.items()):
        result_class, role, modality, channel, fraction, control_class = key
        if result_class in {
                "FULL_MENU_DOWNSAMPLING_REPLICATE",
                "CELL_CONDITIONAL_NULL_REPLICATE",
                "SOURCE_DISJOINT_CONTROL_CANDIDATE"}:
            continue
        row = {
            "result_class": result_class, "role": role,
            "modality": modality, "evidence_channel": channel,
            "fraction": fraction, "control_class": control_class,
            "rows": len(values),
            "available_rows": sum(value.get("status") == "AVAILABLE"
                                  for value in values),
        }
        if result_class == "FULL_MENU_DOWNSAMPLING":
            retention = [finite(value.get("winner_retention_fraction", ""))
                         for value in values]
            category = [finite(value.get("category_retention_fraction", ""))
                        for value in values]
            model = [finite(value.get("model_retention_fraction", ""))
                     for value in values]
            row.update({
                "mean_winner_retention_fraction": mean(retention),
                "mean_category_retention_fraction": mean(category),
                "mean_model_retention_fraction": mean(model),
                "replicates_per_cell": DOWNSAMPLE_REPLICATES,
            })
        elif result_class == "CELL_CONDITIONAL_LOCKED_MODEL_NULL":
            probabilities = [finite(value.get(
                "empirical_upper_tail_probability", "")) for value in values]
            evaluable = [value for value in probabilities if math.isfinite(value)]
            positives = sum(value <= 0.05 for value in evaluable)
            low, high = wilson(positives, len(evaluable))
            row.update({
                "evaluable": len(evaluable), "upper_tail_p_le_0_05": positives,
                "fraction_p_le_0_05": positives / len(evaluable)
                    if evaluable else math.nan,
                "fraction_p_le_0_05_wilson_low": low,
                "fraction_p_le_0_05_wilson_high": high,
                "null_replicates_per_cell": CELL_NULL_REPLICATES,
            })
        elif result_class == "SOURCE_DISJOINT_CONTROL_SUMMARY":
            available = [value for value in values
                         if value.get("status") == "AVAILABLE"]
            fit_eligible = sum(truthy(value.get(
                "fit_interior_eligible", "")) for value in available)
            decoy_wins = sum(truthy(value.get("decoy_won", ""))
                             for value in available)
            low = high = math.nan
            row.update({
                "evaluable": len(available),
                "fit_interior_eligible": fit_eligible,
                "fit_interior_eligible_fraction": fit_eligible / len(available)
                    if available else math.nan,
                "fit_interior_eligible_wilson_low": low,
                "uncertainty": "SEE_INDEPENDENT_PARENT_COMPONENT_SUMMARY",
                "fit_interior_eligible_wilson_high": high,
                "decoy_wins": decoy_wins,
                "calibrated_support_status":
                    "CALIBRATED_SUPPORT_NOT_EVALUATED",
                "interpretation": (
                    "per-channel fit diagnostics only; cross-channel control "
                    "decisions are in control_performance_cluster_summary.tsv"),
            })
        elif result_class == "OBSERVED_FULL_MENU_CANDIDATE":
            optimized = [value for value in values
                         if value.get("engine") == "OPTIMIZED_BATCHED"]
            winners = [value for value in optimized
                       if truthy(value.get("winner", ""))]
            row.update({
                "legal_candidate_rows": len(optimized),
                "winning_rows": len(winners),
                "fit_interior_eligible_winners": sum(
                    value.get("evidence_category") ==
                    "FIT_INTERIOR_ELIGIBLE"
                    for value in winners),
                "weak_addition_compatible_winners": sum(
                    value.get("evidence_category") ==
                    "WEAK_ADDITION_COMPATIBLE"
                    for value in winners),
                "calibrated_support_status":
                    "SEE_TARGET_DECISION_EVIDENCE_FOR_OBSERVED_TARGETS",
                "replacement_or_contributor_boundary_winners": sum(
                    value.get("evidence_category") ==
                    "REPLACEMENT_OR_CONTRIBUTOR_ONLY_BOUNDARY"
                    for value in winners),
            })
        summaries.append(row)
    return summaries


def _cluster_rate_interval(outcomes, clusters, seed_parts, replicates=2000):
    """Equal-cluster bootstrap interval plus sufficient cluster counts."""
    grouped = defaultdict(list)
    for outcome, cluster in zip(outcomes, clusters):
        if cluster and outcome is not None:
            grouped[cluster].append(bool(outcome))
    statistics_by_cluster = {
        key: {"successes": sum(values), "denominator": len(values)}
        for key, values in sorted(grouped.items())
    }
    cluster_rates = [item["successes"] / item["denominator"]
                     for item in statistics_by_cluster.values()
                     if item["denominator"]]
    if not cluster_rates:
        return math.nan, math.nan, math.nan, statistics_by_cluster
    estimate = mean(cluster_rates)
    if len(cluster_rates) < 2:
        return estimate, math.nan, math.nan, statistics_by_cluster
    generator = random.Random(stable_seed("cluster_bootstrap", *seed_parts))
    bootstrap = [mean(generator.choices(cluster_rates, k=len(cluster_rates)))
                 for _ in range(replicates)]
    return (estimate, quantile(bootstrap, 0.025), quantile(bootstrap, 0.975),
            statistics_by_cluster)


def targeted_decision_evidence(rows, root):
    """One human-auditable row per target/comparison with separate components."""
    root = Path(root)
    mapping_path = root / "target_cells_and_matched_comparisons.tsv"
    mapping_rows = list(read_tsv(mapping_path)) if mapping_path.is_file() else []
    manifest_by_cell = {}
    for item in mapping_rows:
        manifest_by_cell[("FROZEN_TARGET", item.get("library", ""),
                          item.get("target_barcode", ""))] = ("target", item)
        if clean(item.get("matched_comparison_barcode", "")):
            manifest_by_cell[("MATCHED_COMPARISON", item.get("library", ""),
                              item["matched_comparison_barcode"])] = \
                ("comparison", item)
    candidate_contract = {}
    for path in sorted((root / "manifests").glob(
            "lib*.*_targeted_manifest.tsv.gz")):
        modality = "RNA" if ".rna_" in path.name else "ATAC"
        for candidate in read_tsv(path):
            candidate_contract[(candidate.get("library", ""), modality,
                                candidate.get("barcode", ""),
                                candidate.get("candidate_id", ""))] = candidate
    calibration_roster_path = \
        root / "primary_calibration_reference_roster.tsv.gz"
    calibration_audit_path = root / "primary_calibration_audit.tsv.gz"
    calibration_rosters = defaultdict(list)
    if calibration_roster_path.is_file():
        for row in read_tsv(calibration_roster_path):
            if row.get("scheme") == "PRIMARY_EXACT_MASK_AND_COVERAGE":
                value = finite(row.get("reference_value", ""))
                if math.isfinite(value):
                    calibration_rosters[(row.get("target_id", ""),
                        row.get("modality", ""),
                        row.get("evidence_channel", ""))].append((
                            row.get("reference_cell_id", ""), value,
                            row.get("reference_stratum", ""),
                            row.get("coarsening", "")))
    calibration_audit = {}
    if calibration_audit_path.is_file():
        for row in read_tsv(calibration_audit_path):
            if row.get("scheme") == "PRIMARY_EXACT_MASK_AND_COVERAGE":
                calibration_audit[(row.get("target_id", ""),
                    row.get("modality", ""),
                    row.get("evidence_channel", ""))] = row
    balance_path = root / "target_comparison_balance.tsv"
    balance_rows = list(read_tsv(balance_path)) if balance_path.is_file() else []
    balance_by_library = {}
    for library in sorted({row.get("library", "") for row in balance_rows}):
        applicable = [row for row in balance_rows
                      if row.get("library") == library and
                      row.get("balance_stage") == "AFTER_MATCHING"]
        balance_by_library[library] = {
            "covariates": len(applicable),
            "all_pass": bool(applicable) and all(
                row.get("balance_status") == "PASS" for row in applicable),
            "maximum_absolute_standardized_difference": max((finite(
                row.get("absolute_standardized_difference", ""))
                for row in applicable), default=math.nan),
        }
    optimized = [row for row in rows
                 if row.get("engine", "OPTIMIZED_BATCHED") ==
                 "OPTIMIZED_BATCHED"]
    observed = defaultdict(list)
    nulls = {}
    downsampling = {}
    for row in optimized:
        cell_key = (row.get("role", ""), row.get("library", ""),
                    row.get("barcode", ""))
        channel_key = cell_key + (row.get("modality", ""),
                                  row.get("evidence_channel", ""))
        if row.get("result_class") == "OBSERVED_FULL_MENU_CANDIDATE":
            observed[channel_key].append(row)
        elif row.get("result_class") == "CELL_CONDITIONAL_LOCKED_MODEL_NULL":
            nulls[channel_key] = row
        elif row.get("result_class") == "FULL_MENU_DOWNSAMPLING":
            downsampling[channel_key + (row.get("fraction", ""),)] = row
    cell_keys = sorted(manifest_by_cell)
    output = []
    for role, library, barcode in cell_keys:
        member_name, manifest = manifest_by_cell.get(
            (role, library, barcode), ("unknown", {}))
        menu_signature_field = "complete_candidate_menu_signature" if \
            member_name == "target" else \
            "comparison_complete_candidate_menu_signature"
        menu_size_field = "complete_candidate_menu_size" if \
            member_name == "target" else \
            "comparison_complete_candidate_menu_size"
        result = {
            "library": library, "barcode": barcode, "role": role,
            "target_id": manifest.get("target_id", ""),
            "locked_source": manifest.get("locked_identity", "") if
                member_name == "target" else manifest.get(
                    "comparison_locked_identity", ""),
            "legacy_proposed_contributor": manifest.get(
                "proposed_contributor", "") if member_name == "target" else
                manifest.get("comparison_proposed_contributor", ""),
            "proposed_candidate_id": manifest.get("proposed_candidate_id", "")
                if member_name == "target" else manifest.get(
                    "comparison_proposed_candidate_id", ""),
            "candidate_physical_pool_state": manifest.get(
                "proposed_physical_pool_state", "") if member_name == "target"
                else manifest.get("comparison_physical_pool_state", ""),
            "candidate_component_only_state": manifest.get(
                "proposed_component_only_state", "") if member_name == "target"
                else manifest.get("comparison_component_only_state", ""),
            "candidate_new_donor_state": manifest.get(
                "proposed_new_donor_state", "") if member_name == "target"
                else manifest.get("comparison_new_donor_state", ""),
            "candidate_genotype_distinguishable": manifest.get(
                "proposed_genotype_distinguishable", "") if
                member_name == "target" else manifest.get(
                    "comparison_genotype_distinguishable", ""),
            "candidate_relationship": manifest.get(
                "proposed_relationship", "") if member_name == "target" else
                manifest.get("comparison_relationship", ""),
            "locked_copy_vector": manifest.get(
                "proposed_locked_copy_vector", "") if member_name == "target"
                else manifest.get("comparison_locked_copy_vector", ""),
            "candidate_copy_vector": manifest.get(
                "proposed_second_copy_vector", "") if member_name == "target"
                else manifest.get("comparison_second_copy_vector", ""),
            "candidate_origin_class": manifest.get(
                "proposed_candidate_origin_class", ""),
            "complete_menu_signature": manifest.get(
                menu_signature_field, ""),
            "complete_menu_size": manifest.get(menu_size_field, ""),
            "comparison_class": manifest.get("comparison_class", "")
                if role == "MATCHED_COMPARISON" else "NOT_APPLICABLE",
            "comparison_match_status": manifest.get("match_status", ""),
            "comparison_match_distance": manifest.get("match_distance", ""),
            "matching_stratum": manifest.get("matching_stratum", ""),
            "strict_site_five": manifest.get("strict_site_five", False)
                if member_name == "target" else False,
            "legacy_dual_molecule_p95_nine": manifest.get(
                "legacy_dual_molecule_p95_nine", False)
                if member_name == "target" else False,
            "legacy_exact_molecule_supported_seven": manifest.get(
                "legacy_exact_molecule_supported_seven", False)
                if member_name == "target" else False,
            "strict_exact_overlap_three": manifest.get(
                "strict_exact_overlap_three", False)
                if member_name == "target" else False,
            "interpretation_guardrail": (
                "component evidence for boss-chat review; this row is not an "
                "automatic physical-doublet call"),
        }
        winners = {}
        categories = {}
        models = {}
        for modality in ("RNA", "ATAC"):
            for channel in ("SITE", "MOLECULE"):
                prefix = f"{modality.lower()}_{channel.lower()}"
                candidates = observed.get(
                    (role, library, barcode, modality, channel), [])
                # Re-rank only the prespecified physical, new-donor,
                # genotype-distinguishable menu.  C++ deliberately emits all
                # legal candidates and a fit-only category; Python owns this
                # observed-target decision layer.
                ranked_available = []
                for item in candidates:
                    contract = candidate_contract.get((
                        library, modality, barcode,
                        item.get("candidate_id", "")), {})
                    eligible = (
                        truthy(contract.get("physical_pool_state", "")) and
                        not truthy(contract.get("component_only_state", "")) and
                        "CONTAINS_NEW_DONOR" in clean(contract.get(
                            "structural_added_state_relationship", "")) and
                        clean(contract.get("locked_copy_vector", "")) !=
                        clean(contract.get("second_copy_vector", "")) and
                        item.get("status") == "AVAILABLE" and
                        math.isfinite(finite(item.get(
                            "delta_log_likelihood", ""))))
                    if eligible:
                        ranked_available.append(item)
                ranked_available.sort(key=lambda item: (
                    finite(item.get("delta_log_likelihood", "")),
                    item.get("candidate_id", "")), reverse=True)
                winner = ranked_available[0] if ranked_available else {}
                runner = ranked_available[1] if len(ranked_available) > 1 else {}
                winners[(modality, channel)] = winner.get("second_state", "")
                categories[(modality, channel)] = winner.get(
                    "evidence_category", "UNAVAILABLE")
                models[(modality, channel)] = winner.get(
                    "preferred_model", "UNAVAILABLE")
                winner_delta = finite(winner.get("delta_log_likelihood", ""))
                runner_delta = finite(runner.get("delta_log_likelihood", ""))
                calibration_key = (manifest.get("target_id", ""), modality,
                                   channel)
                frozen_roster = sorted(calibration_rosters.get(
                    calibration_key, []), key=lambda item: (item[1], item[0])) \
                    if role == "FROZEN_TARGET" else []
                reference_values = [item[1] for item in frozen_roster]
                audit = calibration_audit.get(calibration_key, {}) \
                    if role == "FROZEN_TARGET" else {}
                roster_ids = [item[0] for item in frozen_roster]
                roster_keys = {item[2] for item in frozen_roster}
                roster_coarsening = {item[3] for item in frozen_roster}
                roster_key = next(iter(roster_keys), "")
                roster_coarsening_value = next(
                    iter(roster_coarsening), "")
                audit_winner = finite(audit.get("winner_delta", ""))
                audit_runner = finite(audit.get("runner_delta", ""))
                primary_available = bool(
                    math.isfinite(winner_delta) and
                    math.isfinite(runner_delta) and
                    len(reference_values) >= MIN_CALIBRATION_REFERENCE and
                    len(roster_ids) == len(set(roster_ids)) and
                    len(roster_keys) == 1 and len(roster_coarsening) == 1 and
                    audit.get("status") == "AVAILABLE" and
                    int(audit.get("reference_count", 0) or 0) ==
                    len(reference_values) and
                    audit.get("reference_stratum", "") == roster_key and
                    audit.get("coarsening", "") ==
                    roster_coarsening_value and
                    math.isclose(audit_winner, winner_delta,
                                 rel_tol=1e-12, abs_tol=1e-12) and
                    math.isclose(audit_runner, runner_delta,
                                 rel_tol=1e-12, abs_tol=1e-12))
                winner_count = bisect.bisect_right(
                    reference_values, winner_delta) \
                    if primary_available else None
                runner_count = bisect.bisect_right(
                    reference_values, runner_delta) \
                    if primary_available else None
                winner_cdf = winner_count / len(reference_values) \
                    if winner_count is not None else math.nan
                runner_cdf = runner_count / len(reference_values) \
                    if runner_count is not None else math.nan
                primary_p95 = bool(primary_available and
                    20 * winner_count >= 19 * len(reference_values))
                cdf_margin_pass = bool(primary_available and
                    100 * (winner_count - runner_count) >
                    len(reference_values))
                fit_eligible = bool(winner) and truthy(
                    winner.get("fit_interior_eligible", "")) and \
                    winner.get("evidence_category") == "FIT_INTERIOR_ELIGIBLE"
                supported = bool(fit_eligible and primary_available and
                                 primary_p95 and cdf_margin_pass)
                if supported:
                    final_category = "SUPPORTED_INTERIOR_MIXTURE"
                elif not winner:
                    final_category = "EVIDENCE_UNAVAILABLE"
                elif not fit_eligible:
                    final_category = winner.get(
                        "evidence_category", "EVIDENCE_UNAVAILABLE")
                elif not primary_available:
                    final_category = "PRIMARY_CALIBRATION_UNAVAILABLE"
                elif not primary_p95:
                    final_category = "INTERIOR_FIT_BELOW_PRIMARY_P95"
                else:
                    final_category = "INTERIOR_FIT_WEAK_CDF_SEPARATION"
                categories[(modality, channel)] = final_category
                low_statuses = sorted({item.get("status", "") for item in candidates
                                       if item.get("status", "") not in {
                                           "", "AVAILABLE"}})
                evidence_basis = "SITE_COUNTS" if channel == "SITE" else \
                    manifest.get(
                        f"{member_name}_{modality.lower()}_molecule_evidence_basis",
                        manifest.get(
                            f"target_{modality.lower()}_molecule_evidence_basis",
                            "UNAVAILABLE"))
                result.update({
                    f"{prefix}_evidence_basis": evidence_basis,
                    f"{prefix}_availability": winner.get("status",
                        ";".join(low_statuses) if low_statuses else
                        "UNAVAILABLE"),
                    f"{prefix}_legal_menu_candidates": winner.get(
                        "legal_menu_candidates", len(candidates)),
                    f"{prefix}_winner": winner.get("second_state", ""),
                    f"{prefix}_winner_candidate_id": winner.get(
                        "candidate_id", ""),
                    f"{prefix}_runner_up": runner.get("second_state", ""),
                    f"{prefix}_runner_up_candidate_id": runner.get(
                        "candidate_id", ""),
                    f"{prefix}_winner_delta": winner_delta,
                    f"{prefix}_runner_up_delta": runner_delta,
                    f"{prefix}_winner_margin": winner_delta - runner_delta
                        if math.isfinite(winner_delta) and
                        math.isfinite(runner_delta) else math.nan,
                    f"{prefix}_locked_log_likelihood": winner.get(
                        "locked_log_likelihood", ""),
                    f"{prefix}_interior_log_likelihood": winner.get(
                        "interior_log_likelihood", ""),
                    f"{prefix}_contributor_only_log_likelihood": winner.get(
                        "contributor_only_log_likelihood", ""),
                    f"{prefix}_fitted_fraction": winner.get(
                        "fitted_fraction", ""),
                    f"{prefix}_interval_low": winner.get(
                        "fitted_fraction_profile_low", ""),
                    f"{prefix}_interval_high": winner.get(
                        "fitted_fraction_profile_high", ""),
                    f"{prefix}_preferred_model": winner.get(
                        "preferred_model", "UNAVAILABLE"),
                    f"{prefix}_boundary_choice": winner.get(
                        "preferred_model", "UNAVAILABLE"),
                    f"{prefix}_fit_level_category": winner.get(
                        "evidence_category", "UNAVAILABLE"),
                    f"{prefix}_fit_interior_eligible": fit_eligible,
                    f"{prefix}_primary_calibration_available": primary_available,
                    f"{prefix}_primary_winner_p95": primary_p95,
                    f"{prefix}_reference_cdf_margin_pass": cdf_margin_pass,
                    f"{prefix}_supported_interior_mixture": supported,
                    f"{prefix}_evidence_category": final_category,
                    f"{prefix}_primary_reference_count": len(reference_values),
                    f"{prefix}_primary_winner_count": winner_count,
                    f"{prefix}_primary_runner_count": runner_count,
                    f"{prefix}_primary_winner_cdf": winner_cdf,
                    f"{prefix}_primary_runner_cdf": runner_cdf,
                    f"{prefix}_reference_cdf_margin": winner_cdf-runner_cdf
                        if primary_available else math.nan,
                    f"{prefix}_primary_reference_roster": json.dumps(roster_ids),
                    f"{prefix}_primary_reference_stratum":
                        roster_key,
                    f"{prefix}_primary_reference_coarsening":
                        roster_coarsening_value,
                    f"{prefix}_eligible_runner_set": json.dumps([
                        item.get("candidate_id", "")
                        for item in ranked_available[1:]]),
                    f"{prefix}_tie_rule": "DELTA_THEN_DESCENDING_CANDIDATE_ID",
                    f"{prefix}_evidence_units": winner.get(
                        "evidence_units", ""),
                    f"{prefix}_plain_language_component": {
                        "SUPPORTED_INTERIOR_MIXTURE":
                            "fit and primary calibration support an interior mixture",
                        "INTERIOR_FIT_BELOW_PRIMARY_P95":
                            "interior fit is below the primary calibration threshold",
                        "INTERIOR_FIT_WEAK_CDF_SEPARATION":
                            "interior fit does not separate from the runner-up",
                        "PRIMARY_CALIBRATION_UNAVAILABLE":
                            "fit is eligible but primary calibration is unavailable",
                        "WEAK_ADDITION_COMPATIBLE":
                            "addition-compatible fit fails one or more fit safeguards",
                        "REPLACEMENT_OR_CONTRIBUTOR_ONLY_BOUNDARY":
                            "replacement or contributor-only boundary",
                        "LOCKED_SOURCE_ONLY":
                            "locked source preferred",
                        "CONFLICTING_ASSAY_OR_EVIDENCE":
                            "conflicting assay or evidence result",
                        "LOW_OR_LIMITED_EVIDENCE":
                            "low or limited evidence",
                    }.get(final_category,
                          "evidence unavailable or conflicting"),
                })
                null = nulls.get(
                    (role, library, barcode, modality, channel), {})
                result.update({
                    f"{prefix}_null_status": null.get("status", "UNAVAILABLE"),
                    f"{prefix}_null_empirical_p": null.get(
                        "empirical_upper_tail_probability", ""),
                    f"{prefix}_null_observed_statistic": null.get(
                        "observed_maximum_delta", ""),
                    f"{prefix}_null_median": null.get("null_median", ""),
                    f"{prefix}_null_maximum": null.get("null_maximum", ""),
                    f"{prefix}_null_p95": null.get("null_p95", ""),
                    f"{prefix}_null_p99": null.get("null_p99", ""),
                    f"{prefix}_null_replicates": null.get(
                        "replicate_count", ""),
                })
                result[f"{prefix}_coverage_calibration_status"] = \
                    "AVAILABLE" if primary_available else "UNAVAILABLE"
                result[f"{prefix}_coverage_calibration_scheme"] = \
                    "PRIMARY_EXACT_MASK_AND_COVERAGE"
                for fraction in DOWNSAMPLE_FRACTIONS:
                    fraction_text = f"{fraction:.6g}"
                    downsample = downsampling.get(
                        (role, library, barcode, modality, channel,
                         fraction_text), {})
                    tag = str(fraction).replace(".", "_")
                    for field in ("winner_retention_fraction",
                                  "category_retention_fraction",
                                  "model_retention_fraction",
                                  "fitted_fraction_p025",
                                  "fitted_fraction_p975", "replicate_count",
                                  "successful_replicates",
                                  "unavailable_replicates", "winner_counts",
                                  "category_counts", "model_counts", "status"):
                        result[f"{prefix}_downsample_{tag}_{field}"] = \
                            downsample.get(field, "EVIDENCE_UNAVAILABLE")
        all_balance = balance_by_library.get("ALL_TARGET", {})
        library_balance = balance_by_library.get(library, {})
        matched = manifest.get("match_status") == "MATCHED"
        result.update({
            "site_rna_atac_contributor_consistent": bool(
                winners.get(("RNA", "SITE"))) and
                winners.get(("RNA", "SITE")) ==
                winners.get(("ATAC", "SITE")),
            "molecule_rna_atac_contributor_consistent": bool(
                winners.get(("RNA", "MOLECULE"))) and
                winners.get(("RNA", "MOLECULE")) ==
                winners.get(("ATAC", "MOLECULE")),
            "rna_site_molecule_contributor_consistent": bool(
                winners.get(("RNA", "SITE"))) and
                winners.get(("RNA", "SITE")) ==
                winners.get(("RNA", "MOLECULE")),
            "atac_site_molecule_contributor_consistent": bool(
                winners.get(("ATAC", "SITE"))) and
                winners.get(("ATAC", "SITE")) ==
                winners.get(("ATAC", "MOLECULE")),
            "coverage_conditional_status": (
                "AVAILABLE_MATCHED_AND_BALANCED" if member_name == "target" and
                matched and library_balance.get("all_pass") and
                all_balance.get("all_pass") else
                "UNRESOLVED_IMBALANCE_OR_NO_MATCH" if member_name == "target"
                else "NOT_APPLICABLE_COMPARISON_ROW"),
            "coverage_balance_library_all_covariates_pass":
                library_balance.get("all_pass", False),
            "coverage_balance_overall_all_covariates_pass":
                all_balance.get("all_pass", False),
            "coverage_balance_library_max_abs_standardized_difference":
                library_balance.get(
                    "maximum_absolute_standardized_difference", math.nan),
        })
        conflicts = []
        for channel in ("SITE", "MOLECULE"):
            rna = winners.get(("RNA", channel), "")
            atac = winners.get(("ATAC", channel), "")
            if rna and atac and rna != atac:
                conflicts.append(f"RNA_ATAC_{channel}_CONTRIBUTOR_CONFLICT")
        for modality in ("RNA", "ATAC"):
            site = winners.get((modality, "SITE"), "")
            molecule = winners.get((modality, "MOLECULE"), "")
            if site and molecule and site != molecule:
                conflicts.append(f"{modality}_SITE_MOLECULE_CONTRIBUTOR_CONFLICT")
        if any(value in {"LOW_OR_LIMITED_EVIDENCE", "UNAVAILABLE"}
               for value in categories.values()):
            conflicts.append("LOW_LIMITED_OR_UNAVAILABLE_EVIDENCE")
        nonempty_winners = {value for value in winners.values() if value}
        nonempty_models = {value for value in models.values()
                           if value not in {"", "UNAVAILABLE"}}
        all_four_supported = len(categories) == 4 and all(
            value == "SUPPORTED_INTERIOR_MIXTURE"
            for value in categories.values())
        contributor_concordant = len(nonempty_winners) == 1 and \
            len([value for value in winners.values() if value]) == 4
        model_concordant = len(nonempty_models) == 1 and \
            len([value for value in models.values()
                 if value not in {"", "UNAVAILABLE"}]) == 4
        if len(nonempty_winners) > 1:
            conflicts.append("FOUR_CHANNEL_CONTRIBUTOR_CONFLICT")
        if len(nonempty_models) > 1:
            conflicts.append("FOUR_CHANNEL_MODEL_CONFLICT")
        joint_supported = bool(all_four_supported and contributor_concordant and
                               model_concordant and not conflicts)
        if joint_supported:
            joint_state = "JOINT_FOUR_CHANNEL_SUPPORTED_INTERIOR_MIXTURE"
        elif len(nonempty_winners) > 1 or len(nonempty_models) > 1:
            joint_state = "JOINT_CONTRIBUTOR_OR_MODEL_CONFLICT"
        elif any(value == "PRIMARY_CALIBRATION_UNAVAILABLE"
                 for value in categories.values()):
            joint_state = "JOINT_PRIMARY_CALIBRATION_UNAVAILABLE"
        elif any(value in {"EVIDENCE_UNAVAILABLE", "UNAVAILABLE",
                           "LOW_OR_LIMITED_EVIDENCE"}
                 for value in categories.values()):
            joint_state = "JOINT_EVIDENCE_UNAVAILABLE"
        else:
            joint_state = "JOINT_NOT_SUPPORTED"
        result.update({
            "joint_four_channel_supported": joint_supported,
            "joint_four_channel_state": joint_state,
            "joint_contributor_concordant": contributor_concordant,
            "joint_model_concordant": model_concordant,
            "joint_supported_contributor": next(iter(nonempty_winners), "")
                if joint_supported else "",
            "joint_supported_model": next(iter(nonempty_models), "")
                if joint_supported else "",
            "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
            "calibration_library": CALIBRATION_LIBRARY,
        })
        result["unresolved_conflicts"] = ";".join(conflicts) or "NONE"
        result["plain_language_evidence_category"] = (
            "conflicting assay/evidence result" if conflicts else
            "supported by all four assay/evidence channels" if joint_supported else
            "fit is eligible but primary calibration is unavailable" if any(
                value == "PRIMARY_CALIBRATION_UNAVAILABLE"
                for value in categories.values()) else
            "interior fit is not sufficiently calibrated or separated" if any(
                value in {"INTERIOR_FIT_BELOW_PRIMARY_P95",
                          "INTERIOR_FIT_WEAK_CDF_SEPARATION",
                          "WEAK_ADDITION_COMPATIBLE"}
                for value in categories.values()) else
            "replacement or contributor-only boundary" if any(
                value == "REPLACEMENT_OR_CONTRIBUTOR_ONLY_BOUNDARY"
                for value in categories.values()) else
            "locked source preferred or evidence unavailable")
        output.append(result)
    return output


def targeted_control_performance(rows, root):
    """Control denominators and source-cluster-aware uncertainty summaries."""
    root = Path(root)
    final_path = root / "manifests" / "control_pairs_final.tsv"
    accounting_path = root / "control_capacity_and_reuse.tsv"
    controls = list(read_tsv(final_path)) if final_path.is_file() else []
    accounting = list(read_tsv(accounting_path)) \
        if accounting_path.is_file() else []
    by_id = {row["control_id"]: row for row in controls}
    parent = {}
    def component(key):
        parent.setdefault(key, key)
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key
    for control in controls:
        left = component((control["library"], control["recipient_barcode"]))
        right = component((control["library"], control["source_barcode"]))
        parent[max(left, right)] = min(left, right)
    for control in controls:
        control["parent_component"] = ":".join(component((control["library"], control["recipient_barcode"])))
    physical_new_states = defaultdict(set)
    for path in sorted((root / "manifests").glob(
            "lib*.*_targeted_manifest.tsv.gz")):
        modality = "RNA" if ".rna_" in path.name else "ATAC"
        for candidate in read_tsv(path):
            if truthy(candidate.get("physical_pool_state", "")) and \
                    not truthy(candidate.get("component_only_state", "")) and \
                    "CONTAINS_NEW_DONOR" in clean(candidate.get(
                        "structural_added_state_relationship", "")) and \
                    clean(candidate.get("locked_copy_vector", "")) != clean(
                        candidate.get("second_copy_vector", "")):
                physical_new_states[(candidate.get("library", ""), modality,
                    candidate.get("barcode", ""))].add(
                        candidate.get("second_state", ""))
    result_rows = [row for row in rows
                   if row.get("result_class") ==
                   "SOURCE_DISJOINT_CONTROL_SUMMARY" and
                   row.get("engine", "OPTIMIZED_BATCHED") ==
                   "OPTIMIZED_BATCHED"]
    # A control endpoint is a four-channel decision.  Synthetic controls have
    # no prespecified fraction-specific calibration CDF, so this combines only
    # the serialized C++ fit facts and explicitly does not claim calibrated
    # support.
    control_decisions = {}
    by_control_fraction = defaultdict(list)
    for row in result_rows:
        by_control_fraction[(row.get("control_id", ""),
                             row.get("fraction", ""))].append(row)
    for key, channel_rows in by_control_fraction.items():
        control = by_id.get(key[0], {})
        channel_map = {(row.get("modality", ""),
                        row.get("evidence_channel", "")): row
                       for row in channel_rows}
        expected_channels = {(modality, channel)
                             for modality in ("RNA", "ATAC")
                             for channel in ("SITE", "MOLECULE")}
        complete = set(channel_map) == expected_channels
        ordered = [channel_map[item] for item in sorted(expected_channels)] \
            if complete else []
        fit_eligible = complete and all(
            row.get("status") == "AVAILABLE" and
            truthy(row.get("fit_interior_eligible", "")) and
            row.get("evidence_category") == "FIT_INTERIOR_ELIGIBLE"
            for row in ordered)
        contributors = {clean(row.get("full_menu_winner", ""))
                        for row in ordered} - {"", "NA"}
        models = {clean(row.get("preferred_model", ""))
                  for row in ordered} - {"", "UNAVAILABLE", "NA"}
        contributor_concordant = complete and len(contributors) == 1
        model_concordant = complete and len(models) == 1
        winner = next(iter(contributors), "")
        expected = clean(control.get("expected_contributor", ""))
        control_class = control.get("control_class", "")
        recovery = bool(fit_eligible and contributor_concordant and
                        model_concordant and winner == expected)
        false_positive = bool(
            control_class == "SAME_SOURCE_NULL" and fit_eligible and
            contributor_concordant and model_concordant and winner and
            winner != clean(control.get("recipient_identity", "")) and
            all(winner in physical_new_states[(
                control.get("library", ""), modality,
                control.get("recipient_barcode", ""))]
                for modality in ("RNA", "ATAC")))
        control_decisions[key] = {
            "complete_four_channels": complete,
            "fit_eligible_all_four": fit_eligible,
            "contributor_concordant": contributor_concordant,
            "model_concordant": model_concordant,
            "winner": winner,
            "control_fit_recovery": recovery,
            "control_fit_false_positive": false_positive,
            "calibrated_support_status": "CALIBRATED_SUPPORT_NOT_EVALUATED",
        }
    observed = defaultdict(list)
    for row in result_rows:
        control = by_id.get(row.get("control_id", ""), {})
        recipient_basis = row.get("recipient_evidence_basis", "") or (
            "SITE_COUNTS" if row.get("evidence_channel") == "SITE" else
            control.get(
                f"{row.get('modality', '').lower()}_molecule_evidence_basis",
                "UNAVAILABLE"))
        source_basis = row.get("source_evidence_basis", "") or (
            "SITE_COUNTS" if row.get("evidence_channel") == "SITE" else
            control.get(
                f"source_{row.get('modality', '').lower()}_molecule_evidence_basis",
                "UNAVAILABLE"))
        basis = f"RECIPIENT={recipient_basis};SOURCE={source_basis}"
        row = dict(row)
        row["_basis_contract"] = basis
        observed[(row.get("library", ""), row.get("control_class", ""),
                  row.get("fraction", ""), row.get("modality", ""),
                  row.get("evidence_channel", ""))].append((row, control))
    accounting_by_key = {(row.get("library", ""), row.get("control_class", "")): row
                         for row in accounting if row.get("contributor", "ALL") == "ALL"}
    planned_keys = set(observed)
    for (library, control_class), _accounting in accounting_by_key.items():
        class_controls = [row for row in controls
                          if row.get("library") == library and
                          row.get("control_class") == control_class]
        for fraction in CONTROL_FRACTIONS:
            fraction_text = f"{fraction:.6g}"
            for modality in ("RNA", "ATAC"):
                for channel in ("SITE", "MOLECULE"):
                    planned_keys.add((library, control_class, fraction_text,
                                      modality, channel))
    output = []
    for key in sorted(planned_keys):
        library, control_class, fraction, modality, channel = key
        values = observed.get(key, [])
        available = [(row, control) for row, control in values
                     if row.get("status") == "AVAILABLE"]
        bases = sorted({row.get("_basis_contract", "") for row, _ in values})
        decisions = [(control_decisions.get((row.get("control_id", ""),
                      row.get("fraction", "")), {}), control)
                     for row, control in values]
        complete_decisions = [(decision, control)
                              for decision, control in decisions
                              if truthy(decision.get(
                                  "complete_four_channels", ""))]
        success = [truthy(decision.get("control_fit_recovery", ""))
                   for decision, _ in complete_decisions]
        false_positive = [truthy(decision.get(
            "control_fit_false_positive", ""))
            for decision, _ in complete_decisions] \
            if control_class == "SAME_SOURCE_NULL" else []
        clustered_outcomes = false_positive if control_class == \
            "SAME_SOURCE_NULL" else success
        cell_clusters = [control.get("parent_component", "")
                         for _, control in complete_decisions]
        donor_clusters = [control.get("source_donor_cluster", "")
                          for _, control in complete_decisions]
        cell_est, cell_low, cell_high, cell_stats = _cluster_rate_interval(
            clustered_outcomes, cell_clusters, key + ("SOURCE_CELL",))
        donor_est, donor_low, donor_high, donor_stats = _cluster_rate_interval(
            clustered_outcomes, donor_clusters, key + ("SOURCE_DONOR",))
        decoy_eligible = [(row, control) for row, control in available
                          if truthy(control.get("decoy_comparison_eligible", ""))]
        decoy_outcomes = [
            finite(row.get("expected_contributor_rank", "")) <
            finite(row.get("decoy_rank", ""))
            for row, _ in decoy_eligible
            if math.isfinite(finite(row.get("expected_contributor_rank", ""))) and
            math.isfinite(finite(row.get("decoy_rank", "")))]
        decoy_discriminated = sum(decoy_outcomes)
        decoy_controls = [control for row, control in decoy_eligible
                          if math.isfinite(finite(
                              row.get("expected_contributor_rank", ""))) and
                          math.isfinite(finite(row.get("decoy_rank", "")))]
        decoy_cell_est, decoy_cell_low, decoy_cell_high, decoy_cell_stats = \
            _cluster_rate_interval(
                decoy_outcomes,
                [control.get("parent_component", "")
                 for control in decoy_controls], key + ("DECOY_SOURCE_CELL",))
        plan = accounting_by_key.get((library, control_class), {})
        selected = int(plan.get("selected_pairs", 0) or 0)
        analyzed_ids = {row.get("control_id", "") for row, _ in available}
        result = {
            "library": library, "control_class": control_class,
            "fraction": fraction, "modality": modality,
            "evidence_channel": channel,
            "molecule_evidence_basis": ";".join(bases) or "UNAVAILABLE",
            "requested_pairs": int(plan.get("requested_pairs", 100) or 100),
            "eligible_pairs": int(plan.get("eligible_pairs", 0) or 0),
            "selected_pairs": selected,
            "successfully_analyzed_pairs": len(analyzed_ids),
            "excluded_or_failed_selected_pairs": max(0, selected-len(analyzed_ids)),
            "unique_sources": int(plan.get("unique_sources", 0) or 0),
            "unique_recipients": int(plan.get("unique_recipients", 0) or 0),
            "maximum_source_reuse": int(plan.get("maximum_source_reuse", 0) or 0),
            "source_reuse_distribution": plan.get(
                "source_reuse_distribution", "{}"),
            "available_outcomes": len(available),
            "expected_contributor_recovered": sum(success),
            "expected_recovery_fraction": sum(success)/len(success)
                if success else math.nan,
            "mean_realized_source_fraction": mean([
                finite(row.get("realized_source_fraction", ""))
                for row, _ in available]),
            "mean_planned_source_fraction": mean([
                finite(row.get("planned_source_fraction", ""))
                for row, _ in available]),
            "mean_recipient_units": mean([
                finite(row.get("recipient_units", ""))
                for row, _ in available]),
            "mean_source_units": mean([
                finite(row.get("source_units", ""))
                for row, _ in available]),
            "mean_successful_source_units": mean([
                finite(row.get("successful_source_units", ""))
                for row, _ in available]),
            "same_source_false_positives": sum(false_positive),
            "same_source_false_positive_denominator": len(false_positive),
            "same_source_false_positive_fraction":
                sum(false_positive)/len(false_positive)
                if false_positive else math.nan,
            "control_fit_recovery_definition": (
                "all four RNA/ATAC SITE/MOLECULE results are fit-eligible with "
                "one concordant planted contributor and model"),
            "control_fit_false_positive_definition": (
                "same-source null has four fit-eligible results with one "
                "concordant genetically distinct contributor and model"),
            "calibrated_support_status": "CALIBRATED_SUPPORT_NOT_EVALUATED",
            "dosage_shift_expected_recovery_fraction":
                sum(success)/len(success) if success and control_class ==
                "ALREADY_PRESENT_DONOR_DOSAGE_SHIFT" else math.nan,
            "decoy_eligible_pairs": len(decoy_eligible),
            "planted_source_ranked_ahead_of_decoy": decoy_discriminated,
            "planted_source_ahead_of_decoy_fraction":
                decoy_discriminated/len(decoy_outcomes)
                if decoy_outcomes else math.nan,
            "decoy_parent_component_equal_weight_rate": decoy_cell_est,
            "decoy_parent_component_bootstrap_p025": decoy_cell_low,
            "decoy_parent_component_bootstrap_p975": decoy_cell_high,
            "decoy_parent_component_sufficient_counts": json.dumps(
                decoy_cell_stats, sort_keys=True),
            "parent_component_count": len(cell_stats),
            "clustered_rate_metric": "SAME_SOURCE_FALSE_POSITIVE" if
                control_class == "SAME_SOURCE_NULL" else
                "EXPECTED_CONTRIBUTOR_RECOVERY",
            "parent_component_equal_weight_rate": cell_est,
            "parent_component_bootstrap_p025": cell_low,
            "parent_component_bootstrap_p975": cell_high,
            "parent_component_sufficient_counts": json.dumps(
                cell_stats, sort_keys=True),
            "source_donor_cluster_count": len(donor_stats),
            "source_donor_cluster_equal_weight_rate": donor_est,
            "source_donor_cluster_bootstrap_p025": math.nan,
            "source_donor_cluster_bootstrap_p975": math.nan,
            "source_donor_cluster_sufficient_counts": json.dumps(
                donor_stats, sort_keys=True),
            "cluster_bootstrap_replicates": 2000,
            "cluster_bootstrap_seed": SEED,
            "interpretation": (
                "fit-level control behavior only; calibrated support and a "
                "calibrated false-positive rate are not evaluated"),
        }
        class_controls = [row for row in controls
                          if row.get("library") == library and
                          row.get("control_class") == control_class]
        result["decoy_unmatched_selected_pairs"] = sum(
            control_class == "PLANTED_NEW_SOURCE" and
            not truthy(row.get("decoy_comparison_eligible", ""))
            for row in class_controls)
        for assay in ("rna", "atac"):
            for metric in ("distance_absolute_mismatch",
                           "distance_relative_mismatch",
                           "opportunity_absolute_mismatch",
                           "opportunity_relative_mismatch",
                           "genotype_site_count", "planted_distance",
                           "decoy_distance",
                           "continuous_coverage_mismatch"):
                metric_values = [finite(row.get(f"{assay}_{metric}", ""))
                                 for row in class_controls
                                 if truthy(row.get(
                                     "decoy_comparison_eligible", ""))]
                metric_values = [value for value in metric_values
                                 if math.isfinite(value)]
                result[f"{assay}_{metric}_median"] = median(metric_values)
                result[f"{assay}_{metric}_maximum"] = max(metric_values) \
                    if metric_values else math.nan
        output.append(result)
    return output


def targeted_gather_v4(args):
    root = Path(args.targeted_root).resolve()
    # A validated publication is immutable.  This also guarantees that a
    # manual/retried gather cannot damage a good archive or either bound
    # terminal marker if a later packaging attempt would fail.
    if (root / "TARGETED_VALIDATION_FINISHED").is_file() and \
            (root / "TARGETED_VALIDATION_COMPLETE").is_file():
        benchmark = (root / "BENCHMARK_RENDERED_UNSUBMITTED").is_file()
        prior_problems = _targeted_completion_problems(root, benchmark)
        if not prior_problems:
            print(json.dumps({
                "operational_status": "ALREADY_COMPLETE_VALIDATED",
                "workload_generation_id": args.generation,
                "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
                "calibration_library": CALIBRATION_LIBRARY,
            }, indent=2))
            return 0
    extraction_manifest_path = root / "manifests" / "extraction_tasks.tsv"
    analysis_manifest_path = root / "manifests" / "analysis_tasks.tsv"
    # Reuse the exact same full manifest boundary used by workers before any
    # cache/result artifact is opened.
    _task, extraction_generation = validate_targeted_input_contract(
        extraction_manifest_path, "extraction", 0)
    _task, analysis_generation = validate_targeted_input_contract(
        analysis_manifest_path, "analysis", 0)
    if extraction_generation != args.generation or \
            analysis_generation != args.generation:
        raise RuntimeError("gather generation does not match task manifests")
    extraction_tasks = list(read_tsv(
        extraction_manifest_path))
    validate_targeted_manifest_boundary(extraction_tasks, "extraction")
    tasks = list(read_tsv(analysis_manifest_path))
    validate_targeted_manifest_boundary(tasks, "analysis")
    valid_rows = []
    task_inventory = []
    invalid_extraction = []
    for task in extraction_tasks:
        valid, problems, details = validate_task_artifact(
            task, "extraction", require_marker=True)
        task_inventory.append({
            "task_kind": "extraction", "task_index": task["task_index"],
            "library": task["library"], "modality": task["modality"],
            "action": "NORMALIZED_EXTRACTION", "valid": valid,
            "problems": ";".join(problems), **details,
        })
        if not valid:
            invalid_extraction.append(int(task["task_index"]))
    invalid = []
    for task in tasks:
        valid, problems, details = validate_task_artifact(
            task, "analysis", require_marker=True)
        task_inventory.append({
            "task_kind": "analysis", "task_index": task["task_index"],
            "library": task["library"],
            "modality": task["modality"], "action": task["action"],
            "valid": valid, "problems": ";".join(problems), **details,
        })
        if valid:
            valid_rows.extend(read_tsv(task["analysis_output"]))
        else:
            invalid.append(int(task["task_index"]))
    finalizer_valid, control_finalizer_problems = _validate_control_finalizer(
        root, args.generation, tasks)
    if finalizer_valid:
        control_finalizer_problems = []
    write_tsv(root / "analysis" / "task_artifact_inventory.tsv", task_inventory)
    combined_problems = []
    combined_path = root / "targeted_validation_results.tsv.gz"
    if valid_rows:
        write_tsv(combined_path, valid_rows)
        try:
            header, combined_count, combined_digest = _full_gzip_audit(
                combined_path)
            required = {"schema_version", "workload_generation_id",
                        "scientific_method_version", "calibration_library",
                        "task_index", "result_key", "result_class",
                        "action", "library", "modality"}
            if not required.issubset(header):
                combined_problems.append("COMBINED_SCHEMA_MISSING_FIELDS")
            combined_keys = [row.get("result_key", "")
                             for row in read_tsv(combined_path)]
            if any(not key for key in combined_keys) or \
                    len(combined_keys) != len(set(combined_keys)):
                combined_problems.append("COMBINED_RESULT_KEYS_INVALID")
            if combined_count != len(valid_rows):
                combined_problems.append(
                    f"COMBINED_ROW_COUNT:{combined_count}!={len(valid_rows)}")
            if any(row.get("workload_generation_id") != args.generation
                   for row in valid_rows):
                combined_problems.append("COMBINED_GENERATION_MIXED_OR_STALE")
            if any(row.get("scientific_method_version") !=
                   SCIENTIFIC_METHOD_VERSION or
                   row.get("calibration_library") != "25"
                   for row in valid_rows):
                combined_problems.append(
                    "COMBINED_METHOD_OR_CALIBRATION_MIXED_OR_STALE")
            valid_task_indices = {str(task["task_index"]) for task in tasks}
            if any(row.get("task_index") not in valid_task_indices
                   for row in valid_rows):
                combined_problems.append("COMBINED_FOREIGN_TASK_INDEX")
        except (OSError, EOFError, gzip.BadGzipFile,
                UnicodeDecodeError, csv.Error) as error:
            combined_problems.append(f"COMBINED_OUTPUT_INVALID:{error}")
            combined_count = 0
            combined_digest = ""
    else:
        combined_count = 0
        combined_digest = ""
        combined_problems.append("NO_VALID_ROWS_TO_COMBINE")
    benchmark = root.name == "benchmark_unsubmitted"
    compact_summary = targeted_compact_summary(valid_rows)
    write_tsv(root / "targeted_validation_compact_summary.tsv",
              compact_summary)
    compact_rows = [row for row in valid_rows if row.get("result_class") not in {
        "FULL_MENU_DOWNSAMPLING_REPLICATE",
        "CELL_CONDITIONAL_NULL_REPLICATE",
    }]
    compact_results_path = root / "targeted_validation_compact_results.tsv.gz"
    write_tsv(compact_results_path, compact_rows)
    decision_path = root / "target_decision_evidence.tsv"
    control_summary_path = root / "control_performance_cluster_summary.tsv"
    decision_rows = targeted_decision_evidence(valid_rows, root) \
        if not benchmark else []
    control_rows = targeted_control_performance(valid_rows, root) \
        if not benchmark else []
    if not benchmark:
        write_tsv(decision_path, decision_rows)
        write_tsv(control_summary_path, control_rows)
        expected_cells = 2 * len({(row["library"], row[field])
            for row in read_tsv(root / "target_cells_and_matched_comparisons.tsv")
            for field in ("target_barcode", "matched_comparison_barcode")})
        if len(decision_rows) != expected_cells // 2:
            combined_problems.append(
                f"DECISION_CELL_COUNT:{len(decision_rows)}!={expected_cells // 2}")
        decision_keys = [(row.get("role"), row.get("library"), row.get("barcode"))
                         for row in decision_rows]
        decision_required = {
            "locked_source", "legacy_proposed_contributor",
            "candidate_physical_pool_state", "candidate_new_donor_state",
            "candidate_genotype_distinguishable", "candidate_relationship",
            "locked_copy_vector", "candidate_copy_vector",
            "complete_menu_signature", "complete_menu_size",
            "rna_site_winner", "atac_site_winner", "rna_molecule_winner",
            "atac_molecule_winner", "site_rna_atac_contributor_consistent",
            "molecule_rna_atac_contributor_consistent",
            "rna_site_null_observed_statistic", "rna_site_null_median",
            "rna_site_null_maximum", "rna_site_null_p95",
            "rna_site_null_p99", "rna_site_null_empirical_p",
            "rna_site_downsample_0_25_winner_retention_fraction",
            "atac_molecule_downsample_0_75_model_retention_fraction",
            "coverage_conditional_status", "strict_site_five",
            "rna_site_fit_interior_eligible",
            "rna_site_primary_calibration_available",
            "rna_site_primary_winner_p95",
            "rna_site_reference_cdf_margin_pass",
            "rna_site_supported_interior_mixture",
            "joint_four_channel_supported", "joint_four_channel_state",
            "joint_contributor_concordant", "joint_model_concordant",
            "scientific_method_version", "calibration_library",
            "legacy_dual_molecule_p95_nine",
            "legacy_exact_molecule_supported_seven",
            "strict_exact_overlap_three", "plain_language_evidence_category",
            "unresolved_conflicts", "interpretation_guardrail",
        }
        if decision_rows and not decision_required.issubset(decision_rows[0]):
            combined_problems.append("DECISION_EVIDENCE_SCHEMA_INVALID")
        for decision in decision_rows:
            reconstructed_support = []
            reconstructed_winners = []
            reconstructed_models = []
            for modality in ("rna", "atac"):
                for channel in ("site", "molecule"):
                    prefix = f"{modality}_{channel}"
                    available = truthy(decision.get(
                        f"{prefix}_primary_calibration_available", ""))
                    n_reference = int(decision.get(
                        f"{prefix}_primary_reference_count", 0) or 0)
                    winner_count = decision.get(
                        f"{prefix}_primary_winner_count", "")
                    runner_count = decision.get(
                        f"{prefix}_primary_runner_count", "")
                    try:
                        winner_count = int(winner_count)
                        runner_count = int(runner_count)
                    except (TypeError, ValueError):
                        winner_count = runner_count = -1
                    p95 = bool(available and n_reference >= 20 and
                               20 * winner_count >= 19 * n_reference)
                    margin = bool(available and n_reference >= 20 and
                                  100 * (winner_count-runner_count) >
                                  n_reference)
                    fit = truthy(decision.get(
                        f"{prefix}_fit_interior_eligible", ""))
                    supported = fit and available and p95 and margin
                    if p95 != truthy(decision.get(
                            f"{prefix}_primary_winner_p95", "")) or \
                            margin != truthy(decision.get(
                                f"{prefix}_reference_cdf_margin_pass", "")) or \
                            supported != truthy(decision.get(
                                f"{prefix}_supported_interior_mixture", "")):
                        combined_problems.append(
                            f"DECISION_BOOLEAN_RECONSTRUCTION_FAILED:"
                            f"{decision.get('library')}:{decision.get('barcode')}:"
                            f"{prefix}")
                    reconstructed_support.append(supported)
                    reconstructed_winners.append(decision.get(
                        f"{prefix}_winner", ""))
                    reconstructed_models.append(decision.get(
                        f"{prefix}_preferred_model", ""))
            joint = all(reconstructed_support) and \
                len(set(reconstructed_winners)) == 1 and \
                len(set(reconstructed_models)) == 1
            if joint != truthy(decision.get(
                    "joint_four_channel_supported", "")):
                combined_problems.append(
                    f"JOINT_DECISION_RECONSTRUCTION_FAILED:"
                    f"{decision.get('library')}:{decision.get('barcode')}")
        if len(decision_keys) != len(set(decision_keys)):
            combined_problems.append("DECISION_EVIDENCE_KEYS_NOT_UNIQUE")
        if any(task.get("action") == "CONTROL_BIN" for task in tasks) and \
                not control_rows:
            combined_problems.append("CONTROL_PERFORMANCE_SUMMARY_EMPTY")
        control_keys = [(row.get("library"), row.get("control_class"),
                         row.get("fraction"), row.get("modality"),
                         row.get("evidence_channel"),
                         row.get("molecule_evidence_basis"))
                        for row in control_rows]
        if len(control_keys) != len(set(control_keys)):
            combined_problems.append("CONTROL_PERFORMANCE_KEYS_NOT_UNIQUE")
        expected_control_keys = {
            (f"lib{library}", control_class, f"{fraction:.6g}", modality,
             channel)
            for library in ALLOWED_LIBRARIES
            for control_class in ("PLANTED_NEW_SOURCE", "SAME_SOURCE_NULL",
                                  "ALREADY_PRESENT_DONOR_DOSAGE_SHIFT")
            for fraction in CONTROL_FRACTIONS
            for modality in ("RNA", "ATAC")
            for channel in ("SITE", "MOLECULE")}
        observed_control_keys = {
            (row.get("library"), row.get("control_class"),
             row.get("fraction"), row.get("modality"),
             row.get("evidence_channel")) for row in control_rows}
        if len(control_rows) != len(expected_control_keys) or observed_control_keys != \
                expected_control_keys:
            combined_problems.append(
                "CONTROL_PERFORMANCE_EXPECTED_STRATA_MISMATCH")
    if benchmark:
        benchmark_blueprint = json.loads(
            (root / "workload_blueprint.json").read_text())
        benchmark_resources = _validated_targeted_resources(
            benchmark_blueprint.get("resource_policy", {}))
        extraction_threads = benchmark_resources["extraction_cpus"]
        production_threads = benchmark_resources["analysis_cpus"]
        measurements = []
        for task in extraction_tasks:
            path = root / "analysis" / f"extract_{task['task_index']}.time.txt"
            metrics = {
                "schema_version": "joint_doublet_benchmark_measurement_v4",
                "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
                "calibration_library": CALIBRATION_LIBRARY,
                "workload_generation_id": args.generation,
                "measurement_kind": "CACHE_EXTRACTION_RESOURCE_AND_ACCOUNTING",
                "stage": "EXTRACTION", "task_index": task["task_index"],
                "library": task["library"], "modality": task["modality"],
                "threads": extraction_threads,
                "projected_evidence_records": int(
                    task.get("projected_observation_records", 0) or 0) + int(
                    task.get("projected_molecule_records", 0) or 0),
            }
            if path.is_file():
                metrics.update(_parse_time_v(path))
            interval_path = root / "analysis" / \
                f"extract_{task['task_index']}.interval.tsv"
            if interval_path.is_file():
                metrics.update(_parse_process_interval(interval_path))
            metadata_path = Path(task["cache_prefix"] + ".metadata.json")
            if metadata_path.is_file():
                metadata = json.loads(metadata_path.read_text())
                cache_cells = list(read_tsv(
                    task["cache_prefix"] + ".cells.tsv"))
                menu_counts = Counter(
                    row.get("barcode", "")
                    for row in read_tsv(task["candidate_manifest"]))
                candidate_expanded_records = sum(
                    (int(row.get("observation_count", 0) or 0) +
                     int(row.get("molecule_count", 0) or 0)) *
                    menu_counts.get(row.get("barcode", ""), 0)
                    for row in cache_cells)
                unique_records = int(metadata.get(
                    "observation_records", 0)) + int(metadata.get(
                        "molecule_records", 0))
                metrics.update({
                    "measured_source_bytes": sum(int(metadata.get(field, 0))
                        for field in ("samples_source_bytes", "sites_source_bytes",
                                      "observations_source_bytes",
                                      "molecules_source_bytes")),
                    "measured_selected_cells": int(metadata.get(
                        "selected_cells", 0)),
                    "measured_candidate_rows": int(metadata.get(
                        "candidate_rows", 0)),
                    "measured_site_records": int(metadata.get(
                        "observation_records", 0)),
                    "measured_molecule_site_records": int(metadata.get(
                        "molecule_records", 0)),
                    "measured_linked_units": int(metadata.get(
                        "linked_units", 0)),
                    "measured_unique_evidence_records": unique_records,
                    "measured_candidate_expanded_records":
                        candidate_expanded_records,
                    "measured_candidate_expanded_rows_avoided":
                        candidate_expanded_records - unique_records,
                    "measured_normalization_reduction_factor":
                        candidate_expanded_records / unique_records
                        if unique_records else math.nan,
                    "measured_uncompressed_evidence_bytes":
                        int(metadata.get("observation_records", 0)) *
                        int(metadata.get("observation_record_bytes", 0)) +
                        int(metadata.get("molecule_records", 0)) *
                        int(metadata.get("molecule_record_bytes", 0)),
                    "measured_source_scans": sum(
                        int(value) for value in metadata.get(
                            "source_scans", {}).values()),
                    "written_cache_bytes": sum(
                        Path(task["cache_prefix"] + suffix).stat().st_size
                        for suffix in (".metadata.json", ".cells.tsv",
                                       ".observations.bin", ".molecules.bin",
                                       ".sites.bin", ".genotypes.tsv.gz",
                                       ".samples.tsv")
                        if Path(task["cache_prefix"] + suffix).is_file()),
                })
                metrics["logical_io_bytes"] = (
                    metrics["measured_source_bytes"] +
                    metrics["written_cache_bytes"])
                elapsed = finite(metrics.get("elapsed_seconds", ""))
                metrics["logical_io_mib_per_second"] = (
                    metrics["logical_io_bytes"] / (1 << 20) / elapsed
                    if math.isfinite(elapsed) and elapsed > 0 else math.nan)
            measurements.append(metrics)
        for task in tasks:
            for suffix, threads in (
                    ("threads1", 1),
                    (f"threads{production_threads}", production_threads)):
                path = Path(task["analysis_output"].replace(
                    ".tsv.gz", f".{suffix}.time.txt"))
                metrics = {
                           "schema_version":
                               "joint_doublet_benchmark_measurement_v4",
                           "scientific_method_version":
                               SCIENTIFIC_METHOD_VERSION,
                           "calibration_library": CALIBRATION_LIBRARY,
                           "workload_generation_id": args.generation,
                           "measurement_kind":
                               "THREAD_REFERENCE_EQUIVALENCE_AND_RESOURCE",
                           "stage": "ANALYSIS", "task_index": task["task_index"],
                           "library": task["library"],
                           "modality": task["modality"], "threads": threads,
                           "action": task["action"],
                           "role": task.get("role", ""),
                           "predicted_likelihood_work": task.get(
                               "predicted_likelihood_work", "")}
                if path.is_file():
                    metrics.update(_parse_time_v(path))
                interval_path = Path(task["analysis_output"].replace(
                    ".tsv.gz", f".{suffix}.interval.tsv"))
                if interval_path.is_file():
                    metrics.update(_parse_process_interval(interval_path))
                raw_output = Path(task["analysis_output"].replace(
                    ".tsv.gz", f".{suffix}.tsv.gz"))
                if raw_output.is_file():
                    metrics["output_bytes"] = raw_output.stat().st_size
                    metrics["logical_io_bytes"] = metrics["output_bytes"]
                    elapsed = finite(metrics.get("elapsed_seconds", ""))
                    metrics["logical_io_mib_per_second"] = (
                        metrics["logical_io_bytes"] / (1 << 20) / elapsed
                        if math.isfinite(elapsed) and elapsed > 0 else math.nan)
                    accounting_rows = [row for row in read_tsv(raw_output)
                                       if row.get("result_class") ==
                                       "COMPILED_SHARD_ACCOUNTING" and
                                       row.get("engine") == "OPTIMIZED_BATCHED"]
                    if accounting_rows:
                        for field in (
                                "optimized_fit_calls", "reference_fit_calls",
                                "optimizer_calls", "likelihood_evaluations",
                                "derivative_evaluations", "row_scans",
                                "candidate_owned_raw_evidence_copies",
                                "per_fraction_string_rehashes",
                                "per_fraction_string_sorts",
                                "unchanged_80_pass_refits"):
                            metrics[field] = int(
                                accounting_rows[0].get(field, 0) or 0)
                measurements.append(metrics)
        write_tsv(root / "benchmark_measurements.tsv", measurements)
        measured_elapsed = [finite(row.get("elapsed_seconds", ""))
                            for row in measurements
                            if math.isfinite(finite(row.get(
                                "elapsed_seconds", "")))]
        extraction_measurements = [row for row in measurements
                                   if row.get("stage") == "EXTRACTION"]
        analysis_measurements = [row for row in measurements
                                 if row.get("stage") == "ANALYSIS"]
        measured_unique = sum(int(row.get(
            "measured_unique_evidence_records", 0) or 0)
            for row in extraction_measurements)
        measured_expanded = sum(int(row.get(
            "measured_candidate_expanded_records", 0) or 0)
            for row in extraction_measurements)
        extraction_seconds = sum(finite(row.get("elapsed_seconds", ""), 0)
                                 for row in extraction_measurements)
        analysis_seconds = sum(finite(row.get("elapsed_seconds", ""), 0)
                               for row in analysis_measurements)
        (safe_aggregate_io_mib_per_second,
         maximum_task_io_mib_per_second,
         measured_io_overlap_processes) = _measured_io_capacity(
             extraction_measurements)
        write_tsv(root / "benchmark_accounting.tsv", [{
            "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
            "calibration_library": CALIBRATION_LIBRARY,
            "workload_generation_id": args.generation,
            "measurement_rows": len(measurements),
            "measured_source_bytes": sum(int(row.get(
                "measured_source_bytes", 0) or 0) for row in measurements),
            "measured_unique_evidence_records": sum(int(row.get(
                "measured_unique_evidence_records", 0) or 0)
                for row in measurements),
            "measured_site_records": sum(int(row.get(
                "measured_site_records", 0) or 0) for row in measurements),
            "measured_molecule_site_records": sum(int(row.get(
                "measured_molecule_site_records", 0) or 0)
                for row in measurements),
            "measured_linked_units": sum(int(row.get(
                "measured_linked_units", 0) or 0) for row in measurements),
            "measured_candidate_expanded_records": measured_expanded,
            "measured_candidate_expanded_rows_avoided":
                measured_expanded - measured_unique,
            "measured_normalization_reduction_factor":
                measured_expanded / measured_unique if measured_unique else math.nan,
            "measured_uncompressed_evidence_bytes": sum(int(row.get(
                "measured_uncompressed_evidence_bytes", 0) or 0)
                for row in measurements),
            "measured_written_cache_bytes": sum(int(row.get(
                "written_cache_bytes", 0) or 0) for row in measurements),
            "measured_analysis_output_bytes": sum(int(row.get(
                "output_bytes", 0) or 0) for row in measurements),
            "measured_optimizer_calls": sum(int(row.get(
                "optimizer_calls", 0) or 0) for row in measurements),
            "measured_likelihood_evaluations": sum(int(row.get(
                "likelihood_evaluations", 0) or 0) for row in measurements),
            "measured_task_seconds": sum(measured_elapsed),
            "measured_extraction_task_hours": extraction_seconds / 3600,
            "measured_extraction_core_hours":
                extraction_seconds * extraction_threads / 3600,
            "measured_analysis_task_hours": analysis_seconds / 3600,
            "measured_analysis_core_hours": sum(
                finite(row.get("elapsed_seconds", ""), 0) *
                int(row.get("threads", 0) or 0) / 3600
                for row in analysis_measurements),
            "extraction_task_count": len(extraction_tasks),
            "benchmark_analysis_shards": len(tasks),
            "timed_analysis_executions": len(analysis_measurements),
            "extraction_cpus_per_task": extraction_threads,
            "analysis_cpus_per_task": production_threads,
            "configured_usable_cluster_memory_gib": 512,
            "configured_usable_cluster_cpus": production_threads * 8,
            "measured_safe_aggregate_io_mib_per_second":
                safe_aggregate_io_mib_per_second,
            "measured_maximum_task_io_mib_per_second":
                maximum_task_io_mib_per_second,
            "measured_io_overlap_processes": measured_io_overlap_processes,
            "extraction_memory_limit": benchmark_resources[
                "extraction_memory"],
            "analysis_memory_limit": benchmark_resources[
                "analysis_memory"],
            "extraction_wall_limit": benchmark_resources[
                "extraction_wall_time"],
            "analysis_wall_limit": benchmark_resources[
                "analysis_wall_time"],
            "benchmark_array_concurrency": 4,
            "measured_peak_rss_kb": max(
                (int(row.get("max_rss_kb", 0) or 0) for row in measurements),
                default=0),
            "target_or_comparison_shards": sum(
                row.get("action") == "CELL" for row in tasks),
            "control_shards": sum(
                row.get("action") == "CONTROL_BIN" for row in tasks),
            "control_fractions": ",".join(map(str, CONTROL_FRACTIONS)),
            "channels": "SITE,MOLECULE",
            "downsample_replicates_per_fraction": DOWNSAMPLE_REPLICATES,
            "null_replicates_per_channel": CELL_NULL_REPLICATES,
            "measured_parallel_lower_bound_seconds": max(measured_elapsed)
                if measured_elapsed else math.nan,
            "measured_end_to_end_critical_path_lower_bound_seconds":
                max((finite(row.get("elapsed_seconds", ""), 0)
                     for row in extraction_measurements), default=0) +
                max((finite(row.get("elapsed_seconds", ""), 0)
                     for row in analysis_measurements), default=0),
            "queue_time": "EXCLUDED",
        }])
        combined_problems.extend(_benchmark_measurement_problems(
            root, args.generation, extraction_tasks, tasks, measurements))
    warning_rows = []
    for problem in sorted(set(invalid_extraction)):
        warning_rows.append({
            "scope": "TARGETED_EXTRACTION", "identifier": problem,
            "warning": "INVALID_EXTRACTION_TASK",
            "detail": "artifact failed shared validator",
        })
    for problem in sorted(set(invalid)):
        warning_rows.append({
            "scope": "TARGETED_ANALYSIS", "identifier": problem,
            "warning": "INVALID_ANALYSIS_TASK",
            "detail": "artifact failed shared validator",
        })
    for problem in sorted(set(control_finalizer_problems + combined_problems)):
        warning_rows.append({
            "scope": "TARGETED_GATHER", "identifier": args.generation,
            "warning": "TECHNICAL_VALIDATION_PROBLEM", "detail": problem,
        })
    final_control_path = root / "manifests" / "control_pairs_final.tsv"
    if final_control_path.is_file():
        for control in read_tsv(final_control_path):
            if control.get("control_class") == "PLANTED_NEW_SOURCE" and not \
                    truthy(control.get("decoy_comparison_eligible", "")):
                warning_rows.append({
                    "scope": "CONTROL_DECOY",
                    "identifier": control.get("control_id", ""),
                    "warning": "DECOY_EXCLUDED",
                    "detail": control.get("decoy_exclusion_reason", ""),
                })
    capacity_path = root / "control_capacity_and_reuse.tsv"
    if capacity_path.is_file():
        for row in read_tsv(capacity_path):
            if truthy(row.get("capacity_limited", "")):
                warning_rows.append({
                    "scope": "CONTROL_CAPACITY",
                    "identifier": f"{row.get('library')}:{row.get('control_class')}",
                    "warning": "CONTROL_CLASS_CAPACITY_LIMITED",
                    "detail": row.get("capacity_reason", ""),
                })
    provenance_rows = []
    for task in extraction_tasks + tasks:
        kind = "extraction" if "samples" in task else "analysis"
        record = {
            "schema_version": "joint_doublet_targeted_provenance_v4",
            "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
            "calibration_library": CALIBRATION_LIBRARY,
            "workload_generation_id": args.generation,
            "task_kind": kind, "task_index": task.get("task_index", ""),
            "library": task.get("library", ""),
            "modality": task.get("modality", ""),
            "action": task.get("action", "NORMALIZED_EXTRACTION"),
            "task_contract_sha256": hashlib.sha256(json.dumps(
                task, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
            "candidate_manifest": task.get("candidate_manifest", ""),
            "candidate_manifest_content_digest": task.get(
                "manifest_digest", _manifest_content_digest(
                    task["candidate_manifest"])),
            "artifact_marker": task.get("marker", ""),
        }
        if kind == "extraction":
            metadata_path = Path(task["cache_prefix"] + ".metadata.json")
            if metadata_path.is_file():
                metadata = json.loads(metadata_path.read_text())
                record.update({
                    "cache_generation_id": metadata.get(
                        "workload_generation_id", ""),
                    "samples_path": metadata.get("samples_path", ""),
                    "sites_path": metadata.get("sites_path", ""),
                    "observations_path": metadata.get("observations_path", ""),
                    "molecules_path": metadata.get("molecules_path", ""),
                    "samples_content_digest": metadata.get(
                        "samples_content_digest", ""),
                    "sites_content_digest": metadata.get(
                        "sites_content_digest", ""),
                    "observations_content_digest": metadata.get(
                        "observations_content_digest", ""),
                    "molecules_content_digest": metadata.get(
                        "molecules_content_digest", ""),
                    "observations_binary_content_digest": metadata.get(
                        "observations_binary_content_digest", ""),
                    "molecules_binary_content_digest": metadata.get(
                        "molecules_binary_content_digest", ""),
                    "sites_binary_content_digest": metadata.get(
                        "sites_binary_content_digest", ""),
                    "cells_index_content_digest": metadata.get(
                        "cells_index_content_digest", ""),
                    "samples_dictionary_content_digest": metadata.get(
                        "samples_dictionary_content_digest", ""),
                    "genotypes_dictionary_content_digest": metadata.get(
                        "genotypes_dictionary_content_digest", ""),
                })
        provenance_rows.append(record)
    write_tsv(root / "targeted_provenance.tsv", provenance_rows)
    write_tsv(root / "warnings_and_exclusions.tsv", warning_rows)

    operational = "EXECUTION_COMPLETE" \
        if tasks and not invalid and not invalid_extraction and \
        not control_finalizer_problems and not combined_problems else \
        "TECHNICAL_FAILURE"
    atomic_json(root / "operational_status.json", {
        "schema_version": "joint_doublet_targeted_operational_status_v4",
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "operational_status": operational,
        "scientific_status": "UNASSESSED" if benchmark else
            "COMPLETE_WITH_SCIENTIFIC_WARNINGS",
        "workload_generation_id": args.generation,
        "analysis_manifest_content_sha256": _small_file_sha256(
            root / "manifests" / "analysis_tasks.tsv"),
        "extraction_manifest_content_sha256": _small_file_sha256(
            root / "manifests" / "extraction_tasks.tsv"),
        "combined_content_sha256": combined_digest,
    })
    archive_paths = [
        root / "analysis" / "task_artifact_inventory.tsv",
        root / "analysis" / "decoy_match_diagnostics.tsv",
        root / "manifests",
        root / "targeted_validation_compact_summary.tsv",
        compact_results_path,
        root / "warnings_and_exclusions.tsv",
        root / "targeted_provenance.tsv",
        root / "operational_status.json",
        root / "workload_blueprint.json",
        root / "exact_commands.json",
        root / "workload_accounting.tsv",
        root / "resource_projection.json",
        root / ("BENCHMARK_RENDERED_UNSUBMITTED" if benchmark else
                "TARGETED_WORKLOAD_RENDERED_UNSUBMITTED"),
        root / "logs", root / "markers",
    ]
    if not benchmark:
        archive_paths.extend((decision_path, control_summary_path,
                              root / "control_capacity_and_reuse.tsv",
                              root / "target_comparison_balance.tsv",
                              root / "target_cells_and_matched_comparisons.tsv",
                              root / "primary_calibration_audit.tsv.gz",
                              root / "primary_calibration_reference_roster.tsv.gz",
                              root / "frozen_targets_20260920.tsv",
                              root / "frozen_target_comparisons.tsv",
                              root / "unavailable_cells.tsv",
                              root / "control_parent_eligibility.tsv.gz"))
    if benchmark:
        archive_paths.extend((root / "benchmark_measurements.tsv",
                              root / "benchmark_accounting.tsv"))
        for task in tasks:
            archive_paths.extend((
                Path(task["analysis_output"]),
                Path(str(task["analysis_output"]) +
                     ".field_differences.tsv"),
            ))
        archive_paths.append(root / "frozen_targets_20260920.tsv")
        archive_paths.append(root / "frozen_target_comparisons.tsv")
    final_archive = root / (
        "benchmark_equivalence_performance_return.zip" if benchmark else
        "targeted_validation_compact_return.zip")
    if operational != "EXECUTION_COMPLETE":
        failure = {
            "schema_version": "joint_doublet_targeted_terminal_v4",
            "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
            "calibration_library": CALIBRATION_LIBRARY,
            "operational_status": operational,
            "workload_generation_id": args.generation,
            "combined_problems": combined_problems,
            "invalid_analysis_tasks": invalid,
            "invalid_extraction_tasks": invalid_extraction,
            "control_finalizer_problems": control_finalizer_problems,
        }
        print(json.dumps(failure, indent=2))
        return 2
    candidate_name = f".{final_archive.name}.candidate.{os.getpid()}.zip"
    candidate_archive = compact_zip(
        root, [(path, root, "") for path in archive_paths], candidate_name)
    archive_problems = _targeted_archive_problems(
        candidate_archive, root, args.generation, benchmark)
    if archive_problems:
        combined_problems.extend(archive_problems)
        operational = "TECHNICAL_FAILURE"
        for problem in archive_problems:
            warning_rows.append({
                "scope": "RETURN_ARCHIVE", "identifier": args.generation,
                "warning": "RETURN_ARCHIVE_VALIDATION_FAILED",
                "detail": problem,
            })
        write_tsv(root / "warnings_and_exclusions.tsv", warning_rows)
        atomic_json(root / "operational_status.json", {
            "schema_version": "joint_doublet_targeted_operational_status_v4",
            "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
            "calibration_library": CALIBRATION_LIBRARY,
            "operational_status": operational,
            "scientific_status": "UNASSESSED",
            "workload_generation_id": args.generation,
            "archive_problems": archive_problems,
        })
        candidate_archive.unlink(missing_ok=True)
        print(json.dumps({
            "schema_version": "joint_doublet_targeted_terminal_v4",
            "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
            "calibration_library": CALIBRATION_LIBRARY,
            "operational_status": operational,
            "workload_generation_id": args.generation,
            "archive_problems": archive_problems,
        }, indent=2))
        return 2
    os.replace(candidate_archive, final_archive)
    return_archive = final_archive
    marker = {
        "schema_version": "joint_doublet_targeted_terminal_v4",
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "operational_status": operational,
        "scientific_status": "UNASSESSED" if benchmark else
            "COMPLETE_WITH_SCIENTIFIC_WARNINGS",
        "scientific_summary": (
            "benchmark equivalence/performance only" if benchmark else
            "targeted evidence complete; physical doublets require boss-chat review"),
        "utc": utc_now(), "workload_generation_id": args.generation,
        "required_tasks": len(tasks), "valid_tasks": len(tasks) - len(invalid),
        "invalid_analysis_tasks": invalid,
        "invalid_extraction_tasks": invalid_extraction,
        "control_finalizer_problems": control_finalizer_problems,
        "combined_output": str(combined_path),
        "combined_rows": combined_count,
        "combined_content_sha256": combined_digest,
        "combined_problems": combined_problems,
        "compact_summary_rows": len(compact_summary),
        "compact_results_rows": len(compact_rows),
        "target_decision_rows": len(decision_rows),
        "control_performance_summary_rows": len(control_rows),
        "return_archive": str(return_archive),
        "return_archive_content_sha256": _small_file_sha256(return_archive),
        "return_archive_problems": archive_problems,
        "analysis_manifest_content_sha256": _small_file_sha256(
            root / "manifests" / "analysis_tasks.tsv"),
        "extraction_manifest_content_sha256": _small_file_sha256(
            root / "manifests" / "extraction_tasks.tsv"),
        "accessed_libraries": sorted({
            int(row["library"].removeprefix("lib")) for row in tasks}),
    }
    marker["protected_libraries_accessed"] = sorted(
        set(marker["accessed_libraries"]) & set(PROTECTED_LIBRARIES))
    if benchmark and (root / "benchmark_measurements.tsv").is_file():
        marker["benchmark_measurements_content_sha256"] = _small_file_sha256(
            root / "benchmark_measurements.tsv")
    marker["audited_terminal_files"] = _file_states([combined_path, return_archive])
    atomic_json(root / "TARGETED_VALIDATION_FINISHED", marker)
    atomic_json(root / "TARGETED_VALIDATION_COMPLETE", marker)
    print(json.dumps(marker, indent=2))
    return 0 if operational == "EXECUTION_COMPLETE" else 2


def _promoted_library12_cache_problems(extraction_tasks):
    """Prove both promoted cache modalities match one exact selected/menu panel."""
    problems = []
    tasks = [task for task in extraction_tasks if task.get("library") == "lib12"]
    if len(tasks) != 2 or {task.get("modality") for task in tasks} != {
            "RNA", "ATAC"}:
        return ["PROMOTED_CACHE_REQUIRES_EXACTLY_LIB12_RNA_AND_ATAC"]
    expected_cells_by_modality = {}
    menu_keys_by_modality = {}
    provenance = {}
    for task in tasks:
        valid, task_problems, details = validate_task_artifact(
            task, "extraction", require_marker=True)
        if not valid:
            problems.extend(f"PROMOTED_CACHE_{task['modality']}:{problem}"
                            for problem in task_problems)
            continue
        candidate_rows = list(read_tsv(task["candidate_manifest"]))
        menu_keys = [(row.get("barcode", ""), row.get("candidate_id", ""),
                      row.get("locked_state", ""), row.get("locked_copy_vector", ""),
                      row.get("second_state", ""), row.get("second_copy_vector", ""),
                      row.get("physical_pool_state", ""),
                      row.get("component_only_state", ""),
                      row.get("structural_added_state_relationship", ""))
                     for row in candidate_rows]
        menu_keys_by_modality[task["modality"]] = menu_keys
        expected_cells = {row[0] for row in menu_keys}
        expected_cells_by_modality[task["modality"]] = expected_cells
        cells_path = Path(task["cache_prefix"] + ".cells.tsv")
        cell_rows = list(read_tsv(cells_path))
        if {row.get("barcode", "") for row in cell_rows} != expected_cells:
            problems.append(
                f"PROMOTED_CACHE_{task['modality']}_SELECTED_CELL_SET_MISMATCH")
        metadata_path = Path(task["cache_prefix"] + ".metadata.json")
        metadata = json.loads(metadata_path.read_text())
        if metadata.get("schema_version") != CACHE_SCHEMA or \
                metadata.get("workload_generation_id") != task.get(
                    "cache_generation_id") or metadata.get(
                    "manifest_digest") != task.get("manifest_digest"):
            problems.append(
                f"PROMOTED_CACHE_{task['modality']}_SCHEMA_GENERATION_MENU_MISMATCH")
        if int(metadata.get("selected_cells", -1)) != len(expected_cells) or \
                int(metadata.get("candidate_rows", -1)) != len(menu_keys):
            problems.append(
                f"PROMOTED_CACHE_{task['modality']}_COUNT_MISMATCH")
        if int(metadata.get("observation_records", 0) or 0) <= 0 or \
                int(metadata.get("molecule_records", 0) or 0) <= 0 or \
                int(metadata.get("linked_units", 0) or 0) <= 0:
            problems.append(
                f"PROMOTED_CACHE_{task['modality']}_EVIDENCE_CHANNEL_INCOMPLETE")
        components = [Path(task["cache_prefix"] + suffix) for suffix in (
            ".cells.tsv", ".observations.bin", ".molecules.bin", ".sites.bin",
            ".genotypes.tsv.gz", ".samples.tsv")]
        if metadata.get("atomic_publication") != "metadata published last" or \
                any(path.stat().st_mtime_ns > metadata_path.stat().st_mtime_ns
                    for path in components if path.is_file()):
            problems.append(
                f"PROMOTED_CACHE_{task['modality']}_ATOMIC_PUBLICATION_INVALID")
        provenance[task["modality"]] = tuple(metadata.get(field, "") for field in (
            "samples_content_digest", "sites_content_digest",
            "observations_content_digest", "molecules_content_digest"))
        marker = json.loads(Path(task["marker"]).read_text())
        if marker.get("task_configuration") != dict(task):
            problems.append(
                f"PROMOTED_CACHE_{task['modality']}_MARKER_CONTRACT_MISMATCH")
    if len(expected_cells_by_modality) == 2 and \
            expected_cells_by_modality["RNA"] != expected_cells_by_modality["ATAC"]:
        problems.append("PROMOTED_CACHE_MODALITY_SELECTED_CELL_SET_MISMATCH")
    if len(menu_keys_by_modality) == 2 and \
            menu_keys_by_modality["RNA"] != menu_keys_by_modality["ATAC"]:
        problems.append("PROMOTED_CACHE_MODALITY_COMPLETE_MENU_MISMATCH")
    if len(provenance) == 2 and any(not re.fullmatch(
            r"fnv1a64:[0-9a-f]{16}", value)
            for values in provenance.values() for value in values):
        problems.append("PROMOTED_CACHE_SOURCE_PROVENANCE_INVALID")
    return sorted(set(problems))


def targeted_reproject(args):
    top = Path(args.targeted_root).resolve()
    benchmark = top / "benchmark_unsubmitted"
    shared_completion_problems = _targeted_completion_problems(
        benchmark, True)
    if shared_completion_problems:
        raise RuntimeError(
            "benchmark completion validation failed: " +
            ";".join(shared_completion_problems))
    marker_path = benchmark / "TARGETED_VALIDATION_FINISHED"
    measurements_path = benchmark / "benchmark_measurements.tsv"
    if not marker_path.is_file():
        raise RuntimeError("validated benchmark completion is required for reprojection")
    extraction_manifest = benchmark / "manifests" / "extraction_tasks.tsv"
    analysis_manifest = benchmark / "manifests" / "analysis_tasks.tsv"
    # Validate the complete manifest graph before opening the terminal marker,
    # cache products, measurements, or any other generated artifact.
    validate_targeted_input_contract(extraction_manifest, "extraction", 0)
    validate_targeted_input_contract(analysis_manifest, "analysis", 0)
    benchmark_extraction = list(read_tsv(extraction_manifest))
    benchmark_analysis = list(read_tsv(analysis_manifest))
    validate_targeted_manifest_boundary(benchmark_extraction, "extraction")
    validate_targeted_manifest_boundary(benchmark_analysis, "analysis")
    benchmark_marker = json.loads(marker_path.read_text())
    benchmark_generation = benchmark_analysis[0]["workload_generation_id"] \
        if benchmark_analysis else ""
    benchmark_problems = []
    if benchmark_marker.get("operational_status") != "EXECUTION_COMPLETE" or \
            benchmark_marker.get("workload_generation_id") != benchmark_generation:
        benchmark_problems.append("BENCHMARK_TERMINAL_MARKER_MISMATCH")
    for task in benchmark_extraction:
        valid, problems, _ = validate_task_artifact(task, "extraction", True)
        if not valid:
            benchmark_problems.append(
                f"INVALID_BENCHMARK_EXTRACTION_{task['task_index']}:" +
                ";".join(problems))
    benchmark_problems.extend(
        _promoted_library12_cache_problems(benchmark_extraction))
    for task in benchmark_analysis:
        valid, problems, _ = validate_task_artifact(task, "analysis", True)
        if not valid:
            benchmark_problems.append(
                f"INVALID_BENCHMARK_ANALYSIS_{task['task_index']}:" +
                ";".join(problems))
    finalizer_valid, finalizer_problems = _validate_control_finalizer(
        benchmark, benchmark_generation, benchmark_analysis)
    if not finalizer_valid:
        benchmark_problems.extend(finalizer_problems)
    combined = Path(benchmark_marker.get("combined_output", "")) \
        if benchmark_marker.get("combined_output") else None
    if combined is None or not combined.is_file() or \
            benchmark_marker.get("combined_content_sha256") != \
            _small_file_sha256(combined):
        benchmark_problems.append("BENCHMARK_COMBINED_OUTPUT_INVALID")
    if benchmark_problems:
        raise RuntimeError("benchmark validation failed: " +
                           ";".join(benchmark_problems))
    measurements = list(read_tsv(measurements_path))
    benchmark_problems.extend(_benchmark_measurement_problems(
        benchmark, benchmark_generation, benchmark_extraction,
        benchmark_analysis, measurements))
    if benchmark_marker.get("benchmark_measurements_content_sha256") != \
            _small_file_sha256(measurements_path):
        benchmark_problems.append("BENCHMARK_MEASUREMENT_DIGEST_MISMATCH")
    if benchmark_problems:
        raise RuntimeError("benchmark validation failed: " +
                           ";".join(sorted(set(benchmark_problems))))
    measured = [row for row in measurements
                if math.isfinite(finite(row.get("elapsed_seconds", ""))) and
                math.isfinite(finite(row.get("max_rss_kb", "")))]
    benchmark_blueprint = json.loads(
        (benchmark / "workload_blueprint.json").read_text())
    benchmark_resources = _validated_targeted_resources(
        benchmark_blueprint.get("resource_policy", {}))
    benchmark_production_threads = benchmark_resources["analysis_cpus"]
    analysis_measured = [row for row in measured
                         if row.get("stage") == "ANALYSIS" and
                         int(row.get("threads", 0) or 0) ==
                         benchmark_production_threads and
                         finite(row.get("predicted_likelihood_work", ""), 0) > 0]
    extraction_measured = [row for row in measured
                           if row.get("stage") == "EXTRACTION" and
                           finite(row.get("measured_source_bytes", ""), 0) > 0]
    if len(analysis_measured) != len(benchmark_analysis) or \
            len(extraction_measured) != len(benchmark_extraction):
        raise RuntimeError(
            "benchmark measurement panel is incomplete after strict validation")
    blueprint = json.loads((top / "workload_blueprint.json").read_text())
    previous = blueprint["workload_generation_id"]
    measurement_digest = _small_file_sha256(measurements_path)
    generation = "reprojected_" + hashlib.sha256(
        f"{previous}:{measurement_digest}".encode()).hexdigest()[:16]
    generation_root = top / "generations" / generation
    for name in ("manifests", "slurm_scripts", "logs", "task_scratch",
                 "scores", "cache", "analysis", "markers"):
        (generation_root / name).mkdir(parents=True, exist_ok=True)
    analysis_peak_rss_gib = max(finite(row["max_rss_kb"])
                                for row in analysis_measured) / (1024 * 1024)
    extraction_peak_rss_gib = max(finite(row["max_rss_kb"])
                                  for row in extraction_measured) / (1024 * 1024)
    requested_memory_gib = max(8, int(math.ceil(analysis_peak_rss_gib * 1.5)))
    extraction_memory_gib = max(
        16, int(math.ceil(extraction_peak_rss_gib * 1.5)))
    seconds_per_work = max(
        finite(row["elapsed_seconds"]) /
        finite(row["predicted_likelihood_work"])
        for row in analysis_measured)
    extraction_seconds_per_byte = max(
        finite(row["elapsed_seconds"]) /
        finite(row["measured_source_bytes"])
        for row in extraction_measured)
    target_analysis_seconds = 4 * 3600
    target_bin_work = target_analysis_seconds / seconds_per_work
    benchmark_accounting_rows = list(read_tsv(
        benchmark / "benchmark_accounting.tsv"))
    if len(benchmark_accounting_rows) != 1:
        raise RuntimeError("benchmark accounting must contain exactly one row")
    benchmark_accounting = benchmark_accounting_rows[0]
    concurrency_memory_budget_gib = finite(benchmark_accounting.get(
        "configured_usable_cluster_memory_gib", ""))
    concurrency_cpu_budget = finite(benchmark_accounting.get(
        "configured_usable_cluster_cpus", ""))
    safe_aggregate_io = finite(benchmark_accounting.get(
        "measured_safe_aggregate_io_mib_per_second", ""))
    maximum_task_io = finite(benchmark_accounting.get(
        "measured_maximum_task_io_mib_per_second", ""))
    concurrency, concurrency_caps, limiting_resources = \
        _resource_concurrency_caps(
            concurrency_memory_budget_gib, requested_memory_gib,
            concurrency_cpu_budget, benchmark_production_threads,
            safe_aggregate_io, maximum_task_io, hard_cap=8)
    extraction = [dict(row) for row in blueprint["extraction_tasks"]]
    original_analysis = [dict(row) for row in blueprint["analysis_tasks"]]
    copied_manifests = {}
    for source_text in sorted({row["candidate_manifest"] for row in extraction}):
        source_path = Path(source_text)
        destination = generation_root / "manifests" / source_path.name
        shutil.copyfile(source_path, destination)
        copied_manifests[source_text] = str(destination)
    for row in extraction:
        row["workload_generation_id"] = generation
        row["candidate_manifest"] = copied_manifests[row["candidate_manifest"]]
        row["manifest_digest"] = _manifest_content_digest(
            row["candidate_manifest"])
        if row["library"] == "lib12":
            row["cache_policy"] = "IMMUTABLE_BENCHMARK_REUSE_REQUIRED"
            row["promoted_cache_source_generation"] = row[
                "cache_generation_id"]
            row["promoted_cache_proof_sha256"] = hashlib.sha256(json.dumps({
                "benchmark_generation": benchmark_generation,
                "benchmark_marker_sha256": _small_file_sha256(marker_path),
                "measurement_sha256": measurement_digest,
                "modality": row["modality"],
                "manifest_digest": row["manifest_digest"],
            }, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
            modality = row["modality"].lower()
            row["marker"] = str(generation_root / "markers" /
                f"extract_{row['library']}_{modality}.promoted.complete.json")
        else:
            row["cache_generation_id"] = generation
            modality = row["modality"].lower()
            row["cache_prefix"] = str(generation_root / "cache" /
                f"{row['library']}.{modality}.{generation}")
            row["marker"] = str(generation_root / "markers" /
                f"extract_{row['library']}_{modality}.complete.json")
            row["cache_policy"] = "EXTRACT_ONCE_FOR_GENERATION"
    extraction_by_key = {(row["library"], row["modality"]): row
                         for row in extraction}

    # Cell shards remain one cell/modality.  Control IDs are re-binned from
    # benchmark-calibrated predicted work while retaining the eight-pair cap.
    analysis = [dict(row) for row in original_analysis
                if row["action"] == "CELL"]
    control_items = defaultdict(list)
    control_prototypes = {}
    for row in original_analysis:
        if row["action"] != "CONTROL_BIN":
            continue
        key = (row["library"], row["modality"])
        control_prototypes[key] = row
        identifiers = [value for value in row.get("control_ids", "").split(",")
                       if value]
        per_pair = finite(row.get("predicted_likelihood_work", ""), 0) / \
            max(len(identifiers), 1)
        control_items[key].extend((per_pair, identifier)
                                  for identifier in identifiers)
    for key, items in sorted(control_items.items()):
        bins = []
        for work, identifier in sorted(items, key=lambda item: (-item[0], item[1])):
            choices = [item for item in bins
                       if len(item["ids"]) < 8 and
                       item["work"] + work <= target_bin_work]
            if choices:
                selected = min(choices, key=lambda item: (item["work"],
                                                           len(item["ids"])))
            else:
                selected = {"work": 0.0, "ids": []}
                bins.append(selected)
            selected["ids"].append(identifier)
            selected["work"] += work
        prototype = control_prototypes[key]
        for item in bins:
            row = dict(prototype)
            row["control_ids"] = ",".join(item["ids"])
            row["predicted_likelihood_work"] = item["work"]
            analysis.append(row)
    analysis.sort(key=lambda row: (
        row["action"] != "CELL", row["library"], row["modality"],
        row.get("barcode", ""), row.get("control_ids", "")))
    if len(extraction) >= 999 or len(analysis) >= 999:
        raise RuntimeError("reprojected array size must remain below 999 tasks")
    for index, row in enumerate(analysis):
        row["task_index"] = index
        row["workload_generation_id"] = generation
        extraction_row = extraction_by_key[(row["library"], row["modality"])]
        row["cache_generation_id"] = extraction_row["cache_generation_id"]
        row["cache_prefix"] = extraction_row["cache_prefix"]
        row["candidate_manifest"] = extraction_row["candidate_manifest"]
        stem = "cell" if row["action"] == "CELL" else "control"
        row["analysis_output"] = str(generation_root / "analysis" /
                                     f"{stem}_{index:04d}.tsv.gz")
        row["marker"] = str(generation_root / "markers" /
                            f"analysis_{index:04d}.complete.json")
        row["control_manifest"] = str(
            generation_root / "manifests" / "control_pairs_final.tsv")
    extraction_fields = list(extraction[0]) if extraction else []
    analysis_fields = list(analysis[0]) if analysis else []
    write_tsv(generation_root / "manifests" / "extraction_tasks.tsv",
              extraction, extraction_fields)
    write_tsv(generation_root / "manifests" / "analysis_tasks.tsv",
              analysis, analysis_fields)
    # Publish generation-specific proof markers for the immutable Library-12
    # cache only after the exact promoted task contracts have been frozen.
    for row in extraction:
        if row.get("cache_policy") != "IMMUTABLE_BENCHMARK_REUSE_REQUIRED":
            continue
        valid, problems, details = validate_task_artifact(
            row, "extraction", require_marker=False)
        if not valid:
            raise RuntimeError(
                "promoted Library-12 cache failed generation binding: " +
                ";".join(problems))
        atomic_json(row["marker"], {
            "schema_version": TASK_MARKER_SCHEMA,
            "operational_status": "COMPLETE", "utc": utc_now(),
            **details,
            "promotion_generation_id": generation,
            "promotion_source_generation_id": row[
                "promoted_cache_source_generation"],
            "promoted_cache_proof_sha256": row[
                "promoted_cache_proof_sha256"],
            "immutable_reuse_without_source_rescan": True,
        })
    source = top / "manifests" / "control_pairs_provisional.tsv"
    if not source.is_file():
        raise RuntimeError("reprojection requires the frozen provisional controls")
    shutil.copyfile(source, generation_root / "manifests" /
                    "control_pairs_provisional.tsv")
    for name in ("target_cells_and_matched_comparisons.tsv",
                 "target_comparison_balance.tsv",
                 "control_capacity_and_reuse.tsv",
                 "primary_calibration_audit.tsv.gz",
                 "primary_calibration_reference_roster.tsv.gz",
                 "frozen_targets_20260920.tsv", "frozen_target_comparisons.tsv", "unavailable_cells.tsv", "exact_commands.json"):
        source_path = top / name
        if not source_path.is_file():
            raise RuntimeError(f"reprojection requires frozen workload input {source_path}")
        shutil.copyfile(source_path, generation_root / name)

    projected_analysis_seconds_by_task = [
        finite(row.get("predicted_likelihood_work", ""), 0) * seconds_per_work
        for row in analysis]
    maximum_analysis_seconds = max(projected_analysis_seconds_by_task,
                                   default=0.0)
    analysis_wall_seconds = max(1800, int(math.ceil(
        maximum_analysis_seconds * 1.5 + 600)))
    extraction_source_bytes = []
    for row in extraction:
        source_bytes = sum(Path(row[field]).stat().st_size
                           for field in ("samples", "pileup_sites",
                                         "pileup_observations", "pileup_molecules")
                           if Path(row[field]).is_file())
        extraction_source_bytes.append(source_bytes)
    projected_extraction_seconds = [
        0.0 if row.get("cache_policy") ==
        "IMMUTABLE_BENCHMARK_REUSE_REQUIRED" else
        value * extraction_seconds_per_byte
        for row, value in zip(extraction, extraction_source_bytes)]
    extraction_wall_seconds = max(1800, int(math.ceil(
        max(projected_extraction_seconds, default=0.0) * 1.5 + 600)))

    def slurm_time(seconds):
        seconds = int(math.ceil(seconds / 60.0) * 60)
        days, remainder = divmod(seconds, 86400)
        hours, remainder = divmod(remainder, 3600)
        minutes, seconds = divmod(remainder, 60)
        prefix = f"{days}-" if days else ""
        return f"{prefix}{hours:02d}:{minutes:02d}:{seconds:02d}"

    inherited_resources = _validated_targeted_resources(
        blueprint.get("resource_policy", {}))
    resource_projection = {
        "schema_version": "joint_doublet_resource_projection_v4",
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "array_concurrency": concurrency,
        "analysis_memory": f"{requested_memory_gib}G",
        "analysis_wall_time": slurm_time(analysis_wall_seconds),
        "analysis_cpus": inherited_resources["analysis_cpus"],
        "extraction_memory": f"{extraction_memory_gib}G",
        "extraction_wall_time": slurm_time(extraction_wall_seconds),
        "extraction_cpus": inherited_resources["extraction_cpus"],
        "finalizer_memory": inherited_resources["finalizer_memory"],
        "finalizer_wall_time": inherited_resources["finalizer_wall_time"],
        "finalizer_cpus": inherited_resources["finalizer_cpus"],
        "gather_memory": inherited_resources["gather_memory"],
        "gather_wall_time": inherited_resources["gather_wall_time"],
        "gather_cpus": inherited_resources["gather_cpus"],
        "partition": inherited_resources["partition"],
        "memory_safety_margin": 1.50,
        "concurrency_memory_budget_gib": concurrency_memory_budget_gib,
        "concurrency_cpu_budget": concurrency_cpu_budget,
        "measured_safe_aggregate_io_mib_per_second": safe_aggregate_io,
        "measured_maximum_task_io_mib_per_second": maximum_task_io,
        "concurrency_caps": concurrency_caps,
        "limiting_resources": limiting_resources,
    }
    atomic_json(generation_root / "resource_projection.json",
                resource_projection)
    projected_resources = _validated_targeted_resources(resource_projection)
    for row in extraction:
        if row.get("cache_policy") != "IMMUTABLE_BENCHMARK_REUSE_REQUIRED":
            row["scheduler_memory_bytes"] = projected_resources[
                "extraction_scheduler_memory_bytes"]
            row["launcher_runtime_reserve_bytes"] = projected_resources[
                "extraction_launcher_runtime_reserve_bytes"]
            row["worker_memory_budget_bytes"] = projected_resources[
                "extraction_worker_memory_budget_bytes"]
            row["bounded_memory_limit_bytes"] = row[
                "worker_memory_budget_bytes"]
    for row in analysis:
        row["scheduler_memory_bytes"] = projected_resources[
            "analysis_scheduler_memory_bytes"]
        row["launcher_runtime_reserve_bytes"] = projected_resources[
            "analysis_launcher_runtime_reserve_bytes"]
        row["worker_memory_budget_bytes"] = projected_resources[
            "analysis_worker_memory_budget_bytes"]
    write_tsv(generation_root / "manifests" / "extraction_tasks.tsv",
              extraction, extraction_fields)
    write_tsv(generation_root / "manifests" / "analysis_tasks.tsv",
              analysis, analysis_fields)
    scripts = _write_targeted_scripts(
        generation_root, args.tool_bin_root, generation, extraction, analysis,
        benchmark=False,
        production_threads=resource_projection["analysis_cpus"],
        resource_policy=resource_projection)
    total_predicted_work = sum(finite(row.get("predicted_likelihood_work", ""), 0)
                               for row in analysis)
    projected_analysis_seconds = sum(projected_analysis_seconds_by_task)
    measured_cache_bytes = sum(int(row.get("written_cache_bytes", 0) or 0)
                               for row in extraction_measured)
    measured_records = sum(int(row.get(
        "measured_unique_evidence_records", 0) or 0)
        for row in extraction_measured)
    cache_bytes_per_record = measured_cache_bytes / measured_records \
        if measured_records else math.nan
    projected_records = sum(int(row.get("projected_observation_records", 0) or 0) +
                            int(row.get("projected_molecule_records", 0) or 0)
                            for row in extraction)
    projected_cache_bytes = projected_records * cache_bytes_per_record \
        if math.isfinite(cache_bytes_per_record) else math.nan
    measured_output_bytes = sum(int(row.get("output_bytes", 0) or 0)
                                for row in analysis_measured)
    measured_analysis_work = sum(finite(row.get(
        "predicted_likelihood_work", ""), 0) for row in analysis_measured)
    projected_output_bytes = measured_output_bytes * total_predicted_work / \
        measured_analysis_work if measured_analysis_work else math.nan
    projected_extraction_task_seconds = sum(projected_extraction_seconds)
    measured_optimizer_calls = sum(int(row.get("optimizer_calls", 0) or 0)
                                   for row in analysis_measured)
    measured_likelihood_evaluations = sum(int(row.get(
        "likelihood_evaluations", 0) or 0) for row in analysis_measured)
    projected_optimizer_calls = measured_optimizer_calls * \
        total_predicted_work / measured_analysis_work \
        if measured_analysis_work else math.nan
    projected_likelihood_evaluations = measured_likelihood_evaluations * \
        total_predicted_work / measured_analysis_work \
        if measured_analysis_work else math.nan
    lanes = [0.0] * concurrency
    for seconds in sorted(projected_analysis_seconds_by_task, reverse=True):
        lane = min(range(concurrency), key=lambda index: lanes[index])
        lanes[lane] += seconds
    projected_analysis_critical_seconds = max(lanes, default=0.0)
    prebenchmark_policy = blueprint.get("resource_policy", {})
    accounting = [{
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "workload_generation_id": generation,
        "supersedes": previous,
        "benchmark_measurement_digest": measurement_digest,
        "measured_analysis_peak_rss_gib": analysis_peak_rss_gib,
        "measured_extraction_peak_rss_gib": extraction_peak_rss_gib,
        "memory_safety_margin": 1.50,
        "analysis_requested_memory_gib": requested_memory_gib,
        "extraction_requested_memory_gib": extraction_memory_gib,
        "measured_analysis_seconds_per_predicted_work": seconds_per_work,
        "measured_extraction_seconds_per_source_byte":
            extraction_seconds_per_byte,
        "analysis_wall_limit": resource_projection["analysis_wall_time"],
        "extraction_wall_limit": resource_projection["extraction_wall_time"],
        "target_control_bin_predicted_work": target_bin_work,
        "reprojected_control_bin_tasks": sum(
            row["action"] == "CONTROL_BIN" for row in analysis),
        "reprojected_cell_tasks": sum(row["action"] == "CELL"
                                      for row in analysis),
        "projected_source_bytes": sum(extraction_source_bytes),
        "projected_source_bytes_read_after_cache_reuse": sum(
            value for row, value in zip(extraction, extraction_source_bytes)
            if row.get("cache_policy") !=
            "IMMUTABLE_BENCHMARK_REUSE_REQUIRED"),
        "projected_full_source_scans": sum(
            row.get("cache_policy") != "IMMUTABLE_BENCHMARK_REUSE_REQUIRED"
            for row in extraction),
        "immutable_library12_cache_reuse_tasks": sum(
            row.get("cache_policy") == "IMMUTABLE_BENCHMARK_REUSE_REQUIRED"
            for row in extraction),
        "projected_unique_evidence_records": projected_records,
        "measured_cache_bytes_per_unique_record": cache_bytes_per_record,
        "projected_written_cache_bytes": projected_cache_bytes,
        "projected_analysis_output_bytes": projected_output_bytes,
        "candidate_expanded_rows_avoided": prebenchmark_policy.get(
            "legacy_candidate_expanded_rows_avoided", "UNAVAILABLE"),
        "projected_cache_reduction_factor": prebenchmark_policy.get(
            "projected_reduction_factor", "UNAVAILABLE"),
        "projected_optimizer_calls": projected_optimizer_calls,
        "projected_likelihood_evaluations":
            projected_likelihood_evaluations,
        "target_cells": prebenchmark_policy.get(
            "selected_target_cells", "UNAVAILABLE"),
        "comparison_cells": prebenchmark_policy.get(
            "matched_comparison_cells", "UNAVAILABLE"),
        "control_pairs": prebenchmark_policy.get(
            "control_unique_pairs", "UNAVAILABLE"),
        "control_fractions": ",".join(map(str, CONTROL_FRACTIONS)),
        "channels": "SITE,MOLECULE",
        "downsample_replicates_per_fraction": DOWNSAMPLE_REPLICATES,
        "null_replicates_per_channel": CELL_NULL_REPLICATES,
        "projected_extraction_task_hours":
            projected_extraction_task_seconds / 3600,
        "projected_extraction_core_hours":
            projected_extraction_task_seconds *
            resource_projection["extraction_cpus"] / 3600,
        "projected_analysis_task_hours": projected_analysis_seconds / 3600
            if math.isfinite(projected_analysis_seconds) else math.nan,
        "projected_analysis_core_hours": projected_analysis_seconds *
            resource_projection["analysis_cpus"] / 3600
            if math.isfinite(projected_analysis_seconds) else math.nan,
        "analysis_task_count": len(analysis),
        "extraction_task_count": len(extraction),
        "analysis_cpus_per_task": resource_projection["analysis_cpus"],
        "extraction_cpus_per_task": resource_projection["extraction_cpus"],
        "maximum_projected_analysis_rss_gib": requested_memory_gib,
        "maximum_projected_extraction_rss_gib": extraction_memory_gib,
        "resident_memory_formula":
            "maximum measured benchmark RSS * 1.5 safety margin, rounded up",
        "projected_throttled_analysis_critical_path_hours":
            projected_analysis_critical_seconds / 3600,
        "projected_end_to_end_critical_path_hours": (
            max(projected_extraction_seconds, default=0.0) +
            projected_analysis_critical_seconds
        ) / 3600,
        "queue_time": "EXCLUDED",
        "array_concurrency": concurrency,
        "concurrency_memory_cap": concurrency_caps["MEMORY"],
        "concurrency_cpu_cap": concurrency_caps["CPU"],
        "concurrency_io_cap": concurrency_caps["IO"],
        "concurrency_hard_cap": concurrency_caps["HARD_CAP"],
        "concurrency_limiting_resources": ",".join(limiting_resources),
        "concurrency_basis": (
            "minimum of measured-memory, configured-CPU, measured-I/O, "
            "and declared hard caps"),
        "disk_projection_bytes": projected_cache_bytes + projected_output_bytes,
        "state": "RENDERED_UNSUBMITTED",
    }]
    write_tsv(generation_root / "workload_accounting.tsv", accounting)
    new_blueprint = dict(blueprint)
    new_blueprint.update({
        "workload_generation_id": generation, "state": "RENDERED_UNSUBMITTED",
        "supersedes": previous, "extraction_tasks": extraction,
        "analysis_tasks": analysis, "scripts": scripts,
        "resource_policy": accounting[0], "created_utc": utc_now(),
    })
    atomic_json(generation_root / "workload_blueprint.json", new_blueprint)
    atomic_json(top / "PREBENCHMARK_SUPERSEDED.json", {
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "workload_generation_id": previous, "state": "SUPERSEDED",
        "superseded_by": generation, "utc": utc_now(),
    })
    atomic_json(top / "CURRENT_GENERATION.json", {
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "workload_generation_id": generation,
        "root": str(generation_root), "state": "RENDERED_UNSUBMITTED",
    })
    rendered_marker = {
        "schema_version": WORKLOAD_SCHEMA,
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": CALIBRATION_LIBRARY,
        "status": "RENDERED_UNSUBMITTED",
        "workload_generation_id": generation, "utc": utc_now(),
        "extraction_tasks": len(extraction), "analysis_tasks": len(analysis),
        "manifest_libraries": list(TARGET_LIBRARIES),
        "protected_libraries_accessed": sorted(
            set(TARGET_LIBRARIES) & set(PROTECTED_LIBRARIES)),
        "resource_projection": resource_projection,
    }
    atomic_json(generation_root / "TARGETED_WORKLOAD_RENDERED_UNSUBMITTED",
                rendered_marker)
    reprojected_paths = [
        generation_root / "workload_accounting.tsv",
        generation_root / "workload_blueprint.json",
        generation_root / "resource_projection.json",
        generation_root / "exact_commands.json",
        generation_root / "TARGETED_WORKLOAD_RENDERED_UNSUBMITTED",
        generation_root / "target_cells_and_matched_comparisons.tsv",
        generation_root / "target_comparison_balance.tsv",
        generation_root / "control_capacity_and_reuse.tsv",
        generation_root / "primary_calibration_audit.tsv.gz",
        generation_root / "primary_calibration_reference_roster.tsv.gz",
        generation_root / "frozen_targets_20260920.tsv",
        generation_root / "manifests", generation_root / "slurm_scripts",
    ]
    return_archive = compact_zip(generation_root, [
        (path, generation_root, "") for path in reprojected_paths
    ], "reprojected_workload_return.zip")
    print(json.dumps({"status": "RENDERED_UNSUBMITTED",
                      "generation_root": str(generation_root),
                      "return_archive": str(return_archive),
                      "accounting": accounting[0]}, indent=2))
    return 0


def subset_progression_rows(cells):
    def molecule_preserves_site_genetic(row):
        site_rna = clean(row.get("rna_site_genetic_second_state", ""))
        site_atac = clean(row.get("atac_site_genetic_second_state", ""))
        molecule_rna = clean(row.get("rna_molecule_genetic_second_state", ""))
        molecule_atac = clean(row.get("atac_molecule_genetic_second_state", ""))
        return bool(site_rna) and \
            site_rna == site_atac == molecule_rna == molecule_atac

    base_cohorts = (
        ("ORIGINAL_3183_BOTH_ASSAYS_P95", lambda row:
            truthy(row.get("original_site_both_assays_p95", ""))),
        ("ORIGINAL_2308_EXACT_WINNER_AGREEMENT", lambda row:
            truthy(row.get("original_site_exact_winner_agreement", ""))),
    )
    requirements = (
        ("STARTING_COHORT", lambda row: True),
        ("PHYSICAL_GENOTYPE_VISIBLE_NEW_DONOR_BOTH_P95", lambda row:
            truthy(row.get("genetic_source_both_assays_p95", ""))),
        ("PHYSICAL_NEW_DONOR_EXACT_AGREEMENT", lambda row:
            truthy(row.get("genetic_source_exact_winner_agreement", ""))),
        ("EXCLUDE_REPLACEMENT_BOUNDARY", lambda row:
            row.get("site_existing_evidence_category") !=
                "replacement-like or upper-boundary fit"),
        ("MIXTURE_COMPATIBLE_OR_UNCERTAIN", lambda row:
            row.get("site_existing_evidence_category") in {
                "mixture-compatible in both assays",
                "addition-compatible but uncertain"}),
        ("MOLECULE_EVALUABLE_AND_PRESERVES_SITE_GENETIC_WINNER",
         molecule_preserves_site_genetic),
        ("MOLECULE_MIXTURE_COMPATIBLE_OR_UNCERTAIN", lambda row:
            row.get("molecule_existing_evidence_category") in {
                "mixture-compatible in both assays",
                "addition-compatible but uncertain"}),
    )
    rows = []
    for cohort_name, cohort_predicate in base_cohorts:
        selected = [row for row in cells if cohort_predicate(row)]
        previous_count = None
        for stage, requirement in requirements:
            selected = [row for row in selected if requirement(row)]
            groups = Counter((
                row["library"],
                clean(row.get("rna_site_genetic_second_state", "") or
                      row.get("rna_site_unrestricted_second_state", "")))
                for row in selected)
            rows.append({
                "analysis_type": "DENOMINATOR_PROGRESSION",
                "endpoint": stage, "scope": "ALL_SEVEN",
                "starting_cohort": cohort_name,
                "cells": len(selected),
                "excluded_from_previous": (
                    previous_count - len(selected)
                    if previous_count is not None else ""),
            })
            previous_count = len(selected)
            for (library, contributor), count in sorted(groups.items()):
                rows.append({
                    "analysis_type": "DENOMINATOR_PROGRESSION",
                    "endpoint": stage,
                    "scope": "LIBRARY_AND_PROPOSED_CONTRIBUTOR",
                    "starting_cohort": cohort_name,
                    "library": library,
                    "proposed_contributor": contributor,
                    "cells": count,
                })
    return rows


def molecule_basis_rows(cells):
    groups = Counter()
    for row in cells:
        if not truthy(row.get("original_site_both_assays_p95", "")):
            continue
        contributor = clean(row.get("rna_site_unrestricted_second_state", ""))
        groups[(row["library"], contributor,
                row.get("rna_molecule_evidence_basis_class", "UNAVAILABLE"),
                row.get("atac_molecule_evidence_basis_class", "UNAVAILABLE"),
                truthy(row.get("molecule_retains_original_concordant_winner", "")))] += 1
    return [{
        "analysis_type": "MOLECULE_EVIDENCE_BASIS_SENSITIVITY",
        "endpoint": "ORIGINAL_SITE_FOCAL_CELLS",
        "scope": "LIBRARY_CANDIDATE_LINKAGE_BASIS",
        "library": library, "proposed_contributor": contributor,
        "rna_molecule_basis": rna_basis, "atac_molecule_basis": atac_basis,
        "concordant_winner_retained": retained, "cells": count,
    } for (library, contributor, rna_basis, atac_basis, retained), count
        in sorted(groups.items())]


def calibration_rows(cells, calibration_audit=None):
    rows = []
    for library in (f"lib{value}" for value in ALLOWED_LIBRARIES):
        subset = [row for row in cells if row["library"] == library]
        for endpoint in (
                "original_site_both_assays_p95",
                "global_calibration_library_unstratified_dual_p95",
                "simple_within_library_dual_p95"):
            rows.append({
                "analysis_type": "CALIBRATION_SENSITIVITY",
                "endpoint": endpoint, "scope": "PER_LIBRARY",
                "library": library, "cells": len(subset),
                "positive_cells": sum(truthy(row.get(endpoint, "")) for row in subset),
                "positive_fraction": mean([
                    float(truthy(row.get(endpoint, ""))) for row in subset]),
            })
    primary = [row for row in cells
               if int(row["library"].removeprefix("lib")) in PRIMARY_LIBRARIES]
    external_only = [row for row in primary
                     if truthy(row.get(
                         "global_calibration_library_unstratified_dual_p95", ""))
                     and not truthy(row.get("simple_within_library_dual_p95", ""))]
    within_only = [row for row in primary
                   if truthy(row.get("simple_within_library_dual_p95", ""))
                   and not truthy(row.get(
                       "global_calibration_library_unstratified_dual_p95", ""))]
    both = [row for row in primary
            if truthy(row.get("simple_within_library_dual_p95", "")) and
            truthy(row.get(
                "global_calibration_library_unstratified_dual_p95", ""))]
    calibration_label = f"Library-{CALIBRATION_LIBRARY}"
    rows.extend((
        {
            "analysis_type": "CALIBRATION_DIRECTION_DIAGNOSTIC",
            "endpoint": f"{calibration_label}-positive but within-library-negative",
            "scope": "PRIMARY_SIX_LIBRARIES", "cells": len(external_only),
            "calibration_library": CALIBRATION_LIBRARY,
            "interpretation": "external-reference inflation concern",
        },
        {
            "analysis_type": "CALIBRATION_DIRECTION_DIAGNOSTIC",
            "endpoint": f"within-library-positive but {calibration_label}-negative",
            "scope": "PRIMARY_SIX_LIBRARIES", "cells": len(within_only),
            "calibration_library": CALIBRATION_LIBRARY,
            "interpretation": "reverse-direction discrepancy",
        },
        {
            "analysis_type": "CALIBRATION_DIRECTION_DIAGNOSTIC",
            "endpoint": (
                f"positive under both global {calibration_label} and within library"),
            "scope": "PRIMARY_SIX_LIBRARIES", "cells": len(both),
            "calibration_library": CALIBRATION_LIBRARY,
            "interpretation": "shared positive set",
        },
    ))
    for audit in calibration_audit or []:
        rows.append({
            "analysis_type": "TARGET_CALIBRATION_SCHEME",
            "endpoint": endpoint_label(
                "SITE_GENETIC_SOURCE" if audit["evidence_channel"] == "SITE"
                else "MOLECULE_GENETIC_SOURCE"),
            "scope": "FROZEN_TARGET",
            **audit,
        })
    return rows


def objective_decisions(cells, group_results):
    """Apply prespecified inferential semantics without converting missingness.

    A positive point estimate alone is never a robustness result.  Inferential
    answers require an exchangeable null, a positive block-bootstrap interval,
    and the appropriate multiplicity-adjusted probability.  Descriptive
    addition/molecule summaries stay explicitly non-confirmatory.
    """
    def find_null(endpoint, scope, exclusion="NONE"):
        return next((row for row in group_results
                     if row.get("analysis_type") ==
                        "CANDIDATE_AND_COVERAGE_AWARE_CROSS_ASSAY_NULL" and
                     row.get("endpoint") == endpoint and
                     row.get("scope") == scope and
                     row.get("exclusion", "NONE") == exclusion), None)

    def null_components(row):
        if not row:
            return None
        result = {
            "denominator": int(row.get("denominator", 0) or 0),
            "exchangeable_denominator": int(row.get(
                "exchangeable_only_denominator", 0) or 0),
            "exchangeable_blocks": int(row.get(
                "exchangeable_multirow_blocks", 0) or 0),
            "effect": finite(row.get(
                "exchangeable_only_excess_fraction", "")),
            "interval_low": finite(row.get(
                "exchangeable_only_bootstrap_excess_95pct_low", "")),
            "interval_high": finite(row.get(
                "exchangeable_only_bootstrap_excess_95pct_high", "")),
            "raw_probability": finite(row.get(
                "exchangeable_only_empirical_upper_tail_probability", "")),
            "fixed_singletons": int(row.get("singleton_fixed_cells", 0) or 0),
        }
        result["valid"] = (
            result["denominator"] > 0 and
            result["exchangeable_denominator"] > 0 and
            result["exchangeable_blocks"] > 0 and
            all(math.isfinite(result[field]) for field in (
                "effect", "interval_low", "interval_high",
                "raw_probability")))
        return result

    def wilson(successes, denominator):
        if denominator <= 0:
            return math.nan, math.nan
        z = 1.959963984540054
        proportion = successes / denominator
        scale = 1 + z * z / denominator
        center = (proportion + z * z / (2 * denominator)) / scale
        radius = z * math.sqrt(
            proportion * (1 - proportion) / denominator +
            z * z / (4 * denominator * denominator)) / scale
        return center - radius, center + radius

    genetic = find_null("SITE_GENETIC_SOURCE", "PRIMARY_SIX_LIBRARIES")
    genetic_components = null_components(genetic)
    question1 = "unresolved"
    if genetic_components and genetic_components["valid"]:
        question1 = "yes" if (
            genetic_components["effect"] > 0 and
            genetic_components["interval_low"] > 0 and
            genetic_components["raw_probability"] <= 0.05
        ) else "no"
    sensitivities = [row for row in group_results
                     if row.get("scope") in {
                         "LEAVE_ONE_CANDIDATE_RERANK",
                         "LEAVE_ONE_PATTERN_RERANK"} and
                     row.get("endpoint") == "SITE_GENETIC_SOURCE"]
    leave_one_library = [row for row in group_results
                         if row.get("analysis_type") ==
                            "CANDIDATE_AND_COVERAGE_AWARE_CROSS_ASSAY_NULL" and
                         row.get("endpoint") == "SITE_GENETIC_SOURCE" and
                         row.get("scope") == "LEAVE_ONE_LIBRARY_OUT"]
    critical_sensitivities = sensitivities + leave_one_library
    sensitivity_components = [null_components(row)
                              for row in critical_sensitivities]
    sensitivity_family_size = len(sensitivity_components)
    for row, components in zip(critical_sensitivities,
                               sensitivity_components):
        raw = components["raw_probability"] if components else math.nan
        adjusted = min(1.0, raw * sensitivity_family_size) \
            if math.isfinite(raw) and sensitivity_family_size else math.nan
        row["decision_family_adjustment_rule"] = (
            "Bonferroni across all predeclared leave-one-library, "
            "leave-one-prominent-contributor, and leave-one-prominent-pattern "
            "genetic-source sensitivities")
        row["decision_family_adjusted_probability"] = adjusted
        if components is not None:
            components["adjusted_probability"] = adjusted
    if not sensitivity_components or any(
            not item or not item["valid"] for item in sensitivity_components):
        question2 = "unresolved"
    else:
        robust = all(
            item["effect"] > 0 and item["interval_low"] > 0 and
            item["adjusted_probability"] <= 0.05
            for item in sensitivity_components)
        question2 = "yes" if robust else "no"
    additions = [row for row in cells
                 if int(row["library"].removeprefix("lib")) in PRIMARY_LIBRARIES and
                 truthy(row.get("genetic_source_exact_winner_agreement", "")) and
                 row.get("site_existing_evidence_category") in {
                     "mixture-compatible in both assays",
                     "addition-compatible but uncertain"}]
    addition_libraries = {row["library"] for row in additions}
    question3 = "descriptive only"
    molecule_preserved = [row for row in additions
                          if truthy(row.get(
                              "molecule_genetic_exact_winner_agreement", "")) and
                          clean(row.get("rna_site_genetic_second_state", "")) ==
                          clean(row.get("rna_molecule_genetic_second_state", "")) and
                          row.get("molecule_existing_evidence_category") in {
                              "mixture-compatible in both assays",
                              "addition-compatible but uncertain"}]
    evaluable = [row for row in additions
                 if clean(row.get("rna_molecule_genetic_second_state", "")) and
                 clean(row.get("atac_molecule_genetic_second_state", ""))]
    question4 = "descriptive only" if evaluable else "unresolved"
    molecule_low, molecule_high = wilson(len(molecule_preserved), len(evaluable))
    calibration_evaluable = [row for row in additions
                             if all(clean(row.get(
                                 f"{modality}_site_primary_exact_mask_and_coverage_status",
                                 "")) == "AVAILABLE"
                                 for modality in ("rna", "atac"))]
    calibration_preserved = [row for row in calibration_evaluable
                             if truthy(row.get(
                                 "menu_coverage_calibrated_dual_p95", ""))]
    calibration_libraries = {
        row["library"] for row in calibration_preserved
    }
    calibration_low, calibration_high = wilson(
        len(calibration_preserved), len(calibration_evaluable))
    # Coverage adjustment cannot be promoted until every prespecified
    # after-matching/overlap balance diagnostic passes in the rendered targeted
    # workload.  The no-rescore calibration rows alone are informative but do
    # not satisfy that decision criterion.
    question5 = "unresolved"
    answers = (
        ("Does agreement remain above a candidate- and coverage-aware null?",
         question1, {
             "endpoint": "genetically distinguishable new contributor",
             "criterion": (
                 "valid exchangeable null; positive excess; 95% block-bootstrap "
                 "interval entirely above zero; upper-tail probability <=0.05"),
             **(genetic_components or {"valid": False}),
         }),
        ("Does the excess persist outside one dominant library/candidate pattern?",
         question2, {
             "required_sensitivities": len(critical_sensitivities),
             "valid_sensitivities": sum(bool(item and item["valid"])
                                          for item in sensitivity_components),
             "positive_interval_and_adjusted_significant_sensitivities": sum(
                 bool(item and item["valid"] and item["effect"] > 0 and
                      item["interval_low"] > 0 and
                      item["adjusted_probability"] <= 0.05)
                 for item in sensitivity_components),
             "adjustment_rule": (
                 "Bonferroni across the combined predeclared endpoint-matched "
                 "sensitivity family"),
             "worst_adjusted_probability": max((item.get(
                 "adjusted_probability", math.nan)
                 for item in sensitivity_components if item),
                 default=math.nan),
             "minimum_exchangeable_excess_fraction": min((item["effect"]
                 for item in sensitivity_components
                 if item and math.isfinite(item["effect"])),
                 default=math.nan),
             "minimum_exchangeable_interval_low": min((item["interval_low"]
                 for item in sensitivity_components
                 if item and math.isfinite(item["interval_low"])),
                 default=math.nan),
         } if critical_sensitivities else None),
        ("Is there a subset compatible with addition rather than replacement/boundary behavior?",
         question3, {"cells": len(additions), "libraries": len(addition_libraries),
                     "conclusion": "plausible but narrow; not established physical doublets"}),
        ("Does molecule-balanced evidence preserve that subset?",
         question4, {"preserved": len(molecule_preserved),
                     "molecule_evaluable": len(evaluable),
                     "fraction_retained": len(molecule_preserved) / len(evaluable)
                         if evaluable else math.nan,
                     "fraction_wilson_low": molecule_low,
                     "fraction_wilson_high": molecule_high,
                     "preserved_libraries": sorted({row["library"]
                                                    for row in molecule_preserved}),
                     "criterion": (
                         "same RNA/ATAC molecule contributor, equal to the site "
                         "contributor, with compatible available evidence"),
                     "conclusion": "limited under the legacy reference"}),
        ("Does coverage-conditional calibration preserve it?",
         question5, {
             "coverage_conditioned_targets_evaluable": len(calibration_evaluable),
             "calibration_preserved_addition_cells":
                 len(calibration_preserved),
             "calibration_preserved_libraries":
                 sorted(calibration_libraries),
             "fraction_retained": len(calibration_preserved) /
                 len(calibration_evaluable) if calibration_evaluable else math.nan,
             "fraction_wilson_low": calibration_low,
             "fraction_wilson_high": calibration_high,
             "contributor_concordance_required": True,
             "prespecified_preservation_criterion": (
                 "dual-assay P95 for the same evaluated contributor plus all "
                 "prespecified post-match/overlap balance diagnostics <=0.10"),
             "prespecified_preservation_criterion_met": False,
             "unresolved_reason": (
                 "targeted after-matching and overlap-weighted balance must be "
                 "validated before an adjusted preservation claim"),
         }),
    )
    rows = []
    supporting_tables = {
        1: "group_results.tsv",
        2: "group_results.tsv",
        3: "focal_cells_site_and_molecule.tsv.gz;group_results.tsv",
        4: "focal_cells_site_and_molecule.tsv.gz;group_results.tsv",
        5: "all_cells_reanalysis.tsv.gz;group_results.tsv",
    }
    for index, (question, answer, evidence) in enumerate(answers, 1):
        rows.append({
            "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
            "calibration_library": CALIBRATION_LIBRARY,
            "inference_layer": "LEGACY_DESCRIPTIVE_RECONSTRUCTION",
            "question_number": index, "question": question, "answer": answer,
            "estimate": json.dumps(evidence, sort_keys=True, default=str)
                if evidence is not None else "UNAVAILABLE",
            "supporting_table": supporting_tables[index],
            "interpretation_limit": (
                "historical stored fields only; not scorer-authoritative and "
                "not a validated physical-doublet call set"),
        })
    recommendation = (
        "Cross-assay site concordance is statistically non-random; credible "
        "addition evidence is plausible but narrow; most agreement is "
        "replacement/boundary-like; physical doublets are not yet validated.")
    return rows, recommendation


def slurm_usage():
    job_id = clean(os.environ.get("SLURM_JOB_ID", ""))
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return {
        "job_id": job_id or "NOT_A_SLURM_JOB",
        "user_cpu_seconds": usage.ru_utime,
        "system_cpu_seconds": usage.ru_stime,
        "maximum_resident_set_kb": usage.ru_maxrss,
        "scheduler_accounting_note": (
            "final sacct is intentionally collected by read-only status after exit"),
    }


def compact_zip(root, entries, name):
    """Build an atomic archive from explicit logical-path mappings.

    Each entry is ``(live_source, member_root, archive_prefix)``.  No basename
    fallback exists: external trees must declare their namespace explicitly.
    """
    root = Path(root).resolve()
    archive_path = root / name
    temporary = archive_path.with_name(f".{archive_path.name}.tmp.{os.getpid()}")
    archived_sources = {}

    def normalized_destination(source, member_root, prefix):
        source = Path(source)
        member_root = Path(member_root).resolve()
        if source.is_symlink():
            raise RuntimeError(f"ZIP source symlink rejected: {source}")
        resolved = source.resolve(strict=True)
        try:
            relative = resolved.relative_to(member_root)
        except ValueError as error:
            raise RuntimeError(
                f"ZIP source escapes declared member root: {source}") from error
        destination = Path(prefix) / relative
        if destination.is_absolute() or ".." in destination.parts:
            raise RuntimeError(f"unsafe ZIP destination: {destination}")
        arcname = destination.as_posix()
        if not arcname or arcname.startswith("/"):
            raise RuntimeError(f"invalid ZIP destination: {arcname}")
        return resolved, arcname

    def add_file(archive, path, member_root, prefix):
        resolved, arcname = normalized_destination(path, member_root, prefix)
        prior = archived_sources.get(arcname)
        if prior is not None:
            raise RuntimeError(
                f"ZIP member collision: {arcname}: {prior} versus {resolved}")
        archive.write(resolved, arcname=arcname)
        archived_sources[arcname] = resolved

    try:
        with zipfile.ZipFile(temporary, "w", zipfile.ZIP_DEFLATED,
                             compresslevel=9) as archive:
            for source, member_root, prefix in entries:
                source = Path(source)
                if not source.exists():
                    raise RuntimeError(f"required ZIP source is missing: {source}")
                if source.is_symlink():
                    raise RuntimeError(f"ZIP source symlink rejected: {source}")
                if source.is_file():
                    add_file(archive, source, member_root, prefix)
                elif source.is_dir():
                    for child in sorted(source.rglob("*")):
                        if child.is_symlink():
                            raise RuntimeError(
                                f"ZIP tree symlink rejected: {child}")
                        if child.is_file():
                            add_file(archive, child, member_root, prefix)
        with zipfile.ZipFile(temporary) as archive:
            names = archive.namelist()
            if len(names) != len(set(names)):
                raise RuntimeError("ZIP contains duplicate logical destinations")
            for info in archive.infolist():
                parts = Path(info.filename).parts
                if info.filename.startswith("/") or ".." in parts or \
                        ((info.external_attr >> 16) & 0o170000) == 0o120000:
                    raise RuntimeError(
                        f"unsafe ZIP member rejected: {info.filename}")
            corrupt = archive.testzip()
            if corrupt:
                raise RuntimeError(f"ZIP integrity failure at {corrupt}")
        os.replace(temporary, archive_path)
    except Exception:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise
    return archive_path


def public_row(row):
    return {key: value for key, value in row.items() if not key.startswith("_")}


def run_no_rescore(args):
    global CALIBRATION_LIBRARY
    CALIBRATION_LIBRARY = int(args.calibration_library)
    if CALIBRATION_LIBRARY != 25:
        raise RuntimeError("this pilot requires calibration Library 25")
    if not re.fullmatch(r"no_rescore_[0-9a-f]{16}", clean(args.generation)):
        raise RuntimeError("no-rescore workload generation identifier is invalid")
    phases = PhaseAccounting()
    libraries = tuple(int(value) for value in args.libraries)
    if libraries != ALLOWED_LIBRARIES:
        raise RuntimeError(
            "run requires exactly Libraries " + " ".join(map(str, ALLOWED_LIBRARIES)))
    roots = [args.stage_root, args.existing_stage_root, args.historical_gather_root,
             args.task_search_root, args.targeted_output_root, args.identity_tool,
             args.tool_bin_root, args.frozen_targets]
    if any(not Path(value).is_absolute() for value in roots):
        raise RuntimeError("all no-rescore paths must be absolute")
    if any(f"lib{value}" in "|".join(roots).lower()
           for value in PROTECTED_LIBRARIES):
        raise RuntimeError("a configured path unexpectedly names a protected library")
    stage_root = Path(args.stage_root)
    aggregate_root = stage_root / "aggregate"
    analysis_root = stage_root / "analysis"
    manifest_root = stage_root / "manifests"
    for path in (aggregate_root, analysis_root, manifest_root):
        path.mkdir(parents=True, exist_ok=True)

    frozen_manifest, frozen_rows = load_frozen_targets(args.frozen_targets)
    frozen_keys = {(f"lib{row['library']}", row["barcode"])
                   for row in frozen_rows}
    write_tsv(manifest_root / "joint_doublet_frozen_targets_20260920.tsv", frozen_rows)
    frozen_comparisons = load_frozen_comparisons(
        args.frozen_comparisons or args.frozen_targets, frozen_keys)
    excluded_reference_keys = frozen_keys | {
        (row["library"], row["matched_comparison_barcode"]) for row in frozen_comparisons}
    phases.finish("validate_and_copy_frozen_targets", targets=len(frozen_keys))

    # The upstream task manifest is the source allowlist.  Validate it and the
    # complete RNA/ATAC candidate-menu contents before discover_inputs() is
    # allowed to stat, gzip-probe, or open any ledger/score source.
    task_manifest, task_rows = load_task_rows(
        args.existing_stage_root, libraries)
    validated_candidate_manifests = \
        validate_no_rescore_candidate_manifests(task_rows, libraries)
    phases.finish(
        "validate_source_manifests_before_discovery",
        task_rows=len(task_rows),
        candidate_manifests=len(validated_candidate_manifests))

    (sources, inventory, warnings, task_manifest, task_rows,
     access_audit, accessed_libraries) = discover_inputs(
        args.existing_stage_root, args.task_search_root, libraries,
        task_manifest, task_rows)
    write_tsv(analysis_root / "input_open_audit.tsv", access_audit)
    phases.finish("resolve_completed_inputs", source_records=len(inventory))
    source_manifest = manifest_root / "no_rescore_input_manifest.tsv"
    write_tsv(source_manifest, sources,
              ("library",) + INPUT_ROLES)
    command = [
        sys.executable, args.identity_tool, "joint-aggregate",
        "--input-manifest", str(source_manifest),
        "--output-root", str(aggregate_root),
        "--libraries", *[str(value) for value in libraries],
        "--calibration-library", str(args.calibration_library),
        "--allow-missing",
    ]
    atomic_text(stage_root / "commands_actually_used.sh",
                "#!/bin/bash\n" + " ".join(shlex.quote(value) for value in command) +
                "\n")
    aggregate_run = subprocess.run(command, capture_output=True, text=True, check=False)
    atomic_text(stage_root / "aggregate_command.stdout.log", aggregate_run.stdout)
    atomic_text(stage_root / "aggregate_command.stderr.log", aggregate_run.stderr)
    if aggregate_run.returncode:
        raise RuntimeError(
            "corrected no-rescore aggregation failed: " +
            (aggregate_run.stderr.strip() or aggregate_run.stdout.strip()))
    phases.finish("aggregate_completed_score_tables")

    ledger_path = aggregate_root / "joint_doublet_cell_ledger.tsv.gz"
    candidates_path = aggregate_root / "joint_doublet_candidate_scores.tsv.gz"
    if not ledger_path.is_file() or not candidates_path.is_file():
        raise RuntimeError("aggregation did not produce its ledger and candidate table")
    ledgers, ledger_by_key = load_ledgers(ledger_path)
    (cells, opportunity, opportunity_denominators, _candidate_counts,
     score_rows_by_library, coverage_edges) = \
        stream_candidate_analysis(candidates_path, ledger_by_key)
    phases.finish("stream_and_rerank_candidate_tables",
                  cells=len(cells), candidate_rows=sum(score_rows_by_library.values()))
    reuse_frozen_comparisons(cells, frozen_comparisons)
    planned_parents = control_parent_rows(cells, excluded_reference_keys)
    planned_controls, planned_control_accounting = select_control_pairs(planned_parents)
    excluded_reference_keys.update((row["library"], row[field]) for row in planned_controls
                                   for field in ("recipient_barcode", "source_barcode"))
    changed, calibration_audit, calibration_roster, calibration_checks = \
        apply_calibration_sensitivities(cells, frozen_keys, excluded_reference_keys)
    phases.finish("calibration_reconstruction",
                  calibration_rows=len(calibration_audit),
                  reference_roster_rows=len(calibration_roster))
    group_results = opportunity_rows(opportunity, opportunity_denominators)
    group_results.extend(grouped_evidence_rows(cells))
    group_results.extend(occupancy_descriptive_rows(cells))
    group_results.extend(subset_progression_rows(cells))
    group_results.extend(molecule_basis_rows(cells))
    group_results.extend(calibration_rows(cells, calibration_audit))
    null_rows = run_null_suite(cells, args.permutations)
    group_results.extend(null_rows)
    dominant_rows, candidate_wins, pattern_wins = dominant_sensitivity_rows(
        cells, args.permutations, excluded_reference_keys)
    group_results.extend(dominant_rows)
    coverage_rows, _coverage_matches = coverage_adjustment_rows(
        cells, frozen_keys, calibration_audit)
    group_results.extend(coverage_rows)
    # objective_decisions annotates the complete predeclared sensitivity
    # family with its multiplicity adjustment; run it before group_results is
    # serialized so the audit table contains the exact criteria consumed.
    decisions, recommendation = objective_decisions(cells, group_results)
    phases.finish("null_sensitivity_and_coverage_analysis",
                  group_rows=len(group_results))

    denominator = inventory + denominator_rows(cells, score_rows_by_library)
    write_tsv(analysis_root / "denominator_and_source_inventory.tsv", denominator)
    write_tsv(analysis_root / "coverage_quintile_edges.tsv", coverage_edges)
    write_tsv(analysis_root / "all_cells_reanalysis.tsv.gz",
              (public_row(row) for row in cells))
    focal = [public_row(row) for row in cells
             if truthy(row.get("original_site_both_assays_p95", ""))]
    write_tsv(analysis_root / "focal_cells_site_and_molecule.tsv.gz", focal)
    write_tsv(analysis_root / "group_results.tsv", group_results)
    write_tsv(analysis_root / "changed_cells.tsv.gz", changed)
    write_tsv(analysis_root / "calibration_audit.tsv.gz", calibration_audit)
    write_tsv(analysis_root / "calibration_reference_roster.tsv.gz",
              calibration_roster)
    write_tsv(analysis_root / "calibration_reproduction_checks.tsv",
              calibration_checks)
    targets, _matches, _mapping, _balance = reuse_frozen_comparisons(cells, frozen_comparisons)
    target_evidence = frozen_target_evidence_rows(targets, changed)
    write_tsv(analysis_root / "frozen_59_target_evidence.tsv.gz",
              target_evidence)
    write_tsv(analysis_root / "decision_evidence.tsv", decisions)
    phases.finish("write_corrected_analysis_tables",
                  frozen_target_rows=len(target_evidence))

    expected_files_complete = all(
        row["exists"] and row["gzip_envelope_status"] == "PASS" and
        row["schema_status"] == "PASS"
        for row in inventory)
    molecule_channels = {
        (row["library"], modality)
        for row in cells for modality in ("rna", "atac")
        if clean(row.get(f"{modality}_molecule_unrestricted_second_state", ""))
    }
    missing_molecule_channels = [
        f"lib{library}:{modality.upper()}"
        for library in libraries for modality in ("rna", "atac")
        if (f"lib{library}", modality) not in molecule_channels
    ]
    original_focal = sum(truthy(row.get("original_site_both_assays_p95", ""))
                         for row in cells)
    original_agreement = sum(
        truthy(row.get("original_site_exact_winner_agreement", "")) for row in cells)
    if original_focal != 3183:
        warnings.append({
            "library": "ALL_SEVEN", "warning": "REFERENCE_DENOMINATOR_MISMATCH",
            "detail": f"expected_original_focal=3183;observed={original_focal}",
        })
    if original_agreement != 2308:
        warnings.append({
            "library": "ALL_SEVEN", "warning": "REFERENCE_AGREEMENT_MISMATCH",
            "detail": f"expected_original_agreement=2308;observed={original_agreement}",
        })
    if missing_molecule_channels:
        warnings.append({
            "library": "ALL_SEVEN", "warning": "SITE_AND_MOLECULE_INCOMPLETE",
            "detail": ",".join(missing_molecule_channels),
        })
    for check in calibration_checks:
        if check["reproduction_status"].startswith("DIFFERENCE_"):
            warnings.append({
                "library": "FROZEN_59",
                "warning": "AUDIT_CALIBRATION_REPRODUCTION_DIFFERENCE",
                "detail": json.dumps(check, sort_keys=True),
            })
    for row in null_rows:
        if row.get("fixed_singleton_materiality_warning") == "YES":
            warnings.append({
                "library": row.get("library") or "PRIMARY_SIX",
                "warning": "FIXED_SINGLETON_MATERIALLY_AFFECTS_NULL_EFFECT",
                "detail": f"{row.get('endpoint')}:{row.get('scope')}:{row.get('exclusion')}",
            })
    if args.targeted_input_manifest:
        if not Path(args.targeted_input_manifest).is_absolute():
            raise RuntimeError("targeted input manifest must be absolute")
        seen_inputs = set()
        for row in read_tsv(args.targeted_input_manifest):
            library, modality = library_name(row["library"]), row["modality"].lower()
            if int(library[3:]) not in ALLOWED_LIBRARIES or modality not in {"rna","atac"} or (library,modality) in seen_inputs:
                raise RuntimeError("invalid/duplicate original targeted input row")
            seen_inputs.add((library,modality))
            for field in ("samples","pileup_sites","pileup_observations","pileup_molecules"):
                value = row.get(field, "")
                if not Path(value).is_absolute() or _path_mentions_protected(value):
                    raise RuntimeError("original targeted manifest has an invalid input path")
                task_rows[library][f"{modality}_{field}"] = value
                if field == "samples":
                    task_rows[library][f"{modality}_pileup_samples"] = value
    targeted_accounting, targeted_blueprint = render_targeted_workload(
        args.targeted_output_root, candidates_path, cells, task_rows,
        args.tool_bin_root, frozen_keys, frozen_manifest,
        resource_policy={
            "extraction_cpus": 1,
            "extraction_memory": args.targeted_extraction_memory,
            "analysis_cpus": args.targeted_analysis_cpus,
            "analysis_memory": args.targeted_analysis_memory,
            "finalizer_cpus": args.targeted_finalizer_cpus,
            "finalizer_memory": args.targeted_finalizer_memory,
            "gather_cpus": args.targeted_gather_cpus,
            "gather_memory": args.targeted_gather_memory,
            "extraction_wall_time": args.targeted_time,
            "analysis_wall_time": args.targeted_time,
            "partition": args.targeted_partition,
        },
        model_parameters={
            "rna_error_ref": args.rna_error_ref,
            "rna_error_alt": args.rna_error_alt,
            "atac_error_ref": args.atac_error_ref,
            "atac_error_alt": args.atac_error_alt,
            "min_evidence": args.min_evidence,
            "max_second_fraction": args.max_second_fraction,
        },
        calibration_audit=calibration_audit,
        calibration_roster=calibration_roster,
        frozen_comparisons=frozen_comparisons,
        control_design=(planned_parents, planned_controls, planned_control_accounting))
    phases.finish("render_normalized_targeted_and_benchmark_workloads")

    usage = slurm_usage()
    protected_accessed = sorted(
        set(accessed_libraries) & set(PROTECTED_LIBRARIES))
    targeted_inputs_ready = all(
        row.get("input_status") == "READY"
        for row in targeted_blueprint.get("extraction_tasks", []))
    if not targeted_inputs_ready:
        warnings.append({
            "library": "TARGETED_12_20_29",
            "warning": "TARGETED_EXTRACTION_INPUT_CONTRACT_INCOMPLETE",
            "detail": ";".join(row.get("input_issues", "") for row in
                targeted_blueprint.get("extraction_tasks", [])
                if row.get("input_status") != "READY"),
        })
    warning_counts = Counter(row["warning"] for row in warnings)
    phases.finish("final_accounting", warnings=len(warnings))
    write_tsv(analysis_root / "phase_runtime_accounting.tsv", phases.rows)
    write_tsv(analysis_root / "process_resource_usage.tsv", [usage])
    # This is intentionally the final warning-table write: calibration,
    # matching, balance, control selection, rendering, and package-input
    # discovery have all completed.
    write_tsv(analysis_root / "warnings_and_exclusions.tsv", warnings,
              ("library", "warning", "detail"))
    required_artifacts_complete = (
        expected_files_complete and len(target_evidence) == 59 and
        set(accessed_libraries) == set(libraries) and
        not protected_accessed and targeted_inputs_ready)
    operational_status = "EXECUTION_COMPLETE" if required_artifacts_complete \
        else "TECHNICAL_FAILURE"
    scientific_status = "COMPLETE_WITH_SCIENTIFIC_WARNINGS" \
        if required_artifacts_complete else "UNRESOLVED"
    top_line = "Analysis execution complete; physical doublets are not yet validated." \
        if required_artifacts_complete else \
        "Analysis execution is technically incomplete; biological interpretation is unresolved."
    readme = f"""# Joint RNA/ATAC doublet no-rescore completion

**{top_line}**

Operational status: **{operational_status}**  
Scientific status: **{scientific_status}**

This derivative analysis read only completed ledgers, manifests, and score
tables for Libraries 7, 9, 12, 17, 20, 25, and 29. Libraries 19, 35, and 38
were not accessed. No preparation, pileup generation, RNA/ATAC rescoring,
synthetic-control execution, downsampling execution, or cell-conditional-null
execution occurred in this no-rescore job.

- Ledger cells: {len(cells):,}
- Unique candidate rows: {sum(score_rows_by_library.values()):,}
- Reproduced unrestricted both-assays-P95 cells: {original_focal:,}
- Reproduced exact winner agreements: {original_agreement:,}
- Frozen target cells: {len(target_evidence):,} (32 Library 12, 26 Library 20, 1 Library 29)
- Strict site-mixture targets: {sum(row['strict_site_five'] for row in target_evidence):,}
- Legacy exact molecule-supported targets: {sum(row['legacy_exact_molecule_supported_seven'] for row in target_evidence):,}
- Strict/exact overlap: {sum(row['strict_exact_overlap_three'] for row in target_evidence):,}
- Recommendation for boss-chat review: `{recommendation}`
- Targeted workload and bounded benchmark state: `RENDERED_UNSUBMITTED`
- Warning records: {len(warnings):,}; by type: `{json.dumps(dict(sorted(warning_counts.items())), sort_keys=True)}`
- In-process accounting: `{json.dumps(usage, sort_keys=True)}`

The decision table separates supported, contradicted, unresolved, and
descriptive-only questions. The 59-cell table preserves RNA UMI/gene-linked
molecules and ATAC read-name-based units as distinct evidence bases. The
rendered targeted workload is the next validation stage; no targeted or
benchmark task was submitted by this action.
"""
    atomic_text(analysis_root / "README_COMPLETION.md", readme)
    targeted_root = Path(args.targeted_output_root)
    targeted_paths = [
        targeted_root / "target_cells_and_matched_comparisons.tsv",
        targeted_root / "target_comparison_balance.tsv",
        targeted_root / "frozen_target_comparisons.tsv",
        targeted_root / "unavailable_cells.tsv",
        targeted_root / "primary_calibration_audit.tsv.gz",
        targeted_root / "primary_calibration_reference_roster.tsv.gz",
        targeted_root / "control_parent_eligibility.tsv.gz",
        targeted_root / "control_capacity_and_reuse.tsv",
        targeted_root / "frozen_targets_20260920.tsv",
        targeted_root / "workload_accounting.tsv",
        targeted_root / "resource_projection.json",
        targeted_root / "workload_blueprint.json",
        targeted_root / "CURRENT_GENERATION.json",
        targeted_root / "exact_commands.json",
        targeted_root / "TARGETED_WORKLOAD_RENDERED_UNSUBMITTED",
        targeted_root / "manifests",
        targeted_root / "slurm_scripts",
        targeted_root / "benchmark_unsubmitted" / "manifests",
        targeted_root / "benchmark_unsubmitted" / "slurm_scripts",
        targeted_root / "benchmark_unsubmitted" /
            "BENCHMARK_RENDERED_UNSUBMITTED",
    ]
    targeted_archive = compact_zip(
        stage_root, [(path, targeted_root, "proposal")
                     for path in targeted_paths],
        "proposed_targeted_workload.zip")
    result_paths = [
        analysis_root / "README_COMPLETION.md",
        analysis_root / "denominator_and_source_inventory.tsv",
        analysis_root / "input_open_audit.tsv",
        analysis_root / "coverage_quintile_edges.tsv",
        analysis_root / "all_cells_reanalysis.tsv.gz",
        analysis_root / "focal_cells_site_and_molecule.tsv.gz",
        analysis_root / "frozen_59_target_evidence.tsv.gz",
        analysis_root / "group_results.tsv",
        analysis_root / "calibration_audit.tsv.gz",
        analysis_root / "calibration_reference_roster.tsv.gz",
        analysis_root / "calibration_reproduction_checks.tsv",
        analysis_root / "changed_cells.tsv.gz",
        analysis_root / "decision_evidence.tsv",
        analysis_root / "warnings_and_exclusions.tsv",
        analysis_root / "phase_runtime_accounting.tsv",
        analysis_root / "process_resource_usage.tsv",
        manifest_root / "joint_doublet_frozen_targets_20260920.tsv",
        stage_root / "commands_actually_used.sh",
        targeted_archive,
    ]
    preliminary_archive = compact_zip(
        stage_root,
        [(path, stage_root, "") for path in result_paths] +
        [(path, targeted_root, "proposal") for path in targeted_paths],
        "preliminary_no_rescore_results.zip")
    marker = {
        "schema_version": SCHEMA,
        "scientific_method_version": SCIENTIFIC_METHOD_VERSION,
        "calibration_library": 25,
        "operational_status":
            "ANALYSIS_PROCESS_COMPLETE_AWAITING_ACCOUNTING_VALIDATION"
            if required_artifacts_complete else "TECHNICAL_FAILURE",
        "scientific_status": scientific_status,
        "scientific_summary": "physical doublets not yet validated",
        "utc": utc_now(), "libraries": list(libraries),
        "workload_generation_id": args.generation,
        "action": "NO_RESCORE_REANALYSIS",
        "stage_root": str(stage_root.resolve()),
        "job_id": clean(os.environ.get("SLURM_JOB_ID", "")) or
            "NOT_A_SLURM_JOB",
        "accessed_libraries": accessed_libraries,
        "protected_libraries_accessed": protected_accessed,
        "ledger_cells": len(cells),
        "candidate_rows": sum(score_rows_by_library.values()),
        "original_focal_cells": original_focal,
        "original_exact_agreements": original_agreement,
        "missing_molecule_channels": missing_molecule_channels,
        "warnings": len(warnings), "warning_types": dict(warning_counts),
        "expected_inputs_complete": expected_files_complete,
        "targeted_inputs_ready": targeted_inputs_ready,
        "required_artifacts_complete": required_artifacts_complete,
        "recommendation": recommendation,
        "preliminary_result_archive": str(preliminary_archive),
        "targeted_archive": str(targeted_archive),
        "targeted_accounting": targeted_accounting,
    }
    atomic_json(stage_root / "NO_RESCORE_ANALYSIS_FINISHED_UNVALIDATED", marker)
    print(json.dumps(marker, indent=2))
    return 0 if operational_status == "EXECUTION_COMPLETE" else 2


def write_marker(args):
    atomic_json(args.path, {
        "status": args.status, "detail": args.detail, "utc": utc_now(),
    })
    return 0


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=(
        "run", "targeted-gather", "validate-task",
        "finalize-controls", "compare-equivalence", "targeted-control",
        "targeted-reproject"))
    parser.add_argument("--stage-root", default="")
    parser.add_argument("--existing-stage-root", default="")
    parser.add_argument("--historical-gather-root", default="")
    parser.add_argument("--task-search-root", default="")
    parser.add_argument("--targeted-output-root", default="")
    parser.add_argument("--identity-tool", default="")
    parser.add_argument("--tool-bin-root", default="")
    parser.add_argument("--frozen-targets", default="")
    parser.add_argument("--targeted-input-manifest", default="")
    parser.add_argument("--frozen-comparisons", default="",
                        help="Authoritative frozen target/comparison TSV; never rematch")
    parser.add_argument("--libraries", nargs="*", default=[])
    parser.add_argument("--calibration-library", type=int,
                        default=CALIBRATION_LIBRARY)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--targeted-extraction-cpus", type=int, default=1)
    parser.add_argument("--targeted-extraction-memory", default="64G")
    parser.add_argument("--targeted-analysis-cpus", type=int, default=8)
    parser.add_argument("--targeted-analysis-memory", default="48G")
    parser.add_argument("--targeted-finalizer-cpus", type=int, default=1)
    parser.add_argument("--targeted-finalizer-memory", default="16G")
    parser.add_argument("--targeted-gather-cpus", type=int, default=2)
    parser.add_argument("--targeted-gather-memory", default="64G")
    parser.add_argument("--targeted-time", default="12:00:00")
    parser.add_argument("--targeted-partition", default="compute")
    parser.add_argument("--rna-error-ref", type=float, required=False,
                        default=0.001)
    parser.add_argument("--rna-error-alt", type=float, required=False,
                        default=0.001)
    parser.add_argument("--atac-error-ref", type=float, required=False,
                        default=0.005)
    parser.add_argument("--atac-error-alt", type=float, required=False,
                        default=0.005)
    parser.add_argument("--min-evidence", type=int, default=10)
    parser.add_argument("--max-second-fraction", type=float, default=0.95)
    parser.add_argument("--task-manifest", default="")
    parser.add_argument("--task-index", default="")
    parser.add_argument("--target-manifest", default="")
    parser.add_argument("--control-pairs", default="")
    parser.add_argument("--downsample-replicates", type=int,
                        default=DOWNSAMPLE_REPLICATES)
    parser.add_argument("--null-replicates", type=int,
                        default=CELL_NULL_REPLICATES)
    parser.add_argument("--targeted-root", default="")
    parser.add_argument("--task-kind", choices=("extraction", "analysis"),
                        default="analysis")
    parser.add_argument("--write-marker", action="store_true")
    parser.add_argument("--input-only", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--generation", default="")
    parser.add_argument("--scope", choices=("benchmark", "production"),
                        default="production")
    parser.add_argument("--operation", choices=("launch", "resume", "status"),
                        default="status")
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--left", default="")
    parser.add_argument("--right", default="")
    parser.add_argument("--legacy-reference", default="")
    parser.add_argument("--output", default="")
    parser.add_argument("--marker", default="")
    parser.add_argument("--path", default="")
    parser.add_argument("--status", default="")
    parser.add_argument("--detail", default="")
    return parser


def main():
    args = build_parser().parse_args()
    global SEED
    SEED = int(args.seed)
    try:
        if args.calibration_library != CALIBRATION_LIBRARY:
            raise RuntimeError(
                "this joint-doublet pilot requires calibration Library 25")
        if args.action == "run":
            if args.permutations != 10000:
                raise RuntimeError("no-rescore cross-assay null requires exactly 10000 permutations")
            return run_no_rescore(args)
        if args.action == "targeted-gather":
            return targeted_gather_v4(args)
        if args.action == "validate-task":
            return validate_task_action(args)
        if args.action == "finalize-controls":
            return finalize_controls(args)
        if args.action == "compare-equivalence":
            return compare_equivalence(args)
        if args.action == "targeted-control":
            return targeted_control(args)
        if args.action == "targeted-reproject":
            return targeted_reproject(args)
        raise RuntimeError(f"unsupported action: {args.action}")
    except (OSError, RuntimeError, ValueError, AssertionError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        if args.action == "targeted-control" and args.targeted_root:
            try:
                run_root, _ = _resolve_targeted_root(
                    args.targeted_root, args.scope)
                print(json.dumps({
                    "durably_recorded_submissions_before_failure":
                        _read_submission_ledger(run_root)
                }, indent=2), file=sys.stderr)
            except (OSError, RuntimeError, ValueError, json.JSONDecodeError):
                pass
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
