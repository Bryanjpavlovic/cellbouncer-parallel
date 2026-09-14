#!/usr/bin/env python3
"""Build the canonical per-library cell and ambient ledgers for arm-CNV calling."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

from tetra_arm_common import (
    AMBIENT_SCHEMA,
    CELL_MANIFEST_SCHEMA,
    RELEASE,
    canonical_barcode,
    canonical_pair,
    clean,
    file_record,
    finite_float,
    natural_key,
    open_text,
    parse_headerless_assignments,
    parse_libraries,
    parse_two_or_three_column,
    read_tsv,
    require_file,
    require_outputs_absent,
    truthy,
    write_json_atomic,
    write_tsv_atomic,
)


CELL_FIELDS = [
    "library", "barcode", "production_assignment", "donor_a", "donor_b",
    "donor_pair", "species_a", "species_b", "demux_score",
    "production_assignment_source", "application_state", "application_reason",
    "calibration_group", "uid", "current_ploidy_state",
    "ploidy_evidence_status", "nn_prob_tetraploid", "nn_qc_pass",
    "occupancy_state", "occupancy_evidence_status", "technical_state",
    "nuclear_reconciliation_status", "nuclear_warning_reasons",
    "species_evidence_status", "mitochondrial_evidence_status",
    "mitochondrial_resolution_status", "atac_evidence_status",
    "known_line_relationship", "uid_resolution_status", "metadata_event_status",
    "library_exchange_status",
    "review_required", "review_reasons", "downstream_release_status",
    "downstream_exclusion_reason", "identity_changed", "event_id",
    "event_identity_evidence_disposition", "event_review_scope",
    "cell_exception_reasons", "ambient_evaluation_status",
    "ambient_c", "ambient_c_se", "ambient_arm", "ambient_profile_mode",
    "ambient_profile_status", "ambient_production_minus_original_c",
    "ambient_exact_donor_burden_fields", "ambient_background_shift_fields",
    "model_eligible", "calibration_eligible",
    "eligibility_reasons", "schema_version",
]

AMBIENT_FIELDS = [
    "library", "barcode", "source_label", "scoring_profile_mass",
    "profile_origin", "ambient_arm", "schema_version",
]

FINAL_IDENTITY_SCHEMA = "identity_reconciliation_final_v2_phase3_dispositions"
FINAL_IDENTITY_REQUIRED_FIELDS = {
    "library", "barcode", "production_assignment",
    "production_assignment_source", "review_required",
    "downstream_release_status", "ambient_production_arm",
    "ambient_production_c", "ambient_evaluation_status", "event_id",
    "final_schema_version",
}


def canonical_identity(value: str) -> str:
    text = clean(value)
    pair = canonical_pair(text)
    if pair != ("", ""):
        return "+".join(pair)
    return text


def load_species_map(path: str) -> dict[str, str]:
    if not path:
        return {}
    rows = list(read_tsv(require_file(path, "panel metadata")))
    if not rows:
        return {}
    header = rows[0]
    donor_candidates = [
        "individual", "sample", "donor", "vcf_id", "VCF_ID", "line",
    ]
    species_candidates = ["species", "Species", "species_code"]
    donor_field = next((field for field in donor_candidates if field in header), None)
    species_field = next((field for field in species_candidates if field in header), None)
    if donor_field is None or species_field is None:
        return {}
    result = {}
    for row in rows:
        donor, species = clean(row.get(donor_field)), clean(row.get(species_field))
        if donor and species:
            if donor in result and result[donor] != species:
                raise ValueError(f"conflicting species metadata for donor {donor}")
            result[donor] = species
    return result


def load_groups(path: str, library: int) -> dict[str, str]:
    if not path:
        return {}
    group_path = require_file(path, "cell groups")
    with open_text(group_path) as handle:
        first_line = handle.readline().rstrip("\r\n")
    if not first_line:
        raise ValueError(f"empty cell-group file: {group_path}")
    first_fields = first_line.split("\t")
    normalized_header = {clean(field).lower() for field in first_fields}
    has_header = (
        "barcode" in normalized_header
        and bool({"group", "cluster", "cell_type"} & normalized_header)
    )

    if has_header:
        rows = read_tsv(group_path)
    else:
        def headerless_rows():
            with open_text(group_path) as handle:
                for line_number, line in enumerate(handle, start=1):
                    fields = line.rstrip("\r\n").split("\t")
                    if len(fields) != 2 or not all(clean(value) for value in fields):
                        raise ValueError(
                            f"headerless cell-group row must have two nonempty "
                            f"columns: {group_path}:{line_number}")
                    yield {"barcode": fields[0], "group": fields[1]}
        rows = headerless_rows()

    result = {}
    for row in rows:
        raw_library = clean(row.get("library"))
        if raw_library:
            raw_library = raw_library.lower().removeprefix("lib")
            try:
                if int(raw_library) != library:
                    continue
            except ValueError as exc:
                raise ValueError(f"invalid library in cell-group file: {raw_library}") from exc
        barcode = canonical_barcode(row.get("barcode", ""))
        group = clean(row.get("group") or row.get("cluster") or row.get("cell_type"))
        if not barcode or not group:
            continue
        if barcode in result and result[barcode] != group:
            raise ValueError(f"duplicate/conflicting group for {barcode}")
        result[barcode] = group
    return result


def load_ploidy_nn(path: str, library: int) -> dict[str, dict[str, str]]:
    """Load the optional production PLOIDY_NN table by canonical barcode."""
    if not path:
        return {}
    result = {}
    for row in read_tsv(require_file(path, "PLOIDY_NN calls")):
        if not {"barcode", "library", "prob_tetraploid", "qc_pass"} <= set(row):
            raise ValueError(f"PLOIDY_NN calls have an incompatible schema: {path}")
        raw_library = clean(row.get("library")).lower().removeprefix("lib")
        if raw_library:
            try:
                if int(raw_library) != library:
                    continue
            except ValueError as exc:
                raise ValueError(f"invalid library in PLOIDY_NN table: {raw_library}") from exc
        barcode = canonical_barcode(row.get("barcode", ""))
        if not barcode or barcode in result:
            raise ValueError(f"empty/duplicate PLOIDY_NN barcode: {path}")
        result[barcode] = row
    if not result:
        raise ValueError(f"PLOIDY_NN calls contain no lib{library} rows: {path}")
    return result


def load_source_profile(prefix: str, required: bool) -> dict:
    if not prefix:
        if required:
            raise ValueError("required ambient prefix was not supplied")
        return {}
    rate_path = prefix + ".contam_rate"
    prof_path = prefix + ".contam_prof"
    if not (os.path.isfile(rate_path) and os.path.getsize(rate_path) > 0 and
            os.path.isfile(prof_path) and os.path.getsize(prof_path) > 0):
        if required:
            raise FileNotFoundError(f"incomplete ambient prefix: {prefix}")
        return {}

    rates = parse_two_or_three_column(rate_path)
    global_profile = {}
    with open_text(prof_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\r\n").split()
            if len(fields) not in {2, 3}:
                raise ValueError(f"malformed ambient profile {prof_path}:{line_number}")
            label = fields[0]
            mass = finite_float(fields[1])
            if not math.isfinite(mass) or mass < 0:
                raise ValueError(f"invalid ambient mass {prof_path}:{line_number}")
            if label in global_profile:
                raise ValueError(f"duplicate ambient donor {label}: {prof_path}")
            global_profile[label] = mass
    total = sum(global_profile.values())
    if total <= 0:
        raise ValueError(f"ambient profile has zero total mass: {prof_path}")
    global_profile = {key: value / total for key, value in global_profile.items()}

    cell_profile_path = prefix + ".cell_source_profile.tsv"
    if not (os.path.isfile(cell_profile_path) and os.path.getsize(cell_profile_path) > 0):
        cell_profile_path = ""

    return {
        "prefix": os.path.abspath(prefix),
        "rates": rates,
        "global_profile": global_profile,
        # The long cell-by-source table can contain tens of millions of rows.
        # Keep only its path here; selected rows are validated and streamed below.
        "cell_profile_path": cell_profile_path,
        "mode": "CELL_SPECIFIC" if cell_profile_path else "GLOBAL_PROFILE",
        "records": {
            "contam_rate": file_record(rate_path),
            "contam_prof": file_record(prof_path),
            **({"cell_source_profile": file_record(cell_profile_path)}
               if cell_profile_path else {}),
        },
    }


def scan_cell_source_metadata(profile: dict, selected_barcodes: set[str]) -> dict[str, str]:
    """Validate selected simplexes without materializing their source vectors."""
    path = profile.get("cell_profile_path", "")
    if not path or not selected_barcodes:
        return {}
    required_fields = {
        "barcode", "source_label", "scoring_profile_mass",
        "scoring_profile_status",
    }
    sums = defaultdict(float)
    counts = Counter()
    statuses = {}
    for row in read_tsv(path):
        if not required_fields <= set(row):
            raise ValueError(
                f"cell source profile lacks {sorted(required_fields - set(row))}: {path}")
        barcode = canonical_barcode(row["barcode"])
        if barcode not in selected_barcodes:
            continue
        source = clean(row["source_label"])
        mass = finite_float(row["scoring_profile_mass"])
        if not source or not math.isfinite(mass) or mass < 0:
            raise ValueError(f"invalid cell source-profile row for {barcode}: {path}")
        sums[barcode] += mass
        counts[barcode] += 1
        status = clean(row["scoring_profile_status"]) or "UNKNOWN"
        previous = statuses.get(barcode)
        if previous is not None and previous != status:
            raise ValueError(f"inconsistent source-profile status for {barcode}: {path}")
        statuses[barcode] = status
    for barcode in statuses:
        if counts[barcode] <= 0 or abs(sums[barcode] - 1.0) > 1e-6:
            raise ValueError(
                f"cell-specific ambient simplex does not sum to one for {barcode}: "
                f"{sums[barcode]}")
    return statuses


def iter_ambient_rows(library: int, profile_objects: dict[str, dict],
                      profile_cells: dict[str, set[str]],
                      profile_statuses: dict[str, dict[str, str]],
                      selected_arms: dict[str, str]):
    """Stream exact selected per-cell source masses, then global fallbacks."""
    for prefix in sorted(profile_objects, key=natural_key):
        profile = profile_objects[prefix]
        selected = profile_cells[prefix]
        statuses = profile_statuses[prefix]
        path = profile.get("cell_profile_path", "")
        if path:
            for row in read_tsv(path):
                barcode = canonical_barcode(row.get("barcode", ""))
                if barcode not in selected or barcode not in statuses:
                    continue
                source = clean(row.get("source_label"))
                mass = finite_float(row.get("scoring_profile_mass"))
                yield {
                    "library": library,
                    "barcode": barcode,
                    "source_label": source,
                    "scoring_profile_mass": f"{mass:.17g}",
                    "profile_origin": prefix,
                    "ambient_arm": selected_arms[barcode],
                    "schema_version": AMBIENT_SCHEMA,
                }
        for barcode in sorted(selected - set(statuses), key=natural_key):
            for source, mass in sorted(
                    profile["global_profile"].items(), key=lambda item: natural_key(item[0])):
                yield {
                    "library": library,
                    "barcode": barcode,
                    "source_label": source,
                    "scoring_profile_mass": f"{mass:.17g}",
                    "profile_origin": prefix,
                    "ambient_arm": selected_arms[barcode],
                    "schema_version": AMBIENT_SCHEMA,
                }


def split_ledger(args) -> int:
    libraries = set(parse_libraries(args.libraries))
    input_path = require_file(args.input, "final reconciliation ledger")
    output_root = Path(os.path.abspath(args.output_root))
    output_paths = [
        *(output_root / f"lib{library}.final_cells.tsv.gz"
          for library in sorted(libraries)),
        output_root / "split_ledger_summary.tsv",
        output_root / "split_ledger_contract.json",
    ]
    require_outputs_absent(output_paths)
    output_root.mkdir(parents=True, exist_ok=True)
    grouped = {library: [] for library in libraries}
    header = None
    seen = set()
    with open_text(input_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames
        if not header or not FINAL_IDENTITY_REQUIRED_FIELDS <= set(header):
            raise ValueError("final reconciliation ledger has an incompatible schema")
        for line_number, row in enumerate(reader, start=2):
            raw_library = clean(row.get("library")).lower().removeprefix("lib")
            try:
                library = int(raw_library)
            except ValueError as exc:
                raise ValueError(f"invalid library at {input_path}:{line_number}") from exc
            if library not in libraries:
                continue
            if clean(row.get("final_schema_version")) != FINAL_IDENTITY_SCHEMA:
                raise ValueError(
                    f"unsupported final reconciliation schema at "
                    f"{input_path}:{line_number}: "
                    f"{clean(row.get('final_schema_version'))!r}")
            barcode = canonical_barcode(row.get("barcode", ""))
            key = (library, barcode)
            if not barcode or key in seen:
                raise ValueError(f"empty/duplicate cell key at {input_path}:{line_number}")
            seen.add(key)
            row["library"] = str(library)
            row["barcode"] = barcode
            grouped[library].append(row)

    summary_rows = []
    for library in sorted(libraries):
        rows = grouped[library]
        if not rows:
            raise ValueError(f"final reconciliation ledger has no rows for lib{library}")
        path = output_root / f"lib{library}.final_cells.tsv.gz"
        write_tsv_atomic(str(path), rows, header)
        summary_rows.append({
            "library": library,
            "cells": len(rows),
            "output": str(path),
            "schema_version": "tetra_arm_split_ledger_v1",
        })
    write_tsv_atomic(
        str(output_root / "split_ledger_summary.tsv"), summary_rows,
        ["library", "cells", "output", "schema_version"])
    write_json_atomic(str(output_root / "split_ledger_contract.json"), {
        "schema_version": "tetra_arm_split_ledger_contract_v1",
        "release": RELEASE,
        "input": file_record(input_path),
        "libraries": sorted(libraries),
        "cells": sum(len(value) for value in grouped.values()),
        "status": "PASS",
    })
    return 0


def row_bad_status(value: str) -> bool:
    upper = clean(value).upper()
    return any(token in upper for token in (
        "EXCLUDE", "UNRESOLVED", "MULTIPLET", "REVIEW", "FAIL", "BLOCK",
    ))


def has_event_id(value: str) -> bool:
    event = clean(value)
    return bool(event) and event.upper() not in {"NA", "N/A", "NONE", "."}


def prepare_library(args) -> int:
    library = int(args.library)
    output_dir = Path(os.path.abspath(args.output_dir))
    manifest_path = output_dir / f"lib{library}.cell_manifest.tsv.gz"
    ambient_path = output_dir / f"lib{library}.ambient_sources.tsv.gz"
    qc_path = output_dir / f"lib{library}.prepare_qc.tsv"
    contract_path = output_dir / f"lib{library}.prepare_contract.json"
    require_outputs_absent((manifest_path, ambient_path, qc_path, contract_path))
    ledger_path = require_file(args.ledger, "split final ledger")
    assignments_path = require_file(args.final_assignments, "final assignments")
    assignments = parse_headerless_assignments(assignments_path)
    groups = load_groups(args.cell_groups, library)
    ploidy_nn = load_ploidy_nn(args.ploidy_nn, library)
    species = load_species_map(args.panel_metadata)

    standard = load_source_profile(args.ambient_standard_prefix, required=True)
    arm_a = load_source_profile(args.ambient_arm_a_prefix, required=False)
    arm_c = load_source_profile(args.ambient_arm_c_prefix, required=False)

    ledger_rows = list(read_tsv(ledger_path))
    if not ledger_rows:
        raise ValueError(f"empty split ledger: {ledger_path}")
    has_identity_events = any(has_event_id(row.get("event_id"))
                              for row in ledger_rows)
    expected_production_arm = "C" if has_identity_events else "A"
    expected_ambient_statuses = (
        {"PAIRED_A_B_C_D", "PAIRED_A_B_C_D_NOT_APPLICABLE"}
        if has_identity_events else {"NOT_APPLICABLE_ZERO_EVENT"})
    # A standard original-assignment fit is equivalent to four-arm A only for
    # libraries with no reconciliation event.  In an event-bearing library,
    # silently substituting it would violate the finalized per-cell arm choice.
    arm_a_effective = arm_a or (standard if not has_identity_events else {})
    profiles = {"STANDARD": standard, "A": arm_a_effective, "C": arm_c}
    ledger_by_barcode = {}
    for row in ledger_rows:
        barcode = canonical_barcode(row.get("barcode", ""))
        if barcode in ledger_by_barcode:
            raise ValueError(f"duplicate ledger barcode {barcode}: {ledger_path}")
        ledger_by_barcode[barcode] = row
    if set(ledger_by_barcode) != set(assignments):
        missing = sorted(set(assignments) - set(ledger_by_barcode), key=natural_key)
        extra = sorted(set(ledger_by_barcode) - set(assignments), key=natural_key)
        raise ValueError(
            f"lib{library} final ledger/assignment barcode mismatch: "
            f"missing={missing[:5]} extra={extra[:5]}")

    cells = []
    reason_counts = Counter()
    profile_counts = Counter()
    profile_objects = {}
    profile_cells = defaultdict(set)
    selected_profiles = {}
    selected_arms = {}
    for barcode in sorted(assignments, key=natural_key):
        row = ledger_by_barcode[barcode]
        assignment, _assignment_type, score = assignments[barcode]
        ledger_assignment = clean(row.get("production_assignment"))
        if canonical_identity(assignment) != canonical_identity(ledger_assignment):
            raise ValueError(
                f"lib{library}/{barcode} final assignment mismatch: "
                f"{assignment!r} != {ledger_assignment!r}")
        donor_a, donor_b = canonical_pair(assignment)
        arm = clean(row.get("ambient_production_arm")).upper()
        if arm != expected_production_arm:
            raise ValueError(
                f"lib{library}/{barcode} final ambient production arm is "
                f"{arm!r}; expected uniform Arm {expected_production_arm} for "
                f"{'an event-bearing' if has_identity_events else 'a zero-event'} "
                "library")
        ambient_status = clean(row.get("ambient_evaluation_status")).upper()
        if ambient_status not in expected_ambient_statuses:
            raise ValueError(
                f"lib{library}/{barcode} has incompatible finalized ambient "
                f"status {ambient_status!r}; expected one of "
                f"{sorted(expected_ambient_statuses)}")
        selected = profiles.get(arm)
        if not selected:
            raise FileNotFoundError(
                f"lib{library}/{barcode} requires ambient Arm {arm}, but its bundle is absent")
        rate = selected["rates"].get(barcode)
        if rate is None:
            raise ValueError(
                f"ambient rate is absent for lib{library}/{barcode} in {selected['prefix']}")
        ambient_c, ambient_se = rate
        ledger_c = finite_float(row.get("ambient_production_c"))
        if has_identity_events and not math.isfinite(ledger_c):
            raise ValueError(
                f"final ledger has no finite ambient production rate for "
                f"lib{library}/{barcode}")
        if math.isfinite(ledger_c) and abs(ledger_c - ambient_c) > 1e-5:
            raise ValueError(
                f"ambient-rate mismatch for lib{library}/{barcode}: "
                f"ledger={ledger_c} prefix={ambient_c}")
        if not 0 <= ambient_c < 1:
            raise ValueError(f"ambient rate outside [0,1) for lib{library}/{barcode}")
        if math.isfinite(ambient_se) and ambient_se < 0:
            raise ValueError(
                f"ambient rate standard error is negative for lib{library}/{barcode}")
        prefix = selected["prefix"]
        profile_objects[prefix] = selected
        profile_cells[prefix].add(barcode)
        selected_profiles[barcode] = prefix
        selected_arms[barcode] = arm

        model_reasons = []
        if not donor_a or not donor_b:
            model_reasons.append("NOT_HETEROTYPIC_TWO_DONOR")
        if truthy(row.get("review_required")):
            model_reasons.append("IDENTITY_REVIEW_REQUIRED")
        if row_bad_status(row.get("technical_state")):
            model_reasons.append("TECHNICAL_STATE")
        if row_bad_status(row.get("occupancy_evidence_status")):
            model_reasons.append("OCCUPANCY_UNRESOLVED")
        # The canonical identity finalizer publishes exactly READY for cells
        # released to downstream analyses.  Treat a missing, novel, held, or
        # excluded value as ineligible rather than relying on a blacklist of
        # known failure words.  This also guarantees that a review-required
        # cell cannot enter the arm model if an inconsistent ledger ever pairs
        # it with an apparently benign status string.
        if clean(row.get("downstream_release_status")).upper() != "READY":
            model_reasons.append("DOWNSTREAM_NOT_RELEASED")
        model_eligible = not model_reasons

        calibration_reasons = list(model_reasons)
        ploidy_state = clean(row.get("current_ploidy_state"))
        if ploidy_state and "TETRA" not in ploidy_state.upper():
            calibration_reasons.append("NON_TETRAPLOID_STATE")
        if row_bad_status(row.get("ploidy_evidence_status")):
            calibration_reasons.append("PLOIDY_EVIDENCE_UNRESOLVED")
        if row_bad_status(row.get("nuclear_reconciliation_status")):
            calibration_reasons.append("NUCLEAR_IDENTITY_UNRESOLVED")
        if args.ploidy_nn and barcode not in ploidy_nn:
            raise ValueError(
                f"PLOIDY_NN calls are missing finalized lib{library} barcode "
                f"{barcode}; refusing ledger fallback")
        nn_row = ploidy_nn.get(barcode, {})
        if args.ploidy_nn:
            p_tet = finite_float(nn_row.get("prob_tetraploid"))
            nn_qc = clean(nn_row.get("qc_pass"))
            if not math.isfinite(p_tet) or not 0.0 <= p_tet <= 1.0:
                raise ValueError(
                    f"invalid PLOIDY_NN probability for lib{library}/{barcode}")
            if nn_qc not in {"0", "1"}:
                raise ValueError(
                    f"invalid PLOIDY_NN qc_pass for lib{library}/{barcode}")
        else:
            p_tet = finite_float(
                row.get("nn_prob_tetraploid") or row.get("prob_tetraploid"))
            nn_qc = clean(row.get("nn_qc_pass"))
        if math.isfinite(p_tet) and p_tet < args.min_calibration_tet_probability:
            calibration_reasons.append("LOW_TETRAPLOID_PROBABILITY")
        if nn_qc and not truthy(nn_qc):
            calibration_reasons.append("PLOIDY_NN_QC_FAIL")
        if has_event_id(row.get("event_id")):
            calibration_reasons.append("IDENTITY_EVENT_CELL")
        calibration_eligible = model_eligible and not calibration_reasons
        reasons = sorted(set(model_reasons + calibration_reasons), key=natural_key)
        if not reasons:
            reasons = ["PASS"]
        reason_counts.update(reasons)
        group = groups.get(barcode, f"lib{library}")
        uid = clean(row.get("uid_or_uid_set") or row.get("reconciled_uid") or row.get("uid"))
        cells.append({
            "library": library,
            "barcode": barcode,
            "production_assignment": assignment,
            "donor_a": donor_a,
            "donor_b": donor_b,
            "donor_pair": f"{donor_a}+{donor_b}" if donor_a else "NA",
            "species_a": species.get(donor_a, "NA") if donor_a else "NA",
            "species_b": species.get(donor_b, "NA") if donor_b else "NA",
            "demux_score": "NA" if not math.isfinite(score) else f"{score:.17g}",
            "production_assignment_source": clean(
                row.get("production_assignment_source")) or "NA",
            "application_state": clean(row.get("application_state")) or "NA",
            "application_reason": clean(row.get("application_reason")) or "NA",
            "calibration_group": group,
            "uid": uid or "NA",
            "current_ploidy_state": clean(row.get("current_ploidy_state")) or "NA",
            "ploidy_evidence_status": clean(row.get("ploidy_evidence_status")) or "NA",
            "nn_prob_tetraploid": "NA" if not math.isfinite(p_tet) else f"{p_tet:.17g}",
            "nn_qc_pass": nn_qc or "NA",
            "occupancy_state": clean(row.get("occupancy_state")) or "NA",
            "occupancy_evidence_status": clean(row.get("occupancy_evidence_status")) or "NA",
            "technical_state": clean(row.get("technical_state")) or "NA",
            "nuclear_reconciliation_status": clean(row.get("nuclear_reconciliation_status")) or "NA",
            "nuclear_warning_reasons": clean(row.get("nuclear_warning_reasons")) or "NONE",
            "species_evidence_status": clean(row.get("species_evidence_status")) or "NA",
            "mitochondrial_evidence_status": clean(row.get("mitochondrial_evidence_status")) or "NA",
            "mitochondrial_resolution_status": clean(
                row.get("mitochondrial_resolution_status")) or "NA",
            "atac_evidence_status": clean(row.get("atac_evidence_status")) or "NA",
            "known_line_relationship": clean(row.get("known_line_relationship")) or "NA",
            "uid_resolution_status": clean(row.get("uid_resolution_status")) or "NA",
            "metadata_event_status": clean(row.get("metadata_event_status")) or "NA",
            "library_exchange_status": clean(row.get("library_exchange_status")) or "NA",
            "review_required": clean(row.get("review_required")) or "FALSE",
            "review_reasons": clean(row.get("review_reasons")) or "NONE",
            "downstream_release_status": clean(row.get("downstream_release_status")) or "NA",
            "downstream_exclusion_reason": clean(row.get("downstream_exclusion_reason")) or "NONE",
            "identity_changed": "TRUE" if canonical_identity(
                row.get("demux_original_assignment", "")) != canonical_identity(assignment) else "FALSE",
            "event_id": clean(row.get("event_id")) or "NA",
            "event_identity_evidence_disposition": clean(
                row.get("event_identity_evidence_disposition")) or "NA",
            "event_review_scope": clean(row.get("event_review_scope")) or "NA",
            "cell_exception_reasons": clean(row.get("cell_exception_reasons")) or "NONE",
            "ambient_evaluation_status": clean(row.get("ambient_evaluation_status")) or "NA",
            "ambient_c": f"{ambient_c:.17g}",
            "ambient_c_se": "NA" if not math.isfinite(ambient_se) else f"{ambient_se:.17g}",
            "ambient_arm": arm,
            "ambient_profile_mode": "PENDING",
            "ambient_profile_status": "PENDING",
            "ambient_production_minus_original_c": clean(
                row.get("ambient_production_minus_original_c")) or "NA",
            "ambient_exact_donor_burden_fields": clean(
                row.get("ambient_exact_donor_burden_fields")) or "NA",
            "ambient_background_shift_fields": clean(
                row.get("ambient_background_shift_fields")) or "NA",
            "model_eligible": int(model_eligible),
            "calibration_eligible": int(calibration_eligible),
            "eligibility_reasons": ";".join(reasons),
            "schema_version": CELL_MANIFEST_SCHEMA,
        })

    profile_statuses = {
        prefix: scan_cell_source_metadata(profile_objects[prefix], set(barcodes))
        for prefix, barcodes in profile_cells.items()
    }
    for cell in cells:
        barcode = cell["barcode"]
        prefix = selected_profiles[barcode]
        status = profile_statuses[prefix].get(barcode)
        if status is None:
            cell["ambient_profile_mode"] = "GLOBAL_PROFILE_FALLBACK"
            cell["ambient_profile_status"] = "GLOBAL_PROFILE_FALLBACK"
        else:
            cell["ambient_profile_mode"] = "CELL_SPECIFIC"
            cell["ambient_profile_status"] = status
        profile_counts[cell["ambient_profile_mode"]] += 1

    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv_atomic(str(manifest_path), cells, CELL_FIELDS)
    ambient_rows = write_tsv_atomic(
        str(ambient_path),
        iter_ambient_rows(
            library, profile_objects, profile_cells, profile_statuses, selected_arms),
        AMBIENT_FIELDS)
    write_tsv_atomic(
        str(qc_path),
        [{
            "library": library,
            "cells": len(cells),
            "model_eligible": sum(int(row["model_eligible"]) for row in cells),
            "calibration_eligible": sum(int(row["calibration_eligible"]) for row in cells),
            "ambient_rows": ambient_rows,
            "profile_modes": ";".join(f"{key}={value}" for key, value in sorted(profile_counts.items())),
            "eligibility_reasons": ";".join(f"{key}={value}" for key, value in sorted(reason_counts.items())),
            "status": "PASS",
            "schema_version": "tetra_arm_prepare_qc_v1",
        }],
        ["library", "cells", "model_eligible", "calibration_eligible",
         "ambient_rows", "profile_modes", "eligibility_reasons", "status",
         "schema_version"])
    write_json_atomic(str(contract_path), {
        "schema_version": "tetra_arm_prepare_contract_v1",
        "release": RELEASE,
        "library": library,
        "inputs": {
            "ledger": file_record(ledger_path),
            "final_assignments": file_record(assignments_path),
            "panel_metadata": file_record(args.panel_metadata) if args.panel_metadata else None,
            "cell_groups": file_record(args.cell_groups) if args.cell_groups else None,
            "ploidy_nn": file_record(args.ploidy_nn) if args.ploidy_nn else None,
        },
        "ambient_sources": {
            key: value.get("records") if value else None
            for key, value in profiles.items()
        },
        "ambient_arm_a_standard_equivalent": bool(
            not arm_a and not has_identity_events),
        "outputs": {
            "cell_manifest": str(manifest_path),
            "ambient_sources": str(ambient_path),
        },
        "cells": len(cells),
        "status": "PASS",
    })
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare canonical cells and exact ambient donor mixtures for arm-CNV analysis.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {RELEASE}")
    sub = parser.add_subparsers(dest="command", required=True)

    split = sub.add_parser("split-ledger")
    split.add_argument("--input", required=True)
    split.add_argument("--libraries", nargs="+", required=True)
    split.add_argument("--output-root", required=True)
    split.set_defaults(func=split_ledger)

    prepare = sub.add_parser("prepare-library")
    prepare.add_argument("--library", type=int, required=True)
    prepare.add_argument("--ledger", required=True)
    prepare.add_argument("--final-assignments", required=True)
    prepare.add_argument("--ambient-standard-prefix", required=True)
    prepare.add_argument("--ambient-arm-a-prefix", default="")
    prepare.add_argument("--ambient-arm-c-prefix", default="")
    prepare.add_argument("--panel-metadata", default="")
    prepare.add_argument("--cell-groups", default="")
    prepare.add_argument("--ploidy-nn", default="")
    prepare.add_argument("--min-calibration-tet-probability", type=float, default=0.90)
    prepare.add_argument("--output-dir", required=True)
    prepare.set_defaults(func=prepare_library)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        return args.func(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
