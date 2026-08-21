#!/usr/bin/env python3
"""Validate structural invariants of an identity-reconciliation run."""
from __future__ import annotations

import argparse
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List

from identity_reconciliation_common import (
    DEFAULT_IDENTITY_POLICY, POLICY_VERSION, SCHEMA_VERSION, canonical_genotype, clean, donor_components, parse_library_spec, read_assignments,
    read_tsv, write_tsv,
)

SUMMARY_FIELDS = ["check", "status", "n_failures", "detail"]
FAIL_FIELDS = ["check", "library", "barcode", "detail"]


def semicolon_set(value: str):
    return {clean(x) for x in clean(value).split(";") if clean(x)}


def numeric_float(value):
    try:
        x = float(clean(value))
        return x if math.isfinite(x) else math.nan
    except Exception:
        return math.nan


def numeric_int(value):
    try:
        return int(float(clean(value)))
    except Exception:
        return 0


def enum_value(value):
    return "" if value is None else str(value).strip()


def expected_context_from_metadata(metadata):
    explicit = {}
    components = {}
    composites = {}
    for lib, by_genotype in metadata.items():
        explicit[lib] = set()
        components[lib] = set()
        composites[lib] = defaultdict(set)
        for genotype in by_genotype:
            comps = donor_components(genotype)
            if not comps:
                continue
            components[lib].update(comps)
            if len(comps) == 1:
                explicit[lib].add(comps[0])
            else:
                for donor in set(comps):
                    composites[lib][donor].add(genotype)
    return explicit, components, composites


def expected_singlet_relationship(lib, genotype, explicit, components, composites, global_donors):
    comps = donor_components(genotype)
    if len(comps) != 1:
        return "NOT_SINGLET", set()
    donor = comps[0]
    if donor not in global_donors:
        return "UNRESOLVED", set()
    if donor in explicit.get(lib, set()):
        return "EXPECTED_SINGLET", set()
    if donor in components.get(lib, set()):
        return "EXPECTED_COMPOSITE_COMPONENT_SINGLET", set(composites.get(lib, {}).get(donor, set()))
    return "UNEXPECTED_SINGLET", set()


def expected_pair_roster(library_a, library_b, components):
    a = set(components.get(library_a, set()))
    b = set(components.get(library_b, set()))
    shared = a & b
    a_specific = a - b
    b_specific = b - a
    total = len(a_specific) + len(b_specific)
    if not a_specific and not b_specific:
        relation = "ROSTER_EQUIVALENT_NONDISCRIMINATING"
        discriminability = "NONE"
    else:
        if not shared:
            relation = "DISTINCT_DONOR_ROSTERS"
        elif not a_specific or not b_specific:
            relation = "NESTED_DONOR_ROSTERS"
        else:
            relation = "PARTIALLY_OVERLAPPING_DONOR_ROSTERS"
        discriminability = "WEAK" if total <= 2 else ("MODERATE" if total <= 5 else "STRONG")
    return relation, discriminability, shared, a_specific, b_specific


def expected_directional_strength(n_detected, n_diagnostic):
    if n_diagnostic <= 0 or n_detected <= 0:
        return "NONE"
    coverage = n_detected / n_diagnostic
    if n_detected == 1:
        return "WEAK"
    if n_detected >= 3 and coverage >= 0.75:
        return "STRONG"
    return "MODERATE"


def expected_exchange_classification(
    roster_relation, pair_discriminability, a_to_b_strength, b_to_a_strength,
    a_native_retention, b_native_retention,
):
    if roster_relation == "ROSTER_EQUIVALENT_NONDISCRIMINATING":
        return "NONDISCRIMINATING", "ROSTER_EQUIVALENT_NONDISCRIMINATING", "NONDISCRIMINATING"
    a_to_b = a_to_b_strength != "NONE"
    b_to_a = b_to_a_strength != "NONE"
    if not a_to_b and not b_to_a:
        return "NONE", "NO_LIBRARY_EXCHANGE_SIGNAL", "NONE"
    if a_to_b and b_to_a:
        a_disp = 1.0 - a_native_retention if math.isfinite(a_native_retention) else 0.0
        b_disp = 1.0 - b_native_retention if math.isfinite(b_native_retention) else 0.0
        strong_reciprocal_signature = (
            pair_discriminability in {"MODERATE", "STRONG"}
            and "STRONG" in {a_to_b_strength, b_to_a_strength}
        )
        interpretation = (
            "LIKELY_RECIPROCAL_LIBRARY_SWAP"
            if strong_reciprocal_signature and a_disp >= 0.75 and b_disp >= 0.75
            else "RECIPROCAL_LIBRARY_MIXING"
        )
        confidence = (
            "SUGGESTIVE" if pair_discriminability == "WEAK"
            else "STRONG" if interpretation == "LIKELY_RECIPROCAL_LIBRARY_SWAP" and pair_discriminability == "STRONG"
            else "STRONG" if a_to_b_strength == "STRONG" and b_to_a_strength == "STRONG"
            else "MODERATE"
        )
        return "RECIPROCAL", interpretation, confidence
    if a_to_b:
        reciprocal = "ONE_WAY_A_TO_B"
        strength = a_to_b_strength
        target_retention = b_native_retention
    else:
        reciprocal = "ONE_WAY_B_TO_A"
        strength = b_to_a_strength
        target_retention = a_native_retention
    target_displacement = 1.0 - target_retention if math.isfinite(target_retention) else 0.0
    if strength in {"MODERATE", "STRONG"} and target_displacement >= 0.50:
        interpretation = "LIKELY_PARTIAL_LIBRARY_REPLACEMENT"
    elif strength in {"MODERATE", "STRONG"} and math.isfinite(target_retention) and target_retention >= 0.75:
        interpretation = "LIKELY_CROSS_CONTAMINATION"
    else:
        interpretation = "ONE_WAY_FOREIGN_SIGNATURE"
    confidence = "SUGGESTIVE" if pair_discriminability == "WEAK" or strength == "WEAK" else ("STRONG" if strength == "STRONG" else "MODERATE")
    return reciprocal, interpretation, confidence


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--libraries", nargs="+", required=True)
    p.add_argument("--demux-root", required=True)
    p.add_argument("--metadata-root", required=True)
    p.add_argument("--candidate-root", required=True)
    p.add_argument("--doublet-context-root", default="")
    p.add_argument("--nuclear-root", required=True)
    p.add_argument("--mt-root", required=True)
    p.add_argument("--atac-root", required=True)
    p.add_argument("--decisions-root", required=True)
    p.add_argument("--reports-root", required=True)
    p.add_argument("--evidence-mode", choices=("rna", "rna-atac"), default="rna")
    p.add_argument("--output-root", required=True)
    return p.parse_args()


def main():
    args = parse_args(); libs = parse_library_spec(args.libraries)
    failures: List[dict] = []; summaries = []

    def fail(check, lib="", bc="", detail=""):
        failures.append({"check": check, "library": lib, "barcode": bc, "detail": detail})

    def finish(check, before, detail=""):
        n = len(failures) - before
        summaries.append({"check": check, "status": "PASS" if n == 0 else "FAIL", "n_failures": n, "detail": detail})

    all_path = os.path.join(args.decisions_root, "all_libraries.reconciled_cells.tsv.gz")
    all_rows = read_tsv(all_path)
    by_key = defaultdict(list)
    for row in all_rows: by_key[(clean(row.get("library")), clean(row.get("barcode")))].append(row)

    b = len(failures)
    for key, rows in by_key.items():
        if len(rows) != 1: fail("unique_library_barcode", key[0], key[1], f"rows={len(rows)}")
    finish("unique_library_barcode", b)

    metadata = defaultdict(dict)
    for row in read_tsv(os.path.join(args.metadata_root, "library_expected_genotypes.tsv")):
        metadata[clean(row.get("library"))][canonical_genotype(row.get("canonical_genotype", ""))] = row
    global_line_rows = {
        canonical_genotype(r.get("canonical_genotype", "")): r
        for r in read_tsv(os.path.join(args.metadata_root, "global_biological_lines.tsv"))
        if canonical_genotype(r.get("canonical_genotype", ""))
    }
    global_lines = set(global_line_rows)
    global_donors = {
        clean(r.get("donor_id"))
        for r in read_tsv(os.path.join(args.metadata_root, "global_donors.tsv"))
        if clean(r.get("donor_id"))
    }
    explicit_expected_singlets, expected_donor_components, expected_composites = expected_context_from_metadata(metadata)
    members = defaultdict(set)
    for row in read_tsv(os.path.join(args.metadata_root, "library_uid_members.tsv")):
        members[(clean(row.get("library")), canonical_genotype(row.get("canonical_genotype", "")))].add(clean(row.get("uid")))

    for n in libs:
        lib = f"lib{n}"; prefix = os.path.join(args.demux_root, f"Tet_2025_Multiome-RNA_{n}", "demux_nomito", f"{lib}_demuxed")
        original = read_assignments(prefix + ".assignments")
        lib_rows = read_tsv(os.path.join(args.decisions_root, f"{lib}.reconciled_cells.tsv.gz"))
        lib_map = {clean(r.get("barcode")): r for r in lib_rows}

        b = len(failures)
        if set(original) != set(lib_map):
            for bc in sorted(set(original) ^ set(lib_map))[:100]: fail("all_input_barcodes_once", lib, bc, "barcode set mismatch")
        finish("all_input_barcodes_once", b, lib)

        b = len(failures)
        for bc, orig in original.items():
            row = lib_map.get(bc)
            if row and canonical_genotype(row.get("original_demux_assignment", "")) != canonical_genotype(orig["assignment"]):
                fail("original_assignments_preserved", lib, bc, f"source={orig['assignment']} result={row.get('original_demux_assignment')}")
        finish("original_assignments_preserved", b, lib)

        b = len(failures)
        candidate_rows = read_tsv(os.path.join(args.candidate_root, f"{lib}.identity_candidates.tsv.gz"))
        candidate_set = {(clean(r.get("barcode")), canonical_genotype(r.get("donor_genotype", ""))) for r in candidate_rows}

        b_universe = len(failures)
        for cand in candidate_rows:
            g = canonical_genotype(cand.get("donor_genotype", ""))
            current = canonical_genotype(cand.get("current_donor_genotype", ""))
            if not g or g == current:
                continue
            comps = donor_components(g)
            admiss = clean(cand.get("biological_admissibility"))
            if len(comps) == 1:
                if comps[0] not in global_donors:
                    fail("candidate_global_biological_universe", lib, clean(cand.get("barcode")), f"unknown donor {g}")
                if admiss != "SINGLET_IDENTITY_CANDIDATE":
                    fail("candidate_admissibility", lib, clean(cand.get("barcode")), f"{g}: {admiss}")
            elif len(comps) == 2:
                if g not in global_lines:
                    fail("candidate_global_biological_universe", lib, clean(cand.get("barcode")), f"unmade line {g}")
                if admiss != "BIOLOGICAL_SINGLE_CELL_ALLOWED":
                    fail("candidate_admissibility", lib, clean(cand.get("barcode")), f"{g}: {admiss}")
            else:
                fail("candidate_global_biological_universe", lib, clean(cand.get("barcode")), f"non-single-cell candidate {g}")
        finish("candidate_global_biological_universe", b_universe, lib)

        b = len(failures)
        tech_path = os.path.join(args.candidate_root, f"{lib}.technical_multiplet_candidates.tsv.gz")
        tech_rows = read_tsv(tech_path) if os.path.isfile(tech_path) else []
        for tr in tech_rows:
            comp = canonical_genotype(tr.get("donor_composition", ""))
            comps = donor_components(comp)
            tclass = clean(tr.get("technical_class"))
            if not comps or any(x not in global_donors for x in comps):
                fail("technical_multiplet_known_donors", lib, clean(tr.get("barcode")), comp)
            if tclass == "UNMADE_DONOR_PAIR":
                if len(comps) != 2 or comp in global_lines:
                    fail("technical_multiplet_unmade_pair", lib, clean(tr.get("barcode")), comp)
            if tclass == "KNOWN_TETRAPLOID_PLUS_EXTRA_DONOR":
                current = canonical_genotype(tr.get("current_donor_genotype", ""))
                extra = clean(tr.get("additional_donor"))
                if current not in global_lines or len(donor_components(current)) != 2 or extra not in global_donors or extra in donor_components(current):
                    fail("technical_multiplet_tetraploid_plus_donor", lib, clean(tr.get("barcode")), f"current={current} extra={extra}")
        finish("technical_multiplet_candidate_invariants", b, lib)

        b = len(failures)
        for row in lib_rows:
            applied = clean(row.get("reassignment_applied")).upper() == "TRUE"
            if applied and canonical_genotype(row.get("reconciled_donor_genotype", "")) == canonical_genotype(row.get("original_demux_assignment", "")) and clean(row.get("final_action")) != "RECLASSIFY_PLOIDY":
                fail("applied_change_changes_state", lib, clean(row.get("barcode")), "applied flag without genotype/ploidy change")
            if applied and clean(row.get("final_action")) == "REASSIGN_GENOTYPE":
                key = (clean(row.get("barcode")), canonical_genotype(row.get("reconciled_donor_genotype", "")))
                if key not in candidate_set: fail("changed_genotype_has_candidate", lib, key[0], key[1])
                if clean(row.get("decision_confidence")) != "DECISIVE": fail("only_decisive_autoapplied", lib, key[0], clean(row.get("decision_confidence")))
                proposed_for_threshold = canonical_genotype(row.get("proposed_donor_genotype", ""))
                proposed_components = donor_components(proposed_for_threshold)
                if len(proposed_components) == 1:
                    unexpected_here = clean(row.get("singlet_library_relationship")) == "UNEXPECTED_SINGLET"
                else:
                    unexpected_here = clean(row.get("proposed_library_expected_status")) != "EXPECTED"
                if unexpected_here:
                    line_mass = int(float(clean(row.get("alternative_line_event_mass")) or 0))
                    min_autoapply = int(DEFAULT_IDENTITY_POLICY.get("unexpected_line_autoapply_min_cells", 100))
                    if line_mass < min_autoapply:
                        fail("applied_unexpected_line_autoapply_mass_threshold", lib, key[0], f"{key[1]} mass={line_mass} < {min_autoapply}")

            current = canonical_genotype(row.get("current_refined_assignment", "")) or canonical_genotype(row.get("original_demux_assignment", ""))
            reconciled = canonical_genotype(row.get("reconciled_donor_genotype", ""))
            if not applied and reconciled != current:
                fail("nonapplied_preserves_current_state", lib, clean(row.get("barcode")), f"current={current} reconciled={reconciled}")

            proposed = canonical_genotype(row.get("proposed_donor_genotype", ""))
            if proposed and proposed != current and clean(row.get("proposed_droplet_state")) == "SINGLE_CELL":
                comps = donor_components(proposed)
                if len(comps) == 1 and comps[0] not in global_donors:
                    fail("proposed_state_global_biological_universe", lib, clean(row.get("barcode")), f"unknown donor {proposed}")
                elif len(comps) == 2 and proposed not in global_lines:
                    fail("proposed_state_global_biological_universe", lib, clean(row.get("barcode")), f"unmade line {proposed}")
                elif len(comps) not in {1, 2}:
                    fail("proposed_state_global_biological_universe", lib, clean(row.get("barcode")), f"invalid biological single-cell state {proposed}")
        finish("applied_change_invariants", b, lib)

        b = len(failures)
        for row in lib_rows:
            bc = clean(row.get("barcode"))
            proposed = canonical_genotype(row.get("proposed_donor_genotype", ""))
            expected_relationship, expected_context = expected_singlet_relationship(
                lib, proposed, explicit_expected_singlets, expected_donor_components,
                expected_composites, global_donors,
            )
            observed_relationship = clean(row.get("singlet_library_relationship"))
            if observed_relationship != expected_relationship:
                fail(
                    "singlet_library_relationship",
                    lib, bc,
                    f"{proposed}: observed={observed_relationship} expected={expected_relationship}",
                )
            observed_context = semicolon_set(row.get("expected_composite_context", ""))
            if expected_relationship == "EXPECTED_COMPOSITE_COMPONENT_SINGLET":
                if not expected_context:
                    fail("component_singlet_has_expected_composite_context", lib, bc, proposed)
                if observed_context != expected_context:
                    fail(
                        "component_singlet_expected_composite_context",
                        lib, bc,
                        f"observed={sorted(observed_context)} expected={sorted(expected_context)}",
                    )
                if clean(row.get("library_exchange_evidence_eligible")).upper() == "TRUE":
                    fail("component_singlet_excluded_from_exchange_evidence", lib, bc, proposed)
            elif observed_context:
                fail("noncomponent_singlet_has_no_composite_context", lib, bc, f"{proposed}: {sorted(observed_context)}")
            if expected_relationship == "EXPECTED_SINGLET" and observed_relationship == "UNEXPECTED_SINGLET":
                fail("explicit_expected_singlet_not_unexpected", lib, bc, proposed)
        finish("component_singlet_classification_invariants", b, lib)

        # Independently reconstruct the library-local diploid donor population
        # used by the reconciler's occupancy guard: expected singlet lines,
        # robust singlet events, and robust residual component-singlet populations.
        # This protects
        # against a future regression that labels a repeated A+B donor signal as
        # an intact biological line when D[A]+D[B] remains an equally plausible
        # technical-doublet explanation.
        b = len(failures)
        min_event_cells = int(DEFAULT_IDENTITY_POLICY.get("event_min_cells", 10))
        local_diploid_donors = {
            donor_components(g)[0]
            for g in metadata.get(lib, {})
            if len(donor_components(g)) == 1
        }
        current_component_counts = Counter()
        for row in lib_rows:
            current = canonical_genotype(row.get("current_refined_assignment", "")) or canonical_genotype(row.get("original_demux_assignment", ""))
            current_ploidy = clean(row.get("current_ploidy")).upper()
            if clean(row.get("explicit_multiplet_evidence")).upper() == "TRUE" or current_ploidy in {"T", "TET", "TETRAPLOID", "HOMOTYPIC"}:
                continue
            relationship, _ = expected_singlet_relationship(
                lib, current, explicit_expected_singlets, expected_donor_components,
                expected_composites, global_donors,
            )
            if relationship == "EXPECTED_COMPOSITE_COMPONENT_SINGLET":
                current_component_counts[current] += 1
        for donor, count in current_component_counts.items():
            if count >= min_event_cells:
                local_diploid_donors.add(donor)
        for row in lib_rows:
            proposed = canonical_genotype(row.get("proposed_donor_genotype", ""))
            comps = donor_components(proposed)
            line_mass = int(float(clean(row.get("alternative_line_event_mass")) or 0))
            if len(comps) == 1 and line_mass >= min_event_cells:
                local_diploid_donors.add(comps[0])

        for row in lib_rows:
            bc = clean(row.get("barcode"))
            proposed = canonical_genotype(row.get("proposed_donor_genotype", ""))
            current = canonical_genotype(row.get("current_refined_assignment", "")) or canonical_genotype(row.get("original_demux_assignment", ""))
            comps = donor_components(proposed)
            status = clean(row.get("occupancy_resolution_status"))
            action = clean(row.get("final_action"))
            applied = clean(row.get("reassignment_applied")).upper() == "TRUE"
            expected_here = clean(row.get("proposed_library_expected_status")) == "EXPECTED"
            droplet_flag = clean(row.get("current_droplet_flag"))
            expected_explicit_multiplet = (
                droplet_flag.lower() not in {"", "none", "single", "single_cell", "0", "no"}
            )
            observed_explicit_multiplet = (
                clean(row.get("explicit_multiplet_evidence")).upper() == "TRUE"
            )
            if observed_explicit_multiplet != expected_explicit_multiplet:
                fail(
                    "explicit_multiplet_evidence_definition",
                    lib,
                    bc,
                    f"droplet_flag={droplet_flag or 'NA'} explicit={observed_explicit_multiplet}",
                )

            if proposed != current and len(comps) == 2 and len(set(comps)) == 2 and all(d in local_diploid_donors for d in comps):
                if status != "DONOR_PAIR_CELLULAR_ORIGIN_UNRESOLVED":
                    fail("heterotypic_pair_occupancy_status", lib, bc, f"{proposed}: {status}")
                if applied or action == "REASSIGN_GENOTYPE":
                    fail("heterotypic_pair_occupancy_blocks_reassignment", lib, bc, f"{proposed}: action={action} applied={applied}")

            if (
                len(comps) == 2
                and len(set(comps)) == 2
                and proposed != current
                and expected_explicit_multiplet
                and (applied or action == "REASSIGN_GENOTYPE")
            ):
                fail("explicit_multiplet_evidence_blocks_heterotypic_reassignment", lib, bc, f"{proposed}: action={action} applied={applied}")

            if (
                len(comps) == 2
                and len(set(comps)) == 1
                and proposed != current
                and not expected_here
                and comps[0] in local_diploid_donors
            ):
                if status != "HOMOTET_VS_SAME_DONOR_DOUBLET_UNRESOLVED":
                    fail("unexpected_homotet_occupancy_status", lib, bc, f"{proposed}: {status}")
                if applied or action == "RECLASSIFY_PLOIDY":
                    fail("unexpected_homotet_occupancy_blocks_reclassification", lib, bc, f"{proposed}: action={action} applied={applied}")

            if (
                len(comps) == 2
                and len(set(comps)) == 1
                and proposed != current
                and expected_explicit_multiplet
                and (applied or action == "RECLASSIFY_PLOIDY")
            ):
                fail("explicit_multiplet_evidence_blocks_homotet_reclassification", lib, bc, f"{proposed}: action={action} applied={applied}")
        finish("occupancy_ambiguity_invariants", b, lib)

        b = len(failures)
        multiplet_rows = [r for r in lib_rows if clean(r.get("reconciled_droplet_state")) == "TECHNICAL_MULTIPLET"]
        for row in multiplet_rows:
            g = canonical_genotype(row.get("reconciled_donor_genotype", ""))
            if clean(row.get("original_demux_type")).upper() != "D" or len(donor_components(g)) != 2 or g in global_lines:
                fail("genotype_identifiable_technical_multiplet", lib, clean(row.get("barcode")), f"type={row.get('original_demux_type')} genotype={g}")
        simple_path = os.path.join(args.decisions_root, f"{lib}.reconciled_single_cell.assignments")
        if os.path.isfile(simple_path):
            simple = read_assignments(simple_path)
            multiplets = {clean(r.get("barcode")) for r in multiplet_rows}
            for bc in multiplets & set(simple): fail("multiplets_excluded_from_single_cell_assignments", lib, bc, "technical multiplet present in compatibility assignments")
        finish("multiplets_excluded_from_single_cell_assignments", b, lib)

        b = len(failures)
        if os.path.isfile(simple_path):
            simple = read_assignments(simple_path)
            expected_single = {
                clean(r.get("barcode")): r
                for r in lib_rows
                if clean(r.get("reconciled_droplet_state")) == "SINGLE_CELL"
                and canonical_genotype(r.get("reconciled_donor_genotype", ""))
            }
            if set(simple) != set(expected_single):
                for bc in sorted(set(simple) ^ set(expected_single))[:100]:
                    fail("reconciled_assignments_barcode_set", lib, bc, "assignment/decision barcode set mismatch")
            for bc in sorted(set(simple) & set(expected_single)):
                outrow = simple[bc]; decision = expected_single[bc]
                got_genotype = canonical_genotype(outrow.get("assignment", ""))
                expected_genotype = canonical_genotype(decision.get("reconciled_donor_genotype", ""))
                if got_genotype != expected_genotype:
                    fail("reconciled_assignments_genotype", lib, bc, f"assignment={got_genotype} decision={expected_genotype}")
                ploidy = clean(decision.get("reconciled_biological_ploidy")).upper()
                expected_type = "D" if ploidy == "TETRAPLOID" else "S"
                got_type = clean(outrow.get("type")).upper()
                if got_type != expected_type:
                    fail("reconciled_assignments_ploidy_type", lib, bc, f"ploidy={decision.get('reconciled_biological_ploidy')} type={got_type} expected={expected_type}")
        finish("reconciled_assignments_invariants", b, lib)

        b = len(failures)
        if args.doublet_context_root:
            summary_path = os.path.join(args.doublet_context_root, f"{lib}.doublet_dragon_summary.tsv")
            input_path = os.path.join(args.doublet_context_root, f"{lib}.doublet_dragon_input.assignments")
            if not os.path.isfile(summary_path):
                fail("doublet_dragon_context_present", lib, "", summary_path)
            if os.path.isfile(input_path):
                dd_input = read_assignments(input_path)
                for bc, r in dd_input.items():
                    if clean(r.get("type")).upper() == "D" and canonical_genotype(r.get("assignment", "")) in global_lines:
                        fail("doublet_dragon_excludes_real_biological_lines", lib, bc, r.get("assignment", ""))
        finish("doublet_dragon_context_invariants", b, lib)

        b = len(failures)
        for row in lib_rows:
            bc = clean(row.get("barcode"))
            g = canonical_genotype(row.get("reconciled_donor_genotype", "")); uid = clean(row.get("reconciled_uid")); status = clean(row.get("uid_resolution_status"))
            meta = metadata[lib].get(g, {})
            expected_uid = clean(meta.get("reconciled_uid"))
            if status in {"EXACT_LIBRARY_METADATA_MATCH", "MULTIPLE_EXPECTED_UIDS_SAME_GENOTYPE"}:
                if uid != expected_uid:
                    fail("uid_scope_exact", lib, bc, f"result={uid} metadata={expected_uid}")
                if uid:
                    for token in uid.split("|"):
                        if token not in members[(lib, g)]: fail("uid_member_exists", lib, bc, f"uid={token} genotype={g}")
            elif status in {"EXACT_GLOBAL_METADATA_MATCH", "MULTIPLE_GLOBAL_UIDS_SAME_GENOTYPE", "GLOBAL_LINE_MISSING_UID"}:
                glob = global_line_rows.get(g, {})
                global_uids = clean(glob.get("uid_candidates"))
                if clean(row.get("uid_candidates")) != global_uids:
                    fail("uid_global_scope_exact", lib, bc, f"result={row.get('uid_candidates')} global={global_uids}")

            proposed = canonical_genotype(row.get("proposed_donor_genotype", ""))
            pstatus = clean(row.get("proposed_uid_resolution_status"))
            pcands = clean(row.get("proposed_uid_candidates"))
            if proposed in metadata[lib]:
                pm = metadata[lib][proposed]
                if pstatus != clean(pm.get("uid_resolution_status")):
                    fail("proposed_uid_library_scope", lib, bc, f"{proposed}: {pstatus}")
                if pcands != clean(pm.get("uid_candidates")):
                    fail("proposed_uid_library_candidates", lib, bc, f"{proposed}: {pcands}")
            elif proposed in global_line_rows:
                gm = global_line_rows[proposed]
                expected = clean(gm.get("uid_candidates"))
                if pstatus not in {"EXACT_GLOBAL_METADATA_MATCH", "MULTIPLE_GLOBAL_UIDS_SAME_GENOTYPE", "GLOBAL_LINE_MISSING_UID"}:
                    fail("proposed_uid_global_scope", lib, bc, f"{proposed}: {pstatus}")
                if pcands != expected:
                    fail("proposed_uid_global_candidates", lib, bc, f"{proposed}: {pcands} global={expected}")
            elif len(donor_components(proposed)) == 1 and proposed in global_donors:
                if pstatus != "GLOBAL_DONOR_NO_STANDALONE_LINE_UID":
                    fail("proposed_uid_donor_scope", lib, bc, f"{proposed}: {pstatus}")
        finish("uid_resolution_invariants", b, lib)

        b = len(failures)
        for row in lib_rows:
            bc = clean(row.get("barcode"))
            proposed = canonical_genotype(row.get("proposed_donor_genotype", ""))
            current = canonical_genotype(row.get("current_refined_assignment", "")) or canonical_genotype(row.get("original_demux_assignment", ""))
            proposed_comps = donor_components(proposed)
            current_comps = donor_components(current)
            mt_status = clean(row.get("mt_alternative_status"))
            mt_mode = clean(row.get("mt_verification_mode"))
            mt_prop_component = canonical_genotype(row.get("mt_proposed_component", ""))
            mt_cur_component = canonical_genotype(row.get("mt_current_component", ""))
            all_reasons = [x for x in clean(row.get("decision_reason_codes")).split(",") if x]
            mt_reasons = [x for x in all_reasons if x.startswith("MITOCHONDRIA_")]

            if clean(row.get("mt_haplotype_resolution")) == "MT_HAPLOTYPE_UNRESOLVED" and mt_status == "SUPPORTS_ALTERNATIVE":
                fail("mt_unresolved_not_supportive", lib, bc, "unresolved haplotype called supportive")

            mt_participates = mt_status in {"SUPPORTS_ALTERNATIVE", "SUPPORTS_CURRENT", "CONTRADICTS_ALTERNATIVE"} or bool(mt_reasons)
            singlet_allowed = len(proposed_comps) == 1 and mt_mode in {"", "SINGLET_IDENTITY"}
            disjoint_line_allowed = (
                len(proposed_comps) == 2
                and proposed in global_line_rows
                and bool(current_comps)
                and set(proposed_comps).isdisjoint(set(current_comps))
                and mt_mode == "DISJOINT_REAL_LINE_COMPONENT"
            )
            if mt_participates and not (singlet_allowed or disjoint_line_allowed):
                fail("mt_scope_guardrail", lib, bc, f"current={current} proposed={proposed} mode={mt_mode} status={mt_status}")

            if mt_mode == "DISJOINT_REAL_LINE_COMPONENT":
                if not disjoint_line_allowed:
                    fail("mt_disjoint_line_mode_requires_real_zero_overlap_line", lib, bc, f"current={current} proposed={proposed}")
                if mt_prop_component and mt_prop_component not in set(proposed_comps):
                    fail("mt_proposed_component_not_in_proposed_line", lib, bc, f"{mt_prop_component} not in {proposed}")
                if mt_cur_component and mt_cur_component not in set(current_comps):
                    fail("mt_current_component_not_in_current_line", lib, bc, f"{mt_cur_component} not in {current}")
                if mt_status == "SUPPORTS_ALTERNATIVE" and mt_prop_component not in set(proposed_comps):
                    fail("mt_line_support_requires_proposed_unique_component", lib, bc, f"component={mt_prop_component} proposed={proposed}")
                if mt_status == "SUPPORTS_CURRENT" and mt_cur_component not in set(current_comps):
                    fail("mt_line_current_support_requires_current_unique_component", lib, bc, f"component={mt_cur_component} current={current}")
                line_mt_reason = (
                    "MITOCHONDRIA_SUPPORT_DISJOINT_PROPOSED_LINE_COMPONENT" in mt_reasons
                    or "MITOCHONDRIA_SUPPORT_DISJOINT_CURRENT_LINE_COMPONENT_NONDECISIONAL" in mt_reasons
                )
                if line_mt_reason:
                    # The reconciler only emits these reasons after its full robust
                    # nuclear gate, which includes raw and depth-normalized effect
                    # size, informative depth, fold fraction and evaluable-fold count.
                    # Require the corresponding nuclear reason codes here rather than
                    # reimplementing only a subset of that gate from output columns.
                    if "NUCLEAR_ALTERNATIVE_STRONGLY_SUPPORTED" not in all_reasons or "NUCLEAR_SUPPORT_REPLICATED_ACROSS_SITE_FOLDS" not in all_reasons:
                        fail(
                            "mt_line_support_requires_preexisting_strong_nuclear_proposal",
                            lib, bc,
                            f"reasons={clean(row.get('decision_reason_codes'))}",
                        )

            for field in ("mt_best_identity", "mt_second_identity"):
                ident = canonical_genotype(row.get(field, ""))
                if ident and len(donor_components(ident)) != 1:
                    fail("mt_outputs_donor_only", lib, bc, f"{field}={ident}")

        mt_score_path = os.path.join(args.mt_root, f"{lib}.mt_identity_scores.tsv.gz")
        if os.path.isfile(mt_score_path):
            for mtrow in read_tsv(mt_score_path):
                g = canonical_genotype(mtrow.get("donor_genotype", ""))
                if g and len(donor_components(g)) != 1:
                    fail("mt_score_manifest_donor_only", lib, clean(mtrow.get("barcode")), g)
        finish("mt_guardrail", b, lib)

        b = len(failures)
        if args.evidence_mode == "rna":
            for row in lib_rows:
                if clean(row.get("atac_status")) != "ATAC_NOT_REQUESTED": fail("rna_mode_atac_not_requested", lib, clean(row.get("barcode")), clean(row.get("atac_status")))
            # A stale ATAC artifact from an earlier rna-atac run must not affect an
            # RNA-only result.  Isolation is proven from this run's manifest and cell
            # statuses, not by requiring the filesystem to be globally ATAC-free.
        finish("atac_mode_isolation", b, lib)

    b = len(failures)
    events_path = os.path.join(args.decisions_root, "all_libraries.identity_events.tsv")
    min_event_cells = int(DEFAULT_IDENTITY_POLICY.get("event_min_cells", 10))
    event_rows = read_tsv(events_path) if os.path.isfile(events_path) else []
    if os.path.isfile(events_path):
        for ev in event_rows:
            cls = clean(ev.get("event_class"))
            scope = clean(ev.get("event_scope"))
            ident = canonical_genotype(ev.get("unexpected_component", ""))
            n_imp = int(float(clean(ev.get("n_implicated_cells")) or 0))
            if cls.startswith("LIKELY_") and n_imp < min_event_cells:
                fail("event_mass_threshold", clean(ev.get("library")), "", f"{cls}: n={n_imp} < {min_event_cells}")
            if scope == "EXACT_IDENTITY" and cls.startswith("LIKELY_UNEXPECTED_"):
                comps = donor_components(ident)
                if len(comps) == 2 and ident not in global_lines:
                    fail("exact_identity_event_global_line", clean(ev.get("library")), "", ident)
                if len(comps) == 1 and ident not in global_donors:
                    fail("exact_identity_event_global_donor", clean(ev.get("library")), "", ident)
    finish("event_mass_threshold", b)

    b = len(failures)
    foreign_event_cells = defaultdict(Counter)
    foreign_event_ids = defaultdict(dict)
    foreign_event_libraries = defaultdict(set)
    for ev in event_rows:
        lib = clean(ev.get("library"))
        scope = clean(ev.get("event_scope"))
        ident = canonical_genotype(ev.get("unexpected_component", ""))
        expected_relationship, expected_context = expected_singlet_relationship(
            lib, ident, explicit_expected_singlets, expected_donor_components,
            expected_composites, global_donors,
        )
        observed_relationship = clean(ev.get("singlet_library_relationship"))
        observed_context = semicolon_set(ev.get("expected_composite_context", ""))
        contributes = clean(ev.get("contributes_to_library_exchange_evidence")).upper() == "TRUE"
        n_imp = numeric_int(ev.get("n_implicated_cells"))

        if scope == "EXACT_IDENTITY" and len(donor_components(ident)) == 1:
            if observed_relationship != expected_relationship:
                fail(
                    "event_singlet_library_relationship", lib, "",
                    f"{ident}: observed={observed_relationship} expected={expected_relationship}",
                )
            if expected_relationship == "EXPECTED_COMPOSITE_COMPONENT_SINGLET":
                if observed_context != expected_context:
                    fail(
                        "event_component_singlet_context", lib, "",
                        f"{ident}: observed={sorted(observed_context)} expected={sorted(expected_context)}",
                    )
                if n_imp >= min_event_cells and clean(ev.get("event_class")) != "EXPECTED_COMPOSITE_COMPONENT_SINGLET_POPULATION":
                    fail("component_singlet_event_class", lib, "", f"{ident}: {clean(ev.get('event_class'))}")
                if contributes:
                    fail("component_singlet_event_excluded_from_exchange", lib, "", ident)
            elif observed_context:
                fail("event_noncomponent_has_no_composite_context", lib, "", f"{ident}: {sorted(observed_context)}")

            expected_contributes = expected_relationship == "UNEXPECTED_SINGLET" and n_imp >= min_event_cells
            if contributes != expected_contributes:
                fail(
                    "event_exchange_evidence_eligibility", lib, "",
                    f"{ident}: observed={contributes} expected={expected_contributes}",
                )
            if expected_contributes:
                foreign_event_cells[lib][ident] = max(foreign_event_cells[lib][ident], n_imp)
                foreign_event_ids[lib][ident] = clean(ev.get("event_id"))
                foreign_event_libraries[ident].add(lib)
        elif contributes:
            fail("exchange_evidence_requires_exact_unexpected_singlet", lib, "", ident)
    for ev in event_rows:
        if clean(ev.get("event_scope")) != "EXACT_IDENTITY":
            continue
        ident = canonical_genotype(ev.get("unexpected_component", ""))
        if len(donor_components(ident)) != 1:
            continue
        lib = clean(ev.get("library"))
        relationship, _ = expected_singlet_relationship(
            lib, ident, explicit_expected_singlets, expected_donor_components,
            expected_composites, global_donors,
        )
        observed_recurrence = clean(ev.get("related_library_recurrence")) == "YES"
        if relationship == "UNEXPECTED_SINGLET" and numeric_int(ev.get("n_implicated_cells")) >= min_event_cells:
            expected_recurrence = len(foreign_event_libraries[ident]) > 1
            if observed_recurrence != expected_recurrence:
                fail("unexpected_singlet_recurrence_excludes_component_populations", lib, "", ident)
    finish("component_singlet_event_invariants", b)

    b = len(failures)
    for row in all_rows:
        lib = clean(row.get("library"))
        proposed = canonical_genotype(row.get("proposed_donor_genotype", ""))
        eligible = clean(row.get("library_exchange_evidence_eligible")).upper() == "TRUE"
        relationship, _ = expected_singlet_relationship(
            lib, proposed, explicit_expected_singlets, expected_donor_components,
            expected_composites, global_donors,
        )
        line_mass = numeric_int(row.get("alternative_line_event_mass"))
        expected_eligible = relationship == "UNEXPECTED_SINGLET" and line_mass >= min_event_cells
        if eligible and not expected_eligible:
            fail("cell_exchange_evidence_requires_unexpected_singlet", lib, clean(row.get("barcode")), proposed)
    finish("cell_exchange_evidence_invariants", b)

    b = len(failures)
    exchange_path = os.path.join(args.decisions_root, "all_libraries.library_exchange_events.tsv")
    if not os.path.isfile(exchange_path):
        fail("library_exchange_table_present", "", "", exchange_path)
        exchange_rows = []
    else:
        exchange_rows = read_tsv(exchange_path)

    component_barcodes = defaultdict(lambda: defaultdict(set))
    singlet_cell_counts = defaultdict(Counter)
    for row in all_rows:
        if clean(row.get("reconciled_droplet_state")) != "SINGLE_CELL":
            continue
        lib = clean(row.get("library")); bc = clean(row.get("barcode"))
        genotype = canonical_genotype(row.get("reconciled_donor_genotype", ""))
        comps = donor_components(genotype)
        for donor in set(comps):
            component_barcodes[lib][donor].add(bc)
        if len(comps) == 1:
            singlet_cell_counts[lib][comps[0]] += 1

    expected_pairs = []
    lib_names = [f"lib{n}" for n in libs]
    for i, library_a in enumerate(lib_names):
        for library_b in lib_names[i + 1:]:
            expected_pairs.append((library_a, library_b))
    exchange_by_pair = defaultdict(list)
    for row in exchange_rows:
        exchange_by_pair[(clean(row.get("library_a")), clean(row.get("library_b")))].append(row)
    if set(exchange_by_pair) != set(expected_pairs):
        missing = sorted(set(expected_pairs) - set(exchange_by_pair))
        extra = sorted(set(exchange_by_pair) - set(expected_pairs))
        fail("library_exchange_pair_roster", "", "", f"missing={missing[:10]} extra={extra[:10]}")

    for pair in expected_pairs:
        rows = exchange_by_pair.get(pair, [])
        if len(rows) != 1:
            fail("library_exchange_one_row_per_pair", pair[0], "", f"{pair[1]} rows={len(rows)}")
            continue
        row = rows[0]
        library_a, library_b = pair
        relation, discriminability, shared, a_specific, b_specific = expected_pair_roster(
            library_a, library_b, expected_donor_components,
        )
        if enum_value(row.get("roster_relation")) != relation:
            fail("library_exchange_roster_relation", library_a, "", f"{library_b}: {row.get('roster_relation')} != {relation}")
        if enum_value(row.get("pair_discriminability")) != discriminability:
            fail("library_exchange_pair_discriminability", library_a, "", f"{library_b}: {row.get('pair_discriminability')} != {discriminability}")
        if semicolon_set(row.get("shared_donors", "")) != shared:
            fail("library_exchange_shared_donors", library_a, "", library_b)
        if semicolon_set(row.get("a_specific_donors", "")) != a_specific:
            fail("library_exchange_a_specific_donors", library_a, "", library_b)
        if semicolon_set(row.get("b_specific_donors", "")) != b_specific:
            fail("library_exchange_b_specific_donors", library_a, "", library_b)
        if numeric_int(row.get("n_shared_donors")) != len(shared):
            fail("library_exchange_n_shared", library_a, "", library_b)
        if numeric_int(row.get("n_a_specific_donors")) != len(a_specific) or numeric_int(row.get("n_b_specific_donors")) != len(b_specific):
            fail("library_exchange_n_specific", library_a, "", library_b)

        a_native_detected = {d for d in a_specific if len(component_barcodes[library_a].get(d, set())) >= min_event_cells}
        b_native_detected = {d for d in b_specific if len(component_barcodes[library_b].get(d, set())) >= min_event_cells}
        a_foreign_in_b = {d for d in a_specific if foreign_event_cells[library_b].get(d, 0) >= min_event_cells}
        b_foreign_in_a = {d for d in b_specific if foreign_event_cells[library_a].get(d, 0) >= min_event_cells}

        expected_counts = {
            "a_specific_detected_in_a": len(a_native_detected),
            "a_specific_detected_in_b": len(a_foreign_in_b),
            "b_specific_detected_in_b": len(b_native_detected),
            "b_specific_detected_in_a": len(b_foreign_in_a),
        }
        for field, expected_value in expected_counts.items():
            if numeric_int(row.get(field)) != expected_value:
                fail("library_exchange_detection_count", library_a, "", f"{library_b} {field}={row.get(field)} expected={expected_value}")

        a_native_cells = set()
        for donor in a_specific:
            a_native_cells.update(component_barcodes[library_a].get(donor, set()))
        b_native_cells = set()
        for donor in b_specific:
            b_native_cells.update(component_barcodes[library_b].get(donor, set()))
        if numeric_int(row.get("a_specific_cells_in_a")) != len(a_native_cells):
            fail("library_exchange_a_native_cells", library_a, "", library_b)
        if numeric_int(row.get("b_specific_cells_in_b")) != len(b_native_cells):
            fail("library_exchange_b_native_cells", library_a, "", library_b)
        if numeric_int(row.get("a_specific_cells_in_b")) != sum(foreign_event_cells[library_b].get(d, 0) for d in a_foreign_in_b):
            fail("library_exchange_a_foreign_cells", library_a, "", library_b)
        if numeric_int(row.get("b_specific_cells_in_a")) != sum(foreign_event_cells[library_a].get(d, 0) for d in b_foreign_in_a):
            fail("library_exchange_b_foreign_cells", library_a, "", library_b)

        def expected_fraction(numerator, denominator):
            return numerator / denominator if denominator else math.nan

        expected_fracs = {
            "a_signature_coverage_in_b": expected_fraction(len(a_foreign_in_b), len(a_specific)),
            "b_signature_coverage_in_a": expected_fraction(len(b_foreign_in_a), len(b_specific)),
            "a_native_retention_fraction": expected_fraction(len(a_native_detected), len(a_specific)),
            "b_native_retention_fraction": expected_fraction(len(b_native_detected), len(b_specific)),
        }
        expected_fracs["a_native_displacement_fraction"] = 1.0 - expected_fracs["a_native_retention_fraction"] if math.isfinite(expected_fracs["a_native_retention_fraction"]) else math.nan
        expected_fracs["b_native_displacement_fraction"] = 1.0 - expected_fracs["b_native_retention_fraction"] if math.isfinite(expected_fracs["b_native_retention_fraction"]) else math.nan
        for field, expected_value in expected_fracs.items():
            observed = numeric_float(row.get(field))
            if math.isfinite(expected_value):
                if not math.isfinite(observed) or abs(observed - expected_value) > 1e-8:
                    fail("library_exchange_fraction_consistency", library_a, "", f"{library_b} {field}={row.get(field)} expected={expected_value}")
            elif math.isfinite(observed):
                fail("library_exchange_fraction_consistency", library_a, "", f"{library_b} {field} should be NA")

        a_strength = expected_directional_strength(len(a_foreign_in_b), len(a_specific))
        b_strength = expected_directional_strength(len(b_foreign_in_a), len(b_specific))
        if enum_value(row.get("a_to_b_exchange_strength")) != a_strength or enum_value(row.get("b_to_a_exchange_strength")) != b_strength:
            fail("library_exchange_directional_strength", library_a, "", library_b)
        reciprocal, interpretation, confidence = expected_exchange_classification(
            relation, discriminability, a_strength, b_strength,
            expected_fracs["a_native_retention_fraction"], expected_fracs["b_native_retention_fraction"],
        )
        if enum_value(row.get("reciprocal_exchange_status")) != reciprocal:
            fail("library_exchange_reciprocal_status", library_a, "", library_b)
        if enum_value(row.get("exchange_interpretation")) != interpretation:
            fail("library_exchange_interpretation", library_a, "", f"{library_b}: {row.get('exchange_interpretation')} != {interpretation}")
        if enum_value(row.get("exchange_confidence")) != confidence:
            fail("library_exchange_confidence", library_a, "", library_b)
        if relation == "ROSTER_EQUIVALENT_NONDISCRIMINATING" and enum_value(row.get("exchange_interpretation")) != "ROSTER_EQUIVALENT_NONDISCRIMINATING":
            fail("identical_roster_nonidentifiability", library_a, "", library_b)

        supporting_expected = {
            foreign_event_ids[library_b].get(d, "") for d in a_foreign_in_b
        } | {
            foreign_event_ids[library_a].get(d, "") for d in b_foreign_in_a
        }
        supporting_expected.discard("")
        if semicolon_set(row.get("supporting_event_ids", "")) != supporting_expected:
            fail("library_exchange_supporting_events", library_a, "", library_b)

        for donor in shared:
            if donor in a_foreign_in_b or donor in b_foreign_in_a:
                fail("shared_donor_zero_exchange_evidence", library_a, "", f"{library_b}: {donor}")
    finish("library_exchange_invariants", b)

    b = len(failures)
    donor_evidence_path = os.path.join(args.decisions_root, "all_libraries.library_exchange_donor_evidence.tsv")
    if not os.path.isfile(donor_evidence_path):
        fail("library_exchange_donor_evidence_present", "", "", donor_evidence_path)
    else:
        for row in read_tsv(donor_evidence_path):
            source = clean(row.get("source_library")); target = clean(row.get("target_library"))
            donor = clean(row.get("diagnostic_donor"))
            if donor not in expected_donor_components.get(source, set()) or donor in expected_donor_components.get(target, set()):
                fail("exchange_donor_must_be_source_specific", source, "", f"target={target} donor={donor}")
            if clean(row.get("shared_or_specific")) != "SPECIFIC":
                fail("exchange_donor_shared_excluded", source, "", f"target={target} donor={donor}")
            target_cells = numeric_int(row.get("target_observed_robust_singlet_cells"))
            if target_cells > 0:
                relationship, _ = expected_singlet_relationship(
                    target, donor, explicit_expected_singlets, expected_donor_components,
                    expected_composites, global_donors,
                )
                if relationship != "UNEXPECTED_SINGLET":
                    fail("exchange_donor_target_must_be_unexpected_singlet", target, "", donor)
                if target_cells != foreign_event_cells[target].get(donor, 0):
                    fail("exchange_donor_target_cell_count", target, "", donor)
    finish("library_exchange_donor_evidence_invariants", b)

    b = len(failures)
    if len(all_rows) != sum(len(read_tsv(os.path.join(args.decisions_root, f"lib{n}.reconciled_cells.tsv.gz"))) for n in libs):
        fail("aggregate_row_count", "", "", "all-libraries count does not equal per-library sum")
    finish("aggregate_row_count", b)

    b = len(failures)
    manifest_path = os.path.join(args.decisions_root, "reconciliation_manifest.json")
    try:
        manifest = json.loads(Path(manifest_path).read_text())
        if manifest.get("schema_version") != SCHEMA_VERSION: fail("manifest_schema", "", "", str(manifest.get("schema_version")))
        if manifest.get("policy_version") != POLICY_VERSION: fail("manifest_policy", "", "", str(manifest.get("policy_version")))
        for key in ("contam_dependency", "empty_drops_dependency", "tetra_refine_rerun_dependency"):
            if manifest.get(key): fail("forbidden_dependency", "", "", key)
        if manifest.get("ambient_rna_evaluated") is not False: fail("ambient_not_evaluated", "", "", str(manifest.get("ambient_rna_evaluated")))
        if manifest.get("global_biological_line_source") != "2025_LineMeta":
            fail("global_biological_line_source", "", "", str(manifest.get("global_biological_line_source")))
        if args.doublet_context_root and manifest.get("doublet_dragon_role") != "population_context_only_diploid_resolvable_subset":
            fail("doublet_dragon_role", "", "", str(manifest.get("doublet_dragon_role")))
        if int(manifest.get("event_min_cells", -1)) != int(DEFAULT_IDENTITY_POLICY.get("event_min_cells", 10)):
            fail("event_min_cells_manifest", "", "", str(manifest.get("event_min_cells")))
        if int(manifest.get("unexpected_line_autoapply_min_cells", -1)) != int(DEFAULT_IDENTITY_POLICY.get("unexpected_line_autoapply_min_cells", 100)):
            fail("unexpected_line_autoapply_min_cells_manifest", "", "", str(manifest.get("unexpected_line_autoapply_min_cells")))
        if manifest.get("occupancy_aware_reconciliation") is not True:
            fail("occupancy_aware_reconciliation_manifest", "", "", str(manifest.get("occupancy_aware_reconciliation")))
        if manifest.get("diploid_population_definition") != "expected_singlets_plus_robust_singlet_events_plus_component_derived_singlet_populations":
            fail("diploid_population_definition_manifest", "", "", str(manifest.get("diploid_population_definition")))
        if manifest.get("heterotypic_pair_occupancy_guard") is not True:
            fail("heterotypic_pair_occupancy_guard_manifest", "", "", str(manifest.get("heterotypic_pair_occupancy_guard")))
        if manifest.get("unexpected_homotet_occupancy_guard") is not True:
            fail("unexpected_homotet_occupancy_guard_manifest", "", "", str(manifest.get("unexpected_homotet_occupancy_guard")))
        if manifest.get("component_singlet_classification") is not True:
            fail("component_singlet_classification_manifest", "", "", str(manifest.get("component_singlet_classification")))
        if manifest.get("library_exchange_inference") is not True:
            fail("library_exchange_inference_manifest", "", "", str(manifest.get("library_exchange_inference")))
        if manifest.get("library_exchange_shared_donors_excluded") is not True:
            fail("library_exchange_shared_donors_manifest", "", "", str(manifest.get("library_exchange_shared_donors_excluded")))
        if manifest.get("auto_apply_enabled") is False:
            for row in all_rows:
                if clean(row.get("reassignment_applied")).upper() == "TRUE":
                    fail("auto_apply_disabled_no_rewrites", clean(row.get("library")), clean(row.get("barcode")), clean(row.get("final_action")))
    except Exception as exc:
        fail("manifest_readable", "", "", str(exc))
    finish("manifest_invariants", b)

    out = Path(args.output_root); out.mkdir(parents=True, exist_ok=True)
    write_tsv(str(out / "validation_summary.tsv"), summaries, SUMMARY_FIELDS)
    write_tsv(str(out / "validation_failures.tsv"), failures, FAIL_FIELDS)
    report = ["# Identity reconciliation validation", "", f"Evidence mode: {args.evidence_mode}", f"Failures: {len(failures)}", ""]
    for s in summaries: report.append(f"- {s['status']} {s['check']}: {s['n_failures']} failure(s) {s['detail']}")
    if failures:
        report += ["", "## Failures"] + [f"- {f['check']} {f['library']} {f['barcode']}: {f['detail']}" for f in failures[:200]]
    (out / "validation_report.md").write_text("\n".join(report) + "\n")
    # Structural corruption/invariant violations fail; unresolved science does not.
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
