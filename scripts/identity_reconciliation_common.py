#!/usr/bin/env python3
"""Shared helpers for the CellBouncer identity-reconciliation pipeline.

This module deliberately contains no ambient-RNA/contamination integration.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import re
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Optional, Sequence, Set, Tuple

SCHEMA_VERSION = "identity_reconciliation_v1"
POLICY_VERSION = "identity_reconciliation_policy_v11_component_singlets_and_library_exchange"
NA_TOKENS = {"", ".", "na", "nan", "none", "null", "unavailable"}

DOWNSTREAM_SAFE_ASSIGNMENT_FIELDS = (
    "downstream_safe_assignment",
    "downstream_safe_assignment_source",
    "downstream_safe_change_applied",
    "downstream_assignment_status",
)
DOWNSTREAM_SAFE_ASSIGNMENT_SOURCE_CURRENT = "CURRENT_ASSIGNMENT"
DOWNSTREAM_SAFE_ASSIGNMENT_SOURCE_APPROVED = "APPROVED_RECONCILIATION"
DOWNSTREAM_SAFE_ASSIGNMENT_SOURCE_UNAVAILABLE = "UNAVAILABLE"
DOWNSTREAM_ASSIGNMENT_READY = "RECONCILED_CHANGE_READY"
DOWNSTREAM_ASSIGNMENT_HELD = "CURRENT_RETAINED_HELD_CHANGE"
DOWNSTREAM_ASSIGNMENT_NO_CHANGE = "CURRENT_RETAINED_NO_CHANGE"
DOWNSTREAM_ASSIGNMENT_INVALID_PROPOSAL = "CURRENT_RETAINED_INVALID_PROPOSAL"
DOWNSTREAM_ASSIGNMENT_NO_CURRENT = "NO_VALID_CURRENT_ASSIGNMENT"

ASSIGNMENT_STATUS_FINE = "FINE_NO_CHANGE"
ASSIGNMENT_STATUS_APPLIED = "CHANGE_APPLIED"
ASSIGNMENT_STATUS_REVIEW = "REVIEW_NEEDED"
ASSIGNMENT_STATUSES = (
    ASSIGNMENT_STATUS_FINE,
    ASSIGNMENT_STATUS_APPLIED,
    ASSIGNMENT_STATUS_REVIEW,
)
THREE_STATE_ASSIGNMENT_FIELDS = (
    "assignment_status",
    "current_assignment",
    "proposed_assignment",
    "final_assignment",
    "assignment_change",
    "review_reason",
    "evidence_summary",
)

# One explicit interpretation point for every action emitted by reconciliation.
# These classes describe the upstream diagnostic action; they do not by
# themselves define the primary three-state identity question.  In particular,
# REVIEW actions may describe ploidy, occupancy, event, or other context while
# the production identity remains unchanged.  A primary identity question is
# established separately from the preliminary production assignment below.
RECONCILIATION_ACTION_CLASS = {
    # Blank/NA is the normal state for cells that never entered a cell-level
    # reconciliation decision because no identity change was proposed.
    "": "NO_CHANGE",
    "KEEP": "NO_CHANGE",
    "REASSIGN_GENOTYPE": "CHANGE",
    "RECLASSIFY_PLOIDY": "CHANGE",
    "REVIEW_CELLULAR_ORIGIN": "REVIEW",
    "REVIEW_UNEXPECTED_IDENTITY": "REVIEW",
    "REVIEW_HOMOTET_OCCUPANCY": "REVIEW",
    "KEEP_CURRENT_CONFLICTED": "REVIEW",
    "UNRESOLVED_INSUFFICIENT_EVIDENCE": "REVIEW",
}

DEFAULT_IDENTITY_POLICY = {
    "policy_version": POLICY_VERSION,
    "schema_version": SCHEMA_VERSION,
    # Conservative pilot thresholds. These are intentionally not treated as a
    # calibrated production policy; auto-apply is disabled by default below.
    "nuclear_strong_delta_ll": 100.0,
    "nuclear_suggestive_delta_ll": 20.0,
    "nuclear_strong_depth_normalized_delta": 0.01,
    "nuclear_min_informative_depth": 5000,
    "nuclear_low_information_depth": 500,
    "site_fold_min_evaluable": 4,
    "site_fold_support_fraction": 0.80,
    "mt_support_delta_ll": 4.0,
    "atac_support_delta_ll": 100.0,
    "nn_homotypic_prob_tet": 0.90,
    "event_min_cells": 10,
    # Unexpected non-homotet identities may be reported as coherent events at
    # event_min_cells, but production assignments are only rewritten once the
    # exact canonical identity has this many robust cells in the same library.
    # NN-driven A -> A+A homotet reclassification is exempt from this mass
    # threshold, but the occupancy guards below still apply.
    "unexpected_line_autoapply_min_cells": 100,
    # A repeated A+B donor-composition signal does not by itself distinguish an
    # intact biological T[A+B] cell from a technical M{D[A]|D[B]} droplet when
    # both donor singlets are independently present in the same library.
    "heterotypic_pair_occupancy_guard": True,
    # Likewise, a library-unexpected A+A call is not auto-applied from NN ploidy
    # evidence alone when diploid A is independently present; A+A and a same-
    # donor technical doublet are then observationally confounded here.
    "unexpected_homotet_occupancy_guard": True,
    # A positive quad-pattern score is genotype-composition evidence, not
    # independent proof that multiple cells occupied the droplet.  Only an
    # explicit non-single droplet flag is allowed to trigger the hard
    # technical-multiplet veto.
    "explicit_multiplet_evidence_source": "droplet_flag_only",
    "event_partner_preservation_fraction": 0.80,
    "event_primary_source_fraction": 0.80,
    "full_replacement_fraction": 0.90,
    "partial_replacement_fraction": 0.20,
    "foreign_carryover_max_fraction": 0.10,
    # Explicit opt-in until held-out calibration establishes the false-change rate.
    "auto_apply_decisive": False,
}
DEFAULT_IDENTITY_POLICY_SHA256 = hashlib.sha256(
    json.dumps(DEFAULT_IDENTITY_POLICY, sort_keys=True, separators=(",", ":")).encode("utf-8")
).hexdigest()


def clean(value) -> str:
    if value is None:
        return ""
    s = str(value).strip()
    return "" if s.lower() in NA_TOKENS else s


def natural_key(value: str):
    s = clean(value)
    if re.fullmatch(r"\d+", s):
        return (0, int(s), "")
    return (1, 0, s)


def canonical_uid_set(values: Iterable[str]) -> str:
    toks = set()
    for value in values:
        for token in re.split(r"[|,;]", clean(value)):
            token = clean(token)
            if token:
                toks.add(token)
    return "|".join(sorted(toks, key=natural_key))


def uid_count(value: str) -> int:
    return len([x for x in clean(value).split("|") if x])


def canonical_genotype(value: str, aliases: Optional[Mapping[str, str]] = None) -> str:
    aliases = aliases or {}
    raw = clean(value)
    if not raw:
        return ""
    # Remove state wrappers used by the reconciliation layer.
    m = re.fullmatch(r"(?:D|T|UNKNOWN_SINGLE_CELL)\[(.*)\]", raw)
    if m:
        raw = m.group(1)
    if raw.startswith("M{"):
        return raw
    parts = [clean(x) for x in raw.replace("x", "+").split("+")]
    parts = [aliases.get(x, x) for x in parts if x]
    return "+".join(sorted(parts))


def donor_components(genotype: str) -> List[str]:
    g = canonical_genotype(genotype)
    if not g or g.startswith("M{"):
        return []
    return [x for x in g.split("+") if x]


def reconciliation_action_class(value: object) -> str:
    """Return NO_CHANGE, CHANGE, REVIEW, or UNKNOWN for a cell action."""
    return RECONCILIATION_ACTION_CLASS.get(clean(value).upper(), "UNKNOWN")


def _first_canonical_genotype(*values: object) -> str:
    for value in values:
        genotype = canonical_genotype(value)
        if genotype:
            return genotype
    return ""


def has_real_cell_change_proposal(row: Mapping[str, object]) -> bool:
    """Return whether a distinct primary production change was proposed.

    Candidate nominations and diagnostic action names are not sufficient.  The
    cell must have a change-producing action, a legacy selected production
    change when the action is absent, or an explicit cell-scoped review
    decision.  Event-scoped review metadata never creates a primary cell-level
    identity question.
    """
    current = _first_canonical_genotype(
        row.get("comparison_current_assignment", ""),
        row.get("current_assignment", ""),
        row.get("refined_assignment", ""),
        row.get("demux_original_assignment", ""))
    proposal = _first_canonical_genotype(
        row.get("nominated_proposal", ""),
        row.get("proposed_assignment", ""))
    preliminary = _first_canonical_genotype(
        row.get("preliminary_reconciled_assignment", ""))
    preliminary_applied = clean(
        row.get("preliminary_action_applied", "")).lower() in {
            "1", "true", "yes", "y"}
    action = (
        clean(row.get("preliminary_reconciliation_action", ""))
        or clean(row.get("reconciliation_final_action", ""))).upper()
    action_class = reconciliation_action_class(action)
    disposition = clean(row.get("review_disposition", "")).upper()
    cell_review = (
        clean(row.get("review_record_scope", "")).upper() == "CELL"
        and disposition in {
            "ACCEPT_PROPOSAL", "KEEP_CURRENT", "LEAVE_UNRESOLVED",
            "PENDING"})
    preliminary_selected_change = bool(
        (preliminary and preliminary != current) or preliminary_applied)
    return bool(
        current and proposal and proposal != current
        and (cell_review or action_class == "CHANGE"
             or (not action and preliminary_selected_change)))


def _finite_float(value: object) -> float:
    try:
        parsed = float(clean(value))
        return parsed if math.isfinite(parsed) else math.nan
    except (TypeError, ValueError):
        return math.nan


def _reason_tokens(value: object) -> Set[str]:
    return {
        clean(token).upper()
        for token in re.split(r"[,;]", clean(value))
        if clean(token)
    }


def _three_state_review_reason(row: Mapping[str, object], evidence_mode: str,
                               donor_identity_change: bool,
                               structure_change: bool) -> Tuple[str, List[str]]:
    """Interpret only evidence relevant to the exact current/proposed contrast."""
    reasons: List[str] = []
    nuclear = clean(row.get("nuclear_reconciliation_status")).upper()
    atac = clean(row.get("atac_evidence_status")).upper()
    mito = clean(row.get("mitochondrial_evidence_status")).upper()
    ploidy = clean(row.get("ploidy_evidence_status")).upper()
    occupancy = clean(row.get("occupancy_evidence_status")).upper()
    technical = clean(row.get("technical_state")).upper()

    # A one-donor/two-donor structure change may also change donor identity.
    # Nuclear and ATAC evidence are therefore checked whenever the donor set
    # changes; ploidy/occupancy checks below are additional, not substitutes.
    # Pure A <-> A+A ploidy reclassification is evaluated by the structure
    # evidence because RNA/ATAC donor evidence cannot distinguish copy count.
    if donor_identity_change:
        modalities_disagree = (
            evidence_mode == "rna-atac"
            and ((nuclear == "NUCLEAR_SUPPORTS_PROPOSAL"
                  and atac in {"ATAC_SUPPORTS_CURRENT", "SUPPORTS_CURRENT"})
                 or (nuclear == "NUCLEAR_SUPPORTS_CURRENT"
                     and atac == "ATAC_SUPPORTS_ALTERNATIVE")))
        if modalities_disagree:
            reasons.append("RNA and ATAC disagree")
        elif nuclear == "NUCLEAR_UNAVAILABLE":
            reasons.append("RNA comparison not run")
        elif nuclear != "NUCLEAR_SUPPORTS_PROPOSAL":
            reasons.append("RNA evidence is insufficient")

        if evidence_mode == "rna-atac" and not modalities_disagree:
            if atac != "ATAC_SUPPORTS_ALTERNATIVE":
                reasons.append("ATAC evidence is insufficient")

    decision_tokens = (
        _reason_tokens(row.get("preliminary_decision_reason_codes", ""))
        | _reason_tokens(row.get("decision_reason_codes", "")))
    mt_current_is_explicitly_nondecisional = (
        "MITOCHONDRIA_SUPPORT_DISJOINT_CURRENT_LINE_COMPONENT_NONDECISIONAL"
        in decision_tokens)
    if (donor_identity_change and mito in {
            "SUPPORTS_CURRENT", "CONTRADICTS", "CONTRADICTS_ALTERNATIVE",
            "MITO_SUPPORTS_CURRENT", "MITO_CONTRADICTS"}
            and not mt_current_is_explicitly_nondecisional):
        reasons.append("Mitochondrial evidence conflicts")

    if structure_change and ploidy != "SUPPORTS_PROPOSAL":
        reasons.append("Ploidy or cellular state is unresolved")
    if structure_change and ("UNRESOLVED" in occupancy or "AMBIG" in occupancy):
        reasons.append("Ploidy or cellular state is unresolved")
    if structure_change and (technical not in {"", "NOT_APPLICABLE", "NA"}):
        reasons.append("Technical mixture is unresolved")

    frozen = _finite_float(
        row.get("ambient_frozen_proposal_minus_current_c"))
    refitted = _finite_float(row.get("ambient_assignment_effect_c_minus_b"))
    if (math.isfinite(frozen) and math.isfinite(refitted)
            and frozen != 0 and refitted != 0
            and (frozen < 0) != (refitted < 0)):
        reasons.append("Ambient correction changes the result")

    reasons = list(dict.fromkeys(reasons))
    if len(reasons) > 1:
        return "Multiple evidence sources conflict", reasons
    if reasons:
        return reasons[0], reasons
    return "", []


def derive_three_state_assignment(
        row: Mapping[str, object], evidence_mode: str = "rna") -> Dict[str, str]:
    """Derive the sole user-facing assignment status and safe assignment.

    Optional diagnostics never create a review case until the cell action first
    establishes a real current-versus-proposed identity question.
    """
    current = _first_canonical_genotype(
        row.get("comparison_current_assignment", ""),
        row.get("current_assignment", ""),
        row.get("refined_assignment", ""),
        row.get("demux_original_assignment", ""))
    proposal = _first_canonical_genotype(
        row.get("nominated_proposal", ""),
        row.get("proposed_assignment", ""))
    action = (
        clean(row.get("preliminary_reconciliation_action", ""))
        or clean(row.get("reconciliation_final_action", ""))).upper()
    action_class = reconciliation_action_class(action)
    disposition = clean(row.get("review_disposition", "")).upper()
    cell_disposition = (
        disposition
        if clean(row.get("review_record_scope", "")).upper() == "CELL"
        else "")

    def result(status: str, final: str, reason: str = "NONE",
               details: Optional[Sequence[str]] = None) -> Dict[str, str]:
        # Context/event nominations remain available in nominated_proposal in
        # the verbose ledger.  The compact primary proposal is populated only
        # when there is an applied or unresolved production identity change.
        proposed = (
            proposal if status != ASSIGNMENT_STATUS_FINE
            and proposal and proposal != current else "")
        change = (
            f"{current} -> {proposed}"
            if status != ASSIGNMENT_STATUS_FINE and proposed else "NONE")
        if status == ASSIGNMENT_STATUS_APPLIED:
            summary = (
                "Direct evidence supports the proposed assignment with no "
                "relevant unresolved conflict.")
        elif status == ASSIGNMENT_STATUS_REVIEW:
            detail_text = "; ".join(details or ([reason] if reason != "NONE" else []))
            summary = detail_text or "The proposed identity requires review."
        else:
            summary = "Current assignment retained; no cell-level identity change is recommended."
        release = (
            "HELD_FOR_REVIEW" if status == ASSIGNMENT_STATUS_REVIEW else "READY")
        legacy_status = (
            DOWNSTREAM_ASSIGNMENT_READY
            if status == ASSIGNMENT_STATUS_APPLIED else
            DOWNSTREAM_ASSIGNMENT_HELD
            if status == ASSIGNMENT_STATUS_REVIEW else
            DOWNSTREAM_ASSIGNMENT_NO_CHANGE)
        source = (
            DOWNSTREAM_SAFE_ASSIGNMENT_SOURCE_APPROVED
            if status == ASSIGNMENT_STATUS_APPLIED else
            DOWNSTREAM_SAFE_ASSIGNMENT_SOURCE_CURRENT)
        return {
            "assignment_status": status,
            "current_assignment": current or "NA",
            "proposed_assignment": proposed or "NA",
            "final_assignment": final or current or "NA",
            "assignment_change": change,
            "review_reason": reason if status == ASSIGNMENT_STATUS_REVIEW else "NONE",
            "evidence_summary": summary,
            # Deprecated compatibility aliases.  production_assignment is the
            # safe final assignment, never a held proposal.
            "production_assignment": final or current or "NA",
            "downstream_release_status": release,
            "downstream_safe_assignment": final or current or "NA",
            "downstream_safe_assignment_source": source,
            "downstream_safe_change_applied": (
                "TRUE" if status == ASSIGNMENT_STATUS_APPLIED else "FALSE"),
            "downstream_assignment_status": legacy_status,
            "review_required": (
                "TRUE" if status == ASSIGNMENT_STATUS_REVIEW else "FALSE"),
        }

    if action_class == "UNKNOWN":
        return result(
            ASSIGNMENT_STATUS_REVIEW, current,
            "Unknown reconciliation action", ["Unknown reconciliation action"])

    # Establish the cell-level question before considering review records or
    # evidence.  Event-wide review metadata cannot turn KEEP/context rows into
    # changes, and a proposal equal to current is never a change question.
    real_change_proposal = has_real_cell_change_proposal(row)
    if not real_change_proposal:
        return result(ASSIGNMENT_STATUS_FINE, current)

    if cell_disposition == "KEEP_CURRENT":
        return result(ASSIGNMENT_STATUS_FINE, current)
    if cell_disposition == "LEAVE_UNRESOLVED":
        return result(
            ASSIGNMENT_STATUS_REVIEW, current,
            "Multiple evidence sources conflict",
            ["Explicit review left the identity unresolved"])
    if cell_disposition == "PENDING":
        return result(
            ASSIGNMENT_STATUS_REVIEW, current,
            "Multiple evidence sources conflict",
            ["Explicit cell review remains pending"])
    if cell_disposition == "ACCEPT_PROPOSAL":
        return result(ASSIGNMENT_STATUS_APPLIED, proposal)

    if action_class == "REVIEW":
        preferred = {
            "REVIEW_CELLULAR_ORIGIN": "Technical mixture is unresolved",
            "REVIEW_HOMOTET_OCCUPANCY":
                "Ploidy or cellular state is unresolved",
            "REVIEW_UNEXPECTED_IDENTITY": "RNA evidence is insufficient",
            "KEEP_CURRENT_CONFLICTED": "Multiple evidence sources conflict",
            "UNRESOLVED_INSUFFICIENT_EVIDENCE":
                "RNA evidence is insufficient",
        }.get(action, "Multiple evidence sources conflict")
        return result(ASSIGNMENT_STATUS_REVIEW, current, preferred, [preferred])

    structure_change = (
        len(donor_components(current)) != len(donor_components(proposal)))
    donor_identity_change = (
        set(donor_components(current)) != set(donor_components(proposal)))
    reason, details = _three_state_review_reason(
        row, evidence_mode, donor_identity_change, structure_change)
    if reason:
        return result(ASSIGNMENT_STATUS_REVIEW, current, reason, details)
    return result(ASSIGNMENT_STATUS_APPLIED, proposal)


def derive_downstream_safe_assignment(row: Mapping[str, object]) -> Dict[str, str]:
    """Derive the assignment that may be consumed without review leakage.

    Three-state rows use ``final_assignment`` as authoritative and expose
    ``production_assignment`` only as a deprecated safe alias.  The older
    production/release derivation remains below for resumable v3 checkpoints.
    """
    primary_status = clean(row.get("assignment_status", "")).upper()
    if primary_status in ASSIGNMENT_STATUSES:
        final = canonical_genotype(row.get("final_assignment", ""))
        current = canonical_genotype(
            row.get("current_assignment", "")
            or row.get("comparison_current_assignment", ""))
        if primary_status == ASSIGNMENT_STATUS_APPLIED:
            status = DOWNSTREAM_ASSIGNMENT_READY
            source = DOWNSTREAM_SAFE_ASSIGNMENT_SOURCE_APPROVED
        elif primary_status == ASSIGNMENT_STATUS_REVIEW:
            status = DOWNSTREAM_ASSIGNMENT_HELD
            source = DOWNSTREAM_SAFE_ASSIGNMENT_SOURCE_CURRENT
        else:
            status = DOWNSTREAM_ASSIGNMENT_NO_CHANGE
            source = DOWNSTREAM_SAFE_ASSIGNMENT_SOURCE_CURRENT
        return {
            "downstream_safe_assignment": final or current or "NA",
            "downstream_safe_assignment_source": source,
            "downstream_safe_change_applied": (
                "TRUE" if primary_status == ASSIGNMENT_STATUS_APPLIED else "FALSE"),
            "downstream_assignment_status": status,
        }

    current = canonical_genotype(row.get("comparison_current_assignment", ""))
    production = canonical_genotype(row.get("production_assignment", ""))
    release = clean(row.get("downstream_release_status", "")).upper()

    if not current:
        return {
            "downstream_safe_assignment": "NA",
            "downstream_safe_assignment_source":
                DOWNSTREAM_SAFE_ASSIGNMENT_SOURCE_UNAVAILABLE,
            "downstream_safe_change_applied": "FALSE",
            "downstream_assignment_status": DOWNSTREAM_ASSIGNMENT_NO_CURRENT,
        }
    if not production:
        return {
            "downstream_safe_assignment": current,
            "downstream_safe_assignment_source":
                DOWNSTREAM_SAFE_ASSIGNMENT_SOURCE_CURRENT,
            "downstream_safe_change_applied": "FALSE",
            "downstream_assignment_status":
                DOWNSTREAM_ASSIGNMENT_INVALID_PROPOSAL,
        }
    if current != production and release == "READY":
        return {
            "downstream_safe_assignment": production,
            "downstream_safe_assignment_source":
                DOWNSTREAM_SAFE_ASSIGNMENT_SOURCE_APPROVED,
            "downstream_safe_change_applied": "TRUE",
            "downstream_assignment_status": DOWNSTREAM_ASSIGNMENT_READY,
        }
    if current != production:
        return {
            "downstream_safe_assignment": current,
            "downstream_safe_assignment_source":
                DOWNSTREAM_SAFE_ASSIGNMENT_SOURCE_CURRENT,
            "downstream_safe_change_applied": "FALSE",
            "downstream_assignment_status": DOWNSTREAM_ASSIGNMENT_HELD,
        }
    return {
        "downstream_safe_assignment": current,
        "downstream_safe_assignment_source":
            DOWNSTREAM_SAFE_ASSIGNMENT_SOURCE_CURRENT,
        "downstream_safe_change_applied": "FALSE",
        "downstream_assignment_status": DOWNSTREAM_ASSIGNMENT_NO_CHANGE,
    }


def expected_library_context(expected_genotypes_by_library: Mapping[str, object]):
    """Derive canonical per-library singlet/component/composite metadata views.

    ``expected_genotypes_by_library`` may map each library either to a mapping
    keyed by genotype or to an iterable of genotype strings.  All parsing flows
    through ``canonical_genotype``/``donor_components`` so this does not create a
    second genotype interpretation path.
    """
    explicit_singlets: Dict[str, Set[str]] = {}
    donor_universe: Dict[str, Set[str]] = {}
    composites_by_component: Dict[str, Dict[str, List[str]]] = {}

    for lib, value in expected_genotypes_by_library.items():
        genotypes = value.keys() if isinstance(value, Mapping) else value
        singlets: Set[str] = set()
        donors: Set[str] = set()
        composites = {}
        for raw in genotypes:
            genotype = canonical_genotype(raw)
            comps = donor_components(genotype)
            if not comps:
                continue
            donors.update(comps)
            if len(comps) == 1:
                singlets.add(comps[0])
            else:
                for donor in set(comps):
                    composites.setdefault(donor, set()).add(genotype)
        explicit_singlets[lib] = singlets
        donor_universe[lib] = donors
        composites_by_component[lib] = {
            donor: sorted(values, key=natural_key)
            for donor, values in composites.items()
        }
    return explicit_singlets, donor_universe, composites_by_component


def classify_singlet_library_relationship(
    genotype: str,
    explicit_expected_singlets: Iterable[str],
    expected_donor_components: Iterable[str],
    expected_composites_by_component: Mapping[str, Sequence[str]],
    globally_valid: bool = True,
) -> Tuple[str, str]:
    """Classify one single-donor identity relative to a library's biology."""
    comps = donor_components(genotype)
    if len(comps) != 1:
        return "NOT_SINGLET", ""
    donor = comps[0]
    if not globally_valid:
        return "UNRESOLVED", ""
    singlets = set(explicit_expected_singlets)
    components = set(expected_donor_components)
    if donor in singlets:
        return "EXPECTED_SINGLET", ""
    if donor in components:
        context = sorted(
            {canonical_genotype(x) for x in expected_composites_by_component.get(donor, []) if canonical_genotype(x)},
            key=natural_key,
        )
        return "EXPECTED_COMPOSITE_COMPONENT_SINGLET", ";".join(context)
    return "UNEXPECTED_SINGLET", ""


def pairwise_library_roster(
    library_a: str, library_b: str, expected_donor_components_by_library: Mapping[str, Iterable[str]]
) -> Dict[str, object]:
    """Return transparent donor-roster overlap and discriminability for a pair."""
    a = set(expected_donor_components_by_library.get(library_a, set()))
    b = set(expected_donor_components_by_library.get(library_b, set()))
    shared = a & b
    a_specific = a - b
    b_specific = b - a
    n_distinguishing = len(a_specific) + len(b_specific)
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
        if n_distinguishing <= 2:
            discriminability = "WEAK"
        elif n_distinguishing <= 5:
            discriminability = "MODERATE"
        else:
            discriminability = "STRONG"
    return {
        "roster_relation": relation,
        "pair_discriminability": discriminability,
        "shared": shared,
        "a_specific": a_specific,
        "b_specific": b_specific,
    }


def biological_state(genotype: str, ploidy: str = "") -> str:
    comps = donor_components(genotype)
    if not comps:
        return "UNKNOWN"
    pl = clean(ploidy).upper()
    if len(comps) == 1:
        if pl in {"T", "TET", "TETRAPLOID", "HOMOTYPIC"}:
            return f"T[{comps[0]}+{comps[0]}]"
        return f"D[{comps[0]}]"
    if len(comps) == 2:
        return f"T[{'+'.join(comps)}]"
    return "M{" + "+".join(comps) + "}"


def parse_library(value) -> str:
    s = clean(value)
    m = re.search(r"(?:lib|RNA[_-]?)(\d+)$", s, flags=re.I)
    if m:
        return f"lib{int(m.group(1))}"
    if re.fullmatch(r"\d+", s):
        return f"lib{int(s)}"
    return s


def library_number(value) -> Optional[int]:
    lib = parse_library(value)
    m = re.fullmatch(r"lib(\d+)", lib)
    return int(m.group(1)) if m else None


def parse_library_spec(values: Sequence[str]) -> List[int]:
    out = set()
    for item in values:
        for token in str(item).replace(",", " ").split():
            if "-" in token:
                a, b = token.split("-", 1)
                out.update(range(int(a), int(b) + 1))
            else:
                n = library_number(token)
                out.add(n if n is not None else int(token))
    return sorted(out)


def open_text(path: str, mode: str = "rt"):
    return gzip.open(path, mode, newline="") if str(path).endswith(".gz") else open(path, mode, newline="")


def read_tsv(path: str) -> List[Dict[str, str]]:
    with open_text(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [dict(row) for row in reader]


def iter_tsv(path: str) -> Iterator[Dict[str, str]]:
    with open_text(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            yield dict(row)


def write_tsv(path: str, rows: Iterable[Mapping[str, object]], fields: Sequence[str]) -> None:
    path = str(path)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = path + ".tmp"
    opener = gzip.open if path.endswith(".gz") else open
    kwargs = {"newline": ""}
    if path.endswith(".gz"):
        kwargs["mode"] = "wt"
    else:
        kwargs["mode"] = "w"
    with opener(tmp, **kwargs) as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fields), delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field, "")) for field in fields})
    os.replace(tmp, path)


def write_headerless_tsv(path: str, rows: Iterable[Sequence[object]]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = str(path) + ".tmp"
    with open(tmp, "w", newline="") as fh:
        for row in rows:
            fh.write("\t".join(format_value(x) for x in row) + "\n")
    os.replace(tmp, path)


def format_value(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.10g}"
    s = str(value)
    return s if s else "NA"


def ffloat(value, default: float = math.nan) -> float:
    try:
        x = float(clean(value))
        return x if math.isfinite(x) else default
    except Exception:
        return default


def fint(value, default: int = 0) -> int:
    try:
        return int(float(clean(value)))
    except Exception:
        return default


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def stable_site_fold(tid: int, pos: int, folds: int) -> int:
    if folds <= 1:
        return 0
    payload = f"{tid}:{pos}".encode()
    digest = hashlib.blake2b(payload, digest_size=8, person=b"cellbounc").digest()
    return int.from_bytes(digest, "little") % folds


def read_assignments(path: str) -> Dict[str, Dict[str, str]]:
    """Read the production headerless four-column .assignments convention."""
    out: Dict[str, Dict[str, str]] = {}
    with open_text(path, "rt") as fh:
        for line_no, raw in enumerate(fh, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n\r").split("\t")
            if len(parts) < 2:
                continue
            bc = clean(parts[0])
            out[bc] = {
                "barcode": bc,
                "assignment": canonical_genotype(parts[1]),
                "type": clean(parts[2]) if len(parts) > 2 else "",
                "score": clean(parts[3]) if len(parts) > 3 else "",
                "source_line": str(line_no),
            }
    return out


def read_refined_assignments(path: str) -> Dict[str, Dict[str, str]]:
    if not path or not os.path.isfile(path):
        return {}
    rows = read_tsv(path)
    out = {}
    for row in rows:
        bc = clean(row.get("barcode"))
        if bc:
            out[bc] = row
    return out


def load_alias_table(path: str) -> Dict[str, str]:
    if not path or not os.path.isfile(path):
        return {}
    aliases: Dict[str, str] = {}
    for row in read_tsv(path):
        src = clean(row.get("library_label") or row.get("alias") or row.get("source"))
        dst = clean(row.get("canonical_vcf_id") or row.get("canonical") or row.get("target"))
        if src and dst:
            aliases[src] = dst
        if dst:
            aliases.setdefault(dst, dst)
        for eq in re.split(r"[,;|]", clean(row.get("equivalent_or_reporter_vcf_ids"))):
            eq = clean(eq)
            if eq and dst:
                aliases.setdefault(eq, dst)
    return aliases


def choose_first(row: Mapping[str, object], names: Sequence[str], default: str = "") -> str:
    normalized = {re.sub(r"[^a-z0-9]", "", str(k).lower()): k for k in row}
    for name in names:
        key = normalized.get(re.sub(r"[^a-z0-9]", "", name.lower()))
        if key is not None:
            val = clean(row.get(key))
            if val:
                return val
    return default


def json_dump_atomic(path: str, obj) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, path)
