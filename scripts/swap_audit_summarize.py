#!/usr/bin/env python3
"""Summarize one library's post-hoc CellBouncer swap audit."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Mapping, Tuple

from swap_audit_lib import (
    NA,
    canonical_identity,
    identity_components,
    identity_type,
    median,
    mean,
    normalize_library,
    parse_expected_metadata,
    percentile,
    read_panel_metadata,
    read_refined_assignments,
    refined_assignment_counts,
    composition_distances_raw_identity,
    snp_resolvable_identity,
    species_relation,
    species_relation_booleans,
    species_set,
    write_tsv,
)

DEFAULT_THRESHOLDS = {
    "threshold_min_signal_llr": 10.0,
    "threshold_max_ambiguous_frac": 0.30,
    "threshold_min_cells": 100.0,
    "threshold_clean_overlap_frac": 0.90,
    "threshold_clean_best_frac": 0.30,
    "threshold_full_swap_frac": 0.50,
    "threshold_one_component_swap_frac": 0.50,
    "threshold_coherent_alt_frac": 0.25,
    "threshold_ambiguous_delta_ll": 20.0,
    "threshold_homotypic_frac": 0.50,
    "threshold_species_conflict_frac": 0.20,
    "threshold_species_unexpected_frac": 0.20,
    "threshold_species_disjoint_frac": 0.02,
    "threshold_species_component_missing_frac": 0.20,
    "threshold_min_dosage_concordance": 0.70,
    "threshold_max_low_evidence_frac": 0.30,
    "threshold_ambiguous_gap": 0.0,
    "threshold_clean_concordance": 0.70,
}

REPORT_FIELDS = [
    "library", "preflight_status", "audit_verdict", "audit_flags", "feature_mode",
    "has_call_qc", "has_species_qc", "has_atac_qc", "has_refined_assignments", "has_calibrated_thresholds",
    "missing_optional_features", "warnings", "thresholds_source", "n_cells", "expected_identity",
    "audit_best_identity_unconstrained", "audit_best_identity_constrained", "audit_best_fraction",
    "expected_rank_median", "expected_rank_p90", "delta_ll_best_vs_expected_median",
    "frac_cells_component_overlap_2", "frac_cells_component_overlap_1", "frac_cells_component_overlap_0",
    "median_dosage_concordance_assigned", "median_dosage_concordance_unconstrained",
    "median_dosage_gap_constrained", "median_dosage_gap_unconstrained",
    "frac_cells_neg_gap_constrained", "frac_cells_neg_gap_unconstrained",
    "expected_vs_observed_l1", "expected_vs_observed_cosine", "jensen_shannon_distance",
    "unexpected_identity_fraction", "missing_expected_identity_fraction",
    "refined_n_cells", "refined_n_total_available", "refined_n_overlap_with_audit", "refined_overlap_fraction",
    "refined_n_changed", "refined_n_singlet", "refined_n_heterotypic",
    "refined_n_homotypic", "refined_n_plus_identity", "refined_n_biological_fusion",
    "refined_n_doublet", "refined_n_droplet_flagged",
    "refined_top_identity", "refined_top_fraction",
    "refined_biological_unexpected_fraction", "refined_biological_missing_expected_fraction",
    "species_expected",
    "species_audit_best", "median_species_support_expected", "frac_cells_species_conflict",
    "frac_cells_species_exact_match", "frac_cells_species_component_missing",
    "frac_cells_species_unexpected_extra", "frac_cells_species_partial_overlap_extra",
    "frac_cells_species_disjoint_wrong", "frac_cells_species_unexpected_or_disjoint",
    "unresolved_homotypic_ratio", "atac_qc_mode", "median_atac_dosage_concordance",
    "atac_best_identity", "rna_atac_same_identity_fraction", "rna_atac_same_species_fraction",
    "frac_cells_rna_atac_discordant",
]

SCORE_FIELDS = [
    "barcode", "libname", "raw_demux_assignment", "snp_resolvable_assignment",
    "refined_biological_assignment", "refined_assignment_source",
    "refined_assignment_confidence", "droplet_doublet_flag", "quad_pattern_score",
    "original_assignment", "audit_unconstrained_best", "audit_constrained_best",
    "expected_components", "audit_best_components", "shared_components", "missing_expected_components",
    "unexpected_components", "component_overlap", "expected_rank_unconstrained", "expected_rank_constrained",
    "delta_ll_best_vs_expected", "dosage_concordance_assigned", "dosage_concordance_unconstrained",
    "dosage_gap_constrained", "dosage_gap_unconstrained", "species_support_expected",
    "species_relation", "species_missing_expected_component",
    "species_has_unexpected_component", "species_disjoint_wrong_species",
    "combined_qc_flags", "swap_cell_verdict", "feature_mode", "warnings",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--lib", required=True, help="Library number/name, e.g. 17 or lib17")
    p.add_argument("--capabilities_manifest", required=True)
    p.add_argument("--audit_root", required=True)
    p.add_argument("--expected_metadata", required=True)
    p.add_argument("--panel_metadata", required=True)
    p.add_argument("--thresholds", default=None)
    p.add_argument("--out_prefix", default=None, help="Default: <audit_root>/<lib>/<lib>")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--skip-existing", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--strict", action="store_true")
    p.add_argument("--best-effort", action="store_true")
    return p.parse_args()


def load_thresholds(path: str | None) -> Tuple[Dict[str, float], str, bool]:
    t = dict(DEFAULT_THRESHOLDS)
    if not path:
        return t, "defaults_uncalibrated", False
    p = Path(path)
    if not p.exists():
        return t, "defaults_uncalibrated", False
    with open(p) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2 or parts[0] == "threshold_name":
                continue
            try:
                t[parts[0]] = float(parts[1])
            except ValueError:
                pass
    return t, "calibrated", True


def load_assignments(path: str) -> Dict[str, Dict[str, object]]:
    out: Dict[str, Dict[str, object]] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            bc = parts[0]
            ident = canonical_identity(parts[1], None)
            sd = parts[2] if len(parts) > 2 else ("D" if "+" in ident else "S")
            llr = float(parts[3]) if len(parts) > 3 and parts[3] not in {"", ".", "NA"} else math.nan
            out[bc] = {"assignment": ident, "singlet_doublet": sd, "llr": llr}
    return out


def load_diagnostics(path: str) -> Dict[str, Dict[str, str]]:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return {r["barcode"]: r for r in csv.DictReader(fh, delimiter="\t")}


def load_runner_ups(path: str) -> Dict[str, List[Dict[str, object]]]:
    out: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            bc = r.get("barcode", "")
            if not bc:
                continue
            ident = r.get("identity", "")
            if not ident:
                continue
            try:
                ident = canonical_identity(ident, None)
            except ValueError:
                continue
            rank = int(float(r.get("rank", "0") or 0))
            margin = r.get("min_margin", r.get("llr_vs_winner", "NA"))
            try:
                margin_f = float(margin)
            except ValueError:
                margin_f = math.nan
            out[bc].append({"rank": rank, "identity": ident, "min_margin": margin_f})
    for bc in out:
        out[bc].sort(key=lambda x: int(x["rank"]))
    return out


def load_call_qc(path: str | None) -> Dict[str, Dict[str, str]]:
    if not path or not Path(path).exists():
        return {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return {r["barcode"]: r for r in csv.DictReader(fh, delimiter="\t")}


def load_species_qc(path: str | None) -> Dict[str, str]:
    if not path or not Path(path).exists():
        return {}
    rows = []
    with open(path, "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows[0] if rows else {}


def path_exists(path: str | None) -> bool:
    return bool(path) and Path(str(path)).exists()


def runtime_feature_mode(has_call_qc: bool, has_species_qc: bool, has_atac_qc: bool) -> str:
    if has_call_qc and has_species_qc and has_atac_qc:
        return "CORE_PLUS_CALL_QC_AND_SPECIES_AND_ATAC"
    if has_call_qc and has_species_qc:
        return "CORE_PLUS_CALL_QC_AND_SPECIES"
    if has_call_qc:
        return "CORE_PLUS_CALL_QC"
    if has_species_qc:
        return "CORE_PLUS_SPECIES"
    return "CORE_ONLY"


def safe_float(x, default=math.nan):
    try:
        if x in (None, "", ".", "NA"):
            return default
        return float(x)
    except (TypeError, ValueError):
        return default


def top_identity(assignments: Mapping[str, Mapping[str, object]]) -> Tuple[str, float]:
    c = Counter(str(v["assignment"]) for v in assignments.values())
    if not c:
        return NA, math.nan
    ident, n = c.most_common(1)[0]
    return ident, n / sum(c.values())


def components_csv(ident: str) -> str:
    return "+".join(identity_components(ident))


def runner_expected_rank(runners: List[Mapping[str, object]], expected: str) -> Tuple[float, float]:
    expected_snp = snp_resolvable_identity(expected, None)
    for r in runners:
        try:
            runner_snp = snp_resolvable_identity(str(r["identity"]), None)
        except ValueError:
            continue
        if runner_snp == expected_snp:
            return float(r["rank"]), float(r.get("min_margin", math.nan))
    return math.nan, math.nan


def composition_distances(expected_ids: List[str], observed_assignments: Mapping[str, Mapping[str, object]]) -> Dict[str, float]:
    # SNP evidence cannot distinguish A from A+A.  Collapse only for
    # composition/audit-distance math; raw biological labels remain available in
    # expected_identity, switch_matrix, and tetra_refine outputs.
    obs_counts = Counter(snp_resolvable_identity(str(v["assignment"]), None)
                         for v in observed_assignments.values())
    total = sum(obs_counts.values()) or 1
    obs = {k: v / total for k, v in obs_counts.items()}
    expected_snp_ids = sorted({snp_resolvable_identity(x, None) for x in expected_ids})
    exp = {k: 1.0 / len(expected_snp_ids) for k in expected_snp_ids} if expected_snp_ids else {}
    keys = sorted(set(obs) | set(exp))
    l1 = sum(abs(obs.get(k, 0.0) - exp.get(k, 0.0)) for k in keys)
    dot = sum(obs.get(k, 0.0) * exp.get(k, 0.0) for k in keys)
    norm_o = math.sqrt(sum(obs.get(k, 0.0) ** 2 for k in keys))
    norm_e = math.sqrt(sum(exp.get(k, 0.0) ** 2 for k in keys))
    cosine = 1.0 - dot / (norm_o * norm_e) if norm_o and norm_e else math.nan
    def kl(p, q):
        s = 0.0
        for k in keys:
            if p.get(k, 0.0) > 0 and q.get(k, 0.0) > 0:
                s += p[k] * math.log(p[k] / q[k], 2)
        return s
    m = {k: 0.5 * (obs.get(k, 0.0) + exp.get(k, 0.0)) for k in keys}
    js = math.sqrt(0.5 * kl(obs, m) + 0.5 * kl(exp, m)) if keys else math.nan
    unexpected = sum(v for k, v in obs.items() if k not in exp)
    missing = sum(1 for k in exp if obs.get(k, 0.0) < 0.005) / len(exp) if exp else math.nan
    return {
        "expected_vs_observed_l1": l1,
        "expected_vs_observed_cosine": cosine,
        "jensen_shannon_distance": js,
        "unexpected_identity_fraction": unexpected,
        "missing_expected_identity_fraction": missing,
    }


def verdict_logic(row: Dict[str, object], t: Mapping[str, float]) -> Tuple[str, str]:
    flags: List[str] = []
    n_cells = float(row.get("n_cells", 0) or 0)
    if n_cells < t["threshold_min_cells"]:
        flags.append("LOW_CELL_COUNT")
    if float(row.get("median_original_llr", math.nan)) < t["threshold_min_signal_llr"]:
        flags.append("LOW_SIGNAL_LLR")
    if float(row.get("frac_cells_high_n_close", 0.0) or 0.0) > t["threshold_max_ambiguous_frac"]:
        flags.append("HIGH_AMBIGUOUS_N_CLOSE")
    frac_species_disjoint = float(row.get("frac_cells_species_disjoint_wrong", 0.0) or 0.0)
    frac_species_unexpected = float(row.get("frac_cells_species_unexpected_or_disjoint", 0.0) or 0.0)
    frac_species_component_missing = float(row.get("frac_cells_species_component_missing", 0.0) or 0.0)
    if frac_species_disjoint >= t.get("threshold_species_disjoint_frac", 0.02):
        flags.append("WRONG_SPECIES_SIGNAL")
    elif frac_species_unexpected >= t.get("threshold_species_unexpected_frac", t["threshold_species_conflict_frac"]):
        flags.append("UNEXPECTED_SPECIES_SIGNAL")
    if frac_species_component_missing >= t.get("threshold_species_component_missing_frac", t["threshold_species_conflict_frac"]):
        flags.append("SPECIES_COMPONENT_IMBALANCE")
    if float(row.get("frac_cells_component_overlap_2", 0.0) or 0.0) >= t["threshold_clean_overlap_frac"] and float(row.get("audit_best_fraction", 0.0) or 0.0) >= t["threshold_clean_best_frac"]:
        clean = True
    else:
        clean = False
    if float(row.get("frac_cells_component_overlap_0", 0.0) or 0.0) >= t["threshold_full_swap_frac"]:
        flags.append("LIKELY_FULL_FUSION_SWAP")
    if float(row.get("frac_cells_component_overlap_1", 0.0) or 0.0) >= t["threshold_one_component_swap_frac"] and float(row.get("audit_best_fraction", 0.0) or 0.0) >= t["threshold_coherent_alt_frac"]:
        flags.append("LIKELY_ONE_COMPONENT_SWAP")
    delta = float(row.get("delta_ll_best_vs_expected_median", math.nan))
    if not math.isnan(delta) and abs(delta) < t["threshold_ambiguous_delta_ll"]:
        flags.append("AMBIGUOUS_CLOSE_GENOTYPE")
    if float(row.get("unresolved_homotypic_ratio", 0.0) or 0.0) >= t["threshold_homotypic_frac"]:
        flags.append("UNRESOLVED_HOMOTYPIC_OR_SINGLET")
    if float(row.get("unexpected_identity_fraction", 0.0) or 0.0) > 0.05:
        flags.append("UNEXPECTED_EXTRA_IDENTITY")
    if float(row.get("missing_expected_identity_fraction", 0.0) or 0.0) > 0.05:
        flags.append("EXPECTED_IDENTITY_MISSING")

    if "LOW_CELL_COUNT" in flags or "LOW_SIGNAL_LLR" in flags or "HIGH_AMBIGUOUS_N_CLOSE" in flags:
        verdict = "LOW_SIGNAL"
    elif "WRONG_SPECIES_SIGNAL" in flags:
        verdict = "WRONG_SPECIES_SIGNAL"
    elif "UNEXPECTED_SPECIES_SIGNAL" in flags:
        verdict = "UNEXPECTED_SPECIES_SIGNAL"
    elif "SPECIES_COMPONENT_IMBALANCE" in flags:
        verdict = "SPECIES_COMPONENT_IMBALANCE"
    elif clean and not {"LIKELY_FULL_FUSION_SWAP", "LIKELY_ONE_COMPONENT_SWAP"}.intersection(flags):
        verdict = "OK"
    elif "LIKELY_FULL_FUSION_SWAP" in flags:
        verdict = "LIKELY_FULL_FUSION_SWAP"
    elif "LIKELY_ONE_COMPONENT_SWAP" in flags:
        verdict = "LIKELY_ONE_COMPONENT_SWAP"
    elif "AMBIGUOUS_CLOSE_GENOTYPE" in flags:
        verdict = "AMBIGUOUS_CLOSE_GENOTYPE"
    elif "UNRESOLVED_HOMOTYPIC_OR_SINGLET" in flags:
        verdict = "UNRESOLVED_HOMOTYPIC_OR_SINGLET"
    elif "UNEXPECTED_EXTRA_IDENTITY" in flags:
        verdict = "UNEXPECTED_EXTRA_IDENTITY"
    elif "EXPECTED_IDENTITY_MISSING" in flags:
        verdict = "EXPECTED_IDENTITY_MISSING"
    else:
        verdict = "OK"
    return verdict, ",".join(sorted(set(flags)))


def main() -> int:
    args = parse_args()
    lib = normalize_library(args.lib)
    manifest = json.loads(Path(args.capabilities_manifest).read_text())
    audit_root = Path(args.audit_root).resolve()
    out_prefix = Path(args.out_prefix).resolve() if args.out_prefix else audit_root / lib / lib
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    thresholds_path = args.thresholds or manifest.get("optional_files", {}).get("thresholds")
    thresholds, thresholds_source, has_calibrated = load_thresholds(thresholds_path)
    expected = parse_expected_metadata(args.expected_metadata)
    expected_ids = expected.get(lib, manifest.get("expected_identities", []))
    panel = read_panel_metadata(args.panel_metadata)

    core = manifest["core_files"]
    audit_files = manifest["audit_files"]
    required = [
        core["assignments"], core["diagnostics"],
        audit_files["unconstrained_prefix"] + ".assignments",
        audit_files["unconstrained_prefix"] + ".runner_ups.gz",
        audit_files["constrained_prefix"] + ".assignments",
        audit_files["constrained_prefix"] + ".runner_ups.gz",
    ]
    missing = [p for p in required if not Path(p).exists()]
    if missing:
        row = {f: NA for f in REPORT_FIELDS}
        row.update({
            "library": lib, "preflight_status": "FAIL_MISSING_AUDIT_OUTPUT", "audit_verdict": "NOT_RUN",
            "audit_flags": "FAIL_MISSING_AUDIT_OUTPUT", "feature_mode": manifest.get("feature_mode", "CORE_ONLY"),
            "has_call_qc": False, "has_species_qc": False, "has_atac_qc": False, "has_refined_assignments": False,
            "has_calibrated_thresholds": has_calibrated, "missing_optional_features": ",".join(manifest.get("missing_optional_features", [])),
            "warnings": ",".join(sorted(set(manifest.get("warnings", []) + ["FAIL_MISSING_AUDIT_OUTPUT"]))),
            "thresholds_source": thresholds_source, "expected_identity": ",".join(expected_ids),
        })
        write_tsv(str(out_prefix) + ".swap_report.tsv", [row], REPORT_FIELDS)
        raise SystemExit("ERROR: missing audit outputs:\n" + "\n".join(missing))

    orig = load_assignments(core["assignments"])
    diag = load_diagnostics(core["diagnostics"])
    au = load_assignments(audit_files["unconstrained_prefix"] + ".assignments")
    ac = load_assignments(audit_files["constrained_prefix"] + ".assignments")
    ru_u = load_runner_ups(audit_files["unconstrained_prefix"] + ".runner_ups.gz")
    ru_c = load_runner_ups(audit_files["constrained_prefix"] + ".runner_ups.gz")
    opt_files = manifest.get("optional_files", {})
    call_qc = load_call_qc(opt_files.get("call_qc"))
    species_qc = load_species_qc(opt_files.get("species_qc"))
    refined_assignments = read_refined_assignments(opt_files.get("refined_assignments"))

    stale_optional_warning_tags = {
        "CALL_QC_MISSING_REDUCED_MODE",
        "SPECIES_QC_NOT_YET_COMPUTED",
        "SPECIES_NATIVE_BUNDLE_MISSING",
        "SPECIES_COUNTS_MISSING",
        "SPECIES_CONDF_MISSING",
        "USING_DEFAULT_THRESHOLDS",
    }
    warnings = sorted(set(manifest.get("warnings", [])) - stale_optional_warning_tags)
    has_call_qc = bool(call_qc)
    species_qc_warnings = set(x for x in (species_qc.get("warnings", "") or "").split(",") if x)
    has_species_qc = bool(species_qc) and safe_float(species_qc.get("n_cells_with_species_evidence"), 0.0) > 0 and not ({"NO_SPECIES_INPUTS", "SPECIES_SCORING_DISABLED", "SPECIES_COUNTS_MISSING", "SPECIES_SAMPLES_MISSING", "PANEL_METADATA_MISSING"} & species_qc_warnings)
    has_atac_qc = path_exists(opt_files.get("atac_counts")) and path_exists(opt_files.get("atac_call_qc"))
    has_refined_assignments = bool(refined_assignments)
    if not has_refined_assignments:
        warnings.append("REFINED_ASSIGNMENTS_MISSING_BIOLOGICAL_REPORTING_SKIPPED")
    if not has_call_qc and "CALL_QC_MISSING_REDUCED_MODE" not in warnings:
        warnings.append("CALL_QC_MISSING_REDUCED_MODE")
    if not has_species_qc:
        raw_species_bundle = all(path_exists(opt_files.get(k)) for k in ("species_counts", "species_condf", "species_samples"))
        warnings.append("SPECIES_QC_NOT_YET_COMPUTED" if raw_species_bundle else "SPECIES_NATIVE_BUNDLE_MISSING")
    runtime_missing = []
    if not has_call_qc:
        runtime_missing.append("call_qc")
    if not has_species_qc:
        runtime_missing.append("species_qc")
    if not has_atac_qc:
        runtime_missing.extend(["atac_counts", "atac_call_qc"])
    if not has_refined_assignments:
        runtime_missing.append("refined_assignments")
    if not has_calibrated:
        runtime_missing.append("thresholds")
        warnings.append("USING_DEFAULT_THRESHOLDS")
    warnings = sorted(set(warnings))
    runtime_missing = sorted(set(runtime_missing))
    feature_mode_runtime = runtime_feature_mode(has_call_qc, has_species_qc, has_atac_qc)

    barcodes = sorted(set(orig) & set(au))

    # Refined-assignment summaries must be restricted to the same barcode
    # universe used for the audit.  A stale or broader refined sidecar should not
    # skew lib-level biological composition while per-cell audit rows silently
    # fall back to raw demux assignments.
    refined_for_scored = {bc: refined_assignments[bc] for bc in barcodes if bc in refined_assignments}
    refined_total_available = len(refined_assignments)
    refined_overlap_n = len(refined_for_scored)
    refined_overlap_fraction = (refined_overlap_n / len(barcodes)) if barcodes else math.nan
    if has_refined_assignments and refined_overlap_n == 0:
        warnings.append("REFINED_ASSIGNMENTS_NO_AUDIT_BARCODE_OVERLAP")
    elif has_refined_assignments and refined_overlap_fraction < 0.95:
        warnings.append("REFINED_ASSIGNMENTS_BARCODE_OVERLAP_LOW")
    warnings = sorted(set(warnings))

    score_rows: List[Dict[str, object]] = []
    overlap_counts = Counter()
    ranks_u: List[float] = []
    ranks_c: List[float] = []
    deltas: List[float] = []
    unresolved_homotypic = 0
    species_conflicts = 0
    species_relation_counts = Counter()

    for bc in barcodes:
        o = str(orig[bc]["assignment"])
        ref_row = refined_for_scored.get(bc, {})
        refined_bio = str(ref_row.get("refined_assignment", o))
        refined_source = str(ref_row.get("refined_assignment_source", ref_row.get("ploidy_method", "demux_only")))
        refined_conf = str(ref_row.get("refined_assignment_confidence", ref_row.get("overall_confidence", NA)))
        droplet_flag = str(ref_row.get("droplet_doublet_flag", "NONE"))
        quad_score = str(ref_row.get("quad_pattern_score", NA))
        snp_o = snp_resolvable_identity(o, None)
        ub = str(au[bc]["assignment"])
        cb = str(ac.get(bc, {}).get("assignment", NA))
        ocomp = set(identity_components(o))
        ucomp = set(identity_components(ub)) if ub != NA else set()
        shared = sorted(ocomp & ucomp)
        missing_exp = sorted(ocomp - ucomp)
        unexpected = sorted(ucomp - ocomp)
        overlap = len(shared)
        overlap_counts[overlap] += 1
        if identity_type(o) == "homotypic" or identity_type(ub) == "homotypic":
            unresolved_homotypic += 1
        rank_u, delta_u = runner_expected_rank(ru_u.get(bc, []), o)
        rank_c, _ = runner_expected_rank(ru_c.get(bc, []), o)
        if not math.isnan(rank_u): ranks_u.append(rank_u)
        if not math.isnan(rank_c): ranks_c.append(rank_c)
        if not math.isnan(delta_u): deltas.append(abs(delta_u))
        qc = call_qc.get(bc, {})
        sp_exp = species_set(o, panel)
        sp_best = species_set(ub, panel) if ub != NA else NA
        sp_rel = species_relation(sp_exp, sp_best)
        sp_missing, sp_unexpected, sp_disjoint = species_relation_booleans(sp_exp, sp_best)
        if sp_rel not in {"exact_match", "missing_species_evidence"}:
            species_conflicts += 1
        species_relation_counts[sp_rel] += 1
        flags: List[str] = []
        if overlap == 0:
            flags.append("FULL_SWAP")
            cell_verdict = "FULL_SWAP"
        elif overlap == 1 and len(ocomp) > 1:
            flags.append("ONE_COMPONENT_SWAP")
            cell_verdict = "ONE_COMPONENT_SWAP"
        else:
            cell_verdict = "OK"
        if qc.get("call_qc_flags"):
            flags.extend([x for x in qc.get("call_qc_flags", "").split(",") if x])
        score_rows.append({
            "barcode": bc,
            "libname": lib,
            "raw_demux_assignment": o,
            "snp_resolvable_assignment": snp_o,
            "refined_biological_assignment": refined_bio,
            "refined_assignment_source": refined_source,
            "refined_assignment_confidence": refined_conf,
            "droplet_doublet_flag": droplet_flag,
            "quad_pattern_score": quad_score,
            "original_assignment": o,
            "audit_unconstrained_best": ub,
            "audit_constrained_best": cb,
            "expected_components": "+".join(sorted(ocomp)),
            "audit_best_components": "+".join(sorted(ucomp)),
            "shared_components": "+".join(shared),
            "missing_expected_components": "+".join(missing_exp),
            "unexpected_components": "+".join(unexpected),
            "component_overlap": overlap,
            "expected_rank_unconstrained": rank_u,
            "expected_rank_constrained": rank_c,
            "delta_ll_best_vs_expected": delta_u,
            "dosage_concordance_assigned": qc.get("dosage_concordance", NA),
            "dosage_concordance_unconstrained": NA,
            "dosage_gap_constrained": qc.get("dosage_gap_constrained", NA),
            "dosage_gap_unconstrained": NA,
            "species_support_expected": qc.get("species_support_expected", NA),
            "species_relation": qc.get("species_relation", sp_rel),
            "species_missing_expected_component": qc.get("species_missing_expected_component", str(sp_missing)),
            "species_has_unexpected_component": qc.get("species_has_unexpected_component", str(sp_unexpected)),
            "species_disjoint_wrong_species": qc.get("species_disjoint_wrong_species", str(sp_disjoint)),
            "combined_qc_flags": sorted(set(flags)),
            "swap_cell_verdict": cell_verdict,
            "feature_mode": feature_mode_runtime,
            "warnings": warnings,
        })

    n = len(barcodes)
    top_u, top_u_frac = top_identity({bc: au[bc] for bc in barcodes if bc in au})
    top_c, _ = top_identity({bc: ac[bc] for bc in barcodes if bc in ac})
    llrs = []
    n_close_vals = []
    for bc in barcodes:
        d = diag.get(bc, {})
        try: llrs.append(float(d.get("llr", math.nan)))
        except ValueError: pass
        try: n_close_vals.append(float(d.get("n_close", 0)))
        except ValueError: pass
    comp = composition_distances(expected_ids, orig)
    refined_counts = refined_assignment_counts(refined_for_scored) if refined_for_scored else Counter()
    refined_top_identity, refined_top_fraction = (NA, math.nan)
    if refined_for_scored:
        rc = Counter(str(r.get("refined_assignment", NA)) for r in refined_for_scored.values() if r.get("refined_assignment", NA) != NA)
        if rc:
            refined_top_identity, refined_top_n = rc.most_common(1)[0]
            refined_top_fraction = refined_top_n / sum(rc.values())
    refined_comp = composition_distances_raw_identity(
        expected_ids,
        (str(r.get("refined_assignment", NA)) for r in refined_for_scored.values())
    ) if refined_for_scored else {"unexpected_identity_fraction": math.nan, "missing_expected_identity_fraction": math.nan}
    qcdos = []
    qcgap = []
    qclow = 0
    qcspecies = []
    qc_species_conflicts = 0
    qc_species_den = 0
    qc_species_relation_counts = Counter()
    for r in call_qc.values():
        v = safe_float(r.get("dosage_concordance"))
        if not math.isnan(v): qcdos.append(v)
        v = safe_float(r.get("dosage_gap_constrained"))
        if not math.isnan(v): qcgap.append(v)
        flags = r.get("call_qc_flags", "")
        if "LOW_EVIDENCE" in flags: qclow += 1
        v = safe_float(r.get("species_support_expected"))
        if not math.isnan(v): qcspecies.append(v)
        rel = str(r.get("species_relation", ""))
        if rel:
            qc_species_relation_counts[rel] += 1
        scf = str(r.get("species_conflict_flag", "NA"))
        if scf in {"0", "1"}:
            qc_species_den += 1
            if scf == "1": qc_species_conflicts += 1

    species_expected = ";".join(sorted({species_set(x, panel) for x in expected_ids if species_set(x, panel) != NA})) or NA
    species_audit_best = species_set(top_u, panel) if top_u != NA else NA

    def relation_frac_from_sidecar(key: str, rel_counts: Counter, denom: int) -> float:
        if has_species_qc and species_qc.get(key) not in (None, "", NA):
            return safe_float(species_qc.get(key))
        return (rel_counts.get(key.replace("frac_cells_species_", ""), 0) / denom) if denom else math.nan

    rel_counts = qc_species_relation_counts if qc_species_relation_counts else species_relation_counts
    rel_denom = sum(rel_counts.values()) if rel_counts else (qc_species_den if qc_species_den else n)
    exact_frac = safe_float(species_qc.get("frac_cells_species_exact_match")) if has_species_qc and species_qc.get("frac_cells_species_exact_match") not in (None, "", NA) else (rel_counts.get("exact_match", 0) / rel_denom if rel_denom else math.nan)
    component_missing_frac = safe_float(species_qc.get("frac_cells_species_component_missing")) if has_species_qc and species_qc.get("frac_cells_species_component_missing") not in (None, "", NA) else (rel_counts.get("expected_subset_only_component_missing", 0) / rel_denom if rel_denom else math.nan)
    unexpected_extra_frac = safe_float(species_qc.get("frac_cells_species_unexpected_extra")) if has_species_qc and species_qc.get("frac_cells_species_unexpected_extra") not in (None, "", NA) else (rel_counts.get("expected_superset_with_extra_species", 0) / rel_denom if rel_denom else math.nan)
    partial_overlap_frac = safe_float(species_qc.get("frac_cells_species_partial_overlap_extra")) if has_species_qc and species_qc.get("frac_cells_species_partial_overlap_extra") not in (None, "", NA) else (rel_counts.get("partial_overlap_with_extra_and_missing", 0) / rel_denom if rel_denom else math.nan)
    disjoint_wrong_frac = safe_float(species_qc.get("frac_cells_species_disjoint_wrong")) if has_species_qc and species_qc.get("frac_cells_species_disjoint_wrong") not in (None, "", NA) else (rel_counts.get("disjoint_wrong_species", 0) / rel_denom if rel_denom else math.nan)
    unexpected_or_disjoint_frac = safe_float(species_qc.get("frac_cells_species_unexpected_or_disjoint")) if has_species_qc and species_qc.get("frac_cells_species_unexpected_or_disjoint") not in (None, "", NA) else (sum(rel_counts.get(k, 0) for k in ("expected_superset_with_extra_species", "partial_overlap_with_extra_and_missing", "disjoint_wrong_species")) / rel_denom if rel_denom else math.nan)

    row = {f: NA for f in REPORT_FIELDS}
    row.update({
        "library": lib,
        "preflight_status": "PASS",
        "feature_mode": feature_mode_runtime,
        "has_call_qc": has_call_qc,
        "has_species_qc": has_species_qc,
        "has_atac_qc": has_atac_qc,
        "has_refined_assignments": has_refined_assignments,
        "has_calibrated_thresholds": has_calibrated,
        "missing_optional_features": ",".join(runtime_missing),
        "warnings": warnings,
        "thresholds_source": thresholds_source,
        "n_cells": n,
        "expected_identity": ",".join(expected_ids),
        "audit_best_identity_unconstrained": top_u,
        "audit_best_identity_constrained": top_c,
        "audit_best_fraction": top_u_frac,
        "expected_rank_median": median(ranks_u),
        "expected_rank_p90": percentile(ranks_u, 90),
        "delta_ll_best_vs_expected_median": median(deltas),
        "frac_cells_component_overlap_2": overlap_counts[2] / n if n else math.nan,
        "frac_cells_component_overlap_1": overlap_counts[1] / n if n else math.nan,
        "frac_cells_component_overlap_0": overlap_counts[0] / n if n else math.nan,
        "median_dosage_concordance_assigned": median(qcdos),
        "median_dosage_concordance_unconstrained": NA,
        "median_dosage_gap_constrained": median(qcgap),
        "median_dosage_gap_unconstrained": NA,
        "frac_cells_neg_gap_constrained": sum(1 for v in qcgap if v < 0) / len(qcgap) if qcgap else math.nan,
        "frac_cells_neg_gap_unconstrained": NA,
        "refined_n_cells": refined_counts.get("total", NA) if refined_for_scored else NA,
        "refined_n_total_available": refined_total_available if has_refined_assignments else NA,
        "refined_n_overlap_with_audit": refined_overlap_n if has_refined_assignments else NA,
        "refined_overlap_fraction": refined_overlap_fraction if has_refined_assignments else NA,
        "refined_n_changed": refined_counts.get("changed", NA) if refined_for_scored else NA,
        "refined_n_singlet": refined_counts.get("singlet", NA) if refined_for_scored else NA,
        "refined_n_heterotypic": refined_counts.get("heterotypic", NA) if refined_for_scored else NA,
        "refined_n_homotypic": refined_counts.get("homotypic", NA) if refined_for_scored else NA,
        "refined_n_plus_identity": refined_counts.get("plus_identity", NA) if refined_for_scored else NA,
        "refined_n_biological_fusion": refined_counts.get("biological_fusion", NA) if refined_for_scored else NA,
        "refined_n_doublet": refined_counts.get("doublet", NA) if refined_for_scored else NA,
        "refined_n_droplet_flagged": refined_counts.get("droplet_flagged", NA) if refined_for_scored else NA,
        "refined_top_identity": refined_top_identity,
        "refined_top_fraction": refined_top_fraction,
        "refined_biological_unexpected_fraction": refined_comp.get("unexpected_identity_fraction", math.nan),
        "refined_biological_missing_expected_fraction": refined_comp.get("missing_expected_identity_fraction", math.nan),
        "species_expected": species_expected,
        "species_audit_best": species_audit_best,
        "median_species_support_expected": safe_float(species_qc.get("median_species_support_expected")) if has_species_qc else median(qcspecies),
        "frac_cells_species_conflict": safe_float(species_qc.get("frac_cells_species_conflict")) if has_species_qc else (qc_species_conflicts / qc_species_den if qc_species_den else (species_conflicts / n if n else math.nan)),
        "frac_cells_species_exact_match": exact_frac,
        "frac_cells_species_component_missing": component_missing_frac,
        "frac_cells_species_unexpected_extra": unexpected_extra_frac,
        "frac_cells_species_partial_overlap_extra": partial_overlap_frac,
        "frac_cells_species_disjoint_wrong": disjoint_wrong_frac,
        "frac_cells_species_unexpected_or_disjoint": unexpected_or_disjoint_frac,
        "unresolved_homotypic_ratio": unresolved_homotypic / n if n else math.nan,
        "atac_qc_mode": "ATAC_PRESENT" if has_atac_qc else "ATAC_NOT_PROVIDED",
        "median_atac_dosage_concordance": NA,
        "atac_best_identity": NA,
        "rna_atac_same_identity_fraction": NA,
        "rna_atac_same_species_fraction": NA,
        "frac_cells_rna_atac_discordant": NA,
        "median_original_llr": median(llrs),
        "frac_cells_high_n_close": sum(1 for v in n_close_vals if v > 0) / len(n_close_vals) if n_close_vals else 0.0,
    })
    row.update(comp)
    if has_call_qc:
        row["frac_cells_low_evidence"] = qclow / len(call_qc) if call_qc else math.nan
    verdict, flags = verdict_logic(row, thresholds)
    row["audit_verdict"] = verdict
    row["audit_flags"] = flags

    write_tsv(str(out_prefix) + ".swap_scores.tsv", score_rows, SCORE_FIELDS)
    write_tsv(str(out_prefix) + ".swap_report.tsv", [row], REPORT_FIELDS)

    # Switch matrix and identity fractions.
    row_ids = sorted(set(str(orig[bc]["assignment"]) for bc in barcodes if bc in orig))
    col_ids = sorted(set(str(au[bc]["assignment"]) for bc in barcodes if bc in au))
    matrix_counts = defaultdict(Counter)
    for bc in barcodes:
        matrix_counts[str(orig[bc]["assignment"])][str(au[bc]["assignment"])] += 1
    matrix_rows = []
    for rid in row_ids:
        mr = {"expected_identity": rid}
        for cid in col_ids:
            mr[cid] = matrix_counts[rid][cid]
        matrix_rows.append(mr)
    write_tsv(str(out_prefix) + ".switch_matrix.tsv", matrix_rows, ["expected_identity"] + col_ids)
    frac_rows = []
    total = len(barcodes) or 1
    for ident, cnt in Counter(str(au[bc]["assignment"]) for bc in barcodes if bc in au).most_common():
        frac_rows.append({"identity": ident, "count": cnt, "fraction": cnt / total, "source": "audit_unconstrained"})
    write_tsv(str(out_prefix) + ".identity_fractions.tsv", frac_rows, ["identity", "count", "fraction", "source"])

    print(f"Wrote {out_prefix}.swap_report.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
