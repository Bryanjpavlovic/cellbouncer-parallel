#!/usr/bin/env python3
"""Prepare post-hoc CellBouncer swap-audit inputs.

This script validates per-library core demux outputs, canonicalizes the broad and
constrained identity universes, creates --reuse_counts-compatible audit prefixes,
and writes one capability manifest per library.

Revision note:
- demux_parallel audit thread count -t 80 -> -t 8 (both the unconstrained -i and
  constrained -I commands). The audit runs --reuse_counts, which skips the BAM
  read pass (the only thread-hungry section), so 80 threads were never used.
  Must stay consistent with the POSTHOC --cpus-per-task in the orchestrator
  (also lowered 80 -> 8).
"""
from __future__ import annotations

import argparse
import json
import os
import shlex
import sys
from pathlib import Path
from typing import Dict, List

from swap_audit_lib import (
    build_allowed_from_expected,
    canonical_identity,
    discover_demux_prefixes,
    identity_type,
    now_iso,
    normalize_library,
    parse_expected_metadata,
    parse_library_list,
    read_allowed_identities,
    read_samples,
    safe_symlink,
    sort_identity_key,
    write_tsv,
)

CORE_SUFFIXES = {
    "counts": ".counts",
    "samples": ".samples",
    "condf": ".condf",
    "assignments": ".assignments",
    "diagnostics": ".diagnostics.gz",
    "runner_ups": ".runner_ups.gz",
}

OPTIONAL_SUFFIXES = {
    "call_qc": ".call_qc.tsv.gz",
    "species_counts": ".species_counts",
    "species_condf": ".species_condf",
    "species_samples": ".species_samples",
    "species_assignments": ".species_assignments",
    "species_qc": ".species_qc.tsv",
    "atac_counts": ".atac.counts",
    "atac_call_qc": ".atac.call_qc.tsv.gz",
    "refined_assignments": ".refined_assignments",
}

REPORT_FIELDS = [
    "library", "preflight_status", "audit_verdict", "audit_flags", "feature_mode",
    "has_call_qc", "has_species_qc", "has_atac_qc", "has_calibrated_thresholds",
    "missing_optional_features", "warnings", "thresholds_source", "n_cells",
    "expected_identity", "audit_best_identity_unconstrained", "audit_best_identity_constrained",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--vcf_samples", required=True, help="VCF/demux sample list; one sample ID per line")
    p.add_argument("--expected_pool_metadata", required=True, help="TSV with library_id/lib and cell_identities columns")
    p.add_argument("--allowed_identities", default=None, help="Optional explicit constrained identity universe; one identity per line")
    p.add_argument("--allowed_identities_source", choices=["auto", "file", "expected_pool"], default="auto")
    p.add_argument("--demux_root", required=True, help="Root to recursively search for lib<N>_demuxed.counts prefixes")
    p.add_argument("--audit_root", required=True, help="Output root for audit files")
    p.add_argument("--libraries", default="all", help="Comma list like lib1,lib2 or range like 1-40 or all")
    p.add_argument("--mode", choices=["best-effort", "strict"], default="best-effort")
    p.add_argument("--thresholds", default=None, help="Optional recommended_thresholds.tsv")
    p.add_argument("--panel_metadata", default=None, help="Optional panel_metadata.tsv; used when generating tetra_score_calls commands and manifest provenance")
    p.add_argument("--call_qc_root", default=None, help="Optional root for existing lib<N>.call_qc.tsv.gz files")
    p.add_argument("--refined_assignments_root", default=None, help="Optional root for refined assignment files")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--skip-existing", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def detect_optional(prefix: str, lib: str, audit_dir: Path, args: argparse.Namespace) -> Dict[str, str | None]:
    """Detect optional files without confusing raw species inputs with completed species QC.

    call_qc/species_qc are written under audit_dir by tetra_score_calls.  We store
    their planned paths even before they exist so a manifest generated before
    scoring can still be used after the scoring command finishes.  Runtime
    consumers check existence before claiming a capability is available.
    """
    opt: Dict[str, str | None] = {}
    planned_call_qc = audit_dir / f"{lib}.call_qc.tsv.gz"
    planned_species_qc = audit_dir / f"{lib}.species_qc.tsv"
    for key, suf in OPTIONAL_SUFFIXES.items():
        candidates: List[Path] = []
        if key == "call_qc":
            candidates.extend([planned_call_qc, Path(prefix + suf)])
            if args.call_qc_root:
                candidates.extend([Path(args.call_qc_root) / lib / f"{lib}.call_qc.tsv.gz", Path(args.call_qc_root) / f"{lib}.call_qc.tsv.gz"])
            found = next((str(c.resolve()) for c in candidates if c.exists()), str(planned_call_qc.resolve()))
        elif key == "species_qc":
            candidates.extend([planned_species_qc, Path(prefix + suf)])
            found = next((str(c.resolve()) for c in candidates if c.exists()), str(planned_species_qc.resolve()))
        elif key == "atac_call_qc":
            candidates.extend([audit_dir / f"{lib}.atac.call_qc.tsv.gz", Path(prefix + suf)])
            found = next((str(c.resolve()) for c in candidates if c.exists()), None)
        elif key == "refined_assignments":
            candidates.append(Path(prefix + suf))
            if args.refined_assignments_root:
                candidates.extend([
                    Path(args.refined_assignments_root) / lib / f"{lib}.refined_assignments",
                    Path(args.refined_assignments_root) / f"{lib}.refined_assignments",
                    Path(args.refined_assignments_root) / f"{lib}_demuxed.refined_assignments",
                ])
            found = next((str(c.resolve()) for c in candidates if c.exists()), None)
        else:
            candidates.append(Path(prefix + suf))
            found = next((str(c.resolve()) for c in candidates if c.exists()), None)
        opt[key] = found
    opt["thresholds"] = str(Path(args.thresholds).resolve()) if args.thresholds and Path(args.thresholds).exists() else None
    opt["panel_metadata"] = str(Path(args.panel_metadata).resolve()) if args.panel_metadata and Path(args.panel_metadata).exists() else None
    return opt


def _path_exists(p: str | None) -> bool:
    return bool(p) and Path(str(p)).exists()


def feature_mode_and_warnings(opt: Dict[str, str | None]) -> tuple[str, List[str], List[str]]:
    """Return capabilities based on completed optional products, not raw inputs.

    Raw native species files are inputs for tetra_score_calls.  Species QC is only
    considered present when lib<N>.species_qc.tsv exists.
    """
    missing: List[str] = []
    warnings: List[str] = []
    has_call = _path_exists(opt.get("call_qc"))
    has_species_qc = _path_exists(opt.get("species_qc"))
    has_species_bundle = all(_path_exists(opt.get(k)) for k in ("species_counts", "species_condf", "species_samples"))
    has_atac_counts = _path_exists(opt.get("atac_counts"))
    has_atac_call_qc = _path_exists(opt.get("atac_call_qc"))
    has_atac = has_atac_counts and has_atac_call_qc

    if not has_call:
        missing.append("call_qc")
        warnings.append("CALL_QC_MISSING_REDUCED_MODE")
    if not has_species_qc:
        missing.append("species_qc")
        if has_species_bundle:
            warnings.append("SPECIES_QC_NOT_YET_COMPUTED")
        else:
            for k in ("species_counts", "species_condf", "species_samples"):
                if not _path_exists(opt.get(k)):
                    missing.append(k)
            warnings.append("SPECIES_NATIVE_BUNDLE_MISSING")
    if not has_atac:
        if has_atac_counts and not has_atac_call_qc:
            missing.append("atac_call_qc")
            warnings.append("ATAC_QC_UNAVAILABLE")
        else:
            missing.extend([x for x in ("atac_counts", "atac_call_qc") if not _path_exists(opt.get(x))])
            warnings.append("ATAC_COUNTS_MISSING")
    if not _path_exists(opt.get("thresholds")):
        missing.append("thresholds")
        warnings.append("USING_DEFAULT_THRESHOLDS")

    if has_call and has_species_qc and has_atac:
        mode = "CORE_PLUS_CALL_QC_AND_SPECIES_AND_ATAC"
    elif has_call and has_species_qc:
        mode = "CORE_PLUS_CALL_QC_AND_SPECIES"
    elif has_call:
        mode = "CORE_PLUS_CALL_QC"
    elif has_species_qc:
        mode = "CORE_PLUS_SPECIES"
    else:
        mode = "CORE_ONLY"
    return mode, sorted(set(missing)), sorted(set(warnings))


def write_failed_report(audit_dir: Path, lib: str, reason: str, warnings: List[str]) -> None:
    row = {
        "library": lib,
        "preflight_status": reason,
        "audit_verdict": "NOT_RUN",
        "audit_flags": reason,
        "feature_mode": "NOT_RUN",
        "has_call_qc": False,
        "has_species_qc": False,
        "has_atac_qc": False,
        "has_calibrated_thresholds": False,
        "missing_optional_features": "",
        "warnings": warnings,
        "thresholds_source": "defaults_uncalibrated",
        "n_cells": "NA",
        "expected_identity": "",
        "audit_best_identity_unconstrained": "NA",
        "audit_best_identity_constrained": "NA",
    }
    write_tsv(audit_dir / f"{lib}.swap_report.tsv", [row], REPORT_FIELDS)


def write_command_script(path: Path, command_lines: List[str]) -> None:
    """Write a strict bash command script.

    Per-library command scripts avoid races when multiple POSTHOC jobs call this
    preparer concurrently with the same audit_root.  The root-level scripts are
    still written for manual/batch use, but orchestrated single-library jobs
    should use the per-library scripts written inside audit_root/libN/.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "#!/usr/bin/env bash\nset -Eeuo pipefail\n\n" + "\n".join(command_lines) + "\n"
    path.write_text(body)
    os.chmod(path, 0o755)


def main() -> int:
    args = parse_args()
    audit_root = Path(args.audit_root).resolve()
    audit_root.mkdir(parents=True, exist_ok=True)

    samples = read_samples(args.vcf_samples)
    if not samples:
        raise SystemExit(f"ERROR: sample list is empty: {args.vcf_samples}")
    sample_set = set(samples)
    expected = parse_expected_metadata(args.expected_pool_metadata, samples)

    source = args.allowed_identities_source
    if source == "auto":
        source = "file" if args.allowed_identities else "expected_pool"
    if source == "file":
        if not args.allowed_identities:
            raise SystemExit("ERROR: --allowed_identities_source file requires --allowed_identities")
        allowed, allowed_meta = read_allowed_identities(args.allowed_identities, samples)
        allowed_meta["allowed_identities_source"] = "file"
        allowed_meta["allowed_identities_input"] = str(Path(args.allowed_identities).resolve())
    else:
        allowed, allowed_meta = build_allowed_from_expected(expected, samples, include_all_singlets=True)
        allowed_meta["allowed_identities_source"] = "expected_pool"
        allowed_meta["allowed_identities_input"] = "NA"

    # Broad unconstrained singlets file is generated directly from samples.
    all_singlets_file = audit_root / f"all_{len(samples)}_individuals.txt"
    all_singlets_file.write_text("\n".join(samples) + "\n")

    constrained_file = audit_root / "global_allowed_singlets_homotypics_fusions.txt"
    constrained_file.write_text("\n".join(allowed) + "\n")
    (audit_root / "allowed_identities_manifest.json").write_text(json.dumps(allowed_meta, indent=2, sort_keys=True) + "\n")

    discovered = discover_demux_prefixes(args.demux_root)
    libs = parse_library_list(args.libraries, expected)
    if not libs:
        libs = sorted(discovered.keys(), key=lambda x: int(x[3:]) if x.startswith("lib") and x[3:].isdigit() else 10**9)

    prep_rows = []
    command_lines: List[str] = []
    score_command_lines: List[str] = []
    for lib in libs:
        lib = normalize_library(lib)
        audit_dir = audit_root / lib
        audit_dir.mkdir(parents=True, exist_ok=True)
        prefix = discovered.get(lib)
        if not prefix:
            reason = "FAIL_MISSING_DEMUX_PREFIX"
            write_failed_report(audit_dir, lib, reason, [reason])
            prep_rows.append({"library": lib, "preflight_status": reason, "demux_prefix": ""})
            continue

        print(f"Library {lib}: selected demux prefix {Path(prefix).resolve()}")

        core_files = {k: str(Path(prefix + suf).resolve()) for k, suf in CORE_SUFFIXES.items()}
        missing_core = [k for k, p in core_files.items() if not Path(p).exists()]
        expected_ids = expected.get(lib, [])
        if not expected_ids:
            missing_core.append("expected_pool_metadata")
        unknown_expected = []
        for ident in expected_ids:
            for comp in ident.replace("+", " ").split():
                if comp not in sample_set:
                    unknown_expected.append(comp)
        if unknown_expected:
            missing_core.append("invalid_expected_identity_component:" + ",".join(sorted(set(unknown_expected))))
        if missing_core:
            reason = "FAIL_" + ";".join(missing_core)
            write_failed_report(audit_dir, lib, reason, ["INVALID_EXPECTED_LABEL_LIBRARY_SKIPPED"])
            prep_rows.append({"library": lib, "preflight_status": reason, "demux_prefix": prefix})
            continue

        for audit_name in ("audit_unconstrained", "audit_constrained"):
            audit_prefix = audit_dir / f"{lib}_{audit_name}"
            for key in ("counts", "samples", "condf"):
                safe_symlink(core_files[key], str(audit_prefix) + CORE_SUFFIXES[key], overwrite=args.overwrite)

        opt_files = detect_optional(prefix, lib, audit_dir, args)
        feature_mode, missing_optional, warnings = feature_mode_and_warnings(opt_files)
        if args.mode == "strict" and missing_optional:
            reason = "FAIL_STRICT_MODE_OPTIONAL_MISSING"
            write_failed_report(audit_dir, lib, reason, warnings)
            prep_rows.append({"library": lib, "preflight_status": reason, "demux_prefix": prefix})
            continue

        manifest = {
            "schema_version": "swap_audit_v1",
            "created_at": now_iso(),
            "library": lib,
            "demux_prefix": str(Path(prefix).resolve()),
            "core_files": core_files,
            "optional_files": opt_files,
            "feature_mode": feature_mode,
            "missing_optional_features": missing_optional,
            "warnings": warnings,
            "allowed_identities": allowed_meta,
            "global_files": {"panel_metadata": opt_files.get("panel_metadata")},
            "audit_files": {
                "all_singlets": str(all_singlets_file),
                "constrained_identities": str(constrained_file),
                "unconstrained_prefix": str(audit_dir / f"{lib}_audit_unconstrained"),
                "constrained_prefix": str(audit_dir / f"{lib}_audit_constrained"),
            },
            "expected_identities": expected_ids,
        }
        (audit_dir / f"{lib}.capabilities.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

        lib_audit_commands = [
            "demux_parallel "
            f"-o {audit_dir / (lib + '_audit_unconstrained')} --reuse_counts "
            f"-i {all_singlets_file} -D 0.5 -N 0 --n_runner_ups 250 --close_threshold 20 -t 8",
            "demux_parallel "
            f"-o {audit_dir / (lib + '_audit_constrained')} --reuse_counts "
            f"-I {constrained_file} -D 0.5 -N 0 --n_runner_ups 250 --close_threshold 20 -t 8",
        ]
        command_lines.extend(lib_audit_commands)

        q = shlex.quote
        score_cmd = [
            "tetra_score_calls",
            f"--counts {q(core_files['counts'])}",
            f"--samples {q(core_files['samples'])}",
            f"--assignments {q(core_files['assignments'])}",
            f"--diagnostics {q(core_files['diagnostics'])}",
            f"--runner_ups {q(core_files['runner_ups'])}",
            f"--output {q(str(opt_files['call_qc']))}",
            f"--libname {q(lib)}",
        ]
        if opt_files.get("panel_metadata"):
            score_cmd.append(f"--panel_metadata {q(str(opt_files['panel_metadata']))}")
        if all(opt_files.get(k) and Path(str(opt_files[k])).exists() for k in ("species_counts", "species_condf", "species_samples")):
            score_cmd.extend([
                f"--species_counts {q(str(opt_files['species_counts']))}",
                f"--species_condf {q(str(opt_files['species_condf']))}",
                f"--species_samples {q(str(opt_files['species_samples']))}",
            ])
        lib_score_command = " ".join(score_cmd)
        score_command_lines.append(lib_score_command)

        write_command_script(audit_dir / f"{lib}.run_tetra_score_calls_commands.sh", [lib_score_command])
        write_command_script(audit_dir / f"{lib}.run_audit_demux_parallel_commands.sh", lib_audit_commands)

        prep_rows.append({
            "library": lib,
            "preflight_status": "PASS",
            "demux_prefix": str(Path(prefix).resolve()),
            "feature_mode": feature_mode,
            "missing_optional_features": ",".join(missing_optional),
            "warnings": ",".join(warnings),
        })

    write_tsv(audit_root / "prepare_summary.tsv", prep_rows, ["library", "preflight_status", "demux_prefix", "feature_mode", "missing_optional_features", "warnings"])
    write_command_script(audit_root / "run_audit_demux_parallel_commands.sh", command_lines)
    write_command_script(audit_root / "run_tetra_score_calls_commands.sh", score_command_lines)

    print(f"Wrote {audit_root}")
    print(f"Prepared {sum(1 for r in prep_rows if r.get('preflight_status') == 'PASS')} libraries; {sum(1 for r in prep_rows if r.get('preflight_status') != 'PASS')} failed preflight")
    print(f"Allowed identities: {constrained_file} ({allowed_meta['n_unique_after_canonicalization']} unique)")
    print(f"Audit demux commands: {audit_root / 'run_audit_demux_parallel_commands.sh'}")
    print(f"Per-cell QC commands: {audit_root / 'run_tetra_score_calls_commands.sh'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
