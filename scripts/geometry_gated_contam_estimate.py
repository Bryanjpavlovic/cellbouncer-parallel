#!/usr/bin/env python3
"""Run a frozen CK geometry-gated source-exclusion condition on real libraries.

This production helper is invoked by orchestrate_pipeline.py. It runs the
unchanged tet_contam_estimate binary at a configured base source-exclusion
strength, evaluates the frozen geometry gate from estimator-emitted base
endpoint evidence, lazily runs the configured fallback strength when required
by the gate or by a base per-cell fit failure, and atomically publishes one
selected per-cell result plus an explicit gate audit.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping, Sequence

WORKER_VERSION = "geometry_gated_contam_estimate_V2"
CONTRACT_VERSION = "geometry_gated_contam_estimate_run_contract_V2"

_ENDPOINT_REQUIRED_SUFFIXES = (
    ".contam_rate",
    ".contam_diagnostics.tsv",
    ".allele_ratio",
    ".contam_prof",
    ".profile_fit_diagnostics.tsv",
    ".condf_coverage.tsv",
    ".run_contract.json",
)
_FINAL_SELECTED_SUFFIXES = (
    ".contam_rate",
    ".contam_diagnostics.tsv",
    ".geometry_gate_audit.tsv",
    ".allele_ratio",
    ".contam_prof",
    ".profile_fit_diagnostics.tsv",
    ".condf_coverage.tsv",
)
_OPTIONAL_BASE_SUFFIXES = (
    ".model_fit.tsv",
    ".class_residuals.tsv",
    ".class_residuals.tsv.gz",
    ".decontam.assignments",
    ".pass1.contam_rate",
    ".pass1.contam_prof",
)
_GATE_REQUIRED_FIELDS = {
    "barcode",
    "is_heterotypic",
    "ambient_parent_axis_alpha_scoring",
    "ambient_orthogonal_norm_scoring",
    "excluded_parent_A_mass_raw",
    "excluded_parent_B_mass_raw",
    "parent_axis_geometry_status",
}
_V1_GATE_FIELDS = [
    "geometry_gate_version",
    "geometry_gate_triggered",
    "geometry_gate_parent_axis_alpha",
    "geometry_gate_ambient_orthogonal_norm",
    "geometry_gate_max_raw_excluded_parent_mass",
    "geometry_gate_parent_axis_alpha_threshold",
    "geometry_gate_ambient_orthogonal_norm_threshold",
    "geometry_gate_parent_mass_threshold",
    "source_exclusion_strength_base",
    "source_exclusion_strength_fallback",
    "source_exclusion_strength_selected",
    "geometry_gate_selection_reason",
    "base_c_selected",
    "fallback_c_selected",
    "selected_c",
]
_V2_FIELDS = [
    "geometry_gate_evaluable",
    "base_rate_status",
    "fallback_rate_status",
    "selected_endpoint",
    "selected_rate_status",
    "base_optimizer_status",
    "fallback_optimizer_status",
    "base_profile_validation_status",
    "fallback_profile_validation_status",
    "base_fit_failure_reason",
    "fallback_fit_failure_reason",
    "selected_c_se",
]
_AUDIT_FIELDS = [*_V1_GATE_FIELDS, *_V2_FIELDS]


class GeometryGateError(RuntimeError):
    """Raised when endpoint evidence cannot support deterministic selection."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _fsync_directory(path: Path) -> None:
    try:
        fd = os.open(str(path), os.O_RDONLY)
    except OSError:
        return
    try:
        try:
            os.fsync(fd)
        except OSError:
            pass
    finally:
        os.close(fd)


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with tmp.open("w", encoding="utf-8", newline="") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    finally:
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass


def _atomic_copy(source: Path, destination: Path) -> None:
    if not source.is_file() or source.stat().st_size <= 0:
        raise GeometryGateError(f"required endpoint output is missing: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    try:
        with source.open("rb") as src, tmp.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)
            dst.flush()
            os.fsync(dst.fileno())
        os.replace(tmp, destination)
    finally:
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass


def _option_occurrences(arguments: Sequence[str], names: Iterable[str]) -> list[tuple[int, str, str]]:
    ordered = tuple(names)
    accepted = set(ordered)
    found: list[tuple[int, str, str]] = []
    i = 0
    while i < len(arguments):
        token = str(arguments[i])
        if token in accepted:
            if i + 1 >= len(arguments):
                raise GeometryGateError(f"option {token} is missing its value")
            found.append((i, token, str(arguments[i + 1])))
            i += 2
            continue
        matched = next((name for name in ordered if token.startswith(name + "=")), None)
        if matched is not None:
            found.append((i, matched, token[len(matched) + 1 :]))
        i += 1
    return found


def _option_value(arguments: Sequence[str], names: Iterable[str]) -> str:
    found = _option_occurrences(arguments, names)
    if not found:
        raise GeometryGateError(f"required option is absent: {sorted(set(names))}")
    values = {value for _, _, value in found}
    if len(values) != 1:
        raise GeometryGateError(
            f"conflicting values for {sorted(set(names))}: {sorted(values)}"
        )
    return found[0][2]


def _remove_option(arguments: Sequence[str], names: Iterable[str]) -> list[str]:
    ordered = tuple(names)
    accepted = set(ordered)
    output: list[str] = []
    i = 0
    while i < len(arguments):
        token = str(arguments[i])
        if token in accepted:
            if i + 1 >= len(arguments):
                raise GeometryGateError(f"option {token} is missing its value")
            i += 2
            continue
        if any(token.startswith(name + "=") for name in ordered):
            i += 1
            continue
        output.append(token)
        i += 1
    return output


def _remove_flag(arguments: Sequence[str], names: Iterable[str]) -> list[str]:
    """Remove valueless flags without consuming the following argument."""
    accepted = set(names)
    return [str(token) for token in arguments if str(token) not in accepted]


def _set_option_unique(arguments: Sequence[str], names: Iterable[str], value: str) -> list[str]:
    ordered = tuple(names)
    first_index: int | None = None
    preferred_name = ordered[0]
    found = _option_occurrences(arguments, ordered)
    if found:
        first_index = found[0][0]
        preferred_name = found[0][1]
    stripped = _remove_option(arguments, ordered)
    if first_index is None:
        return [*stripped, preferred_name, value]
    # Count tokens surviving before the first removed occurrence.
    prefix = _remove_option(arguments[:first_index], ordered)
    insert_at = len(prefix)
    return [*stripped[:insert_at], preferred_name, value, *stripped[insert_at:]]


def _assert_no_active_forwarded_inputs(arguments: Sequence[str]) -> None:
    prohibited = {
        "warm_start": ("--warm_start", "-W"),
        "fixed_ambient": ("--fixed_ambient", "-A"),
        "fix_r": ("--fix_r",),
    }
    for label, names in prohibited.items():
        if _option_occurrences(arguments, names):
            raise GeometryGateError(
                f"geometry-gated free-r fitted-profile endpoint cannot use --{label}"
            )


def _endpoint_arguments(
    estimator_arguments: Sequence[str], endpoint_prefix: Path, strength: float,
    condition_key: str,
) -> list[str]:
    _assert_no_active_forwarded_inputs(estimator_arguments)
    condition_values = [value for _, _, value in _option_occurrences(
        estimator_arguments, ("--condition_key",)
    )]
    if any(value != condition_key for value in condition_values):
        raise GeometryGateError(
            "forwarded estimator --condition_key conflicts with helper --condition-key: "
            f"{condition_values!r} != {condition_key!r}"
        )
    arguments = _set_option_unique(
        estimator_arguments, ("-o", "--output_prefix"), str(endpoint_prefix)
    )
    arguments = _set_option_unique(
        arguments, ("--run_contract",), str(endpoint_prefix) + ".run_contract.json"
    )
    arguments = _remove_option(arguments, ("--source_exclusion_strength",))
    arguments = _remove_flag(arguments, ("--loo",))
    arguments = _set_option_unique(arguments, ("--condition_key",), condition_key)
    arguments.extend(["--source_exclusion_strength", format(strength, ".17g")])
    return arguments


def _link_endpoint_inputs(outer_prefix: Path, endpoint_prefix: Path) -> None:
    for suffix in (".counts", ".condf", ".samples", ".assignments"):
        source = Path(str(outer_prefix) + suffix)
        if not source.is_file() or source.stat().st_size <= 0:
            raise GeometryGateError(f"required staged estimator input is missing: {source}")
        destination = Path(str(endpoint_prefix) + suffix)
        try:
            destination.unlink()
        except FileNotFoundError:
            pass
        try:
            os.link(source, destination)
        except OSError:
            shutil.copy2(source, destination)


def _cleanup_prefix(prefix: Path) -> None:
    for path in prefix.parent.glob(prefix.name + ".*"):
        try:
            if path.is_file() or path.is_symlink():
                path.unlink()
        except FileNotFoundError:
            pass


def _cleanup_staging_for_outer(outer_prefix: Path) -> None:
    pattern = f".{outer_prefix.name}.geometry_publish.*"
    for path in outer_prefix.parent.glob(pattern):
        try:
            if path.is_file() or path.is_symlink():
                path.unlink()
        except FileNotFoundError:
            pass


def _cleanup_final_outputs(outer_prefix: Path) -> None:
    for suffix in (*_FINAL_SELECTED_SUFFIXES, *_OPTIONAL_BASE_SUFFIXES, ".run_contract.json"):
        try:
            Path(str(outer_prefix) + suffix).unlink()
        except FileNotFoundError:
            pass


def _validate_endpoint_required_files(prefix: Path) -> None:
    for suffix in _ENDPOINT_REQUIRED_SUFFIXES:
        path = Path(str(prefix) + suffix)
        if not path.is_file():
            raise GeometryGateError(f"required endpoint file is missing: {path}")
        if suffix != ".contam_rate" and path.stat().st_size <= 0:
            raise GeometryGateError(f"required endpoint file is empty: {path}")


def _run_endpoint(
    estimator_binary: Path,
    estimator_arguments: Sequence[str],
    outer_prefix: Path,
    endpoint_prefix: Path,
    strength: float,
    condition_key: str,
) -> None:
    _link_endpoint_inputs(outer_prefix, endpoint_prefix)
    command = [
        str(estimator_binary),
        *_endpoint_arguments(
            estimator_arguments, endpoint_prefix, strength, condition_key
        ),
    ]
    print(
        f"=== geometry-gated endpoint lambda={strength:g}: {endpoint_prefix} ===",
        file=sys.stderr,
        flush=True,
    )
    print(" ".join(command), file=sys.stderr, flush=True)
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        raise GeometryGateError(
            f"tet_contam_estimate endpoint lambda={strength:g} exited "
            f"with status {completed.returncode}"
        )


def _read_tsv_rows(
    path: Path, *, allow_empty: bool = False,
) -> tuple[list[str], list[dict[str, str]], dict[str, dict[str, str]]]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise GeometryGateError(f"endpoint TSV is missing: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        if "barcode" not in fieldnames:
            raise GeometryGateError(f"endpoint TSV lacks barcode column: {path}")
        rows: list[dict[str, str]] = []
        by_barcode: dict[str, dict[str, str]] = {}
        for row in reader:
            normalized = {
                str(key): "" if value is None else str(value)
                for key, value in row.items()
            }
            barcode = normalized.get("barcode", "")
            if not barcode:
                raise GeometryGateError(f"endpoint TSV contains an empty barcode: {path}")
            if barcode in by_barcode:
                raise GeometryGateError(f"duplicate barcode {barcode!r} in {path}")
            rows.append(normalized)
            by_barcode[barcode] = normalized
    if not rows and not allow_empty:
        raise GeometryGateError(f"endpoint TSV has no data rows: {path}")
    return fieldnames, rows, by_barcode


def _diagnostic_value_is_unavailable(raw: object) -> bool:
    value = str(raw).strip().lower()
    return value in {"", "nan", "unavailable", "not_available", "na", "n/a"}


def _diagnostic_fit_failure_reason(row: Mapping[str, str]) -> str:
    prior_reason = str(row.get("prior_training_reason", "")).strip()
    profile_status = str(row.get("profile_validation_status", "")).strip()
    optimizer_status = str(row.get("optimizer_status", "")).strip()
    compile_status = str(row.get("compile_status", "")).strip()
    fallback_status = str(row.get("diagnostic_fallback_status", "")).strip()

    if prior_reason != "fit_failed":
        return ""
    if not _diagnostic_value_is_unavailable(row.get("contam_mle", "")):
        return ""

    # Current estimator writes profile_validation_status from the production MLE
    # optimizer status. Compile failures and free-fit failures are explicit here.
    explicit_statuses = {
        "two_dimensional_boundary_profile_support_too_weak",
        "two_dimensional_boundary_profile_kkt_not_supported",
        "two_dimensional_boundary_profile_too_broad",
        "two_dimensional_boundary_invalid_continuous_profile_interval_inputs",
        "two_dimensional_boundary_nonfinite_profile_likelihood_at_optimum",
        "two_dimensional_boundary_nonfinite_continuous_profile_interval",
        "two_dimensional_boundary_invalid_continuous_profile_interval_invariants",
        "regularized_two_dimensional_boundary_not_selectable",
        "invalid_two_dimensional_candidate",
        "invalid_two_dimensional_profile_inputs",
        "not_a_two_dimensional_c_boundary_candidate",
        "nonfinite_two_dimensional_profile_likelihood",
        "nonfinite_two_dimensional_boundary_likelihood",
        "unidentifiable_flat_profile_c_likelihood",
        "unidentifiable_nearly_flat_profile_c_likelihood",
        "nonfinite_two_dimensional_profile_refinement",
        "invalid_continuous_profile_interval_inputs",
        "nonfinite_profile_likelihood_at_optimum",
        "nonfinite_continuous_profile_interval",
        "invalid_continuous_profile_interval_invariants",
        "two_dimensional_profile_refinement_not_locally_maximal",
        "nonfinite_two_dimensional_projected_gradient",
        "no_observations",
        "regularized_prior_not_log_concave",
        "fixed_r_brent_failed",
        "brent_failed",
        "no_weighted_univariate_observations",
        "regularized_boundary_not_selectable",
        "not_a_boundary_candidate",
        "unidentifiable_flat_c_likelihood",
        "unidentifiable_weak_c_likelihood",
        "nonfinite_boundary_likelihood",
        "boundary_likelihood_support_too_weak",
        "nonfinite_boundary_derivative",
        "boundary_kkt_direction_not_supported",
    }
    if profile_status in explicit_statuses:
        return profile_status
    for status in explicit_statuses:
        if status and status in optimizer_status:
            return status
    if compile_status and compile_status not in {"ok", "native_assignment_rows"}:
        return compile_status
    if fallback_status in {"diagnostic_fixed_r_solver_failed"}:
        return fallback_status
    return ""


def _diagnostic_has_explicit_fit_failure(row: Mapping[str, str]) -> bool:
    return bool(_diagnostic_fit_failure_reason(row))


def _parse_rate_number(raw: str, *, context: str) -> float:
    try:
        value = float(raw)
    except ValueError as exc:
        raise GeometryGateError(f"nonnumeric {context}: {raw!r}") from exc
    if math.isinf(value):
        raise GeometryGateError(f"infinite {context}: {raw!r}")
    return value


def _read_rate_rows(
    path: Path,
    diagnostics_by_barcode: Mapping[str, Mapping[str, str]],
    diagnostic_order: Sequence[str],
) -> dict[str, dict[str, object]]:
    if not path.is_file():
        raise GeometryGateError(f"endpoint contamination-rate file is missing: {path}")
    raw_values: dict[str, tuple[str, str]] = {}
    if path.stat().st_size > 0:
        with path.open("r", encoding="utf-8", newline="") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\r\n")
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) != 3:
                    raise GeometryGateError(
                        f"malformed contamination-rate row {line_number} in {path}: expected 3 fields"
                    )
                barcode, c_raw, se_raw = fields
                if not barcode:
                    raise GeometryGateError(f"empty barcode in contamination-rate row {line_number}: {path}")
                if barcode in raw_values:
                    raise GeometryGateError(f"duplicate barcode {barcode!r} in {path}")
                if barcode not in diagnostics_by_barcode:
                    raise GeometryGateError(
                        f"rate barcode {barcode!r} is absent from matching diagnostics: {path}"
                    )
                raw_values[barcode] = (c_raw, se_raw)

    states: dict[str, dict[str, object]] = {}
    for barcode in diagnostic_order:
        diag = diagnostics_by_barcode[barcode]
        pair = raw_values.get(barcode)
        if pair is None:
            reason = _diagnostic_fit_failure_reason(diag)
            if not reason:
                raise GeometryGateError(
                    f"barcode {barcode!r} is absent from {path} without explicit diagnostic fit failure"
                )
            states[barcode] = {
                "status": "failed", "c": math.nan, "c_se": math.nan,
                "c_text": "nan", "c_se_text": "nan", "failure_reason": reason,
            }
            continue

        c_raw, se_raw = pair
        if c_raw == "nan":
            if se_raw != "nan":
                raise GeometryGateError(
                    f"failed rate must be canonical nan/nan for barcode {barcode!r}: {c_raw!r}, {se_raw!r}"
                )
            reason = _diagnostic_fit_failure_reason(diag)
            if not reason:
                raise GeometryGateError(
                    f"explicit nan rate for barcode {barcode!r} lacks diagnostic fit-failure evidence"
                )
            states[barcode] = {
                "status": "failed", "c": math.nan, "c_se": math.nan,
                "c_text": "nan", "c_se_text": "nan", "failure_reason": reason,
            }
            continue

        c = _parse_rate_number(c_raw, context=f"contamination estimate for {barcode!r} in {path}")
        se = _parse_rate_number(se_raw, context=f"contamination SE for {barcode!r} in {path}")
        if math.isnan(c):
            raise GeometryGateError(
                f"noncanonical NaN contamination estimate for barcode {barcode!r}: {c_raw!r}"
            )
        if not (0.0 <= c <= 1.0):
            raise GeometryGateError(
                f"contamination estimate outside [0,1] for barcode {barcode!r}: {c_raw!r}"
            )
        if math.isnan(se):
            # Only canonical lower-case nan is accepted for uncertainty.
            if se_raw != "nan":
                raise GeometryGateError(
                    f"noncanonical NaN contamination SE for barcode {barcode!r}: {se_raw!r}"
                )
        elif se < 0.0:
            raise GeometryGateError(
                f"negative contamination SE for barcode {barcode!r}: {se_raw!r}"
            )
        if _diagnostic_has_explicit_fit_failure(diag):
            raise GeometryGateError(
                f"finite rate for barcode {barcode!r} conflicts with diagnostics declaring selected fit failure"
            )
        states[barcode] = {
            "status": "finite", "c": c, "c_se": se,
            "c_text": c_raw, "c_se_text": se_raw, "failure_reason": "",
        }
    return states


def _strict_bool(raw: object, field: str) -> bool:
    if isinstance(raw, bool):
        return raw
    value = str(raw).strip().lower()
    if value in {"1", "1.0", "true"}:
        return True
    if value in {"0", "0.0", "false"}:
        return False
    raise GeometryGateError(f"gate diagnostic {field} is not boolean: {raw!r}")


def _float_or_nan(raw: object, field: str) -> float:
    text = str(raw).strip()
    try:
        value = float(text)
    except ValueError as exc:
        raise GeometryGateError(f"gate diagnostic {field} is not numeric: {text!r}") from exc
    if math.isinf(value):
        raise GeometryGateError(f"gate diagnostic {field} is infinite: {text!r}")
    return value


def _fmt_numeric(value: float) -> str:
    return "nan" if not math.isfinite(value) else format(value, ".17g")


def _selection_reason_v1_false(
    *, finite_geometry: bool, is_heterotypic: bool, alpha: float,
    orthogonal: float, parent_mass: float, alpha_threshold: float,
    orthogonal_threshold: float, parent_mass_threshold: float,
) -> str:
    if not finite_geometry:
        if not is_heterotypic:
            return "geometry_not_applicable_nonheterotypic_select_base"
        return "geometry_unavailable_select_base"
    failures: list[str] = []
    if alpha < alpha_threshold:
        failures.append("parent_axis_alpha_below_threshold")
    if orthogonal > orthogonal_threshold:
        failures.append("ambient_orthogonal_norm_above_threshold")
    if parent_mass < parent_mass_threshold:
        failures.append("raw_excluded_parent_mass_below_threshold")
    return "geometry_gate_false_select_base:" + ";".join(
        failures or ["no_threshold_failure"]
    )


def _evaluate_base_gate(
    base_fields: Sequence[str],
    base_rows: Sequence[Mapping[str, str]],
    base_rates: Mapping[str, Mapping[str, object]],
    *, alpha_threshold: float, orthogonal_threshold: float,
    parent_mass_threshold: float,
) -> dict[str, dict[str, object]]:
    missing = sorted(_GATE_REQUIRED_FIELDS - set(base_fields))
    if missing:
        raise GeometryGateError(
            "base endpoint diagnostics lack gate fields: " + ",".join(missing)
        )
    decisions: dict[str, dict[str, object]] = {}
    for row in base_rows:
        barcode = str(row["barcode"])
        state = base_rates[barcode]
        is_heterotypic = _strict_bool(row.get("is_heterotypic", ""), "is_heterotypic")
        alpha = _float_or_nan(row.get("ambient_parent_axis_alpha_scoring", "nan"), "ambient_parent_axis_alpha_scoring")
        orthogonal = _float_or_nan(row.get("ambient_orthogonal_norm_scoring", "nan"), "ambient_orthogonal_norm_scoring")
        raw_a = _float_or_nan(row.get("excluded_parent_A_mass_raw", "nan"), "excluded_parent_A_mass_raw")
        raw_b = _float_or_nan(row.get("excluded_parent_B_mass_raw", "nan"), "excluded_parent_B_mass_raw")
        parent_mass = max(raw_a, raw_b) if math.isfinite(raw_a) and math.isfinite(raw_b) else math.nan
        finite_geometry = all(math.isfinite(v) for v in (alpha, orthogonal, raw_a, raw_b))
        geometry_status = str(row.get("parent_axis_geometry_status", ""))
        evaluable = bool(finite_geometry and (not is_heterotypic or geometry_status == "ok"))

        if state["status"] == "finite" and is_heterotypic:
            if not finite_geometry or geometry_status != "ok":
                raise GeometryGateError(
                    f"heterotypic barcode {barcode!r} has malformed/missing base geometry for successful fit: "
                    f"status={geometry_status!r}, alpha={alpha!r}, orthogonal={orthogonal!r}, "
                    f"raw_parent_A={raw_a!r}, raw_parent_B={raw_b!r}"
                )
            evaluable = True

        # Preserve V1 behavior: gate formula is evaluated whenever all numeric
        # geometry is finite. Nonheterotypic rows commonly carry not-applicable
        # values and therefore remain base-selected.
        triggered = bool(
            finite_geometry
            and alpha >= alpha_threshold
            and orthogonal <= orthogonal_threshold
            and parent_mass >= parent_mass_threshold
        )
        if state["status"] == "failed" and not evaluable:
            triggered = False

        if triggered:
            v1_reason = "geometry_gate_triggered_select_fallback"
        else:
            v1_reason = _selection_reason_v1_false(
                finite_geometry=finite_geometry,
                is_heterotypic=is_heterotypic,
                alpha=alpha,
                orthogonal=orthogonal,
                parent_mass=parent_mass,
                alpha_threshold=alpha_threshold,
                orthogonal_threshold=orthogonal_threshold,
                parent_mass_threshold=parent_mass_threshold,
            )
        decisions[barcode] = {
            "triggered": triggered,
            "evaluable": evaluable,
            "v1_reason": v1_reason,
            "alpha": alpha,
            "orthogonal": orthogonal,
            "parent_mass": parent_mass,
            "is_heterotypic": is_heterotypic,
        }
    return decisions


def _read_profile_diagnostics(path: Path) -> tuple[list[str], dict[str, str]]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise GeometryGateError(f"profile-fit diagnostics missing: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = [
            {str(k): "" if v is None else str(v) for k, v in row.items()}
            for row in reader
        ]
    if not fields:
        raise GeometryGateError(f"profile-fit diagnostics lacks header: {path}")
    if len(rows) != 1:
        raise GeometryGateError(
            f"profile-fit diagnostics must contain exactly one data row: {path} has {len(rows)}"
        )
    return fields, rows[0]


def _path_value_empty(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, Mapping):
        return str(value.get("path", "")).strip() == ""
    return str(value).strip() == ""


def _basis_is_no_input(value: object) -> bool:
    return str(value if value is not None else "").strip().lower() in {
        "", "none", "unspecified"
    }


def _load_endpoint_contract(
    path: Path, expected_strength: float, condition_key: str,
) -> dict[str, object]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise GeometryGateError(f"endpoint run contract is missing: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise GeometryGateError(f"cannot read endpoint run contract {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise GeometryGateError(f"endpoint run contract is not an object: {path}")

    exact = {
        "production_contract_pass": True,
        "run_class": "production",
        "panel_mode": "interindividual",
        "condition_key": condition_key,
        "strict_condf": True,
        "assignments_basis": "library",
        "expected_lines_basis": "library",
        "ambient_candidates_basis": "library",
        "source_exclusion_explicit": True,
        "truth_assisted": False,
    }
    for field, expected in exact.items():
        if field not in value:
            raise GeometryGateError(f"endpoint run contract lacks {field}: {path}")
        if value[field] != expected:
            raise GeometryGateError(
                f"endpoint run contract {path} has {field}={value[field]!r}, expected {expected!r}"
            )
    try:
        observed_strength = float(value["source_exclusion_strength"])
    except (KeyError, TypeError, ValueError) as exc:
        raise GeometryGateError(
            f"endpoint run contract has invalid source-exclusion strength: {path}"
        ) from exc
    if not math.isfinite(observed_strength) or not math.isclose(
        observed_strength, expected_strength, rel_tol=0.0, abs_tol=1e-12
    ):
        raise GeometryGateError(
            f"endpoint run-contract strength mismatch: {observed_strength} != {expected_strength}: {path}"
        )

    for field in ("warm_start_basis", "fixed_ambient_basis", "fix_r_basis", "fixed_r_basis"):
        if field in value and not _basis_is_no_input(value[field]):
            raise GeometryGateError(f"endpoint contract has active {field}: {value[field]!r}")
    for field in ("warm_start", "fixed_ambient", "fix_r"):
        if field in value and not _path_value_empty(value[field]):
            raise GeometryGateError(f"endpoint contract has active {field} input: {value[field]!r}")
    if value.get("fixed_r_enabled", False) is not False:
        raise GeometryGateError(f"endpoint contract records fixed_r_enabled=true: {path}")
    if value.get("fixed_ambient_enabled", False) is not False:
        raise GeometryGateError(f"endpoint contract records fixed_ambient_enabled=true: {path}")
    return value


def _validate_contract_invariants(base: Mapping[str, object], fallback: Mapping[str, object]) -> None:
    fields = (
        "tool", "tool_version", "run_class", "panel_mode", "strict_condf",
        "assignments_basis", "expected_lines_basis", "ambient_candidates_basis",
        "condition_key", "synthetic_id", "source_exclusion_explicit",
        "production_contract_pass", "truth_assisted",
    )
    for field in fields:
        if base.get(field) != fallback.get(field):
            raise GeometryGateError(
                f"base/fallback run-contract invariant mismatch for {field}: "
                f"{base.get(field)!r} != {fallback.get(field)!r}"
            )


def _contract_summary(contract: Mapping[str, object]) -> dict[str, object]:
    keep = (
        "contract_version", "tool", "tool_version", "run_class", "panel_mode",
        "strict_condf", "assignments_basis", "expected_lines_basis",
        "ambient_candidates_basis", "warm_start_basis", "fixed_ambient_basis",
        "fix_r_basis", "condition_key", "synthetic_id", "source_exclusion_strength",
        "source_exclusion_explicit", "production_contract_pass",
        "production_contract_reason", "truth_assisted",
    )
    return {key: contract.get(key) for key in keep if key in contract}


def _read_allele_rows(
    path: Path, diagnostic_barcodes: set[str],
) -> tuple[list[str], dict[str, dict[str, str]]]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise GeometryGateError(f"endpoint allele-ratio file is missing or empty: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if "barcode" not in fields:
            raise GeometryGateError(f"endpoint allele-ratio file lacks barcode header: {path}")
        by_barcode: dict[str, dict[str, str]] = {}
        for row in reader:
            normalized = {str(k): "" if v is None else str(v) for k, v in row.items()}
            barcode = normalized.get("barcode", "")
            if not barcode:
                raise GeometryGateError(f"endpoint allele-ratio contains empty barcode: {path}")
            if barcode in by_barcode:
                raise GeometryGateError(f"duplicate barcode {barcode!r} in {path}")
            if barcode not in diagnostic_barcodes:
                raise GeometryGateError(
                    f"allele-ratio barcode {barcode!r} absent from matching diagnostics: {path}"
                )
            by_barcode[barcode] = normalized
    return fields, by_barcode


def _write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with tmp.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(fieldnames), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for row in rows:
                writer.writerow({name: row.get(name, "") for name in fieldnames})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    finally:
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass


def _build_selected_rows(
    base_fields: Sequence[str],
    base_rows: Sequence[Mapping[str, str]],
    base_rates: Mapping[str, Mapping[str, object]],
    fallback_by_barcode: Mapping[str, Mapping[str, str]],
    fallback_rates: Mapping[str, Mapping[str, object]],
    decisions: Mapping[str, Mapping[str, object]],
    *, fallback_evaluated: bool, gate_version: str, base_strength: float,
    fallback_strength: float, alpha_threshold: float, orthogonal_threshold: float,
    parent_mass_threshold: float,
) -> tuple[list[str], list[dict[str, str]], list[dict[str, str]], dict[str, str], Counter[str], Counter[str]]:
    final_fields = [f for f in base_fields if f not in _AUDIT_FIELDS] + _AUDIT_FIELDS
    selected_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []
    selected_endpoint: dict[str, str] = {}
    reason_counts: Counter[str] = Counter()
    summary_counts: Counter[str] = Counter()

    for base in base_rows:
        barcode = str(base["barcode"])
        base_state = base_rates[barcode]
        decision = decisions[barcode]
        fallback_state = fallback_rates.get(barcode)
        fallback_row = fallback_by_barcode.get(barcode)
        if fallback_evaluated and (fallback_state is None or fallback_row is None):
            raise GeometryGateError(f"fallback endpoint evidence incomplete for barcode {barcode!r}")

        if base_state["status"] == "failed":
            if not fallback_evaluated:
                raise GeometryGateError(f"base failed but fallback was not evaluated for barcode {barcode!r}")
            if fallback_state["status"] == "finite":
                endpoint = "fallback"
                selected_status = "finite"
                reason = "base_fit_failed_fallback_rescue"
                summary_counts["selected_fallback_rescue_cells"] += 1
            else:
                endpoint = "none"
                selected_status = "na"
                reason = "base_and_fallback_fit_failed_na"
                summary_counts["selected_na_cells"] += 1
        elif bool(decision["triggered"]):
            if not fallback_evaluated:
                raise GeometryGateError(f"geometry gate triggered but fallback was not evaluated for barcode {barcode!r}")
            if fallback_state["status"] == "finite":
                endpoint = "fallback"
                selected_status = "finite"
                reason = "geometry_gate_triggered_select_fallback"
                summary_counts["selected_fallback_geometry_cells"] += 1
            else:
                endpoint = "none"
                selected_status = "na"
                reason = "geometry_gate_triggered_fallback_fit_failed_na"
                summary_counts["selected_na_cells"] += 1
        else:
            endpoint = "base"
            selected_status = "finite"
            reason = str(decision["v1_reason"])
            summary_counts["selected_base_cells"] += 1

        if endpoint == "base":
            selected_state = base_state
            chosen = dict(base)
        elif endpoint == "fallback":
            selected_state = fallback_state
            chosen = dict(fallback_row)
        else:
            # If fallback was attempted and failed, its diagnostic row is the
            # most recent endpoint evidence while both failures are retained.
            selected_state = {"c_text": "nan", "c_se_text": "nan"}
            chosen = dict(fallback_row if fallback_evaluated else base)

        fallback_status = (
            str(fallback_state["status"]) if fallback_evaluated else "not_evaluated"
        )
        fallback_c = (
            str(fallback_state["c_text"])
            if fallback_evaluated and fallback_state["status"] == "finite"
            else "nan"
        )
        base_c = str(base_state["c_text"]) if base_state["status"] == "finite" else "nan"
        selected_c = str(selected_state["c_text"]) if selected_status == "finite" else "nan"
        selected_se = str(selected_state["c_se_text"]) if selected_status == "finite" else "nan"

        gate_values = {
            "geometry_gate_version": gate_version,
            "geometry_gate_triggered": "1" if decision["triggered"] else "0",
            "geometry_gate_parent_axis_alpha": _fmt_numeric(float(decision["alpha"])),
            "geometry_gate_ambient_orthogonal_norm": _fmt_numeric(float(decision["orthogonal"])),
            "geometry_gate_max_raw_excluded_parent_mass": _fmt_numeric(float(decision["parent_mass"])),
            "geometry_gate_parent_axis_alpha_threshold": format(alpha_threshold, ".17g"),
            "geometry_gate_ambient_orthogonal_norm_threshold": format(orthogonal_threshold, ".17g"),
            "geometry_gate_parent_mass_threshold": format(parent_mass_threshold, ".17g"),
            "source_exclusion_strength_base": format(base_strength, ".17g"),
            "source_exclusion_strength_fallback": format(fallback_strength, ".17g"),
            "source_exclusion_strength_selected": (
                format(base_strength, ".17g") if endpoint == "base"
                else format(fallback_strength, ".17g") if endpoint == "fallback"
                else "nan"
            ),
            "geometry_gate_selection_reason": reason,
            "base_c_selected": base_c,
            "fallback_c_selected": fallback_c,
            "selected_c": selected_c,
            "geometry_gate_evaluable": "1" if decision["evaluable"] else "0",
            "base_rate_status": str(base_state["status"]),
            "fallback_rate_status": fallback_status,
            "selected_endpoint": endpoint,
            "selected_rate_status": selected_status,
            "base_optimizer_status": str(base.get("optimizer_status", "")),
            "fallback_optimizer_status": str(fallback_row.get("optimizer_status", "")) if fallback_evaluated else "",
            "base_profile_validation_status": str(base.get("profile_validation_status", "")),
            "fallback_profile_validation_status": str(fallback_row.get("profile_validation_status", "")) if fallback_evaluated else "",
            "base_fit_failure_reason": str(base_state.get("failure_reason", "")),
            "fallback_fit_failure_reason": str(fallback_state.get("failure_reason", "")) if fallback_evaluated else "",
            "selected_c_se": selected_se,
        }
        for field in _AUDIT_FIELDS:
            chosen.pop(field, None)
        chosen.update(gate_values)
        selected_rows.append(chosen)
        audit_rows.append({"barcode": barcode, **gate_values})
        selected_endpoint[barcode] = endpoint
        reason_counts[reason] += 1

    return final_fields, selected_rows, audit_rows, selected_endpoint, reason_counts, summary_counts


def _write_selected_rate(
    path: Path, order: Sequence[str], selected_rows_by_barcode: Mapping[str, Mapping[str, str]]
) -> None:
    lines = [
        f"{barcode}\t{selected_rows_by_barcode[barcode]['selected_c']}\t{selected_rows_by_barcode[barcode]['selected_c_se']}\n"
        for barcode in order
    ]
    _atomic_write_text(path, "".join(lines))


def _write_selected_allele(
    path: Path,
    order: Sequence[str],
    selected_endpoint: Mapping[str, str],
    selected_rate_status: Mapping[str, str],
    base_fields: Sequence[str], base_rows: Mapping[str, Mapping[str, str]],
    fallback_fields: Sequence[str] | None, fallback_rows: Mapping[str, Mapping[str, str]],
) -> None:
    if fallback_fields is not None and list(base_fields) != list(fallback_fields):
        raise GeometryGateError("base and fallback allele-ratio schemas differ")
    rows: list[Mapping[str, str]] = []
    for barcode in order:
        if selected_rate_status[barcode] != "finite":
            continue
        endpoint = selected_endpoint[barcode]
        row = base_rows.get(barcode) if endpoint == "base" else fallback_rows.get(barcode) if endpoint == "fallback" else None
        if row is not None:
            rows.append(row)
    _write_tsv(path, base_fields, rows)


def _validate_selected_staging(
    staging_prefix: Path, authoritative_order: Sequence[str], final_contract_path: Path,
) -> None:
    for suffix in _FINAL_SELECTED_SUFFIXES:
        path = Path(str(staging_prefix) + suffix)
        if not path.is_file() or path.stat().st_size <= 0:
            raise GeometryGateError(f"required staged selected output is missing: {path}")

    rate_order: list[str] = []
    with Path(str(staging_prefix) + ".contam_rate").open("r", encoding="utf-8") as handle:
        for raw in handle:
            if raw.strip():
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 3:
                    raise GeometryGateError("staged selected contamination rate is malformed")
                rate_order.append(fields[0])
    _, diag_rows, _ = _read_tsv_rows(Path(str(staging_prefix) + ".contam_diagnostics.tsv"))
    _, audit_rows, _ = _read_tsv_rows(Path(str(staging_prefix) + ".geometry_gate_audit.tsv"))
    diag_order = [row["barcode"] for row in diag_rows]
    audit_order = [row["barcode"] for row in audit_rows]
    expected = list(authoritative_order)
    if rate_order != expected or diag_order != expected or audit_order != expected:
        raise GeometryGateError("selected rate/diagnostic/audit roster or order mismatch")

    try:
        parsed = json.loads(final_contract_path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise GeometryGateError(f"invalid staged final run contract: {exc}") from exc
    if not isinstance(parsed, dict) or parsed.get("production_contract_pass") is not True:
        raise GeometryGateError("staged final run contract is not a passing JSON object")


def _publish_staging(
    staging_prefix: Path, outer_prefix: Path, staged_suffixes: Sequence[str],
    staged_contract: Path,
) -> None:
    moved: list[Path] = []
    try:
        for suffix in staged_suffixes:
            source = Path(str(staging_prefix) + suffix)
            destination = Path(str(outer_prefix) + suffix)
            os.replace(source, destination)
            moved.append(destination)
        _fsync_directory(outer_prefix.parent)
        destination_contract = Path(str(outer_prefix) + ".run_contract.json")
        os.replace(staged_contract, destination_contract)
        _fsync_directory(outer_prefix.parent)
    except Exception as exc:
        for path in moved:
            try:
                path.unlink()
            except FileNotFoundError:
                pass
        try:
            Path(str(outer_prefix) + ".run_contract.json").unlink()
        except FileNotFoundError:
            pass
        raise GeometryGateError(f"atomic publication failed: {exc}") from exc


def _write_final_contract(
    path: Path, *, outer_prefix: Path, condition_key: str, gate_version: str,
    base_strength: float, fallback_strength: float, alpha_threshold: float,
    orthogonal_threshold: float, parent_mass_threshold: float,
    total_cells: int, evaluable_cells: int, triggered_cells: int,
    base_rates: Mapping[str, Mapping[str, object]],
    fallback_rates: Mapping[str, Mapping[str, object]],
    fallback_required: bool, fallback_evaluated: bool,
    summary_counts: Mapping[str, int], reason_counts: Mapping[str, int],
    base_contract_path: Path, fallback_contract_path: Path | None,
    base_contract: Mapping[str, object], fallback_contract: Mapping[str, object] | None,
) -> None:
    base_finite = sum(1 for state in base_rates.values() if state["status"] == "finite")
    base_failed = total_cells - base_finite
    fallback_finite = sum(1 for state in fallback_rates.values() if state["status"] == "finite") if fallback_evaluated else 0
    fallback_failed = sum(1 for state in fallback_rates.values() if state["status"] == "failed") if fallback_evaluated else 0
    selected_na = int(summary_counts.get("selected_na_cells", 0))
    contract: dict[str, object] = {
        "contract_version": CONTRACT_VERSION,
        "tool": "geometry_gated_contam_estimate",
        "tool_version": WORKER_VERSION,
        "condition_key": condition_key,
        "output_prefix": str(outer_prefix),
        "run_class": "validation_candidate",
        "panel_mode": "interindividual",
        "geometry_gate_version": gate_version,
        "base_source_exclusion_strength": base_strength,
        "fallback_source_exclusion_strength": fallback_strength,
        "source_exclusion_strength_selected": "per_cell_geometry_gate",
        "geometry_gate_parent_axis_alpha_threshold": alpha_threshold,
        "geometry_gate_ambient_orthogonal_norm_threshold": orthogonal_threshold,
        "geometry_gate_parent_mass_threshold": parent_mass_threshold,
        "published_profile_basis": "base_endpoint_fitted_profile",
        "rate_roster_basis": "base_endpoint_contam_diagnostics",
        "failed_rate_encoding": "nan",
        "per_cell_fit_failures_nonfatal": True,
        "fallback_endpoint_required_by_gate_or_rescue": bool(fallback_required),
        "fallback_endpoint_evaluated": bool(fallback_evaluated),
        "geometry_gate_total_cells": int(total_cells),
        "geometry_gate_evaluable_cells": int(evaluable_cells),
        "geometry_gate_triggered_cells": int(triggered_cells),
        "geometry_gate_triggered_fraction": float(triggered_cells / total_cells) if total_cells else 0.0,
        "base_rate_finite_cells": int(base_finite),
        "base_rate_failed_cells": int(base_failed),
        "fallback_rate_finite_cells": int(fallback_finite),
        "fallback_rate_failed_cells": int(fallback_failed),
        "selected_base_cells": int(summary_counts.get("selected_base_cells", 0)),
        "selected_fallback_geometry_cells": int(summary_counts.get("selected_fallback_geometry_cells", 0)),
        "selected_fallback_rescue_cells": int(summary_counts.get("selected_fallback_rescue_cells", 0)),
        "selected_na_cells": selected_na,
        "selected_na_fraction": float(selected_na / total_cells) if total_cells else 0.0,
        "geometry_gate_selection_reason_counts": {
            str(key): int(value) for key, value in sorted(reason_counts.items())
        },
        "base_endpoint_run_contract_sha256": _sha256(base_contract_path),
        "fallback_endpoint_run_contract_sha256": _sha256(fallback_contract_path) if fallback_contract_path else None,
        "base_endpoint_contract_summary": _contract_summary(base_contract),
        "fallback_endpoint_contract_summary": _contract_summary(fallback_contract) if fallback_contract is not None else None,
        "production_contract_pass": True,
        "production_contract_reason": (
            "geometry_gate_outputs_validated_with_nonfatal_cell_fit_failures"
            if selected_na > 0
            else "geometry_gate_outputs_validated"
        ),
    }
    _atomic_write_text(path, json.dumps(contract, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _copy_base_outputs_to_staging(base_prefix: Path, staging_prefix: Path) -> list[str]:
    staged = []
    for suffix in (".contam_prof", ".profile_fit_diagnostics.tsv", ".condf_coverage.tsv"):
        _atomic_copy(Path(str(base_prefix) + suffix), Path(str(staging_prefix) + suffix))
        staged.append(suffix)
    for suffix in _OPTIONAL_BASE_SUFFIXES:
        source = Path(str(base_prefix) + suffix)
        if source.is_file() and source.stat().st_size > 0:
            _atomic_copy(source, Path(str(staging_prefix) + suffix))
            staged.append(suffix)
    return staged


def run(args: argparse.Namespace, estimator_arguments: Sequence[str]) -> None:
    estimator_binary = Path(args.estimator_binary)
    if not estimator_binary.is_file() or not os.access(estimator_binary, os.X_OK):
        raise GeometryGateError(
            f"estimator binary is missing or not executable: {estimator_binary}"
        )
    _assert_no_active_forwarded_inputs(estimator_arguments)

    outer_prefix = Path(_option_value(estimator_arguments, ("-o", "--output_prefix")))
    requested_contract = Path(_option_value(estimator_arguments, ("--run_contract",)))
    expected_final_contract = Path(str(outer_prefix) + ".run_contract.json")
    if requested_contract != expected_final_contract:
        # The orchestrator already uses the canonical path; direct callers may
        # supply an equivalent path spelling but not redirect the completeness marker.
        raise GeometryGateError(
            f"geometry helper requires --run_contract {expected_final_contract}, got {requested_contract}"
        )

    outer_prefix.parent.mkdir(parents=True, exist_ok=True)
    base_prefix = Path(str(outer_prefix) + ".geometry_base_endpoint")
    fallback_prefix = Path(str(outer_prefix) + ".geometry_fallback_endpoint")
    staging_prefix = outer_prefix.with_name(
        f".{outer_prefix.name}.geometry_publish.{os.getpid()}"
    )
    staged_contract = Path(str(staging_prefix) + ".run_contract.json")

    # Retry-safe startup cleanup. Never touch upstream input suffixes.
    _cleanup_prefix(base_prefix)
    _cleanup_prefix(fallback_prefix)
    _cleanup_staging_for_outer(outer_prefix)
    _cleanup_final_outputs(outer_prefix)

    success = False
    try:
        _run_endpoint(
            estimator_binary, estimator_arguments, outer_prefix, base_prefix,
            args.base_strength, args.condition_key,
        )
        _validate_endpoint_required_files(base_prefix)

        base_diag_path = Path(str(base_prefix) + ".contam_diagnostics.tsv")
        base_fields, base_rows, base_by = _read_tsv_rows(base_diag_path)
        base_order = [str(row["barcode"]) for row in base_rows]
        base_rates = _read_rate_rows(
            Path(str(base_prefix) + ".contam_rate"), base_by, base_order
        )
        decisions = _evaluate_base_gate(
            base_fields, base_rows, base_rates,
            alpha_threshold=args.parent_axis_alpha_threshold,
            orthogonal_threshold=args.ambient_orthogonal_norm_threshold,
            parent_mass_threshold=args.parent_mass_threshold,
        )
        fallback_required = any(
            bool(decisions[bc]["triggered"]) or base_rates[bc]["status"] == "failed"
            for bc in base_order
        )

        base_contract_path = Path(str(base_prefix) + ".run_contract.json")
        base_contract = _load_endpoint_contract(
            base_contract_path, args.base_strength, args.condition_key
        )
        base_profile_fields, _ = _read_profile_diagnostics(
            Path(str(base_prefix) + ".profile_fit_diagnostics.tsv")
        )
        base_condf = Path(str(base_prefix) + ".condf_coverage.tsv")
        if not base_condf.is_file() or base_condf.stat().st_size <= 0:
            raise GeometryGateError(f"base endpoint condf coverage missing: {base_condf}")
        base_allele_fields, base_allele = _read_allele_rows(
            Path(str(base_prefix) + ".allele_ratio"), set(base_order)
        )

        fallback_fields: list[str] | None = None
        fallback_by: dict[str, dict[str, str]] = {}
        fallback_rates: dict[str, dict[str, object]] = {}
        fallback_contract_path: Path | None = None
        fallback_contract: dict[str, object] | None = None
        fallback_allele_fields: list[str] | None = None
        fallback_allele: dict[str, dict[str, str]] = {}

        if fallback_required:
            _run_endpoint(
                estimator_binary, estimator_arguments, outer_prefix, fallback_prefix,
                args.fallback_strength, args.condition_key,
            )
            _validate_endpoint_required_files(fallback_prefix)
            fallback_fields, fallback_rows, fallback_by = _read_tsv_rows(
                Path(str(fallback_prefix) + ".contam_diagnostics.tsv")
            )
            if fallback_fields != base_fields:
                raise GeometryGateError("base and fallback diagnostic schemas differ")
            fallback_order = [str(row["barcode"]) for row in fallback_rows]
            if set(fallback_order) != set(base_order) or len(fallback_order) != len(base_order):
                raise GeometryGateError("base and fallback diagnostic barcode rosters differ")
            fallback_rates = _read_rate_rows(
                Path(str(fallback_prefix) + ".contam_rate"), fallback_by, base_order
            )
            fallback_contract_path = Path(str(fallback_prefix) + ".run_contract.json")
            fallback_contract = _load_endpoint_contract(
                fallback_contract_path, args.fallback_strength, args.condition_key
            )
            _validate_contract_invariants(base_contract, fallback_contract)
            fallback_profile_fields, _ = _read_profile_diagnostics(
                Path(str(fallback_prefix) + ".profile_fit_diagnostics.tsv")
            )
            if fallback_profile_fields != base_profile_fields:
                raise GeometryGateError("base and fallback profile-fit diagnostic schemas differ")
            fallback_condf = Path(str(fallback_prefix) + ".condf_coverage.tsv")
            if not fallback_condf.is_file() or fallback_condf.stat().st_size <= 0:
                raise GeometryGateError(f"fallback endpoint condf coverage missing: {fallback_condf}")
            if base_condf.read_bytes() != fallback_condf.read_bytes():
                raise GeometryGateError("base and fallback .condf_coverage.tsv files differ")
            fallback_allele_fields, fallback_allele = _read_allele_rows(
                Path(str(fallback_prefix) + ".allele_ratio"), set(base_order)
            )
            if fallback_allele_fields != base_allele_fields:
                raise GeometryGateError("base and fallback allele-ratio schemas differ")

        final_fields, selected_rows, audit_rows, selected_endpoint, reason_counts, summary_counts = _build_selected_rows(
            base_fields, base_rows, base_rates, fallback_by, fallback_rates, decisions,
            fallback_evaluated=fallback_required,
            gate_version=args.gate_version,
            base_strength=args.base_strength,
            fallback_strength=args.fallback_strength,
            alpha_threshold=args.parent_axis_alpha_threshold,
            orthogonal_threshold=args.ambient_orthogonal_norm_threshold,
            parent_mass_threshold=args.parent_mass_threshold,
        )
        selected_by = {str(row["barcode"]): row for row in selected_rows}
        selected_rate_status = {
            bc: str(selected_by[bc]["selected_rate_status"]) for bc in base_order
        }

        _write_selected_rate(
            Path(str(staging_prefix) + ".contam_rate"), base_order, selected_by
        )
        _write_tsv(
            Path(str(staging_prefix) + ".contam_diagnostics.tsv"), final_fields, selected_rows
        )
        _write_tsv(
            Path(str(staging_prefix) + ".geometry_gate_audit.tsv"),
            ["barcode", *_AUDIT_FIELDS], audit_rows,
        )
        _write_selected_allele(
            Path(str(staging_prefix) + ".allele_ratio"), base_order,
            selected_endpoint, selected_rate_status,
            base_allele_fields, base_allele,
            fallback_allele_fields if fallback_required else None,
            fallback_allele,
        )
        staged_suffixes = [
            ".contam_rate", ".contam_diagnostics.tsv", ".geometry_gate_audit.tsv",
            ".allele_ratio",
        ]
        staged_suffixes.extend(_copy_base_outputs_to_staging(base_prefix, staging_prefix))

        _write_final_contract(
            staged_contract,
            outer_prefix=outer_prefix,
            condition_key=args.condition_key,
            gate_version=args.gate_version,
            base_strength=args.base_strength,
            fallback_strength=args.fallback_strength,
            alpha_threshold=args.parent_axis_alpha_threshold,
            orthogonal_threshold=args.ambient_orthogonal_norm_threshold,
            parent_mass_threshold=args.parent_mass_threshold,
            total_cells=len(base_order),
            evaluable_cells=sum(bool(decisions[bc]["evaluable"]) for bc in base_order),
            triggered_cells=sum(bool(decisions[bc]["triggered"]) for bc in base_order),
            base_rates=base_rates,
            fallback_rates=fallback_rates,
            fallback_required=fallback_required,
            fallback_evaluated=fallback_required,
            summary_counts=summary_counts,
            reason_counts=reason_counts,
            base_contract_path=base_contract_path,
            fallback_contract_path=fallback_contract_path,
            base_contract=base_contract,
            fallback_contract=fallback_contract,
        )
        _validate_selected_staging(staging_prefix, base_order, staged_contract)
        _publish_staging(staging_prefix, outer_prefix, staged_suffixes, staged_contract)
        success = True
        print(
            f"Geometry gate completed {len(base_order)} cells; "
            f"fallback evaluated={fallback_required}; "
            f"selected NA={summary_counts.get('selected_na_cells', 0)}.",
            file=sys.stderr, flush=True,
        )
    finally:
        _cleanup_prefix(staging_prefix)
        if success:
            _cleanup_prefix(base_prefix)
            _cleanup_prefix(fallback_prefix)
        else:
            _cleanup_final_outputs(outer_prefix)
            print(
                "Geometry-gated run failed; retained endpoint evidence for diagnosis:\n"
                f"  base: {base_prefix}\n"
                f"  fallback: {fallback_prefix}",
                file=sys.stderr, flush=True,
            )


def parse_args(argv: Sequence[str] | None = None) -> tuple[argparse.Namespace, list[str]]:
    parser = argparse.ArgumentParser(
        description="Internal production CK geometry-gated endpoint selector"
    )
    parser.add_argument("--estimator-binary", required=True)
    parser.add_argument("--condition-key", required=True)
    parser.add_argument("--gate-version", required=True)
    parser.add_argument("--base-strength", type=float, required=True)
    parser.add_argument("--fallback-strength", type=float, required=True)
    parser.add_argument("--parent-axis-alpha-threshold", type=float, required=True)
    parser.add_argument("--ambient-orthogonal-norm-threshold", type=float, required=True)
    parser.add_argument("--parent-mass-threshold", type=float, required=True)
    parser.add_argument("estimator_arguments", nargs=argparse.REMAINDER)
    parsed = parser.parse_args(argv)
    forwarded = list(parsed.estimator_arguments)
    if forwarded and forwarded[0] == "--":
        forwarded = forwarded[1:]
    if not forwarded:
        parser.error("estimator arguments are required after --")
    for name in (
        "base_strength", "fallback_strength", "parent_axis_alpha_threshold",
        "ambient_orthogonal_norm_threshold", "parent_mass_threshold",
    ):
        value = float(getattr(parsed, name))
        if not math.isfinite(value):
            parser.error(f"{name.replace('_', '-')} must be finite")
    if not (0.0 <= parsed.base_strength <= 1.0):
        parser.error("base-strength must be in [0,1]")
    if not (0.0 <= parsed.fallback_strength <= 1.0):
        parser.error("fallback-strength must be in [0,1]")
    if math.isclose(parsed.base_strength, parsed.fallback_strength):
        parser.error("base-strength and fallback-strength must differ")
    return parsed, forwarded


def main(argv: Sequence[str] | None = None) -> int:
    try:
        args, estimator_arguments = parse_args(argv)
        run(args, estimator_arguments)
    except GeometryGateError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
