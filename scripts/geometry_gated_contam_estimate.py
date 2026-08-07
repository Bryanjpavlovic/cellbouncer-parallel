#!/usr/bin/env python3
"""Run the frozen CK geometry-gated source-exclusion condition on real libraries.

This production helper is invoked by orchestrate_pipeline.py. It runs the
unchanged tet_contam_estimate binary at lambda=0.25, evaluates the frozen
CK_GEOMETRY_GATE_V1 from estimator-emitted base-endpoint geometry, lazily runs
lambda=1.0 only when at least one barcode triggers, and publishes one selected
per-cell result plus an explicit gate audit.
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

WORKER_VERSION = "geometry_gated_contam_estimate_V1"


class GeometryGateError(RuntimeError):
    """Raised when endpoint evidence cannot support deterministic selection."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)


def _atomic_copy(source: Path, destination: Path) -> None:
    if not source.is_file() or source.stat().st_size <= 0:
        raise GeometryGateError(f"required endpoint output is missing: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    with source.open("rb") as src, tmp.open("wb") as dst:
        shutil.copyfileobj(src, dst, length=1024 * 1024)
        dst.flush()
        os.fsync(dst.fileno())
    os.replace(tmp, destination)


def _option_value(arguments: Sequence[str], names: Iterable[str]) -> str:
    ordered = tuple(names)
    accepted = set(ordered)
    for index, token in enumerate(arguments):
        if token in accepted:
            if index + 1 >= len(arguments):
                raise GeometryGateError(f"option {token} is missing its value")
            return str(arguments[index + 1])
        for name in ordered:
            prefix = name + "="
            if token.startswith(prefix):
                return token[len(prefix):]
    raise GeometryGateError(f"required option is absent: {sorted(accepted)}")


def _replace_option(arguments: Sequence[str], names: Iterable[str], value: str) -> list[str]:
    ordered = tuple(names)
    accepted = set(ordered)
    output: list[str] = []
    replaced = False
    index = 0
    while index < len(arguments):
        token = str(arguments[index])
        if token in accepted:
            if index + 1 >= len(arguments):
                raise GeometryGateError(f"option {token} is missing its value")
            output.extend([token, value])
            replaced = True
            index += 2
            continue
        matched = next((name for name in ordered if token.startswith(name + "=")), None)
        if matched is not None:
            output.append(f"{matched}={value}")
            replaced = True
            index += 1
            continue
        output.append(token)
        index += 1
    if not replaced:
        output.extend([ordered[0], value])
    return output


def _remove_option(arguments: Sequence[str], names: Iterable[str]) -> list[str]:
    ordered = tuple(names)
    accepted = set(ordered)
    output: list[str] = []
    index = 0
    while index < len(arguments):
        token = str(arguments[index])
        if token in accepted:
            if index + 1 >= len(arguments):
                raise GeometryGateError(f"option {token} is missing its value")
            index += 2
            continue
        if any(token.startswith(name + "=") for name in ordered):
            index += 1
            continue
        output.append(token)
        index += 1
    return output


def _endpoint_arguments(
    estimator_arguments: Sequence[str], endpoint_prefix: Path, strength: float,
) -> list[str]:
    arguments = _replace_option(
        estimator_arguments, ("-o", "--output_prefix"), str(endpoint_prefix)
    )
    arguments = _replace_option(
        arguments, ("--run_contract",), str(endpoint_prefix) + ".run_contract.json"
    )
    arguments = _remove_option(arguments, ("--source_exclusion_strength", "--loo"))
    arguments.extend(["--source_exclusion_strength", format(strength, ".17g")])
    return arguments


def _link_endpoint_inputs(outer_prefix: Path, endpoint_prefix: Path) -> None:
    for suffix in (".counts", ".condf", ".samples", ".assignments"):
        source = Path(str(outer_prefix) + suffix)
        if not source.is_file() or source.stat().st_size <= 0:
            raise GeometryGateError(f"required staged estimator input is missing: {source}")
        destination = Path(str(endpoint_prefix) + suffix)
        try:
            os.link(source, destination)
        except OSError:
            shutil.copy2(source, destination)


def _cleanup_endpoint(prefix: Path) -> None:
    for path in prefix.parent.glob(prefix.name + ".*"):
        try:
            if path.is_file() or path.is_symlink():
                path.unlink()
        except FileNotFoundError:
            pass


def _run_endpoint(
    estimator_binary: Path,
    estimator_arguments: Sequence[str],
    outer_prefix: Path,
    endpoint_prefix: Path,
    strength: float,
) -> None:
    _link_endpoint_inputs(outer_prefix, endpoint_prefix)
    command = [
        str(estimator_binary),
        *_endpoint_arguments(estimator_arguments, endpoint_prefix, strength),
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


def _parse_rate_number(raw: str, *, allow_nan: bool, context: str) -> float:
    try:
        value = float(raw)
    except ValueError as exc:
        raise GeometryGateError(f"nonnumeric {context}: {raw!r}") from exc
    if math.isinf(value):
        raise GeometryGateError(f"infinite {context}: {raw!r}")
    if math.isnan(value) and not allow_nan:
        raise GeometryGateError(f"NaN {context}: {raw!r}")
    return value


def _read_rate_rows(path: Path) -> tuple[list[str], dict[str, tuple[str, str]]]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise GeometryGateError(f"endpoint contamination-rate file is missing: {path}")
    order: list[str] = []
    values: dict[str, tuple[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise GeometryGateError(
                    f"malformed contamination-rate row {line_number} in {path}"
                )
            barcode = fields[0]
            c_value = fields[1]
            c_se = fields[2] if len(fields) >= 3 else "NaN"
            if barcode in values:
                raise GeometryGateError(f"duplicate barcode {barcode!r} in {path}")
            _parse_rate_number(
                c_value, allow_nan=False,
                context=f"contamination estimate at row {line_number} in {path}",
            )
            _parse_rate_number(
                c_se, allow_nan=True,
                context=f"contamination SE at row {line_number} in {path}",
            )
            order.append(barcode)
            values[barcode] = (c_value, c_se)
    if not values:
        raise GeometryGateError(f"endpoint contamination-rate file has no rows: {path}")
    return order, values


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


def _finite_float(row: Mapping[str, str], field: str) -> float:
    raw = str(row.get(field, ""))
    try:
        value = float(raw)
    except ValueError as exc:
        raise GeometryGateError(f"gate diagnostic {field} is not numeric: {raw!r}") from exc
    return value


def _strict_bool(raw: object, field: str) -> bool:
    if isinstance(raw, bool):
        return raw
    value = str(raw).strip().lower()
    if value in {"1", "1.0", "true"}:
        return True
    if value in {"0", "0.0", "false"}:
        return False
    raise GeometryGateError(f"gate diagnostic {field} is not boolean: {raw!r}")


def _selection_reason(
    *,
    triggered: bool,
    finite_geometry: bool,
    is_heterotypic: bool,
    alpha: float,
    orthogonal: float,
    parent_mass: float,
    alpha_threshold: float,
    orthogonal_threshold: float,
    parent_mass_threshold: float,
) -> str:
    if triggered:
        return "geometry_gate_triggered_select_fallback"
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
    base_rows: Sequence[Mapping[str, str]],
    *,
    alpha_threshold: float,
    orthogonal_threshold: float,
    parent_mass_threshold: float,
) -> dict[str, dict[str, object]]:
    required = {
        "barcode",
        "is_heterotypic",
        "ambient_parent_axis_alpha_scoring",
        "ambient_orthogonal_norm_scoring",
        "excluded_parent_A_mass_raw",
        "excluded_parent_B_mass_raw",
        "parent_axis_geometry_status",
    }
    if base_rows:
        missing = sorted(required - set(base_rows[0]))
        if missing:
            raise GeometryGateError(
                "base endpoint diagnostics lack gate fields: " + ",".join(missing)
            )

    decisions: dict[str, dict[str, object]] = {}
    for base in base_rows:
        barcode = str(base["barcode"])
        is_heterotypic = _strict_bool(base.get("is_heterotypic", ""), "is_heterotypic")
        alpha = _finite_float(base, "ambient_parent_axis_alpha_scoring")
        orthogonal = _finite_float(base, "ambient_orthogonal_norm_scoring")
        raw_a = _finite_float(base, "excluded_parent_A_mass_raw")
        raw_b = _finite_float(base, "excluded_parent_B_mass_raw")
        if not (math.isfinite(raw_a) and math.isfinite(raw_b)):
            raise GeometryGateError(
                f"barcode {barcode!r} has malformed raw excluded-parent masses: "
                f"A={raw_a!r}, B={raw_b!r}"
            )
        parent_mass = max(raw_a, raw_b)
        finite_geometry = all(math.isfinite(v) for v in (alpha, orthogonal, raw_a, raw_b))
        geometry_status = str(base.get("parent_axis_geometry_status", ""))
        if is_heterotypic and (not finite_geometry or geometry_status != "ok"):
            raise GeometryGateError(
                f"heterotypic barcode {barcode!r} has malformed/missing base geometry: "
                f"status={geometry_status!r}, alpha={alpha!r}, orthogonal={orthogonal!r}, "
                f"raw_parent_A={raw_a!r}, raw_parent_B={raw_b!r}"
            )
        triggered = bool(
            finite_geometry
            and alpha >= alpha_threshold
            and orthogonal <= orthogonal_threshold
            and parent_mass >= parent_mass_threshold
        )
        decisions[barcode] = {
            "triggered": triggered,
            "reason": _selection_reason(
                triggered=triggered,
                finite_geometry=finite_geometry,
                is_heterotypic=is_heterotypic,
                alpha=alpha,
                orthogonal=orthogonal,
                parent_mass=parent_mass,
                alpha_threshold=alpha_threshold,
                orthogonal_threshold=orthogonal_threshold,
                parent_mass_threshold=parent_mass_threshold,
            ),
            "alpha": alpha,
            "orthogonal": orthogonal,
            "parent_mass": parent_mass,
            "is_heterotypic": is_heterotypic,
        }
    return decisions


def _build_gate_selection(
    base_rows: Sequence[Mapping[str, str]],
    fallback_by_barcode: Mapping[str, Mapping[str, str]],
    base_rates: Mapping[str, tuple[str, str]],
    fallback_rates: Mapping[str, tuple[str, str]],
    *,
    gate_version: str,
    alpha_threshold: float,
    orthogonal_threshold: float,
    parent_mass_threshold: float,
    base_strength: float,
    fallback_strength: float,
    decisions: Mapping[str, Mapping[str, object]],
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, bool], Counter[str]]:
    selected_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []
    triggered_by_barcode: dict[str, bool] = {}
    reason_counts: Counter[str] = Counter()
    gate_fields = [
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

    for base in base_rows:
        barcode = str(base["barcode"])
        decision = decisions.get(barcode)
        if decision is None or barcode not in base_rates:
            raise GeometryGateError(f"base gate evidence is incomplete for barcode {barcode!r}")
        triggered = bool(decision["triggered"])
        fallback = fallback_by_barcode.get(barcode)
        fallback_rate = fallback_rates.get(barcode)
        if triggered and (fallback is None or fallback_rate is None):
            raise GeometryGateError(
                f"triggered barcode {barcode!r} lacks the fallback endpoint"
            )
        selected_rate = fallback_rate if triggered else base_rates[barcode]
        selected = dict(fallback if triggered else base)
        for field in gate_fields:
            selected.pop(field, None)

        gate_values = {
            "geometry_gate_version": gate_version,
            "geometry_gate_triggered": "1" if triggered else "0",
            "geometry_gate_parent_axis_alpha": format(float(decision["alpha"]), ".17g"),
            "geometry_gate_ambient_orthogonal_norm": format(float(decision["orthogonal"]), ".17g"),
            "geometry_gate_max_raw_excluded_parent_mass": format(float(decision["parent_mass"]), ".17g"),
            "geometry_gate_parent_axis_alpha_threshold": format(alpha_threshold, ".17g"),
            "geometry_gate_ambient_orthogonal_norm_threshold": format(orthogonal_threshold, ".17g"),
            "geometry_gate_parent_mass_threshold": format(parent_mass_threshold, ".17g"),
            "source_exclusion_strength_base": format(base_strength, ".17g"),
            "source_exclusion_strength_fallback": format(fallback_strength, ".17g"),
            "source_exclusion_strength_selected": format(
                fallback_strength if triggered else base_strength, ".17g"
            ),
            "geometry_gate_selection_reason": str(decision["reason"]),
            "base_c_selected": base_rates[barcode][0],
            "fallback_c_selected": fallback_rate[0] if fallback_rate is not None else "",
            "selected_c": selected_rate[0],
        }
        selected.update(gate_values)
        selected_rows.append(selected)
        audit_rows.append({"barcode": barcode, **gate_values})
        triggered_by_barcode[barcode] = triggered
        reason_counts[str(decision["reason"])] += 1

    return selected_rows, audit_rows, triggered_by_barcode, reason_counts


def _write_rate_selection(
    path: Path,
    order: Sequence[str],
    base_rates: Mapping[str, tuple[str, str]],
    fallback_rates: Mapping[str, tuple[str, str]],
    triggered: Mapping[str, bool],
) -> None:
    lines: list[str] = []
    for barcode in order:
        selected = (
            fallback_rates[barcode]
            if triggered.get(barcode, False)
            else base_rates[barcode]
        )
        lines.append(f"{barcode}\t{selected[0]}\t{selected[1]}\n")
    _atomic_write_text(path, "".join(lines))


def _write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fieldnames),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)


def _merge_barcode_tsv(
    base_path: Path,
    fallback_path: Path,
    output_path: Path,
    triggered: Mapping[str, bool],
) -> None:
    base_fields, base_rows, _ = _read_tsv_rows(base_path, allow_empty=True)
    fallback_fields, fallback_rows, fallback_by = _read_tsv_rows(
        fallback_path, allow_empty=True
    )
    if base_fields != fallback_fields:
        raise GeometryGateError(f"endpoint TSV schemas differ: {base_path.name}")
    base_roster = {str(row["barcode"]) for row in base_rows}
    if base_roster != {str(row["barcode"]) for row in fallback_rows}:
        raise GeometryGateError(f"endpoint TSV barcode rosters differ: {base_path.name}")
    selected = [
        dict(fallback_by[str(row["barcode"])])
        if triggered.get(str(row["barcode"]), False)
        else dict(row)
        for row in base_rows
    ]
    _write_tsv(output_path, base_fields, selected)


def _load_endpoint_contract(path: Path, expected_strength: float) -> dict[str, object]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise GeometryGateError(f"endpoint run contract is missing: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise GeometryGateError(f"cannot read endpoint run contract {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise GeometryGateError(f"endpoint run contract is not an object: {path}")
    if value.get("production_contract_pass") is False:
        raise GeometryGateError(f"endpoint production contract failed: {path}")
    for field, expected in (
        ("run_class", "production"),
        ("assignments_basis", "library"),
        ("expected_lines_basis", "library"),
        ("ambient_candidates_basis", "library"),
    ):
        observed = value.get(field)
        if observed is not None and str(observed) != expected:
            raise GeometryGateError(
                f"endpoint run contract {path} has {field}={observed!r}, "
                f"expected {expected!r}"
            )
    if value.get("strict_condf") is False:
        raise GeometryGateError(f"endpoint run contract did not enable strict_condf: {path}")
    if value.get("truth_assisted") is True:
        raise GeometryGateError(f"endpoint run contract is truth-assisted: {path}")
    if "source_exclusion_strength" in value:
        try:
            observed_strength = float(value["source_exclusion_strength"])
        except (TypeError, ValueError) as exc:
            raise GeometryGateError(
                f"endpoint run contract has invalid source-exclusion strength: {path}"
            ) from exc
        if not math.isclose(
            observed_strength, expected_strength, rel_tol=0.0, abs_tol=1e-12
        ):
            raise GeometryGateError(
                f"endpoint run-contract strength mismatch: {observed_strength} "
                f"!= {expected_strength}: {path}"
            )
    return value


def _write_final_contract(
    path: Path,
    *,
    outer_prefix: Path,
    condition_key: str,
    gate_version: str,
    base_strength: float,
    fallback_strength: float,
    alpha_threshold: float,
    orthogonal_threshold: float,
    parent_mass_threshold: float,
    total_cells: int,
    triggered_cells: int,
    reason_counts: Mapping[str, int],
    base_contract_path: Path,
    fallback_contract_path: Path | None,
) -> None:
    base_contract = _load_endpoint_contract(base_contract_path, base_strength)
    fallback_contract = (
        _load_endpoint_contract(fallback_contract_path, fallback_strength)
        if fallback_contract_path is not None
        else None
    )
    contract: dict[str, object] = {
        "contract_version": "geometry_gated_contam_estimate_run_contract_V1",
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
        "geometry_gate_total_cells": int(total_cells),
        "geometry_gate_triggered_cells": int(triggered_cells),
        "geometry_gate_triggered_fraction": (
            float(triggered_cells / total_cells) if total_cells else 0.0
        ),
        "geometry_gate_selection_reason_counts": {
            str(key): int(value) for key, value in sorted(reason_counts.items())
        },
        "fallback_endpoint_evaluated": fallback_contract_path is not None,
        "published_profile_basis": "base_endpoint_fitted_profile",
        "production_contract_pass": True,
        "production_contract_reason": "geometry_gate_outputs_validated",
        "base_endpoint_run_contract_sha256": (
            _sha256(base_contract_path) if base_contract_path.is_file() else None
        ),
        "fallback_endpoint_run_contract_sha256": (
            _sha256(fallback_contract_path)
            if fallback_contract_path is not None and fallback_contract_path.is_file()
            else None
        ),
        "base_endpoint_contract": base_contract,
        "fallback_endpoint_contract": fallback_contract,
    }
    _atomic_write_text(
        path,
        json.dumps(contract, indent=2, sort_keys=True, allow_nan=False) + "\n",
    )


def _copy_optional_base_outputs(base_prefix: Path, outer_prefix: Path) -> None:
    for suffix in (
        ".contam_prof",
        ".profile_fit_diagnostics.tsv",
        ".condf_coverage.tsv",
        ".model_fit.tsv",
        ".class_residuals.tsv",
        ".class_residuals.tsv.gz",
        ".decontam.assignments",
        ".pass1.contam_rate",
        ".pass1.contam_prof",
    ):
        source = Path(str(base_prefix) + suffix)
        if source.is_file() and source.stat().st_size > 0:
            _atomic_copy(source, Path(str(outer_prefix) + suffix))


def run(args: argparse.Namespace, estimator_arguments: Sequence[str]) -> None:
    estimator_binary = Path(args.estimator_binary)
    if not estimator_binary.is_file() or not os.access(estimator_binary, os.X_OK):
        raise GeometryGateError(
            f"estimator binary is missing or not executable: {estimator_binary}"
        )

    outer_prefix = Path(_option_value(estimator_arguments, ("-o", "--output_prefix")))
    final_run_contract = Path(
        _option_value(estimator_arguments, ("--run_contract",))
    )
    outer_prefix.parent.mkdir(parents=True, exist_ok=True)
    base_prefix = Path(str(outer_prefix) + ".geometry_base_endpoint")
    fallback_prefix = Path(str(outer_prefix) + ".geometry_fallback_endpoint")
    for prefix in (base_prefix, fallback_prefix):
        _cleanup_endpoint(prefix)

    _run_endpoint(
        estimator_binary,
        estimator_arguments,
        outer_prefix,
        base_prefix,
        args.base_strength,
    )

    base_rate_order, base_rates = _read_rate_rows(
        Path(str(base_prefix) + ".contam_rate")
    )
    base_fields, base_diagnostic_rows, _ = _read_tsv_rows(
        Path(str(base_prefix) + ".contam_diagnostics.tsv")
    )
    if set(base_rates) != {str(row["barcode"]) for row in base_diagnostic_rows}:
        raise GeometryGateError(
            "base contamination-rate and diagnostic barcode rosters differ"
        )

    decisions = _evaluate_base_gate(
        base_diagnostic_rows,
        alpha_threshold=args.parent_axis_alpha_threshold,
        orthogonal_threshold=args.ambient_orthogonal_norm_threshold,
        parent_mass_threshold=args.parent_mass_threshold,
    )
    fallback_needed = any(bool(item["triggered"]) for item in decisions.values())

    fallback_rates: dict[str, tuple[str, str]] = {}
    fallback_diagnostic_by: dict[str, dict[str, str]] = {}
    if fallback_needed:
        _run_endpoint(
            estimator_binary,
            estimator_arguments,
            outer_prefix,
            fallback_prefix,
            args.fallback_strength,
        )
        fallback_order, fallback_rates = _read_rate_rows(
            Path(str(fallback_prefix) + ".contam_rate")
        )
        if base_rate_order != fallback_order:
            raise GeometryGateError(
                "base and fallback contamination-rate barcode order differs"
            )
        fallback_fields, fallback_rows, fallback_diagnostic_by = _read_tsv_rows(
            Path(str(fallback_prefix) + ".contam_diagnostics.tsv")
        )
        if base_fields != fallback_fields:
            raise GeometryGateError("base and fallback diagnostic schemas differ")
        if {str(row["barcode"]) for row in fallback_rows} != set(base_rates):
            raise GeometryGateError(
                "base and fallback diagnostic barcode rosters differ"
            )

    selected_diagnostics, audit_rows, triggered, reason_counts = _build_gate_selection(
        base_diagnostic_rows,
        fallback_diagnostic_by,
        base_rates,
        fallback_rates,
        gate_version=args.gate_version,
        alpha_threshold=args.parent_axis_alpha_threshold,
        orthogonal_threshold=args.ambient_orthogonal_norm_threshold,
        parent_mass_threshold=args.parent_mass_threshold,
        base_strength=args.base_strength,
        fallback_strength=args.fallback_strength,
        decisions=decisions,
    )

    appended_fields = [
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
    final_fields = [f for f in base_fields if f not in appended_fields] + appended_fields

    _write_rate_selection(
        Path(str(outer_prefix) + ".contam_rate"),
        base_rate_order,
        base_rates,
        fallback_rates,
        triggered,
    )
    _write_tsv(
        Path(str(outer_prefix) + ".contam_diagnostics.tsv"),
        final_fields,
        selected_diagnostics,
    )
    _write_tsv(
        Path(str(outer_prefix) + ".geometry_gate_audit.tsv"),
        ["barcode", *appended_fields],
        audit_rows,
    )

    base_allele = Path(str(base_prefix) + ".allele_ratio")
    final_allele = Path(str(outer_prefix) + ".allele_ratio")
    if fallback_needed:
        _merge_barcode_tsv(
            base_allele,
            Path(str(fallback_prefix) + ".allele_ratio"),
            final_allele,
            triggered,
        )
    else:
        _atomic_copy(base_allele, final_allele)

    _copy_optional_base_outputs(base_prefix, outer_prefix)
    for required_suffix in (
        ".contam_prof",
        ".profile_fit_diagnostics.tsv",
        ".condf_coverage.tsv",
    ):
        required = Path(str(outer_prefix) + required_suffix)
        if not required.is_file() or required.stat().st_size <= 0:
            raise GeometryGateError(f"required published output is missing: {required}")

    base_contract_path = Path(str(base_prefix) + ".run_contract.json")
    fallback_contract_path = (
        Path(str(fallback_prefix) + ".run_contract.json")
        if fallback_needed
        else None
    )
    _write_final_contract(
        final_run_contract,
        outer_prefix=outer_prefix,
        condition_key=args.condition_key,
        gate_version=args.gate_version,
        base_strength=args.base_strength,
        fallback_strength=args.fallback_strength,
        alpha_threshold=args.parent_axis_alpha_threshold,
        orthogonal_threshold=args.ambient_orthogonal_norm_threshold,
        parent_mass_threshold=args.parent_mass_threshold,
        total_cells=len(triggered),
        triggered_cells=sum(triggered.values()),
        reason_counts=reason_counts,
        base_contract_path=base_contract_path,
        fallback_contract_path=fallback_contract_path,
    )

    print(
        f"Geometry gate selected fallback for {sum(triggered.values())}/"
        f"{len(triggered)} cells.",
        file=sys.stderr,
        flush=True,
    )
    _cleanup_endpoint(base_prefix)
    _cleanup_endpoint(fallback_prefix)


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
        "base_strength",
        "fallback_strength",
        "parent_axis_alpha_threshold",
        "ambient_orthogonal_norm_threshold",
        "parent_mass_threshold",
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
