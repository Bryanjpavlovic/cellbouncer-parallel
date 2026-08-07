#!/usr/bin/env python3
"""Structural identifiability pre-flight for tetraploid ambient RNA profiles.

The utility uses only production-available design inputs: CellBouncer conditional
expectation matrices, source/sample ordering, receiver classes, optional species
matrices/maps, and optional design-row weights.  It does not estimate contamination
and never consumes synthetic truth or a fitted ambient profile.

The module exposes :func:`run_preflight` for pipeline integration and a CLI for
standalone use.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

SCHEMA_VERSION = "tet_ambient_identifiability_preflight_V2"


def _json_safe(value):
    if isinstance(value, Mapping):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, (float, np.floating)) and not math.isfinite(float(value)):
        return None
    if isinstance(value, np.integer):
        return int(value)
    return value


@dataclass(frozen=True)
class ReceiverClass:
    name: str
    parent_1: str
    parent_2: str


def _sha256(path: Optional[Path]) -> str:
    if path is None:
        return ""
    p = Path(path)
    if not p.is_file():
        return ""
    h = hashlib.sha256()
    with p.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return "sha256:" + h.hexdigest()


def _fingerprint(path: Optional[Path]) -> Dict[str, object]:
    if path is None:
        return {"path": "", "exists": False, "size_bytes": 0, "sha256": ""}
    p = Path(path).resolve()
    return {
        "path": str(p),
        "exists": p.is_file(),
        "size_bytes": int(p.stat().st_size) if p.is_file() else 0,
        "sha256": _sha256(p),
    }


def _atomic_write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(_json_safe(dict(payload)), indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def _write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    os.replace(tmp, path)


def _read_samples(path: Path) -> List[str]:
    samples: List[str] = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            value = line.strip().split("\t", 1)[0]
            if value and not value.startswith("#"):
                samples.append(value)
    if not samples:
        raise ValueError(f"no sample labels found in {path}")
    if len(set(samples)) != len(samples):
        duplicates = sorted({x for x in samples if samples.count(x) > 1})
        raise ValueError(f"duplicate sample labels in {path}: {duplicates[:10]}")
    return samples


def _parse_source_labels(raw: Optional[str | Sequence[str]], samples: Sequence[str]) -> List[str]:
    if raw is None:
        labels = list(samples)
    elif isinstance(raw, str):
        labels = [x.strip() for x in raw.split(",") if x.strip()]
    else:
        labels = [str(x).strip() for x in raw if str(x).strip()]
    if not labels:
        raise ValueError("source label selection is empty")
    if len(set(labels)) != len(labels):
        duplicates = sorted({x for x in labels if labels.count(x) > 1})
        raise ValueError(f"duplicate --source-labels entries: {duplicates}")
    missing = [label for label in labels if label not in set(samples)]
    if missing:
        raise ValueError(f"source labels absent from samples: {missing}")
    return labels


def _read_species_map(path: Optional[Path], samples: Sequence[str]) -> Tuple[Dict[str, str], List[str]]:
    result = {sample: "UNMAPPED" for sample in samples}
    if path is not None:
        with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for row in reader:
                if not row or row[0].startswith("#"):
                    continue
                if row[0].strip().lower() in {"sample", "donor", "source_label"}:
                    continue
                if len(row) < 2:
                    raise ValueError(f"species map row needs donor and species: {row}")
                donor, species = row[0].strip(), row[1].strip()
                if donor not in result:
                    continue
                if not species:
                    raise ValueError(f"empty species label for donor {donor!r}")
                result[donor] = species
    missing = sorted(sample for sample, species in result.items() if not species or species == "UNMAPPED")
    return result, missing


def _read_receiver_classes(path: Path) -> List[ReceiverClass]:
    rows: List[ReceiverClass] = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if row[0].strip().lower() in {"fusion_class", "receiver_class", "identity"}:
                continue
            if len(row) < 3:
                raise ValueError("receiver-class rows require fusion_class, parent_1, parent_2")
            rows.append(ReceiverClass(row[0].strip(), row[1].strip(), row[2].strip()))
    if not rows:
        raise ValueError(f"no receiver classes found in {path}")
    names = [row.name for row in rows]
    if len(set(names)) != len(names):
        raise ValueError("receiver class names must be unique")
    return rows


def _all_pairs(samples: Sequence[str], include_homotypic: bool) -> List[ReceiverClass]:
    rows: List[ReceiverClass] = []
    for i, left in enumerate(samples):
        start = i if include_homotypic else i + 1
        for j in range(start, len(samples)):
            right = samples[j]
            rows.append(ReceiverClass(f"{left}+{right}", left, right))
    return rows


def _read_condf(
    path: Path,
    samples: Sequence[str],
    required_labels: Sequence[str],
) -> Tuple[List[Tuple[int, int]], np.ndarray, Dict[str, object]]:
    sample_index = {label: index for index, label in enumerate(samples)}
    required_indices = {sample_index[label]: label for label in required_labels}
    values: Dict[Tuple[int, int], Dict[int, float]] = defaultdict(dict)
    total_lines = 0
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            total_lines += 1
            fields = line.split()
            if len(fields) != 4:
                raise ValueError(f"{path}:{lineno}: expected four columns")
            idx1, genotype_type, idx2 = map(int, fields[:3])
            frac = float(fields[3])
            if idx2 < 0 or idx2 >= len(samples):
                continue
            if not math.isfinite(frac):
                raise ValueError(f"{path}:{lineno}: non-finite conditional fraction")
            values[(idx1, genotype_type)][idx2] = frac
    row_keys = sorted(values)
    if not row_keys:
        raise ValueError(f"no usable conditional fractions found in {path}")
    matrix = np.full((len(row_keys), len(samples)), np.nan, dtype=float)
    for row_index, key in enumerate(row_keys):
        for sample_idx, value in values[key].items():
            matrix[row_index, sample_idx] = value
    complete = np.ones(len(row_keys), dtype=bool)
    missing_by_source: Dict[str, int] = {}
    for index, label in required_indices.items():
        missing = ~np.isfinite(matrix[:, index])
        missing_by_source[label] = int(missing.sum())
        complete &= ~missing
    retained = int(complete.sum())
    if retained == 0:
        raise ValueError("no condf rows contain every required source and receiver-parent column")
    retained_keys = [row_keys[i] for i in np.where(complete)[0]]
    stats = {
        "noncomment_input_lines": int(total_lines),
        "unique_condf_rows": int(len(row_keys)),
        "retained_condf_rows": retained,
        "dropped_condf_rows": int(len(row_keys) - retained),
        "dropped_condf_fraction": float((len(row_keys) - retained) / len(row_keys)),
        "missing_rows_by_required_label": missing_by_source,
    }
    return retained_keys, matrix[complete], stats


def _read_weights(path: Optional[Path], row_keys: Sequence[Tuple[int, int]]) -> Tuple[np.ndarray, Dict[str, object]]:
    if path is None:
        return np.ones(len(row_keys), dtype=float), {"mode": "uniform", "path": ""}
    lines = [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]
    if not lines:
        raise ValueError(f"weights file is empty: {path}")
    split = [line.split() for line in lines]
    # Legacy one-value-per-retained-row mode.
    if all(len(row) == 1 for row in split):
        raw = [float(row[0]) for row in split]
        if len(raw) != len(row_keys):
            raise ValueError("legacy weights must contain one value per retained condf row")
        weights = np.asarray(raw, dtype=float)
        mode = "legacy_retained_row_order"
    else:
        keyed: Dict[Tuple[int, int], float] = {}
        for row in split:
            if row[0].lower() in {"idx1", "row_idx"}:
                continue
            if len(row) < 3:
                raise ValueError("keyed weights require idx1, genotype_type, weight")
            key = (int(row[0]), int(row[1]))
            if key in keyed:
                raise ValueError(f"duplicate keyed weight row: {key}")
            keyed[key] = float(row[2])
        missing = [key for key in row_keys if key not in keyed]
        if missing:
            raise ValueError(f"keyed weights missing {len(missing)} retained condf rows; examples={missing[:5]}")
        weights = np.asarray([keyed[key] for key in row_keys], dtype=float)
        mode = "keyed_idx1_genotype_type"
    if np.any(~np.isfinite(weights)) or np.any(weights < 0) or float(weights.sum()) <= 0:
        raise ValueError("weights must be finite, non-negative, and have positive sum")
    weights = weights / float(np.mean(weights))
    return weights, {"mode": mode, "path": str(Path(path).resolve()), "fingerprint": _fingerprint(path)}


def _weighted_norm(vector: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sqrt(np.sum(weights * vector * vector)))


def _residualize_affine(source: np.ndarray, parent_1: np.ndarray, parent_2: np.ndarray, weights: np.ndarray) -> np.ndarray:
    direction = parent_1 - parent_2
    centered = source - parent_2
    denom = float(np.sum(weights * direction * direction))
    if denom <= np.finfo(float).eps:
        return centered
    coefficient = float(np.sum(weights * centered * direction) / denom)
    return centered - coefficient * direction


def _effective_rank(matrix: np.ndarray, rank_tol: float) -> Tuple[int, np.ndarray, float]:
    if matrix.size == 0:
        return 0, np.array([], dtype=float), math.inf
    singular = np.linalg.svd(matrix, compute_uv=False)
    if singular.size == 0 or singular[0] <= 0:
        return 0, singular, math.inf
    threshold = float(rank_tol) * float(singular[0])
    rank = int(np.sum(singular > threshold))
    nonzero = singular[singular > threshold]
    condition = float(nonzero[0] / nonzero[-1]) if nonzero.size else math.inf
    return rank, singular, condition


def _pair_equivalent(
    left: Sequence[np.ndarray], right: Sequence[np.ndarray], weights: np.ndarray,
    cosine_threshold: float, norm_floor: float, distance_threshold: float,
) -> bool:
    stacked_left = np.concatenate(left)
    stacked_right = np.concatenate(right)
    stacked_weights = np.tile(weights, len(left))
    left_norm = _weighted_norm(stacked_left, stacked_weights)
    right_norm = _weighted_norm(stacked_right, stacked_weights)
    distance = _weighted_norm(stacked_left - stacked_right, stacked_weights)
    scale = max(left_norm, right_norm, norm_floor)
    if distance / scale > distance_threshold:
        return False
    if left_norm < norm_floor and right_norm < norm_floor:
        return True
    for lvec, rvec in zip(left, right):
        ln = _weighted_norm(lvec, weights)
        rn = _weighted_norm(rvec, weights)
        if ln < norm_floor and rn < norm_floor:
            continue
        if min(ln, rn) < norm_floor:
            return False
        cosine = float(np.sum(weights * lvec * rvec) / (ln * rn))
        if cosine < cosine_threshold:
            return False
    return True


def _complete_link_groups(labels: Sequence[str], equivalent: np.ndarray) -> Dict[str, str]:
    """Deterministic all-pairs groups; avoids transitive chaining of near matches."""
    groups: List[List[int]] = []
    for idx, _label in sorted(enumerate(labels), key=lambda item: item[1]):
        placed = False
        for group in groups:
            if all(bool(equivalent[idx, member]) for member in group):
                group.append(idx)
                placed = True
                break
        if not placed:
            groups.append([idx])
    assignments: Dict[str, str] = {}
    for number, group in enumerate(groups, 1):
        members = sorted(labels[i] for i in group)
        group_name = f"EQ{number:03d}" if len(members) > 1 else f"UNIQUE_{members[0]}"
        for member in members:
            assignments[member] = group_name
    return assignments


def _analyze_matrix(
    *, matrix: np.ndarray, samples: Sequence[str], source_labels: Sequence[str], receivers: Sequence[ReceiverClass],
    weights: np.ndarray, species_map: Mapping[str, str], cosine_threshold: float, norm_floor: float,
    distance_threshold: float, rank_tolerance: float,
) -> Tuple[List[Dict[str, object]], Dict[str, str], Dict[str, object]]:
    sample_index = {label: index for index, label in enumerate(samples)}
    for receiver in receivers:
        if receiver.parent_1 not in sample_index or receiver.parent_2 not in sample_index:
            raise ValueError(f"receiver parents absent from samples: {receiver}")
    residuals: Dict[str, Dict[str, np.ndarray]] = defaultdict(dict)
    source_rows: List[Dict[str, object]] = []
    for receiver in receivers:
        p1 = matrix[:, sample_index[receiver.parent_1]]
        p2 = matrix[:, sample_index[receiver.parent_2]]
        for source in source_labels:
            residual = _residualize_affine(matrix[:, sample_index[source]], p1, p2, weights)
            residuals[source][receiver.name] = residual
            norm = _weighted_norm(residual, weights)
            source_rows.append({
                "receiver_class": receiver.name,
                "parent_1": receiver.parent_1,
                "parent_2": receiver.parent_2,
                "source_label": source,
                "source_species": species_map.get(source, "UNMAPPED"),
                "residual_norm": f"{norm:.12g}",
                "informative": int(norm >= norm_floor),
            })
    equivalent = np.eye(len(source_labels), dtype=bool)
    for i, left in enumerate(source_labels):
        for j in range(i + 1, len(source_labels)):
            right = source_labels[j]
            eq = _pair_equivalent(
                [residuals[left][receiver.name] for receiver in receivers],
                [residuals[right][receiver.name] for receiver in receivers],
                weights, cosine_threshold, norm_floor, distance_threshold,
            )
            equivalent[i, j] = equivalent[j, i] = eq
    groups = _complete_link_groups(source_labels, equivalent)
    for row in source_rows:
        row["equivalence_group"] = groups[str(row["source_label"])]
    blocks = []
    for receiver in receivers:
        block = np.column_stack([residuals[source][receiver.name] for source in source_labels])
        blocks.append(np.sqrt(weights)[:, None] * block)
    stacked = np.vstack(blocks) if blocks else np.empty((0, len(source_labels)))
    rank, singular, condition = _effective_rank(stacked, rank_tolerance)
    summary = {
        "effective_rank": rank,
        "maximum_simplex_contrast_rank": max(0, len(source_labels) - 1),
        "singular_values": singular.tolist(),
        "condition_number": condition,
        "equivalence_group_count": len(set(groups.values())),
        "equivalence_map": groups,
    }
    return source_rows, groups, summary


def run_preflight(
    *, condf: Path, samples: Path, output_prefix: Path,
    receiver_classes: Optional[Path] = None, all_pairs: bool = False, include_homotypic: bool = False,
    species_map: Optional[Path] = None, source_labels: Optional[str | Sequence[str]] = None,
    weights: Optional[Path] = None, species_condf: Optional[Path] = None,
    species_samples: Optional[Path] = None, cosine_threshold: float = 0.999,
    norm_floor: float = 1e-8, distance_threshold: float = 1e-3,
    rank_tolerance: float = 1e-8, max_dropped_row_fraction: float = 0.25,
) -> Dict[str, object]:
    condf = Path(condf).resolve(); samples = Path(samples).resolve(); output_prefix = Path(output_prefix).resolve()
    if bool(receiver_classes) == bool(all_pairs):
        raise ValueError("exactly one of receiver_classes or all_pairs must be supplied")
    donor_samples = _read_samples(samples)
    selected_sources = _parse_source_labels(source_labels, donor_samples)
    receivers = _all_pairs(donor_samples, include_homotypic) if all_pairs else _read_receiver_classes(Path(receiver_classes))
    required = sorted(set(selected_sources) | {r.parent_1 for r in receivers} | {r.parent_2 for r in receivers})
    row_keys, full_matrix, condf_stats = _read_condf(condf, donor_samples, required)
    if condf_stats["dropped_condf_fraction"] > float(max_dropped_row_fraction):
        raise ValueError(
            f"condf completeness gate failed: dropped fraction {condf_stats['dropped_condf_fraction']:.6g} "
            f"> {max_dropped_row_fraction:.6g}"
        )
    donor_weights, weight_meta = _read_weights(weights, row_keys)
    donor_species_map, missing_species_all = _read_species_map(species_map, donor_samples)
    donor_rows, equivalence, donor_summary = _analyze_matrix(
        matrix=full_matrix, samples=donor_samples, source_labels=selected_sources,
        receivers=receivers, weights=donor_weights, species_map=donor_species_map,
        cosine_threshold=cosine_threshold, norm_floor=norm_floor,
        distance_threshold=distance_threshold, rank_tolerance=rank_tolerance,
    )

    # Prefer the exact species matrix used by CellBouncer. Fall back to donor
    # centroids only as a diagnostic and label that basis explicitly.
    species_summary: Dict[str, object]
    species_basis = "unavailable"
    species_labels: List[str] = []
    mapped_receiver_count = 0
    unmappable_species_receivers: List[str] = []
    if species_condf is not None and species_samples is not None and Path(species_condf).is_file() and Path(species_samples).is_file():
        sp_samples = _read_samples(Path(species_samples))
        mapped_receivers: List[ReceiverClass] = []
        for receiver in receivers:
            p1 = donor_species_map.get(receiver.parent_1, "UNMAPPED")
            p2 = donor_species_map.get(receiver.parent_2, "UNMAPPED")
            if p1 == "UNMAPPED" or p2 == "UNMAPPED" or p1 not in sp_samples or p2 not in sp_samples:
                unmappable_species_receivers.append(receiver.name)
                continue
            mapped_receivers.append(ReceiverClass(receiver.name, p1, p2))
        # Receiver names can repeat after species projection; deduplicate exact triples.
        mapped_receivers = list({(r.name, r.parent_1, r.parent_2): r for r in mapped_receivers}.values())
        mapped_receiver_count = len(mapped_receivers)
        species_labels = list(sp_samples)
        sp_required = sorted(set(species_labels) | {r.parent_1 for r in mapped_receivers} | {r.parent_2 for r in mapped_receivers})
        if mapped_receivers:
            sp_keys, sp_matrix, sp_stats = _read_condf(Path(species_condf), sp_samples, sp_required)
            sp_weights = np.ones(len(sp_keys), dtype=float)
            _sp_rows, _sp_groups, species_summary = _analyze_matrix(
                matrix=sp_matrix, samples=sp_samples, source_labels=species_labels,
                receivers=mapped_receivers, weights=sp_weights,
                species_map={x: x for x in sp_samples}, cosine_threshold=cosine_threshold,
                norm_floor=norm_floor, distance_threshold=distance_threshold,
                rank_tolerance=rank_tolerance,
            )
            species_summary["condf_completeness"] = sp_stats
            species_basis = "cellbouncer_species_condf"
        else:
            species_summary = {"effective_rank": 0, "maximum_simplex_contrast_rank": max(0, len(sp_samples)-1), "condition_number": math.inf, "singular_values": [], "equivalence_group_count": len(sp_samples), "equivalence_map": {x: f"UNIQUE_{x}" for x in sp_samples}}
            species_basis = "cellbouncer_species_condf_no_mappable_receiver_classes"
    else:
        species_labels = sorted({donor_species_map[s] for s in selected_sources if donor_species_map.get(s, "UNMAPPED") != "UNMAPPED"})
        species_summary = {"effective_rank": 0, "maximum_simplex_contrast_rank": max(0, len(species_labels)-1), "condition_number": math.inf, "singular_values": [], "equivalence_group_count": len(species_labels), "equivalence_map": {x: f"UNIQUE_{x}" for x in species_labels}}
        species_basis = "not_computed_without_species_condf"

    group_count = int(donor_summary["equivalence_group_count"])
    donor_rank = int(donor_summary["effective_rank"])
    species_rank = int(species_summary.get("effective_rank", 0))
    selected_missing_species = sorted(source for source in selected_sources if donor_species_map.get(source, "UNMAPPED") == "UNMAPPED")
    species_official_available = bool(
        not selected_missing_species
        and species_basis == "cellbouncer_species_condf"
        and mapped_receiver_count == len(receivers)
        and not unmappable_species_receivers
    )
    if donor_rank >= max(0, len(selected_sources) - 1) and group_count == len(selected_sources):
        recommended = "donor_resolved"
    elif donor_rank >= max(0, group_count - 1):
        recommended = "equivalence_group_resolved"
    elif species_official_available and species_rank >= max(0, len(species_labels) - 1):
        recommended = "species_resolved"
    else:
        recommended = "observable_profile_only"

    donor_classes: Dict[str, List[str]] = defaultdict(list)
    donor_partners: Dict[str, set[str]] = defaultdict(set)
    for receiver in receivers:
        donor_classes[receiver.parent_1].append(receiver.name)
        donor_classes[receiver.parent_2].append(receiver.name)
        donor_partners[receiver.parent_1].add(receiver.parent_2)
        donor_partners[receiver.parent_2].add(receiver.parent_1)
    fusion_rows = []
    for donor in sorted(donor_classes):
        classes = sorted(set(donor_classes[donor])); partners = sorted(donor_partners[donor])
        fusion_rows.append({
            "donor": donor, "fusion_class_count": len(classes), "fusion_classes": ",".join(classes),
            "partner_count": len(partners), "partners": ",".join(partners),
            "potential_cross_fusion_anchor": int(len(classes) >= 2 and len(partners) >= 2),
        })

    source_path = Path(str(output_prefix) + ".source_identifiability.tsv")
    eq_path = Path(str(output_prefix) + ".equivalence_groups.tsv")
    fusion_path = Path(str(output_prefix) + ".fusion_graph.tsv")
    order_path = Path(str(output_prefix) + ".source_order.tsv")
    summary_path = Path(str(output_prefix) + ".identifiability_summary.json")
    _write_tsv(source_path, ["receiver_class", "parent_1", "parent_2", "source_label", "source_species", "residual_norm", "informative", "equivalence_group"], donor_rows)
    _write_tsv(eq_path, ["source_label", "source_species", "equivalence_group"], ({"source_label": source, "source_species": donor_species_map.get(source, "UNMAPPED"), "equivalence_group": equivalence[source]} for source in selected_sources))
    _write_tsv(fusion_path, ["donor", "fusion_class_count", "fusion_classes", "partner_count", "partners", "potential_cross_fusion_anchor"], fusion_rows)
    _write_tsv(order_path, ["source_index", "source_label", "source_species", "samples_line_number"], ({"source_index": donor_samples.index(source), "source_label": source, "source_species": donor_species_map.get(source, "UNMAPPED"), "samples_line_number": donor_samples.index(source)+1} for source in selected_sources))

    inputs = {
        "condf": _fingerprint(condf), "samples": _fingerprint(samples),
        "receiver_classes": _fingerprint(Path(receiver_classes) if receiver_classes else None),
        "species_map": _fingerprint(species_map), "weights": _fingerprint(weights),
        "species_condf": _fingerprint(species_condf), "species_samples": _fingerprint(species_samples),
        "script": _fingerprint(Path(__file__).resolve()),
    }
    summary: Dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "receiver_class_count": len(receivers),
        "source_count": len(selected_sources),
        "species_count": len(species_labels),
        "source_labels": list(selected_sources),
        "species_labels": species_labels,
        "donor_effective_rank": donor_rank,
        "donor_maximum_simplex_contrast_rank": donor_summary["maximum_simplex_contrast_rank"],
        "donor_singular_values": donor_summary["singular_values"],
        "donor_condition_number": donor_summary["condition_number"],
        "species_effective_rank": species_rank,
        "species_maximum_simplex_contrast_rank": species_summary.get("maximum_simplex_contrast_rank", max(0, len(species_labels)-1)),
        "species_singular_values": species_summary.get("singular_values", []),
        "species_condition_number": species_summary.get("condition_number", math.inf),
        "species_design_basis": species_basis,
        "species_receiver_class_count": mapped_receiver_count,
        "species_unmappable_receiver_classes": sorted(set(unmappable_species_receivers)),
        "species_official_available": species_official_available,
        "unmapped_selected_sources": selected_missing_species,
        "unmapped_all_samples": missing_species_all,
        "equivalence_group_count": group_count,
        "recommended_resolution": recommended,
        "condf_completeness": condf_stats,
        "weights": weight_meta,
        "thresholds": {
            "cosine": cosine_threshold, "norm_floor": norm_floor,
            "distance": distance_threshold, "rank_tolerance": rank_tolerance,
            "max_dropped_row_fraction": max_dropped_row_fraction,
        },
        "potential_cross_fusion_anchor_donors": [row["donor"] for row in fusion_rows if row["potential_cross_fusion_anchor"]],
        "inputs": inputs,
        "outputs": {
            "source_identifiability": str(source_path), "equivalence_groups": str(eq_path),
            "fusion_graph": str(fusion_path), "source_order": str(order_path),
        },
        "output_fingerprints": {
            "source_identifiability": _fingerprint(source_path),
            "equivalence_groups": _fingerprint(eq_path),
            "fusion_graph": _fingerprint(fusion_path),
            "source_order": _fingerprint(order_path),
        },
    }
    signature_payload = _json_safe(dict(summary))
    summary = _json_safe(summary)
    summary["contract_signature"] = "sha256:" + hashlib.sha256(json.dumps(signature_payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")).hexdigest()
    _atomic_write_json(summary_path, summary)
    print(json.dumps(_json_safe(summary), indent=2, sort_keys=True, allow_nan=False))
    return summary


def validate_summary(path: Path) -> Tuple[bool, str, Dict[str, object]]:
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return False, "missing_summary", {}
    try:
        summary = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return False, f"invalid_json:{exc}", {}
    if summary.get("schema_version") != SCHEMA_VERSION:
        return False, "stale_schema", summary
    saved = str(summary.get("contract_signature", ""))
    unsigned = dict(summary); unsigned.pop("contract_signature", None)
    calc = "sha256:" + hashlib.sha256(json.dumps(_json_safe(unsigned), sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")).hexdigest()
    if not saved or saved != calc:
        return False, "contract_signature_mismatch", summary
    for name, record in dict(summary.get("inputs", {})).items():
        if not isinstance(record, Mapping) or not record.get("path"):
            continue
        p = Path(str(record["path"]))
        if not p.is_file() or _sha256(p) != str(record.get("sha256", "")):
            return False, f"stale_input:{name}", summary
    output_fingerprints = dict(summary.get("output_fingerprints", {}))
    for name, raw in dict(summary.get("outputs", {})).items():
        p = Path(str(raw))
        if not p.is_file() or p.stat().st_size == 0:
            return False, f"missing_output:{name}", summary
        recorded = output_fingerprints.get(name, {})
        if isinstance(recorded, Mapping) and recorded.get("sha256"):
            if _sha256(p) != str(recorded.get("sha256")):
                return False, f"modified_output:{name}", summary
    return True, "ok", summary


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--condf", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    receiver_group = parser.add_mutually_exclusive_group(required=True)
    receiver_group.add_argument("--receiver-classes", type=Path)
    receiver_group.add_argument("--all-pairs", action="store_true")
    parser.add_argument("--include-homotypic", action="store_true")
    parser.add_argument("--species-map", type=Path)
    parser.add_argument("--species-condf", type=Path)
    parser.add_argument("--species-samples", type=Path)
    parser.add_argument("--source-labels", help="Comma-separated subset of source labels")
    parser.add_argument("--weights", type=Path, help="Keyed idx1/genotype_type/weight table, or legacy one value per retained row")
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--cosine-threshold", type=float, default=0.999)
    parser.add_argument("--norm-floor", type=float, default=1e-8)
    parser.add_argument("--distance-threshold", type=float, default=1e-3)
    parser.add_argument("--rank-tolerance", type=float, default=1e-8)
    parser.add_argument("--max-dropped-row-fraction", type=float, default=0.25)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_parser().parse_args(argv)
    run_preflight(
        condf=args.condf, samples=args.samples, output_prefix=args.output_prefix,
        receiver_classes=args.receiver_classes, all_pairs=args.all_pairs,
        include_homotypic=args.include_homotypic, species_map=args.species_map,
        source_labels=args.source_labels, weights=args.weights,
        species_condf=args.species_condf, species_samples=args.species_samples,
        cosine_threshold=args.cosine_threshold, norm_floor=args.norm_floor,
        distance_threshold=args.distance_threshold, rank_tolerance=args.rank_tolerance,
        max_dropped_row_fraction=args.max_dropped_row_fraction,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
