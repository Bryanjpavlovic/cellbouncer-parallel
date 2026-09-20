#!/usr/bin/env python3
"""Shared, dependency-light helpers for the Tetraploid Arm CNV workflow."""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import re
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Sequence, Tuple


RELEASE = "2.5.0"
CELL_MANIFEST_SCHEMA = "tetra_arm_cell_manifest_v1"
AMBIENT_SCHEMA = "tetra_arm_ambient_sources_v1"
ASE_SCHEMA = "tetra_arm_ase_evidence_v2"
EXPRESSION_SCHEMA = "tetra_arm_expression_evidence_v1"
CALL_SCHEMA = "tetra_arm_cnv_calls_v2"


def clean(value) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.upper() in {"NA", "N/A", "NONE", "NULL", "."} else text


def finite_float(value, default=math.nan) -> float:
    text = clean(value)
    if not text:
        return default
    try:
        result = float(text)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def truthy(value) -> bool:
    return clean(value).upper() in {"1", "TRUE", "T", "YES", "Y", "PASS"}


def natural_key(value) -> Tuple:
    return tuple(
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", str(value))
    )


def canonical_barcode(value: str) -> str:
    """Normalize CellBouncer/10x display wrappers, never the DNA sequence itself."""
    text = clean(value)
    text = re.sub(r"^(?:lib)?\d+_", "", text, flags=re.I)
    return re.sub(r"(?:-\d+)+$", "", text)


def donor_components(identity: str) -> List[str]:
    """Return biological donor tokens without treating technical M{...} as ASE."""
    text = clean(identity)
    if not text or text.startswith("M{"):
        return []
    return [part.strip() for part in text.split("+") if part.strip()]


def canonical_pair(identity: str) -> Tuple[str, str]:
    components = donor_components(identity)
    # The finalized identity grammar for a heterotypic tetraploid is exactly
    # two distinct donor tokens.  Do not silently collapse multiplicity or a
    # malformed three-component identity with set().
    if len(components) != 2 or components[0] == components[1]:
        return "", ""
    unique = sorted(components, key=natural_key)
    return unique[0], unique[1]


def parse_libraries(values: Sequence[str], minimum: int = 1,
                    maximum: int = 40) -> List[int]:
    libraries = set()
    for raw in values:
        for token in str(raw).split(","):
            token = token.strip()
            if not token:
                continue
            match = re.fullmatch(r"(?:lib)?(\d+)(?:-(\d+))?", token, re.I)
            if not match:
                raise ValueError(f"invalid library selection: {token!r}")
            start = int(match.group(1))
            end = int(match.group(2) or start)
            if end < start:
                raise ValueError(f"descending library range is not allowed: {token}")
            for library in range(start, end + 1):
                if library < minimum or library > maximum:
                    raise ValueError(
                        f"library {library} is outside {minimum}-{maximum}")
                libraries.add(library)
    if not libraries:
        raise ValueError("at least one library must be selected")
    return sorted(libraries)


def open_text(path: str, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def read_tsv(path: str) -> Iterator[Dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"empty or headerless TSV: {path}")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"duplicate columns in TSV header: {path}")
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"malformed TSV row {path}:{line_number}")
            yield {str(k): str(v) for k, v in row.items()}


@contextmanager
def atomic_text(path: str, gzip_output: bool | None = None):
    destination = Path(path)
    if os.path.lexists(destination):
        raise FileExistsError(
            f"refusing to replace an existing output: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    suffix = ".gz" if (gzip_output if gzip_output is not None
                        else destination.suffix == ".gz") else ""
    fd, temporary = tempfile.mkstemp(
        prefix=destination.name + ".tmp.", suffix=suffix,
        dir=str(destination.parent))
    os.close(fd)
    try:
        if suffix:
            handle = gzip.open(temporary, "wt", encoding="utf-8", newline="")
        else:
            handle = open(temporary, "w", encoding="utf-8", newline="")
        with handle:
            yield handle
            handle.flush()
        if os.path.getsize(temporary) == 0:
            raise ValueError(f"refusing to publish empty output: {path}")
        # The temporary file lives beside its destination, so a hard link is
        # an atomic, no-clobber publish. Unlike os.replace(), this also closes
        # the race between the initial existence check and publication.
        os.link(temporary, destination)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def require_outputs_absent(paths: Iterable[os.PathLike | str]) -> None:
    """Fail before a multi-file producer publishes any member of its bundle."""
    conflicts = [
        os.path.abspath(os.fspath(path)) for path in paths
        if os.path.lexists(path)
    ]
    if conflicts:
        raise FileExistsError(
            "refusing to replace existing output(s): " + ", ".join(conflicts))


def write_tsv_atomic(path: str, rows: Iterable[Mapping[str, object]],
                     fields: Sequence[str]) -> int:
    count = 0
    with atomic_text(path) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t",
            lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
            count += 1
    return count


def write_json_atomic(path: str, payload: Mapping) -> None:
    with atomic_text(path, gzip_output=False) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def file_record(path: str) -> Dict[str, object]:
    return {"path": os.path.abspath(path)}


def require_file(path: str, label: str = "input") -> str:
    target = os.path.abspath(path)
    if not os.path.isfile(target) or os.path.getsize(target) <= 0:
        raise FileNotFoundError(f"{label} is missing or empty: {target}")
    return target


def parse_headerless_assignments(path: str) -> Dict[str, Tuple[str, str, float]]:
    assignments: Dict[str, Tuple[str, str, float]] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\r\n").split()
            if not fields:
                continue
            if fields[0].lower() == "barcode":
                continue
            if len(fields) != 4:
                raise ValueError(
                    f"expected four assignment columns at {path}:{line_number}")
            barcode = canonical_barcode(fields[0])
            if barcode in assignments:
                raise ValueError(f"duplicate assignment barcode {barcode}: {path}")
            score = finite_float(fields[3])
            assignments[barcode] = (fields[1], fields[2], score)
    if not assignments:
        raise ValueError(f"assignment file contains no records: {path}")
    return assignments


def parse_two_or_three_column(path: str) -> Dict[str, Tuple[float, float]]:
    result: Dict[str, Tuple[float, float]] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\r\n").split()
            if not fields:
                continue
            if len(fields) not in {2, 3}:
                raise ValueError(
                    f"expected two or three columns at {path}:{line_number}")
            key = canonical_barcode(fields[0])
            if key in result:
                raise ValueError(f"duplicate key {key}: {path}")
            value = finite_float(fields[1])
            error = finite_float(fields[2]) if len(fields) == 3 else math.nan
            if not math.isfinite(value):
                raise ValueError(f"non-finite value at {path}:{line_number}")
            result[key] = (value, error)
    return result


def ensure_unique_canonical(values: Iterable[str], label: str) -> Dict[str, str]:
    result: Dict[str, str] = {}
    for original in values:
        canonical = canonical_barcode(original)
        previous = result.get(canonical)
        if previous is not None and previous != original:
            raise ValueError(
                f"canonical barcode collision in {label}: {previous!r}, {original!r}")
        result[canonical] = original
    return result
