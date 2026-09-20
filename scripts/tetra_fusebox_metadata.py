#!/usr/bin/env python3
"""Build reconciled per-cell metadata and barcode lists for Fusebox ASE.

The input is the canonical finalized identity-reconciliation ledger.  Only
cells released downstream (``downstream_release_status == READY``) are
published.  The output metadata index matches the Tet2025 counts-H5AD index:

    <16-base barcode>-Tet_2025_Multiome-RNA_<library>

One raw 16-base barcode list is also written per library for ``count_ase -B``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Iterable, Iterator


RELEASE = "1.1.2"
MISSING = {"", ".", "NA", "N/A", "NONE", "NULL", "NAN"}
SPECIES_ORDER = {
    "Human": 0,
    "Chimp": 1,
    "Bonobo": 2,
    "Chinobo": 3,
    "Orangutan": 4,
}
KNOWN_SPECIES = frozenset(SPECIES_ORDER)


def clean(value: object) -> str:
    text = "" if value is None else str(value).strip()
    return "" if text.upper() in MISSING else text


def open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.endswith(".gz") else open(
        path, "r", encoding="utf-8", newline="")


def parse_libraries(values: Iterable[str]) -> list[int]:
    result: set[int] = set()
    for value in values:
        for token in re.split(r"[\s,]+", value.strip()):
            if not token:
                continue
            token = token.lower().removeprefix("lib")
            if "-" in token:
                left, right = token.split("-", 1)
                start, end = int(left), int(right)
                if start > end:
                    raise ValueError(f"descending library range: {token}")
                result.update(range(start, end + 1))
            else:
                result.add(int(token))
    if not result or min(result) < 1:
        raise ValueError("at least one positive library number is required")
    return sorted(result)


def library_number(value: object) -> int:
    text = clean(value).lower().removeprefix("lib")
    return int(text)


def canonical_species(value: object) -> str:
    raw = clean(value)
    key = re.sub(r"[^a-z0-9]+", "", raw.lower())
    aliases = {
        "h": "Human", "human": "Human", "homo": "Human",
        "c": "Chimp", "chimp": "Chimp", "chimpanzee": "Chimp",
        "pantroglodytes": "Chimp",
        "b": "Bonobo", "bonobo": "Bonobo", "panpaniscus": "Bonobo",
        "o": "Orangutan", "orang": "Orangutan", "orangutan": "Orangutan",
        "pongo": "Orangutan",
        "chinobo": "Chinobo", "chinobomcherry": "Chinobo",
        "chimpbonobo": "Chinobo", "hy": "Chinobo", "hybrid": "Chinobo",
    }
    return aliases.get(key, raw.replace(" ", "_"))


def normalized_header(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", clean(value).lower())


def infer_species_from_donor(value: object) -> str:
    """Resolve project donor IDs when panel metadata has no species column."""
    donor = clean(value)
    key = re.sub(r"[^a-z0-9]+", "", donor.lower())
    if not key:
        return ""
    if key.startswith("chinobo"):
        return "Chinobo"
    if key == "congoa4b":
        return "Bonobo"
    if key == "jos3c1":
        return "Orangutan"
    if key in {"kolf", "h1", "h9"} or re.fullmatch(r"h[0-9]+", key):
        return "Human"
    if re.fullmatch(r"c[0-9]+", key):
        return "Chimp"
    if re.fullmatch(r"b[0-9]+", key):
        return "Bonobo"
    if re.fullmatch(r"o[0-9]+", key):
        return "Orangutan"
    return ""


def load_species_map(path: str) -> dict[str, str]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if not fields:
            raise ValueError("panel metadata has no header")
        normalized_fields = {normalized_header(name): name for name in fields}
        species_field = next((normalized_fields[name] for name in (
            "species", "speciescode", "sp"
        ) if name in normalized_fields), None)
        if species_field:
            donor_field = next((normalized_fields[name] for name in (
                "individual", "individ", "vcfid", "donor", "sample", "line", "id"
            ) if name in normalized_fields), None)
        else:
            # CellBouncer's production panel-metadata reader has always used
            # the first two columns positionally: donor ID, then species.
            donor_field = fields[0]
        if donor_field is None:
            raise ValueError(
                "panel metadata must contain a donor column "
                "(individual/INDIVID/VCF_ID/donor/sample/line/id)")
        positional_species_field = (
            fields[1] if species_field is None and len(fields) >= 2 else None)
        result: dict[str, str] = {}
        unresolved: list[str] = []
        for row_number, row in enumerate(reader, start=2):
            donor = clean(row.get(donor_field))
            if not donor:
                continue
            species = ""
            if species_field:
                candidate = canonical_species(row.get(species_field))
                if candidate in KNOWN_SPECIES:
                    species = candidate
            elif positional_species_field:
                candidate = canonical_species(row.get(positional_species_field))
                if candidate in KNOWN_SPECIES:
                    species = candidate
            if not species:
                species = infer_species_from_donor(donor)
            if not species:
                unresolved.append(donor)
                continue
            previous = result.get(donor)
            if previous and previous != species:
                raise ValueError(
                    f"conflicting species for donor {donor!r} at row {row_number}: "
                    f"{previous!r} versus {species!r}")
            result[donor] = species
    if not result:
        detail = ",".join(sorted(set(unresolved))[:10])
        raise ValueError(
            "panel metadata yielded no donor-to-species mappings"
            + (f"; unresolved donor examples: {detail}" if detail else ""))
    return result


def assignment_components(value: object) -> list[str]:
    text = clean(value)
    if not text:
        return []
    return [part.strip() for part in text.replace(",", "+").split("+") if part.strip()]


def species_identity(assignment: str, donor_species: dict[str, str]) -> tuple[str, list[str]]:
    donors = assignment_components(assignment)
    missing = [donor for donor in donors if donor not in donor_species]
    species = {donor_species[donor] for donor in donors if donor in donor_species}
    ordered = sorted(species, key=lambda value: (SPECIES_ORDER.get(value, 99), value))
    return "_".join(ordered) if ordered else "Unknown", missing


def atomic_writer(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent,
        prefix=path.name + ".tmp.", delete=False)
    return handle


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {RELEASE}")
    parser.add_argument("--ledger", required=True)
    parser.add_argument("--panel-metadata", required=True)
    parser.add_argument("--libraries", nargs="+", required=True)
    parser.add_argument("--output-meta", required=True)
    parser.add_argument("--barcode-dir", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument(
        "--library-prefix", default="Tet_2025_Multiome-RNA_",
        help="H5AD obs-name library suffix prefix")
    args = parser.parse_args(argv)

    libraries = parse_libraries(args.libraries)
    requested = set(libraries)
    donor_species = load_species_map(args.panel_metadata)
    output_meta = Path(args.output_meta).resolve()
    barcode_dir = Path(args.barcode_dir).resolve()
    summary_path = Path(args.summary).resolve()
    for destination in (output_meta, summary_path):
        if destination.exists():
            raise FileExistsError(f"refusing to overwrite output: {destination}")
    for library in libraries:
        destination = barcode_dir / f"lib{library}.barcodes.tsv"
        if destination.exists():
            raise FileExistsError(f"refusing to overwrite output: {destination}")
    barcode_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    barcodes: dict[int, list[str]] = {library: [] for library in libraries}
    seen: set[tuple[int, str]] = set()
    unresolved: Counter[str] = Counter()
    status_counts: Counter[str] = Counter()

    with open_text(args.ledger) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        required = {
            "library", "barcode", "assignment_status", "final_assignment",
            "downstream_release_status",
        }
        if not required <= fields:
            raise ValueError(
                "identity ledger lacks required columns: "
                + ",".join(sorted(required - fields)))
        for row_number, row in enumerate(reader, start=2):
            try:
                library = library_number(row.get("library"))
            except (TypeError, ValueError):
                continue
            if library not in requested:
                continue
            if clean(row.get("downstream_release_status")).upper() != "READY":
                continue
            barcode = clean(row.get("barcode"))
            if not re.fullmatch(r"[ACGT]{16}", barcode):
                raise ValueError(
                    f"invalid selected barcode at ledger row {row_number}: {barcode!r}")
            key = (library, barcode)
            if key in seen:
                raise ValueError(f"duplicate selected cell in ledger: lib{library}/{barcode}")
            seen.add(key)
            assignment = clean(row.get("final_assignment"))
            if not assignment:
                continue
            donors = assignment_components(assignment)
            if not 1 <= len(donors) <= 2:
                raise ValueError(
                    "Fusebox requires one resolved donor or one resolved donor pair; "
                    f"ledger row {row_number} has {assignment!r}")
            species, missing = species_identity(assignment, donor_species)
            unresolved.update(missing)
            status = clean(row.get("assignment_status"))
            status_counts[status or "UNKNOWN"] += 1
            uid = clean(
                row.get("uid_or_uid_set") or row.get("reconciled_uid") or row.get("uid"))
            obs_name = f"{barcode}-{args.library_prefix}{library}"
            rows.append({
                "obs_name": obs_name,
                "library": str(library),
                "barcode": barcode,
                "composition": assignment,
                "donor1": donors[0],
                "donor2": donors[1] if len(donors) == 2 else "NA",
                "donor_count": str(len(donors)),
                "unique_donor_count": str(len(set(donors))),
                "species": species,
                "uid": uid or "NA",
                "assignment_status": status or "NA",
                "final_assignment": assignment,
                "production_assignment_source": clean(
                    row.get("production_assignment_source")) or "NA",
                "review_required": clean(row.get("review_required")) or "NA",
                "downstream_release_status": "READY",
                "event_id": clean(row.get("event_id")) or "NA",
                "current_ploidy_state": clean(row.get("current_ploidy_state")) or "NA",
                "nn_prob_tetraploid": clean(row.get("nn_prob_tetraploid")) or "NA",
                "uid_resolution_status": clean(row.get("uid_resolution_status")) or "NA",
            })
            barcodes[library].append(barcode)

    missing_libraries = [library for library in libraries if not barcodes[library]]
    if missing_libraries:
        raise ValueError(
            "no READY reconciled cells for "
            + ",".join(f"lib{library}" for library in missing_libraries))
    if not rows:
        raise ValueError("no READY reconciled cells selected")
    if unresolved:
        raise ValueError(
            "READY reconciled assignments contain donors absent from panel metadata: "
            + ",".join(
                f"{donor}:{count}" for donor, count in unresolved.most_common()))

    fieldnames = list(rows[0])
    meta_tmp = atomic_writer(output_meta)
    try:
        with meta_tmp:
            writer = csv.DictWriter(meta_tmp, fieldnames=fieldnames, delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(sorted(rows, key=lambda row: (int(row["library"]), row["barcode"])))
        os.replace(meta_tmp.name, output_meta)
    except Exception:
        try:
            os.unlink(meta_tmp.name)
        except OSError:
            pass
        raise

    for library in libraries:
        destination = barcode_dir / f"lib{library}.barcodes.tsv"
        tmp = atomic_writer(destination)
        try:
            with tmp:
                for barcode in sorted(barcodes[library]):
                    tmp.write(barcode + "\n")
            os.replace(tmp.name, destination)
        except Exception:
            try:
                os.unlink(tmp.name)
            except OSError:
                pass
            raise

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_payload = {
        "schema_version": "tetra_fusebox_metadata_v1",
        "release": RELEASE,
        "ledger": str(Path(args.ledger).resolve()),
        "panel_metadata": str(Path(args.panel_metadata).resolve()),
        "metadata": str(output_meta),
        "barcode_dir": str(barcode_dir),
        "libraries": libraries,
        "ready_cells": len(rows),
        "library_cells": {
            f"lib{library}": len(barcodes[library]) for library in libraries
        },
        "assignment_status_counts": dict(sorted(status_counts.items())),
        "unresolved_donor_counts": dict(sorted(unresolved.items())),
    }
    summary_tmp = atomic_writer(summary_path)
    try:
        with summary_tmp:
            json.dump(summary_payload, summary_tmp, indent=2, sort_keys=True)
            summary_tmp.write("\n")
        os.replace(summary_tmp.name, summary_path)
    except Exception:
        try:
            os.unlink(summary_tmp.name)
        except OSError:
            pass
        raise

    print(f"tetra_fusebox_metadata.py {RELEASE}")
    print(f"READY_CELLS\t{len(rows)}")
    for library in libraries:
        print(f"LIBRARY_CELLS\tlib{library}\t{len(barcodes[library])}")
    for status, count in sorted(status_counts.items()):
        print(f"ASSIGNMENT_STATUS\t{status}\t{count}")
    print(f"METADATA\t{output_meta}")
    print(f"BARCODE_DIR\t{barcode_dir}")
    print(f"SUMMARY\t{summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, csv.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
