#!/usr/bin/env python3
"""Workers for Tet 2025 ATAC clustering and RNA-ATAC embedding.

This file is invoked by SLURM scripts rendered by
``orchestrate_joint_atac_rna.py``.  Heavy imports are intentionally deferred so
that --help and unit tests work outside the cluster module environment.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import json
import math
import os
from pathlib import Path
import re
import sys
from typing import Any, Iterable, Sequence


VERSION = "1.0.3"
BARCODE_RE = re.compile(r"^[ACGT]{16}$")

CLASS_COLORS = {
    "AlloTet-Human-Chimp": "#FF9999",
    "AlloTet-Human-Orang": "#FFB366",
    "AlloTet-Human-Bonobo": "#FFCC99",
    "AlloTet-Chimp-Orang": "#99CCFF",
    "AlloTet-Chimp-Bonobo": "#CC99FF",
    "AlloTet-Bonobo-Orang": "#FF99CC",
    "AlloTet-Bonobo-Hy": "#FFCCFF",
    "AlloTet-Human-Hy": "#99FFCC",
    "AlloTet-Chimp-Hy": "#CCFFCC",
    "AlloTet-Orang-Hy": "#FFFFCC",
    "AutoTet-Human": "#66B2FF",
    "AutoTet-Chimp": "#FF6666",
    "AutoTet-Orang": "#FFB84D",
    "AutoTet-Bonobo": "#B366FF",
    "AutoTet-Hy": "#66FFB2",
    "Dip-Chimp": "#FF3333",
    "Dip-Human": "#3366FF",
    "Dip-Orang": "#FF8C00",
    "Dip-Bonobo": "#9933FF",
    "Dip-Hy": "#33FF99",
    "Unknown": "#CCCCCC",
}

SPECIES_NAMES = {
    "H": "Human",
    "C": "Chimp",
    "B": "Bonobo",
    "O": "Orang",
    "Hy": "Hy",
}


def log(message: str) -> None:
    print(message, flush=True)


def warn(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr, flush=True)


def canonical_barcode(value: object) -> str:
    raw = str(value).strip().upper()
    if not raw:
        return ""
    core = raw.split("-", 1)[0]
    return core if BARCODE_RE.fullmatch(core) else ""


def canonical_library(value: object) -> str:
    text = str(value).strip()
    patterns = (
        r"(?i)(?:^|[^a-z0-9])lib(?:rary)?[_ -]?(\d+)(?:$|[^0-9])",
        r"(?i)multiome-(?:atac|rna)_(\d+)(?:$|[^0-9])",
        r"^(\d+)$",
    )
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return f"lib{int(match.group(1))}"
    raise ValueError(f"invalid library label: {value}")


def natural_key(value: object) -> list[object]:
    return [
        int(piece) if piece.isdigit() else piece.lower()
        for piece in re.split(r"(\d+)", str(value))
    ]


def observation_frame(data_or_obs, source: str):
    """Materialize observations from pandas, AnnData, or backed SnapATAC2 data.

    SnapATAC2's backed ``.obs`` object is a PyDataFrameElem proxy rather than a
    pandas DataFrame. A full slice materializes that proxy as a Polars-like
    table, which can then be converted safely with ``to_pandas``.
    """
    import pandas as pd

    if isinstance(data_or_obs, pd.DataFrame):
        return data_or_obs.copy()

    table = getattr(data_or_obs, "obs", data_or_obs)
    names = None
    if hasattr(data_or_obs, "obs_names"):
        names = list(map(str, data_or_obs.obs_names))

    if isinstance(table, pd.DataFrame):
        frame = table.copy()
    else:
        try:
            materialized = table[:]
        except Exception as exc:
            raise TypeError(
                f"{source}: could not materialize the observation table: "
                f"{type(exc).__name__}: {exc}"
            ) from exc

        if isinstance(materialized, pd.DataFrame):
            frame = materialized.copy()
        else:
            to_pandas = getattr(materialized, "to_pandas", None)
            if not callable(to_pandas):
                raise TypeError(
                    f"{source}: unsupported materialized observation-table type "
                    f"{type(materialized).__name__}"
                )
            try:
                frame = to_pandas()
            except Exception as exc:
                raise TypeError(
                    f"{source}: could not convert the materialized observation "
                    f"table to pandas: {type(exc).__name__}: {exc}"
                ) from exc

    if names is not None:
        if len(frame) != len(names):
            raise ValueError(
                f"{source}: observation table has {len(frame)} rows but "
                f"obs_names has {len(names)} entries"
            )
        frame.index = names
    frame.index = frame.index.astype(str)
    return frame


def remove_null_h5ad_encodings(path: Path) -> list[str]:
    """Remove AnnData null-encoded metadata unsupported by anndata-rs.

    AnnData 0.12+ writes Python ``None`` values as HDF5 objects with
    ``encoding-type='null'``. SnapATAC2 2.10's Rust-backed AnnDataSet reader
    cannot open files containing those objects. Null objects carry no matrix or
    annotation values, so deleting them preserves the actual RNA data.
    """
    import h5py

    null_paths: list[str] = []
    with h5py.File(path, "r+") as handle:
        def collect(name: str, item) -> None:
            encoding = item.attrs.get("encoding-type")
            if hasattr(encoding, "item"):
                encoding = encoding.item()
            if isinstance(encoding, bytes):
                encoding = encoding.decode("utf-8", errors="replace")
            if encoding == "null":
                null_paths.append(name)

        handle.visititems(collect)
        for name in sorted(null_paths, key=lambda value: value.count("/"), reverse=True):
            if name in handle:
                del handle[name]
        if null_paths:
            handle.flush()
    return null_paths


def numeric_column_median(frame, column: str) -> float:
    """Return a numeric column median, or NaN when the column is unavailable."""
    import pandas as pd

    if column not in frame.columns:
        return math.nan
    values = pd.to_numeric(frame[column], errors="coerce")
    if values.notna().sum() == 0:
        return math.nan
    return float(values.median())


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", errors="replace")
    return path.open(mode.replace("t", ""), encoding="utf-8", errors="replace")


def unlink_if_exists(path: Path) -> None:
    if path.exists() and path.is_file():
        path.unlink()


def write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def write_frame_atomic(frame, path: Path, *, index: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    compression = "gzip" if str(path).endswith(".gz") else None
    frame.to_csv(
        temporary,
        sep="\t",
        index=index,
        na_rep="NA",
        compression=compression,
    )
    os.replace(temporary, path)


def read_task_table(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "library",
        "library_number",
        "fragments",
        "rna_matrix_dir",
        "library_output_dir",
    }
    if not rows:
        raise ValueError(f"task table is empty: {path}")
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"task table lacks columns: {', '.join(sorted(missing))}")
    rows.sort(key=lambda row: int(row["library_number"]))
    return rows


def load_chrom_sizes(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 2:
                fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected at least two columns")
            chrom = fields[0]
            try:
                size = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid chromosome size") from exc
            if size <= 0:
                raise ValueError(f"{path}:{line_number}: chromosome size must be positive")
            if chrom in result and result[chrom] != size:
                raise ValueError(f"{path}:{line_number}: conflicting size for {chrom}")
            result[chrom] = size
    if not result:
        raise ValueError(f"chromosome-size table has no records: {path}")
    return result


def read_rna_barcodes(path: Path) -> list[str]:
    barcodes: list[str] = []
    raw_by_core: dict[str, str] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            raw = line.strip()
            if not raw:
                continue
            core = canonical_barcode(raw)
            if not core:
                raise ValueError(f"{path}:{line_number}: invalid 16-base barcode: {raw}")
            previous = raw_by_core.setdefault(core, raw)
            if previous != raw:
                raise ValueError(
                    f"{path}:{line_number}: canonical barcode collision: {previous} and {raw}"
                )
            barcodes.append(core)
    if len(set(barcodes)) != len(barcodes):
        raise ValueError(f"duplicate canonical barcodes in {path}")
    return barcodes


def inspect_fragment_columns(path: Path, sample_records: int = 1000) -> dict[str, Any]:
    """Inspect enough records to distinguish native 5+ column and 4-column files.

    SnapATAC2 requires a fragment count in column five. Some upstream fragment
    writers emit only chromosome, start, end, and barcode. The source files are
    therefore sampled before import instead of assuming either representation.
    """
    column_counts: Counter[int] = Counter()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 4:
                raise ValueError(
                    f"{path}:{line_number}: fragment record has {len(fields)} columns; "
                    "expected at least chromosome, start, end, and barcode"
                )
            column_counts[len(fields)] += 1
            if sum(column_counts.values()) >= sample_records:
                break
    if not column_counts:
        raise ValueError(f"fragment file has no data records: {path}")
    return {
        "sampled_records": int(sum(column_counts.values())),
        "sampled_column_counts": {
            str(columns): int(count)
            for columns, count in sorted(column_counts.items())
        },
        "sample_contains_four_columns": 4 in column_counts,
    }


def prepare_fragment_input(
    source: Path,
    temp_dir: Path,
    library: str,
    whitelist_cores: set[str],
) -> tuple[Path, dict[str, Any], bool]:
    """Return a SnapATAC2-compatible fragment path without modifying source.

    Native files with five or more columns are returned directly. If the
    sample contains four-column records, a temporary gzip stream is written in
    which a count of one is appended. During that otherwise unavoidable copy,
    records are restricted to RNA-called barcodes to reduce temporary storage.
    """
    inspection = inspect_fragment_columns(source)
    if not inspection["sample_contains_four_columns"]:
        return (
            source,
            {
                **inspection,
                "fragment_input_mode": "native_5plus_columns",
                "temporary_copy": False,
            },
            False,
        )

    temp_dir.mkdir(parents=True, exist_ok=True)
    prepared = temp_dir / f"{library}.snapatac2_fragments.tsv.gz"
    partial = prepared.with_name(prepared.name + ".tmp")
    unlink_if_exists(partial)
    total_records = 0
    kept_records = 0
    dropped_nonwhitelist = 0
    converted_four_column = 0
    preserved_five_plus = 0
    try:
        with open_text(source) as source_handle, gzip.open(
            partial, "wt", encoding="utf-8", compresslevel=1
        ) as output_handle:
            for line_number, line in enumerate(source_handle, 1):
                if not line.strip():
                    continue
                if line.startswith("#"):
                    output_handle.write(line)
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 4:
                    raise ValueError(
                        f"{source}:{line_number}: fragment record has "
                        f"{len(fields)} columns; expected at least four"
                    )
                total_records += 1
                core = canonical_barcode(fields[3])
                if not core:
                    raise ValueError(
                        f"{source}:{line_number}: invalid 16-base fragment barcode: "
                        f"{fields[3]}"
                    )
                if core not in whitelist_cores:
                    dropped_nonwhitelist += 1
                    continue
                if len(fields) == 4:
                    fields.append("1")
                    converted_four_column += 1
                else:
                    preserved_five_plus += 1
                output_handle.write("\t".join(fields) + "\n")
                kept_records += 1
        if kept_records == 0:
            raise RuntimeError(
                f"{library}: no fragment records matched the RNA barcode whitelist"
            )
        os.replace(partial, prepared)
    except Exception:
        unlink_if_exists(partial)
        raise

    return (
        prepared,
        {
            **inspection,
            "fragment_input_mode": "temporary_four_to_five_column_normalization",
            "temporary_copy": True,
            "source_records_scanned": total_records,
            "records_kept": kept_records,
            "records_dropped_nonwhitelist": dropped_nonwhitelist,
            "four_column_records_converted": converted_four_column,
            "five_plus_column_records_preserved": preserved_five_plus,
        },
        True,
    )


def normalized_barcodes(values: Iterable[object], source: str) -> list[str]:
    raw_values = list(values)
    result = [canonical_barcode(value) for value in raw_values]
    invalid = [str(value) for value, core in zip(raw_values, result) if not core]
    if invalid:
        raise ValueError(f"{source}: invalid barcode examples: {invalid[:5]}")
    if len(set(result)) != len(result):
        counts = Counter(result)
        duplicates = [key for key, count in counts.items() if count > 1]
        raise ValueError(f"{source}: canonical barcode collisions: {duplicates[:5]}")
    return result


def library_paths(output_dir: Path) -> dict[str, Path]:
    return {
        "atac": output_dir / "atac_tiles.h5ad",
        "rna": output_dir / "rna_log1p.h5ad",
        "qc": output_dir / "cell_qc.tsv.gz",
        "fragment_sizes": output_dir / "fragment_size_distribution.tsv",
        "summary": output_dir / "library_summary.tsv",
        "complete": output_dir / "import.complete.json",
    }


def import_library(args: argparse.Namespace) -> int:
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import snapatac2 as snap

    library = canonical_library(args.library)
    fragments = args.fragments.resolve()
    rna_dir = args.rna_matrix_dir.resolve()
    output_dir = args.output_dir.resolve()
    temp_dir = args.temp_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)
    paths = library_paths(output_dir)
    if (
        not args.force
        and paths["complete"].is_file()
        and paths["atac"].is_file()
        and paths["rna"].is_file()
    ):
        log(f"{library}: completed outputs already exist; skipping")
        return 0

    required = [
        fragments,
        rna_dir / "barcodes.tsv.gz",
        rna_dir / "features.tsv.gz",
        rna_dir / "matrix.mtx.gz",
        args.chrom_sizes.resolve(),
    ]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("required inputs missing or empty: " + ", ".join(missing))

    chrom_sizes = load_chrom_sizes(args.chrom_sizes.resolve())
    rna_barcodes = read_rna_barcodes(rna_dir / "barcodes.tsv.gz")
    # Project fragments contain the extracted 16-mer.  The -1 alternatives make
    # the importer robust to a conventional 10X suffix without changing the
    # canonical downstream join key.
    whitelist = rna_barcodes + [barcode + "-1" for barcode in rna_barcodes]
    excluded = [item for item in args.exclude_chroms.split(",") if item]

    fragment_input, fragment_input_info, remove_fragment_input = prepare_fragment_input(
        fragments,
        temp_dir,
        library,
        set(rna_barcodes),
    )
    log(
        f"{library}: fragment input mode is "
        f"{fragment_input_info['fragment_input_mode']}"
    )

    atac_tmp = output_dir / "atac_tiles.tmp.h5ad"
    rna_tmp = output_dir / "rna_log1p.tmp.h5ad"
    unlink_if_exists(atac_tmp)
    unlink_if_exists(rna_tmp)

    log(f"{library}: importing {fragment_input}")
    atac = snap.pp.import_fragments(
        fragment_input,
        chrom_sizes=chrom_sizes,
        is_paired=True,
        file=atac_tmp,
        min_num_fragments=1,
        sorted_by_barcode=False,
        whitelist=whitelist,
        chrM=["chrM", "M", "MT"],
        tempdir=temp_dir,
        n_jobs=args.cpus,
    )
    if remove_fragment_input:
        unlink_if_exists(fragment_input)
    imported_barcodes = normalized_barcodes(atac.obs_names, f"{library} imported ATAC")

    tsse_available = False
    tsse_computed = False
    tsse_message = "not attempted"
    if args.gtf.is_file() and args.gtf.stat().st_size > 0:
        try:
            log(f"{library}: computing TSS enrichment from {args.gtf}")
            snap.metrics.tsse(atac, args.gtf.resolve(), n_jobs=args.cpus)
            tsse_computed = True
        except Exception as exc:  # optional QC metric; fragment-count analysis can proceed
            tsse_message = f"failed: {type(exc).__name__}: {exc}"
            warn(f"{library}: TSS enrichment {tsse_message}; continuing without a TSSe filter")
    else:
        tsse_message = f"annotation missing: {args.gtf}"
        warn(f"{library}: {tsse_message}; continuing without a TSSe filter")

    qc = observation_frame(atac, f"{library} imported ATAC observations")
    if tsse_computed:
        tsse_available = "tsse" in qc.columns
        tsse_message = "computed" if tsse_available else "metric returned without tsse"
        if not tsse_available:
            warn(f"{library}: TSS enrichment returned without a tsse column")
    qc.index = imported_barcodes
    qc.index.name = "barcode"
    qc.insert(0, "library", library)
    qc.insert(1, "cell_id", [f"{library}:{barcode}" for barcode in imported_barcodes])
    n_fragment = pd.to_numeric(qc.get("n_fragment"), errors="coerce")
    qc["pass_min_fragments"] = n_fragment >= args.min_fragments
    if tsse_available:
        tsse = pd.to_numeric(qc["tsse"], errors="coerce")
        qc["pass_min_tsse"] = tsse >= args.min_tsse
    else:
        qc["pass_min_tsse"] = True
    qc["pass_qc"] = qc["pass_min_fragments"] & qc["pass_min_tsse"]
    write_frame_atomic(qc.reset_index(), paths["qc"], index=False)

    log(
        f"{library}: filtering at n_fragment >= {args.min_fragments}"
        + (f" and TSSe >= {args.min_tsse}" if tsse_available else "")
    )
    snap.pp.filter_cells(
        atac,
        min_counts=args.min_fragments,
        min_tsse=args.min_tsse if tsse_available else None,
        inplace=True,
        n_jobs=args.cpus,
    )
    filtered_barcodes = normalized_barcodes(atac.obs_names, f"{library} filtered ATAC")
    if not filtered_barcodes:
        raise RuntimeError(f"{library}: no cells passed ATAC QC")
    cell_ids = [f"{library}:{barcode}" for barcode in filtered_barcodes]
    atac.obs_names = cell_ids
    atac.obs["library"] = [library] * len(cell_ids)
    atac.obs["barcode"] = filtered_barcodes
    atac.obs["cell_id"] = cell_ids

    fragment_size_status = "not computed"
    try:
        snap.metrics.frag_size_distr(
            atac,
            max_recorded_size=800,
            add_key="frag_size_distr",
            inplace=True,
            n_jobs=args.cpus,
        )
        distribution = np.asarray(atac.uns["frag_size_distr"], dtype=np.float64)
        fragment_sizes = pd.DataFrame(
            {
                "library": library,
                "fragment_size": np.arange(distribution.size, dtype=int),
                "fragment_count": distribution,
                "bin_note": [
                    ">800" if index == 0 else str(index)
                    for index in range(distribution.size)
                ],
            }
        )
        write_frame_atomic(fragment_sizes, paths["fragment_sizes"], index=False)
        fragment_size_status = "computed"
    except Exception as exc:
        fragment_size_status = f"failed: {type(exc).__name__}: {exc}"
        warn(f"{library}: fragment-size distribution {fragment_size_status}")

    log(f"{library}: generating {args.tile_size}-bp tile matrix")
    snap.pp.add_tile_matrix(
        atac,
        bin_size=args.tile_size,
        exclude_chroms=excluded,
        inplace=True,
        n_jobs=args.cpus,
    )
    if atac.n_vars == 0:
        raise RuntimeError(f"{library}: tile matrix has no features")

    log(f"{library}: loading matched RNA matrix")
    rna = sc.read_10x_mtx(
        rna_dir,
        var_names="gene_ids",
        make_unique=True,
        cache=False,
        gex_only=True,
    )
    rna_cores = normalized_barcodes(rna.obs_names, f"{library} RNA matrix")
    rna_lookup = {barcode: index for index, barcode in enumerate(rna_cores)}
    absent = [barcode for barcode in filtered_barcodes if barcode not in rna_lookup]
    if absent:
        raise RuntimeError(
            f"{library}: {len(absent)} ATAC cells are absent from the RNA matrix; "
            f"examples: {absent[:5]}"
        )
    rna = rna[[rna_lookup[barcode] for barcode in filtered_barcodes], :].copy()
    rna.obs_names = cell_ids
    rna.obs["library"] = library
    rna.obs["barcode"] = filtered_barcodes
    rna.obs["cell_id"] = cell_ids
    rna.layers["counts"] = rna.X.copy()
    sc.pp.calculate_qc_metrics(rna, percent_top=None, log1p=False, inplace=True)
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    log1p_metadata = rna.uns.get("log1p")
    if isinstance(log1p_metadata, dict) and log1p_metadata.get("base") is None:
        # Natural log is equivalent to base e, but a numeric value avoids the
        # AnnData 0.12 null encoding that anndata-rs cannot read.
        log1p_metadata["base"] = math.e
    rna.uns["normalization"] = "normalize_total(target_sum=1e4), log1p"
    rna.write_h5ad(rna_tmp, compression="gzip")
    removed_nulls = remove_null_h5ad_encodings(rna_tmp)
    if removed_nulls:
        log(
            f"{library}: removed unsupported null H5AD metadata: "
            + ", ".join("/" + value for value in removed_nulls)
        )

    # Close the backed ATAC file before its atomic rename.
    atac.close()
    os.replace(atac_tmp, paths["atac"])
    os.replace(rna_tmp, paths["rna"])

    passing = qc.loc[qc["pass_qc"]]
    summary = pd.DataFrame(
        [
            {
                "library": library,
                "rna_whitelist_cells": len(rna_barcodes),
                "atac_cells_with_fragments": len(imported_barcodes),
                "atac_cells_passing_qc": len(filtered_barcodes),
                "fraction_whitelist_passing_qc": len(filtered_barcodes) / len(rna_barcodes),
                "median_fragments_passing": numeric_column_median(
                    passing, "n_fragment"
                ),
                "median_tsse_passing": numeric_column_median(
                    passing, "tsse"
                )
                if tsse_available
                else math.nan,
                "tsse_status": tsse_message,
                "fragment_size_status": fragment_size_status,
                "fragment_input_mode": fragment_input_info["fragment_input_mode"],
                "fragment_input_details": json.dumps(
                    fragment_input_info, sort_keys=True, separators=(",", ":")
                ),
                "tile_size": args.tile_size,
                "atac_h5ad": str(paths["atac"]),
                "rna_h5ad": str(paths["rna"]),
            }
        ]
    )
    write_frame_atomic(summary, paths["summary"], index=False)
    write_json_atomic(
        paths["complete"],
        {
            "library": library,
            "status": "complete",
            "atac_cells_passing_qc": len(filtered_barcodes),
            "tsse_status": tsse_message,
            "fragment_input": fragment_input_info,
        },
    )
    log(f"{library}: complete with {len(filtered_barcodes):,} matched cells")
    return 0


def fallback_species(donor: str) -> str:
    donor = donor.strip()
    if donor.startswith("Chinobo"):
        return "Hy"
    if donor.startswith("Congo"):
        return "B"
    if donor.startswith("JOS"):
        return "O"
    if donor in {"KOLF", "H1", "H9"} or donor.startswith("H"):
        return "H"
    if donor.startswith("C"):
        return "C"
    return ""


def species_from_spsx(value: object) -> str:
    label = str(value).strip()
    if label.startswith("Hy"):
        return "Hy"
    return label[:1] if label[:1] in {"H", "C", "B", "O"} else ""


def split_assignment(value: object) -> list[str]:
    text = str(value).strip()
    if not text or text.lower() in {"nan", "na", "none", "unassigned"}:
        return []
    return [piece.strip() for piece in re.split(r"[+,]", text) if piece.strip()]


def classify_assignment(value: object, donor_species: dict[str, str]) -> tuple[str, str]:
    donors = split_assignment(value)
    if not donors:
        return "Unknown", "Unknown"
    species = [donor_species.get(donor, fallback_species(donor)) for donor in donors]
    if any(not item for item in species):
        return "Unknown", "Unknown"
    if len(donors) == 1:
        return f"Dip-{SPECIES_NAMES[species[0]]}", "Diploid"
    unique_species = set(species)
    if len(unique_species) == 1:
        return f"AutoTet-{SPECIES_NAMES[species[0]]}", "HomotypicTetraploid"
    order = ["H", "C", "B", "O", "Hy"]
    pair = [item for item in order if item in unique_species]
    pair_names = "-".join(SPECIES_NAMES[item] for item in pair)
    return f"AlloTet-{pair_names}", "HeterotypicTetraploid"


def read_xlsx_sheet_stdlib(path: Path, sheet_name: str):
    """Read a simple XLSX worksheet without requiring openpyxl."""
    import pandas as pd
    import xml.etree.ElementTree as et
    import zipfile

    main_ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    rel_ns = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    package_rel_ns = "http://schemas.openxmlformats.org/package/2006/relationships"

    def cell_text(cell, shared_strings: list[str]) -> str:
        cell_type = cell.get("t", "")
        if cell_type == "inlineStr":
            return "".join(
                node.text or "" for node in cell.findall(f".//{{{main_ns}}}t")
            )
        value = cell.find(f"{{{main_ns}}}v")
        text = "" if value is None or value.text is None else value.text
        if cell_type == "s" and text:
            return shared_strings[int(text)]
        if cell_type == "b":
            return "TRUE" if text == "1" else "FALSE"
        return text

    with zipfile.ZipFile(path) as archive:
        workbook = et.fromstring(archive.read("xl/workbook.xml"))
        relationship_id = None
        for sheet in workbook.findall(f".//{{{main_ns}}}sheet"):
            if sheet.get("name") == sheet_name:
                relationship_id = sheet.get(f"{{{rel_ns}}}id")
                break
        if not relationship_id:
            raise ValueError(f"worksheet is absent: {sheet_name}")

        relationships = et.fromstring(
            archive.read("xl/_rels/workbook.xml.rels")
        )
        target = None
        for relationship in relationships.findall(
            f"{{{package_rel_ns}}}Relationship"
        ):
            if relationship.get("Id") == relationship_id:
                target = relationship.get("Target")
                break
        if not target:
            raise ValueError(f"worksheet relationship is absent: {sheet_name}")
        target = target.lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target

        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared = et.fromstring(archive.read("xl/sharedStrings.xml"))
            shared_strings = [
                "".join(
                    node.text or ""
                    for node in item.findall(f".//{{{main_ns}}}t")
                )
                for item in shared.findall(f"{{{main_ns}}}si")
            ]

        worksheet = et.fromstring(archive.read(target))
        rows: list[dict[int, str]] = []
        for row in worksheet.findall(f".//{{{main_ns}}}row"):
            values: dict[int, str] = {}
            for cell in row.findall(f"{{{main_ns}}}c"):
                reference = cell.get("r", "")
                match = re.match(r"([A-Z]+)", reference)
                if not match:
                    continue
                index = 0
                for character in match.group(1):
                    index = index * 26 + ord(character) - ord("A") + 1
                values[index - 1] = cell_text(cell, shared_strings)
            if values:
                rows.append(values)
        if not rows:
            return pd.DataFrame()

        header_row = rows[0]
        width = max(max(row) for row in rows) + 1
        headers = [header_row.get(index, "").strip() for index in range(width)]
        records = [
            {
                headers[index]: row.get(index, "")
                for index in range(width)
                if headers[index]
            }
            for row in rows[1:]
        ]
        return pd.DataFrame.from_records(records)


def load_workbook_annotations(path: Path):
    import pandas as pd

    donor_species: dict[str, str] = {}
    library_metadata: dict[str, dict[str, str]] = {}
    if not path.is_file() or path.stat().st_size == 0:
        warn(f"workbook unavailable; class/batch annotations will be limited: {path}")
        return donor_species, library_metadata
    try:
        try:
            conversion = pd.read_excel(path, sheet_name="convert_extra", dtype=str)
            lanes = pd.read_excel(path, sheet_name="LaneQuickGuide", dtype=str)
        except ImportError:
            warn(
                "openpyxl is unavailable; reading workbook annotations with "
                "the built-in XLSX reader"
            )
            conversion = read_xlsx_sheet_stdlib(path, "convert_extra")
            lanes = read_xlsx_sheet_stdlib(path, "LaneQuickGuide")
        for _, row in conversion.iterrows():
            donor = str(row.get("VCF_ID", "")).strip()
            species = species_from_spsx(row.get("spsx", ""))
            if donor and donor.lower() != "nan" and species:
                donor_species[donor] = species

        for _, row in lanes.iterrows():
            raw_number = str(row.get("TotalLane", "")).strip()
            match = re.search(r"\d+", raw_number)
            if not match:
                continue
            library = f"lib{int(match.group())}"
            library_metadata[library] = {
                "diff_batch": str(row.get("DiffBatch", "Unknown")).strip(),
                "library_pool_role": str(
                    row.get("Treatment_Patterning", "Unknown")
                ).strip(),
            }
    except Exception as exc:
        warn(
            f"could not read optional workbook annotations from {path}: "
            f"{type(exc).__name__}: {exc}"
        )
    return donor_species, library_metadata


def canonical_cell_ids_from_obs(data_or_obs, source: str) -> list[str]:
    """Build stable library:16-mer IDs while preserving the current row order."""
    obs = observation_frame(data_or_obs, source)
    names = list(map(str, obs.index))
    if "library" in obs and "barcode" in obs:
        libraries = [canonical_library(value) for value in obs["library"]]
        barcodes = [canonical_barcode(value) for value in obs["barcode"]]
    else:
        libraries = [canonical_library(value.split(":", 1)[0]) for value in names]
        barcodes = [canonical_barcode(value.rsplit(":", 1)[-1]) for value in names]
    invalid = [names[index] for index, barcode in enumerate(barcodes) if not barcode]
    if invalid:
        raise ValueError(f"{source}: invalid cell barcode examples: {invalid[:5]}")
    cell_ids = [
        f"{library}:{barcode}" for library, barcode in zip(libraries, barcodes)
    ]
    if len(set(cell_ids)) != len(cell_ids):
        counts = Counter(cell_ids)
        duplicates = [key for key, count in counts.items() if count > 1]
        raise ValueError(f"{source}: duplicate canonical cell IDs: {duplicates[:5]}")
    return cell_ids


def annotated_obs(data_or_obs, identity_table: Path, workbook: Path):
    import pandas as pd

    obs = observation_frame(data_or_obs, "ATAC cohort observations")
    cell_ids = canonical_cell_ids_from_obs(obs, "ATAC cohort observations")
    obs["library"] = [value.split(":", 1)[0] for value in cell_ids]
    obs["barcode"] = [value.rsplit(":", 1)[-1] for value in cell_ids]
    obs["cell_id"] = cell_ids
    obs.index = obs["cell_id"]
    obs.index.name = "cell_id"

    identity_fields = [
        "assignment_status",
        "current_assignment",
        "proposed_assignment",
        "final_assignment",
        "assignment_change",
        "review_reason",
        "evidence_summary",
    ]
    joined = None
    if identity_table.is_file() and identity_table.stat().st_size > 0:
        try:
            identity = pd.read_csv(identity_table, sep="\t", dtype=str)
            required = {"library", "barcode", "final_assignment"}
            if not required.issubset(identity.columns):
                raise ValueError(
                    "missing required columns " + ", ".join(sorted(required - set(identity.columns)))
                )
            identity["library"] = identity["library"].map(canonical_library)
            identity["barcode"] = identity["barcode"].map(canonical_barcode)
            identity["cell_id"] = identity["library"] + ":" + identity["barcode"]
            if identity["cell_id"].duplicated().any():
                raise ValueError("duplicate library+barcode keys")
            available = [field for field in identity_fields if field in identity.columns]
            joined = identity.set_index("cell_id")[available].reindex(obs.index)
        except Exception as exc:
            warn(
                f"could not attach optional reconciled identity table {identity_table}: "
                f"{type(exc).__name__}: {exc}"
            )
    else:
        warn(f"optional reconciled identity table is unavailable: {identity_table}")
    for field in identity_fields:
        if joined is not None and field in joined:
            obs[field] = joined[field].fillna("Unassigned").astype(str)
        else:
            obs[field] = "Unassigned"

    donor_species, library_metadata = load_workbook_annotations(workbook)
    obs["diff_batch"] = [
        library_metadata.get(library, {}).get("diff_batch", "Unknown")
        for library in obs["library"]
    ]
    obs["library_pool_role"] = [
        library_metadata.get(library, {}).get("library_pool_role", "Unknown")
        for library in obs["library"]
    ]
    classifications = [
        classify_assignment(value, donor_species) for value in obs["final_assignment"]
    ]
    obs["identity_class"] = [item[0] for item in classifications]
    obs["biological_ploidy_class"] = [item[1] for item in classifications]
    return obs


def attach_obs(data, obs) -> None:
    data.obs_names = list(obs.index.astype(str))
    for column in obs.columns:
        values = obs[column]
        if str(values.dtype) == "category":
            values = values.astype(str)
        data.obs[column] = values.to_numpy()


def choose_batch_key(obs, requested: str) -> str | None:
    if requested == "none":
        return None
    if requested not in obs.columns:
        warn(f"batch column {requested!r} is unavailable; Harmony will be skipped")
        return None
    values = obs[requested].astype(str)
    informative = values[~values.isin(["", "Unknown", "NA", "nan"])]
    if informative.nunique() < 2:
        warn(f"batch column {requested!r} has fewer than two usable values")
        return None
    return requested


def plotting_modules():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns

    sns.set_theme(style="white", context="notebook")
    return np, pd, plt, sns


def plot_indices(n_cells: int, maximum: int, seed: int = 1729):
    import numpy as np

    if n_cells <= maximum:
        return np.arange(n_cells, dtype=int)
    return np.sort(np.random.default_rng(seed).choice(n_cells, size=maximum, replace=False))


def categorical_colors(categories: Sequence[str], palette: dict[str, str] | None = None):
    np, pd, plt, sns = plotting_modules()
    ordered = sorted(set(map(str, categories)), key=natural_key)
    if palette:
        colors = {category: palette.get(category, "#BDBDBD") for category in ordered}
    else:
        cmap = plt.get_cmap("turbo" if len(ordered) > 20 else "tab20")
        colors = {
            category: cmap(index / max(1, len(ordered) - 1))
            for index, category in enumerate(ordered)
        }
    return ordered, colors


def scatter_categorical(
    coordinates,
    labels,
    path: Path,
    title: str,
    maximum: int,
    *,
    palette: dict[str, str] | None = None,
    point_size: float = 1.0,
) -> None:
    np, pd, plt, sns = plotting_modules()
    coordinates = np.asarray(coordinates)
    labels = np.asarray(labels, dtype=str)
    selected = plot_indices(coordinates.shape[0], maximum)
    x = coordinates[selected]
    y = labels[selected]
    categories, colors = categorical_colors(y, palette)
    figure, axis = plt.subplots(figsize=(10.5, 8.5))
    for category in categories:
        mask = y == category
        axis.scatter(
            x[mask, 0],
            x[mask, 1],
            s=point_size,
            alpha=0.55,
            linewidths=0,
            rasterized=True,
            color=colors[category],
            label=category,
        )
    axis.set_title(title)
    axis.set_xlabel("UMAP 1")
    axis.set_ylabel("UMAP 2")
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_aspect("equal", adjustable="datalim")
    if len(categories) <= 50:
        columns = 2 if len(categories) > 16 else 1
        axis.legend(
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
            markerscale=5,
            fontsize=7 if len(categories) > 25 else 8,
            ncol=columns,
        )
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def scatter_continuous(
    coordinates,
    values,
    path: Path,
    title: str,
    maximum: int,
    *,
    label: str,
    log1p: bool = False,
) -> None:
    np, pd, plt, sns = plotting_modules()
    coordinates = np.asarray(coordinates)
    values = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy()
    selected = plot_indices(coordinates.shape[0], maximum)
    x = coordinates[selected]
    y = values[selected]
    if log1p:
        y = np.log1p(np.maximum(y, 0))
    finite = np.isfinite(y)
    figure, axis = plt.subplots(figsize=(9.5, 8.0))
    scatter = axis.scatter(
        x[finite, 0],
        x[finite, 1],
        c=y[finite],
        cmap="viridis",
        s=1.0,
        alpha=0.65,
        linewidths=0,
        rasterized=True,
    )
    axis.set_title(title)
    axis.set_xlabel("UMAP 1")
    axis.set_ylabel("UMAP 2")
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_aspect("equal", adjustable="datalim")
    colorbar = figure.colorbar(scatter, ax=axis, fraction=0.04, pad=0.02)
    colorbar.set_label(("log1p " if log1p else "") + label)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def compare_embeddings(
    left,
    right,
    labels,
    path: Path,
    left_title: str,
    right_title: str,
    maximum: int,
) -> None:
    np, pd, plt, sns = plotting_modules()
    left = np.asarray(left)
    right = np.asarray(right)
    labels = np.asarray(labels, dtype=str)
    selected = plot_indices(left.shape[0], maximum)
    categories, colors = categorical_colors(labels[selected])
    figure, axes = plt.subplots(1, 2, figsize=(16, 7.2))
    for axis, coords, title in zip(axes, (left, right), (left_title, right_title)):
        coords = coords[selected]
        sample_labels = labels[selected]
        for category in categories:
            mask = sample_labels == category
            axis.scatter(
                coords[mask, 0],
                coords[mask, 1],
                s=0.8,
                alpha=0.5,
                linewidths=0,
                rasterized=True,
                color=colors[category],
            )
        axis.set_title(title)
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
        axis.set_aspect("equal", adjustable="datalim")
    figure.suptitle("Raw versus batch-corrected embedding, colored by library")
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def composition_heatmap(
    obs,
    row_key: str,
    column_key: str,
    path: Path,
    title: str,
    *,
    normalize: str = "index",
) -> None:
    np, pd, plt, sns = plotting_modules()
    table = pd.crosstab(obs[row_key].astype(str), obs[column_key].astype(str))
    table = table.loc[
        sorted(table.index, key=natural_key), sorted(table.columns, key=natural_key)
    ]
    if normalize == "index":
        values = table.div(table.sum(axis=1).replace(0, np.nan), axis=0)
        colorbar_label = "Fraction within row"
    elif normalize == "columns":
        values = table.div(table.sum(axis=0).replace(0, np.nan), axis=1)
        colorbar_label = "Fraction within column"
    else:
        values = table
        colorbar_label = "Cells"
    width = min(24, max(8, 0.34 * len(values.columns) + 5))
    height = min(24, max(6, 0.32 * len(values.index) + 3))
    figure, axis = plt.subplots(figsize=(width, height))
    sns.heatmap(values, cmap="mako", ax=axis, cbar_kws={"label": colorbar_label})
    axis.set_title(title)
    axis.set_xlabel(column_key.replace("_", " ").title())
    axis.set_ylabel(row_key.replace("_", " ").title())
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def qc_by_library(obs, path: Path, maximum_per_library: int = 3000) -> None:
    np, pd, plt, sns = plotting_modules()
    frames = []
    rng = np.random.default_rng(1729)
    for library, frame in obs.groupby("library", observed=True):
        if len(frame) > maximum_per_library:
            frame = frame.iloc[
                np.sort(rng.choice(len(frame), size=maximum_per_library, replace=False))
            ]
        frames.append(frame)
    sampled = pd.concat(frames, axis=0)
    order = sorted(sampled["library"].astype(str).unique(), key=natural_key)
    figure, axes = plt.subplots(2, 1, figsize=(18, 10), sharex=True)
    fragments = pd.to_numeric(sampled.get("n_fragment"), errors="coerce")
    plot_frame = sampled.assign(log10_fragments=np.log10(np.maximum(fragments, 1)))
    sns.boxplot(
        data=plot_frame,
        x="library",
        y="log10_fragments",
        order=order,
        showfliers=False,
        color="#4C78A8",
        ax=axes[0],
    )
    axes[0].set_ylabel("log10 unique fragments")
    axes[0].set_xlabel("")
    if "tsse" in sampled:
        sampled = sampled.assign(tsse_numeric=pd.to_numeric(sampled["tsse"], errors="coerce"))
        sns.boxplot(
            data=sampled,
            x="library",
            y="tsse_numeric",
            order=order,
            showfliers=False,
            color="#F58518",
            ax=axes[1],
        )
        axes[1].set_ylabel("TSS enrichment")
    else:
        axes[1].text(0.5, 0.5, "TSS enrichment unavailable", ha="center", va="center")
        axes[1].set_yticks([])
    axes[1].set_xlabel("Library")
    for axis in axes:
        axis.tick_params(axis="x", rotation=90, labelsize=7)
    figure.suptitle("ATAC QC distributions by library")
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def fragment_size_plot(rows: Sequence[dict[str, str]], path: Path) -> None:
    np, pd, plt, sns = plotting_modules()
    frames = []
    for row in rows:
        source = library_paths(Path(row["library_output_dir"]))["fragment_sizes"]
        if source.is_file() and source.stat().st_size > 0:
            frame = pd.read_csv(source, sep="\t")
            frame = frame[pd.to_numeric(frame["fragment_size"], errors="coerce") > 0]
            total = pd.to_numeric(frame["fragment_count"], errors="coerce").sum()
            if total > 0:
                frame["fraction"] = pd.to_numeric(
                    frame["fragment_count"], errors="coerce"
                ) / total
                frames.append(frame)
    if not frames:
        raise RuntimeError("no per-library fragment-size distributions are available")
    combined = pd.concat(frames, ignore_index=True)
    pivot = combined.pivot_table(
        index="fragment_size", columns="library", values="fraction", fill_value=0
    ).sort_index()
    figure, axis = plt.subplots(figsize=(13, 6.5))
    for library in sorted(pivot.columns, key=natural_key):
        axis.plot(pivot.index, pivot[library], color="#7A7A7A", alpha=0.18, linewidth=0.7)
    axis.plot(
        pivot.index,
        pivot.median(axis=1),
        color="#D62728",
        linewidth=2.0,
        label="Across-library median",
    )
    axis.set_xlim(1, min(800, int(pivot.index.max())))
    axis.set_xlabel("Fragment length (bp)")
    axis.set_ylabel("Fraction of recorded fragments")
    axis.set_title("ATAC fragment-size periodicity across libraries")
    axis.legend(frameon=False)
    sns.despine(ax=axis)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def safe_plot(description: str, function, *args, **kwargs) -> None:
    try:
        function(*args, **kwargs)
    except Exception as exc:
        warn(f"plot {description!r} failed: {type(exc).__name__}: {exc}")


def clean_obs_for_h5ad(obs):
    import pandas as pd

    result = pd.DataFrame(obs).copy()
    result.index = result.index.astype(str)
    for column in result.columns:
        if result[column].dtype == object:
            result[column] = result[column].fillna("NA").astype(str)
    return result


def write_light_h5ad(path: Path, obs, embeddings: dict[str, Any], metadata: dict[str, Any]) -> None:
    import anndata as ad
    import scipy.sparse as sp

    obs = clean_obs_for_h5ad(obs)
    result = ad.AnnData(X=sp.csr_matrix((len(obs), 0), dtype="float32"), obs=obs)
    for key, value in embeddings.items():
        result.obsm[key] = value
    result.uns.update(metadata)
    temporary = path.with_name(path.stem + ".tmp.h5ad")
    unlink_if_exists(temporary)
    result.write_h5ad(temporary, compression="gzip")
    os.replace(temporary, path)


def export_cell_table(obs, embeddings: dict[str, Any], path: Path) -> None:
    import numpy as np
    import pandas as pd

    frame = clean_obs_for_h5ad(obs)
    for key, coordinates in embeddings.items():
        array = np.asarray(coordinates)
        if array.ndim != 2 or array.shape[1] < 2 or "umap" not in key.lower():
            continue
        label = key.removeprefix("X_")
        frame[f"{label}_1"] = array[:, 0]
        frame[f"{label}_2"] = array[:, 1]
    write_frame_atomic(frame.reset_index(), path, index=False)


def write_cluster_tables(obs, output_dir: Path, cluster_keys: Sequence[str]) -> None:
    import pandas as pd

    output_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for cluster_key in cluster_keys:
        if cluster_key not in obs:
            continue
        counts = (
            obs.groupby(cluster_key, observed=True)
            .size()
            .rename("n_cells")
            .reset_index()
        )
        counts.insert(0, "cluster_field", cluster_key)
        summaries.append(counts)
        for annotation in (
            "library",
            "diff_batch",
            "library_pool_role",
            "identity_class",
            "assignment_status",
        ):
            if annotation not in obs:
                continue
            table = pd.crosstab(
                obs[cluster_key].astype(str), obs[annotation].astype(str)
            )
            table.index.name = cluster_key
            write_frame_atomic(
                table.reset_index(),
                output_dir / f"{cluster_key}_by_{annotation}.tsv",
                index=False,
            )
    if summaries:
        write_frame_atomic(
            pd.concat(summaries, ignore_index=True),
            output_dir / "cluster_counts.tsv",
            index=False,
        )


def completed_component_paths(rows: Sequence[dict[str, str]], kind: str) -> list[tuple[str, Path]]:
    key = "atac" if kind == "atac" else "rna"
    result: list[tuple[str, Path]] = []
    failures: list[str] = []
    for row in rows:
        output_dir = Path(row["library_output_dir"])
        paths = library_paths(output_dir)
        path = paths[key]
        marker = paths["complete"]
        if not marker.is_file() or not path.is_file() or path.stat().st_size == 0:
            failures.append(f"{row['library']}: {path}")
        else:
            result.append((canonical_library(row["library"]), path.resolve()))
    if failures:
        raise FileNotFoundError(
            f"incomplete {kind.upper()} component outputs:\n  " + "\n  ".join(failures)
        )
    return result


def embedding_array(data, key: str):
    import numpy as np

    if key not in data.obsm:
        raise KeyError(f"embedding is absent: {key}")
    return np.asarray(data.obsm[key])


def run_harmony_or_raw(data, batch_key: str | None, source: str, target: str, cpus: int) -> str:
    import numpy as np
    import snapatac2 as snap

    if batch_key is None:
        data.obsm[target] = embedding_array(data, source)
        return "skipped"
    source_embedding = embedding_array(data, source)
    try:
        corrected = snap.pp.harmony(
            data,
            batch=batch_key,
            use_rep=source,
            inplace=False,
            n_jobs=cpus,
            max_iter_harmony=20,
        )
        corrected = np.asarray(corrected)
        # SnapATAC2 2.9/2.10 can return Harmony's native PCs-by-cells
        # orientation. Backed AnnDataSet.obsm requires cells-by-PCs.
        if corrected.shape == source_embedding.T.shape:
            corrected = corrected.T
        if corrected.shape != source_embedding.shape:
            raise RuntimeError(
                "Harmony returned shape "
                f"{corrected.shape}; expected {source_embedding.shape}"
            )
        data.obsm[target] = corrected
        return f"computed using {batch_key}"
    except Exception as exc:
        warn(
            f"Harmony failed for {source} using {batch_key}: "
            f"{type(exc).__name__}: {exc}; using the raw representation"
        )
        data.obsm[target] = embedding_array(data, source)
        return f"failed; raw fallback: {type(exc).__name__}: {exc}"


def cluster_from_representation(
    data,
    *,
    representation: str,
    cluster_key: str,
    umap_key: str,
    neighbors: int,
    resolution: float,
    extra_resolutions: bool = False,
) -> None:
    import snapatac2 as snap

    snap.pp.knn(
        data,
        n_neighbors=neighbors,
        use_rep=representation,
        method="hora",
        inplace=True,
        random_state=0,
    )
    snap.tl.leiden(
        data,
        resolution=resolution,
        key_added=cluster_key,
        random_state=0,
        inplace=True,
    )
    if extra_resolutions:
        for label, value in (("low", resolution / 2.0), ("high", resolution * 1.5)):
            snap.tl.leiden(
                data,
                resolution=value,
                key_added=f"{cluster_key}_{label}",
                random_state=0,
                inplace=True,
            )
    snap.tl.umap(
        data,
        use_rep=representation,
        key_added=umap_key.removeprefix("X_"),
        random_state=0,
        inplace=True,
        n_neighbors=neighbors,
    )


def cluster_atac(args: argparse.Namespace) -> int:
    import numpy as np
    import pandas as pd
    import snapatac2 as snap

    analysis_root = args.analysis_root.resolve()
    output_dir = analysis_root / "atac"
    table_dir = output_dir / "tables"
    figure_dir = args.figure_root.resolve() / "atac"
    output_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    dataset_path = output_dir / "atac_dataset.h5ads"
    light_path = output_dir / "atac_embedding.h5ad"
    complete = output_dir / "cluster.complete.json"
    if not args.force and complete.is_file() and dataset_path.is_file() and light_path.is_file():
        log("ATAC cohort clustering is already complete; skipping")
        return 0

    rows = read_task_table(args.task_table.resolve())
    components = completed_component_paths(rows, "atac")
    temporary = output_dir / "atac_dataset.tmp.h5ads"
    unlink_if_exists(temporary)
    if args.force or not complete.is_file():
        unlink_if_exists(dataset_path)

    log(f"Creating backed ATAC dataset from {len(components)} libraries")
    data = snap.AnnDataSet(
        adatas=[(key, path) for key, path in components],
        filename=temporary,
        add_key="sample",
        use_absolute_path=True,
    )
    try:
        names = list(map(str, data.obs_names))
        if len(set(names)) != len(names):
            raise RuntimeError("ATAC dataset cell IDs are not unique")
        obs = annotated_obs(data, args.identity_table.resolve(), args.workbook.resolve())
        if len(obs) != len(names):
            raise RuntimeError("ATAC annotations changed the number of observations")
        attach_obs(data, obs)

        log(f"Selecting {args.n_features:,} informative ATAC tiles")
        snap.pp.select_features(
            data,
            n_features=args.n_features,
            inplace=True,
            n_jobs=args.cpus,
            verbose=True,
        )
        selected = np.asarray(data.var["selected"], dtype=bool)
        if not selected.any():
            raise RuntimeError("ATAC feature selection returned zero tiles")

        log("Computing ATAC spectral embedding")
        snap.tl.spectral(
            data,
            n_comps=args.n_components,
            features="selected",
            random_state=0,
            distance_metric="cosine",
            inplace=True,
            num_threads=args.cpus,
        )
        cluster_from_representation(
            data,
            representation="X_spectral",
            cluster_key="leiden_atac_raw",
            umap_key="X_atac_raw_umap",
            neighbors=args.neighbors,
            resolution=args.leiden_resolution,
        )

        obs = observation_frame(data, "ATAC observations before Harmony")
        batch_key = choose_batch_key(obs, args.batch_key)
        harmony_status = run_harmony_or_raw(
            data,
            batch_key=batch_key,
            source="X_spectral",
            target="X_atac_harmony",
            cpus=args.cpus,
        )
        cluster_from_representation(
            data,
            representation="X_atac_harmony",
            cluster_key="leiden_atac",
            umap_key="X_atac_umap",
            neighbors=args.neighbors,
            resolution=args.leiden_resolution,
            extra_resolutions=True,
        )

        obs = observation_frame(data, "clustered ATAC observations")
        embeddings = {
            "X_atac_spectral": embedding_array(data, "X_spectral"),
            "X_atac_harmony": embedding_array(data, "X_atac_harmony"),
            "X_atac_raw_umap": embedding_array(data, "X_atac_raw_umap"),
            "X_atac_umap": embedding_array(data, "X_atac_umap"),
        }
        write_light_h5ad(
            light_path,
            obs,
            embeddings,
            {
                "workflow": "Tet 2025 ATAC clustering",
                "snapatac2_version": snap.__version__,
                "selected_atac_features": int(selected.sum()),
                "harmony_status": harmony_status,
                "batch_key": batch_key or "none",
            },
        )
        export_cell_table(obs, embeddings, table_dir / "atac_cells.tsv.gz")
        write_cluster_tables(
            obs,
            table_dir,
            ["leiden_atac_raw", "leiden_atac", "leiden_atac_low", "leiden_atac_high"],
        )

        summary_rows = []
        for library, frame in obs.groupby("library", observed=True):
            summary_rows.append(
                {
                    "library": library,
                    "n_cells": len(frame),
                    "median_n_fragment": numeric_column_median(
                        frame, "n_fragment"
                    ),
                    "median_tsse": numeric_column_median(frame, "tsse"),
                    "identity_annotated_fraction": (
                        frame["final_assignment"].astype(str) != "Unassigned"
                    ).mean(),
                }
            )
        summary = pd.DataFrame(summary_rows)
        summary["_library_number"] = summary["library"].str.extract(r"(\d+)").astype(int)
        summary = summary.sort_values("_library_number").drop(columns="_library_number")
        write_frame_atomic(summary, table_dir / "atac_library_summary.tsv", index=False)

        raw_umap = embeddings["X_atac_raw_umap"]
        corrected_umap = embeddings["X_atac_umap"]
        safe_plot(
            "ATAC clusters",
            scatter_categorical,
            corrected_umap,
            obs["leiden_atac"],
            figure_dir / "atac_umap_leiden.png",
            "ATAC clusters (Harmony representation)",
            args.max_plot_cells,
        )
        safe_plot(
            "raw ATAC clusters",
            scatter_categorical,
            raw_umap,
            obs["leiden_atac_raw"],
            figure_dir / "atac_raw_umap_leiden.png",
            "ATAC clusters (raw spectral representation)",
            args.max_plot_cells,
        )
        for field, title, palette in (
            ("library", "ATAC UMAP by library", None),
            ("diff_batch", "ATAC UMAP by differentiation batch", None),
            ("library_pool_role", "ATAC UMAP by library pool role", None),
            ("identity_class", "ATAC UMAP by reconciled biological class", CLASS_COLORS),
            ("assignment_status", "ATAC UMAP by identity-reconciliation status", None),
        ):
            safe_plot(
                field,
                scatter_categorical,
                corrected_umap,
                obs[field],
                figure_dir / f"atac_umap_{field}.png",
                title,
                args.max_plot_cells,
                palette=palette,
            )
        if "n_fragment" in obs:
            safe_plot(
                "fragments",
                scatter_continuous,
                corrected_umap,
                obs["n_fragment"],
                figure_dir / "atac_umap_fragments.png",
                "Unique ATAC fragments",
                args.max_plot_cells,
                label="unique fragments",
                log1p=True,
            )
        if "tsse" in obs:
            safe_plot(
                "TSS enrichment",
                scatter_continuous,
                corrected_umap,
                obs["tsse"],
                figure_dir / "atac_umap_tsse.png",
                "ATAC TSS enrichment",
                args.max_plot_cells,
                label="TSSe",
            )
        safe_plot(
            "raw versus corrected",
            compare_embeddings,
            raw_umap,
            corrected_umap,
            obs["library"],
            figure_dir / "atac_raw_vs_harmony_library.png",
            "Raw spectral UMAP",
            "Harmony UMAP",
            args.max_plot_cells,
        )
        safe_plot(
            "cluster by identity class",
            composition_heatmap,
            obs,
            "leiden_atac",
            "identity_class",
            figure_dir / "atac_cluster_by_identity_class.png",
            "ATAC cluster composition by biological class",
        )
        safe_plot(
            "cluster by library",
            composition_heatmap,
            obs,
            "leiden_atac",
            "library",
            figure_dir / "atac_cluster_by_library.png",
            "ATAC cluster composition by library",
        )
        safe_plot(
            "QC by library",
            qc_by_library,
            obs,
            figure_dir / "atac_qc_by_library.png",
        )
        safe_plot(
            "fragment-size periodicity",
            fragment_size_plot,
            rows,
            figure_dir / "atac_fragment_size_distribution.png",
        )
    finally:
        data.close()

    os.replace(temporary, dataset_path)
    write_json_atomic(
        complete,
        {
            "status": "complete",
            "n_libraries": len(components),
            "n_cells": int(len(obs)),
            "dataset": str(dataset_path),
            "embedding": str(light_path),
        },
    )
    log(f"ATAC clustering complete: {len(obs):,} cells")
    return 0


def select_rna_features(component_paths: Sequence[Path], n_features: int):
    """Select HVGs from log-normalized components with streaming moments."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    gene_names = None
    gene_symbols = None
    sums = None
    sums_of_squares = None
    n_cells = 0
    for path in component_paths:
        data = ad.read_h5ad(path, backed="r")
        try:
            current_names = np.asarray(data.var_names.astype(str))
            if gene_names is None:
                gene_names = current_names
                gene_symbols = (
                    data.var["gene_symbols"].astype(str).to_numpy()
                    if "gene_symbols" in data.var
                    else current_names.copy()
                )
                sums = np.zeros(len(gene_names), dtype=np.float64)
                sums_of_squares = np.zeros(len(gene_names), dtype=np.float64)
            elif not np.array_equal(gene_names, current_names):
                raise RuntimeError(f"RNA feature order differs in {path}")
            for item in data.chunked_X(2048):
                matrix = item[0] if isinstance(item, tuple) else item
                if sp.issparse(matrix):
                    sums += np.asarray(matrix.sum(axis=0)).ravel()
                    squared = matrix.copy()
                    squared.data = np.square(squared.data, dtype=np.float64)
                    sums_of_squares += np.asarray(squared.sum(axis=0)).ravel()
                else:
                    array = np.asarray(matrix, dtype=np.float64)
                    sums += array.sum(axis=0)
                    sums_of_squares += np.square(array).sum(axis=0)
                n_cells += matrix.shape[0]
        finally:
            data.file.close()

    if gene_names is None or sums is None or sums_of_squares is None or n_cells < 2:
        raise RuntimeError("RNA feature selection received no usable cells")
    means = sums / n_cells
    variances = np.maximum(
        (sums_of_squares - np.square(sums) / n_cells) / (n_cells - 1), 0
    )
    dispersion = np.log1p(variances / np.maximum(means, 1e-12))
    valid = np.isfinite(dispersion) & (means > 0)
    symbol_upper = np.char.upper(np.asarray(gene_symbols, dtype=str))
    valid &= ~np.char.startswith(symbol_upper, "MT-")

    score = np.full(len(gene_names), -np.inf, dtype=np.float64)
    valid_indices = np.flatnonzero(valid)
    if valid_indices.size:
        ranks = pd.Series(means[valid_indices]).rank(method="first")
        bin_count = min(20, max(1, valid_indices.size))
        bins = pd.qcut(ranks, q=bin_count, labels=False, duplicates="drop")
        valid_dispersion = dispersion[valid_indices]
        for bin_value in sorted(pd.unique(bins)):
            member = np.asarray(bins == bin_value)
            values = valid_dispersion[member]
            standard_deviation = values.std(ddof=1) if values.size > 1 else 0.0
            if standard_deviation > 0:
                standardized = (values - values.mean()) / standard_deviation
            else:
                standardized = values - values.mean()
            score[valid_indices[member]] = standardized

    count = min(n_features, int(np.isfinite(score).sum()))
    if count == 0:
        raise RuntimeError("RNA feature selection produced zero usable genes")
    selected_indices = np.argsort(score)[-count:]
    selected = np.zeros(len(gene_names), dtype=bool)
    selected[selected_indices] = True
    statistics = pd.DataFrame(
        {
            "gene_id": gene_names,
            "gene_symbol": gene_symbols,
            "mean_log_normalized_expression": means,
            "variance_log_normalized_expression": variances,
            "dispersion_score": dispersion,
            "standardized_dispersion": score,
            "selected": selected,
        }
    )
    return selected, statistics, n_cells


def concordance_table(obs, first: str, second: str):
    import pandas as pd

    table = pd.crosstab(obs[first].astype(str), obs[second].astype(str))
    table.index.name = first
    return table


def modality_comparison_figure(obs, embeddings: dict[str, Any], path: Path, maximum: int) -> None:
    np, pd, plt, sns = plotting_modules()
    selected = plot_indices(len(obs), maximum)
    labels = obs["identity_class"].astype(str).to_numpy()[selected]
    categories, colors = categorical_colors(labels, CLASS_COLORS)
    panels = [
        ("X_atac_umap", "ATAC"),
        ("X_rna_umap", "RNA"),
        ("X_joint_umap", "Joint"),
    ]
    figure, axes = plt.subplots(1, 3, figsize=(20, 6.4))
    for axis, (key, title) in zip(axes, panels):
        coords = np.asarray(embeddings[key])[selected]
        for category in categories:
            mask = labels == category
            axis.scatter(
                coords[mask, 0],
                coords[mask, 1],
                s=0.8,
                alpha=0.55,
                linewidths=0,
                rasterized=True,
                color=colors[category],
                label=category,
            )
        axis.set_title(f"{title} embedding")
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
        axis.set_aspect("equal", adjustable="datalim")
    handles, legend_labels = axes[-1].get_legend_handles_labels()
    if len(legend_labels) <= 30:
        figure.legend(
            handles,
            legend_labels,
            loc="center left",
            bbox_to_anchor=(0.995, 0.5),
            frameon=False,
            fontsize=8,
        )
    figure.suptitle("Matched-cell modality comparison, colored by reconciled biological class")
    figure.tight_layout(rect=(0, 0, 0.98, 0.96))
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def joint_embed(args: argparse.Namespace) -> int:
    import numpy as np
    import pandas as pd
    import snapatac2 as snap

    analysis_root = args.analysis_root.resolve()
    atac_path = analysis_root / "atac" / "atac_dataset.h5ads"
    output_dir = analysis_root / "joint"
    table_dir = output_dir / "tables"
    figure_dir = args.figure_root.resolve() / "joint"
    output_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    rna_path = output_dir / "rna_dataset.h5ads"
    light_path = output_dir / "joint_embedding.h5ad"
    complete = output_dir / "joint.complete.json"
    if not args.force and complete.is_file() and rna_path.is_file() and light_path.is_file():
        log("RNA-ATAC joint embedding is already complete; skipping")
        return 0
    if not atac_path.is_file() or atac_path.stat().st_size == 0:
        raise FileNotFoundError(f"ATAC dataset is missing: {atac_path}")

    rows = read_task_table(args.task_table.resolve())
    components = completed_component_paths(rows, "rna")
    for key, path in components:
        removed_nulls = remove_null_h5ad_encodings(path)
        if removed_nulls:
            log(
                f"{key}: repaired SnapATAC2-incompatible null H5AD metadata: "
                + ", ".join("/" + value for value in removed_nulls)
            )
    rna_temporary = output_dir / "rna_dataset.tmp.h5ads"
    unlink_if_exists(rna_temporary)
    unlink_if_exists(rna_path)

    atac = snap.read_dataset(atac_path)
    rna = snap.AnnDataSet(
        adatas=[(key, path) for key, path in components],
        filename=rna_temporary,
        add_key="sample",
        use_absolute_path=True,
    )
    obs = None
    try:
        atac_names = np.asarray(
            canonical_cell_ids_from_obs(atac, "ATAC joint-stage observations"),
            dtype=str,
        )
        rna_names = np.asarray(
            canonical_cell_ids_from_obs(rna, "RNA joint-stage observations"),
            dtype=str,
        )
        if atac_names.shape != rna_names.shape or not np.array_equal(atac_names, rna_names):
            mismatch = np.flatnonzero(atac_names != rna_names) if atac_names.shape == rna_names.shape else []
            example = int(mismatch[0]) if len(mismatch) else "shape"
            raise RuntimeError(
                "ATAC and RNA observations are not in identical order; "
                f"first mismatch={example}, ATAC={atac_names.shape}, RNA={rna_names.shape}"
            )

        atac.obs_names = list(atac_names)
        rna.obs_names = list(rna_names)

        atac_obs = observation_frame(atac, "ATAC joint-stage observations")
        for column in atac_obs.columns:
            values = atac_obs[column]
            if str(values.dtype) == "category":
                values = values.astype(str)
            rna.obs[column] = values.to_numpy()
        rna.obs_names = list(rna_names)

        log(f"Selecting {args.rna_features:,} RNA features using streaming moments")
        selected, hvg_statistics, n_hvg_cells = select_rna_features(
            [path for _, path in components], args.rna_features
        )
        if n_hvg_cells != len(rna_names):
            raise RuntimeError(
                f"RNA feature-statistic cell count {n_hvg_cells} != dataset count {len(rna_names)}"
            )
        rna.var["selected"] = selected
        write_frame_atomic(
            hvg_statistics,
            table_dir / "rna_feature_statistics.tsv.gz",
            index=False,
        )

        log("Computing RNA spectral embedding")
        snap.tl.spectral(
            rna,
            n_comps=args.n_components,
            features="selected",
            random_state=0,
            distance_metric="cosine",
            inplace=True,
            num_threads=args.cpus,
        )
        cluster_from_representation(
            rna,
            representation="X_spectral",
            cluster_key="leiden_rna_raw",
            umap_key="X_rna_raw_umap",
            neighbors=args.neighbors,
            resolution=args.leiden_resolution,
        )
        batch_key = choose_batch_key(atac_obs, args.batch_key)
        rna_harmony_status = run_harmony_or_raw(
            rna,
            batch_key=batch_key,
            source="X_spectral",
            target="X_rna_harmony",
            cpus=args.cpus,
        )
        cluster_from_representation(
            rna,
            representation="X_rna_harmony",
            cluster_key="leiden_rna",
            umap_key="X_rna_umap",
            neighbors=args.neighbors,
            resolution=args.leiden_resolution,
        )

        if "selected" not in atac.var:
            raise RuntimeError("ATAC dataset lacks selected features from the clustering stage")
        log("Computing matched-cell RNA-ATAC multi-spectral embedding")
        eigenvalues, joint_embedding = snap.tl.multi_spectral(
            [rna, atac],
            n_comps=args.n_components,
            features=["selected", "selected"],
            weights=None,
            random_state=0,
            weighted_by_sd=True,
        )
        atac.obsm["X_joint"] = joint_embedding
        cluster_from_representation(
            atac,
            representation="X_joint",
            cluster_key="leiden_joint_raw",
            umap_key="X_joint_raw_umap",
            neighbors=args.neighbors,
            resolution=args.leiden_resolution,
        )
        joint_harmony_status = run_harmony_or_raw(
            atac,
            batch_key=batch_key,
            source="X_joint",
            target="X_joint_harmony",
            cpus=args.cpus,
        )
        cluster_from_representation(
            atac,
            representation="X_joint_harmony",
            cluster_key="leiden_joint",
            umap_key="X_joint_umap",
            neighbors=args.neighbors,
            resolution=args.leiden_resolution,
            extra_resolutions=True,
        )

        # Keep RNA results next to the ATAC and joint coordinates in the
        # canonical cohort dataset and in a compact, portable h5ad.
        atac.obsm["X_rna_spectral"] = embedding_array(rna, "X_spectral")
        atac.obsm["X_rna_harmony"] = embedding_array(rna, "X_rna_harmony")
        atac.obsm["X_rna_raw_umap"] = embedding_array(rna, "X_rna_raw_umap")
        atac.obsm["X_rna_umap"] = embedding_array(rna, "X_rna_umap")
        rna_obs = observation_frame(rna, "clustered RNA observations")
        atac.obs["leiden_rna_raw"] = rna_obs["leiden_rna_raw"].astype(str).to_numpy()
        atac.obs["leiden_rna"] = rna_obs["leiden_rna"].astype(str).to_numpy()

        obs = observation_frame(atac, "clustered joint observations")
        embeddings = {
            "X_atac_spectral": embedding_array(atac, "X_spectral"),
            "X_atac_harmony": embedding_array(atac, "X_atac_harmony"),
            "X_atac_raw_umap": embedding_array(atac, "X_atac_raw_umap"),
            "X_atac_umap": embedding_array(atac, "X_atac_umap"),
            "X_rna_spectral": embedding_array(rna, "X_spectral"),
            "X_rna_harmony": embedding_array(rna, "X_rna_harmony"),
            "X_rna_raw_umap": embedding_array(rna, "X_rna_raw_umap"),
            "X_rna_umap": embedding_array(rna, "X_rna_umap"),
            "X_joint": embedding_array(atac, "X_joint"),
            "X_joint_harmony": embedding_array(atac, "X_joint_harmony"),
            "X_joint_raw_umap": embedding_array(atac, "X_joint_raw_umap"),
            "X_joint_umap": embedding_array(atac, "X_joint_umap"),
        }
        write_light_h5ad(
            light_path,
            obs,
            embeddings,
            {
                "workflow": "Tet 2025 matched RNA-ATAC embedding",
                "snapatac2_version": snap.__version__,
                "selected_rna_features": int(selected.sum()),
                "rna_harmony_status": rna_harmony_status,
                "joint_harmony_status": joint_harmony_status,
                "batch_key": batch_key or "none",
                "joint_spectral_eigenvalues": np.asarray(eigenvalues),
            },
        )
        export_cell_table(obs, embeddings, table_dir / "joint_cells.tsv.gz")
        write_cluster_tables(
            obs,
            table_dir,
            [
                "leiden_atac",
                "leiden_rna_raw",
                "leiden_rna",
                "leiden_joint_raw",
                "leiden_joint",
                "leiden_joint_low",
                "leiden_joint_high",
            ],
        )
        for first, second in (
            ("leiden_joint", "leiden_atac"),
            ("leiden_joint", "leiden_rna"),
            ("leiden_atac", "leiden_rna"),
        ):
            table = concordance_table(obs, first, second)
            write_frame_atomic(
                table.reset_index(),
                table_dir / f"{first}_by_{second}.tsv",
                index=False,
            )

        joint_umap = embeddings["X_joint_umap"]
        for field, title, palette in (
            ("leiden_joint", "Joint RNA-ATAC clusters", None),
            ("leiden_atac", "Joint UMAP colored by ATAC cluster", None),
            ("leiden_rna", "Joint UMAP colored by RNA cluster", None),
            ("library", "Joint UMAP by library", None),
            ("diff_batch", "Joint UMAP by differentiation batch", None),
            ("library_pool_role", "Joint UMAP by library pool role", None),
            ("identity_class", "Joint UMAP by reconciled biological class", CLASS_COLORS),
            ("assignment_status", "Joint UMAP by identity-reconciliation status", None),
        ):
            safe_plot(
                f"joint {field}",
                scatter_categorical,
                joint_umap,
                obs[field],
                figure_dir / f"joint_umap_{field}.png",
                title,
                args.max_plot_cells,
                palette=palette,
            )
        safe_plot(
            "modality comparison",
            modality_comparison_figure,
            obs,
            embeddings,
            figure_dir / "atac_rna_joint_comparison.png",
            args.max_plot_cells,
        )
        safe_plot(
            "joint raw versus corrected",
            compare_embeddings,
            embeddings["X_joint_raw_umap"],
            embeddings["X_joint_umap"],
            obs["library"],
            figure_dir / "joint_raw_vs_harmony_library.png",
            "Raw joint UMAP",
            "Harmony joint UMAP",
            args.max_plot_cells,
        )
        for second, filename, title in (
            ("leiden_atac", "joint_by_atac_cluster.png", "Joint versus ATAC clusters"),
            ("leiden_rna", "joint_by_rna_cluster.png", "Joint versus RNA clusters"),
            ("identity_class", "joint_by_identity_class.png", "Joint clusters by biological class"),
        ):
            safe_plot(
                f"joint composition {second}",
                composition_heatmap,
                obs,
                "leiden_joint",
                second,
                figure_dir / filename,
                title,
            )
    finally:
        rna.close()
        atac.close()

    os.replace(rna_temporary, rna_path)
    write_json_atomic(
        complete,
        {
            "status": "complete",
            "n_libraries": len(components),
            "n_cells": int(len(obs)),
            "rna_dataset": str(rna_path),
            "joint_embedding": str(light_path),
        },
    )
    log(f"Joint RNA-ATAC embedding complete: {len(obs):,} cells")
    return 0


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def common_cohort_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--task-table", type=Path, required=True)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--figure-root", type=Path, required=True)
    parser.add_argument("--batch-key", choices=("library", "diff_batch", "none"), default="library")
    parser.add_argument("--cpus", type=positive_int, default=8)
    parser.add_argument("--n-components", type=positive_int, default=30)
    parser.add_argument("--neighbors", type=positive_int, default=50)
    parser.add_argument("--leiden-resolution", type=float, default=1.0)
    parser.add_argument("--max-plot-cells", type=positive_int, default=200000)
    parser.add_argument("--force", action="store_true")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=VERSION)
    subparsers = parser.add_subparsers(dest="command", required=True)

    importer = subparsers.add_parser(
        "import-library",
        help="Import one ATAC fragment file, apply matched-cell QC, and prepare RNA.",
    )
    importer.add_argument("--library", required=True)
    importer.add_argument("--fragments", type=Path, required=True)
    importer.add_argument("--rna-matrix-dir", type=Path, required=True)
    importer.add_argument("--chrom-sizes", type=Path, required=True)
    importer.add_argument("--gtf", type=Path, required=True)
    importer.add_argument("--output-dir", type=Path, required=True)
    importer.add_argument("--temp-dir", type=Path, required=True)
    importer.add_argument("--cpus", type=positive_int, default=8)
    importer.add_argument("--min-fragments", type=positive_int, default=1000)
    importer.add_argument("--min-tsse", type=float, default=5.0)
    importer.add_argument("--tile-size", type=positive_int, default=5000)
    importer.add_argument(
        "--exclude-chroms",
        default="chrM,M,MT,chrY,Y,human_chimp_bonoborefChr1218",
    )
    importer.add_argument("--force", action="store_true")
    importer.set_defaults(function=import_library)

    atac = subparsers.add_parser(
        "cluster-atac", help="Create a backed cohort ATAC dataset and cluster it."
    )
    common_cohort_arguments(atac)
    atac.add_argument("--identity-table", type=Path, required=True)
    atac.add_argument("--workbook", type=Path, required=True)
    atac.add_argument("--n-features", type=positive_int, default=50000)
    atac.set_defaults(function=cluster_atac)

    joint = subparsers.add_parser(
        "joint-embed", help="Compute RNA-specific and matched RNA-ATAC embeddings."
    )
    common_cohort_arguments(joint)
    joint.add_argument("--rna-features", type=positive_int, default=3000)
    joint.set_defaults(function=joint_embed)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if hasattr(args, "min_tsse") and args.min_tsse < 0:
        raise ValueError("--min-tsse must be non-negative")
    if hasattr(args, "leiden_resolution") and args.leiden_resolution <= 0:
        raise ValueError("--leiden-resolution must be positive")
    return int(args.function(args))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, KeyError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        raise SystemExit(1)
