#!/usr/bin/env python3
"""Export the established Tet 2025 RNA analysis for matched-cell SCENIC+."""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path


DEFAULT_RNA_H5AD = Path(
    "/media/b/Rocket2_Plus/2025_Tet_testCluster/000_000_JupNotebooks/"
    "H5AD/02_3P_all40_processed_full_with_counts.h5ad"
)
DEFAULT_ANCHORED_CELLS = Path(
    "/media/b/Rocket2_Plus/2025_Tet_testCluster/000_000_JupNotebooks/"
    "joint_atac_rna_batch2_batch3_pools_v1/tables/rna_anchored_cells.tsv.gz"
)
DEFAULT_CLUSTER_INPUT = Path(
    "/home/b/Desktop/IndigoBeeGFS/tetraploid_multiome_cis_trans/ATAC/analysis/"
    "scenicplus_batch2_batch3_pools_v1/input"
)
DEFAULT_LIBRARIES = (15, 16, 23, 24, 31, 32, 39, 40)
BARCODE_RE = re.compile(r"^[ACGTN]{16}$", re.IGNORECASE)


def log(message: str) -> None:
    print(message, flush=True)


def require_file(path: Path, description: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"{description} is missing or empty: {path}")
    return path


def canonical_library(value: object) -> str:
    text = str(value).strip()
    for pattern in (
        r"(?i)(?:^|[^a-z0-9])lib(?:rary)?[_ -]?(\d+)(?:$|[^0-9])",
        r"(?i)multiome-(?:atac|rna)_(\d+)(?:$|[^0-9])",
        r"^(\d+)$",
    ):
        match = re.search(pattern, text)
        if match:
            return f"lib{int(match.group(1))}"
    raise ValueError(f"cannot parse library from {value!r}")


def canonical_barcode(value: object) -> str:
    barcode = str(value).strip().upper().split("-", 1)[0]
    if not BARCODE_RE.fullmatch(barcode):
        raise ValueError(f"invalid 16-bp barcode: {value!r}")
    return barcode


def truthy(series):
    import pandas as pd

    if pd.api.types.is_bool_dtype(series.dtype):
        return series.fillna(False).astype(bool)
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin({"1", "true", "t", "yes", "y"})
    )


def clean_label(value: object) -> str:
    label = re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip()).strip("_")
    return label or "unassigned"


def atomic_write_h5ad(adata, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + f".tmp.{os.getpid()}")
    try:
        adata.write_h5ad(temporary, compression="gzip")
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_write_table(frame, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + f".tmp.{os.getpid()}")
    try:
        frame.to_csv(temporary, sep="\t", index=False, compression="gzip")
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Subset the completed RNA H5AD to exact paired RNA/ATAC cells, retain "
            "the established RNA labels, and expose integer counts through .raw."
        )
    )
    parser.add_argument("--rna-h5ad", type=Path, default=DEFAULT_RNA_H5AD)
    parser.add_argument(
        "--anchored-cell-table", type=Path, default=DEFAULT_ANCHORED_CELLS
    )
    parser.add_argument(
        "--output-h5ad",
        type=Path,
        default=DEFAULT_CLUSTER_INPUT / "rna_for_scenicplus.h5ad",
    )
    parser.add_argument(
        "--output-cell-metadata",
        type=Path,
        default=DEFAULT_CLUSTER_INPUT / "cell_metadata.tsv.gz",
    )
    parser.add_argument(
        "--libraries",
        type=int,
        nargs="+",
        default=list(DEFAULT_LIBRARIES),
    )
    parser.add_argument("--force", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
        import scipy.sparse as sp

        rna_path = require_file(args.rna_h5ad, "processed RNA H5AD")
        cell_table_path = require_file(
            args.anchored_cell_table, "RNA-anchored cell table"
        )
        output_h5ad = args.output_h5ad.expanduser().resolve()
        output_metadata = args.output_cell_metadata.expanduser().resolve()
        if not args.force and (output_h5ad.exists() or output_metadata.exists()):
            raise FileExistsError(
                "one or more outputs already exist; use --force only when replacing "
                f"both exports: {output_h5ad}, {output_metadata}"
            )

        log(f"Reading matched-cell table: {cell_table_path}")
        cells = pd.read_csv(cell_table_path, sep="\t", low_memory=False)
        required = {
            "rna_cell_id",
            "rna_reference_matched",
            "rna_leiden_reference",
            "library",
            "barcode",
        }
        missing = required.difference(cells.columns)
        if missing:
            raise KeyError(
                "RNA-anchored cell table is missing: " + ", ".join(sorted(missing))
            )
        cells = cells.loc[truthy(cells["rna_reference_matched"])].copy()
        cells["sample_id"] = cells["library"].map(canonical_library)
        selected_libraries = {f"lib{int(number)}" for number in args.libraries}
        cells = cells.loc[cells["sample_id"].isin(selected_libraries)].copy()
        cells["barcode"] = cells["barcode"].map(canonical_barcode)
        cells["rna_cell_id"] = cells["rna_cell_id"].astype(str)
        cells["rna_leiden_reference"] = (
            cells["rna_leiden_reference"].astype(str).str.strip()
        )
        invalid_group = cells["rna_leiden_reference"].isin({"", "nan", "NA"})
        if invalid_group.any():
            raise ValueError(
                f"{int(invalid_group.sum()):,} matched cells lack an RNA Leiden label"
            )
        cells["cell_id"] = cells["barcode"] + "___" + cells["sample_id"]
        cells["peak_group"] = [
            f"{sample}__rna_{clean_label(cluster)}"
            for sample, cluster in zip(
                cells["sample_id"], cells["rna_leiden_reference"]
            )
        ]
        if cells["rna_cell_id"].duplicated().any():
            raise ValueError("matched-cell table has duplicate rna_cell_id values")
        if cells["cell_id"].duplicated().any():
            raise ValueError("constructed SCENIC+ cell identifiers are not unique")
        present = set(cells["sample_id"])
        absent_libraries = sorted(selected_libraries.difference(present))
        if absent_libraries:
            raise ValueError(
                "selected libraries have no exactly matched cells: "
                + ", ".join(absent_libraries)
            )

        log(f"Reading processed RNA reference in backed mode: {rna_path}")
        source = ad.read_h5ad(rna_path, backed="r")
        try:
            if not source.obs_names.is_unique:
                raise ValueError("processed RNA observation names are not unique")
            if not source.var_names.is_unique:
                raise ValueError(
                    "processed RNA gene names are not unique; they cannot be changed "
                    "without breaking the hg38 gene annotation join"
                )
            positions = source.obs_names.get_indexer(cells["rna_cell_id"])
            if (positions < 0).any():
                raise ValueError(
                    f"{int((positions < 0).sum()):,} matched RNA cells are absent "
                    "from the processed RNA H5AD"
                )
            subset = source[positions, :].to_memory()
        finally:
            source.file.close()

        if "counts" not in subset.layers:
            raise KeyError("processed RNA H5AD does not contain layers['counts']")
        counts = subset.layers["counts"].copy()
        values = counts.data if sp.issparse(counts) else np.asarray(counts).ravel()
        if values.size and (
            np.nanmin(values) < 0
            or not np.allclose(values, np.rint(values), rtol=0, atol=1e-6)
        ):
            raise ValueError("layers['counts'] is not a nonnegative integer-count matrix")

        # The current SCENIC+ CLI consumes adata.raw by default.  Keep the
        # processed matrix in X, but put the integer counts in raw explicitly.
        exported_obs = cells.set_index("cell_id", drop=False).copy()
        exported_obs.index = exported_obs.index.astype(str)
        exported_obs.index.name = None
        for column in exported_obs.columns:
            series = exported_obs[column]
            if isinstance(series.dtype, pd.CategoricalDtype):
                series = series.astype(object)
            if pd.api.types.is_object_dtype(series.dtype) or pd.api.types.is_string_dtype(
                series.dtype
            ):
                exported_obs[column] = (
                    series.astype(object).where(series.notna(), "NA").map(str)
                )
        subset.obs = exported_obs
        subset.obs_names = subset.obs.index.astype(str)
        subset.obs_names.name = None
        subset.var_names = subset.var_names.astype(str)
        raw = ad.AnnData(
            X=counts,
            obs=subset.obs.copy(),
            var=subset.var.copy(),
        )
        raw.obs_names.name = None
        subset.raw = raw
        subset.uns.clear()
        subset.obsm.clear()
        subset.varm.clear()
        subset.obsp.clear()

        metadata = cells.copy()
        for column in metadata.columns:
            if isinstance(metadata[column].dtype, pd.CategoricalDtype):
                metadata[column] = metadata[column].astype(str)
        metadata = metadata.replace({np.nan: "NA"})

        if args.force:
            output_h5ad.unlink(missing_ok=True)
            output_metadata.unlink(missing_ok=True)
        log(f"Writing SCENIC+ RNA H5AD: {output_h5ad}")
        atomic_write_h5ad(subset, output_h5ad)
        log(f"Writing exact paired-cell metadata: {output_metadata}")
        atomic_write_table(metadata, output_metadata)
        log(
            f"Exported {subset.n_obs:,} cells, {subset.n_vars:,} genes, "
            f"{metadata['rna_leiden_reference'].nunique():,} established RNA clusters"
        )
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
