#!/usr/bin/env python3
"""Build a newest-remap normalized H5AD for the existing ploidy NN.

The missing historical generate_h5ad script is not assumed.  Instead this
builder proves the exact transformation used by a retained reference H5AD from
its raw ``counts`` layer and normalized ``X``.  It then applies only that
validated cell-wise transform to the selected newest filtered MEX matrices.
No cells, expression values, or calls are copied from the reference H5AD.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import shlex
import stat
import subprocess
import sys
import tempfile
from typing import Sequence


RELEASE = "2.4.0"
LIBRARY_PREFIX = "Tet_2025_Multiome-RNA_"
DEFAULT_MAPPING_ROOT = (
    "/mnt/beegfs/tet2025_mapping_staging/"
    "rna3_all40_full_saturation_20260830_v2/rna3/mapping_output")
DEFAULT_REFERENCE_H5AD = (
    "/mnt/beegfs/tetmultiome_rna_mapped/mapping_output/h5a5_outs/"
    "unfiltered_normed_tetmultiome_rna.h5ad")


def parse_libraries(values: Sequence[str]) -> list[int]:
    result = set()
    for raw in values:
        for token in str(raw).split(","):
            token = token.strip().lower().removeprefix("lib")
            if not token:
                continue
            if "-" in token:
                first_text, last_text = token.split("-", 1)
                first, last = int(first_text), int(last_text)
            else:
                first = last = int(token)
            if first < 1 or last > 40 or last < first:
                raise ValueError(f"invalid library range: {raw}")
            result.update(range(first, last + 1))
    if not result:
        raise ValueError("at least one library is required")
    return sorted(result)


def require_file(path: str, label: str) -> str:
    path = os.path.abspath(path)
    if not os.path.isfile(path) or os.path.getsize(path) <= 0:
        raise ValueError(f"{label} is missing or empty: {path}")
    return path


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: str, include_sha256: bool = False) -> dict[str, object]:
    path = require_file(path, "input")
    info = os.stat(path)
    result: dict[str, object] = {
        "path": path,
        "realpath": os.path.realpath(path),
        "size": info.st_size,
        "mtime_ns": info.st_mtime_ns,
        "device": info.st_dev,
        "inode": info.st_ino,
    }
    if include_sha256:
        result["sha256"] = sha256_file(path)
    return result


def read_feature_rows(path: str) -> tuple[list[tuple[str, ...]], str]:
    """Read the full decompressed 10x feature identity and fingerprint it."""
    rows: list[tuple[str, ...]] = []
    digest = hashlib.sha256()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.endswith("\n"):
                raw += "\n"
            digest.update(raw.encode("utf-8"))
            fields = tuple(raw.rstrip("\r\n").split("\t"))
            if len(fields) < 2 or not fields[0] or not fields[1]:
                raise ValueError(
                    f"malformed 10x feature row {line_number}: {path}")
            rows.append(fields)
    if not rows:
        raise ValueError(f"10x feature table is empty: {path}")
    return rows, digest.hexdigest()


def atomic_text(path: str, payload: str) -> None:
    destination = os.path.abspath(path)
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=os.path.basename(destination) + ".tmp.",
        dir=os.path.dirname(destination), text=True)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        if os.path.lexists(destination):
            with open(destination, "r", encoding="utf-8") as handle:
                if handle.read() != payload:
                    raise ValueError(
                        f"refusing to replace differing file: {destination}; "
                        "select a new output path")
        else:
            os.link(temporary, destination)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def reference_library_number(value: object) -> int:
    text = str(value).strip()
    match = re.fullmatch(r"(?:lib)?(\d+)", text, flags=re.I)
    if not match:
        match = re.fullmatch(r"Tet_2025_Multiome-RNA_(\d+)", text)
    if not match:
        raise ValueError(f"invalid reference H5AD library label: {value!r}")
    library = int(match.group(1))
    if not 1 <= library <= 40:
        raise ValueError(f"reference H5AD library outside 1-40: {value!r}")
    return library


def h5_encoding_type(node: object) -> str:
    """Return a normalized AnnData HDF5 encoding type."""
    value = getattr(node, "attrs", {}).get("encoding-type", "")
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    return str(value).strip().lower()


def h5_scalar(value: object) -> object:
    """Decode a scalar read through h5py without stringifying byte literals."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    item = getattr(value, "item", None)
    if callable(item):
        value = item()
        if isinstance(value, bytes):
            return value.decode("utf-8")
    return value


def read_h5ad_obs_column(handle: object, column: str,
                         np: object) -> list[object]:
    """Read one H5AD obs column, including categorical encodings."""
    if "obs" not in handle or column not in handle["obs"]:
        raise ValueError(f"reference H5AD has no {column} obs column")
    node = handle["obs"][column]
    if hasattr(node, "shape"):
        values = np.asarray(node[...])
        if values.ndim != 1:
            raise ValueError(
                f"reference H5AD obs/{column} is not one-dimensional")
        return [h5_scalar(value) for value in values]

    encoding = h5_encoding_type(node)
    if encoding == "categorical":
        if "codes" not in node or "categories" not in node:
            raise ValueError(
                f"reference H5AD obs/{column} categorical encoding is incomplete")
        codes = np.asarray(node["codes"][...])
        categories = np.asarray(node["categories"][...])
        if codes.ndim != 1 or categories.ndim != 1:
            raise ValueError(
                f"reference H5AD obs/{column} categorical arrays are invalid")
        if codes.dtype.kind not in "iu":
            raise ValueError(
                f"reference H5AD obs/{column} categorical codes are not integers")
        if np.any(codes < 0) or np.any(codes >= len(categories)):
            raise ValueError(
                f"reference H5AD obs/{column} contains missing/invalid categories")
        decoded = [h5_scalar(value) for value in categories]
        return [decoded[int(code)] for code in codes]

    if encoding in {"nullable-integer", "nullable-boolean", "nullable-string"}:
        if "values" not in node or "mask" not in node:
            raise ValueError(
                f"reference H5AD obs/{column} nullable encoding is incomplete")
        values = np.asarray(node["values"][...])
        mask = np.asarray(node["mask"][...], dtype=bool)
        if values.ndim != 1 or mask.shape != values.shape or np.any(mask):
            raise ValueError(
                f"reference H5AD obs/{column} contains missing/invalid values")
        return [h5_scalar(value) for value in values]

    raise ValueError(
        f"unsupported reference H5AD obs/{column} encoding: {encoding or 'unknown'}")


def describe_h5ad_matrix(handle: object, location: str,
                         np: object) -> dict[str, object]:
    """Describe a dense or row-sliceable CSR H5AD matrix without loading it."""
    location = location.strip("/")
    if location not in handle:
        raise ValueError(f"reference H5AD has no {location} matrix")
    node = handle[location]
    encoding = h5_encoding_type(node)
    if hasattr(node, "shape"):
        shape = tuple(int(value) for value in node.shape)
        if len(shape) != 2:
            raise ValueError(f"reference H5AD {location} is not two-dimensional")
        if node.dtype.kind not in "biuf":
            raise ValueError(
                f"reference H5AD {location} has non-real-numeric dtype")
        if encoding in {"csc", "csc_matrix"}:
            raise ValueError(
                f"reference H5AD {location} uses CSC storage; rewrite it as "
                "CSR before bounded row-wise normalization validation")
        if encoding in {"csr", "csr_matrix"}:
            raise ValueError(
                f"reference H5AD {location} declares CSR but is stored as a dataset")
        return {
            "kind": "dense", "encoding": "DENSE", "node": node,
            "shape": shape,
        }


    if encoding in {"csc", "csc_matrix"}:
        raise ValueError(
            f"reference H5AD {location} uses backed CSC storage; rewrite it as "
            "CSR before bounded row-wise normalization validation")
    if encoding not in {"csr", "csr_matrix"}:
        raise ValueError(
            f"unsupported reference H5AD {location} encoding: "
            f"{encoding or 'unknown'}")
    if not all(name in node for name in ("data", "indices", "indptr")):
        raise ValueError(
            f"reference H5AD {location} CSR encoding is incomplete")
    raw_shape = node.attrs.get("shape")
    if raw_shape is None:
        raise ValueError(
            f"reference H5AD {location} CSR encoding has no shape")
    shape = tuple(int(value) for value in np.asarray(raw_shape).tolist())
    if len(shape) != 2 or shape[0] < 1 or shape[1] < 1:
        raise ValueError(
            f"reference H5AD {location} CSR shape is invalid: {shape}")
    data = node["data"]
    indices = node["indices"]
    indptr = node["indptr"]
    if (data.ndim != 1 or indices.ndim != 1 or indptr.ndim != 1 or
            len(data) != len(indices) or len(indptr) != shape[0] + 1):
        raise ValueError(
            f"reference H5AD {location} CSR arrays do not match its shape")
    if data.dtype.kind not in "biuf" or indices.dtype.kind not in "iu" or \
            indptr.dtype.kind not in "iu":
        raise ValueError(
            f"reference H5AD {location} CSR arrays have invalid dtypes")
    first_pointer = int(indptr[0])
    last_pointer = int(indptr[-1])
    if first_pointer != 0 or last_pointer != len(data):
        raise ValueError(
            f"reference H5AD {location} CSR pointers are invalid")
    return {
        "kind": "csr", "encoding": "CSR", "node": node,
        "shape": shape,
    }


def read_h5ad_matrix_rows(descriptor: dict[str, object], start: int,
                          stop: int, np: object,
                          sparse: object) -> object:
    """Read a bounded row interval from a direct HDF5 matrix descriptor."""
    rows, columns = descriptor["shape"]
    if start < 0 or stop < start or stop > rows:
        raise ValueError(f"invalid H5AD matrix row interval: {start}:{stop}")
    node = descriptor["node"]
    if descriptor["kind"] == "dense":
        return np.asarray(node[start:stop, :], dtype=np.float64)

    pointers = np.asarray(node["indptr"][start:stop + 1], dtype=np.int64)
    if len(pointers) != stop - start + 1 or np.any(np.diff(pointers) < 0):
        raise ValueError(
            f"invalid H5AD CSR pointers in rows {start}:{stop}")
    first = int(pointers[0])
    last = int(pointers[-1])
    if first < 0 or last < first or last > len(node["data"]):
        raise ValueError(
            f"out-of-range H5AD CSR pointers in rows {start}:{stop}")
    indices = np.asarray(node["indices"][first:last], dtype=np.int64)
    if np.any(indices < 0) or np.any(indices >= columns):
        raise ValueError(
            f"out-of-range H5AD CSR column index in rows {start}:{stop}")
    data = np.asarray(node["data"][first:last], dtype=np.float64)
    return sparse.csr_matrix(
        (data, indices, pointers - first), shape=(stop - start, columns))


def prove_reference_transform(path: str, libraries: Sequence[int],
                              chunk_cells: int,
                              tolerance: float) -> dict[str, object]:
    """Verify per-library log normalization across every reference cell."""
    import h5py
    import numpy as np
    import scipy.sparse

    with h5py.File(path, "r") as reference:
        counts_store = describe_h5ad_matrix(
            reference, "layers/counts", np)
        x_store = describe_h5ad_matrix(reference, "X", np)
        if counts_store["shape"] != x_store["shape"]:
            raise ValueError("reference counts and X shapes differ")
        total_cells, total_genes = counts_store["shape"]
        if total_cells < 1 or total_genes < 1:
            raise ValueError("reference H5AD is empty")
        library_values = read_h5ad_obs_column(reference, "library", np)
        if len(library_values) != total_cells:
            raise ValueError(
                "reference H5AD obs/library length differs from matrix rows")
        library_ids = np.asarray([
            reference_library_number(value)
            for value in library_values
        ], dtype=np.int16)

        targets = np.full(total_cells, np.nan, dtype=np.float64)
        maximum_error = 0.0
        for start in range(0, total_cells, chunk_cells):
            stop = min(start + chunk_cells, total_cells)
            raw_block = read_h5ad_matrix_rows(
                counts_store, start, stop, np, scipy.sparse)
            x_block = read_h5ad_matrix_rows(
                x_store, start, stop, np, scipy.sparse)
            if scipy.sparse.issparse(raw_block):
                raw_array = raw_block.tocsr().astype(np.float64)
                raw_array.sum_duplicates()
                raw_array.eliminate_zeros()
                raw_array.sort_indices()
                raw_values = raw_array.data
            else:
                raw_dense = np.asarray(raw_block, dtype=np.float64)
                raw_values = raw_dense
                raw_array = scipy.sparse.csr_matrix(raw_dense)
                raw_array.eliminate_zeros()
                raw_array.sort_indices()
            if (not np.isfinite(raw_values).all() or np.any(raw_values < 0) or
                    np.any(np.abs(raw_values - np.rint(raw_values)) > 1e-5)):
                raise ValueError(
                    "reference counts contain non-finite, negative, or "
                    "non-integer values")
            if scipy.sparse.issparse(x_block):
                x_array = x_block.tocsr().astype(np.float64)
                x_array.sum_duplicates()
                x_array.eliminate_zeros()
                x_array.sort_indices()
                x_values = x_array.data
            else:
                x_dense = np.asarray(x_block, dtype=np.float64)
                x_values = x_dense
                x_array = scipy.sparse.csr_matrix(x_dense)
                x_array.eliminate_zeros()
                x_array.sort_indices()
            if (not np.isfinite(x_values).all() or np.any(x_values < 0)):
                raise ValueError("reference X contains invalid normalized values")
            if (not np.array_equal(raw_array.indptr, x_array.indptr) or
                    not np.array_equal(raw_array.indices, x_array.indices)):
                raise ValueError(
                    "reference counts and X have different nonzero support in "
                    f"rows {start}:{stop}")
            row = np.repeat(
                np.arange(stop - start, dtype=np.int64),
                np.diff(raw_array.indptr))
            totals = np.bincount(
                row, weights=raw_array.data, minlength=stop - start)
            with np.errstate(over="ignore", invalid="ignore"):
                inferred = np.bincount(
                    row, weights=np.expm1(x_array.data),
                    minlength=stop - start)
            nonempty = totals > 0
            if np.any(~np.isfinite(inferred[nonempty])) or np.any(
                    inferred[nonempty] <= 0):
                raise ValueError(
                    f"reference normalization targets are invalid in rows {start}:{stop}")
            targets[start + np.flatnonzero(nonempty)] = inferred[nonempty]
            if len(row):
                expected = np.log1p(
                    raw_array.data * inferred[row] / totals[row])
                error = (float(np.max(np.abs(expected - x_array.data)))
                         if len(expected) else 0.0)
                if not math.isfinite(error):
                    raise ValueError(
                        "reference normalization proof produced non-finite error")
                maximum_error = max(maximum_error, error)

    target_by_library: dict[str, float] = {}
    checked_by_library: dict[str, int] = {}
    nonempty_by_library: dict[str, int] = {}
    spread_by_library: dict[str, float] = {}
    for library in sorted(set(int(value) for value in library_ids)):
        selected = library_ids == library
        values = targets[selected & np.isfinite(targets)]
        checked_by_library[str(library)] = int(np.count_nonzero(selected))
        nonempty_by_library[str(library)] = int(len(values))
        if len(values) < 2:
            if library in libraries:
                raise ValueError(
                    f"reference H5AD has too few nonempty lib{library} cells")
            continue
        center = float(np.median(values))
        spread = float(np.max(np.abs(values / center - 1.0)))
        spread_by_library[str(library)] = spread
        if (not math.isfinite(center) or center <= 0 or
                not math.isfinite(spread) or spread > tolerance):
            raise ValueError(
                "reference H5AD does not prove one normalize_total target "
                f"within lib{library}: spread={spread:.6g}, tolerance={tolerance:.6g}")
        rounded = float(round(center))
        target_by_library[str(library)] = (
            rounded if abs(center - rounded) / center <= tolerance else center)
    missing = [library for library in libraries
               if str(library) not in target_by_library]
    if missing:
        raise ValueError(
            "reference H5AD lacks validated normalization targets for: "
            + ",".join(f"lib{value}" for value in missing))
    if maximum_error > tolerance:
        raise ValueError(
            "reference H5AD does not match normalize_total + natural-log1p: "
            f"max_error={maximum_error:.6g}, tolerance={tolerance:.6g}")
    return {
        "transform": (
            "X=ln(1+counts*target_sum_by_library[library]/cell_total_counts)"),
        "normalization_scope": "per_library",
        "target_sum_by_library": target_by_library,
        "reference_scan_mode": "all_cells_bounded_chunks",
        "reference_total_cells": total_cells,
        "reference_cells_checked": total_cells,
        "reference_nonempty_cells": int(np.count_nonzero(np.isfinite(targets))),
        "reference_zero_count_cells": int(np.count_nonzero(~np.isfinite(targets))),
        "reference_cells_checked_by_library": checked_by_library,
        "reference_nonempty_cells_by_library": nonempty_by_library,
        "reference_chunk_cells": chunk_cells,
        "maximum_target_relative_spread_by_library": spread_by_library,
        "maximum_absolute_X_error": maximum_error,
        "counts_encoding": counts_store["encoding"],
        "X_encoding": x_store["encoding"],
        "tolerance": tolerance,
    }


def canonical_barcode(value: str) -> str:
    text = str(value).strip()
    while text.rsplit("-", 1)[-1].isdigit() and "-" in text:
        text = text.rsplit("-", 1)[0]
    return text


def build(args: argparse.Namespace) -> int:
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import scipy.sparse

    libraries = parse_libraries(args.libraries)
    mapping_root = os.path.abspath(args.mapping_input_root)
    reference_h5ad = require_file(
        args.normalization_reference_h5ad, "normalization reference H5AD")
    output = os.path.abspath(args.output)
    contract_path = output + ".contract.json"
    if os.path.lexists(output) or os.path.lexists(contract_path):
        raise ValueError(
            f"output already exists: {output}; choose a new versioned output path")
    if not os.path.isdir(mapping_root):
        raise ValueError(f"mapping input root is not a directory: {mapping_root}")

    normalization = prove_reference_transform(
        reference_h5ad, libraries, args.reference_chunk_cells,
        args.normalization_tolerance)
    target_sum_by_library = {
        int(key): float(value) for key, value in
        dict(normalization["target_sum_by_library"]).items()
    }
    matrices = []
    observations = []
    reference_var = None
    reference_feature_rows = None
    reference_feature_digest = None
    input_records = []
    cells_by_library: dict[str, int] = {}
    for library in libraries:
        mex = os.path.join(
            mapping_root, f"{LIBRARY_PREFIX}{library}", "filtered")
        barcodes = require_file(os.path.join(mex, "barcodes.tsv.gz"), "MEX barcodes")
        features = require_file(os.path.join(mex, "features.tsv.gz"), "MEX features")
        matrix = require_file(os.path.join(mex, "matrix.mtx.gz"), "MEX matrix")
        input_records.append({
            "library": library,
            "barcodes": file_record(barcodes, True),
            "features": file_record(features, True),
            "matrix": file_record(matrix, True),
        })
        feature_rows, feature_digest = read_feature_rows(features)
        if reference_feature_rows is None:
            reference_feature_rows = feature_rows
            reference_feature_digest = feature_digest
        elif feature_rows != reference_feature_rows:
            raise ValueError(
                f"lib{library} full 10x feature table differs from the first "
                "selected library")
        current = sc.read_10x_mtx(
            mex, var_names="gene_symbols", make_unique=True,
            cache=False, gex_only=True)
        current.X = scipy.sparse.csr_matrix(current.X, dtype=np.float32)
        if current.n_obs < 1 or current.n_vars < 1:
            raise ValueError(f"lib{library} filtered MEX is empty")
        if (not np.isfinite(current.X.data).all() or
                np.any(current.X.data < 0) or
                np.any(np.abs(current.X.data - np.rint(current.X.data)) > 1e-5)):
            raise ValueError(
                f"lib{library} MEX contains non-finite, negative, or "
                "non-integer count values")
        if np.any(np.asarray(current.X.sum(axis=1)).ravel() <= 0):
            raise ValueError(f"lib{library} filtered MEX contains a zero-count cell")
        current_barcodes = [canonical_barcode(value) for value in current.obs_names]
        if not all(current_barcodes) or len(current_barcodes) != len(set(current_barcodes)):
            raise ValueError(f"lib{library} MEX has empty/duplicate canonical barcodes")
        current.obs = pd.DataFrame({
            "library": library,
            "barcode": current_barcodes,
        }, index=[f"{barcode}-lib{library}" for barcode in current_barcodes])
        if reference_var is None:
            reference_var = current.var.copy()
        elif list(current.var_names) != list(reference_var.index):
            raise ValueError(
                f"lib{library} feature order differs from the first selected library")
        matrices.append(current.X)
        observations.append(current.obs)
        cells_by_library[str(library)] = current.n_obs

    if reference_var is None:
        raise ValueError("no MEX matrices were loaded")
    combined = ad.AnnData(
        X=scipy.sparse.vstack(matrices, format="csr", dtype=np.float32),
        obs=pd.concat(observations, axis=0),
        var=reference_var)
    if not combined.obs_names.is_unique or not combined.var_names.is_unique:
        raise ValueError("constructed H5AD names are not unique")
    combined.var["mt"] = combined.var_names.str.startswith("MT-")
    combined.var["ribo"] = combined.var_names.str.match(r"^(RPS|RPL)")
    combined.layers["counts"] = combined.X.copy()
    sc.pp.calculate_qc_metrics(
        combined, qc_vars=["mt", "ribo"], percent_top=None,
        log1p=False, inplace=True)
    raw_totals = np.asarray(combined.X.sum(axis=1)).ravel().astype(np.float64)
    row_targets = np.asarray([
        target_sum_by_library[int(value)]
        for value in combined.obs["library"]
    ], dtype=np.float64)
    normalized = combined.X.astype(np.float64).multiply(
        (row_targets / raw_totals)[:, None]).tocsr()
    scaled_totals = np.asarray(normalized.sum(axis=1)).ravel()
    if not np.allclose(
            scaled_totals, row_targets, rtol=args.normalization_tolerance,
            atol=args.normalization_tolerance):
        raise ValueError("failed to apply validated per-library normalization targets")
    np.log1p(normalized.data, out=normalized.data)
    normalized.data = normalized.data.astype(np.float32)
    combined.X = normalized
    combined.uns["log1p"] = {"base": None}
    if (not np.isfinite(combined.X.data).all() or
            np.any(combined.X.data < 0)):
        raise ValueError("normalized H5AD contains invalid expression values")
    combined.uns["tetra_arm_ploidy_h5ad"] = {
        "schema_version": "tetra_arm_ploidy_h5ad_v3",
        "mapping_input_root": mapping_root,
        "libraries": libraries,
        "normalization_reference_h5ad": reference_h5ad,
        "normalization": normalization,
    }

    os.makedirs(os.path.dirname(output), exist_ok=True)
    temporary = output + f".tmp.{os.getpid()}"
    try:
        combined.write_h5ad(temporary, compression="gzip")
        with open(temporary, "rb") as handle:
            os.fsync(handle.fileno())
        os.link(temporary, output)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass

    output_record = file_record(output, True)
    contract = {
        "schema_version": "tetra_arm_ploidy_h5ad_v3",
        "release": RELEASE,
        "status": "PASS",
        "mapping_input_root": mapping_root,
        "libraries": libraries,
        "cells": combined.n_obs,
        "genes": combined.n_vars,
        "cells_by_library": cells_by_library,
        "feature_name_source": "features.tsv.gz column 2 gene symbols",
        "full_feature_table_sha256": reference_feature_digest,
        "expression_source": "newest filtered MEX only",
        "normalization_reference": file_record(reference_h5ad, False),
        "normalization": normalization,
        "inputs": input_records,
        "output": output_record,
    }
    atomic_text(contract_path, json.dumps(contract, indent=2, sort_keys=True) + "\n")
    print(f"PASS: {output}")
    print(f"Contract: {contract_path}")
    print(f"Cells: {combined.n_obs}; genes: {combined.n_vars}")
    print("Validated target sums: " + ", ".join(
        f"lib{library}={target_sum_by_library[library]:g}"
        for library in libraries))
    return 0


def submit(args: argparse.Namespace) -> int:
    output = os.path.abspath(args.output)
    work_root = os.path.dirname(output)
    script_dir = os.path.join(work_root, "slurm_scripts")
    log_dir = os.path.join(work_root, "logs")
    os.makedirs(script_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    script_token = hashlib.sha256(output.encode("utf-8")).hexdigest()[:12]
    script_path = os.path.join(
        script_dir, f"build_ploidy_h5ad.{script_token}.sbatch")
    forwarded = [
        "--mapping-input-root", os.path.abspath(args.mapping_input_root),
        "--libraries", ",".join(str(value) for value in parse_libraries(args.libraries)),
        "--output", output,
        "--normalization-reference-h5ad",
        os.path.abspath(args.normalization_reference_h5ad),
        "--reference-chunk-cells", str(args.reference_chunk_cells),
        "--normalization-tolerance", str(args.normalization_tolerance),
    ]
    command = " ".join(shlex.quote(value) for value in (
        ["python3", os.path.abspath(__file__)] + forwarded))
    payload = f"""#!/bin/bash
# Generated by tetra_arm_build_ploidy_h5ad.py {RELEASE}.
#SBATCH --job-name=tetarm_ploidy_h5ad
#SBATCH --output={log_dir}/build_ploidy_h5ad_%j.out
#SBATCH --error={log_dir}/build_ploidy_h5ad_%j.err
#SBATCH --partition={args.partition}
#SBATCH --nodes=1
#SBATCH --cpus-per-task={args.cpus}
#SBATCH --mem={args.memory}
#SBATCH --time={args.time}
set -euo pipefail
module purge
module load miniforge/3
module load genomics-base/latest
module list 2>&1
command -v python3
python3 - <<'PY'
import anndata
import h5py
import numpy
import pandas
import scanpy
import scipy
PY
date
hostname
{command}
date
"""
    atomic_text(script_path, payload)
    os.chmod(script_path, os.stat(script_path).st_mode | stat.S_IXUSR)
    syntax = subprocess.run(
        ["bash", "-n", script_path], capture_output=True, text=True,
        check=False)
    if syntax.returncode != 0:
        raise ValueError(
            f"generated sbatch failed bash -n: {syntax.stderr.strip()}")
    result = subprocess.run(
        ["sbatch", "--parsable", "--chdir", work_root, script_path],
        capture_output=True, text=True, check=False)
    job_id = result.stdout.strip().split(";", 1)[0]
    if result.returncode != 0 or not job_id.isdigit():
        detail = result.stderr.strip() or result.stdout.strip()
        raise ValueError(f"sbatch failed: {detail}")
    print(f"Submitted ploidy H5AD build job {job_id}")
    print(f"Script: {script_path}")
    print(f"Logs: {log_dir}")
    print(f"Output after completion: {output}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate the historical normalization transform and build a "
            "newest-remap H5AD for PLOIDY_NN."))
    parser.add_argument("--version", action="version", version=f"%(prog)s {RELEASE}")
    parser.add_argument("--mapping-input-root", default=DEFAULT_MAPPING_ROOT)
    parser.add_argument("--libraries", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--normalization-reference-h5ad", default=DEFAULT_REFERENCE_H5AD,
        help=("Retained H5AD used only to prove its counts-to-X transform; "
              "no reference cells or values are copied"))
    parser.add_argument(
        "--reference-chunk-cells", type=int, default=256,
        help="Reference rows read per chunk; every cell is still checked")
    parser.add_argument("--normalization-tolerance", type=float, default=1e-4)
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--partition", default="compute")
    parser.add_argument("--cpus", type=int, default=8)
    parser.add_argument("--memory", default="128G")
    parser.add_argument("--time", default="1-00:00:00")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if not 1 <= args.reference_chunk_cells <= 4096:
            raise ValueError("--reference-chunk-cells must be between 1 and 4096")
        if not 0 < args.normalization_tolerance <= 0.01:
            raise ValueError("--normalization-tolerance must be in (0,0.01]")
        if args.cpus < 1 or args.cpus > 256:
            raise ValueError("--cpus must be between 1 and 256")
        if not args.memory or any(character.isspace() for character in args.memory):
            raise ValueError("invalid --memory")
        if not args.time or any(character.isspace() for character in args.time):
            raise ValueError("invalid --time")
        parse_libraries(args.libraries)
        return submit(args) if args.submit else build(args)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
