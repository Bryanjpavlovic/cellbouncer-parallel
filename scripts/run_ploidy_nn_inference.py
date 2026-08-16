#!/usr/bin/env python3
"""
run_ploidy_nn_inference.py
Created in conversation: https://claude.ai/chat/<this-conversation-id>

Cluster CPU inference wrapper for the trained ploidy NN. Runs the model
(ploidy_nn_weights.pt + _scaler.npz) on a cluster h5ad and writes per-library
TSVs in the format consumed by `tetra_refine --external_ploidy`.

V1_R5 vs V1_R4: usage examples updated for the current NoMito retrained model.

V1_R3 vs V1_R2: doc/path touch-ups only (V1_R4 of generate_h5ad now produces
both filtered and unfiltered h5ads in /h5a5_outs/, so the example commands
reference real paths).

V1_R2 fixes four bugs from V1_R1:
  1. Encodes batch from a hardcoded LIBRARY_TO_BATCH dict (training used
     adata.obs['batch'], which cluster h5ads don't have).
  2. Writes the column name 'ploidy_probability' (tetra_refine expects this
     literal header), not 'prob_tetraploid'.
  3. Writes ploidy_probability = max(p_tet, 1 - p_tet) (confidence in the
     call, not raw P(tetraploid)).
  4. Aligns HVGs by name from the checkpoint, not from h5ad's own HVG mask.

Adds Option C support: predicts on every cell present in the input h5ad
(filtered or unfiltered), computes a per-cell qc_pass flag using the same
QC thresholds as generate_h5ad_V1_R4.py, and emits qc_pass as an
informational column. Use --qc_only to filter out qc_pass=0 cells before
inference (matches the original training distribution).

Usage:
    # Filtered cells only (matches training distribution)
    python run_ploidy_nn_inference.py \\
        --h5ad /mnt/beegfs/tetmultiome_rna_mapped/mapping_output/h5a5_outs/filtered_normed_tetmultiome_rna.h5ad \\
        --weights /mnt/beegfs/tetmultiome_rna_mapped/ploidy_classifier/retrain_nomito_20260814/model/ploidy_nn_weights.pt \\
        --output_dir /mnt/beegfs/tetmultiome_rna_mapped/ploidy_classifier/retrain_nomito_20260814/ploidy_calls_filtered \\
        --lib_range 1-40 \\
        --qc_only

    # All cells from unfiltered h5ad (Option C; preserves coverage)
    python run_ploidy_nn_inference.py \\
        --h5ad /mnt/beegfs/tetmultiome_rna_mapped/mapping_output/h5a5_outs/unfiltered_normed_tetmultiome_rna.h5ad \\
        --weights /mnt/beegfs/tetmultiome_rna_mapped/ploidy_classifier/retrain_nomito_20260814/model/ploidy_nn_weights.pt \\
        --output_dir /mnt/beegfs/tetmultiome_rna_mapped/ploidy_classifier/retrain_nomito_20260814/ploidy_calls_unfiltered \\
        --lib_range 1-40

Revision History:
V1_R1 (2026-04-30): Initial cluster CPU inference wrapper.
V1_R2 (2026-04-30): Fix batch encoding column, output column names,
  probability semantics. Add qc_pass column and --qc_only flag.
V1_R3 (2026-04-30): Doc/path touch-ups; reflect generate_h5ad_V1_R4 outputs
  in h5a5_outs/ and updated example commands. No code changes.
V1_R4 (2026-05-22): Integration packaging for the current CellBouncer plan.
  Keep the tetra_refine-compatible schema, make library parsing more robust,
  and mark the built-in library-to-batch map as verified against
  Library_conversions.xlsx / LaneQuickGuide.
"""

import argparse
import os
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

# ============================================================================
# CONFIGURATION
# ============================================================================
# Library-to-batch mapping. The trained model uses adata.obs['batch'] for
# batch encoding, but cluster h5ads only have adata.obs['library']. This dict
# replicates the training-time mapping. To populate this, on the desktop run:
#   adata = sc.read_h5ad('/media/b/.../2025_opt_umap_Leiden.h5ad', backed='r')
#   print(adata.obs[['library','batch']].drop_duplicates().sort_values('library'))
# Then fill in the dict below. Wrapper hard-fails at startup if any library
# in the input h5ad is missing from this dict.

LIBRARY_TO_BATCH = {
    # Verified against Library_conversions.xlsx / LaneQuickGuide: libs 1-8
    # Batch1, libs 9-24 Batch2, libs 25-40 Batch3.
    1: "Batch1", 2: "Batch1", 3: "Batch1", 4: "Batch1",
    5: "Batch1", 6: "Batch1", 7: "Batch1", 8: "Batch1",
    9: "Batch2", 10: "Batch2", 11: "Batch2", 12: "Batch2",
    13: "Batch2", 14: "Batch2", 15: "Batch2", 16: "Batch2",
    17: "Batch2", 18: "Batch2", 19: "Batch2", 20: "Batch2",
    21: "Batch2", 22: "Batch2", 23: "Batch2", 24: "Batch2",
    25: "Batch3", 26: "Batch3", 27: "Batch3", 28: "Batch3",
    29: "Batch3", 30: "Batch3", 31: "Batch3", 32: "Batch3",
    33: "Batch3", 34: "Batch3", 35: "Batch3", 36: "Batch3",
    37: "Batch3", 38: "Batch3", 39: "Batch3", 40: "Batch3",
}

# QC thresholds from generate_h5ad_V1_R4.py - used to compute qc_pass flag
QC_MIN_GENES = 200
QC_MAX_PCT_MT = 20.0
QC_MAX_PCT_RIBO = 20.0
QC_MAD_THRESHOLD = 5.0  # for total_counts, total_counts_mt, total_counts_ribo

# ============================================================================
# Lib range parsing
# ============================================================================

def parse_lib_range(s):
    out = []
    for tok in s.split(","):
        tok = tok.strip()
        if "-" in tok:
            lo, hi = tok.split("-", 1)
            out.extend(range(int(lo), int(hi) + 1))
        else:
            out.append(int(tok))
    return sorted(set(out))


def extract_lib_num(library_str):
    """Extract integer library number from common project encodings.

    Handles values such as 19, "19", "lib19", and
    "Tet_2025_Multiome-RNA_19".
    """
    s = str(library_str).strip()
    if s.lower().startswith("lib"):
        s = s[3:]
    else:
        s = s.split("_")[-1]
    return int(s)


# ============================================================================
# Model definition (must match training architecture exactly)
# ============================================================================

def build_model(n_features, n_batches, batch_embed_dim=16):
    import torch
    import torch.nn as nn

    class PloidyNet(nn.Module):
        def __init__(self):
            super().__init__()
            self.batch_embed = nn.Embedding(n_batches, batch_embed_dim)
            input_dim = n_features + batch_embed_dim
            self.net = nn.Sequential(
                nn.Linear(input_dim, 512),
                nn.BatchNorm1d(512),
                nn.ReLU(),
                nn.Dropout(0.3),
                nn.Linear(512, 128),
                nn.BatchNorm1d(128),
                nn.ReLU(),
                nn.Dropout(0.3),
                nn.Linear(128, 32),
                nn.BatchNorm1d(32),
                nn.ReLU(),
                nn.Dropout(0.2),
                nn.Linear(32, 2),
            )

        def forward(self, x, batch_ids):
            b_embed = self.batch_embed(batch_ids)
            combined = torch.cat([x, b_embed], dim=1)
            return self.net(combined)

    return PloidyNet()


# ============================================================================
# Loading
# ============================================================================

def load_model_and_scaler(weights_path):
    import torch

    print(f"Loading checkpoint: {weights_path}")
    checkpoint = torch.load(weights_path, map_location="cpu", weights_only=False)

    n_features = checkpoint["n_features"]
    n_batches = checkpoint["n_batches"]
    hvg_genes = checkpoint["hvg_genes"]
    batch_to_idx = checkpoint["batch_to_idx"]

    model = build_model(n_features, n_batches)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    scaler_path = weights_path.replace(".pt", "_scaler.npz")
    if not os.path.exists(scaler_path):
        print(f"ERROR: Scaler file missing: {scaler_path}", file=sys.stderr)
        sys.exit(1)
    sc_data = np.load(scaler_path)
    scaler_mean = sc_data["mean"]
    scaler_scale = sc_data["scale"]

    print(f"  Features: {n_features}, Batches: {n_batches}, HVGs: {len(hvg_genes)}")
    print(f"  batch_to_idx keys: {sorted(batch_to_idx.keys())}")
    return model, scaler_mean, scaler_scale, hvg_genes, batch_to_idx, n_batches


def compute_qc_pass(adata):
    """Compute per-cell qc_pass using same thresholds as generate_h5ad_V1_R4.py.

    Returns an int8 array of length adata.n_obs: 1 = passes all QC, 0 = fails any.
    Computes QC metrics on the fly if not already in adata.obs.
    """
    import scanpy as sc
    obs = adata.obs

    needed = ["n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo",
              "total_counts", "total_counts_mt", "total_counts_ribo"]
    if any(c not in obs.columns for c in needed):
        print("  Computing QC metrics on the fly...")
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.match("^(RPS|RPL)")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                                   log1p=False, inplace=True)
        sc.pp.calculate_qc_metrics(adata, qc_vars=["ribo"], percent_top=None,
                                   log1p=False, inplace=True)
        obs = adata.obs

    # Hard thresholds
    mask = (obs["n_genes_by_counts"] >= QC_MIN_GENES) & \
           (obs["pct_counts_mt"] < QC_MAX_PCT_MT) & \
           (obs["pct_counts_ribo"] < QC_MAX_PCT_RIBO)

    # MAD outlier filter on three count metrics
    for metric in ("total_counts", "total_counts_mt", "total_counts_ribo"):
        v = obs[metric].values
        med = np.median(v)
        mad = np.median(np.abs(v - med))
        lo = med - QC_MAD_THRESHOLD * mad
        hi = med + QC_MAD_THRESHOLD * mad
        mask &= (v >= lo) & (v <= hi)

    return mask.astype(np.int8)


def load_h5ad_and_align(h5ad_path, hvg_genes, batch_to_idx,
                        lib_filter, qc_only):
    """
    Load h5ad, optionally filter by library, optionally drop qc_pass=0 cells,
    align HVG columns by name, encode batches via LIBRARY_TO_BATCH.
    Returns: X_dense_float32, batch_ids_int, qc_pass_int, obs_df
    """
    import scanpy as sc
    import scipy.sparse

    print(f"Loading h5ad: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    if "library" not in adata.obs.columns:
        print("ERROR: h5ad has no 'library' column", file=sys.stderr)
        sys.exit(1)

    # Filter to requested libraries
    lib_nums_all = adata.obs["library"].astype(str).apply(extract_lib_num).values
    keep = np.isin(lib_nums_all, list(lib_filter))
    if keep.sum() == 0:
        print(f"ERROR: No cells in h5ad match lib_range {sorted(lib_filter)}",
              file=sys.stderr)
        sys.exit(1)
    adata = adata[keep].copy()
    print(f"  After lib filter: {adata.n_obs:,} cells")

    # Validate every library is in the LIBRARY_TO_BATCH dict
    libs_present = set(adata.obs["library"].astype(str).apply(extract_lib_num).unique())
    missing_in_dict = libs_present - set(LIBRARY_TO_BATCH.keys())
    if missing_in_dict:
        print(f"ERROR: Libraries missing from LIBRARY_TO_BATCH dict: "
              f"{sorted(missing_in_dict)}", file=sys.stderr)
        print("       Edit the dict at the top of this script.", file=sys.stderr)
        sys.exit(1)

    # Compute qc_pass for every cell
    qc_pass = compute_qc_pass(adata)
    n_pass = int(qc_pass.sum())
    n_fail = int(adata.n_obs - n_pass)
    print(f"  QC: {n_pass:,} pass, {n_fail:,} fail")

    if qc_only:
        adata = adata[qc_pass.astype(bool)].copy()
        qc_pass = np.ones(adata.n_obs, dtype=np.int8)
        print(f"  --qc_only: dropped failing cells, {adata.n_obs:,} cells remain")

    # Align HVGs by name from the checkpoint
    var_names = list(adata.var_names)
    name_to_col = {n: i for i, n in enumerate(var_names)}
    missing = [g for g in hvg_genes if g not in name_to_col]
    if missing:
        print(f"ERROR: {len(missing)} HVGs from training are missing from "
              f"this h5ad", file=sys.stderr)
        print(f"  First 5 missing: {missing[:5]}", file=sys.stderr)
        sys.exit(1)
    cols = np.array([name_to_col[g] for g in hvg_genes], dtype=np.int64)
    print(f"  HVG alignment OK: {len(hvg_genes)} genes")

    # Densify the HVG sub-matrix
    X_full = adata.X
    if scipy.sparse.issparse(X_full):
        X = X_full[:, cols].toarray().astype(np.float32)
    else:
        X = X_full[:, cols].astype(np.float32)
    print(f"  HVG matrix: {X.shape}")

    # Encode batch via LIBRARY_TO_BATCH lookup, then via checkpoint's batch_to_idx
    lib_nums_keep = adata.obs["library"].astype(str).apply(extract_lib_num).values
    batch_strings = np.array([LIBRARY_TO_BATCH[L] for L in lib_nums_keep])

    unknown_batches = set(batch_strings) - set(batch_to_idx.keys())
    if unknown_batches:
        print(f"ERROR: Batch strings not in checkpoint batch_to_idx: "
              f"{unknown_batches}", file=sys.stderr)
        print(f"       checkpoint batch_to_idx keys: "
              f"{sorted(batch_to_idx.keys())}", file=sys.stderr)
        sys.exit(1)
    batch_ids = np.array([batch_to_idx[b] for b in batch_strings], dtype=np.int64)

    # Build a lean obs frame for output construction
    obs_keep = adata.obs[["library"]].copy()
    obs_keep["lib_num"] = lib_nums_keep
    obs_keep["qc_pass"] = qc_pass
    obs_keep.index = adata.obs.index

    return X, batch_ids, obs_keep


def predict(model, X_scaled, batch_ids, batch_size=2048):
    import torch
    import torch.nn.functional as F

    n = len(X_scaled)
    probs = np.zeros((n, 2), dtype=np.float32)

    with torch.no_grad():
        for start in range(0, n, batch_size):
            end = min(start + batch_size, n)
            xb = torch.from_numpy(X_scaled[start:end])
            bb = torch.from_numpy(batch_ids[start:end])
            logits = model(xb, bb)
            p = F.softmax(logits, dim=1).numpy()
            probs[start:end] = p

    return probs


# ============================================================================
# Output
# ============================================================================

def write_per_library(probs, obs_df, output_dir, libs):
    """
    Output column layout (header line):
      barcode  library  ploidy_call  ploidy_probability  classification_group
      prob_tetraploid  qc_pass

    The first 5 columns are read by tetra_refine. Columns 6-7 are
    informational. ploidy_probability is the confidence in the call (always
    >= 0.5), prob_tetraploid is the raw P(tet).
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    obs_index = obs_df.index.astype(str).values
    barcodes = np.array([s.split("-")[0] for s in obs_index])
    p_tet = probs[:, 1]
    calls = np.where(p_tet >= 0.5, "tetraploid", "diploid")
    # confidence in the call: always >= 0.5
    ploidy_prob = np.where(p_tet >= 0.5, p_tet, 1.0 - p_tet).astype(float)

    df = pd.DataFrame({
        "barcode": barcodes,
        "library": obs_df["lib_num"].values.astype(int),
        "ploidy_call": calls,
        "ploidy_probability": ploidy_prob,
        "classification_group": "group3",
        "prob_tetraploid": p_tet.astype(float),
        "qc_pass": obs_df["qc_pass"].values.astype(int),
    })

    # Combined output (all requested libraries)
    combined_path = os.path.join(output_dir, "all_libs.ploidy_calls_nn.tsv")
    df.to_csv(combined_path, sep="\t", index=False)
    print(f"Wrote combined: {combined_path} ({len(df):,} rows)")

    # Per-library output
    written = []
    for lib in libs:
        sub = df[df["library"] == lib]
        if len(sub) == 0:
            print(f"  lib{lib}: no cells in h5ad, skipping")
            continue
        path = os.path.join(output_dir, f"lib{lib}.ploidy_calls_nn.tsv")
        sub.to_csv(path, sep="\t", index=False)
        n_tet = int((sub["ploidy_call"] == "tetraploid").sum())
        n_dip = int((sub["ploidy_call"] == "diploid").sum())
        n_qc_pass = int(sub["qc_pass"].sum())
        n_qc_fail = len(sub) - n_qc_pass
        print(f"  lib{lib}: tet={n_tet} dip={n_dip} "
              f"(qc_pass={n_qc_pass} qc_fail={n_qc_fail}) -> {path}")
        written.append(lib)

    return written


# ============================================================================
# Main
# ============================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True,
                    help="Path to cluster h5ad (filtered or unfiltered).")
    ap.add_argument("--weights", required=True,
                    help="Path to ploidy_nn_weights.pt; sibling _scaler.npz "
                         "must exist.")
    ap.add_argument("--output_dir", required=True,
                    help="Where to write per-library and combined TSVs.")
    ap.add_argument("--lib_range", default="1-40",
                    help="Libraries to predict (e.g. '13' or '1-40' or "
                         "'3,5,12-20').")
    ap.add_argument("--qc_only", action="store_true",
                    help="Drop cells with qc_pass=0 before inference. "
                         "Default: predict on all cells, propagate qc_pass.")
    ap.add_argument("--force", action="store_true",
                    help="Overwrite existing per-library TSVs.")
    args = ap.parse_args()

    libs = parse_lib_range(args.lib_range)
    print("=" * 70)
    print("  ploidy NN inference V1_R5 (CPU)")
    print("=" * 70)
    print(f"  Started: {datetime.now()}")
    print(f"  Libraries: {libs}")
    print(f"  h5ad: {args.h5ad}")
    print(f"  Weights: {args.weights}")
    print(f"  Output: {args.output_dir}")
    print(f"  Mode: {'QC-only' if args.qc_only else 'all cells (qc_pass propagated)'}")
    print()

    if not args.force:
        all_done = all(os.path.exists(os.path.join(
            args.output_dir, f"lib{L}.ploidy_calls_nn.tsv")) for L in libs)
        if all_done:
            print(f"All {len(libs)} per-library outputs already exist. "
                  f"Use --force to overwrite. Nothing to do.")
            return

    t0 = time.time()
    model, sc_mean, sc_scale, hvg_genes, batch_to_idx, n_batches = \
        load_model_and_scaler(args.weights)

    X, batch_ids, obs_df = load_h5ad_and_align(
        args.h5ad, hvg_genes, batch_to_idx,
        lib_filter=set(libs), qc_only=args.qc_only,
    )

    X_scaled = ((X - sc_mean) / sc_scale).astype(np.float32)

    print("Running inference (CPU)...")
    probs = predict(model, X_scaled, batch_ids)
    print(f"  Inference done in {time.time() - t0:.1f}s")

    write_per_library(probs, obs_df, args.output_dir, libs)
    print(f"\nTotal wall-clock: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
