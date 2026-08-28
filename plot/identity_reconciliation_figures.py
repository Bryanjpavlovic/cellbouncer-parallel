#!/usr/bin/env python3
"""
identity_reconciliation_figures.py

Version: V2_R12 (V2 = this conversation, second design; R12 = twelfth revision)
Provenance: accessory figure generator for the Tet2025 demux -> tetra_refine ->
identity-reconciliation chain. Deployed to the CellBouncer scripts folder
(/nvme/software/packages/cellbouncer/dev/bin/) and callable from
orchestrate_tetraploid.py, or run standalone. Stable pipeline basename by
contract: version lives in this header only, never in the filename.

DESIGN CHANGE FROM V1
---------------------
V1 built delta-only figures (what changed) and three broke on a real lib19-only
run: all-library heatmaps collapsed to one stretched row when a single library
was present, and count-vector similarity was empty on one library. V2 replaces
that framing. The organizing question is now, for EVERY library: which expected
lines were present, which present-but-unexpected, which expected-but-absent,
and how does whole-library demux confidence compare to the cells reconciliation
flagged as wrong. This is the visual language of compare_demux_results_V5,
tetra_refine_layout/plots, tetra_score_distributions, and selection_audit_plots.

Every figure is robust to any library subset, including a single library.

FIGURE FAMILIES
---------------
Per-library pages (one file per library):
  LP  library_page          6-panel: outcome donut, expected-detection bar,
                             all identities by cell count colored by presence
                             status, and tables of unexpected-present and
                             expected-absent identities.
  LS  library_score_page     per-library demux score distributions, whole
                             library vs reassigned/held cells. Needs the demux
                             .diagnostics.gz join; skips gracefully if absent.

All-library aggregates (span every library present; degrade to any subset):
  A1  roster_status_matrix   library x identity presence/absence/unexpected
  A2  library_outcome_bars   per-library outcome composition + non-kept counts
  A3  expected_recovery_bars per-library found/absent/unexpected identity counts
  A4  score_separation       library median demux LLR vs flagged-cell median
  A5  unexpected_catalogue   residual component singlets, true unexpected identities, and expected-missing identities
  A6  full_identity_catalogue every final identity with v11-aware biological relationship
  A7  library_exchange        directed distinguishing-roster exchange coverage plus candidate summary

Ambient swap-test figures (standalone mode, driven by --ambient-swap-spec):
  AS  fixed-profile assignment scatter, fixed-profile KDE by stratum, and
      joint-profile-versus-fixed-profile delta. The same three-figure family is
      written for all candidates across all libraries, for all candidates
      together within each library, and separately for each proposed identity.
      Figure files use descriptive scope/identity/condition prefixes in a flat
      per-library directory; identities and conditions never create directories.
      Plotting and summary ownership lives here; the orchestrator only schedules
      estimator jobs and supplies a compact JSON input manifest.

DATA CONTRACT (under --reconciliation-root)
-------------------------------------------
  decisions/all_libraries.reconciled_cells.tsv.gz
  decisions/all_libraries.identity_events.tsv   (optional)
  decisions/all_libraries.library_exchange_events.tsv  (v11; optional but used by A7)
  metadata/library_expected_genotypes.tsv
Optional per-library score sources (LS, A4):
  <demux-root>/Tet_2025_Multiome-RNA_{N}/demux_nomito/lib{N}_demuxed.diagnostics.gz

Ambient swap-test mode requires a JSON manifest written by
orchestrate_tetraploid.py. Each condition supplies four matched contamination
rate files: G/H use one frozen empty-drop profile, while J/K independently fit
the profile with original/proposed assignments. Candidate groups are discovered
upstream from reconciliation evidence; this script never accepts a handwritten
candidate roster.

If reconciliation column names differ from the documented v11 schema, run
harvest_reconciliation_schema.py and pass --schema <harvest.json>.

Usage (cluster)
---------------
  module purge
  module load miniforge/3
  module load genomics-base/latest
  python3 /nvme/software/packages/cellbouncer/dev/bin/identity_reconciliation_figures.py \
      --reconciliation-root /mnt/beegfs/tetmultiome_rna_mapped/ploidy_classifier/retrain_nomito_20260814/identity_reconciliation \
      --demux-root /mnt/beegfs/tetmultiome_rna_mapped/mapping_output \
      --output-dir /mnt/beegfs/tetmultiome_rna_mapped/ploidy_classifier/retrain_nomito_20260814/identity_reconciliation_figures \
      --figures all

Ambient swap-test figures are normally scheduled by AMBIENT_SWAP_TEST. They
can also be regenerated from its saved manifest without rerunning estimators:
  python3 /nvme/software/packages/cellbouncer/dev/bin/identity_reconciliation_figures.py \
      --ambient-swap-spec <ambient_swap_test_*.worker.json> \
      --output-dir <AMBIENT_SWAP_TEST_output_directory>

Revision history at bottom.

"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import shutil
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap, BoundaryNorm


AMBIENT_PLOT_MAX_POINTS = 50000


# =============================================================================
# Style (plotting_style_guide_V5_R1: Dialect A, semantic colors, manual legends)
# =============================================================================

PAL = {
    "primary": "#2E86AB", "secondary": "#F18F01", "tertiary": "#C73E1D",
    "composite": "#6C5B7B", "good": "#2A9D8F", "warn": "#F4A261",
    "bad": "#E76F51", "accent": "#C73E1D", "neutral": "#888888",
    "bg_light": "#F8F9FA", "grid": "#E0E0E0", "text_dark": "#2B2D42",
    "footnote": "#444444",
}

STATUS_COLORS = {
    "expected_present": "#2A9D8F",
    "component_present": "#2E86AB",
    "expected_absent": "#FFD93D",
    "unexpected_present": "#E76F51",
    "absent": "#EDEDED",
}
STATUS_PRETTY = {
    "expected_present": "Expected exact identity & present",
    "component_present": "Residual expected-composite component singlet",
    "expected_absent": "Expected exact identity, absent",
    "unexpected_present": "Truly unexpected identity, present",
}
INCORRECT_COLOR = "#FF6B6B"

DECISION_COLORS = {
    "KEEP": PAL["neutral"],
    "APPLIED": PAL["good"],
    "REVIEW_UNEXPECTED_IDENTITY": PAL["bad"],
    "REVIEW_CELLULAR_ORIGIN": PAL["warn"],
    "REVIEW_HOMOTET_OCCUPANCY": "#B5651D",
    "OTHER_REVIEW": PAL["composite"],
}
DECISION_ORDER = ["KEEP", "APPLIED", "REVIEW_UNEXPECTED_IDENTITY",
                  "REVIEW_CELLULAR_ORIGIN", "REVIEW_HOMOTET_OCCUPANCY",
                  "OTHER_REVIEW"]
DECISION_PRETTY = {
    "KEEP": "Kept",
    "APPLIED": "Applied change",
    "REVIEW_UNEXPECTED_IDENTITY": "Held: unexpected identity",
    "REVIEW_CELLULAR_ORIGIN": "Held: cellular origin",
    "REVIEW_HOMOTET_OCCUPANCY": "Held: homotet occupancy",
    "OTHER_REVIEW": "Held: other review",
}

BATCH_OF_LIB = {**{n: 1 for n in range(1, 9)},
                **{n: 2 for n in range(9, 25)},
                **{n: 3 for n in range(25, 41)}}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 10, "axes.linewidth": 0.8, "axes.edgecolor": "#333333",
    "axes.labelsize": 11, "axes.titlesize": 12.5, "xtick.labelsize": 9,
    "ytick.labelsize": 9, "legend.fontsize": 8.5, "figure.dpi": 150,
    "savefig.dpi": 200,
})


def style_ax(ax, title=None, xlabel=None, ylabel=None, grid_axis="y",
             title_pad=14):
    if title:
        ax.set_title(title, fontweight="bold", pad=title_pad,
                     color=PAL["text_dark"])
    if xlabel:
        ax.set_xlabel(xlabel, color=PAL["text_dark"])
    if ylabel:
        ax.set_ylabel(ylabel, color=PAL["text_dark"])
    ax.set_facecolor(PAL["bg_light"])
    if grid_axis:
        ax.grid(axis=grid_axis, color=PAL["grid"], linewidth=0.5, alpha=0.7)
    ax.tick_params(colors=PAL["text_dark"])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def help_text(ax, lines, y=-0.30, x=0.0):
    ax.text(x, y, "\n".join(lines), transform=ax.transAxes, ha="left",
            va="top", fontsize=8.2, style="italic", color=PAL["footnote"],
            linespacing=1.42, clip_on=False)


def wrap_identity(label, max_chars=16):
    parts = [p.strip() for p in str(label).split("+") if p.strip()]
    if len(parts) < 2 or len(str(label)) <= max_chars:
        return str(label)
    out, cur = [], parts[0]
    for p in parts[1:]:
        cand = f"{cur}+{p}"
        if len(cand) <= max_chars:
            cur = cand
        else:
            out.append(cur + "+")
            cur = p
    out.append(cur)
    return "\n".join(out)


def save_fig(fig, out_dir, name, dpi=200):
    path = os.path.join(out_dir, name)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {path}")
    return path


# =============================================================================
# Column binding
# =============================================================================

DEFAULT_COLUMNS = {
    "library": "library", "barcode": "barcode",
    "original_demux_assignment": "original_demux_assignment",
    "original_demux_type": "original_demux_type",
    "current_refined_assignment": "current_refined_assignment",
    "reconciled_donor_genotype": "reconciled_donor_genotype",
    "proposed_donor_genotype": "proposed_donor_genotype",
    "occupancy_resolution_status": "occupancy_resolution_status",
    "explicit_multiplet_evidence": "explicit_multiplet_evidence",
    "nn_prob_tetraploid": "nn_prob_tetraploid",
    "nuclear_delta_ll": "nuclear_delta_ll",
    "nuclear_informative_depth": "nuclear_informative_depth",
    "final_action": "final_action",
    "reassignment_applied": "reassignment_applied",
    "decision_reason_codes": "decision_reason_codes",
    "event_id": "event_id", "event_class": "event_class",
    "singlet_library_relationship": "singlet_library_relationship",
    "expected_composite_context": "expected_composite_context",
    "library_exchange_evidence_eligible": "library_exchange_evidence_eligible",
    "policy_version": "policy_version",
}


def bind_columns(header, schema_json=None):
    mapping = dict(DEFAULT_COLUMNS)
    if schema_json and os.path.isfile(schema_json):
        try:
            j = json.load(open(schema_json))
            real = set(j.get("reconciled_cells_columns", []))
            if "library" not in real and j.get("library_column"):
                mapping["library"] = j["library_column"]
        except Exception as e:
            print(f"  WARNING: could not read schema JSON ({e}); using defaults")
    present = {k: v for k, v in mapping.items() if v in header}
    missing = [k for k in mapping if k not in present]
    if missing:
        print(f"  NOTE: columns not in reconciled table (dependent figures "
              f"degrade): {', '.join(missing)}")
    return present


# =============================================================================
# Derived helpers
# =============================================================================

def normalize_library(series):
    def one(v):
        t = str(v).strip()
        for pre in ("Tet_2025_Multiome-RNA_", "Tet_2025_Multiome-ATAC_"):
            if t.startswith(pre):
                t = t[len(pre):]
        if t.lower().startswith("lib"):
            t = t[3:]
        t = t.split("_")[0].split(".")[0]
        try:
            return int(t)
        except ValueError:
            return -1
    return series.map(one)


def component_set(genotype):
    t = str(genotype).strip()
    if not t or t in ("NA", ".", "nan", "None"):
        return []
    return [p for p in t.split("+") if p.strip()]


def component_count(genotype):
    return len(component_set(genotype))


def is_homotet(g):
    parts = component_set(g)
    return len(parts) == 2 and parts[0] == parts[1]


def to_bool(series):
    return series.astype(str).str.strip().str.upper().isin(
        {"TRUE", "1", "T", "YES"})


def decision_group_of(final_action, applied):
    if bool(applied):
        return "APPLIED"
    fa = str(final_action).strip()
    if fa == "KEEP":
        return "KEEP"
    if fa == "REVIEW_UNEXPECTED_IDENTITY":
        return "REVIEW_UNEXPECTED_IDENTITY"
    if fa == "REVIEW_CELLULAR_ORIGIN":
        return "REVIEW_CELLULAR_ORIGIN"
    if fa == "REVIEW_HOMOTET_OCCUPANCY":
        return "REVIEW_HOMOTET_OCCUPANCY"
    return "OTHER_REVIEW"


# =============================================================================
# Data loading
# =============================================================================

def load_data(root, libraries=None, schema_json=None, line_metadata=None):
    dec = os.path.join(root, "decisions")
    meta = os.path.join(root, "metadata")
    combined_cells = os.path.join(dec, "all_libraries.reconciled_cells.tsv.gz")
    events_path = os.path.join(dec, "all_libraries.identity_events.tsv")
    expected_path = os.path.join(meta, "library_expected_genotypes.tsv")
    if not os.path.isfile(expected_path):
        sys.exit(f"ERROR: required input missing: {expected_path}")

    # Resolve the reconciled-cells source. Prefer the canonical all-libraries
    # aggregate when it truly spans multiple libraries. If a prior subset run
    # left that aggregate containing only one library, fall back to the durable
    # per-library files. If the caller restricted --libraries, read only those.
    import glob as _glob
    perlib = {}
    for pth in _glob.glob(os.path.join(dec, "lib*.reconciled_cells.tsv.gz")):
        base = os.path.basename(pth)
        ln = normalize_library(pd.Series([base.split(".")[0]]))[0]
        if ln >= 0:
            perlib[ln] = pth

    def _lib_count(path):
        try:
            col = pd.read_csv(path, sep="\t", usecols=[0], dtype=str)
            return col.iloc[:, 0].map(lambda v: normalize_library(
                pd.Series([v]))[0]).nunique()
        except Exception:
            return 0

    use_combined = os.path.isfile(combined_cells) and \
        _lib_count(combined_cells) > 1
    if use_combined:
        sources = {"__combined__": combined_cells}
        print(f"Reading consolidated table {combined_cells}")
    elif perlib:
        wanted = set(libraries) if libraries else set(perlib)
        sources = {ln: p for ln, p in sorted(perlib.items()) if ln in wanted}
        if not sources:
            sys.exit(f"ERROR: no per-library reconciled files match the "
                     f"requested libraries {sorted(wanted)} in {dec}")
        print(f"Consolidated table absent or single-library; gathering "
              f"{len(sources)} per-library file(s) from {dec}")
    elif os.path.isfile(combined_cells):
        sources = {"__combined__": combined_cells}
        print(f"Reading {combined_cells} (single library)")
    else:
        sys.exit(f"ERROR: no reconciled-cells input found in {dec} "
                 f"(neither all_libraries.reconciled_cells.tsv.gz nor "
                 f"lib*.reconciled_cells.tsv.gz)")

    # Bind columns from the first source's header, then read every source with
    # the same column set and concatenate.
    first_path = next(iter(sources.values()))
    header = pd.read_csv(first_path, sep="\t", nrows=0).columns.tolist()
    cmap = bind_columns(header, schema_json)
    read_cols = list(dict.fromkeys(cmap.values()))
    frames = []
    for key, pth in sources.items():
        h = pd.read_csv(pth, sep="\t", nrows=0).columns.tolist()
        cols = [c for c in read_cols if c in h]
        frames.append(pd.read_csv(pth, sep="\t", usecols=cols, dtype=str,
                                  low_memory=False))
    cells = pd.concat(frames, ignore_index=True)
    cells = cells.rename(columns={v: k for k, v in cmap.items()})

    cells["lib_num"] = normalize_library(cells["library"])
    dropped = int((cells["lib_num"] < 0).sum())
    if dropped:
        print(f"  WARNING: {dropped} rows had unparseable library labels; "
              f"excluded")
        cells = cells[cells["lib_num"] >= 0].copy()

    for col in ("nn_prob_tetraploid", "nuclear_delta_ll",
                "nuclear_informative_depth"):
        if col in cells.columns:
            cells[col] = pd.to_numeric(cells[col], errors="coerce")
    cells["reassignment_applied"] = to_bool(
        cells.get("reassignment_applied",
                  pd.Series("FALSE", index=cells.index)))
    if "explicit_multiplet_evidence" in cells.columns:
        cells["explicit_multiplet_evidence"] = to_bool(
            cells["explicit_multiplet_evidence"])
    else:
        cells["explicit_multiplet_evidence"] = False
    if "final_action" not in cells.columns:
        cells["final_action"] = np.where(cells["reassignment_applied"],
                                         "REASSIGN_GENOTYPE", "KEEP")
    cells["decision_group"] = [
        decision_group_of(a, b) for a, b in
        zip(cells["final_action"], cells["reassignment_applied"])]
    if "reconciled_donor_genotype" not in cells.columns:
        cells["reconciled_donor_genotype"] = cells.get(
            "current_refined_assignment", cells["original_demux_assignment"])
    cells["final_component_count"] = cells[
        "reconciled_donor_genotype"].map(component_count)
    if "proposed_donor_genotype" in cells.columns:
        cells["proposed_component_count"] = cells[
            "proposed_donor_genotype"].map(component_count)
        cells["is_homotet_proposal"] = cells[
            "proposed_donor_genotype"].map(is_homotet)
    else:
        cells["proposed_component_count"] = 0
        cells["is_homotet_proposal"] = False

    events = None
    ev_frames = []
    combined_events = events_path
    ev_perlib = _glob.glob(os.path.join(dec, "lib*.identity_events.tsv"))
    if os.path.isfile(combined_events):
        try:
            ecol = pd.read_csv(combined_events, sep="\t", usecols=[0],
                               dtype=str)
            multi = ecol.iloc[:, 0].map(lambda v: normalize_library(
                pd.Series([v]))[0]).nunique() > 1
        except Exception:
            multi = False
        if multi or not ev_perlib:
            ev_frames.append(pd.read_csv(combined_events, sep="\t", dtype=str))
    if not ev_frames and ev_perlib:
        for pth in sorted(ev_perlib):
            try:
                ev_frames.append(pd.read_csv(pth, sep="\t", dtype=str))
            except Exception:
                pass
    if ev_frames:
        events = pd.concat(ev_frames, ignore_index=True)
        if "library" in events.columns:
            events["lib_num"] = normalize_library(events["library"])
        for col in ("n_implicated_cells", "fraction_library_implicated",
                    "fraction_primary_source_displaced"):
            if col in events.columns:
                events[col] = pd.to_numeric(events[col], errors="coerce")

    expected = pd.read_csv(expected_path, sep="\t", dtype=str)
    expected["lib_num"] = normalize_library(expected["library"])
    gcol = "canonical_genotype" if "canonical_genotype" in expected.columns \
        else expected.columns[1]
    roster = {int(k): set(v) for k, v in expected.groupby("lib_num")[gcol]}
    donor_roster = {}
    component_context = {}
    for lib, gs in roster.items():
        donors = set()
        by_component = {}
        for g in gs:
            comps = component_set(g)
            donors.update(comps)
            if len(comps) >= 2:
                for donor in set(comps):
                    by_component.setdefault(donor, set()).add(str(g))
        donor_roster[lib] = donors
        component_context[lib] = {
            donor: sorted(values) for donor, values in by_component.items()
        }
    n_expected_libs = expected.groupby(gcol)["lib_num"].nunique().to_dict()

    # v11 library-exchange outputs are pairwise library-level evidence. Keep
    # them separate from cell/event tables so a plotting decision can never
    # accidentally reinterpret shared donors as foreign identities.
    exchange = None
    exchange_path = os.path.join(dec, "all_libraries.library_exchange_events.tsv")
    if os.path.isfile(exchange_path):
        exchange = pd.read_csv(exchange_path, sep="\t", dtype=str)
        for col in ("a_signature_coverage_in_b", "b_signature_coverage_in_a",
                    "a_native_retention_fraction", "b_native_retention_fraction",
                    "a_native_displacement_fraction", "b_native_displacement_fraction"):
            if col in exchange.columns:
                exchange[col] = pd.to_numeric(exchange[col], errors="coerce")
        if "library_a" in exchange.columns:
            exchange["lib_a_num"] = normalize_library(exchange["library_a"])
        if "library_b" in exchange.columns:
            exchange["lib_b_num"] = normalize_library(exchange["library_b"])

    if libraries:
        keep = set(libraries)
        cells = cells[cells["lib_num"].isin(keep)].copy()
        if events is not None and "lib_num" in events.columns:
            events = events[events["lib_num"].isin(keep)].copy()
        if exchange is not None and {"lib_a_num", "lib_b_num"}.issubset(exchange.columns):
            exchange = exchange[
                exchange["lib_a_num"].isin(keep) & exchange["lib_b_num"].isin(keep)
            ].copy()

    def final_relationship(lib, genotype):
        """Classify the FINAL reconciled identity using v11 biology.

        Exact expected genotypes are expected. A one-donor genotype that is
        not explicitly expected but is a donor component of an expected
        composite is a residual/component-derived singlet, not foreign. Only
        identities absent from the complete expected biology are unexpected.
        """
        gs = str(genotype).strip()
        if not gs or gs in ("NA", ".", "nan"):
            return "absent"
        lib = int(lib)
        if gs in roster.get(lib, set()):
            return "expected"
        comps = component_set(gs)
        if len(comps) == 1 and comps[0] in donor_roster.get(lib, set()):
            return "component"
        return "unexpected"

    cells["final_identity_relationship"] = [
        final_relationship(l, g) for l, g in
        zip(cells["lib_num"], cells["reconciled_donor_genotype"])]
    cells["is_component_derived_final"] = cells[
        "final_identity_relationship"].eq("component")
    cells["is_unexpected_final"] = cells[
        "final_identity_relationship"].eq("unexpected")

    libs_present = sorted(cells["lib_num"].unique())
    policy = cells["policy_version"].dropna().unique().tolist() \
        if "policy_version" in cells.columns else []
    print(f"Loaded {len(cells):,} cells across {len(libs_present)} "
          f"librar{'y' if len(libs_present) == 1 else 'ies'} "
          f"{libs_present if len(libs_present) <= 8 else ''}; "
          f"{int(cells['reassignment_applied'].sum()):,} applied; "
          f"policy={policy}")
    meta_map = load_line_metadata(line_metadata)
    if meta_map:
        print(f"Loaded line metadata for {len(meta_map)} genotype keys "
              f"(sort by CorrectedFZGRP then UID)")
    return {"cells": cells, "events": events, "exchange": exchange,
            "expected": expected, "roster": roster,
            "donor_roster": donor_roster, "component_context": component_context,
            "n_expected_libs": n_expected_libs, "libs_present": libs_present,
            "policy": policy, "meta_map": meta_map}


def load_line_metadata(xlsx_path):
    """Read Library_conversions.xlsx and return a genotype-string -> metadata
    map used to sort identities by CorrectedFZGRP then UID, split diploid from
    tetraploid.

    Verified against the workbook: the 2025_LineMeta sheet's WGS_Key column
    holds genotypes in the same A+B form as the reconciled genotypes
    (C40210+H27322, H20961+H20961) and single tokens for diploids (H20961).
    Each WGS_Key maps to exactly one CorrectedFZGRP (0/188 keys ambiguous), so
    FZGRP reliably groups fusion partners. Partner orientation varies
    (H20961+C3624 and C3624+H20961 both exist under FZGRP 11), so lookups try
    the key as written and the reversed two-partner order. Diploids carry
    Tet.Class == 'Diploid' and a NaN FZGRP; they are grouped separately and
    ordered by their own key.

    Returns dict: genotype_str -> {fzgrp: float|inf, uid: int|inf,
    is_diploid: bool, matched: bool}. Unmatched genotypes get inf/inf so they
    sort last within their diploid/tetraploid block, deterministically by name.
    """
    if not xlsx_path or not os.path.isfile(xlsx_path):
        return None
    try:
        lm = pd.read_excel(xlsx_path, sheet_name="2025_LineMeta")
    except Exception as e:
        print(f"  WARNING: could not read line metadata ({e}); identities will "
              f"sort by cell mass instead of FZGRP/UID")
        return None
    lm["WGS_Key"] = lm["WGS_Key"].astype(str).str.strip()
    key_to_meta = {}
    for _, r in lm.iterrows():
        k = r["WGS_Key"]
        if k in ("", "nan", "None"):
            continue
        try:
            fz = float(r.get("CorrectedFZGRP")) if pd.notna(
                r.get("CorrectedFZGRP")) else float("inf")
        except (TypeError, ValueError):
            fz = float("inf")
        try:
            uid = int(r.get("UID")) if pd.notna(r.get("UID")) else 10 ** 9
        except (TypeError, ValueError):
            uid = 10 ** 9
        is_dip = str(r.get("Tet.Class")).strip() == "Diploid" or \
            component_count(k) < 2
        # keep the smallest UID seen for a key so the row order is stable
        prev = key_to_meta.get(k)
        if prev is None or uid < prev["uid"]:
            key_to_meta[k] = {"fzgrp": fz, "uid": uid, "is_diploid": is_dip,
                              "matched": True}
    return key_to_meta


def genotype_meta(meta_map, genotype):
    """Resolve one reconciled genotype string to its metadata, trying the key
    as written then the reversed two-partner order. Falls back to a component-
    count diploid judgement when the genotype is not in the workbook."""
    g = str(genotype).strip()
    default = {"fzgrp": float("inf"), "uid": 10 ** 9,
               "is_diploid": component_count(g) < 2, "matched": False}
    if not meta_map:
        return default
    if g in meta_map:
        return meta_map[g]
    parts = [p for p in g.split("+") if p.strip()]
    if len(parts) == 2:
        rev = "+".join(parts[::-1])
        if rev in meta_map:
            return meta_map[rev]
    return default


def identity_sort_key(meta_map, genotype):
    """Sort key implementing: diploids grouped separately from tetraploids,
    then by CorrectedFZGRP, then UID, then genotype string. Diploids sort as a
    block before tetraploids; unmatched lines fall to the end of their block via
    inf FZGRP/UID but stay deterministic by name."""
    m = genotype_meta(meta_map, genotype)
    ploidy_block = 0 if m["is_diploid"] else 1
    return (ploidy_block, m["fzgrp"], m["uid"], str(genotype))


def load_library_scores(lib, args):
    if not args.demux_root:
        return None
    diag = os.path.join(args.demux_root, f"Tet_2025_Multiome-RNA_{lib}",
                        "demux_nomito", f"lib{lib}_demuxed.diagnostics.gz")
    if not os.path.isfile(diag):
        return None
    try:
        df = pd.read_csv(diag, sep="\t", dtype=str, low_memory=False)
    except Exception as e:
        print(f"    lib{lib}: could not read diagnostics ({e})")
        return None
    if "barcode" not in df.columns:
        return None
    # This diagnostics file belongs to exactly one library. Tag it so any
    # downstream merge keys on (lib_num, barcode), never barcode alone:
    # barcodes are NOT unique across libraries (~5% collide), so a barcode-only
    # join would cross-contaminate the moment more than one library is in play.
    df["lib_num"] = lib
    # Guard against a barcode appearing twice in one diagnostics file, which
    # would duplicate cells on an inner join.
    before = len(df)
    df = df.drop_duplicates(subset="barcode", keep="first")
    if len(df) < before:
        print(f"    lib{lib}: dropped {before - len(df)} duplicate-barcode "
              f"diagnostics rows before join")
    for col in ("llr_vs_runner_up", "min_margin", "total_depth",
                "margin_softmax_score", "margin_entropy"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def fig_footer(fig, data, extra=""):
    n = len(data["cells"])
    applied = int(data["cells"]["reassignment_applied"].sum())
    libs = len(data["libs_present"])
    policies = data.get("policy", [])
    policy_txt = ",".join(
        str(x).replace("identity_reconciliation_policy_", "")
        for x in policies
    ) if policies else "unknown"
    txt = (f"n = {n:,} reconciled cells across {libs} librar"
           f"{'y' if libs == 1 else 'ies'}; {applied:,} applied genotype "
           f"changes; policy {policy_txt}. {extra}")
    fig.text(0.01, 0.004, txt, fontsize=8.0, style="italic",
             color=PAL["footnote"], ha="left", va="bottom")


def library_status_table(data, lib):
    """Final identity presence classified with v11 biological expectedness.

    A one-donor final identity that is not an exact roster genotype but is a
    component of an expected composite is a residual/component-derived
    singlet. It is intentionally separate from true foreign/unexpected
    identities.
    """
    cells = data["cells"]
    roster = data["roster"].get(lib, set())
    donor_roster = data["donor_roster"].get(lib, set())
    sub = cells[cells["lib_num"] == lib]
    present_counts = sub["reconciled_donor_genotype"].value_counts().to_dict()
    present = {g: c for g, c in present_counts.items()
               if str(g).strip() not in ("", "NA", ".", "nan")}
    table = {}
    for g in roster:
        table[g] = ("expected_present", present[g]) if present.get(g, 0) > 0 \
            else ("expected_absent", 0)
    for g, c in present.items():
        if g in roster:
            continue
        comps = component_set(g)
        if len(comps) == 1 and comps[0] in donor_roster:
            table[g] = ("component_present", c)
        else:
            table[g] = ("unexpected_present", c)
    return table


def draw_batch_rules(ax, libs, axis="x"):
    for j in range(len(libs) - 1):
        if BATCH_OF_LIB.get(libs[j], 0) != BATCH_OF_LIB.get(libs[j + 1], 0):
            if axis == "x":
                ax.axvline(j + 0.5, color=PAL["text_dark"], linewidth=1.0)
            else:
                ax.axhline(j + 0.5, color=PAL["text_dark"], linewidth=1.0)


# =============================================================================
# LP  Per-library everything page
# =============================================================================

def figure_lp(data, out_dir, args):
    return [_library_page(data, lib, out_dir, args)
            for lib in data["libs_present"]]


def _library_page(data, lib, out_dir, args):
    cells = data["cells"]
    sub = cells[cells["lib_num"] == lib]
    n = len(sub)
    status = library_status_table(data, lib)
    roster = data["roster"].get(lib, set())
    n_present = sum(1 for s, _ in status.values() if s == "expected_present")
    n_component = sum(1 for s, _ in status.values() if s == "component_present")
    n_absent = sum(1 for s, _ in status.values() if s == "expected_absent")
    n_unexp = sum(1 for s, _ in status.values() if s == "unexpected_present")
    grp = sub["decision_group"]
    n_applied = int((grp == "APPLIED").sum())
    held_classes = [g for g in DECISION_ORDER if g not in ("KEEP", "APPLIED")]
    n_held = int(grp.isin(held_classes).sum())
    n_kept = int((grp == "KEEP").sum())
    batch = BATCH_OF_LIB.get(lib, 0)

    # Figure height scales with the identity count so the bar panel never
    # crushes rows together (lib5-class libraries have ~50 identities). Rows get
    # a fixed vertical pitch in inches; two-line wrapped labels need extra room,
    # so the pitch accounts for the longest label's line count.
    m = len(status)
    max_lines = max((wrap_identity(k).count("\n") + 1 for k in status), default=1)
    row_in = 0.24 + (0.14 if max_lines > 1 else 0.0)
    bar_h = max(3.0, m * row_in)
    summary_h = 3.4
    table_h = 3.0
    total_h = summary_h + bar_h + table_h
    fig = plt.figure(figsize=(16, total_h))
    gs = GridSpec(3, 2, height_ratios=[summary_h, bar_h, table_h],
                  hspace=0.45 * (10.0 / total_h) + 0.18, wspace=0.24,
                  left=0.07, right=0.965,
                  top=1 - 1.05 / total_h, bottom=0.7 / total_h)

    ax1 = fig.add_subplot(gs[0, 0])
    parts = [(n_kept, "Kept", PAL["neutral"]),
             (n_applied, "Applied", PAL["good"]),
             (n_held, "Held", PAL["warn"])]
    keep = [(p, l, c) for p, l, c in parts if p > 0]
    if keep:
        vals, labs, cols = zip(*keep)
        total = sum(vals)
        # No wedge labels or autopct text (they collide on tiny slices); put
        # the counts and percents in a legend to the side instead.
        ax1.pie(vals, colors=cols, startangle=90,
                wedgeprops=dict(width=0.42, edgecolor="white"))
        leg_handles = [Patch(facecolor=c, edgecolor="white",
                             label=f"{l}: {p:,} ({p / total * 100:.0f}%)")
                       for p, l, c in keep]
        ax1.legend(handles=leg_handles, loc="center left",
                   bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=9,
                   handlelength=1.2)
    ax1.set_title(f"Cell outcomes  (n={n:,})", fontweight="bold",
                  fontsize=11.5, pad=10, color=PAL["text_dark"])

    ax2 = fig.add_subplot(gs[0, 1])
    cats = ["Expected exact\n& present", "Residual component\nsinglet",
            "Expected exact,\nabsent", "Truly unexpected,\npresent"]
    vals = [n_present, n_component, n_absent, n_unexp]
    cols = [STATUS_COLORS["expected_present"], STATUS_COLORS["component_present"],
            STATUS_COLORS["expected_absent"], STATUS_COLORS["unexpected_present"]]
    bars = ax2.bar(cats, vals, color=cols, edgecolor="black", linewidth=0.5)
    for b, v in zip(bars, vals):
        ax2.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.05,
                 str(v), ha="center", va="bottom", fontweight="bold",
                 fontsize=10)
    ax2.set_ylim(0, max(vals + [1]) * 1.22)
    style_ax(ax2, title="Identity relationship to expected library biology", ylabel="Identity count")

    ax3 = fig.add_subplot(gs[1, :])
    rows = sorted(status.items(), key=lambda kv: kv[1][1])
    ident = [wrap_identity(k) for k, _ in rows]
    counts = [v[1] for _, v in rows]
    scolor = [STATUS_COLORS[v[0]] for _, v in rows]
    y = np.arange(m)
    ax3.barh(y, counts, color=scolor, edgecolor="black", linewidth=0.4,
             height=0.78)
    ax3.set_yticks(y)
    # With dynamic figure height each row has a fixed ~0.24in, so a constant
    # readable label size works without collisions.
    fs = 7.5 if m <= 45 else 6.5
    ax3.set_yticklabels(ident, fontsize=fs)
    ax3.set_ylim(-0.7, m - 0.3)
    cmax = max(counts) if counts else 1
    ax3.set_xlim(0, cmax * 1.20)
    for yy, c in zip(y, counts):
        if c > 0:
            ax3.text(c + cmax * 0.008, yy, str(int(c)), va="center",
                     fontsize=fs, color=PAL["footnote"])
    style_ax(ax3, title="Every roster / present identity by cell count",
             xlabel="Cells", grid_axis="x")
    handles = [Patch(facecolor=STATUS_COLORS[k], edgecolor="black",
                     label=STATUS_PRETTY[k]) for k in
               ("expected_present", "component_present", "expected_absent",
                "unexpected_present")]
    ax3.legend(handles=handles, loc="lower right", frameon=True,
               facecolor="white", edgecolor="none", framealpha=0.85,
               fontsize=8.5)

    ax5 = fig.add_subplot(gs[2, 0])
    ax5.axis("off")
    ax5.set_title("Non-roster identities present", fontweight="bold",
                  fontsize=10.5, pad=6, color=PAL["text_dark"])
    non_roster = []
    for k, v in status.items():
        if v[0] == "component_present":
            non_roster.append((k, "Residual component", v[1]))
        elif v[0] == "unexpected_present":
            non_roster.append((k, "Unexpected", v[1]))
    non_roster.sort(key=lambda x: (-x[2], x[0]))
    if non_roster:
        rows_t = [[k, cls, f"{v:,}"] for k, cls, v in non_roster[:18]]
        tab = ax5.table(cellText=rows_t,
                        colLabels=["Identity", "Relationship", "Cells"],
                        loc="upper center", cellLoc="left",
                        colColours=[PAL["bg_light"]] * 3)
        tab.auto_set_font_size(False)
        tab.set_fontsize(max(6, min(8, 130 / max(len(rows_t), 1))))
        tab.scale(1.0, 1.18)
        if len(non_roster) > 18:
            ax5.text(0.5, 0.0, f"+ {len(non_roster) - 18} more", ha="center",
                     transform=ax5.transAxes, fontsize=8,
                     color=PAL["footnote"])
    else:
        ax5.text(0.5, 0.55, "None: all present identities are exact expected lines",
                 ha="center", va="center", fontsize=10.5, color=PAL["good"],
                 fontweight="bold")

    ax6 = fig.add_subplot(gs[2, 1])
    ax6.axis("off")
    ax6.set_title("Expected identities absent", fontweight="bold",
                  fontsize=10.5, pad=6, color=PAL["text_dark"])
    absent = sorted([k for k, v in status.items()
                     if v[0] == "expected_absent"])
    if absent:
        rows_t = [[k] for k in absent[:22]]
        tab = ax6.table(cellText=rows_t,
                        colLabels=["Expected but not present"],
                        loc="upper center", cellLoc="left",
                        colColours=[STATUS_COLORS["expected_absent"]])
        tab.auto_set_font_size(False)
        tab.set_fontsize(max(6, min(8, 130 / max(len(rows_t), 1))))
        tab.scale(1.0, 1.18)
        if len(absent) > 22:
            ax6.text(0.5, 0.0, f"+ {len(absent) - 22} more", ha="center",
                     transform=ax6.transAxes, fontsize=8,
                     color=PAL["footnote"])
    else:
        ax6.text(0.5, 0.55, "None: every expected identity present",
                 ha="center", va="center", fontsize=10.5, color=PAL["good"],
                 fontweight="bold")

    fig.suptitle(f"lib{lib} identity reconciliation summary  (batch {batch})",
                 fontweight="bold", fontsize=16, y=1 - 0.32 / total_h)
    fig.text(0.5, 1 - 0.68 / total_h,
             f"{len(roster)} expected exact identities: {n_present} present, "
             f"{n_absent} absent; {n_component} residual component singlets; "
             f"{n_unexp} truly unexpected present; {n_applied:,} applied, "
             f"{n_held:,} held",
             ha="center", fontsize=10.5, color=PAL["footnote"])
    fig.text(0.07, 0.22 / total_h,
             "Present means the identity is the final reconciled genotype of at "
             "least one cell. Residual component singlets are not exact roster "
             "genotypes but their donor is represented in an expected composite "
             "line, so v11 does not treat them as foreign. Truly unexpected "
             "identities are absent from the complete expected library biology. "
             "Outcomes use final_action / reassignment_applied.",
             fontsize=8.2, style="italic", color=PAL["footnote"],
             linespacing=1.4)
    return save_fig(fig, out_dir, f"LP_lib{lib:02d}_page.png", args.dpi)


# =============================================================================
# LS  Per-library score page
# =============================================================================

SCORE_SPECS = [
    ("llr_vs_runner_up", "Demux LLR vs runner-up", "log"),
    ("min_margin", "Demux min margin", "log"),
    ("total_depth", "Demux total depth", "log"),
    ("margin_softmax_score", "Demux posterior", "linear"),
    ("margin_entropy", "Demux entropy (bits)", "linear"),
]


def figure_ls(data, out_dir, args):
    written = []
    for lib in data["libs_present"]:
        path = _library_score_page(data, lib, out_dir, args)
        if path:
            written.append(path)
    # Global page pooling every library.
    gpath = _global_score_page(data, out_dir, args)
    if gpath:
        written.append(gpath)
    if not written:
        print("  LS produced no pages: demux diagnostics not found under "
              "--demux-root. Confirm the path with harvest_reconciliation_"
              "schema.py or pass the correct --demux-root.")
    return written


def _grid_dims(n):
    cols = int(np.ceil(np.sqrt(n)))
    rows = int(np.ceil(n / cols))
    return rows, cols


def _render_score_grid(merged, wrong, present, title, subtitle, out_name,
                       out_dir, data, dpi):
    """Render a near-square grid of score panels, each overlaying the whole set
    (grey) and the reconciliation-flagged subset (red). Grid layout keeps log
    axes from sitting horizontally adjacent at cramped spacing, which is what
    made the old single-column version look poor."""
    from matplotlib.ticker import (MaxNLocator, LogLocator, NullFormatter,
                                   NullLocator)
    n = len(present)
    # Cap at 2 columns: wide log/linear panels need horizontal room, and 3-wide
    # rows put panel-edge tick labels too close regardless of scale ordering.
    cols = min(2, n)
    rows = int(np.ceil(n / cols))
    # Interleave log-scaled and linear panels so two log axes are never placed
    # horizontally adjacent in the grid; adjacent log exponent labels at panel
    # edges are the source of tick collisions.
    logs = [p for p in present if p[2] == "log"]
    lins = [p for p in present if p[2] != "log"]
    ordered, li, ki = [], 0, 0
    take_log = True
    while li < len(logs) or ki < len(lins):
        if take_log and li < len(logs):
            ordered.append(logs[li]); li += 1
        elif ki < len(lins):
            ordered.append(lins[ki]); ki += 1
        elif li < len(logs):
            ordered.append(logs[li]); li += 1
        take_log = not take_log
    present = ordered
    fig = plt.figure(figsize=(cols * 5.4, rows * 3.3 + 1.6))
    gs = GridSpec(rows, cols, left=0.09, right=0.965,
                  top=1 - 1.7 / (rows * 3.3 + 1.6),
                  bottom=1.35 / (rows * 3.3 + 1.6), wspace=0.62, hspace=0.58)
    for k, (col, lab, scale) in enumerate(present):
        r, c = divmod(k, cols)
        ax = fig.add_subplot(gs[r, c])
        allv = pd.to_numeric(merged[col], errors="coerce").dropna()
        wrongv = pd.to_numeric(wrong[col], errors="coerce").dropna()
        if scale == "log":
            allv = allv[allv > 0]
            wrongv = wrongv[wrongv > 0]
            if len(allv):
                lo, hi = np.log10(allv.min()), np.log10(allv.max() + 1e-9)
                bins = np.logspace(lo, hi, 34)
                ax.set_xscale("log")
                ax.set_xlim(10 ** (lo - 0.3), 10 ** (hi + 0.7))
                # only label decades actually inside the data range, so the
                # right padding never carries a spurious edge exponent that
                # meets the neighbouring linear panel
                decades = [10.0 ** d for d in
                           range(int(np.floor(lo)), int(np.floor(hi)) + 1)]
                if decades:
                    ax.set_xticks(decades)
                else:
                    ax.xaxis.set_major_locator(LogLocator(numticks=2))
                ax.xaxis.set_minor_locator(NullLocator())
                ax.xaxis.set_minor_formatter(NullFormatter())
            else:
                bins = 30
        else:
            lo = float(min(allv.min() if len(allv) else 0,
                           wrongv.min() if len(wrongv) else 0))
            hi = float(max(allv.max() if len(allv) else 1,
                           wrongv.max() if len(wrongv) else 1))
            span = (hi - lo) or 1.0
            bins = np.linspace(lo, hi + 1e-9, 34)
            # left padding plus lower-tick prune so no label sits at the left
            # edge where it would meet a log panel's right exponent
            ax.set_xlim(lo - span * 0.22, hi + span * 0.08)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4, prune="lower"))
        if len(allv):
            ax.hist(allv, bins=bins, color=PAL["neutral"], alpha=0.5,
                    density=True, label="Whole library", edgecolor="white",
                    linewidth=0.2)
        if len(wrongv):
            ax.hist(wrongv, bins=bins, color=INCORRECT_COLOR, alpha=0.62,
                    density=True, label="Reassigned / held", edgecolor="white",
                    linewidth=0.2)
            amed = float(allv.median()) if len(allv) else None
            wmed = float(wrongv.median())
            ax.axvline(wmed, color=PAL["accent"], linewidth=1.3,
                       label="Flagged median")
            if amed is not None:
                ax.axvline(amed, color=PAL["text_dark"], linewidth=1.0,
                           linestyle="--", label="Library median")
        style_ax(ax, xlabel=lab, ylabel="Density", grid_axis="y")
        ax.set_title(f"{lab}  (n={len(allv):,})", fontweight="bold",
                     fontsize=10.5, pad=8, color=PAL["text_dark"])
    # blank any unused cells
    for k in range(n, rows * cols):
        r, c = divmod(k, cols)
        fig.add_subplot(gs[r, c]).axis("off")
    # one shared legend in the first panel
    handles = [Patch(facecolor=PAL["neutral"], alpha=0.5, label="Whole set"),
               Patch(facecolor=INCORRECT_COLOR, alpha=0.62,
                     label="Reassigned / held"),
               Line2D([0], [0], color=PAL["accent"], linewidth=1.3,
                      label="Flagged median"),
               Line2D([0], [0], color=PAL["text_dark"], linewidth=1.0,
                      linestyle="--", label="Set median")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=9.5, bbox_to_anchor=(0.5, 0.03))
    figh = rows * 3.3 + 1.6
    fig.suptitle(title, fontweight="bold", fontsize=15, y=1 - 0.5 / figh)
    fig.text(0.5, 1 - 1.05 / figh, subtitle, ha="center",
             fontsize=10.5, color=PAL["footnote"])
    fig_footer(fig, data)
    return save_fig(fig, out_dir, out_name, dpi)


def _library_score_page(data, lib, out_dir, args):
    cells = data["cells"]
    scores = load_library_scores(lib, args)
    if scores is None:
        return None
    sub = cells[cells["lib_num"] == lib].copy()
    # Key the join on (lib_num, barcode), never barcode alone: barcodes are not
    # unique across libraries. Both sides carry lib_num, so a foreign-library
    # barcode collision cannot match here.
    n_sub = len(sub)
    merged = sub.merge(scores, on=["lib_num", "barcode"], how="inner",
                       suffixes=("", "_diag"))
    if len(merged) > n_sub:
        raise RuntimeError(
            f"lib{lib}: score join inflated rows {n_sub} -> {len(merged)}; "
            f"duplicate keys present, aborting rather than emit a corrupted "
            f"figure")
    print(f"    lib{lib}: {len(merged):,}/{n_sub:,} cells joined to demux "
          f"diagnostics")
    wrong = merged[merged["decision_group"] != "KEEP"]
    present = [(c, lab, sc) for c, lab, sc in SCORE_SPECS
               if c in merged.columns]
    if not present:
        print(f"    lib{lib}: diagnostics lack expected score columns; skipping")
        return None
    return _render_score_grid(
        merged, wrong, present,
        title=f"lib{lib} demux confidence: library vs flagged cells",
        subtitle=f"{len(merged):,} cells; {len(wrong):,} flagged "
                 f"(reassigned or held). Flagged cells overlapping the strong "
                 f"bulk mean confident demux calls reconciliation still "
                 f"corrected.",
        out_name=f"LS_lib{lib:02d}_scores.png", out_dir=out_dir, data=data,
        dpi=args.dpi)


def _global_score_page(data, out_dir, args):
    cells = data["cells"]
    frames = []
    for lib in data["libs_present"]:
        scores = load_library_scores(lib, args)
        if scores is None:
            continue
        sub = cells[cells["lib_num"] == lib].copy()
        m = sub.merge(scores, on=["lib_num", "barcode"], how="inner",
                      suffixes=("", "_diag"))
        if len(m):
            frames.append(m)
    if not frames:
        return None
    merged = pd.concat(frames, ignore_index=True)
    wrong = merged[merged["decision_group"] != "KEEP"]
    present = [(c, lab, sc) for c, lab, sc in SCORE_SPECS
               if c in merged.columns]
    if not present:
        return None
    n_libs = merged["lib_num"].nunique()
    return _render_score_grid(
        merged, wrong, present,
        title="All libraries: demux confidence, whole set vs flagged cells",
        subtitle=f"{len(merged):,} cells across {n_libs} libraries; "
                 f"{len(wrong):,} flagged (reassigned or held)",
        out_name="LS_all_libraries_scores.png", out_dir=out_dir, data=data,
        dpi=args.dpi)


# =============================================================================
# A1  Roster-status matrix
# =============================================================================

def figure_a1(data, out_dir, args):
    libs = data["libs_present"]
    per_lib = {lib: library_status_table(data, lib) for lib in libs}
    identities = set()
    for tab in per_lib.values():
        identities.update(tab)

    def ident_key(g):
        libcount = sum(1 for tab in per_lib.values() if g in tab)
        return (-libcount, str(g))
    identities = sorted(identities, key=ident_key)
    cap = 90
    shown = identities[:cap]
    n_ident, n_lib = len(shown), len(libs)

    code = {"absent": 0, "expected_present": 1, "component_present": 2,
            "expected_absent": 3, "unexpected_present": 4}
    grid = np.zeros((n_ident, n_lib))
    for j, lib in enumerate(libs):
        tab = per_lib[lib]
        for i, g in enumerate(shown):
            grid[i, j] = code[tab.get(g, ("absent", 0))[0]]

    cmap = ListedColormap([STATUS_COLORS["absent"],
                           STATUS_COLORS["expected_present"],
                           STATUS_COLORS["component_present"],
                           STATUS_COLORS["expected_absent"],
                           STATUS_COLORS["unexpected_present"]])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], cmap.N)

    row_h = max(0.17, min(0.30, 22.0 / max(n_ident, 1)))
    fig_h = max(7.5, n_ident * row_h + 3.4)
    fig_w = max(9.0, n_lib * 0.36 + 4.5)
    fig = plt.figure(figsize=(fig_w, fig_h))
    top_frac = 1 - 1.7 / fig_h
    bot_frac = 2.2 / fig_h
    gs = GridSpec(1, 1, left=0.30, right=0.965, top=top_frac, bottom=bot_frac)
    ax = fig.add_subplot(gs[0, 0])
    ax.imshow(grid, cmap=cmap, norm=norm, aspect="auto",
              interpolation="nearest")
    ax.set_xticks(range(n_lib))
    ax.set_xticklabels([str(l) for l in libs], fontsize=7, rotation=90)
    ax.set_yticks(range(n_ident))
    ax.set_yticklabels([wrap_identity(g, 22) for g in shown],
                       fontsize=max(4.5, min(7, 340 / max(n_ident, 1))))
    draw_batch_rules(ax, libs, "x")
    ax.set_xlabel("Library", color=PAL["text_dark"])
    ax.set_title("A1. Identity presence across every library",
                 fontweight="bold", pad=12, color=PAL["text_dark"])
    handles = [Patch(facecolor=STATUS_COLORS[k], edgecolor="#999999",
                     label=STATUS_PRETTY[k]) for k in
               ("expected_present", "component_present", "expected_absent",
                "unexpected_present")]
    handles.append(Patch(facecolor=STATUS_COLORS["absent"],
                         edgecolor="#999999", label="Not in roster, absent"))
    ax.legend(handles=handles, loc="upper center",
              bbox_to_anchor=(0.5, -1.3 / (fig_h * (top_frac - bot_frac))),
              ncol=4, frameon=False, fontsize=8.5)
    tail = "" if len(identities) <= cap else \
        f"  Showing {cap} of {len(identities)} identities (rarest omitted)."
    fig.suptitle("Which lines were present, unexpected, or missing everywhere?",
                 fontweight="bold", fontsize=15, y=0.997)
    fig.text(0.5, 1 - 1.0 / fig_h,
             "Teal marks exact expected lines, blue residual component singlets, "
             "amber expected-but-absent lines, and coral truly unexpected identities."
             + tail, ha="center", fontsize=10,
             color=PAL["footnote"])
    fig_footer(fig, data)
    return [save_fig(fig, out_dir, "A1_roster_status_matrix.png", args.dpi)]


# =============================================================================
# A2  Per-library outcome composition
# =============================================================================

def figure_a2(data, out_dir, args):
    cells = data["cells"]
    libs = data["libs_present"]
    frac = {g: [] for g in DECISION_ORDER}
    nonkeep = []
    for lib in libs:
        sub = cells[cells["lib_num"] == lib]
        vc = sub["decision_group"].value_counts()
        for g in DECISION_ORDER:
            frac[g].append(vc.get(g, 0) / max(len(sub), 1))
        nonkeep.append(int(len(sub) - vc.get("KEEP", 0)))

    fig = plt.figure(figsize=(max(9.0, len(libs) * 0.36 + 3), 8.6))
    gs = GridSpec(2, 1, height_ratios=[3, 1.15], left=0.09, right=0.975,
                  top=0.85, bottom=0.245, hspace=0.34)
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(libs))
    bottom = np.zeros(len(libs))
    for g in DECISION_ORDER:
        vals = np.array(frac[g])
        if vals.sum() == 0:
            continue
        ax.bar(x, vals, bottom=bottom, color=DECISION_COLORS[g],
               edgecolor="white", linewidth=0.3, label=DECISION_PRETTY[g],
               width=0.86)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([str(l) for l in libs], fontsize=7, rotation=90)
    ax.set_ylim(0, 1)
    draw_batch_rules(ax, libs, "x")
    style_ax(ax, title="A2. Reconciliation outcome composition per library",
             ylabel="Fraction of cells", grid_axis="y")
    handles = [Patch(facecolor=DECISION_COLORS[g], edgecolor="white",
                     label=DECISION_PRETTY[g]) for g in DECISION_ORDER
               if np.array(frac[g]).sum() > 0]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, 0.135), ncol=5, frameon=False,
               fontsize=8.5)

    ax2 = fig.add_subplot(gs[1, 0])
    ax2.bar(x, nonkeep, color=PAL["composite"], width=0.86)
    ax2.set_xticks(x)
    ax2.set_xticklabels([str(l) for l in libs], fontsize=7, rotation=90)
    ax2.set_yscale("symlog")
    draw_batch_rules(ax2, libs, "x")
    style_ax(ax2, ylabel="Non-kept cells", xlabel="Library", grid_axis="y")

    fig.suptitle("Where did reconciliation change or hold cells?",
                 fontweight="bold", fontsize=15, y=0.95)
    fig.text(0.5, 0.90, "Outcome fractions (top) and absolute non-kept counts "
             "(bottom, symlog); batch rules separate 1/2/3", ha="center",
             fontsize=10, color=PAL["footnote"])
    fig_footer(fig, data)
    return [save_fig(fig, out_dir, "A2_library_outcome_bars.png", args.dpi)]


# =============================================================================
# A3  Expected-recovery bars
# =============================================================================

def figure_a3(data, out_dir, args):
    libs = data["libs_present"]
    found, component, absent, unexpected = [], [], [], []
    for lib in libs:
        tab = library_status_table(data, lib)
        found.append(sum(1 for s, _ in tab.values()
                         if s == "expected_present"))
        component.append(sum(1 for s, _ in tab.values()
                             if s == "component_present"))
        absent.append(sum(1 for s, _ in tab.values()
                          if s == "expected_absent"))
        unexpected.append(sum(1 for s, _ in tab.values()
                              if s == "unexpected_present"))

    fig = plt.figure(figsize=(max(9.0, len(libs) * 0.38 + 3), 7.8))
    gs = GridSpec(1, 1, left=0.09, right=0.975, top=0.82, bottom=0.30)
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(libs))
    w = 0.20
    ax.bar(x - 1.5 * w, found, width=w, color=STATUS_COLORS["expected_present"],
           label="Expected exact & present", edgecolor="black", linewidth=0.3)
    ax.bar(x - 0.5 * w, component, width=w,
           color=STATUS_COLORS["component_present"],
           label="Residual component singlet", edgecolor="black", linewidth=0.3)
    ax.bar(x + 0.5 * w, absent, width=w, color=STATUS_COLORS["expected_absent"],
           label="Expected exact, absent", edgecolor="black", linewidth=0.3)
    ax.bar(x + 1.5 * w, unexpected, width=w,
           color=STATUS_COLORS["unexpected_present"],
           label="Truly unexpected, present", edgecolor="black", linewidth=0.3)
    ax.set_xticks(x)
    ax.set_xticklabels([str(l) for l in libs], fontsize=7, rotation=90)
    draw_batch_rules(ax, libs, "x")
    style_ax(ax, title="A3. Expected-line recovery per library",
             xlabel="Library", ylabel="Identity count", grid_axis="y")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=4,
              frameon=False, fontsize=8.5)
    help_text(ax, [
        "For every library: exact expected lines found/missing, residual donor "
        "singlets that are components of expected composites, and genuinely "
        "foreign identities. Component singlets are not counted as unexpected "
        "under the v11 biological expectedness rules.",
    ], y=-0.22)
    fig.suptitle("Did each library contain the lines it was supposed to?",
                 fontweight="bold", fontsize=15, y=0.95)
    fig_footer(fig, data)
    return [save_fig(fig, out_dir, "A3_expected_recovery_bars.png", args.dpi)]


# =============================================================================
# A4  Score-separation summary
# =============================================================================

def figure_a4(data, out_dir, args):
    cells = data["cells"]
    libs = data["libs_present"]
    rows = []
    for lib in libs:
        scores = load_library_scores(lib, args)
        if scores is None or "llr_vs_runner_up" not in scores.columns:
            continue
        sub = cells[cells["lib_num"] == lib].copy()
        m = sub.merge(scores[["lib_num", "barcode", "llr_vs_runner_up"]],
                      on=["lib_num", "barcode"], how="inner")
        if len(m) == 0:
            continue
        allv = pd.to_numeric(m["llr_vs_runner_up"], errors="coerce").dropna()
        wrong = m[m["decision_group"] != "KEEP"]
        wrongv = pd.to_numeric(wrong["llr_vs_runner_up"],
                               errors="coerce").dropna()
        rows.append({"lib": lib,
                     "all_med": float(allv.median()) if len(allv) else np.nan,
                     "wrong_med": float(wrongv.median()) if len(wrongv)
                     else np.nan,
                     "n_wrong": len(wrongv)})
    if not rows:
        print("  A4 skipped: no demux diagnostics under --demux-root")
        return []
    df = pd.DataFrame(rows).sort_values("lib").reset_index(drop=True)

    fig = plt.figure(figsize=(max(9.0, len(df) * 0.42 + 3), 6.8))
    gs = GridSpec(1, 1, left=0.09, right=0.975, top=0.81, bottom=0.30)
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(df))
    ax.plot(x, df["all_med"], marker="o", color=PAL["neutral"],
            label="Whole-library median LLR", linewidth=1.4, zorder=2)
    wm = df["wrong_med"].to_numpy()
    good = ~np.isnan(wm)
    ax.scatter(x[good], wm[good], color=INCORRECT_COLOR, s=48, zorder=3,
               label="Flagged-cell median LLR", edgecolor="white",
               linewidth=0.6)
    for xi, r in zip(x, df.itertuples()):
        if not np.isnan(r.wrong_med):
            ax.plot([xi, xi], [r.all_med, r.wrong_med], color="#CCCCCC",
                    linewidth=0.8, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels([str(int(l)) for l in df["lib"]], fontsize=7,
                       rotation=90)
    ax.set_yscale("symlog")
    style_ax(ax, title="A4. Demux confidence: library median vs flagged cells",
             xlabel="Library", ylabel="Median demux LLR (symlog)",
             grid_axis="y")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2,
              frameon=False, fontsize=9)
    help_text(ax, [
        "For each library with demux diagnostics: whole-library median demux "
        "LLR (grey) vs the median LLR of cells reconciliation reassigned or "
        "held (red). A red point near or above the grey line means the flagged "
        "cells were confidently demuxed yet wrong on identity, which is the "
        "case reconciliation exists to catch.",
    ], y=-0.235)
    fig.suptitle("Were the flagged cells low-confidence, or confidently wrong?",
                 fontweight="bold", fontsize=15, y=0.95)
    fig_footer(fig, data)
    return [save_fig(fig, out_dir, "A4_score_separation_summary.png",
                     args.dpi)]


# =============================================================================
# A5  Unexpected-present catalogue
# =============================================================================

def figure_a5(data, out_dir, args):
    cells = data["cells"]
    libs = data["libs_present"]
    recs = []
    for lib in libs:
        tab = library_status_table(data, lib)
        sub = cells[cells["lib_num"] == lib]
        # v11 residual component singlets are present-but-not-exact-roster, but
        # they are biologically represented in an expected composite and must
        # never be plotted as foreign/unexpected.
        component = sub[sub["is_component_derived_final"]]
        for g, grp in component.groupby("reconciled_donor_genotype"):
            recs.append({"lib": lib, "identity": str(g), "cells": len(grp),
                         "hetero": False, "kind": "component"})
        # genuinely unexpected-present bubbles (sized by cells)
        unexp = sub[sub["is_unexpected_final"]]
        for g, grp in unexp.groupby("reconciled_donor_genotype"):
            recs.append({"lib": lib, "identity": str(g), "cells": len(grp),
                         "hetero": component_count(g) >= 2, "kind": "unexpected"})
        # expected-but-absent markers (open circle, no size)
        for g, (status_val, _) in tab.items():
            if status_val == "expected_absent":
                recs.append({"lib": lib, "identity": str(g), "cells": 0,
                             "hetero": component_count(g) >= 2,
                             "kind": "missing"})
    if not recs:
        fig = plt.figure(figsize=(9, 4))
        ax = fig.add_subplot(111)
        ax.axis("off")
        ax.text(0.5, 0.5, "No residual-component, unexpected, or missing identities "
                "in the selected libraries.", ha="center", va="center", fontsize=12,
                color=PAL["good"], fontweight="bold")
        fig.suptitle("Residual / unexpected / missing identity catalogue",
                     fontweight="bold", fontsize=15)
        fig_footer(fig, data)
        return [save_fig(fig, out_dir, "A5_unexpected_catalogue.png",
                         args.dpi)]
    df = pd.DataFrame(recs)
    meta_map = data.get("meta_map")
    # Sort identities by CorrectedFZGRP then UID, diploids as a separate block
    # before tetraploids (verified from Library_conversions.xlsx 2025_LineMeta).
    idents = sorted(dict.fromkeys(df["identity"]),
                    key=lambda g: identity_sort_key(meta_map, g))
    id_idx = {g: i for i, g in enumerate(idents)}
    m = len(idents)

    fig_h = max(6.0, m * 0.36 + 3.2)
    fig = plt.figure(figsize=(max(9.0, len(libs) * 0.32 + 5), fig_h))
    gs = GridSpec(1, 1, left=0.30, right=0.90, top=1 - 1.5 / fig_h,
                  bottom=1.6 / fig_h)
    ax = fig.add_subplot(gs[0, 0])
    lib_idx = {l: j for j, l in enumerate(libs)}
    for r in df.itertuples():
        x = lib_idx[r.lib]
        yy = id_idx[r.identity]
        if r.kind == "missing":
            # expected but not found: open circle, dashed red-black outline,
            # fixed size (there is no cell count to encode)
            ax.scatter(x, yy, s=95, facecolors="none",
                       edgecolors=STATUS_COLORS["expected_absent"],
                       linewidths=1.6, marker="o", zorder=3)
            ax.scatter(x, yy, s=95, facecolors="none", edgecolors="#B00020",
                       linewidths=0.7, linestyle=(0, (2, 2)), marker="o",
                       zorder=4)
        else:
            size = 30 + np.sqrt(r.cells) * 22
            if r.kind == "component":
                color = STATUS_COLORS["component_present"]
            else:
                color = PAL["composite"] if r.hetero else STATUS_COLORS["unexpected_present"]
            ax.scatter(x, yy, s=size, color=color, alpha=0.82,
                       edgecolor="white", linewidth=0.6)
            if r.cells >= 30:
                ax.text(x, yy, str(int(r.cells)), ha="center", va="center",
                        fontsize=5.6, color="white")
    ax.set_xticks(range(len(libs)))
    ax.set_xticklabels([str(l) for l in libs], fontsize=7, rotation=90,
                       fontweight="bold")
    ax.set_yticks(range(m))
    ax.set_yticklabels([wrap_identity(g, 22) for g in idents],
                       fontsize=max(5, min(7.5, 300 / max(m, 1))),
                       fontweight="bold")
    ax.set_ylim(-0.7, m - 0.3)
    ax.set_xlim(-0.7, len(libs) - 0.3)
    draw_batch_rules(ax, libs, "x")
    prev_block, prev_fz = None, None
    for i, g in enumerate(idents):
        mk = genotype_meta(meta_map, g)
        block = 0 if mk["is_diploid"] else 1
        if prev_block is not None and block != prev_block:
            ax.axhline(i - 0.5, color=PAL["text_dark"], linewidth=1.3)
        elif prev_fz is not None and mk["fzgrp"] != prev_fz:
            ax.axhline(i - 0.5, color=PAL["grid"], linewidth=0.6)
        prev_block, prev_fz = block, mk["fzgrp"]
    style_ax(ax, xlabel="Library", grid_axis=None)
    ax.xaxis.label.set_fontweight("bold")
    handles = [
        Line2D([0], [0], marker="o", linestyle="",
               color=STATUS_COLORS["component_present"],
               label="Residual expected-composite component singlet", markersize=8),
        Line2D([0], [0], marker="o", linestyle="",
               color=STATUS_COLORS["unexpected_present"],
               label="Truly unexpected singlet", markersize=8),
        Line2D([0], [0], marker="o", linestyle="", color=PAL["composite"],
               label="Unexpected heterotypic line", markersize=8),
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor="none",
               markeredgecolor="#B00020", label="Expected, not found",
               markersize=9, markeredgewidth=1.4),
        Line2D([0], [0], marker="o", linestyle="", color=PAL["neutral"],
               label="size = cells", markersize=11, alpha=0.5),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              frameon=False, borderaxespad=0, fontsize=8.5)
    fig.suptitle("Residual components, unexpected identities, and missing lines",
                 fontweight="bold", fontsize=15, y=1 - 0.4 / fig_h)
    fig.text(0.30, 1 - 1.05 / fig_h,
             "Blue bubbles are residual/component-derived singlets already represented "
             "inside expected composite lines; coral/purple bubbles are truly unexpected "
             "singlets or heterotypic lines (size = cells). Open red rings are exact "
             "expected lines found in no cell. Rows sort by CorrectedFZGRP then UID.",
             ha="left",
             fontsize=9.5, color=PAL["footnote"])
    fig_footer(fig, data)
    return [save_fig(fig, out_dir, "A5_unexpected_catalogue.png", args.dpi)]


# =============================================================================
# A6  Full-category identity catalogue (expected + unexpected, singlet + doublet)
# =============================================================================

# Five v11 presence categories: exact expected singlet/fusion, residual
# component-derived singlet, and truly unexpected singlet/fusion.
A6_CATEGORY_COLORS = {
    "expected_singlet": "#2A9D8F",      # teal
    "component_singlet": "#2E86AB",     # blue
    "expected_doublet": "#1F6F78",      # deep teal
    "unexpected_singlet": "#E76F51",    # coral
    "unexpected_doublet": "#6C5B7B",    # purple
}
A6_CATEGORY_PRETTY = {
    "expected_singlet": "Expected singlet",
    "component_singlet": "Residual expected-composite component singlet",
    "expected_doublet": "Expected doublet / fusion",
    "unexpected_singlet": "Truly unexpected singlet",
    "unexpected_doublet": "Unexpected doublet / fusion",
}


def _a6_category(relationship, genotype):
    hetero = component_count(genotype) >= 2
    if relationship == "component":
        return "component_singlet"
    if relationship == "unexpected":
        return "unexpected_doublet" if hetero else "unexpected_singlet"
    return "expected_doublet" if hetero else "expected_singlet"


def figure_a6(data, out_dir, args):
    cells = data["cells"]
    libs = data["libs_present"]
    meta_map = data.get("meta_map")
    recs = []
    for lib in libs:
        sub = cells[cells["lib_num"] == lib]
        for g, grp in sub.groupby("reconciled_donor_genotype"):
            gs = str(g).strip()
            if not gs or gs in ("NA", ".", "nan"):
                continue
            relationship = str(grp["final_identity_relationship"].iloc[0])
            recs.append({"lib": lib, "identity": gs, "cells": len(grp),
                         "category": _a6_category(relationship, gs),
                         "missing": False})
        # expected-but-absent: open ring, no cell count. This makes A6 a full
        # catalogue so the A1 status matrix is no longer needed.
        tab = library_status_table(data, lib)
        for g, (status_val, _) in tab.items():
            if status_val == "expected_absent":
                recs.append({"lib": lib, "identity": str(g), "cells": 0,
                             "category": "missing", "missing": True})
    if not recs:
        fig = plt.figure(figsize=(9, 4))
        ax = fig.add_subplot(111)
        ax.axis("off")
        ax.text(0.5, 0.5, "No identities in the selected libraries.",
                ha="center", va="center", fontsize=12, color=PAL["good"],
                fontweight="bold")
        fig.suptitle("Full identity catalogue", fontweight="bold", fontsize=15)
        fig_footer(fig, data)
        return [save_fig(fig, out_dir, "A6_full_identity_catalogue.png",
                         args.dpi)]
    df = pd.DataFrame(recs)
    idents = sorted(dict.fromkeys(df["identity"]),
                    key=lambda g: identity_sort_key(meta_map, g))
    id_idx = {g: i for i, g in enumerate(idents)}
    m = len(idents)
    lib_idx = {l: j for j, l in enumerate(libs)}

    fig_h = max(6.5, m * 0.34 + 3.4)
    fig = plt.figure(figsize=(max(10.0, len(libs) * 0.34 + 5.5), fig_h))
    gs = GridSpec(1, 1, left=0.32, right=0.86, top=1 - 1.5 / fig_h,
                  bottom=1.7 / fig_h)
    ax = fig.add_subplot(gs[0, 0])
    # Alternating row shading so the eye can track a genotype across 40 columns.
    for i in range(m):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color="#00000008", zorder=0)
    for r in df.itertuples():
        x = lib_idx[r.lib]
        yy = id_idx[r.identity]
        if r.missing:
            ax.scatter(x, yy, s=70, facecolors="none",
                       edgecolors=STATUS_COLORS["expected_absent"],
                       linewidths=1.5, marker="o", zorder=3)
            ax.scatter(x, yy, s=70, facecolors="none", edgecolors="#B00020",
                       linewidths=0.7, linestyle=(0, (2, 2)), marker="o",
                       zorder=4)
        else:
            size = 26 + np.sqrt(r.cells) * 20
            ax.scatter(x, yy, s=size, color=A6_CATEGORY_COLORS[r.category],
                       alpha=0.88, edgecolor="white", linewidth=0.6, zorder=2)
            if r.cells >= 60:
                ax.text(x, yy, str(int(r.cells)), ha="center", va="center",
                        fontsize=5.4, color="white", zorder=3)
    ax.set_xticks(range(len(libs)))
    ax.set_xticklabels([str(l) for l in libs], fontsize=7.5, rotation=90,
                       fontweight="bold")
    ax.set_yticks(range(m))
    ax.set_yticklabels([wrap_identity(g, 22) for g in idents],
                       fontsize=max(4.8, min(7.5, 300 / max(m, 1))),
                       fontweight="bold")
    ax.set_ylim(-0.7, m - 0.3)
    ax.set_xlim(-0.7, len(libs) - 0.3)
    draw_batch_rules(ax, libs, "x")
    prev_block, prev_fz = None, None
    for i, g in enumerate(idents):
        mk = genotype_meta(meta_map, g)
        block = 0 if mk["is_diploid"] else 1
        if prev_block is not None and block != prev_block:
            ax.axhline(i - 0.5, color=PAL["text_dark"], linewidth=1.4, zorder=1)
        elif prev_fz is not None and mk["fzgrp"] != prev_fz:
            ax.axhline(i - 0.5, color=PAL["grid"], linewidth=0.6, zorder=1)
        prev_block, prev_fz = block, mk["fzgrp"]
    style_ax(ax, xlabel="Library", grid_axis=None)
    ax.xaxis.label.set_fontweight("bold")
    handles = [Line2D([0], [0], marker="o", linestyle="",
                      color=A6_CATEGORY_COLORS[k], label=A6_CATEGORY_PRETTY[k],
                      markersize=8) for k in
               ("expected_singlet", "component_singlet", "expected_doublet",
                "unexpected_singlet", "unexpected_doublet")]
    handles.append(Line2D([0], [0], marker="o", linestyle="",
                          markerfacecolor="none", markeredgecolor="#B00020",
                          label="Expected, not found", markersize=9,
                          markeredgewidth=1.4))
    handles.append(Line2D([0], [0], marker="o", linestyle="",
                          color=PAL["neutral"], label="size = cells",
                          markersize=11, alpha=0.5))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              frameon=False, borderaxespad=0, fontsize=8.5)
    fig.suptitle("Every final identity in every library, by category",
                 fontweight="bold", fontsize=15, y=1 - 0.4 / fig_h)
    fig.text(0.32, 1 - 1.05 / fig_h,
             "All reconciled identities with v11 biological expectedness: exact "
             "expected lines, residual component-derived singlets, and genuinely "
             "unexpected singlets/fusions (filled, size = cells), plus exact "
             "expected-but-not-found lines as open red rings. Rows sort by "
             "CorrectedFZGRP then UID.",
             ha="left", fontsize=9.5, color=PAL["footnote"])
    fig_footer(fig, data)
    return [save_fig(fig, out_dir, "A6_full_identity_catalogue.png", args.dpi)]


# =============================================================================
# A7  Library-exchange matrix (v11 distinguishing-roster evidence)
# =============================================================================

EXCHANGE_INTERPRETATION_PRETTY = {
    "ONE_WAY_FOREIGN_SIGNATURE": "One-way foreign signature",
    "LIKELY_CROSS_CONTAMINATION": "Likely cross-contamination",
    "LIKELY_PARTIAL_LIBRARY_REPLACEMENT": "Likely partial replacement",
    "RECIPROCAL_LIBRARY_MIXING": "Reciprocal library mixing",
    "LIKELY_RECIPROCAL_LIBRARY_SWAP": "Likely reciprocal library swap",
}


def figure_a7(data, out_dir, args):
    exchange = data.get("exchange")
    libs = data["libs_present"]
    if exchange is None or len(exchange) == 0:
        print("  A7 skipped: v11 all_libraries.library_exchange_events.tsv not found")
        return []
    if not {"lib_a_num", "lib_b_num", "a_signature_coverage_in_b",
            "b_signature_coverage_in_a"}.issubset(exchange.columns):
        print("  A7 skipped: library-exchange table lacks required v11 columns")
        return []

    lib_idx = {lib: i for i, lib in enumerate(libs)}
    n = len(libs)
    coverage = np.full((n, n), np.nan, dtype=float)
    equivalent = np.zeros((n, n), dtype=bool)
    np.fill_diagonal(equivalent, True)

    for row in exchange.itertuples(index=False):
        a = int(getattr(row, "lib_a_num"))
        b = int(getattr(row, "lib_b_num"))
        if a not in lib_idx or b not in lib_idx:
            continue
        ia, ib = lib_idx[a], lib_idx[b]
        relation = str(getattr(row, "roster_relation", ""))
        if relation == "ROSTER_EQUIVALENT_NONDISCRIMINATING":
            equivalent[ia, ib] = True
            equivalent[ib, ia] = True
            continue
        # Row = source library whose distinguishing donor signature is being
        # tracked; column = target library where that signature was observed.
        av = getattr(row, "a_signature_coverage_in_b")
        bv = getattr(row, "b_signature_coverage_in_a")
        if pd.notna(av):
            coverage[ia, ib] = float(av)
        if pd.notna(bv):
            coverage[ib, ia] = float(bv)

    fig = plt.figure(figsize=(17.0, 10.5))
    gs = GridSpec(1, 2, width_ratios=[1.45, 1.0], left=0.07, right=0.98,
                  top=0.86, bottom=0.16, wspace=0.22)
    ax = fig.add_subplot(gs[0, 0])
    cmap = plt.get_cmap("YlOrRd").copy()
    cmap.set_bad(STATUS_COLORS["absent"])
    im = ax.imshow(coverage, vmin=0, vmax=1, cmap=cmap, aspect="equal",
                   interpolation="nearest")
    ax.set_xticks(range(n))
    ax.set_xticklabels([str(x) for x in libs], fontsize=7, rotation=90,
                       fontweight="bold")
    ax.set_yticks(range(n))
    ax.set_yticklabels([str(x) for x in libs], fontsize=7,
                       fontweight="bold")
    draw_batch_rules(ax, libs, "x")
    draw_batch_rules(ax, libs, "y")
    ax.set_xlabel("Target library where source-specific singlets were observed",
                  fontweight="bold", color=PAL["text_dark"])
    ax.set_ylabel("Source library defining the distinguishing donor signature",
                  fontweight="bold", color=PAL["text_dark"])
    ax.set_title("A7. Directed distinguishing-roster signature coverage",
                 fontweight="bold", pad=12, color=PAL["text_dark"])
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_label("Fraction of source-specific donors detected in target",
                   fontsize=9)

    # Mark roster-equivalent/non-identifiable pairs without pretending their
    # missing directed coverage is a biological zero.
    if n <= 24:
        for i in range(n):
            for j in range(n):
                if i != j and equivalent[i, j]:
                    ax.text(j, i, "ND", ha="center", va="center", fontsize=5.5,
                            color=PAL["footnote"], fontweight="bold")

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.axis("off")
    candidates = exchange[~exchange["exchange_interpretation"].isin(
        ["NO_LIBRARY_EXCHANGE_SIGNAL", "ROSTER_EQUIVALENT_NONDISCRIMINATING"]
    )].copy()
    if len(candidates):
        conf_rank = {"STRONG": 0, "MODERATE": 1, "SUGGESTIVE": 2,
                     "NONDISCRIMINATING": 3, "NONE": 4}
        interp_rank = {"LIKELY_RECIPROCAL_LIBRARY_SWAP": 0,
                       "RECIPROCAL_LIBRARY_MIXING": 1,
                       "LIKELY_PARTIAL_LIBRARY_REPLACEMENT": 2,
                       "LIKELY_CROSS_CONTAMINATION": 3,
                       "ONE_WAY_FOREIGN_SIGNATURE": 4}
        candidates["_conf"] = candidates["exchange_confidence"].map(
            conf_rank).fillna(9)
        candidates["_interp"] = candidates["exchange_interpretation"].map(
            interp_rank).fillna(9)
        candidates["_coverage"] = candidates[
            ["a_signature_coverage_in_b", "b_signature_coverage_in_a"]
        ].max(axis=1, skipna=True).fillna(0)
        candidates = candidates.sort_values(
            ["_interp", "_conf", "_coverage"], ascending=[True, True, False])
        rows = []
        for r in candidates.head(14).itertuples(index=False):
            pair = f"{int(r.lib_a_num)} ↔ {int(r.lib_b_num)}"
            interp = EXCHANGE_INTERPRETATION_PRETTY.get(
                str(r.exchange_interpretation), str(r.exchange_interpretation))
            acov = r.a_signature_coverage_in_b
            bcov = r.b_signature_coverage_in_a
            directional = (
                f"{float(acov):.0%} / {float(bcov):.0%}"
                if pd.notna(acov) and pd.notna(bcov)
                else f"{float(acov):.0%} / NA" if pd.notna(acov)
                else f"NA / {float(bcov):.0%}" if pd.notna(bcov)
                else "NA / NA"
            )
            rows.append([pair, interp, str(r.exchange_confidence), directional])
        tab = ax2.table(cellText=rows,
                        colLabels=["Library pair", "Interpretation",
                                   "Confidence", "A→B / B→A"],
                        loc="upper center", cellLoc="left",
                        colWidths=[0.18, 0.43, 0.18, 0.21],
                        colColours=[PAL["bg_light"]] * 4)
        tab.auto_set_font_size(False)
        tab.set_fontsize(7.5)
        tab.scale(1.0, 1.35)
        ax2.set_title("Exchange candidates", fontweight="bold", fontsize=12,
                      pad=12, color=PAL["text_dark"])
        if len(candidates) > 14:
            ax2.text(0.5, 0.03, f"+ {len(candidates) - 14} additional candidate pair(s)",
                     transform=ax2.transAxes, ha="center", va="bottom",
                     fontsize=8, color=PAL["footnote"])
    else:
        ax2.set_title("Exchange candidates", fontweight="bold", fontsize=12,
                      pad=12, color=PAL["text_dark"])
        ax2.text(0.5, 0.55, "No library-exchange candidates in this selection",
                 transform=ax2.transAxes, ha="center", va="center",
                 fontsize=11, color=PAL["good"], fontweight="bold")

    fig.suptitle("Do distinguishing donor signatures point to another library?",
                 fontweight="bold", fontsize=15, y=0.96)
    fig.text(0.5, 0.905,
             "Shared expected donors contribute zero evidence. Heatmap cells show only "
             "the fraction of source-library-specific donors recovered as robust truly "
             "unexpected singlets in the target; roster-equivalent pairs are non-identifiable.",
             ha="center", fontsize=9.5, color=PAL["footnote"])
    fig_footer(fig, data,
               extra="A7 uses v11 library_exchange_events; shared donors are excluded by construction.")
    return [save_fig(fig, out_dir, "A7_library_exchange_matrix.png", args.dpi)]


# =============================================================================
# AS  Ambient swap-test figures (G/H fixed profile; J/K joint profile)
# =============================================================================

AMBIENT_STRATUM_ORDER = ("changed_target", "scrutinized_other", "background")
AMBIENT_STRATUM_PRETTY = {
    "changed_target": "changed target",
    "scrutinized_other": "scrutinized other",
    "background": "background",
}
AMBIENT_STRATUM_COLORS = {
    "changed_target": "#CC6677",
    "scrutinized_other": "#4477AA",
    "background": "#BBBBBB",
}


def _ambient_slug(value):
    text = str(value).strip().replace("+", "_plus_")
    text = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in text)
    text = "_".join(part for part in text.split("_") if part)
    return text or "unnamed"


def _ambient_read_rate(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise RuntimeError(f"ambient swap-test rate file missing or empty: {path}")
    values = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            fields = raw.split()
            if len(fields) < 2 or fields[0].lower() in {"barcode", "cell"}:
                continue
            try:
                value = float(fields[1])
            except ValueError:
                continue
            if np.isfinite(value):
                values[fields[0]] = value
    if not values:
        raise RuntimeError(f"ambient swap-test rate file has no finite rows: {path}")
    return values


def _ambient_read_rows(path, required):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise RuntimeError(f"ambient swap-test table missing or empty: {path}")
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or ())
        if missing:
            raise RuntimeError(
                f"ambient swap-test table {path} lacks {sorted(missing)}")
        rows = list(reader)
    return rows


def _ambient_plot_sample(rows, maximum=AMBIENT_PLOT_MAX_POINTS):
    """Deterministically thin plotted points without changing summaries."""
    if len(rows) <= maximum:
        return rows
    positions = np.linspace(0, len(rows) - 1, maximum, dtype=int)
    return [rows[int(pos)] for pos in positions]


def _ambient_kde_curve(values, grid):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.zeros_like(grid)
    if values.size > AMBIENT_PLOT_MAX_POINTS:
        positions = np.linspace(
            0, values.size - 1, AMBIENT_PLOT_MAX_POINTS, dtype=int)
        values = np.sort(values)[positions]
    spread = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
    bandwidth = max(0.005, 1.06 * spread * values.size ** (-0.2))
    density = np.zeros_like(grid)
    normalizer = bandwidth * np.sqrt(2.0 * np.pi) * values.size
    for start in range(0, values.size, 2000):
        chunk = values[start:start + 2000]
        scaled = (grid[:, None] - chunk[None, :]) / bandwidth
        density += np.exp(-0.5 * scaled * scaled).sum(axis=1)
    return density / normalizer


def _ambient_summarize(rows, left, right):
    left_values = np.asarray([row[left] for row in rows], dtype=float)
    right_values = np.asarray([row[right] for row in rows], dtype=float)
    delta = right_values - left_values
    return {
        "n_paired": int(len(rows)),
        "median_left_c": float(np.median(left_values)),
        "median_right_c": float(np.median(right_values)),
        "median_delta_c": float(np.median(delta)),
        "q25_delta_c": float(np.quantile(delta, 0.25)),
        "q75_delta_c": float(np.quantile(delta, 0.75)),
        "fraction_right_lower": float(np.mean(delta < 0)),
        "fraction_equal": float(np.mean(delta == 0)),
    }


def _ambient_save_pair(fig, out_dir, stem, dpi, manifest, metadata):
    os.makedirs(out_dir, exist_ok=True)
    written = []
    for suffix in ("pdf", "png"):
        path = os.path.join(out_dir, f"{stem}.{suffix}")
        fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
        print(f"  wrote {path}")
        manifest.append({**metadata, "format": suffix, "path": path})
        written.append(path)
    plt.close(fig)
    return written


def _ambient_scope_rows(records, candidate_identity=None):
    """Return one clean three-stratum view from a library/condition slice.

    The combined scope uses every selected candidate cell and every source
    control. An identity scope keeps only that identity's candidate cells,
    derives its source controls from their original identities, and excludes
    other candidate populations instead of mislabeling them as background.
    """
    if candidate_identity is None:
        return [(row, row["combined_stratum"]) for row in records]
    target = [
        row for row in records
        if row["candidate_identity"] == candidate_identity]
    source_identities = {row["original_identity"] for row in target}
    scoped = [(row, "changed_target") for row in target]
    for row in records:
        if row["combined_stratum"] == "changed_target":
            continue
        stratum = (
            "scrutinized_other"
            if row["original_identity"] in source_identities
            else "background")
        scoped.append((row, stratum))
    return scoped


def _ambient_render_three_figure_set(
        scoped, out_dir, title, condition, dpi, manifest, metadata,
        shared_joint_profile, filename_prefix):
    scoped_rows = [row for row, _ in scoped]
    if not scoped_rows:
        return [], []
    by_stratum = {
        stratum: [row for row, observed in scoped if observed == stratum]
        for stratum in AMBIENT_STRATUM_ORDER
    }
    by_stratum = {key: value for key, value in by_stratum.items() if value}
    if "changed_target" not in by_stratum:
        return [], []

    summaries = []
    for stratum, rows in list(by_stratum.items()) + [("all", scoped_rows)]:
        for comparison, left, right in (
                ("joint_profile", "joint_original_c", "joint_proposed_c"),
                ("fixed_empty", "fixed_original_c", "fixed_proposed_c")):
            summaries.append({
                **metadata,
                "stratum": stratum,
                "comparison": comparison,
                **_ambient_summarize(rows, left, right),
            })

    written = []
    fig, ax = plt.subplots(figsize=(7.4, 6.8))
    for stratum in AMBIENT_STRATUM_ORDER[::-1]:
        rows = by_stratum.get(stratum, [])
        if not rows:
            continue
        plotted = _ambient_plot_sample(rows)
        ax.scatter(
            [row["fixed_original_c"] for row in plotted],
            [row["fixed_proposed_c"] for row in plotted],
            s=7 if stratum == "background" else 17,
            alpha=0.18 if stratum == "background" else 0.62,
            color=AMBIENT_STRATUM_COLORS[stratum],
            label=f"{AMBIENT_STRATUM_PRETTY[stratum]} (n={len(rows):,})")
    all_fixed = [
        row[key] for row in scoped_rows
        for key in ("fixed_original_c", "fixed_proposed_c")]
    upper = min(1.0, max(0.10, float(np.quantile(all_fixed, 0.997)) * 1.04))
    ax.plot([0, upper], [0, upper], color="#333333", linewidth=0.8,
            linestyle="--")
    ax.set_xlim(0, upper)
    ax.set_ylim(0, upper)
    style_ax(
        ax, title=f"{title}\n{condition}",
        xlabel="Original assignments, frozen empty-drop profile (c)",
        ylabel="Proposed assignments, same frozen profile (c)",
        grid_axis=None)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    written.extend(_ambient_save_pair(
        fig, out_dir,
        f"{filename_prefix}__fixed_empty_assignment_scatter", dpi, manifest,
        {**metadata, "figure_kind": "fixed_empty_assignment_scatter"}))

    all_values = [
        row[key] for row in scoped_rows
        for key in ("fixed_original_c", "fixed_proposed_c")]
    kde_upper = min(
        1.0, max(0.10, float(np.quantile(all_values, 0.995)) * 1.10))
    grid = np.linspace(0.0, kde_upper, 320)
    kde_panels = [
        (value, by_stratum[value])
        for value in AMBIENT_STRATUM_ORDER if value in by_stratum]
    kde_panels.append(("all", scoped_rows))
    fig, axes = plt.subplots(
        len(kde_panels), 1, sharex=True, squeeze=False,
        figsize=(9.2, 2.55 * len(kde_panels) + 1.7))
    for ax, (stratum, rows) in zip(axes.ravel(), kde_panels):
        left = [row["fixed_original_c"] for row in rows]
        right = [row["fixed_proposed_c"] for row in rows]
        ax.plot(grid, _ambient_kde_curve(left, grid), color="#4477AA",
                linewidth=1.8, label="Original assignments / fixed profile")
        ax.plot(grid, _ambient_kde_curve(right, grid), color="#CC6677",
                linewidth=1.8, label="Proposed assignments / same profile")
        panel_title = (
            "All paired cells in scope"
            if stratum == "all" else AMBIENT_STRATUM_PRETTY[stratum])
        style_ax(
            ax,
            title=f"{panel_title} (n={len(rows):,})",
            ylabel="Density", grid_axis="both", title_pad=9)
    axes[0, 0].legend(frameon=False, fontsize=8, ncol=2)
    axes[-1, 0].set_xlabel("Estimated contamination fraction (c)")
    fig.suptitle(f"{title}: assignment-only fixed-profile KDE\n{condition}",
                 fontsize=12.5, fontweight="bold", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    written.extend(_ambient_save_pair(
        fig, out_dir,
        f"{filename_prefix}__fixed_empty_kde_by_stratum", dpi, manifest,
        {**metadata, "figure_kind": "fixed_empty_kde_by_stratum"}))

    strata = [value for value in AMBIENT_STRATUM_ORDER if value in by_stratum]
    x = np.arange(len(strata))
    width = 0.36
    joint = [
        _ambient_summarize(
            by_stratum[value], "joint_original_c", "joint_proposed_c")
        ["median_delta_c"] for value in strata]
    fixed = [
        _ambient_summarize(
            by_stratum[value], "fixed_original_c", "fixed_proposed_c")
        ["median_delta_c"] for value in strata]
    fig, ax = plt.subplots(figsize=(8.8, 5.4))
    ax.bar(x - width / 2, joint, width, color="#CCBB44",
           label="Joint-profile original → proposed")
    ax.bar(x + width / 2, fixed, width, color="#228833",
           label="Fixed-profile original → proposed")
    ax.axhline(0, color="#333333", linewidth=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels([AMBIENT_STRATUM_PRETTY[value] for value in strata])
    style_ax(
        ax, title=f"{title}: joint versus assignment-only contrast\n{condition}",
        ylabel="Median paired Δc (proposed − original)", grid_axis="y")
    ax.legend(frameon=False, fontsize=8)
    if shared_joint_profile:
        help_text(ax, [
            "The per-identity cell subset is isolated, while the joint-profile fit",
            "is shared by all candidates in this library; no redundant per-candidate refit."
        ], y=-0.23)
        fig.subplots_adjust(bottom=0.20)
    else:
        fig.tight_layout()
    written.extend(_ambient_save_pair(
        fig, out_dir,
        f"{filename_prefix}__joint_vs_fixed_delta", dpi, manifest,
        {**metadata, "figure_kind": "joint_vs_fixed_delta"}))
    return written, summaries


def render_ambient_swap_figures(spec_path, output_dir, dpi=200):
    """Render full-run ambient swap analysis from an orchestrator manifest."""
    with open(spec_path, "r", encoding="utf-8") as handle:
        spec = json.load(handle)
    condition_short_names = spec.get("condition_short_names", {})
    libraries = spec.get("libraries", [])
    if not libraries:
        raise RuntimeError("ambient swap-test manifest has no applicable libraries")
    os.makedirs(output_dir, exist_ok=True)
    data_dir = os.path.join(output_dir, "data")
    figure_root = os.path.join(output_dir, "figures")
    os.makedirs(data_dir, exist_ok=True)
    # This directory is wholly owned by this rendering mode. Rebuild it so an
    # older nested by-identity layout cannot remain beside the current figures.
    if os.path.lexists(figure_root):
        if os.path.islink(figure_root) or os.path.isfile(figure_root):
            raise RuntimeError(
                f"ambient swap figure root is not a directory: {figure_root}")
        shutil.rmtree(figure_root)
    os.makedirs(figure_root, exist_ok=True)

    records = []
    candidate_audit = []
    discovery_audit = []
    for library_spec in libraries:
        library = int(library_spec["library"])
        candidates = _ambient_read_rows(
            library_spec["candidate_cells"],
            {"barcode", "original_identity", "candidate_identity", "event_id"})
        candidate_by_barcode = {}
        for row in candidates:
            barcode = row["barcode"]
            if barcode in candidate_by_barcode:
                raise RuntimeError(
                    f"lib{library}: duplicate candidate barcode {barcode}")
            candidate_by_barcode[barcode] = row
        strata_rows = _ambient_read_rows(
            library_spec["cell_strata"],
            {"barcode", "original_identity", "candidate_identity", "stratum"})
        strata = {row["barcode"]: row for row in strata_rows}
        if len(strata) != len(strata_rows):
            raise RuntimeError(f"lib{library}: duplicate cell-strata barcode")
        candidate_audit.extend({"library": f"lib{library}", **row}
                               for row in candidates)
        event_discovery = library_spec.get("event_discovery")
        if event_discovery:
            discovery_audit.extend(_ambient_read_rows(
                event_discovery,
                {"library", "candidate_identity", "event_id"}))

        for condition, paths in sorted(library_spec["conditions"].items()):
            required_arms = {"G", "H", "J", "K"}
            missing_arms = required_arms - set(paths)
            if missing_arms:
                raise RuntimeError(
                    f"lib{library} {condition}: manifest lacks arms "
                    f"{sorted(missing_arms)}")
            rates = {arm: _ambient_read_rate(paths[arm])
                     for arm in sorted(required_arms)}
            common = sorted(set.intersection(
                *(set(values) for values in rates.values())))
            if not common:
                raise RuntimeError(
                    f"lib{library} {condition}: G/H/J/K have no paired cells")
            for barcode in common:
                candidate = candidate_by_barcode.get(barcode)
                stratum_row = strata.get(barcode, {})
                combined = stratum_row.get("stratum", "background")
                combined = {
                    "candidate_cell": "changed_target",
                    "source_control": "scrutinized_other",
                    "background": "background",
                }.get(combined, "background")
                records.append({
                    "library": f"lib{library}",
                    "library_num": library,
                    "condition": condition,
                    "barcode": barcode,
                    "original_identity": stratum_row.get(
                        "original_identity",
                        candidate.get("original_identity", "NA")
                        if candidate else "NA"),
                    "candidate_identity": (
                        candidate.get("candidate_identity", "NA")
                        if candidate else "NA"),
                    "event_id": (
                        candidate.get("event_id", "NA")
                        if candidate else "NA"),
                    "combined_stratum": combined,
                    "fixed_original_c": rates["G"][barcode],
                    "fixed_proposed_c": rates["H"][barcode],
                    "joint_original_c": rates["J"][barcode],
                    "joint_proposed_c": rates["K"][barcode],
                })

    if not records:
        raise RuntimeError("ambient swap-test manifest produced no paired records")
    records_df = pd.DataFrame(records)
    records_df["fixed_delta_c"] = (
        records_df["fixed_proposed_c"] - records_df["fixed_original_c"])
    records_df["joint_delta_c"] = (
        records_df["joint_proposed_c"] - records_df["joint_original_c"])
    records_df.to_csv(
        os.path.join(data_dir, "ambient_swap_cell_contrasts.tsv.gz"),
        sep="\t", index=False, compression="gzip")
    pd.DataFrame(candidate_audit).to_csv(
        os.path.join(data_dir, "ambient_swap_candidate_cells.tsv"),
        sep="\t", index=False)
    if discovery_audit:
        pd.DataFrame(discovery_audit).to_csv(
            os.path.join(data_dir, "ambient_swap_event_discovery_audit.tsv"),
            sep="\t", index=False)
    if spec.get("applicability"):
        pd.DataFrame(spec["applicability"]).to_csv(
            os.path.join(data_dir, "ambient_swap_library_applicability.tsv"),
            sep="\t", index=False)

    manifest = []
    summaries = []
    written = []
    conditions = sorted(records_df["condition"].unique())
    for condition in conditions:
        condition_label = condition_short_names.get(condition, condition)
        condition_slug = _ambient_slug(condition_label)
        condition_records = [
            row for row in records if row["condition"] == condition]
        scoped = _ambient_scope_rows(condition_records)
        out_dir = os.path.join(figure_root, "all_libraries_combined")
        files, rows = _ambient_render_three_figure_set(
            scoped, out_dir, "All applicable libraries | combined candidates",
            condition_label, dpi, manifest,
            {"scope": "all_libraries_combined", "library": "all",
             "candidate_identity": "all", "condition": condition,
             "condition_short_name": condition_label},
            shared_joint_profile=False,
            filename_prefix=f"all_libraries__{condition_slug}")
        written.extend(files)
        summaries.extend(rows)

    for library in sorted(records_df["library_num"].unique()):
        library_name = f"lib{int(library)}"
        lib_records = [row for row in records if row["library_num"] == library]
        candidates = sorted({
            row["candidate_identity"] for row in lib_records
            if row["candidate_identity"] not in {"", "NA"}})
        for condition in sorted({row["condition"] for row in lib_records}):
            condition_label = condition_short_names.get(condition, condition)
            condition_slug = _ambient_slug(condition_label)
            condition_records = [
                row for row in lib_records if row["condition"] == condition]
            library_dir = os.path.join(
                figure_root, "by_library", library_name)
            files, rows = _ambient_render_three_figure_set(
                _ambient_scope_rows(condition_records), library_dir,
                f"{library_name} | combined candidates", condition_label, dpi,
                manifest,
                {"scope": "library_combined", "library": library_name,
                 "candidate_identity": "all", "condition": condition,
                 "condition_short_name": condition_label},
                shared_joint_profile=False,
                filename_prefix=(
                    f"{library_name}__all_candidates__{condition_slug}"))
            written.extend(files)
            summaries.extend(rows)
            for identity in candidates:
                files, rows = _ambient_render_three_figure_set(
                    _ambient_scope_rows(condition_records, identity),
                    library_dir, f"{library_name} | {identity}",
                    condition_label,
                    dpi, manifest,
                    {"scope": "identity", "library": library_name,
                     "candidate_identity": identity,
                     "condition": condition,
                     "condition_short_name": condition_label},
                    shared_joint_profile=(len(candidates) > 1),
                    filename_prefix=(
                        f"{library_name}__{_ambient_slug(identity)}__"
                        f"{condition_slug}"))
                written.extend(files)
                summaries.extend(rows)

    pd.DataFrame(summaries).to_csv(
        os.path.join(data_dir, "ambient_swap_group_summaries.tsv"),
        sep="\t", index=False)
    pd.DataFrame(manifest).to_csv(
        os.path.join(data_dir, "ambient_swap_figure_manifest.tsv"),
        sep="\t", index=False)
    summary = {
        "schema_version": 1,
        "producer": "identity_reconciliation_figures.py",
        "input_spec": os.path.abspath(spec_path),
        "identity_event_path": spec.get("identity_event_path", "NA"),
        "reconciliation_manifest": spec.get(
            "reconciliation_manifest", "NA"),
        "discovery_id": spec.get("discovery_id", "NA"),
        "applicable_libraries": len({row["library"] for row in records}),
        "not_applicable_libraries": sum(
            row.get("applicable") != "yes"
            for row in spec.get("applicability", ())),
        "conditions": conditions,
        "condition_short_names": {
            condition: condition_short_names.get(condition, condition)
            for condition in conditions
        },
        "candidate_identities": sorted({
            row["candidate_identity"] for row in records
            if row["candidate_identity"] not in {"", "NA"}}),
        "paired_cells": len(records),
        "figure_files": len(written),
        "figure_layout": (
            "flat_per_library_with_scope_identity_condition_filenames"),
        "individual_candidate_refits": False,
        "interpretation": (
            "G/H hold one candidate-roster empty-drop profile fixed and differ "
            "only in original versus proposed assignments. J/K use the same "
            "receiver and donor rosters but fit their profiles independently. "
            "Per-identity figures isolate that candidate's cells and source "
            "controls; J/K remain the efficient library-wide combined fit."),
    }
    with open(os.path.join(output_dir, "ambient_swap_figure_summary.json"),
              "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"Ambient swap figures complete: {len(written)} file(s) under {figure_root}")
    return written


# =============================================================================
# Registry
# =============================================================================

FIGURES = {
    "LP": (figure_lp, "Per-library everything page (one file per library)"),
    "LS": (figure_ls, "Per-library + global demux score pages"),
    "A2": (figure_a2, "Per-library outcome composition, all libs"),
    "A3": (figure_a3, "Expected-line recovery per library, all libs"),
    "A4": (figure_a4, "Demux confidence: library vs flagged cells, all libs"),
    "A5": (figure_a5, "Residual component + unexpected + missing catalogue"),
    "A6": (figure_a6, "Full v11-aware identity catalogue + missing"),
    "A7": (figure_a7, "Library-exchange matrix from distinguishing-roster evidence"),
}


# =============================================================================
# CLI
# =============================================================================

def parse_lib_spec(tokens):
    out = []
    for tok in tokens:
        if "-" in tok:
            a, b = tok.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(tok))
    return sorted(set(out))


def parse_args():
    ap = argparse.ArgumentParser(
        description="Identity-reconciliation figure generator (v11 outputs)")
    ap.add_argument("--reconciliation-root",
                    default="/mnt/beegfs/tetmultiome_rna_mapped/"
                            "ploidy_classifier/retrain_nomito_20260814/"
                            "identity_reconciliation")
    ap.add_argument("--demux-root",
                    default="/mnt/beegfs/tetmultiome_rna_mapped/mapping_output",
                    help="Root with Tet_2025_Multiome-RNA_{N}/demux_nomito for "
                         "LS/A4 score pages; optional")
    ap.add_argument("--audit-root",
                    default="/mnt/beegfs/tetmultiome_rna_mapped/"
                            "ploidy_classifier/retrain_nomito_20260814/"
                            "swap_audit")
    ap.add_argument("--refined-root",
                    default="/mnt/beegfs/tetmultiome_rna_mapped/"
                            "ploidy_classifier/retrain_nomito_20260814/"
                            "tetra_refine")
    ap.add_argument("--schema", default=None,
                    help="Optional harvest JSON to bind real column names")
    ap.add_argument("--line-metadata",
                    default="/mnt/beegfs/tetmultiome_rna_mapped/mapping_output/"
                            "Library_conversions.xlsx",
                    help="Library_conversions.xlsx; identities are sorted by "
                         "CorrectedFZGRP then UID (diploids separate from "
                         "tetraploids) when present")
    ap.add_argument("--output-dir", default=None)
    ap.add_argument("--figures", nargs="+", default=["all"])
    ap.add_argument("--libraries", nargs="+", default=None)
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--list", action="store_true")
    ap.add_argument(
        "--ambient-swap-spec", default=None, metavar="JSON",
        help=("Render the standard G/H fixed-profile and J/K joint-profile "
              "swap-test figure sets from an orchestrator manifest. This mode "
              "does not load the normal reconciliation figure inputs."))
    return ap.parse_args()


def main():
    args = parse_args()
    if args.list:
        for tag, (_, desc) in FIGURES.items():
            print(f"  {tag:3s} {desc}")
        print("  AS  Ambient swap-test combined/per-library/per-identity figure sets")
        return 0
    if not args.output_dir:
        sys.exit("ERROR: --output-dir is required to render figures")
    if args.ambient_swap_spec:
        render_ambient_swap_figures(
            os.path.abspath(args.ambient_swap_spec),
            os.path.abspath(args.output_dir), args.dpi)
        return 0
    tags = list(FIGURES) if args.figures == ["all"] else args.figures
    unknown = [t for t in tags if t not in FIGURES]
    if unknown:
        sys.exit(f"ERROR: unknown figure tag(s): {', '.join(unknown)}; "
                 f"available: {', '.join(FIGURES)}")
    libraries = parse_lib_spec(args.libraries) if args.libraries else None
    os.makedirs(args.output_dir, exist_ok=True)
    data = load_data(os.path.abspath(args.reconciliation_root), libraries,
                     args.schema, args.line_metadata)
    written = []
    for tag in tags:
        fn, desc = FIGURES[tag]
        print(f"[{tag}] {desc}")
        written.extend(fn(data, args.output_dir, args))
    print(f"Done: {len(written)} figure file(s) in {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())


# =============================================================================
# Revision history
# =============================================================================
# V2_R12 - Added an all-paired-cells panel to every existing fixed-profile KDE
#          while preserving the three stratum panels, figure paths, manifests,
#          and numerical summaries. The pooled panel compares the same paired
#          G/H assignment-only contamination estimates represented by the
#          mutually exclusive strata in its plotting scope.
# V2_R11 - Flattened ambient swap-test output. All-library plots now share one
#          all_libraries_combined directory, and each library has one directory
#          containing both combined-candidate and per-identity plots. Scope,
#          identity, short condition, and plot kind are encoded in filenames.
#          Rendering replaces its owned figures directory so the retired deep
#          by_identity/<identity>/<condition> tree cannot survive regeneration.
# V2_R10 - Ambient swap-test output folders and figure titles use short
#          condition names supplied by the orchestrator while preserving the
#          canonical condition key in every analysis table and summary.
# V2_R9 - Added production ambient swap-test figure ownership. The script reads
#         orchestrator G/H fixed-profile and J/K independently fitted-profile
#         manifests and writes the three standard diagnostics as all-library,
#         per-library combined-candidate, and per-library per-identity sets.
#         Per-identity views isolate the candidate and its source controls while
#         retaining one efficient library-wide J/K fit. Removed the embedded
#         synthetic test-slice generator from the deployed production script.
# V2_R8 - v11 biological-expectedness update. Exact expected genotypes, residual
#         expected-composite component singlets, and truly unexpected identities
#         are now separate plotting categories. LP/A3/A5/A6 no longer color a
#         donor such as C3651 as foreign merely because only C3651+partner is
#         explicitly expected. Footer policy text is read from the data instead
#         of hardcoding v10. Loader reads the v11 library_exchange_events table.
#         New A7 shows directed distinguishing-roster signature coverage and a
#         concise exchange-candidate table; roster-equivalent pairs remain
#         non-identifiable and shared donors contribute no plotted evidence.
# V2_R7 - A5 now shows expected-but-missing identities as open red dashed
#         rings alongside the unexpected-present bubbles; A6 does the same and
#         becomes the full catalogue (all four present categories plus missing
#         rings), retiring A1 from the default set. A6/A5 axis labels bolded;
#         A6 gained alternating row shading. LP pages now scale their height to
#         the identity count (fixes the crushed/overlapping labels on
#         many-line libraries like lib5), with extra row height when labels
#         wrap to two lines. LS score pages rebuilt as a clean interleaved
#         2-column grid matching the reference quality, with a new global
#         all-library LS page (LS_all_libraries_scores.png). Every figure gates
#         to zero for the full run, single-library, and two-library subsets.
# V2_R6 - Identity sorting keyed to Library_conversions.xlsx (2025_LineMeta):
#         rows now sort by CorrectedFZGRP then UID with diploids as a block
#         above tetraploids, so fusion partners cluster (verified: WGS_Key maps
#         to one CorrectedFZGRP, reversed A+B orientations resolve to the same
#         group, diploids are single tokens with Tet.Class Diploid). A5 gained
#         a diploid/tetraploid divider and FZGRP group hairlines. New figure A6
#         (full_identity_catalogue) mirrors A5 but shows ALL categories as
#         distinct colors: expected singlet, expected doublet/fusion,
#         unexpected singlet, unexpected doublet/fusion. New --line-metadata
#         arg (defaults to the cluster xlsx); when the workbook is absent,
#         sorting falls back to a component-count diploid/tetraploid split then
#         genotype name, so figures still render.
# V2_R5 - Loader now gathers the 40 per-library lib{N}.reconciled_cells.tsv.gz
#         files. The consolidated all_libraries.reconciled_cells.tsv.gz in this
#         run is misnamed and holds only one library; the loader detects that
#         (uses the combined file only when it truly spans >1 library) and
#         otherwise concatenates every per-library file, so all 40 libraries
#         load. --libraries reads only the requested per-library files rather
#         than loading all then filtering. Per-library identity_events files
#         are gathered the same way.
# V2_R4 - CRITICAL join fix: barcodes are NOT unique across libraries (~5%
#         collide), so LS/A4 now key the demux join on (lib_num, barcode), not
#         barcode alone. Per-library diagnostics are tagged with lib_num and
#         deduped on barcode; the LS merge asserts the row count cannot inflate
#         and aborts rather than emit a corrupted figure. Verified by injecting
#         a cross-library barcode collision and confirming isolation.
# V2_R3 - Hardcoded against the verified production schema (harvest V1_R2).
#         Reconciled and demux barcodes are identical raw 16-mers (overlap
#         8042/8042 on lib19), so the join is a plain barcode equality; removed
#         the -1 suffix stripping, the match_rate computation, and the
#         skip-if-thin fallback (no runtime hedging). Added the fifth real
#         final_action class REVIEW_UNEXPECTED_IDENTITY (441 cells on lib19)
#         as its own decision group, color, and legend entry; held counts now
#         derive from DECISION_ORDER so a new class cannot be silently dropped.
# V2_R2 - LS score pages now print the per-library demux join rate during the
#         run and warn loudly (⚠️) when a library joins under 90%, so the
#         harvester barcode spot-check is no longer needed before rendering.
# V2_R1 - Full redesign after V1 delta-figures broke on a real lib19-only run.
#   - New organizing question: for every library, which lines were present,
#     present-but-unexpected, expected-but-absent, and how whole-library demux
#     confidence compares to the flagged cells.
#   - Per-library "everything" pages (LP) modeled on compare_demux_results_V5:
#     outcome donut, roster-detection bar, all-identity cell-count bars colored
#     by presence status, unexpected/absent tables.
#   - Per-library demux score pages (LS) modeled on tetra_score_distributions:
#     whole library vs reassigned/held cells, real barcode join (trailing -1
#     stripped), graceful skip when diagnostics absent.
#   - Five all-library aggregates (A1 status matrix, A2 outcome bars, A3
#     recovery, A4 score separation, A5 unexpected catalogue) robust to any
#     library subset; they no longer collapse to one stretched row.
#   - Library-label normalization handles Tet_2025_Multiome-RNA_{N}, lib{N},
#     {N}, {0N}; unparseable rows dropped with a warning.
#   - Optional --schema binds real column names from the harvest JSON; missing
#     columns degrade instead of crashing.
